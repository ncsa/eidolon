//! NEAT-aware false-negative attribution.
//!
//! Each FN that survives exact-match classification + equivalence sweep is
//! tagged with one or more reasons drawn from the simulator's own
//! configuration:
//!
//!   - [`Reason::OutsideSimulatedContigs`] — the FN's contig wasn't in the
//!     simulator's `contigs_simulated` list. Reported alone; the BED checks
//!     are skipped because they presuppose simulation.
//!   - [`Reason::OutsideMutationBed`]      — `mutation_bed` was configured
//!     and the FN's position falls outside its regions.
//!   - [`Reason::OutsideTargetBed`]        — `target_bed` was configured
//!     and the FN's position falls outside its regions.
//!   - [`Reason::Unknown`]                 — none of the above; the
//!     simulator can't explain it.
//!
//! Port of NEAT 4's `compare_vcfs/attribution.py` adapted to eidolon's
//! existing BedRecord/HashMap shape.
use eidolon_core::structs::{bed_record::BedRecord, variants::Variant};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Reason {
    OutsideSimulatedContigs,
    OutsideMutationBed,
    OutsideTargetBed,
    Unknown,
}

impl Reason {
    pub fn as_str(self) -> &'static str {
        match self {
            Reason::OutsideSimulatedContigs => "outside_simulated_contigs",
            Reason::OutsideMutationBed => "outside_mutation_bed",
            Reason::OutsideTargetBed => "outside_target_bed",
            Reason::Unknown => "unknown",
        }
    }
}

/// Tag one FN. Mirrors NEAT 4's `attribute_fn` semantics: if the contig
/// wasn't simulated, that root cause stands alone; otherwise the configured
/// BEDs are checked and `Unknown` is the fallback.
pub fn attribute_fn(
    chrom: &str,
    pos: usize,
    contigs_simulated: Option<&HashSet<String>>,
    mutation_intervals: Option<&HashMap<String, Vec<BedRecord>>>,
    target_intervals: Option<&HashMap<String, Vec<BedRecord>>>,
) -> Vec<Reason> {
    if let Some(set) = contigs_simulated
        && !set.contains(chrom)
    {
        return vec![Reason::OutsideSimulatedContigs];
    }

    let mut reasons = Vec::new();
    if let Some(beds) = mutation_intervals
        && !position_in_intervals(chrom, pos, beds)
    {
        reasons.push(Reason::OutsideMutationBed);
    }
    if let Some(beds) = target_intervals
        && !position_in_intervals(chrom, pos, beds)
    {
        reasons.push(Reason::OutsideTargetBed);
    }
    if reasons.is_empty() {
        reasons.push(Reason::Unknown);
    }
    reasons
}

/// Tag every FN. Returns `(record, reasons)` pairs alongside the rolled-up
/// `{reason → count}` summary for the report.
pub fn attribute_fns(
    fns: &[(&str, &Variant)],
    contigs_simulated: Option<&HashSet<String>>,
    mutation_intervals: Option<&HashMap<String, Vec<BedRecord>>>,
    target_intervals: Option<&HashMap<String, Vec<BedRecord>>>,
) -> AttributionResult {
    let mut per_fn = Vec::with_capacity(fns.len());
    let mut counts: BTreeMap<Reason, usize> = BTreeMap::new();
    for (chrom, v) in fns {
        let reasons = attribute_fn(
            chrom,
            v.location,
            contigs_simulated,
            mutation_intervals,
            target_intervals,
        );
        for r in &reasons {
            *counts.entry(*r).or_insert(0) += 1;
        }
        per_fn.push(((*chrom).to_string(), (*v).clone(), reasons));
    }
    AttributionResult { per_fn, counts }
}

pub struct AttributionResult {
    pub per_fn: Vec<(String, Variant, Vec<Reason>)>,
    pub counts: BTreeMap<Reason, usize>,
}

/// `BedRecord::contains_vcf_pos` is O(N) over a contig's regions. For our small in-
/// repo cases this is fine; if we ever hit pathological BEDs we can swap in
/// a sorted-interval lookup later. Returns true when no regions are defined
/// for the contig only when the BED is empty for that contig — same default
/// as NEAT 4's "missing-contig means outside".
fn position_in_intervals(
    chrom: &str,
    pos: usize,
    intervals: &HashMap<String, Vec<BedRecord>>,
) -> bool {
    match intervals.get(chrom) {
        Some(records) => records.iter().any(|r| r.contains_vcf_pos(chrom, pos)),
        None => false,
    }
}

// ── chrom aliases ────────────────────────────────────────────────────────────

/// Load a two-column TSV mapping `bed_chrom -> reference_chrom`. Empty lines
/// and `#`-prefixed comments are skipped. Returns an empty map if `path` is
/// `None`.
pub fn load_chrom_aliases(
    path: Option<&std::path::Path>,
) -> Result<HashMap<String, String>, std::io::Error> {
    use std::fs;
    use std::io::{BufRead, BufReader};
    let mut out = HashMap::new();
    let Some(p) = path else { return Ok(out) };
    let file = fs::File::open(p)?;
    for (lineno, raw) in BufReader::new(file).lines().enumerate() {
        let line = raw?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 2 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "{}:{}: expected two tab-separated columns, got {} field(s)",
                    p.display(),
                    lineno + 1,
                    parts.len()
                ),
            ));
        }
        out.insert(parts[0].to_string(), parts[1].to_string());
    }
    Ok(out)
}

/// Apply the alias map to every key of `beds` in place. Used to normalize a
/// BED's chrom names to the reference's canonical names before attribution.
pub fn apply_aliases_to_beds(
    beds: &mut HashMap<String, Vec<BedRecord>>,
    aliases: &HashMap<String, String>,
) {
    if aliases.is_empty() {
        return;
    }
    let keys: Vec<String> = beds.keys().cloned().collect();
    for k in keys {
        if let Some(canonical) = aliases.get(&k)
            && canonical != &k
        {
            let v = beds.remove(&k).unwrap();
            // Rewrite the BedRecord's internal contig name so subsequent
            // `record.contains_vcf_pos(chrom, pos)` calls use the canonical name.
            // BedRecord::new_bed_record returns Result, but our start/end
            // already pass the start < end check (otherwise the original
            // BedRecord wouldn't exist), so unwrap is safe.
            let rewritten: Vec<BedRecord> = v
                .into_iter()
                .map(|r| BedRecord::new_bed_record(canonical.clone(), r.start, r.end).unwrap())
                .collect();
            beds.entry(canonical.clone()).or_default().extend(rewritten);
        }
    }
}

// ── chrom-name-mismatch warning ──────────────────────────────────────────────

#[derive(Debug, Clone, Serialize)]
pub struct ChromNamingWarning {
    pub bed_label: String,
    /// Every BED chrom (post-alias) that did NOT match any reference contig.
    /// On a clean run this list is empty and no warning is produced.
    pub unmatched_chroms: Vec<String>,
    /// Up to five reference contigs, alphabetical, for context in the warning.
    pub reference_chroms_sample: Vec<String>,
    pub message: String,
}

/// Detect BED rows whose chrom doesn't match any contig in the reference set.
///
/// Fires when **any** BED chrom (post-alias) is absent from the reference's
/// contigs — catches single-typo cases like `chr_MP` in an otherwise correct
/// H1N1 BED. NEAT 4's original semantics (warn only on *zero* overlap) made
/// the single-typo case silent because the seven correct rows would silence
/// the eighth mistyped one.
pub fn detect_chrom_naming_mismatch(
    bed_label: &str,
    beds: &HashMap<String, Vec<BedRecord>>,
    reference_chroms: &HashSet<String>,
) -> Option<ChromNamingWarning> {
    if beds.is_empty() || reference_chroms.is_empty() {
        return None;
    }
    let mut unmatched: Vec<String> = beds
        .keys()
        .filter(|c| !reference_chroms.contains(*c))
        .cloned()
        .collect();
    if unmatched.is_empty() {
        return None;
    }
    unmatched.sort();

    let mut ref_sample: Vec<String> = reference_chroms.iter().cloned().collect();
    ref_sample.sort();
    ref_sample.truncate(5);

    let n_unmatched = unmatched.len();
    let bed_chroms_count = beds.len();
    let preview = unmatched
        .iter()
        .take(5)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");
    let message = if n_unmatched == bed_chroms_count {
        format!(
            "{bed_label} chrom names don't overlap with the reference's \
             ({n_unmatched}/{bed_chroms_count} unmatched: {preview}). \
             Attribution against this BED will report every FN as \
             'outside_{bed_label}'. Consider passing a chrom_aliases TSV."
        )
    } else {
        format!(
            "{bed_label} has {n_unmatched}/{bed_chroms_count} unmatched chrom name(s): \
             {preview}. These rows will not contribute to attribution and any \
             variants in them will be misattributed."
        )
    };
    Some(ChromNamingWarning {
        bed_label: bed_label.to_string(),
        unmatched_chroms: unmatched,
        reference_chroms_sample: ref_sample,
        message,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::structs::variants::{AlternateType, Provenance};
    use eidolon_core::structs::{
        nucleotides::Nucleotide,
        variants::{Genotype, Variant, VariantType},
    };

    fn snp(loc: usize) -> Variant {
        Variant {
            variant_type: VariantType::SNP,
            location: loc,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Literal(vec![Nucleotide::G]),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: Some("PASS".to_string()),
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        }
    }

    fn bed(chrom: &str, start: usize, end: usize) -> HashMap<String, Vec<BedRecord>> {
        let mut h = HashMap::new();
        h.insert(
            chrom.to_string(),
            vec![BedRecord::new_bed_record(chrom.to_string(), start, end).unwrap()],
        );
        h
    }

    #[test]
    fn outside_simulated_contigs_is_reported_alone() {
        let simulated: HashSet<String> = ["chr1".to_string()].into_iter().collect();
        let mutation = bed("chr1", 0, 1_000_000);
        let target = bed("chr1", 0, 1_000_000);
        let reasons = attribute_fn(
            "chrZ",
            500,
            Some(&simulated),
            Some(&mutation),
            Some(&target),
        );
        // The BED checks are skipped — would have said outside both,
        // but the unsimulated-contig short-circuit wins.
        assert_eq!(reasons, vec![Reason::OutsideSimulatedContigs]);
    }

    #[test]
    fn outside_mutation_bed_alone() {
        let simulated: HashSet<String> = ["chr1".to_string()].into_iter().collect();
        let mutation = bed("chr1", 0, 100);
        let target = bed("chr1", 0, 1_000_000);
        let reasons = attribute_fn(
            "chr1",
            500,
            Some(&simulated),
            Some(&mutation),
            Some(&target),
        );
        assert_eq!(reasons, vec![Reason::OutsideMutationBed]);
    }

    #[test]
    fn outside_both_beds_reports_both() {
        let simulated: HashSet<String> = ["chr1".to_string()].into_iter().collect();
        let mutation = bed("chr1", 0, 100);
        let target = bed("chr1", 0, 100);
        let reasons = attribute_fn(
            "chr1",
            500,
            Some(&simulated),
            Some(&mutation),
            Some(&target),
        );
        assert_eq!(
            reasons,
            vec![Reason::OutsideMutationBed, Reason::OutsideTargetBed]
        );
    }

    #[test]
    fn fully_in_bounds_yields_unknown() {
        let simulated: HashSet<String> = ["chr1".to_string()].into_iter().collect();
        let mutation = bed("chr1", 0, 1_000_000);
        let target = bed("chr1", 0, 1_000_000);
        let reasons = attribute_fn(
            "chr1",
            500,
            Some(&simulated),
            Some(&mutation),
            Some(&target),
        );
        assert_eq!(reasons, vec![Reason::Unknown]);
    }

    /// #770: FN `pos` is the 1-based VCF POS. `chr1 100 200` covers
    /// POS 101..=200 for both the mutation and the target BED.
    #[test]
    fn bed_edges_use_one_based_pos() {
        let b = bed("chr1", 100, 200);
        for (pos, inside) in [(100, false), (101, true), (200, true), (201, false)] {
            let want = if inside {
                vec![Reason::Unknown]
            } else {
                vec![Reason::OutsideMutationBed, Reason::OutsideTargetBed]
            };
            assert_eq!(
                attribute_fn("chr1", pos, None, Some(&b), Some(&b)),
                want,
                "POS {pos}"
            );
        }
    }

    #[test]
    fn no_beds_configured_yields_unknown() {
        let reasons = attribute_fn("chr1", 500, None, None, None);
        assert_eq!(reasons, vec![Reason::Unknown]);
    }

    #[test]
    fn attribute_fns_rolls_up_counts() {
        let v_inside = snp(50);
        let v_outside = snp(500);
        let mutation = bed("chr1", 0, 100);
        let fns: Vec<(&str, &Variant)> = vec![("chr1", &v_inside), ("chr1", &v_outside)];
        let result = attribute_fns(&fns, None, Some(&mutation), None);
        assert_eq!(*result.counts.get(&Reason::Unknown).unwrap(), 1);
        assert_eq!(*result.counts.get(&Reason::OutsideMutationBed).unwrap(), 1);
    }

    #[test]
    fn chrom_aliases_loader_parses_two_columns() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "# header comment").unwrap();
        writeln!(f, "1\tchr1").unwrap();
        writeln!(f, "X\tchrX").unwrap();
        writeln!(f).unwrap();
        writeln!(f, "MT\tchrM").unwrap();
        let aliases = load_chrom_aliases(Some(f.path())).unwrap();
        assert_eq!(aliases.get("1"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("X"), Some(&"chrX".to_string()));
        assert_eq!(aliases.get("MT"), Some(&"chrM".to_string()));
        assert_eq!(aliases.len(), 3);
    }

    #[test]
    fn apply_aliases_remaps_bed_keys_and_record_contigs() {
        let mut beds: HashMap<String, Vec<BedRecord>> = HashMap::new();
        beds.insert(
            "1".to_string(),
            vec![BedRecord::new_bed_record("1".to_string(), 0, 100).unwrap()],
        );
        let aliases: HashMap<String, String> = [("1".to_string(), "chr1".to_string())]
            .into_iter()
            .collect();
        apply_aliases_to_beds(&mut beds, &aliases);
        assert!(beds.contains_key("chr1"));
        assert!(!beds.contains_key("1"));
        // Inner record's contig field is also remapped so `contains()` works.
        assert!(beds["chr1"].iter().any(|r| r.contains("chr1", 50)));
    }

    #[test]
    fn naming_mismatch_warning_when_zero_overlap() {
        let beds = bed("1", 0, 100);
        let reference: HashSet<String> = ["chr1".to_string(), "chr2".to_string()]
            .into_iter()
            .collect();
        let w = detect_chrom_naming_mismatch("target_bed", &beds, &reference).unwrap();
        assert_eq!(w.bed_label, "target_bed");
        assert_eq!(w.unmatched_chroms, vec!["1".to_string()]);
        // All-unmatched message variant mentions the aliases TSV.
        assert!(w.message.contains("chrom_aliases"));
    }

    /// Single typo: BED has 8 H1N1 contigs, one of which is mistyped. The
    /// old "any overlap silences" semantics swallowed this; the new
    /// "any unmatched fires" semantics catches it.
    #[test]
    fn naming_mismatch_warning_when_one_chrom_is_mistyped() {
        let mut beds = bed("H1N1_HA", 0, 100);
        for good in ["H1N1_NA", "H1N1_NP", "H1N1_PA"] {
            beds.insert(
                good.to_string(),
                vec![BedRecord::new_bed_record(good.to_string(), 0, 100).unwrap()],
            );
        }
        beds.insert(
            "chr_MP".to_string(),
            vec![BedRecord::new_bed_record("chr_MP".to_string(), 0, 100).unwrap()],
        );
        let reference: HashSet<String> = ["H1N1_HA", "H1N1_MP", "H1N1_NA", "H1N1_NP", "H1N1_PA"]
            .into_iter()
            .map(String::from)
            .collect();
        let w = detect_chrom_naming_mismatch("target_bed", &beds, &reference).unwrap();
        assert_eq!(w.unmatched_chroms, vec!["chr_MP".to_string()]);
        // Partial-mismatch message variant doesn't recommend aliases — it's
        // probably a typo, not a convention difference.
        assert!(w.message.contains("1/5 unmatched"));
    }

    #[test]
    fn naming_mismatch_silent_when_all_match() {
        let mut beds = bed("chr1", 0, 100);
        beds.insert(
            "chr2".to_string(),
            vec![BedRecord::new_bed_record("chr2".to_string(), 0, 100).unwrap()],
        );
        let reference: HashSet<String> = ["chr1".to_string(), "chr2".to_string()]
            .into_iter()
            .collect();
        assert!(detect_chrom_naming_mismatch("target_bed", &beds, &reference).is_none());
    }
}
