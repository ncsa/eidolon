//! Driver for `compare-vcfs`.
//!
//! Algorithm:
//!   1. Read golden + called VCFs via `eidolon_core::file_tools::vcf_tools::read_vcf`.
//!   2. For each VCF, drop variants the user has filtered out: outside the
//!      simulated contigs, outside the target BED, hom-ref, non-PASS — each
//!      counted into a [`SkipCounts`] for the report.
//!   3. Per contig, key surviving variants by `(position, ref, alt)` into a
//!      `HashSet<VariantKey>`. TP = golden ∩ called; FN = golden \ called;
//!      FP = called \ golden.
//!   4. Unless `fast: true`, stream the reference FASTA and for each contig
//!      with FPs, run [`equivalence::sweep`] to catch denotation-different
//!      alternates. FNs consumed by the sweep are promoted to TP; FPs
//!      consumed are removed from FP.
//!   5. Tag every surviving FN with one or more NEAT-aware attribution
//!      reasons (outside_simulated_contigs / outside_mutation_bed /
//!      outside_target_bed / unknown) via [`attribute_fns`].
//!   6. Aggregate per-contig and total counts, compute precision/recall/F1,
//!      write `comparison_summary.{json,txt}`, `FN_with_reasons.vcf`, and
//!      optionally `FP.vcf`.
use crate::compare_vcfs::{
    errors::CompareVcfsError,
    utils::{
        attribution::{
            ChromNamingWarning, apply_aliases_to_beds, attribute_fns, detect_chrom_naming_mismatch,
            load_chrom_aliases,
        },
        config::RunConfiguration,
        equivalence,
        report::{
            ComparisonSummary, ContigCounts, Metrics, SCHEMA_VERSION, SkipCounts, SummaryInputs,
            SummaryOutputs, write_json, write_txt,
        },
        vcf_writer::{write_fn_with_reasons, write_fp_vcf},
    },
};
use eidolon_core::{
    file_tools::{bed_reader::read_bed, fasta_stream::FastaStream, vcf_tools::read_vcf},
    structs::{bed_record::BedRecord, nucleotides::Nucleotide, variants::Variant},
};
use log::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;

/// Comparison key: the three columns that decide a true positive. Hashing on
/// `Variant` directly is wrong because its derived `Hash` includes
/// genotype/filter/info, none of which the truth-vs-called comparison cares about.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct VariantKey {
    location: usize,
    reference: Vec<Nucleotide>,
    alternate: Vec<Nucleotide>,
}

impl From<&Variant> for VariantKey {
    fn from(v: &Variant) -> Self {
        // filter_vcf removes symbolic ALTs (`<DEL>`, `<DUP>`, ...) into the
        // `symbolic` skip bucket before any VariantKey is built, so by the
        // time we get here the ALT must be literal. Assert it in debug builds
        // so a future regression fails loudly with context instead of a bare
        // `unwrap()` panic in release.
        debug_assert!(
            v.alternate.is_literal(),
            "symbolic ALT reached VariantKey::from at location {}",
            v.location
        );
        VariantKey {
            location: v.location,
            reference: v.reference.clone(),
            alternate: v.alternate.as_literal().unwrap().to_vec(),
        }
    }
}

/// Per-contig comparison buckets retained between exact-match classification
/// and the optional equivalence sweep. TP records are not retained — only
/// the count, since the sweep cannot promote to or out of TP.
struct ContigBuckets {
    tp_count: usize,
    /// FNs on this contig that the equivalence sweep promoted to TP. Already
    /// included in `tp_count`; tracked separately so per-contig and total
    /// rollups can surface how many matches were rescued by denotation-aware
    /// comparison.
    equivalents_promoted: usize,
    fns: Vec<Variant>,
    fps: Vec<Variant>,
}

pub fn runner(config: &RunConfiguration) -> Result<(), CompareVcfsError> {
    fs::create_dir_all(&config.output_dir)?;

    // Load chrom aliases once and apply them to every BED at load time so
    // downstream `contains(chrom, pos)` calls use reference-canonical names.
    let aliases = load_chrom_aliases(config.chrom_aliases.as_deref())?;

    let mut target_bed = match &config.target_bed {
        Some(p) => {
            info!("Reading target BED: {}", p.display());
            Some(read_bed(p, false)?)
        }
        None => None,
    };
    if let Some(t) = target_bed.as_mut() {
        apply_aliases_to_beds(t, &aliases);
    }
    let mut mutation_bed = match &config.mutation_bed {
        Some(p) => {
            info!("Reading mutation BED: {}", p.display());
            Some(read_bed(p, false)?)
        }
        None => None,
    };
    if let Some(m) = mutation_bed.as_mut() {
        apply_aliases_to_beds(m, &aliases);
    }
    let simulated_set = config
        .contigs_simulated
        .as_ref()
        .map(|v| v.iter().cloned().collect::<HashSet<_>>());

    info!("Reading golden VCF: {}", config.golden_vcf.display());
    let golden_raw = read_vcf(config.golden_vcf.clone())?;
    info!("Reading called VCF: {}", config.called_vcf.display());
    let called_raw = read_vcf(config.called_vcf.clone())?;

    let (golden_by_contig, skipped_golden) = filter_vcf(
        &golden_raw,
        target_bed.as_ref(),
        simulated_set.as_ref(),
        config.include_homs,
        config.include_filtered,
    );
    let (called_by_contig, skipped_called) = filter_vcf(
        &called_raw,
        target_bed.as_ref(),
        simulated_set.as_ref(),
        config.include_homs,
        config.include_filtered,
    );

    let mut buckets = classify(&golden_by_contig, &called_by_contig);
    let raw_tp: usize = buckets.values().map(|b| b.tp_count).sum();
    let raw_fn: usize = buckets.values().map(|b| b.fns.len()).sum();
    let raw_fp: usize = buckets.values().map(|b| b.fps.len()).sum();
    info!("Exact-match classification: TP={raw_tp} FN={raw_fn} FP={raw_fp}");

    if !config.fast {
        run_equivalence_sweep(&mut buckets, &config.reference, config.equivalence_window)?;
    } else {
        info!("Skipping equivalence sweep (fast mode)");
    }

    let (per_contig, totals) = aggregate(&buckets);
    let metrics = Metrics::from_counts(totals.tp, totals.fn_, totals.fp);
    info!(
        "Final classification: TP={} (promoted={}) FN={} FP={}",
        totals.tp, totals.equivalents_promoted, totals.fn_, totals.fp
    );

    // Attribute every remaining FN (post-sweep). Contigs in the bucket map
    // line up with the surviving FNs after `swap_remove` in
    // `run_equivalence_sweep`.
    let mut fn_pairs: Vec<(&str, &Variant)> = Vec::new();
    for (chrom, b) in &buckets {
        for v in &b.fns {
            fn_pairs.push((chrom.as_str(), v));
        }
    }
    let attribution_result = attribute_fns(
        &fn_pairs,
        simulated_set.as_ref(),
        mutation_bed.as_ref(),
        target_bed.as_ref(),
    );
    if !attribution_result.counts.is_empty() {
        let display: Vec<String> = attribution_result
            .counts
            .iter()
            .map(|(r, c)| format!("{}={c}", r.as_str()))
            .collect();
        info!("FN attribution: {}", display.join(", "));
    }

    // Detect BED-vs-reference chrom-naming mismatches. The reference set
    // should be the FASTA's full contig list, not the VCFs' — a BED can
    // legitimately cover contigs that have no variants in the golden VCF
    // (e.g., simulator ran on 8 contigs but only put variants in 1).
    // Falls back to the VCF contig union if FASTA scanning fails for any
    // reason, so the warning still has a chance to fire.
    let reference_chroms: HashSet<String> =
        match eidolon_core::file_tools::fasta_stream::scan_fasta_lengths(&config.reference) {
            Ok(pairs) => pairs.into_iter().map(|(name, _)| name).collect(),
            Err(e) => {
                warn!(
                    "Could not scan reference FASTA contigs ({e}); falling back to VCF \
                 contig union for chrom-naming-mismatch detection. Some warnings \
                 may be missed."
                );
                derive_reference_chroms(&golden_raw, &called_raw)
            }
        };
    // Diagnostic logging for the chrom-naming-mismatch check: emitted at
    // debug level so users can `--log-level debug` to see exactly what the
    // detector is comparing when a warning fails to fire (or fires
    // unexpectedly).
    if log_enabled!(log::Level::Debug) {
        let mut ref_sorted: Vec<&String> = reference_chroms.iter().collect();
        ref_sorted.sort();
        debug!(
            "Reference contigs ({}): [{}]",
            reference_chroms.len(),
            ref_sorted
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
        if let Some(t) = target_bed.as_ref() {
            let mut t_sorted: Vec<&String> = t.keys().collect();
            t_sorted.sort();
            debug!(
                "target_bed contigs ({}): [{}]",
                t.len(),
                t_sorted
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        if let Some(m) = mutation_bed.as_ref() {
            let mut m_sorted: Vec<&String> = m.keys().collect();
            m_sorted.sort();
            debug!(
                "mutation_bed contigs ({}): [{}]",
                m.len(),
                m_sorted
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }
    let warnings = collect_naming_warnings(
        target_bed.as_ref(),
        mutation_bed.as_ref(),
        &reference_chroms,
    );
    debug!(
        "chrom-naming-mismatch check produced {} warning(s)",
        warnings.len()
    );

    // Write per-record VCF artifacts before the JSON/TXT report so the
    // report's `outputs` section can name them. FN_with_reasons.vcf is
    // always written; FP.vcf is gated on the config flag.
    // Each artifact inherits declarations from the VCF ITS records came from:
    // FNs are golden-VCF records, FPs are called-VCF records (#444).
    let fn_vcf_path = write_fn_with_reasons(
        &attribution_result,
        &config.output_dir,
        config.overwrite_output,
        &config.golden_vcf,
    )?;
    info!("Wrote {}", fn_vcf_path.display());

    let fp_vcf_path = if config.write_fp_vcf {
        let mut fp_pairs: Vec<(String, &Variant)> = Vec::new();
        for (chrom, b) in &buckets {
            for v in &b.fps {
                fp_pairs.push((chrom.clone(), v));
            }
        }
        let p = write_fp_vcf(
            &fp_pairs,
            &config.output_dir,
            config.overwrite_output,
            &config.called_vcf,
        )?;
        info!("Wrote {}", p.display());
        Some(p)
    } else {
        None
    };

    let json_path = config.output_dir.join("comparison_summary.json");
    let txt_path = config.output_dir.join("comparison_summary.txt");

    let summary = ComparisonSummary {
        schema_version: SCHEMA_VERSION,
        inputs: SummaryInputs {
            golden_vcf: config.golden_vcf.clone(),
            called_vcf: config.called_vcf.clone(),
            reference: config.reference.clone(),
            target_bed: config.target_bed.clone(),
            contigs_simulated: config.contigs_simulated.clone(),
            include_homs: config.include_homs,
            include_filtered: config.include_filtered,
            equivalence_window: config.equivalence_window,
            fast: config.fast,
            mutation_bed: config.mutation_bed.clone(),
            chrom_aliases: config.chrom_aliases.clone(),
        },
        totals,
        metrics,
        per_contig,
        skipped_golden,
        skipped_called,
        fn_attribution: attribution_result.counts,
        warnings,
        outputs: SummaryOutputs {
            comparison_summary_json: json_path.clone(),
            comparison_summary_txt: txt_path.clone(),
            fn_with_reasons_vcf: fn_vcf_path.clone(),
            fp_vcf: fp_vcf_path.clone(),
        },
    };

    let written_json = write_json(&summary, &config.output_dir, config.overwrite_output)?;
    let written_txt = write_txt(&summary, &config.output_dir, config.overwrite_output)?;
    info!("Wrote {}", written_json.display());
    info!("Wrote {}", written_txt.display());
    Ok(())
}

/// Stream the reference FASTA and run the equivalence sweep on each contig
/// with at least one FP. Mutates `buckets` in place: consumed FNs are
/// removed from the bucket's `fns` list (and tracked separately as a count
/// of promoted equivalents); consumed FPs are removed from `fps`.
fn run_equivalence_sweep(
    buckets: &mut BTreeMap<String, ContigBuckets>,
    reference_path: &std::path::Path,
    window_bp: usize,
) -> Result<(), CompareVcfsError> {
    let need_reference = buckets.values().any(|b| !b.fps.is_empty());
    if !need_reference {
        return Ok(());
    }

    let stream = FastaStream::open(&reference_path.to_path_buf())
        .map_err(|e| CompareVcfsError::ConfigurationError(format!("reference: {e}")))?;
    for result in stream {
        let (contig, sequence) =
            result.map_err(|e| CompareVcfsError::ConfigurationError(format!("FASTA: {e}")))?;
        let Some(bucket) = buckets.get_mut(&contig) else {
            continue;
        };
        if bucket.fps.is_empty() {
            continue;
        }
        let ref_bytes = sequence.as_bytes();
        let sweep = equivalence::sweep(&bucket.fps, &bucket.fns, ref_bytes, window_bp);
        if sweep.equivalents_promoted == 0 {
            continue;
        }
        let promoted = sweep.equivalents_promoted;
        // Sort indices descending so each `swap_remove` doesn't disturb the next.
        let mut fp_idx = sweep.fp_consumed;
        fp_idx.sort_unstable_by(|a, b| b.cmp(a));
        for i in fp_idx {
            bucket.fps.swap_remove(i);
        }
        let mut fn_idx = sweep.fn_consumed;
        fn_idx.sort_unstable_by(|a, b| b.cmp(a));
        for i in fn_idx {
            bucket.fns.swap_remove(i);
        }
        bucket.tp_count += promoted;
        bucket.equivalents_promoted += promoted;
        debug!("{contig}: {promoted} equivalents promoted by sweep");
    }
    Ok(())
}

/// Drops disqualified variants and tallies the reasons. The returned map keeps
/// only contigs whose surviving variant list is non-empty so downstream
/// classification doesn't iterate empty buckets.
fn filter_vcf(
    raw: &HashMap<String, Vec<Variant>>,
    target_bed: Option<&HashMap<String, Vec<BedRecord>>>,
    simulated_set: Option<&HashSet<String>>,
    include_homs: bool,
    include_filtered: bool,
) -> (HashMap<String, Vec<Variant>>, SkipCounts) {
    let mut kept: HashMap<String, Vec<Variant>> = HashMap::new();
    let mut skipped = SkipCounts::default();
    for (chrom, variants) in raw {
        if let Some(sim) = simulated_set
            && !sim.contains(chrom)
        {
            skipped.outside_simulated_contigs += variants.len();
            continue;
        }
        for v in variants {
            // Symbolic / structural ALTs (`<DEL>`, `<DUP>`, `<CNV>`, ...) have
            // no comparable byte sequence — compare-vcfs is built on byte-wise
            // REF/ALT identity (and an equivalence sweep that splices literal
            // bases). Skip them before anything calls `.as_literal()`.
            if v.alternate.is_symbolic() {
                skipped.symbolic += 1;
                continue;
            }
            if !include_homs && v.reference.as_slice() == v.alternate.as_literal().unwrap() {
                skipped.homozygous_ref += 1;
                continue;
            }
            if !include_filtered && !filter_is_pass(v) {
                skipped.filtered += 1;
                continue;
            }
            if let Some(bed) = target_bed {
                let regions = bed.get(chrom);
                let in_target = regions
                    .map(|rs| rs.iter().any(|r| r.contains(chrom, v.location)))
                    .unwrap_or(false);
                if !in_target {
                    skipped.outside_target_bed += 1;
                    continue;
                }
            }
            kept.entry(chrom.clone()).or_default().push(v.clone());
        }
    }
    (kept, skipped)
}

fn filter_is_pass(v: &Variant) -> bool {
    match v.filter.as_deref() {
        None | Some("PASS") | Some(".") | Some("") => true,
        Some(_) => false,
    }
}

/// Exact-match classification. Splits each contig's variants into
/// TP (count only), FN (kept as `Variant` so equivalence sweep can apply
/// them), and FP (same).
fn classify(
    golden: &HashMap<String, Vec<Variant>>,
    called: &HashMap<String, Vec<Variant>>,
) -> BTreeMap<String, ContigBuckets> {
    let mut out: BTreeMap<String, ContigBuckets> = BTreeMap::new();
    let all_contigs: HashSet<&String> = golden.keys().chain(called.keys()).collect();
    for chrom in all_contigs {
        let g_keys: HashSet<VariantKey> = golden
            .get(chrom)
            .map(|vs| vs.iter().map(VariantKey::from).collect())
            .unwrap_or_default();
        let c_keys: HashSet<VariantKey> = called
            .get(chrom)
            .map(|vs| vs.iter().map(VariantKey::from).collect())
            .unwrap_or_default();
        let tp = g_keys.intersection(&c_keys).count();
        let fns: Vec<Variant> = golden
            .get(chrom)
            .into_iter()
            .flatten()
            .filter(|v| !c_keys.contains(&VariantKey::from(*v)))
            .cloned()
            .collect();
        let fps: Vec<Variant> = called
            .get(chrom)
            .into_iter()
            .flatten()
            .filter(|v| !g_keys.contains(&VariantKey::from(*v)))
            .cloned()
            .collect();
        out.insert(
            chrom.clone(),
            ContigBuckets {
                tp_count: tp,
                equivalents_promoted: 0,
                fns,
                fps,
            },
        );
    }
    out
}

/// Roll `ContigBuckets` (post-sweep) into the report's `ContigCounts` shape.
fn aggregate(
    buckets: &BTreeMap<String, ContigBuckets>,
) -> (BTreeMap<String, ContigCounts>, ContigCounts) {
    let mut per_contig: BTreeMap<String, ContigCounts> = BTreeMap::new();
    let mut totals = ContigCounts {
        tp: 0,
        fn_: 0,
        fp: 0,
        equivalents_promoted: 0,
    };
    for (chrom, b) in buckets {
        let cc = ContigCounts {
            tp: b.tp_count,
            fn_: b.fns.len(),
            fp: b.fps.len(),
            equivalents_promoted: b.equivalents_promoted,
        };
        totals.tp += cc.tp;
        totals.fn_ += cc.fn_;
        totals.fp += cc.fp;
        totals.equivalents_promoted += cc.equivalents_promoted;
        per_contig.insert(chrom.clone(), cc);
    }
    (per_contig, totals)
}

/// Compose the reference's "known contig set" from the two VCFs. Used only
/// to detect chrom-naming mismatches against BEDs — if every BED-named
/// contig is absent from both VCFs, the BED is almost certainly using a
/// different naming convention than the reference.
fn derive_reference_chroms(
    golden: &HashMap<String, Vec<Variant>>,
    called: &HashMap<String, Vec<Variant>>,
) -> HashSet<String> {
    let mut out: HashSet<String> = HashSet::new();
    out.extend(golden.keys().cloned());
    out.extend(called.keys().cloned());
    out
}

fn collect_naming_warnings(
    target_bed: Option<&HashMap<String, Vec<eidolon_core::structs::bed_record::BedRecord>>>,
    mutation_bed: Option<&HashMap<String, Vec<eidolon_core::structs::bed_record::BedRecord>>>,
    reference_chroms: &HashSet<String>,
) -> Vec<ChromNamingWarning> {
    let mut out = Vec::new();
    if let Some(t) = target_bed
        && let Some(w) = detect_chrom_naming_mismatch("target_bed", t, reference_chroms)
    {
        out.push(w);
    }
    if let Some(m) = mutation_bed
        && let Some(w) = detect_chrom_naming_mismatch("mutation_bed", m, reference_chroms)
    {
        out.push(w);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::structs::variants::{AlternateType, Genotype, Provenance, VariantType};

    fn snp(loc: usize, ref_b: Nucleotide, alt_b: Nucleotide, filter: Option<&str>) -> Variant {
        Variant {
            variant_type: VariantType::SNP,
            location: loc,
            reference: vec![ref_b],
            alternate: AlternateType::Literal(vec![alt_b]),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: filter.map(str::to_string),
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        }
    }

    #[test]
    fn classify_all_tp() {
        let v = snp(10, Nucleotide::A, Nucleotide::C, Some("PASS"));
        let mut g = HashMap::new();
        g.insert("chr1".into(), vec![v.clone()]);
        let mut c = HashMap::new();
        c.insert("chr1".into(), vec![v]);
        let buckets = classify(&g, &c);
        let (_per, tot) = aggregate(&buckets);
        assert_eq!(tot.tp, 1);
        assert_eq!(tot.fn_, 0);
        assert_eq!(tot.fp, 0);
        assert_eq!(buckets["chr1"].tp_count, 1);
    }

    #[test]
    fn classify_mixed() {
        // golden has (10,A,C) and (20,G,T); called has (10,A,C) and (30,T,A).
        // → 1 TP, 1 FN, 1 FP.
        let mut g = HashMap::new();
        g.insert(
            "chr1".into(),
            vec![
                snp(10, Nucleotide::A, Nucleotide::C, Some("PASS")),
                snp(20, Nucleotide::G, Nucleotide::T, Some("PASS")),
            ],
        );
        let mut c = HashMap::new();
        c.insert(
            "chr1".into(),
            vec![
                snp(10, Nucleotide::A, Nucleotide::C, Some("PASS")),
                snp(30, Nucleotide::T, Nucleotide::A, Some("PASS")),
            ],
        );
        let buckets = classify(&g, &c);
        let (_per, tot) = aggregate(&buckets);
        assert_eq!(tot.tp, 1);
        assert_eq!(tot.fn_, 1);
        assert_eq!(tot.fp, 1);
        assert_eq!(buckets["chr1"].fns.len(), 1);
        assert_eq!(buckets["chr1"].fps.len(), 1);
    }

    #[test]
    fn classify_contig_only_in_one_side() {
        // Golden has variants on chr1 only; called has variants on chr2 only.
        // Per-contig: chr1 → 1 FN, chr2 → 1 FP. Totals: 0 TP, 1 FN, 1 FP.
        let mut g = HashMap::new();
        g.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("PASS"))],
        );
        let mut c = HashMap::new();
        c.insert(
            "chr2".into(),
            vec![snp(50, Nucleotide::G, Nucleotide::A, Some("PASS"))],
        );
        let buckets = classify(&g, &c);
        let (per, tot) = aggregate(&buckets);
        assert_eq!(tot.tp, 0);
        assert_eq!(tot.fn_, 1);
        assert_eq!(tot.fp, 1);
        assert_eq!(per["chr1"].fn_, 1);
        assert_eq!(per["chr2"].fp, 1);
    }

    #[test]
    fn filter_skips_homref_when_not_included() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::A, Some("PASS"))],
        );
        let (kept, skipped) = filter_vcf(&raw, None, None, false, false);
        assert!(kept.is_empty());
        assert_eq!(skipped.homozygous_ref, 1);
    }

    #[test]
    fn filter_skips_non_pass_when_not_included() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("LowQual"))],
        );
        let (kept, skipped) = filter_vcf(&raw, None, None, false, false);
        assert!(kept.is_empty());
        assert_eq!(skipped.filtered, 1);
    }

    #[test]
    fn filter_includes_non_pass_when_flag_set() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("LowQual"))],
        );
        let (kept, skipped) = filter_vcf(&raw, None, None, false, true);
        assert_eq!(kept.get("chr1").unwrap().len(), 1);
        assert_eq!(skipped.filtered, 0);
    }

    #[test]
    fn filter_drops_unsimulated_contigs() {
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![snp(10, Nucleotide::A, Nucleotide::C, Some("PASS"))],
        );
        raw.insert(
            "chr2".into(),
            vec![snp(20, Nucleotide::G, Nucleotide::T, Some("PASS"))],
        );
        let sim: HashSet<String> = ["chr1".to_string()].into_iter().collect();
        let (kept, skipped) = filter_vcf(&raw, None, Some(&sim), false, false);
        assert!(kept.contains_key("chr1"));
        assert!(!kept.contains_key("chr2"));
        assert_eq!(skipped.outside_simulated_contigs, 1);
    }

    #[test]
    fn filter_skips_symbolic_alts_without_panicking() {
        use eidolon_core::structs::variants::{SvData, SvType};
        let sv = Variant {
            variant_type: VariantType::Complex,
            location: 100,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            genotype_str: "1/1".to_string(),
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: Some("PASS".to_string()),
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };
        let mut raw = HashMap::new();
        raw.insert(
            "chr1".into(),
            vec![sv, snp(50, Nucleotide::A, Nucleotide::C, Some("PASS"))],
        );
        // filter_vcf must route the symbolic record into `skipped.symbolic`
        // rather than calling `.as_literal().unwrap()` on it — that path
        // previously panicked when a round-tripped <DEL> reached the
        // byte-wise REF/ALT comparison.
        let (kept, skipped) = filter_vcf(&raw, None, None, false, false);
        assert_eq!(kept.get("chr1").map(|v| v.len()), Some(1));
        assert_eq!(kept["chr1"][0].location, 50);
        assert_eq!(skipped.symbolic, 1);
    }
}
