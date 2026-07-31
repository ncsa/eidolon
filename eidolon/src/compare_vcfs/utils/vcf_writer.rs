//! VCF output writers for `compare-vcfs`.
//!
//! Two artifacts:
//!   - `FN_with_reasons.vcf` — every surviving FN annotated with a
//!     `EIDOLON_REASON` INFO tag listing the attribution reasons.
//!   - `FP.vcf` — every surviving FP, as-is. Optional; gated on
//!     `write_fp_vcf` in the config because some pipelines accumulate
//!     large FP sets and the user may not want the extra artifact.
//!
//! Both writers reuse the original Variant's per-column data (ID, QUAL, FILTER,
//! INFO, FORMAT, SAMPLE), and therefore **inherit the declaration lines from the
//! VCF those records came from** (#444).
//!
//! An earlier version emitted a fixed minimal header, on the reasoning that
//! re-reading the input header was avoidable and "downstream tools that need
//! `##contig` can re-create it from the reference FASTA". That does not hold: a
//! record's INFO column is passed through verbatim, so the artifact carried tags
//! its own header never declared, and no downstream tool can invent a declaration
//! for a tag it has never seen. bcftools tolerates that streaming VCF→VCF but
//! **hard-fails on BCF translation** — which `bcftools concat | sort` does
//! internally — so it surfaced far downstream as an opaque failure.
//!
//! A fixed declaration set cannot work either, because the two artifacts have
//! different provenance: FN records come from the **golden** VCF, while FP records
//! come from the **called** VCF, whose INFO is whatever the caller emitted
//! (Mutect2's `TLOD`/`MBQ`, Strelka's `SomaticEVS`, …). Each artifact therefore
//! inherits from its own source.
use crate::compare_vcfs::errors::CompareVcfsError;
use crate::compare_vcfs::utils::attribution::{AttributionResult, Reason};
use eidolon_core::file_tools::file_io::{is_gzipped_file, read_gzip_lines, read_lines};
use eidolon_core::structs::{
    nucleotides::sequence_array_to_string,
    variants::{AlternateType, Variant},
};
use std::collections::HashSet;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Declaration lines carried over from the VCF a set of records came from.
struct SourceHeader {
    /// `##INFO` / `##FILTER` / `##FORMAT` / `##ALT` / `##contig` lines, in source order.
    meta: Vec<String>,
    /// Declared `##contig` IDs, so undeclared record contigs can be backfilled.
    contig_ids: HashSet<String>,
    /// True if the source declared a `GT` FORMAT key. When it didn't, we add one,
    /// because `format_record` synthesizes a GT column for records with no FORMAT.
    has_gt: bool,
}

/// Read the `##`-meta lines up to (not including) `#CHROM`.
fn read_meta_lines(vcf: &PathBuf) -> Result<Vec<String>, CompareVcfsError> {
    let mut out = Vec::new();
    // The two readers have different concrete types; box them rather than making
    // this generic (matches the idiom in gen_seq_error_model's runner).
    let lines: Box<dyn Iterator<Item = std::io::Result<String>>> = if is_gzipped_file(vcf)? {
        Box::new(read_gzip_lines(vcf)?)
    } else {
        Box::new(read_lines(vcf)?)
    };
    for line in lines {
        let line = line?;
        if line.starts_with("#CHROM") {
            break;
        }
        out.push(line);
    }
    Ok(out)
}

/// Extract the `ID=` value from a structured header line (`##INFO=<ID=FOO,...>`).
fn header_line_id(line: &str) -> Option<&str> {
    let open = line.find('<')?;
    line[open + 1..]
        .trim_end_matches('>')
        .split(',')
        .find_map(|f| f.strip_prefix("ID="))
}

/// Collect the declarations an artifact must carry, from the VCF its records came from.
fn source_header(vcf: &PathBuf) -> Result<SourceHeader, CompareVcfsError> {
    const KEEP: [&str; 5] = ["##INFO=", "##FILTER=", "##FORMAT=", "##ALT=", "##contig="];
    let mut meta = Vec::new();
    let mut contig_ids = HashSet::new();
    let mut has_gt = false;
    for line in read_meta_lines(vcf)? {
        if !KEEP.iter().any(|k| line.starts_with(k)) {
            continue;
        }
        if let Some(id) = header_line_id(&line) {
            if line.starts_with("##contig=") {
                contig_ids.insert(id.to_string());
            } else if line.starts_with("##FORMAT=") && id == "GT" {
                has_gt = true;
            }
        }
        meta.push(line);
    }
    Ok(SourceHeader {
        meta,
        contig_ids,
        has_gt,
    })
}

/// Write the header: fileformat, inherited declarations, a backfilled `##contig`
/// for any record contig the source didn't declare, then the column line.
///
/// The backfill emits an ID-only contig (length is optional in VCF 4.2). It exists
/// so an artifact is self-describing even when the source header was itself
/// incomplete — a truth VCF from a third-party tool need not declare contigs.
fn write_header<W: Write>(
    out: &mut W,
    src: &SourceHeader,
    record_contigs: &[&str],
    extra_info: Option<&str>,
) -> Result<(), CompareVcfsError> {
    out.write_all(b"##fileformat=VCFv4.2\n")?;
    out.write_all(b"##source=eidolon compare-vcfs\n")?;
    for line in &src.meta {
        out.write_all(line.as_bytes())?;
        out.write_all(b"\n")?;
    }
    if !src.has_gt {
        out.write_all(GT_FORMAT.as_bytes())?;
    }
    if let Some(info) = extra_info {
        out.write_all(info.as_bytes())?;
    }
    let mut seen: HashSet<&str> = HashSet::new();
    for c in record_contigs {
        if src.contig_ids.contains(*c) || !seen.insert(*c) {
            continue;
        }
        writeln!(out, "##contig=<ID={c}>")?;
    }
    out.write_all(COLUMN_HEADER.as_bytes())?;
    Ok(())
}

const GT_FORMAT: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
const EIDOLON_REASON_INFO: &str = "##INFO=<ID=EIDOLON_REASON,Number=.,Type=String,\
Description=\"Comma-separated NEAT-aware false-negative attribution reasons\">\n";
const COLUMN_HEADER: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

/// Write FN records with EIDOLON_REASON annotations to
/// `<output_dir>/FN_with_reasons.vcf`. Returns the written path.
///
/// `golden_vcf` is the VCF the FN records came from; its declaration lines are
/// inherited so every INFO/FILTER/FORMAT tag the records carry is declared.
///
/// If `attribution.per_fn` is empty, the file is still written (header
/// only) so downstream tools can rely on the artifact existing.
pub fn write_fn_with_reasons(
    attribution: &AttributionResult,
    output_dir: &Path,
    overwrite_output: bool,
    golden_vcf: &PathBuf,
) -> Result<PathBuf, CompareVcfsError> {
    let path = output_dir.join("FN_with_reasons.vcf");
    check_overwrite(&path, overwrite_output)?;
    let src = source_header(golden_vcf)?;
    let contigs: Vec<&str> = attribution
        .per_fn
        .iter()
        .map(|(chrom, _, _)| chrom.as_str())
        .collect();
    let mut file = fs::File::create(&path)?;
    write_header(&mut file, &src, &contigs, Some(EIDOLON_REASON_INFO))?;
    for (chrom, variant, reasons) in &attribution.per_fn {
        let line = format_record(chrom, variant, Some(reasons));
        file.write_all(line.as_bytes())?;
        file.write_all(b"\n")?;
    }
    Ok(path)
}

/// Write FP records to `<output_dir>/FP.vcf`. Returns the written path.
///
/// `called_vcf` is the VCF the FP records came from — *not* the golden VCF. An FP's
/// INFO is whatever the caller emitted, so its declarations must come from there.
pub fn write_fp_vcf(
    fps_by_contig: &[(String, &Variant)],
    output_dir: &Path,
    overwrite_output: bool,
    called_vcf: &PathBuf,
) -> Result<PathBuf, CompareVcfsError> {
    let path = output_dir.join("FP.vcf");
    check_overwrite(&path, overwrite_output)?;
    let src = source_header(called_vcf)?;
    let contigs: Vec<&str> = fps_by_contig.iter().map(|(c, _)| c.as_str()).collect();
    let mut file = fs::File::create(&path)?;
    write_header(&mut file, &src, &contigs, None)?;
    for (chrom, v) in fps_by_contig {
        let line = format_record(chrom, v, None);
        file.write_all(line.as_bytes())?;
        file.write_all(b"\n")?;
    }
    Ok(path)
}

fn check_overwrite(path: &Path, allow: bool) -> Result<(), CompareVcfsError> {
    if path.is_file() && !allow {
        return Err(CompareVcfsError::OverwriteFileError(
            path.display().to_string(),
        ));
    }
    Ok(())
}

/// Serialize one Variant as a tab-separated VCF record. If `reasons` is
/// `Some`, append `EIDOLON_REASON=<csv>` to the INFO column.
fn format_record(chrom: &str, v: &Variant, reasons: Option<&[Reason]>) -> String {
    let id = v.id.as_deref().unwrap_or(".");
    let ref_str = sequence_array_to_string(&v.reference);
    let alt_str = match &v.alternate {
        AlternateType::Literal(bases) => sequence_array_to_string(bases),
        AlternateType::Symbolic(sv) => sv.raw_alt.clone(),
    };
    let qual = match v.quality_score {
        Some(q) => q.to_string(),
        None => ".".to_string(),
    };
    let filter = v.filter.as_deref().unwrap_or(".");
    let info_base = v.info.as_deref().unwrap_or(".");
    let info = match reasons {
        Some(rs) if !rs.is_empty() => {
            let reason_csv = rs.iter().map(|r| r.as_str()).collect::<Vec<_>>().join(",");
            // Preserve an existing INFO blob if present; replace bare "." since
            // it's the "no info" sentinel.
            if info_base == "." || info_base.is_empty() {
                format!("EIDOLON_REASON={reason_csv}")
            } else {
                format!("{info_base};EIDOLON_REASON={reason_csv}")
            }
        }
        _ => info_base.to_string(),
    };
    let (format_col, sample_col) = if v.format.is_empty() {
        // Synthesize a minimal GT column from genotype_str so the line is
        // well-formed (we declared GT in the header).
        ("GT".to_string(), v.genotype_str.clone())
    } else {
        (v.format.join(":"), v.sample.join(":"))
    };

    format!(
        "{chrom}\t{pos}\t{id}\t{ref_str}\t{alt_str}\t{qual}\t{filter}\t{info}\t{format_col}\t{sample_col}",
        pos = v.location,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::structs::variants::{AlternateType, Provenance};
    use eidolon_core::structs::{
        nucleotides::Nucleotide,
        variants::{Genotype, VariantType},
    };

    fn snp_with(loc: usize, info: Option<&str>) -> Variant {
        Variant {
            variant_type: VariantType::SNP,
            location: loc,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Literal(vec![Nucleotide::C]),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: Some(60),
            filter: Some("PASS".to_string()),
            info: info.map(str::to_string),
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        }
    }

    #[test]
    fn format_record_adds_neat_reason_when_info_is_dot() {
        let v = snp_with(100, Some("."));
        let line = format_record("chr1", &v, Some(&[Reason::Unknown]));
        assert!(line.contains("\tEIDOLON_REASON=unknown\t"));
        // Original "." should be replaced, not preserved as ".;EIDOLON_REASON=...".
        assert!(!line.contains(".;EIDOLON_REASON"));
    }

    #[test]
    fn format_record_preserves_existing_info_and_appends_reason() {
        let v = snp_with(100, Some("DP=30;AF=0.5"));
        let line = format_record(
            "chr1",
            &v,
            Some(&[Reason::OutsideMutationBed, Reason::OutsideTargetBed]),
        );
        assert!(
            line.contains("DP=30;AF=0.5;EIDOLON_REASON=outside_mutation_bed,outside_target_bed"),
            "got: {line}"
        );
    }

    #[test]
    fn format_record_without_reasons_keeps_info_unchanged() {
        let v = snp_with(100, Some("DP=30"));
        let line = format_record("chr1", &v, None);
        assert!(line.contains("\tDP=30\tGT\t0/1"));
        assert!(!line.contains("EIDOLON_REASON"));
    }

    /// A source VCF whose header declares tags the artifact must inherit. `DP` stands
    /// in for a caller-specific tag we could never have enumerated ahead of time.
    fn source_vcf(dir: &std::path::Path, name: &str) -> PathBuf {
        let p = dir.join(name);
        std::fs::write(
            &p,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=1000>\n\
             ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
             ##FILTER=<ID=LowQual,Description=\"low\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .unwrap();
        p
    }

    #[test]
    fn format_record_synthesizes_gt_column_when_format_empty() {
        let v = snp_with(100, None);
        let line = format_record("chr1", &v, None);
        // Trailing two tab-separated fields are FORMAT and SAMPLE.
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts[8], "GT");
        assert_eq!(parts[9], "0/1");
    }

    /// #444: the artifact must declare every INFO tag its records carry. A record's
    /// INFO is passed through verbatim from the source VCF, so a fixed header left
    /// undeclared tags behind — fine streaming VCF->VCF, a hard failure on BCF
    /// translation (which `bcftools concat | sort` does internally).
    #[test]
    fn artifact_declares_every_info_tag_its_records_carry() {
        let dir = tempfile::tempdir().unwrap();
        let mut v = snp_with(100, None);
        // A caller-specific tag we could not have enumerated in a fixed header.
        v.info = Some("DP=30".to_string());
        let attribution = AttributionResult {
            per_fn: vec![("chr1".to_string(), v, vec![Reason::Unknown])],
            counts: Default::default(),
        };
        let src = source_vcf(dir.path(), "golden.vcf");
        let path = write_fn_with_reasons(&attribution, dir.path(), false, &src).unwrap();
        let body = std::fs::read_to_string(&path).unwrap();

        let declared: HashSet<&str> = body
            .lines()
            .filter(|l| l.starts_with("##INFO="))
            .filter_map(header_line_id)
            .collect();
        let record = body
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("no record written");
        let info = record.split('\t').nth(7).unwrap();
        for field in info.split(';') {
            let key = field.split('=').next().unwrap();
            assert!(
                declared.contains(key),
                "INFO key `{key}` used but not declared (declared: {declared:?}) in:\n{body}"
            );
        }
        // Specifically: inherited from the source, not invented by us.
        assert!(body.contains("##INFO=<ID=DP,"), "DP not inherited:\n{body}");
        assert!(body.contains("##INFO=<ID=EIDOLON_REASON,"));
    }

    /// FILTER and contig declarations are inherited too, and a record contig the
    /// source failed to declare is backfilled so the artifact is self-describing.
    #[test]
    fn artifact_inherits_filter_and_backfills_undeclared_contigs() {
        let dir = tempfile::tempdir().unwrap();
        let mut v = snp_with(100, None);
        v.filter = Some("LowQual".to_string());
        let other = snp_with(200, None);
        let attribution = AttributionResult {
            per_fn: vec![
                ("chr1".to_string(), v, vec![Reason::Unknown]),
                // chrUn is absent from the source header -> must be backfilled.
                ("chrUn".to_string(), other, vec![Reason::Unknown]),
            ],
            counts: Default::default(),
        };
        let src = source_vcf(dir.path(), "golden.vcf");
        let path = write_fn_with_reasons(&attribution, dir.path(), false, &src).unwrap();
        let body = std::fs::read_to_string(&path).unwrap();
        assert!(
            body.contains("##FILTER=<ID=LowQual,"),
            "FILTER not inherited"
        );
        assert!(
            body.contains("##contig=<ID=chr1,length=1000>"),
            "contig not inherited"
        );
        assert!(
            body.contains("##contig=<ID=chrUn>"),
            "undeclared contig not backfilled"
        );
        // GT was declared by the source, so we must not emit a duplicate.
        assert_eq!(
            body.matches("##FORMAT=<ID=GT,").count(),
            1,
            "duplicate GT decl"
        );
    }

    /// FP records come from the CALLED vcf, so they must inherit its declarations —
    /// using the golden VCF's header here would leave caller tags undeclared.
    #[test]
    fn fp_artifact_inherits_from_the_called_vcf() {
        let dir = tempfile::tempdir().unwrap();
        let called = dir.path().join("called.vcf");
        std::fs::write(
            &called,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=1000>\n\
             ##INFO=<ID=TLOD,Number=1,Type=Float,Description=\"mutect\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .unwrap();
        let mut v = snp_with(100, None);
        v.info = Some("TLOD=12.5".to_string());
        let fps = vec![("chr1".to_string(), &v)];
        let path = write_fp_vcf(&fps, dir.path(), false, &called).unwrap();
        let body = std::fs::read_to_string(&path).unwrap();
        assert!(
            body.contains("##INFO=<ID=TLOD,"),
            "caller tag not inherited:\n{body}"
        );
    }

    #[test]
    fn fn_with_reasons_writes_header_and_records() {
        let dir = tempfile::tempdir().unwrap();
        let v = snp_with(100, None);
        let attribution = AttributionResult {
            per_fn: vec![("chr1".to_string(), v, vec![Reason::Unknown])],
            counts: Default::default(),
        };
        let src = source_vcf(dir.path(), "golden.vcf");
        let path = write_fn_with_reasons(&attribution, dir.path(), false, &src).unwrap();
        let body = std::fs::read_to_string(&path).unwrap();
        assert!(body.contains("##fileformat=VCFv4.2"));
        assert!(body.contains("##INFO=<ID=EIDOLON_REASON"));
        assert!(body.contains("#CHROM\tPOS"));
        assert!(body.contains("EIDOLON_REASON=unknown"));
    }

    #[test]
    fn fn_with_reasons_writes_header_only_when_no_fns() {
        let dir = tempfile::tempdir().unwrap();
        let attribution = AttributionResult {
            per_fn: Vec::new(),
            counts: Default::default(),
        };
        let src = source_vcf(dir.path(), "golden.vcf");
        let path = write_fn_with_reasons(&attribution, dir.path(), false, &src).unwrap();
        let body = std::fs::read_to_string(&path).unwrap();
        assert!(body.contains("##fileformat=VCFv4.2"));
        // Last line should be the column header — no data rows past it.
        let after_header: Vec<&str> = body
            .lines()
            .skip_while(|l| l.starts_with('#') && !l.starts_with("#CHROM"))
            .skip(1)
            .collect();
        assert!(
            after_header.iter().all(|l| l.is_empty()),
            "expected no data rows, found: {after_header:?}"
        );
    }

    #[test]
    fn fp_vcf_writes_records_without_neat_reason_info() {
        let dir = tempfile::tempdir().unwrap();
        let v = snp_with(200, None);
        let fps: Vec<(String, &Variant)> = vec![("chr2".to_string(), &v)];
        let src = source_vcf(dir.path(), "called.vcf");
        let path = write_fp_vcf(&fps, dir.path(), false, &src).unwrap();
        let body = std::fs::read_to_string(&path).unwrap();
        assert!(!body.contains("EIDOLON_REASON"));
        assert!(body.contains("chr2\t200"));
    }

    #[test]
    fn overwrite_refused_when_existing() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("FN_with_reasons.vcf"), "stale").unwrap();
        let attribution = AttributionResult {
            per_fn: Vec::new(),
            counts: Default::default(),
        };
        let src = source_vcf(dir.path(), "golden.vcf");
        let err = write_fn_with_reasons(&attribution, dir.path(), false, &src).unwrap_err();
        assert!(matches!(err, CompareVcfsError::OverwriteFileError(_)));
    }
}
