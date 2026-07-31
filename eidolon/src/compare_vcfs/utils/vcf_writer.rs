//! VCF output writers for `compare-vcfs`.
//!
//! Two artifacts:
//!   - `FN_with_reasons.vcf` — every surviving FN annotated with a
//!     `EIDOLON_REASON` INFO tag listing the attribution reasons.
//!   - `FP.vcf` — every surviving FP, as-is. Optional; gated on
//!     `write_fp_vcf` in the config because some pipelines accumulate
//!     large FP sets and the user may not want the extra artifact.
//!
//! Both writers emit minimal VCFv4.2 headers and reuse the original
//! Variant's per-column data (ID, QUAL, FILTER, INFO, FORMAT, SAMPLE).
//! We do **not** preserve the source VCF's `##contig=` lines because
//! that would require re-reading the input header; downstream tools that
//! need it can re-create from the reference FASTA.
use crate::compare_vcfs::errors::CompareVcfsError;
use crate::compare_vcfs::utils::attribution::{AttributionResult, Reason};
use eidolon_core::structs::{
    nucleotides::sequence_array_to_string,
    variants::{AlternateType, Variant},
};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

const HEADER_BASE: &str = "\
##fileformat=VCFv4.2
##source=eidolon compare-vcfs
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
";
const EIDOLON_REASON_INFO: &str = "##INFO=<ID=EIDOLON_REASON,Number=.,Type=String,\
Description=\"Comma-separated NEAT-aware false-negative attribution reasons\">\n";
const COLUMN_HEADER: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

/// Write FN records with EIDOLON_REASON annotations to
/// `<output_dir>/FN_with_reasons.vcf`. Returns the written path.
///
/// If `attribution.per_fn` is empty, the file is still written (header
/// only) so downstream tools can rely on the artifact existing.
pub fn write_fn_with_reasons(
    attribution: &AttributionResult,
    output_dir: &Path,
    overwrite_output: bool,
) -> Result<PathBuf, CompareVcfsError> {
    let path = output_dir.join("FN_with_reasons.vcf");
    check_overwrite(&path, overwrite_output)?;
    let mut file = fs::File::create(&path)?;
    file.write_all(HEADER_BASE.as_bytes())?;
    file.write_all(EIDOLON_REASON_INFO.as_bytes())?;
    file.write_all(COLUMN_HEADER.as_bytes())?;
    for (chrom, variant, reasons) in &attribution.per_fn {
        let line = format_record(chrom, variant, Some(reasons));
        file.write_all(line.as_bytes())?;
        file.write_all(b"\n")?;
    }
    Ok(path)
}

/// Write FP records to `<output_dir>/FP.vcf`. Returns the written path.
pub fn write_fp_vcf(
    fps_by_contig: &[(String, &Variant)],
    output_dir: &Path,
    overwrite_output: bool,
) -> Result<PathBuf, CompareVcfsError> {
    let path = output_dir.join("FP.vcf");
    check_overwrite(&path, overwrite_output)?;
    let mut file = fs::File::create(&path)?;
    file.write_all(HEADER_BASE.as_bytes())?;
    file.write_all(COLUMN_HEADER.as_bytes())?;
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

    #[test]
    fn format_record_synthesizes_gt_column_when_format_empty() {
        let v = snp_with(100, None);
        let line = format_record("chr1", &v, None);
        // Trailing two tab-separated fields are FORMAT and SAMPLE.
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts[8], "GT");
        assert_eq!(parts[9], "0/1");
    }

    #[test]
    fn fn_with_reasons_writes_header_and_records() {
        let dir = tempfile::tempdir().unwrap();
        let v = snp_with(100, None);
        let attribution = AttributionResult {
            per_fn: vec![("chr1".to_string(), v, vec![Reason::Unknown])],
            counts: Default::default(),
        };
        let path = write_fn_with_reasons(&attribution, dir.path(), false).unwrap();
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
        let path = write_fn_with_reasons(&attribution, dir.path(), false).unwrap();
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
        let path = write_fp_vcf(&fps, dir.path(), false).unwrap();
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
        let err = write_fn_with_reasons(&attribution, dir.path(), false).unwrap_err();
        assert!(matches!(err, CompareVcfsError::OverwriteFileError(_)));
    }
}
