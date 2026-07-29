//! Integration tests for `eidolon compare-vcfs`: FN_with_reasons.vcf and optional FP.vcf writers.

mod common;

use common::eidolon;
use std::io::Write;
use std::path::Path;

fn write_fasta(path: &Path, contig: &str, sequence: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">{contig}").unwrap();
    writeln!(f, "{sequence}").unwrap();
}

fn write_vcf(path: &Path, body_lines: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        f,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .unwrap();
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .unwrap();
    for line in body_lines {
        writeln!(f, "{line}").unwrap();
    }
}

fn write_text(path: &Path, content: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    write!(f, "{content}").unwrap();
}

fn load_summary(dir: &Path) -> serde_json::Value {
    let raw = std::fs::read_to_string(dir.join("comparison_summary.json")).unwrap();
    serde_json::from_str(&raw).unwrap()
}

/// FN_with_reasons.vcf is written unconditionally and contains the surviving
/// FN records, each annotated with a EIDOLON_REASON INFO tag.
#[test]
fn fn_with_reasons_vcf_is_written_with_neat_reason_tag() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    // Golden has two SNPs, called VCF is empty → 2 FNs.
    write_vcf(
        &golden,
        &[
            "chr1\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1",
            "chr1\t500\t.\tA\tG\t60\tPASS\t.\tGT\t1/1",
        ],
    );
    write_vcf(&called, &[]);

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
        ),
    );

    eidolon()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let fn_vcf = out.join("FN_with_reasons.vcf");
    assert!(fn_vcf.is_file(), "FN_with_reasons.vcf was not written");
    let body = std::fs::read_to_string(&fn_vcf).unwrap();
    assert!(body.contains("##INFO=<ID=EIDOLON_REASON"));
    assert!(body.contains("#CHROM\tPOS\tID\tREF\tALT"));
    // Both FNs should have EIDOLON_REASON=unknown (no mutation/target BED configured).
    let data_lines: Vec<&str> = body
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(
        data_lines.len(),
        2,
        "expected 2 FN rows, got: {data_lines:?}"
    );
    assert!(
        data_lines
            .iter()
            .all(|l| l.contains("EIDOLON_REASON=unknown"))
    );

    // outputs section in JSON should point at this file.
    let summary = load_summary(&out);
    assert!(
        summary["outputs"]["fn_with_reasons_vcf"]
            .as_str()
            .unwrap()
            .ends_with("FN_with_reasons.vcf")
    );
    // FP.vcf should NOT be present unless explicitly opted into.
    assert!(summary["outputs"].get("fp_vcf").is_none());
    assert!(!out.join("FP.vcf").exists());
}

/// `write_fp_vcf: true` produces FP.vcf with the called-only records, and
/// the JSON report's outputs section references it.
#[test]
fn fp_vcf_written_when_opted_in() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    write_vcf(&golden, &[]);
    write_vcf(
        &called,
        &[
            "chr1\t200\t.\tA\tT\t60\tPASS\tDP=42\tGT\t0/1",
            "chr1\t800\t.\tA\tG\t60\tPASS\t.\tGT\t1/1",
        ],
    );

    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\nwrite_fp_vcf: true\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
        ),
    );

    eidolon()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let fp_vcf = out.join("FP.vcf");
    assert!(fp_vcf.is_file(), "FP.vcf was not written");
    let body = std::fs::read_to_string(&fp_vcf).unwrap();
    // Two FP records expected.
    let data_lines: Vec<&str> = body
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(data_lines.len(), 2);
    // INFO column is preserved (one record has DP=42).
    assert!(data_lines.iter().any(|l| l.contains("DP=42")));
    // No EIDOLON_REASON tag on FP records.
    assert!(!body.contains("EIDOLON_REASON"));

    let summary = load_summary(&out);
    assert!(
        summary["outputs"]["fp_vcf"]
            .as_str()
            .unwrap()
            .ends_with("FP.vcf")
    );
}

/// FN_with_reasons.vcf is well-formed enough to round-trip back through
/// the existing read_vcf path — i.e., another eidolon pipeline could consume
/// it as input.
#[test]
fn fn_with_reasons_vcf_round_trips_via_existing_reader() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    write_vcf(&golden, &["chr1\t100\t.\tA\tC\t60\tPASS\tDP=20\tGT\t0/1"]);
    write_vcf(&called, &[]);
    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\noverwrite_output: true\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
        ),
    );
    eidolon()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    // Read back via the same read_vcf path the rest of eidolon uses.
    let fn_vcf = out.join("FN_with_reasons.vcf");
    let parsed = ::eidolon_core::file_tools::vcf_tools::read_vcf(fn_vcf).unwrap();
    let chr1_variants = parsed.get("chr1").expect("chr1 variants present");
    assert_eq!(chr1_variants.len(), 1);
    let v = &chr1_variants[0];
    assert_eq!(v.location, 100);
    // The reader stores the INFO blob; check our annotation made it through.
    let info = v.info.as_deref().unwrap_or("");
    assert!(
        info.contains("DP=20") && info.contains("EIDOLON_REASON=unknown"),
        "info field was: {info:?}"
    );
}

/// TXT report should list every output path in the Outputs section.
#[test]
fn txt_report_lists_every_output_path() {
    let tmp = tempfile::tempdir().unwrap();
    let fa = tmp.path().join("ref.fa");
    let golden = tmp.path().join("golden.vcf");
    let called = tmp.path().join("called.vcf");
    let out = tmp.path().join("report");
    write_fasta(&fa, "chr1", &"A".repeat(1000));
    write_vcf(&golden, &["chr1\t100\t.\tA\tC\t60\tPASS\t.\tGT\t0/1"]);
    write_vcf(&called, &["chr1\t300\t.\tA\tT\t60\tPASS\t.\tGT\t0/1"]);
    let yaml = tmp.path().join("cfg.yml");
    write_text(
        &yaml,
        &format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {}\noutput_dir: {}\n\
             overwrite_output: true\nwrite_fp_vcf: true\n",
            golden.display(),
            called.display(),
            fa.display(),
            out.display(),
        ),
    );
    eidolon()
        .args(["compare-vcfs", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let txt = std::fs::read_to_string(out.join("comparison_summary.txt")).unwrap();
    assert!(txt.contains("Outputs"));
    assert!(txt.contains("comparison_summary.json"));
    assert!(txt.contains("comparison_summary.txt"));
    assert!(txt.contains("FN_with_reasons.vcf"));
    assert!(txt.contains("FP.vcf"));
}
