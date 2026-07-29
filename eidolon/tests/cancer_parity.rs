//! Parity test: native `eidolon gen-cancer-reads` vs the legacy
//! `tools/cancer_simulate.sh`.
//!
//! The native subcommand (#239) is a port of the shell orchestrator. Before the
//! script can be retired we need proof the two produce equivalent output. This
//! test runs both over the same H1N1 reference, with identical purity / coverage
//! / seed, and asserts:
//!
//!   1. The merged FASTQs carry the same record **multiset** (both tag `N_`/`T_`
//!      and concatenate normal→tumor; gen-reads output *order* is not stable —
//!      see `determinism.rs` — so we compare sorted records, not lines).
//!   2. The per-pass golden VCF bodies match (same gen-reads driven identically).
//!   3. The origin-tagged truth VCF carries the same
//!      `(chrom, pos, ref, alt, EIDOLON_ORIGIN)` set — i.e. the native Rust merge
//!      classifies germline / somatic / shared exactly as the bcftools+awk merge.
//!
//! The bash path needs `bash` (+ `awk`/`gzip`, universally present on the dev/CI
//! Linux image). The truth-VCF merge additionally needs `bcftools`/`bgzip`; if
//! those are absent the script skips that step, and so do we (with a note).
//! When `bash` itself is missing the whole test skips rather than fails, so a
//! minimal environment doesn't break the suite.

mod common;

use common::{eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

/// `true` if `name` resolves on PATH.
fn tool_available(name: &str) -> bool {
    Command::new("sh")
        .arg("-c")
        .arg(format!("command -v {name}"))
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Group FASTQ lines into 4-line records and sort — turns output order (which is
/// not deterministic across runs) into a stable multiset comparison.
fn sorted_records(lines: Vec<String>) -> Vec<[String; 4]> {
    assert!(
        lines.len().is_multiple_of(4),
        "FASTQ line count must be a multiple of 4"
    );
    let mut records: Vec<[String; 4]> = lines
        .chunks(4)
        .map(|c| [c[0].clone(), c[1].clone(), c[2].clone(), c[3].clone()])
        .collect();
    records.sort();
    records
}

/// Sorted multiset of a gzipped VCF's body (non-`#`) lines.
fn vcf_body_sorted(path: &Path) -> Vec<String> {
    let mut body: Vec<String> = read_gzip_fastq_lines(path)
        .into_iter()
        .filter(|l| !l.starts_with('#'))
        .collect();
    body.sort();
    body
}

/// Set of `(chrom, pos, ref, alt, origin)` for every body record. `origin` is the
/// `EIDOLON_ORIGIN=` value parsed out of the INFO column (col 8).
fn vcf_origin_set(path: &Path) -> BTreeSet<(String, String, String, String, String)> {
    read_gzip_fastq_lines(path)
        .into_iter()
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            if c.len() < 8 {
                return None;
            }
            let origin = c[7]
                .split(';')
                .find_map(|f| f.strip_prefix("EIDOLON_ORIGIN="))
                .unwrap_or("?")
                .to_string();
            Some((c[0].into(), c[1].into(), c[3].into(), c[4].into(), origin))
        })
        .collect()
}

#[test]
fn native_gen_cancer_reads_matches_cancer_simulate_sh() {
    if !tool_available("bash") {
        eprintln!("skipping cancer parity test: `bash` not available");
        return;
    }

    let script = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("tools")
        .join("cancer_simulate.sh");
    assert!(
        script.is_file(),
        "cancer_simulate.sh not found at {script:?}"
    );

    // The same eidolon binary drives both paths (bash shells out to it via --eidolon-bin).
    let bin = assert_cmd::cargo::cargo_bin("eidolon");

    let (_d, work) = fresh_workdir();
    let native_dir = work.join("native");
    let bash_dir = work.join("bash");
    std::fs::create_dir_all(&native_dir).unwrap();
    std::fs::create_dir_all(&bash_dir).unwrap();

    // Shared parameters. purity 0.6 → normal 12x / tumor 18x (asymmetric, so a
    // swapped-pass bug in either path would surface). Defaults match across the
    // two: tumor_mutation_rate 1e-5, sv_rate_scale 0.0, model-fitted normal rate.
    const TOTAL_COV: &str = "30";
    const PURITY: &str = "0.6";
    const READ_LEN: &str = "70";
    const FRAG_MEAN: &str = "250";
    const FRAG_SD: &str = "30";
    const SEED: &str = "parity-seed";

    // ── native ──────────────────────────────────────────────────────────────
    let yaml = native_dir.join("cancer.yml");
    std::fs::write(
        &yaml,
        format!(
            "reference: {ref}\n\
             output_dir: {out}\n\
             output_prefix: sample\n\
             total_coverage: {TOTAL_COV}\n\
             purity: {PURITY}\n\
             read_len: {READ_LEN}\n\
             paired_ended: true\n\
             fragment_mean: {FRAG_MEAN}\n\
             fragment_st_dev: {FRAG_SD}\n\
             rng_seed: {SEED}\n\
             overwrite_output: true\n",
            ref = h1n1_reference().display(),
            out = native_dir.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-cancer-reads", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    // ── bash ────────────────────────────────────────────────────────────────
    let status = Command::new("bash")
        .arg(&script)
        .arg("--reference")
        .arg(h1n1_reference())
        .arg("--output-prefix")
        .arg(bash_dir.join("sample"))
        .args([
            "--total-coverage",
            TOTAL_COV,
            "--purity",
            PURITY,
            "--read-len",
            READ_LEN,
            "--paired-ended",
            "--fragment-mean",
            FRAG_MEAN,
            "--fragment-st-dev",
            FRAG_SD,
            "--rng-seed",
            SEED,
            "--eidolon-bin",
        ])
        .arg(&bin)
        .status()
        .expect("failed to spawn cancer_simulate.sh");
    assert!(status.success(), "cancer_simulate.sh exited non-zero");

    // ── 1. merged FASTQs: same record multiset ───────────────────────────────
    for r in ["r1", "r2"] {
        let nat = sorted_records(read_gzip_fastq_lines(
            &native_dir.join(format!("sample_merged_{r}.fastq.gz")),
        ));
        let bash = sorted_records(read_gzip_fastq_lines(
            &bash_dir.join(format!("sample_merged_{r}.fastq.gz")),
        ));
        assert!(!nat.is_empty(), "native produced no {r} reads");
        assert_eq!(
            nat, bash,
            "merged {r} FASTQ record multiset differs between native and bash"
        );
    }

    // ── 2. per-pass golden VCF bodies match ──────────────────────────────────
    for pass in ["normal", "tumor"] {
        let nat = vcf_body_sorted(&native_dir.join(format!("sample_{pass}.vcf.gz")));
        let bash = vcf_body_sorted(&bash_dir.join(format!("sample_{pass}.vcf.gz")));
        assert!(!nat.is_empty(), "native {pass} golden VCF has no records");
        assert_eq!(nat, bash, "{pass} golden VCF body differs native vs bash");
    }

    // ── 3. origin-tagged truth VCF: same classification set ──────────────────
    // The bash truth merge is gated on bcftools+bgzip; skip the comparison if the
    // script couldn't produce it (it warns and leaves the per-pass VCFs only).
    let bash_truth = bash_dir.join("sample_merged_truth.vcf.gz");
    let native_truth = native_dir.join("sample_merged_truth.vcf.gz");
    if bash_truth.is_file() {
        let nat = vcf_origin_set(&native_truth);
        let bash = vcf_origin_set(&bash_truth);
        assert!(!nat.is_empty(), "native truth VCF has no records");
        assert_eq!(
            nat, bash,
            "origin-tagged truth VCF (chrom,pos,ref,alt,EIDOLON_ORIGIN) set differs \
             between the native merge and bcftools+awk"
        );
    } else {
        eprintln!(
            "note: cancer_simulate.sh produced no truth VCF (bcftools/bgzip absent) \
             — skipping truth-VCF parity check"
        );
    }
}
