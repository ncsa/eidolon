//! A bad config is user input, and user input never panics: the run stops with the
//! subcommand's own error, its message on stderr, and exit status 1.
//!
//! Every case here used to panic: Rust's panic exit status is 101, and the user saw
//! `thread 'main' panicked at …` plus a backtrace hint, sometimes after the real reason had
//! been logged and sometimes with nothing more than "Error generating configuration!".
//! Each case asserts the specific reason reaches stderr, not just that the run failed.

mod common;
use common::{eidolon, fresh_workdir, h1n1_reference};
use std::path::Path;

/// Runs `subcommand -c <config>` and returns (exit code, stderr).
fn run(work: &Path, subcommand: &str, config: &str) -> (Option<i32>, String) {
    let cfg = work.join("config.yml");
    std::fs::write(&cfg, config).unwrap();
    let out = eidolon()
        .current_dir(work)
        .arg("--log-dest")
        .arg(work.join("run.log"))
        .args([subcommand, "-c"])
        .arg(&cfg)
        .output()
        .unwrap();
    (
        out.status.code(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

fn assert_clean_error(code: Option<i32>, stderr: &str, reason: &str) {
    assert!(
        !stderr.contains("panicked"),
        "a config error must not panic:\n{stderr}"
    );
    assert_eq!(code, Some(1), "a config error exits 1, not 101:\n{stderr}");
    assert!(
        stderr.contains(reason),
        "stderr must name the reason ({reason}):\n{stderr}"
    );
}

#[test]
fn gen_reads_a_config_check_failure_names_its_reason() {
    let (_tmp, work) = fresh_workdir();
    // Paired-end with no fragment-length source: a check in check_and_log_config.
    let cfg = format!(
        "reference: {}\noutput_dir: {}\npaired_ended: true\nproduce_fastq: true\n",
        h1n1_reference().display(),
        work.display()
    );
    let (code, err) = run(&work, "gen-reads", &cfg);
    assert_clean_error(code, &err, "no fragment-length source");
}

#[test]
fn gen_reads_an_existing_output_is_refused_without_a_panic() {
    let (_tmp, work) = fresh_workdir();
    let cfg = format!(
        "reference: {}\noutput_dir: {}\noutput_prefix: x\nread_len: 50\ncoverage: 1\n\
         produce_fastq: true\noverwrite_output: false\n",
        h1n1_reference().display(),
        work.display()
    );
    std::fs::write(work.join("x_r1.fastq.gz"), "").unwrap();
    let (code, err) = run(&work, "gen-reads", &cfg);
    assert_clean_error(code, &err, "x_r1.fastq.gz");
}

#[test]
fn gen_reads_a_missing_config_file_names_gen_reads() {
    let (_tmp, work) = fresh_workdir();
    let out = eidolon()
        .arg("--log-dest")
        .arg(work.join("run.log"))
        .args(["gen-reads", "-c"])
        .arg(work.join("missing.yml"))
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert_clean_error(out.status.code(), &err, "missing.yml");
    // It used to be reported as a filter-reads error for every subcommand.
    assert!(!err.contains("filtering reads"), "{err}");
}

#[test]
fn gen_mut_model_a_missing_reference_names_it() {
    let (_tmp, work) = fresh_workdir();
    let cfg = format!(
        "reference: {}/nope.fa\nvcf_file: {}/nope.vcf\noutput_file: {}/m.json.gz\n",
        work.display(),
        work.display(),
        work.display()
    );
    let (code, err) = run(&work, "gen-mut-model", &cfg);
    assert_clean_error(code, &err, "nope.fa");
}

#[test]
fn gen_mut_model_a_missing_required_key_names_it() {
    let (_tmp, work) = fresh_workdir();
    let cfg = format!(
        "vcf_file: {}/x.vcf\noutput_file: {}/m.json.gz\n",
        work.display(),
        work.display()
    );
    let (code, err) = run(&work, "gen-mut-model", &cfg);
    assert_clean_error(code, &err, "reference");
}

#[test]
fn filter_reads_an_unknown_extension_names_it() {
    let (_tmp, work) = fresh_workdir();
    std::fs::write(work.join("regions.bed"), "H1N1_HA\t0\t100\n").unwrap();
    std::fs::write(work.join("reads.bam"), "x").unwrap();
    let cfg = format!(
        "bed_file: {d}/regions.bed\nfiles_to_filter:\n  - {d}/reads.bam\nfilter_key: .\n\
         overwrite_output: true\n",
        d = work.display()
    );
    let (code, err) = run(&work, "filter-reads", &cfg);
    assert_clean_error(code, &err, "reads.bam");
}

#[test]
fn filter_reads_an_unknown_key_names_it() {
    let (_tmp, work) = fresh_workdir();
    std::fs::write(work.join("regions.bed"), "H1N1_HA\t0\t100\n").unwrap();
    let cfg = format!(
        "bed_file: {d}/regions.bed\nfiles_to_filter: []\nfilter_key: .\noverwrite_output: true\n\
         filter_keyy: x\n",
        d = work.display()
    );
    let (code, err) = run(&work, "filter-reads", &cfg);
    assert_clean_error(code, &err, "filter_keyy");
}
