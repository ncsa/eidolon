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
        "reference: {}\noutput_dir: {}\noutput_basename: x\nread_len: 50\ncoverage: 1\n\
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

/// v5.0.0 keeps one name per key (#879). A removed name stops the run and names its
/// replacement, rather than being read under the old name or rejected without a hint.
#[test]
fn gen_reads_output_filename_is_refused_naming_output_basename() {
    let (_tmp, work) = fresh_workdir();
    let cfg = format!(
        "reference: {}\noutput_dir: {}\noutput_filename: x\nread_len: 50\ncoverage: 1\n",
        h1n1_reference().display(),
        work.display()
    );
    let (code, err) = run(&work, "gen-reads", &cfg);
    assert_clean_error(
        code,
        &err,
        "`output_filename` (removed in v5.0.0; use `output_basename`)",
    );
}

#[test]
fn gen_cancer_reads_removed_keys_are_refused_naming_their_replacements() {
    for (removed, replacement) in [
        ("output_filename", "output_basename"),
        ("coverage", "total_coverage"),
        ("rng_seed_root", "rng_seed"),
    ] {
        let (_tmp, work) = fresh_workdir();
        let cfg = format!(
            "reference: {}\noutput_dir: {}\npurity: 0.5\n{removed}: x\n",
            h1n1_reference().display(),
            work.display()
        );
        let (code, err) = run(&work, "gen-cancer-reads", &cfg);
        assert_clean_error(
            code,
            &err,
            &format!("`{removed}` (removed in v5.0.0; use `{replacement}`"),
        );
    }
}

/// `output_prefix` was renamed `output_basename` in v5.0.0; the old name must say so.
#[test]
fn output_prefix_is_refused_naming_output_basename() {
    for sub in ["gen-reads", "gen-cancer-reads"] {
        let (_tmp, work) = fresh_workdir();
        let cfg = format!(
            "reference: {}\noutput_dir: {}\n{extra}output_prefix: x\n",
            h1n1_reference().display(),
            work.display(),
            extra = if sub == "gen-cancer-reads" {
                "purity: 0.5\n"
            } else {
                ""
            }
        );
        let (code, err) = run(&work, sub, &cfg);
        assert_clean_error(
            code,
            &err,
            "`output_prefix` (renamed `output_basename` in v5.0.0)",
        );
    }
}

/// `output_basename` is a file name: a path separator would silently replace `output_dir`
/// (absolute) or fail with a bare I/O error (relative subdirectory).
#[test]
fn output_basename_with_a_path_separator_is_refused() {
    for sub in ["gen-reads", "gen-cancer-reads"] {
        for bad in ["sub/run", "/abs/run", "a\\\\b"] {
            let (_tmp, work) = fresh_workdir();
            let cfg = format!(
                "reference: {}\noutput_dir: {}\n{extra}output_basename: \"{bad}\"\n",
                h1n1_reference().display(),
                work.display(),
                extra = if sub == "gen-cancer-reads" {
                    "purity: 0.5\n"
                } else {
                    ""
                }
            );
            let (code, err) = run(&work, sub, &cfg);
            assert_clean_error(code, &err, "contains a path separator");
            assert!(err.contains("output_dir"), "{sub} {bad}: {err}");
        }
    }
}

/// filter-reads no longer reads the pre-v3.0.0 `RNEAT_generated_` prefix (#879). Such a read
/// is an error naming the prefix, never an ordinary line: dropping it silently would leave an
/// empty filtered FASTQ at exit 0.
#[test]
fn filter_reads_a_pre_v3_read_name_is_refused() {
    let (_tmp, work) = fresh_workdir();
    std::fs::write(work.join("regions.bed"), "chr1\t0\t2000\n").unwrap();
    std::fs::write(
        work.join("reads.fastq"),
        "@RNEAT_generated_chr1_0000001000_0000002000_0000000000000000/1\nACGT\n+\nIIII\n",
    )
    .unwrap();
    let cfg = format!(
        "bed_file: {d}/regions.bed\nfiles_to_filter:\n  - {d}/reads.fastq\nfilter_key: .\n\
         overwrite_output: true\n",
        d = work.display()
    );
    let (code, err) = run(&work, "filter-reads", &cfg);
    assert_clean_error(code, &err, "RNEAT_generated_");
    assert_clean_error(code, &err, "EIDOLON_generated_");
}
