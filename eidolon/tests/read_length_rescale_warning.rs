//! #742 — generating at a read length the quality model was not fitted at must say so.
//!
//! CORRECTNESS CRITERION. `quality_index_remap` stretches or compresses a model fitted at N bp
//! onto reads of M bp. That is NEAT2's behavior, ported deliberately, and NEAT2 warned about it.
//! If this works, a run with M != N logs one warning naming both lengths, and a run with M == N
//! logs none. It is FALSIFIED by a mismatched run with no warning, a warning that omits either
//! length, a warning repeated per read or per contig, or any warning on a matched run.
//!
//! Fixture: H1N1. Nothing here depends on genome size.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::path::Path;

/// Run gen-reads and return everything it logged, stdout and stderr together.
fn run_log(work: &Path, read_len: usize, model: Option<&Path>, tag: &str) -> String {
    run_log_extra(work, read_len, model, tag, true, "")
}

/// As `run_log`, with extra YAML lines appended to the config.
fn run_log_extra(
    work: &Path,
    read_len: usize,
    model: Option<&Path>,
    tag: &str,
    paired: bool,
    extra: &str,
) -> String {
    let mut cfg = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), tag);
    cfg.paired_ended = paired;
    cfg.read_len = read_len;
    cfg.coverage = 2;
    cfg.fragment_mean = Some(400.0);
    cfg.fragment_st_dev = Some(30.0);
    cfg.sequence_error_model = model.map(Path::to_path_buf);
    let yaml = cfg.write_yaml();
    let mut body = std::fs::read_to_string(yaml.path()).unwrap();
    body.push_str(extra);
    std::fs::write(yaml.path(), body).unwrap();
    let out = eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .output()
        .unwrap();
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(out.status.success(), "gen-reads {tag} failed:\n{log}");
    log
}

fn rescale_warnings(log: &str) -> Vec<&str> {
    log.lines()
        .filter(|l| l.contains("WARN") && l.contains("rescal"))
        .collect()
}

#[test]
fn generating_at_another_read_length_warns_once_naming_both_lengths() {
    let (_g, work) = fresh_workdir();
    // The shipped default is fitted at 250 bp.
    let log = run_log(&work, 100, None, "mismatch");
    let warnings = rescale_warnings(&log);
    assert_eq!(
        warnings.len(),
        1,
        "expected exactly one rescale warning, got {}:\n{log}",
        warnings.len()
    );
    let w = warnings[0];
    assert!(
        w.contains("250") && w.contains("100"),
        "the warning must name the model's length (250) and the run's (100): {w}"
    );
}

/// Must not fire: a model fitted at the generated length.
#[test]
fn generating_at_the_fitted_read_length_does_not_warn() {
    let (_g, work) = fresh_workdir();
    let fq = work.join("train.fastq");
    let seq: String = "ACGT".chars().cycle().take(100).collect();
    let qual: String = "I".repeat(100);
    let mut body = String::new();
    for i in 0..200 {
        body.push_str(&format!("@r{i}\n{seq}\n+\n{qual}\n"));
    }
    std::fs::write(&fq, body).unwrap();
    let model = work.join("model100.json.gz");
    let cfg = work.join("fit.yml");
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n",
            fq.display(),
            model.display()
        ),
    )
    .unwrap();
    let fit = eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&cfg)
        .output()
        .unwrap();
    assert!(
        fit.status.success(),
        "{}",
        String::from_utf8_lossy(&fit.stderr)
    );

    let log = run_log(&work, 100, Some(&model), "matched");
    let warnings = rescale_warnings(&log);
    assert!(
        warnings.is_empty(),
        "no rescale warning expected:\n{}",
        warnings.join("\n")
    );
}

/// Long-read mode rescales every read to its own length, so a warning naming one read length
/// would be wrong; it must say the profile is rescaled per read, still once.
#[test]
fn long_read_mode_warns_once_that_each_read_is_rescaled() {
    let (_g, work) = fresh_workdir();
    let log = run_log_extra(&work, 250, None, "long", false, "long_reads: true\n");
    let warnings = rescale_warnings(&log);
    assert_eq!(warnings.len(), 1, "expected one rescale warning:\n{log}");
    assert!(
        warnings[0].contains("Long-read") && warnings[0].contains("250"),
        "the long-read warning must say so and name the model's length: {}",
        warnings[0]
    );
}
