//! Differential test: `eidolon validate` must agree with the tools that consume our files.
//!
//! This is the feature's correctness criterion, stated before it was implemented: for
//! every corpus file, our accept/reject verdict must match the tool's. Any disagreement
//! falsifies it. Without this the validator would only be an opinion about what is valid.
//!
//! The expected verdicts live in `observed_tool_behaviour.tsv` and were CAPTURED by
//! running the real tools, not written from documentation. Verdicts were identical
//! across samtools 1.18 / 1.19.2 / 1.22.1 and bcftools 1.17 / 1.19 / 1.22 (1.22 being
//! what Delta runs), so this asserts on verdicts rather than message text — messages are
//! what drift between versions.
//!
//! `validate_matches_live_tools` re-runs the actual binaries and is `#[ignore]`d, since
//! CI has neither installed. Run it wherever they are:
//!     cargo test --test validate_differential -- --ignored

mod common;

use common::eidolon;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn corpus_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_data/validation_corpus")
}

/// file -> "reject" if ANY recorded operation rejects it, else "accept".
fn recorded_verdicts() -> BTreeMap<String, String> {
    let tsv = std::fs::read_to_string(corpus_root().join("observed_tool_behaviour.tsv"))
        .expect("corpus manifest missing");
    let mut out: BTreeMap<String, String> = BTreeMap::new();
    for line in tsv.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        // A file is rejected if any single operation rejects it — `bcftools view`
        // tolerating unsorted records does not make an unsorted file usable, because
        // `tabix` still refuses to index it.
        let entry = out
            .entry(cols[0].to_string())
            .or_insert_with(|| "accept".to_string());
        if cols[2] == "reject" {
            *entry = "reject".to_string();
        }
    }
    assert!(!out.is_empty(), "manifest parsed to nothing");
    out
}

/// Run `eidolon validate` and return `(verdict, combined_output)`. A non-zero exit means
/// at least one ERROR-severity finding, which is our reject.
fn run_validate(args: &[&str]) -> (String, String) {
    let out = eidolon().arg("validate").args(args).output().unwrap();
    let mut text = String::from_utf8_lossy(&out.stdout).into_owned();
    text.push_str(&String::from_utf8_lossy(&out.stderr));
    let verdict = if out.status.success() {
        "accept".to_string()
    } else {
        "reject".to_string()
    };
    (verdict, text)
}

#[test]
fn validate_agrees_with_recorded_tool_verdicts() {
    let root = corpus_root();
    let mut checked = 0;
    let mut disagreements = Vec::new();
    for (file, expected) in recorded_verdicts() {
        let path = root.join(&file);
        assert!(path.exists(), "corpus file missing: {file}");
        let (ours, _) = run_validate(&[path.to_str().unwrap()]);
        if ours != expected {
            disagreements.push(format!("  {file}: eidolon={ours}, tool={expected}"));
        }
        checked += 1;
    }
    assert!(
        disagreements.is_empty(),
        "eidolon validate disagrees with the tools on {} of {checked} file(s):\n{}",
        disagreements.len(),
        disagreements.join("\n")
    );
    // Coverage of the harness's own input, not just the metric: a manifest that parsed
    // to two rows would otherwise "pass" while testing almost nothing.
    assert!(
        checked >= 14,
        "only {checked} corpus file(s) checked — the manifest is not being read"
    );
}

#[test]
fn every_rejection_names_the_tool_that_rejects_it() {
    // An ERROR that does not say what breaks is an unfalsifiable opinion. Each rejected
    // corpus file must name a tool and quote its message.
    let root = corpus_root();
    let mut seen = 0;
    for (file, expected) in recorded_verdicts() {
        if expected != "reject" {
            continue;
        }
        let (_, text) = run_validate(&[root.join(&file).to_str().unwrap()]);
        assert!(
            text.contains("rejects this:"),
            "{file}: rejected without naming the tool that rejects it:\n{text}"
        );
        assert!(
            text.contains("samtools") || text.contains("bcftools") || text.contains("tabix"),
            "{file}: citation names no known tool:\n{text}"
        );
        seen += 1;
    }
    assert!(seen >= 8, "only {seen} rejections examined");
}

#[test]
fn the_silent_data_loss_case_says_no_tool_will_catch_it() {
    // The most dangerous category and the reason this subcommand exists: bcftools
    // converts a type-mismatched INFO value to `.` without a word, which is exactly how
    // AF=AF=0.3000 propagated and then vanished. The warning must say so, because a
    // reader who assumes "no error means fine" is the person this is written for.
    let (verdict, text) = run_validate(&[corpus_root()
        .join("vcf/malformed_val.vcf")
        .to_str()
        .unwrap()]);
    assert_eq!(
        verdict, "accept",
        "must not reject what every tool accepts — that would break verdict parity"
    );
    assert!(
        text.contains("Type=Integer") && text.contains("AF=0.3"),
        "the finding must name the tag's declared type and the offending value:\n{text}"
    );
    assert!(
        text.contains("No tool rejects this"),
        "must state that nothing downstream will catch it:\n{text}"
    );
}

#[test]
fn well_formed_files_produce_no_findings() {
    // The must-NOT-fire case. A validator that rejected everything would satisfy every
    // rejection assertion above.
    for good in ["fastq/good.fq", "vcf/good.vcf"] {
        let (verdict, text) = run_validate(&[corpus_root().join(good).to_str().unwrap()]);
        assert_eq!(verdict, "accept", "{good} was rejected:\n{text}");
        assert!(
            !text.contains("ERROR") && !text.contains("WARNING"),
            "{good} is well-formed but produced findings:\n{text}"
        );
    }
}

#[test]
fn format_override_beats_a_misleading_extension() {
    let tmp = tempfile::tempdir().unwrap();
    let odd = tmp.path().join("actually_a_vcf.fq");
    std::fs::copy(corpus_root().join("vcf/unsorted.vcf"), &odd).unwrap();
    let (verdict, text) = run_validate(&["--format", "vcf", odd.to_str().unwrap()]);
    assert_eq!(
        verdict, "reject",
        "explicit --format vcf should apply:\n{text}"
    );
    assert!(
        text.contains("unsorted"),
        "should have found the sortedness error:\n{text}"
    );
}

#[test]
fn unknown_extension_is_an_error_not_a_silent_pass() {
    let tmp = tempfile::tempdir().unwrap();
    let odd = tmp.path().join("mystery.dat");
    std::fs::write(&odd, "not a recognised format\n").unwrap();
    let (verdict, text) = run_validate(&[odd.to_str().unwrap()]);
    assert_eq!(
        verdict, "reject",
        "an unrecognised file must fail loudly rather than report OK:\n{text}"
    );
}

/// Re-run the real tools and confirm the recorded verdicts have not drifted with a new
/// release. Ignored because CI has neither samtools nor bcftools.
#[test]
#[ignore]
fn validate_matches_live_tools() {
    use std::process::Command;
    let root = corpus_root();
    let mut mismatches = Vec::new();
    for (file, expected) in recorded_verdicts() {
        let path = root.join(&file);
        let live = if file.starts_with("bam/") {
            // BAM needs all three operations: they disagree with each other, and a
            // file is only usable if every one of them accepts it.
            let ok = ["quickcheck", "view", "index"].iter().all(|op| {
                let mut c = Command::new("samtools");
                c.arg(op);
                if *op == "view" {
                    c.args(["-o", "/dev/null"]);
                }
                let _ = std::fs::remove_file(path.with_extension("bam.bai"));
                c.arg(&path)
                    .output()
                    .expect("samtools not on PATH — this test needs it")
                    .status
                    .success()
            });
            if ok { "accept" } else { "reject" }
        } else if file.starts_with("fastq/") {
            let out = Command::new("samtools")
                .args(["import", "-0"])
                .arg(&path)
                .output()
                .expect("samtools not on PATH — this test needs it");
            if out.status.success() {
                "accept"
            } else {
                "reject"
            }
        } else {
            let out = Command::new("bcftools")
                .args(["view", "-o", "/dev/null"])
                .arg(&path)
                .output()
                .expect("bcftools not on PATH — this test needs it");
            if out.status.success() {
                "accept"
            } else {
                "reject"
            }
        };
        // The manifest folds several operations together, so only the accepts are
        // directly comparable to the single operation re-run here.
        if expected == "accept" && live != "accept" {
            mismatches.push(format!("  {file}: recorded accept, live {live}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "recorded tool behaviour has drifted:\n{}",
        mismatches.join("\n")
    );
}
