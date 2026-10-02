//! #496 — a config key a subcommand does not read must stop the run, not be ignored.
//!
//! CORRECTNESS CRITERION. Before #496, `tumor_mutation_model:` (for `tumor_model:`) ran to
//! completion and simulated the tumor with the germline model. If this works, every
//! subcommand that takes `-c` exits non-zero and names the unknown key, and each shipped
//! template is accepted as written. It is FALSIFIED by a template + bogus key that runs on
//! or fails without naming the key, or by any shipped template being rejected for a key.
//!
//! The templates carry placeholder paths, so the unchanged runs fail later, on a missing
//! file. The must-not-fire check is only that the failure is not the unknown-key one.

mod common;

use common::eidolon;
use std::path::Path;

const CASES: &[(&str, &str)] = &[
    ("gen-reads", "gen_reads_template.yml"),
    ("gen-cancer-reads", "gen_cancer_reads_template.yml"),
    ("gen-mut-model", "gen_mut_model_template.yml"),
    ("gen-seq-error-model", "gen_seq_error_model_template.yml"),
    (
        "gen-frag-length-model",
        "gen_frag_length_model_template.yml",
    ),
    ("gen-gc-bias-model", "gen_gc_bias_model.yml"),
    ("gen-bam-models", "gen_bam_models.yml"),
    ("compare-vcfs", "compare_vcfs_template.yml"),
    ("filter-reads", "filter_reads_template.yml"),
];

fn template(name: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("template_config")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Run `subcommand -c <yaml>` in a scratch directory; return (success, stdout + stderr).
fn run(subcommand: &str, yaml: &str) -> (bool, String) {
    let work = tempfile::tempdir().unwrap();
    let cfg = work.path().join("config.yml");
    std::fs::write(&cfg, yaml).unwrap();
    let out = eidolon()
        .current_dir(work.path())
        .args([subcommand, "-c"])
        .arg(&cfg)
        .output()
        .unwrap();
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    (out.status.success(), log)
}

#[test]
fn every_subcommand_rejects_an_unknown_key_by_name() {
    for (subcommand, file) in CASES {
        let yaml = format!("{}\ntotally_bogus_key: 42\n", template(file));
        let (ok, log) = run(subcommand, &yaml);
        assert!(!ok, "{subcommand} accepted an unknown key:\n{log}");
        assert!(
            log.contains(&format!("unknown config key(s) in {subcommand}")),
            "{subcommand} failed, but not on the unknown key:\n{log}"
        );
        assert!(log.contains("`totally_bogus_key`"), "{subcommand}:\n{log}");
    }
}

#[test]
fn every_shipped_template_is_accepted_as_written() {
    for (subcommand, file) in CASES {
        let (_, log) = run(subcommand, &template(file));
        assert!(
            !log.contains("unknown config key"),
            "{subcommand} rejected a key in template_config/{file}:\n{log}"
        );
    }
}

#[test]
fn the_issue_496_typo_names_the_key() {
    let yaml = template("gen_cancer_reads_template.yml")
        + "\ntumor_mutation_model: tools/cosmic_v104_pancancer_model.json.gz\n";
    let (ok, log) = run("gen-cancer-reads", &yaml);
    assert!(!ok, "the #496 typo ran:\n{log}");
    assert!(
        log.contains("unknown config key(s) in gen-cancer-reads: `tumor_mutation_model`"),
        "{log}"
    );
}
