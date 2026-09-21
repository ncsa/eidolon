//! The SHIPPED DEFAULT, end to end: does a plain paired run actually draw from both of its
//! degraded populations?
//!
//! CORRECTNESS CRITERION. The default carries a degraded population per mate, each with a
//! `read_fraction` recorded in the model file. The class draw is per read, so the fraction of
//! generated reads whose tail collapses must equal that number — and R2's must exceed R1's,
//! because the library it was fitted from is measurably worse on the second mate. Both
//! expectations are read out of the model at runtime, so this is a known answer computable
//! without reference to the generator.
//!
//! WHY THIS EXISTS SEPARATELY FROM THE UNIT TESTS. `the_shipped_default_is_the_hg002_fit`
//! proves the fields DESERIALIZE. A model can carry a perfectly good degraded tensor that
//! nothing ever samples: the per-read class draw sits inside the `Some` branch of
//! `generate_quality_scores`, and a model whose degradation loaded but was never drawn from
//! would pass every metadata assertion while emitting reads from the healthy arm alone. The
//! only way to tell those apart is to generate reads and count.
//!
//! It is FALSIFIED if either mate's collapsed-tail rate lands near zero (the degraded arm is
//! not being sampled), if the two mates come out equal (one model is serving both), or if the
//! rates drift away from what the model declares.
//!
//! NO MODELS ARE CONFIGURED. That is the point — this is what a user gets out of the box.
//!
//! Fixture: ecoli, and `read_len` 250 to match the model's own fitted length, so this measures
//! the model rather than the read-length rescale (#742, measured separately).

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, read_gzip_fastq_lines};
use eidolon_core::models::sequencing_error_model::SequencingErrorModel;
use std::path::{Path, PathBuf};

const READ_LEN: usize = 250;
/// The same window the model was fitted with (`degradation_tail_window`).
const TAIL: usize = 50;

fn ecoli() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
}

/// Percentage of reads whose last `TAIL` bases average below `cut`, and the denominator.
fn collapsed_tail(path: &Path, cut: f64) -> (f64, usize) {
    let lines = read_gzip_fastq_lines(path);
    let mut n = 0usize;
    let mut low = 0usize;
    for q in lines.iter().skip(3).step_by(4) {
        if q.len() != READ_LEN {
            continue;
        }
        n += 1;
        let tail: u32 = q.bytes().rev().take(TAIL).map(|b| (b - 33) as u32).sum();
        if f64::from(tail) / TAIL as f64 > cut {
            continue;
        }
        low += 1;
    }
    assert!(n > 5_000, "only {n} full-length reads; too few to rate");
    (100.0 * low as f64 / n as f64, n)
}

#[test]
fn a_plain_paired_run_samples_both_degraded_populations() {
    // What the model DECLARES. Read at runtime so a future reship is checked against its own
    // numbers rather than against these.
    let model = SequencingErrorModel::default().expect("the shipped default must load");
    let r1_declared = model
        .quality_score_model()
        .degradation
        .as_ref()
        .expect("the shipped default carries an R1 degraded population")
        .read_fraction
        * 100.0;
    let r2_declared = model
        .quality_score_model_r2()
        .expect("the shipped default carries an R2 mate")
        .degradation
        .as_ref()
        .expect("R2 carries its own degraded population")
        .read_fraction
        * 100.0;

    let (_g, work) = fresh_workdir();
    let mut cfg = GenReadsConfig::new(ecoli(), work.clone(), "plain");
    cfg.read_len = READ_LEN;
    cfg.coverage = 2;
    cfg.paired_ended = true;
    // A fragment must hold both mates, or this measures adapter readthrough instead.
    cfg.fragment_mean = Some(600.0);
    cfg.fragment_st_dev = Some(60.0);
    cfg.mutation_rate = Some(0.0);
    cfg.rng_seed = "default-model-end-to-end".to_string();
    // Deliberately no sequence_error_model and no quality_score_model: the shipped defaults.
    let yaml = cfg.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let (r1_25, n1) = collapsed_tail(&work.join("plain_r1.fastq.gz"), 25.0);
    let (r2_25, n2) = collapsed_tail(&work.join("plain_r2.fastq.gz"), 25.0);
    let (r1_20, _) = collapsed_tail(&work.join("plain_r1.fastq.gz"), 20.0);
    let (r2_20, _) = collapsed_tail(&work.join("plain_r2.fastq.gz"), 20.0);
    eprintln!(
        "declared R1 {r1_declared:.2}% R2 {r2_declared:.2}%  |  generated Q<25 R1 {r1_25:.2}% \
         ({n1} reads) R2 {r2_25:.2}% ({n2} reads)  |  Q<20 R1 {r1_20:.2}% R2 {r2_20:.2}%"
    );

    // 1. THE DEGRADED ARM IS SAMPLED, not merely deserialized. A model whose degraded tensor
    //    loaded but was never drawn from lands near the single-population rate, which was
    //    measured at 0.54% on this library. Either mate near zero means the class draw is not
    //    happening.
    assert!(
        r1_25 > 5.0,
        "R1 collapsed-tail rate is {r1_25:.2}%. The degraded population is in the model but \
         does not appear to be sampled — a single-population fit of this library reads 0.54%."
    );
    assert!(
        r2_25 > 10.0,
        "R2 collapsed-tail rate is {r2_25:.2}%, far below its declared {r2_declared:.2}%"
    );

    // 2. KNOWN ANSWER: the rate is the model's own declared fraction, because the class draw
    //    is per read at exactly that probability.
    assert!(
        (r1_25 - r1_declared).abs() < 3.0,
        "R1 generated {r1_25:.2}% against a declared {r1_declared:.2}%"
    );
    assert!(
        (r2_25 - r2_declared).abs() < 3.0,
        "R2 generated {r2_25:.2}% against a declared {r2_declared:.2}%"
    );

    // 3. THE DECISION ASSERTION: the mates stay apart. One model serving both, or the two
    //    wired together, passes every aggregate measure over the pair.
    assert!(
        r2_25 > r1_25 * 1.5,
        "R2 must be materially worse than R1: {r2_25:.2}% against {r1_25:.2}%. Equal rates \
         mean one quality model is serving both mates."
    );

    // 4. The DEEP tail is the discriminating threshold (#694): a single-population fit of this
    //    library reads 0.00% here, so any material rate proves the degraded arm reached it.
    assert!(
        r1_20 > 1.0 && r2_20 > 1.0,
        "the deep tail is where a single-population model reads 0.00%: R1 {r1_20:.2}%, \
         R2 {r2_20:.2}%"
    );
    assert!(
        r2_20 > r1_20,
        "R2's deep tail must exceed R1's: {r2_20:.2}% against {r1_20:.2}%"
    );
}
