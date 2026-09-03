//! The #661 homopolymer indel-error curve must reach the generated reads.
//!
//! The unit tests in `sequencing_error_model` prove the curve scales a probability, and
//! those in `fastq_tools` prove a run length is measured correctly. Neither can show the
//! two are wired together: a model carrying a perfect curve is worth nothing if
//! `generate_read` never passes it any context. That gap is exactly #405's shape — unit
//! tests green while the chain they described was broken — so it gets its own test here.
//!
//! **Why these references.** Both arms are 0% GC. A poly-A reference and an alternating
//! `ATAT` reference differ ONLY in homopolymer structure, so the GC-bias model cannot
//! produce the contrast being measured. An earlier design contrasting poly-A against
//! cycled `ACGT` would have confounded 0% GC against 50% GC, and CLAUDE.md records two
//! separate defects in one day from coverage artifacts masquerading as a signal.
//!
//! **The known answer is a ratio.** Per aligned base,
//!
//! ```text
//! indel errors = P(error | quality) x indel_probability x E[curve over the reference]
//! ```
//!
//! `P(error)` and `indel_probability` are identical across arms and cancel, leaving the
//! ratio of the two references' mean curve values — computable from the published table
//! and the reference sequence alone, with no reference to the implementation. For poly-A
//! against `ATAT` that is 39.20 / 0.64 = 61.25x.

mod common;

use common::eidolon;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::{Path, PathBuf};

/// The curve from Delta job 21674484, restated here so the test does not import the value
/// it is checking. A copy that drifts from the shipped constant is the point: these must
/// be changed together, deliberately.
const CURVE: [f64; 10] = [
    0.64, 0.76, 0.82, 1.11, 1.58, 1.84, 5.64, 12.16, 24.24, 39.20,
];

/// Mean curve value over a sequence, computed by an independent, deliberately naive
/// run-length scan. This is the "correct output computable independently of the code under
/// test" the repo requires — it shares no code with `homopolymer_run_at`.
fn expected_mean_curve(seq: &[u8]) -> f64 {
    let mut total = 0.0;
    for i in 0..seq.len() {
        // Walk out from i in both directions. O(n * run), fine for a fixture.
        let mut run = 1usize;
        let mut j = i;
        while j > 0 && seq[j - 1] == seq[i] {
            run += 1;
            j -= 1;
        }
        let mut k = i;
        while k + 1 < seq.len() && seq[k + 1] == seq[i] {
            run += 1;
            k += 1;
        }
        total += CURVE[run.min(CURVE.len()) - 1];
    }
    total / seq.len() as f64
}

fn write_reference(dir: &Path, name: &str, seq: &[u8]) -> PathBuf {
    let path = dir.join(format!("{name}.fa"));
    let hdr = format!(">{name}\n");
    let body = String::from_utf8(seq.to_vec()).unwrap();
    std::fs::write(&path, format!("{hdr}{body}\n")).unwrap();
    std::fs::write(
        dir.join(format!("{name}.fa.fai")),
        format!(
            "{name}\t{}\t{}\t{}\t{}\n",
            seq.len(),
            hdr.len(),
            seq.len(),
            seq.len() + 1
        ),
    )
    .unwrap();
    path
}

/// Indel CIGAR ops and query-consuming bases across a produced BAM.
///
/// This is a golden BAM written by the simulator itself, not an aligner's opinion, so an
/// `I` or `D` op is a sequencing-error indel by construction — provided mutations are off,
/// which the config below enforces with `mutation_rate: 0.0`.
fn indel_ops_and_bases(path: &Path) -> (usize, usize, usize) {
    let file = std::fs::File::open(path).expect("bam not produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let (mut indels, mut bases, mut records) = (0usize, 0usize, 0usize);
    for r in reader.records() {
        let rec = r.unwrap();
        records += 1;
        for op in rec.cigar().iter() {
            let op = op.unwrap();
            match op.kind() {
                Kind::Insertion | Kind::Deletion => indels += 1,
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => bases += op.len(),
                _ => {}
            }
        }
    }
    (indels, bases, records)
}

/// Simulate over `seq` and return indel ops per aligned base.
fn indel_rate(dir: &Path, name: &str, seq: &[u8]) -> (f64, usize, usize) {
    let reference = write_reference(dir, name, seq);
    let cfg_text = format!(
        "reference: {ref}\nread_len: 100\ncoverage: 60\nploidy: 1\npaired_ended: false\n\
         mutation_rate: 0.0\n\
         produce_bam: true\nproduce_fastq: false\nproduce_vcf: false\n\
         overwrite_output: true\noutput_dir: {out}\noutput_filename: {name}\n\
         rng_seed: indel context fidelity\nnum_threads: 1\n",
        ref = reference.display(),
        out = dir.display(),
    );
    let cfg = dir.join(format!("{name}.yml"));
    std::fs::write(&cfg, cfg_text).unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    let (indels, bases, records) = indel_ops_and_bases(&dir.join(format!("{name}.bam")));
    assert!(
        records > 1000,
        "{name}: only {records} BAM records — too few to measure anything"
    );
    assert!(
        bases > 100_000,
        "{name}: only {bases} aligned bases — denominator too small to trust"
    );
    (indels as f64 / bases as f64, indels, bases)
}

const REF_LEN: usize = 120_000;

fn poly_a() -> Vec<u8> {
    vec![b'A'; REF_LEN]
}

fn alternating_at() -> Vec<u8> {
    (0..REF_LEN)
        .map(|i| if i % 2 == 0 { b'A' } else { b'T' })
        .collect()
}

/// Alternating 40 bp poly-A tracts and 40 bp `ATAT`, so run length VARIES with position
/// rather than being uniform across the whole reference.
fn mixed_tracts() -> Vec<u8> {
    (0..REF_LEN)
        .map(|i| {
            // Inside a tract every base is A; outside, A on even positions only.
            let in_tract = (i / 40) % 2 == 0;
            if in_tract || i % 2 == 0 { b'A' } else { b'T' }
        })
        .collect()
}

/// Run length at every position, by the same independent scan `expected_mean_curve` uses.
fn run_lengths(seq: &[u8]) -> Vec<usize> {
    let mut out = vec![0usize; seq.len()];
    let mut i = 0;
    while i < seq.len() {
        let mut j = i;
        while j < seq.len() && seq[j] == seq[i] {
            j += 1;
        }
        for slot in out.iter_mut().take(j).skip(i) {
            *slot = j - i;
        }
        i = j;
    }
    out
}

/// Reference positions of every indel CIGAR op in a BAM.
///
/// An `I` op consumes no reference, so it is attributed to the base it follows — the
/// position where `generate_read` actually drew the error.
fn indel_reference_positions(path: &Path) -> Vec<usize> {
    let file = std::fs::File::open(path).expect("bam not produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let mut positions = Vec::new();
    for r in reader.records() {
        let rec = r.unwrap();
        let Some(Ok(start)) = rec.alignment_start() else {
            continue;
        };
        let mut ref_pos = start.get() - 1;
        for op in rec.cigar().iter() {
            let op = op.unwrap();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => ref_pos += op.len(),
                Kind::Deletion => {
                    positions.push(ref_pos);
                    ref_pos += op.len();
                }
                Kind::Insertion => positions.push(ref_pos.saturating_sub(1)),
                Kind::Skip => ref_pos += op.len(),
                _ => {}
            }
        }
    }
    positions
}

#[test]
fn a_homopolymer_reference_gets_the_curves_indel_enrichment_in_its_reads() {
    let tmp = tempfile::tempdir().unwrap();

    let (rate_hp, n_hp, b_hp) = indel_rate(tmp.path(), "polya", &poly_a());
    let (rate_ctl, n_ctl, b_ctl) = indel_rate(tmp.path(), "atat", &alternating_at());

    // Known answer from the table alone: every poly-A base sits in a run at or past the
    // cap (39.20x); every ATAT base sits in a run of 1 (0.64x).
    let expected_hp = expected_mean_curve(&poly_a());
    let expected_ctl = expected_mean_curve(&alternating_at());
    assert!(
        (expected_hp - 39.20).abs() < 1e-9 && (expected_ctl - 0.64).abs() < 1e-9,
        "fixture is not the one this test reasons about: {expected_hp} / {expected_ctl}"
    );
    let expected_ratio = expected_hp / expected_ctl; // 61.25

    eprintln!(
        "[indelctx] polyA {n_hp} indels / {b_hp} bases = {rate_hp:.3e}\n\
         [indelctx] ATAT  {n_ctl} indels / {b_ctl} bases = {rate_ctl:.3e}\n\
         [indelctx] observed ratio {:.2}x, expected {expected_ratio:.2}x",
        rate_hp / rate_ctl
    );

    // Non-vacuity: a zero control makes the ratio meaningless, and a zero test arm would
    // pass any "not enriched" phrasing. Both must actually have events.
    assert!(
        n_ctl > 30,
        "control produced {n_ctl} indel errors — denominator too small for a ratio"
    );
    assert!(n_hp > 500, "homopolymer arm produced only {n_hp} indels");

    let observed_ratio = rate_hp / rate_ctl;
    // Tolerance is Poisson, not fitted: the control arm's count sets the noise floor, so
    // the band is +/- 4 sigma on the smaller count, widened to 40% only if that is tighter.
    let sigma = (1.0 / n_ctl as f64).sqrt() + (1.0 / n_hp as f64).sqrt();
    let tol = (4.0 * sigma).max(0.40);
    assert!(
        (observed_ratio - expected_ratio).abs() / expected_ratio < tol,
        "indel enrichment {observed_ratio:.2}x is not the curve's {expected_ratio:.2}x \
         (tolerance {:.0}%). Either the curve is not reaching generate_read, or it is \
         being applied with the wrong run length.",
        tol * 100.0
    );
}

#[test]
fn a_reference_with_varying_run_lengths_lands_where_its_composition_predicts() {
    // The uniform contrast above cannot catch a wrong POSITION: on a reference where every
    // base has the same run length, any index returns the same answer. This arm's run
    // length varies every 20 bp, so a curve applied at the wrong coordinate — or applied
    // as a constant — misses the predicted value.
    let tmp = tempfile::tempdir().unwrap();

    let (rate_mixed, n_mixed, _) = indel_rate(tmp.path(), "mixed", &mixed_tracts());
    let (rate_ctl, n_ctl, _) = indel_rate(tmp.path(), "atat2", &alternating_at());

    let expected_ratio =
        expected_mean_curve(&mixed_tracts()) / expected_mean_curve(&alternating_at());
    let observed_ratio = rate_mixed / rate_ctl;
    eprintln!(
        "[indelctx] mixed {n_mixed} indels, ratio to control {observed_ratio:.2}x, \
         expected {expected_ratio:.2}x"
    );

    assert!(n_ctl > 30 && n_mixed > 200, "too few events to compare");
    // The prediction must be a real discriminator, not a value both hypotheses satisfy.
    assert!(
        expected_ratio > 2.0,
        "fixture is too weak to distinguish anything: expected ratio only {expected_ratio:.2}"
    );
    let sigma = (1.0 / n_ctl as f64).sqrt() + (1.0 / n_mixed as f64).sqrt();
    let tol = (4.0 * sigma).max(0.40);
    assert!(
        (observed_ratio - expected_ratio).abs() / expected_ratio < tol,
        "mixed-composition arm landed at {observed_ratio:.2}x, not the {expected_ratio:.2}x \
         its sequence predicts (tolerance {:.0}%)",
        tol * 100.0
    );
}

/// Build a sequencing-error model file, then rewrite its curve to `curve`.
///
/// Patching the serialized model is what makes the must-not-fire below a real experiment:
/// it isolates the curve as the single variable on an otherwise identical run. It also
/// exercises the field's round trip through a real model file rather than through a
/// struct built in memory.
fn model_with_curve(dir: &Path, tag: &str, curve: Option<[f64; 10]>) -> PathBuf {
    use flate2::{Compression, read::GzDecoder, write::GzEncoder};
    use std::io::{Read, Write};

    // Train at a uniform Phred 8 (~0.158 per-base error). The shared synthetic-FASTQ
    // helper writes Q33-Q42, which yields so few errors that the comparison below runs on
    // counts in the tens — a denominator too small to distinguish 0.64 from anything.
    let fq = dir.join(format!("{tag}_train.fastq.gz"));
    {
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;
        let seq: String = "ACGT".chars().cycle().take(100).collect();
        let qual: String = std::iter::repeat_n(char::from(b'!' + 8), 100).collect();
        let mut enc = GzEncoder::new(std::fs::File::create(&fq).unwrap(), Compression::default());
        for i in 0..500 {
            writeln!(enc, "@read{i}\n{seq}\n+\n{qual}").unwrap();
        }
        enc.finish().unwrap();
    }
    let model = dir.join(format!("{tag}_model.json.gz"));
    let cfg = dir.join(format!("{tag}_build.yml"));
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\n",
            fq.display(),
            model.display()
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&cfg)
        .assert()
        .success();

    let Some(curve) = curve else { return model };

    let mut raw = String::new();
    GzDecoder::new(std::fs::File::open(&model).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let mut json: serde_json::Value = serde_json::from_str(&raw).unwrap();
    // Assert the field is actually there before replacing it. A silent no-op patch and a
    // real one look identical in the output, and CLAUDE.md records that exact trap.
    assert!(
        json.get("indel_context_curve").is_some(),
        "built model has no indel_context_curve field — patch would be a silent no-op"
    );
    json["indel_context_curve"] = serde_json::json!(curve.to_vec());
    let patched = dir.join(format!("{tag}_patched.json.gz"));
    let mut enc = GzEncoder::new(
        std::fs::File::create(&patched).unwrap(),
        Compression::default(),
    );
    enc.write_all(serde_json::to_string(&json).unwrap().as_bytes())
        .unwrap();
    enc.finish().unwrap();
    patched
}

/// Simulate over `seq` with an explicit model file.
fn indel_rate_with_model(dir: &Path, name: &str, seq: &[u8], model: &Path) -> (f64, usize) {
    let reference = write_reference(dir, name, seq);
    let cfg = dir.join(format!("{name}.yml"));
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 100\ncoverage: 60\nploidy: 1\npaired_ended: false\n\
             mutation_rate: 0.0\nsequence_error_model: {model}\n\
             produce_bam: true\nproduce_fastq: false\nproduce_vcf: false\n\
             overwrite_output: true\noutput_dir: {out}\noutput_filename: {name}\n\
             rng_seed: indel context fidelity\nnum_threads: 1\n",
            ref = reference.display(),
            model = model.display(),
            out = dir.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    let (indels, bases, _) = indel_ops_and_bases(&dir.join(format!("{name}.bam")));
    (indels as f64 / bases as f64, indels)
}

#[test]
fn flattening_the_curve_removes_the_effect_and_nothing_else() {
    // The must-not-fire, as a controlled experiment rather than an inequality that holds
    // for any input. Same reference, same seed, same everything — only the curve differs:
    // the shipped one against an all-1.0 curve, which is by definition "no context
    // effect". On `ATAT` every base is a run of 1, so the shipped curve must land at
    // exactly its first entry, 0.64x, relative to flat.
    //
    // This is also the test that would catch the curve being applied where it should not
    // be: if context leaked into the flat run, the two arms would converge on 1.0.
    let tmp = tempfile::tempdir().unwrap();
    let flat = model_with_curve(tmp.path(), "flat", Some([1.0; 10]));
    let shipped = model_with_curve(tmp.path(), "shipped", None);

    let (rate_flat, n_flat) =
        indel_rate_with_model(tmp.path(), "atat_flat", &alternating_at(), &flat);
    let (rate_curve, n_curve) =
        indel_rate_with_model(tmp.path(), "atat_curve", &alternating_at(), &shipped);

    let observed = rate_curve / rate_flat;
    eprintln!(
        "[indelctx] ATAT flat curve {n_flat} indels ({rate_flat:.3e}), \
         shipped curve {n_curve} indels ({rate_curve:.3e}) -> {observed:.3}x, expected 0.640x"
    );
    assert!(
        n_flat > 100 && n_curve > 30,
        "too few indels to compare (flat {n_flat}, curve {n_curve})"
    );
    let sigma = (1.0 / n_flat as f64).sqrt() + (1.0 / n_curve as f64).sqrt();
    let tol = (4.0 * sigma).max(0.30);
    assert!(
        (observed - 0.64).abs() / 0.64 < tol,
        "a run-1 reference must be suppressed to exactly the curve's first entry (0.64x) \
         against a flat curve; got {observed:.3}x (tolerance {:.0}%)",
        tol * 100.0
    );
    // Direction, stated separately: suppression, not enrichment. A run through #378's
    // VARIANT curve would push this above 1.0 instead.
    assert!(
        rate_curve < rate_flat,
        "run-1 context must suppress indel errors, not enrich them"
    );
}

#[test]
fn indel_errors_land_in_the_homopolymers_not_merely_at_the_right_overall_rate() {
    // THE positional test, and the one the ticket's acceptance criterion actually asks
    // for: "simulated candidate sites show the same homopolymer enrichment".
    //
    // Every other test here measures a rate averaged over the whole reference, and a mean
    // is blind to placement. This was not hypothetical — a mutant that measured the run
    // length at a FIXED index instead of the current base survived all three of them,
    // because `sequence` is the fragment, so a fixed index still samples a different
    // reference position per fragment and the average comes out unchanged. Only the
    // distribution over positions can tell those apart.
    //
    // The statistic is the one `indel_context_summarise.awk:49-52` computes, so a
    // simulated run and a real one are scored the same way:
    //
    //     enrichment(r) = (share of indels at run r) / (share of reference bases at run r)
    //
    // Since indels at run r are proportional to curve(r) x bases(r), that reduces to
    // curve(r) / E[curve] — predictable from the table and the reference alone.
    let tmp = tempfile::tempdir().unwrap();
    let seq = mixed_tracts();
    let reference = write_reference(tmp.path(), "posmix", &seq);
    let cfg = tmp.path().join("posmix.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 100\ncoverage: 60\nploidy: 1\npaired_ended: false\n\
             mutation_rate: 0.0\n\
             produce_bam: true\nproduce_fastq: false\nproduce_vcf: false\n\
             overwrite_output: true\noutput_dir: {out}\noutput_filename: posmix\n\
             rng_seed: indel context position\nnum_threads: 1\n",
            ref = reference.display(),
            out = tmp.path().display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();

    let runs = run_lengths(&seq);
    let positions = indel_reference_positions(&tmp.path().join("posmix.bam"));
    assert!(
        positions.len() > 2000,
        "only {} indel ops — too few to stratify",
        positions.len()
    );

    // Stratify into "in a long homopolymer" vs "not", at the run length where the measured
    // curve turns sharply (7 -> 5.64x, the first entry well clear of 1.0).
    const LONG: usize = 7;
    let bases_long = runs.iter().filter(|&&r| r >= LONG).count();
    let bases_share = bases_long as f64 / runs.len() as f64;
    let indels_long = positions
        .iter()
        .filter(|&&p| p < runs.len() && runs[p] >= LONG)
        .count();
    let indels_share = indels_long as f64 / positions.len() as f64;

    // Prediction, computed from the table and the sequence, not from the simulator.
    let mean_curve = expected_mean_curve(&seq);
    let expected_enrichment = {
        let mut weighted = 0.0;
        for (i, &r) in runs.iter().enumerate() {
            let _ = i;
            if r >= LONG {
                weighted += CURVE[r.min(CURVE.len()) - 1];
            }
        }
        (weighted / runs.len() as f64) / mean_curve / bases_share
    };
    let observed_enrichment = indels_share / bases_share;

    eprintln!(
        "[indelctx] positional: {indels_long}/{} indels in runs >= {LONG}; \
         background share {bases_share:.4}; enrichment {observed_enrichment:.2}x, \
         expected {expected_enrichment:.2}x",
        positions.len()
    );

    // Non-vacuity, both sides. A background share of 0 or 1 makes the ratio undefined or
    // trivially 1, and either would let this pass while measuring nothing.
    assert!(
        bases_share > 0.05 && bases_share < 0.95,
        "fixture background share {bases_share:.4} cannot support an enrichment ratio"
    );
    assert!(
        expected_enrichment > 1.5,
        "fixture predicts only {expected_enrichment:.2}x — too weak to discriminate"
    );

    let sigma = (1.0 / indels_long as f64).sqrt() + (1.0 / positions.len() as f64).sqrt();
    let tol = (4.0 * sigma).max(0.20);
    assert!(
        (observed_enrichment - expected_enrichment).abs() / expected_enrichment < tol,
        "indels are not landing where the curve puts them: {observed_enrichment:.2}x \
         against a predicted {expected_enrichment:.2}x (tolerance {:.0}%). The overall \
         rate can be right while the placement is wrong — that is what this test exists \
         to separate.",
        tol * 100.0
    );
}
