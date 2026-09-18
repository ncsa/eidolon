//! #723 — a model fitted from two mate FASTQs must draw R1 reads from the R1 fit and R2 reads
//! from the R2 fit.
//!
//! CORRECTNESS CRITERION. Fit R1 from an all-Q40 library and R2 from an all-Q20 one. If this
//! works, generated R1 carries Q40 and generated R2 carries Q20 — computable from the fixtures
//! without reference to the implementation. It is FALSIFIED if the two output FASTQs carry the
//! same quality distribution (one model still serving both mates), or if swapping the two input
//! FASTQs does not swap the two outputs (the mates wired backwards, which every aggregate
//! measure over the pair would miss).
//!
//! WHY THE CRITERION IS THE PAIR AND NOT THE RATE. Measuring "R2 is worse than R1" would pass on
//! a model that simply made every read worse. Two disjoint, planted quality levels and an
//! assertion on *which* file carries which is what distinguishes a per-mate model from a
//! louder one.
//!
//! Fixture: H1N1 reference. This is a quality-model test — nothing here depends on genome size,
//! window statistics, or contig count, so ecoli's 4.6 Mb would only cost runtime.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Q40 = 'I', Q20 = '5' under Phred+33.
const R1_QUAL: char = 'I';
const R2_QUAL: char = '5';
const READ_LEN: usize = 60;

/// A FASTQ whose every base carries `qual`. One quality value means the fitted model has one
/// option and its output is that value, whatever the transition machinery does with it.
fn write_flat_fastq(path: &Path, n_reads: usize, qual: char) {
    let mut f = std::fs::File::create(path).unwrap();
    let seq = "ACGT".chars().cycle().take(READ_LEN).collect::<String>();
    let quals: String = std::iter::repeat_n(qual, READ_LEN).collect();
    for i in 0..n_reads {
        writeln!(f, "@r{i}\n{seq}\n+\n{quals}").unwrap();
    }
}

/// Fit a model from one or two mate FASTQs. `r2` absent is today's single-model behavior.
fn fit(work: &Path, r1: &Path, r2: Option<&Path>, tag: &str) -> PathBuf {
    let model = work.join(format!("model_{tag}.json.gz"));
    let cfg = work.join(format!("fit_{tag}.yml"));
    let r2_line = match r2 {
        Some(p) => format!("fastq_file_r2: {}\n", p.display()),
        None => String::new(),
    };
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n{r2_line}",
            r1.display(),
            model.display()
        ),
    )
    .unwrap();
    let out = eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&cfg)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "fitting {tag} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    model
}

/// Generate paired reads from `model` and return (R1 lines, R2 lines).
fn generate(work: &Path, model: &Path, tag: &str) -> (Vec<String>, Vec<String>) {
    let mut cfg = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), &format!("gen_{tag}"));
    cfg.paired_ended = true;
    cfg.read_len = READ_LEN;
    cfg.coverage = 6;
    cfg.rng_seed = "723".to_string();
    cfg.sequence_error_model = Some(model.to_path_buf());
    let yaml = cfg.write_yaml();
    let out = eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "generation {tag} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    (
        read_gzip_fastq_lines(&work.join(format!("gen_{tag}_r1.fastq.gz"))),
        read_gzip_fastq_lines(&work.join(format!("gen_{tag}_r2.fastq.gz"))),
    )
}

/// Mean quality over every base of every record, and the number of records it covered.
/// The denominator is returned rather than assumed: a mean over zero reads is not a result.
fn mean_quality(lines: &[String]) -> (f64, usize) {
    let mut sum = 0u64;
    let mut bases = 0u64;
    let mut records = 0usize;
    for q in lines.iter().skip(3).step_by(4) {
        records += 1;
        for b in q.bytes() {
            sum += (b - 33) as u64;
            bases += 1;
        }
    }
    assert!(bases > 0, "no quality bases found in {records} records");
    (sum as f64 / bases as f64, records)
}

#[test]
fn each_mate_is_drawn_from_its_own_fit() {
    let (_g, work) = fresh_workdir();
    let r1_fq = work.join("train_r1.fastq");
    let r2_fq = work.join("train_r2.fastq");
    write_flat_fastq(&r1_fq, 200, R1_QUAL);
    write_flat_fastq(&r2_fq, 200, R2_QUAL);

    let model = fit(&work, &r1_fq, Some(&r2_fq), "pair");
    let (out1, out2) = generate(&work, &model, "pair");
    let (m1, n1) = mean_quality(&out1);
    let (m2, n2) = mean_quality(&out2);

    eprintln!("R1 mean Q{m1:.3} over {n1} reads;  R2 mean Q{m2:.3} over {n2} reads");
    assert!(n1 > 50 && n2 > 50, "too few reads to rate: {n1} / {n2}");

    // Known answer: the fixtures carry exactly one quality value each.
    assert!(
        (m1 - 40.0).abs() < 0.5,
        "R1 was fitted from an all-Q40 library and must generate Q40; got Q{m1:.3}"
    );
    assert!(
        (m2 - 20.0).abs() < 0.5,
        "R2 was fitted from an all-Q20 library and must generate Q20; got Q{m2:.3}. \
         Q40 here means one model is still serving both mates."
    );
}

/// THE DECISION TEST. Swapping the two input FASTQs must swap the two outputs. Wiring the mates
/// backwards passes every aggregate measure over the pair — this is the only assertion that
/// distinguishes it.
#[test]
fn swapping_the_two_inputs_swaps_the_two_outputs() {
    let (_g, work) = fresh_workdir();
    let hi = work.join("hi.fastq");
    let lo = work.join("lo.fastq");
    write_flat_fastq(&hi, 200, R1_QUAL);
    write_flat_fastq(&lo, 200, R2_QUAL);

    let forward = fit(&work, &hi, Some(&lo), "fwd");
    let reverse = fit(&work, &lo, Some(&hi), "rev");
    let (f1, f2) = generate(&work, &forward, "fwd");
    let (r1, r2) = generate(&work, &reverse, "rev");

    let (fm1, _) = mean_quality(&f1);
    let (fm2, _) = mean_quality(&f2);
    let (rm1, _) = mean_quality(&r1);
    let (rm2, _) = mean_quality(&r2);
    eprintln!("forward R1 Q{fm1:.2} R2 Q{fm2:.2}  |  reversed R1 Q{rm1:.2} R2 Q{rm2:.2}");

    assert!(
        fm1 > fm2,
        "forward: R1 was fitted high and R2 low, so R1 must be the better file"
    );
    assert!(
        rm1 < rm2,
        "reversed: the inputs were swapped, so R2 must now be the better file. \
         R1 still better means the mates are wired backwards."
    );
}

/// MUST NOT FIRE: a fit with no second FASTQ produces a model with no R2 population, and paired
/// generation from it is byte-identical to a build that predates this feature. The frozen
/// baselines are the other half of this guarantee; this pins the paired path directly.
#[test]
fn a_single_fastq_fit_still_serves_both_mates() {
    let (_g, work) = fresh_workdir();
    let only = work.join("only.fastq");
    write_flat_fastq(&only, 200, R1_QUAL);

    let model = fit(&work, &only, None, "single");
    let (out1, out2) = generate(&work, &model, "single");
    let (m1, n1) = mean_quality(&out1);
    let (m2, n2) = mean_quality(&out2);
    eprintln!("single-fit R1 Q{m1:.3} ({n1} reads), R2 Q{m2:.3} ({n2} reads)");

    assert!(
        (m1 - 40.0).abs() < 0.5 && (m2 - 40.0).abs() < 0.5,
        "with one fitted library both mates must carry it: R1 Q{m1:.3}, R2 Q{m2:.3}"
    );
}
