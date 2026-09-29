//! The BAM-fitted substitution matrix, round-tripped through a real aligner (#752).
//!
//! CORRECTNESS CRITERION. `gen-seq-error-model`'s `bam_file:` exists so that simulated reads
//! carry the fitted library's substitution spectrum. If the fit works, planting a known matrix,
//! simulating reads with it, aligning them with a real aligner and fitting from the alignment
//! must give the planted matrix back, cell by cell. It is FALSIFIED by any cell off by more
//! than the tolerance.
//!
//! THE FRAME. Both paired writers generate R2 forward, draw its errors, and only then
//! reverse-complement the record (`fastq_tools.rs`, `reverse_complement_record`), so
//! substitution errors are applied in REFERENCE orientation. The BAM counter reads mismatches
//! in reference orientation too. The two must agree, and this test is what says so.
//!
//! WHY THE PLANTED MATRIX IS ASYMMETRIC. A matrix that equals its own reverse complement
//! (A->C = T->G, and so on) comes back unchanged whichever frame the counter uses, so it cannot
//! tell them apart. Here A->C is 0.80 and T->G 0.10. Measured: a counter that complements
//! reverse-strand reads (read orientation, i.e. disagreeing with generation) returns 0.454
//! for A->C, T->G off by 0.354, and this test fails.
//!
//! The denominator is asserted: the fit must see enough mismatches for the tolerance to mean
//! something, and simulation runs with mutation_rate 0, so every mismatch is a sequencing
//! error.
//!
//! Fixture: ecoli (4.6 Mb, one contig), which aligns unambiguously. Needs bwa-mem2, so it is
//! `#[ignore]` and runs in release-gates.

mod common;

use common::gate2::bwa_mem2;
use common::{GenReadsConfig, eidolon, fresh_workdir};
use std::path::{Path, PathBuf};
use std::process::Command;

const READ_LEN: usize = 150;

/// Planted substitution matrix, rows A, C, G, T; the diagonal is ignored.
const PLANTED: [[f64; 4]; 4] = [
    [0.0, 0.80, 0.10, 0.10],
    [0.70, 0.0, 0.20, 0.10],
    [0.10, 0.30, 0.0, 0.60],
    [0.30, 0.60, 0.10, 0.0],
];

fn ecoli() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
}

/// A FASTQ whose every base is Q20 (1% error), for the quality half of both models.
fn write_training_fastq(path: &Path) {
    let seq: String = "ACGT".chars().cycle().take(READ_LEN).collect();
    let qual = "5".repeat(READ_LEN);
    let body: String = (0..200)
        .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
        .collect();
    std::fs::write(path, body).unwrap();
}

fn fit(work: &Path, tag: &str, extra: &str) -> (PathBuf, String) {
    let fq = work.join("train.fastq");
    let model = work.join(format!("{tag}.json.gz"));
    let cfg = work.join(format!("{tag}.yml"));
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n{extra}",
            fq.display(),
            model.display()
        ),
    )
    .unwrap();
    let out = eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&cfg)
        .output()
        .unwrap();
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(out.status.success(), "fit {tag} failed:\n{log}");
    (model, log)
}

/// Row-by-row substitution probabilities out of a model file, rows and columns A, C, G, T.
fn matrix(model: &Path) -> [[f64; 4]; 4] {
    let file = std::fs::File::open(model).unwrap();
    let v: serde_json::Value = serde_json::from_reader(flate2::read::GzDecoder::new(file)).unwrap();
    let mut m = [[0.0; 4]; 4];
    for (r, key) in ["a", "c", "g", "t"].iter().enumerate() {
        let cum: Vec<f64> = v["transition_distros"][key]["weights"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        for c in 0..4 {
            m[r][c] = if c == 0 { cum[0] } else { cum[c] - cum[c - 1] };
        }
    }
    m
}

/// bwa-mem2's SAM, rewritten as the BAM `gen-seq-error-model` reads.
fn sam_to_bam(sam: &Path, bam: &Path) {
    use noodles::{bam as nbam, sam as nsam};
    let mut reader = nsam::io::reader::Builder::default()
        .build_from_path(sam)
        .unwrap();
    let header = reader.read_header().unwrap();
    let mut writer = nbam::io::Writer::new(std::fs::File::create(bam).unwrap());
    writer.write_header(&header).unwrap();
    for record in reader.records() {
        let record = record.unwrap();
        use noodles::sam::alignment::io::Write as _;
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

#[test]
#[ignore = "needs bwa-mem2; runs in release-gates"]
fn a_planted_substitution_matrix_survives_simulation_alignment_and_refit() {
    let bwa = bwa_mem2();
    let (_g, work) = fresh_workdir();
    write_training_fastq(&work.join("train.fastq"));

    let tsv = work.join("planted.tsv");
    let rows: String = PLANTED
        .iter()
        .map(|r| format!("{}\t{}\t{}\t{}\n", r[0], r[1], r[2], r[3]))
        .collect();
    std::fs::write(&tsv, format!("A\tC\tG\tT\n{rows}")).unwrap();
    let (planted_model, _) = fit(
        &work,
        "planted",
        &format!("transition_matrix_file: {}\n", tsv.display()),
    );

    let mut cfg = GenReadsConfig::new(ecoli(), work.clone(), "sim");
    cfg.paired_ended = true;
    cfg.read_len = READ_LEN;
    cfg.coverage = 5;
    cfg.fragment_mean = Some(350.0);
    cfg.fragment_st_dev = Some(30.0);
    cfg.mutation_rate = Some(0.0);
    cfg.rng_seed = "round-trip".to_string();
    cfg.sequence_error_model = Some(planted_model);
    let yaml = cfg.write_yaml();
    let out = eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "gen-reads failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let reference = work.join("ecoli.fa");
    std::fs::copy(ecoli(), &reference).unwrap();
    let idx = Command::new(&bwa)
        .arg("index")
        .arg(&reference)
        .output()
        .unwrap();
    assert!(
        idx.status.success(),
        "{}",
        String::from_utf8_lossy(&idx.stderr)
    );
    let aln = Command::new(&bwa)
        .args(["mem", "-t", "4"])
        .arg(&reference)
        .arg(work.join("sim_r1.fastq.gz"))
        .arg(work.join("sim_r2.fastq.gz"))
        .output()
        .unwrap();
    assert!(
        aln.status.success(),
        "{}",
        String::from_utf8_lossy(&aln.stderr)
    );
    let sam = work.join("sim.sam");
    std::fs::write(&sam, &aln.stdout).unwrap();
    let bam = work.join("sim.bam");
    sam_to_bam(&sam, &bam);

    let (fitted_model, log) = fit(&work, "fitted", &format!("bam_file: {}\n", bam.display()));

    // Denominator: ~0.8M bases at 1% error per contig pass at 5x is ~230k substitution
    // errors; require well over the count the tolerance needs.
    let observed: usize = log
        .lines()
        .find_map(|l| {
            l.split("Observed ")
                .nth(1)
                .and_then(|r| r.split(' ').next())
                .and_then(|n| n.parse().ok())
        })
        .expect("the fit did not report its mismatch count");
    assert!(
        observed > 50_000,
        "only {observed} mismatches: too few to test a 0.03 tolerance"
    );

    let fitted = matrix(&fitted_model);
    let bases = ['A', 'C', 'G', 'T'];
    let mut worst = (0.0f64, String::new());
    for r in 0..4 {
        eprintln!(
            "{}: planted {:?}  fitted [{:.3}, {:.3}, {:.3}, {:.3}]",
            bases[r], PLANTED[r], fitted[r][0], fitted[r][1], fitted[r][2], fitted[r][3]
        );
        for c in 0..4 {
            if r == c {
                continue;
            }
            let d = (fitted[r][c] - PLANTED[r][c]).abs();
            if d > worst.0 {
                worst = (d, format!("{}->{}", bases[r], bases[c]));
            }
        }
    }
    eprintln!(
        "mismatches fitted: {observed}; worst cell {} off by {:.3}",
        worst.1, worst.0
    );
    assert!(
        worst.0 < 0.03,
        "the fitted matrix does not reproduce the planted one: {} is off by {:.3} \
         (a counter in a different frame from generation averages each cell with its \
         reverse complement)",
        worst.1,
        worst.0
    );
}
