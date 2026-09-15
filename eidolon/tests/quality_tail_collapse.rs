//! #694: eidolon's quality model cannot represent degraded reads.
//!
//! Real Illumina data is bimodal. Most reads hold high quality to the end; a minority collapse
//! partway through and STAY collapsed, because a cluster that loses sync does not recover.
//! `QualityScoreModel` is a first-order Markov chain over (position, previous score), which has
//! no per-read latent variable to carry that state, so long low excursions are exponentially
//! unlikely no matter what it is fitted on.
//!
//! Measured on HG002 (GIAB 2x250), reads whose last 50 bases average below a threshold:
//!
//! | source | tail Q<25 | tail Q<20 |
//! |---|---|---|
//! | real, library-wide | 9.80% | 4.45% |
//! | simulated, model fitted from those same reads | 0.36% | **0.00%** |
//!
//! THIS FILE REPRODUCES THAT LOCALLY, in seconds rather than a Delta campaign, so the fix can
//! be iterated on. The fixture below is calibrated against BOTH measured rates -- one target is
//! not enough, and calibrating against only Q<25 hid the fact that an earlier fixture's decay
//! was far too deep (nearly all its collapsed reads fell below Q20, where real data splits
//! about 2.18 to 1 between the thresholds).
//!
//! HOW CLOSE THIS GETS, measured at matched 250 bp read length:
//!
//! | threshold | fixture | generated | local | real HG002 |
//! |---|---|---|---|---|
//! | tail Q<25 | 9.92% | 2.68% | 3.7x short | ~27x short |
//! | tail Q<20 | 4.43% | 0.06% | 74x short  | 0.00% generated |
//!
//! So the DEEP population -- the one that ends up well below Q20 -- is what the chain
//! destroys, and that reproduces sharply. The milder Q<25 population survives here far better
//! than on real data, so this fixture UNDERSTATES the defect at that threshold. A fix that
//! satisfies this file has not thereby been shown to close #694 on HG002; only a refit there
//! shows that.
//!
//! Read length must match the model's. Generating 151 bp reads from a 250 bp model routes
//! through `quality_index_remap` and produced an apparent 10.8x gap that was an artifact of
//! the remap, not of the model form.
//!
//! THE TARGET, deliberately not written as an `#[ignore]`d test. In this repo `#[ignore]`
//! means "runs in release-gates" (`cargo test -- --ignored`), not "skipped" -- that job exists
//! because an ignored test nothing ever runs is not a gate (#582). A permanently-red ignored
//! test would therefore break it. The tripwire is `quality_model_currently_loses_the_degraded_
//! population` below: it pins the CURRENT numbers and fails the moment they change, which is
//! when whoever changed them should replace it with the assertion stated here.
//!
//! A model fitted on reads where a known fraction carry a collapsed tail must generate reads
//! carrying that same fraction, within a factor of two at both thresholds. Against a ~27x
//! shortfall, landing anywhere between half and double the fitted rate is a complete success.
//!
//! WHAT IS GROUNDED AND WHAT IS NOT. The healthy per-position profile is #694's own measured
//! table. The decay floor is Q17, which #692 measured independently as the mean quality of the
//! soft-clipped portion of real reads -- a grid search over the fixture landed on the same
//! value, from a different direction. The degraded fraction is fitted to the two rates above.
//!
//! The decay SHAPE, the onset distribution, and the transient-dip statistics are NOT measured;
//! they are chosen to reproduce the rates. `scripts/delta/measure_quality_degradation.sh`
//! exists to replace them with real numbers, and this fixture should be revisited when it has
//! been run. The transient dips matter more than their arbitrariness suggests: healthy reads
//! that visit low scores and RECOVER are what teach the chain to pull back toward the position
//! mean, and without them the chain reproduces degradation perfectly well.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Deterministic PRNG. Written out rather than pulled in so the fixture is reproducible from
/// this file alone and adds no dependency to a crate that ships a lean tree.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(6364136223846793005).wrapping_add(1))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn unit(&mut self) -> f64 {
        (self.next_u64() % 1_000_000) as f64 / 1_000_000.0
    }
    fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + (self.next_u64() as usize) % (hi - lo + 1)
    }
    /// Box-Muller. The tails matter here, so a sum-of-uniforms approximation will not do.
    fn gauss(&mut self, mean: f64, sd: f64) -> f64 {
        let u1 = self.unit().max(1e-9);
        let u2 = self.unit();
        mean + sd * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// #694's measured mean-quality-by-position table for real HG002 R1, interpolated.
/// Note the low 5' start (Q28.4): real reads begin poorly, and eidolon's do not.
const PROFILE: &[(usize, f64)] = &[
    (1, 28.4),
    (41, 34.9),
    (81, 36.1),
    (121, 36.8),
    (161, 36.3),
    (201, 35.2),
    (241, 33.3),
    (250, 32.8),
];

fn healthy_mean(pos: usize) -> f64 {
    if pos <= PROFILE[0].0 {
        return PROFILE[0].1;
    }
    for w in PROFILE.windows(2) {
        let ((x0, y0), (x1, y1)) = (w[0], w[1]);
        if pos <= x1 {
            return y0 + (y1 - y0) * (pos - x0) as f64 / (x1 - x0) as f64;
        }
    }
    PROFILE[PROFILE.len() - 1].1
}

/// Fixture parameters. Calibrated as a set; changing one without re-checking both rates in
/// `the_fixture_reproduces_the_measured_degradation_rates` invalidates the calibration.
// Calibrated against the LIBRARY-WIDE rates, not the head-sample: 9.80% of reads with tail
// Q<25 and 4.45% with Q<20, measured over 2,007,617 reads drawn at stride 65 across all
// 130,495,089 records of HG002 R1 (jobs 22080254 / 22080262). The head-sample figures this was
// first built on -- 12.22% / 5.61% -- were one flowcell tile and ~25% high.
//
// DECAY_FLOOR sits at the Q17 that #692 measured for clipped portions; 16.9 rather than 17.0
// is well inside that measurement, and the ratio between the two thresholds is what pins it.
const READ_LEN: usize = 250;
const DEGRADED_FRACTION: f64 = 0.129;
const DECAY_FLOOR: f64 = 16.9;
const DIP_RATE: f64 = 0.02;
const N_FIXTURE_READS: usize = 15_000;

fn encode(q: f64) -> u8 {
    (q.round().clamp(2.0, 40.0) as u8) + 33
}

/// Write a FASTQ whose degraded reads collapse at a per-read onset and never recover, and whose
/// healthy reads dip briefly and do recover.
fn write_fixture(path: &Path, n_reads: usize, seed: u64) {
    let mut rng = Lcg::new(seed);
    let file = std::fs::File::create(path).unwrap();
    let mut gz = GzEncoder::new(file, Compression::default());
    for i in 0..n_reads {
        let degraded = rng.unit() < DEGRADED_FRACTION;
        let onset = rng.range(READ_LEN * 15 / 100, READ_LEN - 5);
        let mut quals = Vec::with_capacity(READ_LEN);
        let mut in_dip = 0usize;
        for pos in 1..=READ_LEN {
            let base = healthy_mean(pos);
            let q = if degraded && pos >= onset {
                let t = (pos - onset) as f64 / ((READ_LEN - onset).max(1)) as f64;
                rng.gauss(base - (base - DECAY_FLOOR) * t, 3.0)
            } else if in_dip > 0 {
                in_dip -= 1;
                rng.gauss(15.0, 4.0)
            } else if rng.unit() < DIP_RATE {
                in_dip = rng.range(1, 4);
                rng.gauss(15.0, 4.0)
            } else {
                rng.gauss(base, 2.5)
            };
            quals.push(encode(q));
        }
        let seq: String = (0..READ_LEN)
            .map(|_| ['A', 'C', 'G', 'T'][(rng.next_u64() % 4) as usize])
            .collect();
        write!(
            gz,
            "@r{i}\n{seq}\n+\n{}\n",
            String::from_utf8(quals).unwrap()
        )
        .unwrap();
    }
    gz.finish().unwrap();
}

/// Percentage of reads whose last `window` bases average below 25 and below 20.
/// Returns (pct_lt_25, pct_lt_20, n_reads) -- the denominator is reported because a rate over an
/// unknown denominator is not a result.
fn tail_collapse_rates(path: &Path, window: usize) -> (f64, f64, usize) {
    // MultiGzDecoder, not GzDecoder: eidolon writes FASTQ in chunks, so the output is a
    // CONCATENATED gzip stream. GzDecoder stops at the end of the first member and yields a
    // clean, short, entirely plausible read set -- 456 reads of 37,000 here, which is how this
    // was caught. eidolon-core's own `read_gzip_lines` uses MultiGzDecoder for the same reason.
    let reader = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
    let (mut n, mut lt25, mut lt20) = (0usize, 0usize, 0usize);
    for (i, line) in reader.lines().enumerate() {
        if i % 4 != 3 {
            continue;
        }
        let line = line.unwrap();
        if line.len() < window {
            continue;
        }
        let tail: f64 = line.as_bytes()[line.len() - window..]
            .iter()
            .map(|&b| (b - 33) as f64)
            .sum::<f64>()
            / window as f64;
        n += 1;
        if tail < 25.0 {
            lt25 += 1;
        }
        if tail < 20.0 {
            lt20 += 1;
        }
    }
    assert!(n > 0, "no reads measured in {}", path.display());
    (
        100.0 * lt25 as f64 / n as f64,
        100.0 * lt20 as f64 / n as f64,
        n,
    )
}

/// Fit a model on the fixture, generate reads with it, and return the generated R1 path.
fn fit_and_generate(work: &Path, fixture: &Path) -> PathBuf {
    let model = work.join("model.json.gz");
    let cfg = work.join("fit.yml");
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n",
            fixture.display(),
            model.display()
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-seq-error-model", "-c", cfg.to_str().unwrap()])
        .assert()
        .success();

    let reference = PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ));
    let mut reads_cfg = GenReadsConfig::new(reference, work.to_path_buf(), "gen");
    reads_cfg.read_len = READ_LEN;
    reads_cfg.coverage = 2;
    reads_cfg.paired_ended = false;
    reads_cfg.produce_fastq = true;
    reads_cfg.rng_seed = "q694".to_string();
    reads_cfg.sequence_error_model = Some(model);
    let gen_cfg = reads_cfg.write_yaml();
    eidolon()
        .args(["gen-reads", "-c", gen_cfg.path().to_str().unwrap()])
        .assert()
        .success();

    work.join("gen_r1.fastq.gz")
}

/// The fixture itself must match what was measured on HG002, or nothing downstream means
/// anything. Both rates, because one does not constrain the shape: an earlier fixture hit
/// Q<25 while putting nearly every collapsed read below Q20, where real data splits 2.18 to 1.
#[test]
fn the_fixture_reproduces_the_measured_degradation_rates() {
    let (_g, work) = fresh_workdir();
    let fixture = work.join("fixture.fastq.gz");
    write_fixture(&fixture, N_FIXTURE_READS, 694);

    let (lt25, lt20, n) = tail_collapse_rates(&fixture, 50);
    assert_eq!(n, N_FIXTURE_READS, "every read must be measured");
    assert!(
        (8.5..11.5).contains(&lt25),
        "fixture tail Q<25 is {lt25:.2}%, outside the band around HG002's measured 9.80%"
    );
    assert!(
        (3.5..5.5).contains(&lt20),
        "fixture tail Q<20 is {lt20:.2}%, outside the band around HG002's measured 4.45%"
    );
    let ratio = lt25 / lt20;
    assert!(
        (1.7..2.8).contains(&ratio),
        "fixture splits {ratio:.2} to 1 between the thresholds; HG002 splits 2.20 to 1. A \
         fixture matching one rate but not the ratio has the wrong decay depth."
    );
}

/// CHARACTERIZATION of the #694 defect as it stands. This passes today and is expected to FAIL
/// when the two-component model lands -- at which point it should be replaced by the target
/// assertion in the file header. It exists so the gap cannot silently change size meanwhile.
#[test]
fn quality_model_currently_loses_the_degraded_population() {
    let (_g, work) = fresh_workdir();
    let fixture = work.join("fixture.fastq.gz");
    write_fixture(&fixture, N_FIXTURE_READS, 694);
    let (in25, in20, _) = tail_collapse_rates(&fixture, 50);

    let generated = fit_and_generate(&work, &fixture);
    let (out25, out20, n_out) = tail_collapse_rates(&generated, 50);

    assert!(
        n_out > 5_000,
        "only {n_out} reads generated; too few to rate"
    );
    eprintln!(
        "MEASURED  fixture Q<25 {in25:.2}% Q<20 {in20:.2}%  |  generated Q<25 {out25:.2}% Q<20 {out20:.2}%  n={n_out}"
    );
    // Thresholds set from measurement, not guessed. Measured on this fixture:
    //   Q<25  9.92% -> 2.68%   (3.7x short; real HG002 is ~27x)
    //   Q<20  4.43% -> 0.06%   (74x short; real HG002 reaches 0.00%)
    // The DEEP population is what the chain destroys, and that reproduces sharply here. The
    // milder Q<25 population survives locally far better than it does on real data, so this
    // fixture understates the defect at that threshold -- see the file header.
    assert!(
        out25 < in25 * 0.5,
        "#694 says the degraded population is lost. Fixture had {in25:.2}% of reads with a \
         collapsed tail; generation produced {out25:.2}%. If this now reads close to the \
         input, the model has been fixed -- replace this with the target assertion stated in \
         the file header, and delete this test."
    );
    assert!(
        out20 < in20 * 0.1,
        "the Q<20 population all but disappears (real HG002: 0.00%); got {out20:.2}% against \
         a fixture carrying {in20:.2}%"
    );
}
