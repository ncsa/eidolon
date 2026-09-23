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
//! | tail Q<25 | 9.84% | 1.80% | 5.5x short | ~27x short |
//! | tail Q<20 | 4.37% | 0.07% | 62x short  | 0.00% generated |
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
//! WHAT IS GROUNDED AND WHAT IS NOT, after job 22093395 measured the library:
//!
//! MEASURED -- the per-position profile, the mean depth of a low run (Q15.5), the mean run
//! length (2.19 bases), and the head-of-read low-base rate of each population separately
//! (0.398% healthy, 8.186% collapsed). Those last two are what make the degraded population a
//! PROPENSITY rather than an onset, and they were measured, not chosen.
//!
//! TUNED -- the degraded fraction, the dip rate at the 3' end, the spread of per-read severity,
//! and how much runs lengthen toward the end. Four parameters against five constraints (both
//! thresholds, their ratio, and the two head-window statistics), so there is less freedom here
//! than the count suggests, but it is still fitting.
//!
//! NOT REPRODUCED -- `longest run in a collapsed read` reads 17.5 bases against a measured
//! 20.47, and low bases inside the tail window 12.3% against 16.6%. Both err in the same
//! direction: the fixture is slightly LESS extreme than the library. Left as is rather than
//! tuned away, because the residual is honest and chasing it would be fitting one library's
//! noise.

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
    /// Box-Muller. The tails matter here, so a sum-of-uniforms approximation will not do.
    fn gauss(&mut self, mean: f64, sd: f64) -> f64 {
        let u1 = self.unit().max(1e-9);
        let u2 = self.unit();
        mean + sd * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// Mean quality by position for real HG002 R1, interpolated. MEASURED over 2,007,617 reads
/// sampled at stride 65 across the whole library (job 22093395) -- not the table in #694, which
/// does not reproduce: it records Q28.4 at position 1 and a mid-read peak, where the library
/// declines monotonically from position 21.
///
/// The last point matters more than its single entry suggests. Quality falls from Q32.4 at 241
/// to **Q24.2 at 250** -- a sharp final-cycle drop sitting inside the 50-base tail window, and
/// the reason healthy reads contribute as many low bases to that window as they do.
const PROFILE: &[(usize, f64)] = &[
    (1, 34.5),
    (21, 39.3),
    (41, 39.1),
    (61, 39.0),
    (81, 38.7),
    (101, 38.5),
    (121, 38.4),
    (141, 38.0),
    (161, 37.5),
    (181, 36.7),
    (201, 35.6),
    (221, 34.1),
    (241, 32.4),
    (250, 24.2),
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
/// 15,000 reads gives a binomial standard error of about 0.24 percentage points on a ~9.8%
/// rate, so the +/-1.5 point assertion bands below sit roughly six sigma out and will not flake.
/// Larger costs test time for no added confidence.
const N_FIXTURE_READS: usize = 15_000;

/// UNMEASURED, and flagged as such. The per-base spread around the profile (2.5) and around a
/// dip's depth (4.0) are not in any measurement taken so far -- the script reports means, not
/// variances. They are plausible rather than known, and are the first thing to revisit if the
/// fixture's shape is ever questioned.
///
/// One known divergence from the library, left in deliberately. Real HG002 R1 has 31 DISCRETE
/// quality levels -- [2, 10, 11, ..., 40], with no Q3-Q9 and no Q17 (recorded in #694) -- and
/// this fixture emits a dense axis clamped to [2, 40]. It therefore produces scores the
/// instrument never emits. That is immaterial to the tail-collapse statistics this fixture is
/// calibrated on, but it WILL matter to a mixture model fitted from it, since the fitted
/// `quality_score_options` would not match a real one. #694 makes the same point about NovaSeq
/// emitting four bins.
const BASELINE_SD: f64 = 2.5;
const DIP_DEPTH_SD: f64 = 4.0;

/// MEASURED, not tuned. Mean depth of a low run across the library.
const DIP_DEPTH: f64 = 15.5;
/// Geometric continuation probability for a low run. Two values, because one cannot fit the
/// measured shape: 80.5% of real runs are 1-2 bases, yet the distribution has a heavy tail
/// (max 204) and collapsed reads carry a longest run averaging 20.47 bases. A single geometric
/// giving the right mean puts far too much weight at 3-5 and none at 26+. Runs lengthening
/// where the read is already failing produces both.
const DIP_CONTINUE: f64 = 0.42;
/// ...rising to this at the 3' end of a degraded read.
const DIP_CONTINUE_TAIL: f64 = 0.92;

/// Per-base probability that a HEALTHY read starts a dip. Derived from the measured 0.398% of
/// low bases in a healthy read's first 100 positions via f = rL/(1+rL), and VERIFIED against
/// the fixture: it reproduces 0.386%.
const HEALTHY_DIP_RATE: f64 = 0.00183;
/// The same for a DEGRADED read at the start. TUNED, not derived, and the distinction matters.
///
/// The same arithmetic applied to the measured 8.186% gives 0.0407, and that value produces
/// 11.357% in the fixture -- 39% over target. The derivation is wrong because the 8.186% is
/// measured on reads SELECTED for collapsing, which are the high-severity tail of the degraded
/// population, not a fair sample of it. 0.026 is the rate that makes the selected subset match.
/// This value was a placeholder carrying a derivation that does not hold; it is now fitted to
/// the statistic it is supposed to reproduce (fixture: 8.617%).
const DEGRADED_DIP_HEAD: f64 = 0.026;
/// ...ramping to this by the 3' end. Tuned, with the two tail thresholds as the constraint.
const DEGRADED_DIP_TAIL: f64 = 0.36;
/// Spread of per-read severity within the degraded population, so some cross Q20 and some
/// only cross Q25. Without it every degraded read collapses to the same depth and the
/// 2.20 ratio between the thresholds cannot be reproduced.
const SEVERITY_SPREAD: f64 = 0.50;
const DEGRADED_FRACTION: f64 = 0.125;

fn encode(q: f64) -> u8 {
    (q.round().clamp(2.0, 40.0) as u8) + 33
}

/// Where to dump the generated fixture for inspection, if set. The fixture is calibrated
/// against several statistics now, not one, so being able to run the SAME script that measured
/// the real library over it is what keeps the two comparable:
///
///   QTC_DUMP_FIXTURE=/tmp/fx.fastq.gz cargo test --test quality_tail_collapse the_fixture
///   HEAD_WINDOW=100 FASTQ=/tmp/fx.fastq.gz MAX_READS=0 OUT=/tmp/fx.txt \
///       bash scripts/delta/measure_quality_degradation.sh
fn dump_path() -> Option<PathBuf> {
    std::env::var_os("QTC_DUMP_FIXTURE").map(PathBuf::from)
}

/// Write a FASTQ whose degraded reads are noisier THROUGHOUT and worst at the 3' end, and whose
/// healthy reads dip rarely and briefly.
fn write_fixture(path: &Path, n_reads: usize, seed: u64) {
    write_fixture_at(path, n_reads, seed, DEGRADED_FRACTION)
}

/// As `write_fixture`, with the degraded fraction chosen by the caller. The two mates of a
/// real library carry different fractions -- HG002 measures R2 at 2.34x R1 on collapsed-tail
/// rate -- and a four-population test needs them to differ or it cannot tell the mates apart.
fn write_fixture_at(path: &Path, n_reads: usize, seed: u64, degraded_fraction: f64) {
    let mut rng = Lcg::new(seed);
    let file = std::fs::File::create(path).unwrap();
    let mut gz = GzEncoder::new(file, Compression::default());
    for i in 0..n_reads {
        let degraded = rng.unit() < degraded_fraction;
        // Per-read severity. A degraded read is not "collapsed" or not -- reads differ in how
        // noisy they are, which is what the head-window measurement showed (#694): reads that
        // end badly carry 20x the low-base rate of a healthy read a hundred positions earlier.
        let severity = if degraded {
            (1.0 + rng.gauss(0.0, SEVERITY_SPREAD)).max(0.05)
        } else {
            1.0
        };
        let mut quals = Vec::with_capacity(READ_LEN);
        let mut in_dip = 0usize;
        for pos in 1..=READ_LEN {
            let base = healthy_mean(pos);
            // The dip rate ramps along the read rather than switching at an onset. A degraded
            // read is noisier everywhere and worst at the 3' end; a healthy one is flat.
            let t = (pos - 1) as f64 / (READ_LEN - 1) as f64;
            let dip_rate = if degraded {
                severity * (DEGRADED_DIP_HEAD + (DEGRADED_DIP_TAIL - DEGRADED_DIP_HEAD) * t * t)
            } else {
                HEALTHY_DIP_RATE
            };
            let q = if in_dip > 0 {
                in_dip -= 1;
                rng.gauss(DIP_DEPTH, DIP_DEPTH_SD)
            } else if rng.unit() < dip_rate {
                // Geometric run length: mean 2.19 bases, most runs 1-2, a long tail.
                let cont = if degraded {
                    DIP_CONTINUE + (DIP_CONTINUE_TAIL - DIP_CONTINUE) * t * t
                } else {
                    DIP_CONTINUE
                };
                in_dip = 0;
                while rng.unit() < cont && in_dip < 120 {
                    in_dip += 1;
                }
                rng.gauss(DIP_DEPTH, DIP_DEPTH_SD)
            } else {
                rng.gauss(base, BASELINE_SD)
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

/// Head-of-read statistics split by whether the read ends with a collapsed tail: returns
/// (collapsed mean Q, healthy mean Q, collapsed low-base %, healthy low-base %).
///
/// This is what distinguishes a PROPENSITY fixture from an ONSET one, and the tail rates alone
/// cannot: the previous onset-shaped fixture hit the same two thresholds while its collapsed
/// reads were indistinguishable from healthy ones at position 1-100.
fn head_stats(path: &Path, head_win: usize, tail_win: usize) -> (f64, f64, f64, f64) {
    let reader = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
    let (mut cq, mut hq, mut cn, mut hn) = (0.0f64, 0.0f64, 0usize, 0usize);
    let (mut clow, mut hlow) = (0usize, 0usize);
    for (i, line) in reader.lines().enumerate() {
        if i % 4 != 3 {
            continue;
        }
        let line = line.unwrap();
        if line.len() < tail_win.max(head_win) {
            continue;
        }
        let q: Vec<f64> = line.as_bytes().iter().map(|&b| (b - 33) as f64).collect();
        let tail: f64 = q[q.len() - tail_win..].iter().sum::<f64>() / tail_win as f64;
        let head = &q[..head_win];
        let mean = head.iter().sum::<f64>() / head_win as f64;
        let low = head.iter().filter(|&&x| x <= 25.0).count();
        if tail < 25.0 {
            cq += mean;
            clow += low;
            cn += 1;
        } else {
            hq += mean;
            hlow += low;
            hn += 1;
        }
    }
    assert!(cn > 0 && hn > 0, "both populations must be present");
    (
        cq / cn as f64,
        hq / hn as f64,
        100.0 * clow as f64 / (cn * head_win) as f64,
        100.0 * hlow as f64 / (hn * head_win) as f64,
    )
}

/// Fit a model on the fixture, generate reads with it, and return the generated R1 path.
fn fit_and_generate(work: &Path, fixture: &Path) -> PathBuf {
    fit_and_generate_with(work, fixture, false)
}

/// `degradation` fits the two-population model (#694) rather than one.
fn fit_and_generate_with(work: &Path, fixture: &Path, degradation: bool) -> PathBuf {
    let model = work.join("model.json.gz");
    let cfg = work.join("fit.yml");
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\nfit_quality_degradation: {degradation}\n",
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

/// Fit BOTH mates into one model with the degraded population on -- four tensors, R1/R2 x
/// healthy/degraded -- then generate PAIRED reads and return both output FASTQs.
fn fit_pair_and_generate(work: &Path, r1: &Path, r2: &Path) -> (PathBuf, PathBuf) {
    let model = work.join("model_pair.json.gz");
    let cfg = work.join("fit_pair.yml");
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\nfastq_file_r2: {}\noutput_file: {}\noverwrite_output: true\n\
             max_reads: 0\nqual_offset: 33\nfit_quality_degradation: true\n",
            r1.display(),
            r2.display(),
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
    let mut reads_cfg = GenReadsConfig::new(reference, work.to_path_buf(), "genpair");
    reads_cfg.read_len = READ_LEN;
    reads_cfg.coverage = 2;
    reads_cfg.paired_ended = true;
    // A fragment must hold both mates. The helper's 200 bp default is shorter than one
    // 250 bp read, which would make this a test of adapter readthrough instead.
    reads_cfg.fragment_mean = Some(600.0);
    reads_cfg.fragment_st_dev = Some(60.0);
    reads_cfg.produce_fastq = true;
    reads_cfg.rng_seed = "q694pair".to_string();
    reads_cfg.sequence_error_model = Some(model);
    let gen_cfg = reads_cfg.write_yaml();
    eidolon()
        .args(["gen-reads", "-c", gen_cfg.path().to_str().unwrap()])
        .assert()
        .success();

    (
        work.join("genpair_r1.fastq.gz"),
        work.join("genpair_r2.fastq.gz"),
    )
}

/// The fixture itself must match what was measured on HG002, or nothing downstream means
/// anything. Both rates, because one does not constrain the shape: an earlier fixture hit
/// Q<25 while putting nearly every collapsed read below Q20, where real data splits 2.18 to 1.
#[test]
fn the_fixture_reproduces_the_measured_degradation_rates() {
    let (_g, work) = fresh_workdir();
    let fixture = dump_path().unwrap_or_else(|| work.join("fixture.fastq.gz"));
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

/// The fixture must implement a PROPENSITY, not an onset: reads that end with a collapsed tail
/// must already be measurably noisier at position 1-100, a hundred positions before the decline
/// appears in the mean profile.
///
/// Measured on HG002 R1 (job 22093395): collapsed reads run Q35.58 against Q39.11 for healthy
/// ones, carrying 8.186% low bases against 0.398% -- a 20.6x ratio. An onset model predicts no
/// difference at all, and the fixture this file used to carry reported -0.08 Q and 0.88x.
///
/// Without this test the tail rates alone would accept either shape, because they did.
#[test]
fn the_fixture_degraded_population_is_a_propensity_not_an_onset() {
    let (_g, work) = fresh_workdir();
    let fixture = work.join("fixture.fastq.gz");
    write_fixture(&fixture, N_FIXTURE_READS, 694);

    let (c_mean, h_mean, c_low, h_low) = head_stats(&fixture, 100, 50);
    let ratio = c_low / h_low;
    assert!(
        c_mean < h_mean - 1.5,
        "collapsed reads must START worse: collapsed Q{c_mean:.2}, healthy Q{h_mean:.2} \
         (real HG002 differs by -3.53). An onset fixture shows no difference here."
    );
    assert!(
        (8.0..45.0).contains(&ratio),
        "collapsed reads carry {c_low:.3}% low bases in their first 100 positions against \
         {h_low:.3}% for healthy, a {ratio:.1}x ratio; HG002 measures 20.6x. Near 1x means the \
         fixture has reverted to an onset shape."
    );
}

/// THE TARGET for #694, end to end: fit the fixture with the two-population model and the
/// generated reads must carry the degraded population the fixture planted.
///
/// This is the whole chain -- fixture, fitter, model file, generation -- rather than any one
/// layer. The single-population fit of the same fixture is the characterization test below, and
/// it reads about 1.8% against a planted 9.8%.
#[test]
fn a_two_population_fit_reproduces_the_degraded_population() {
    let (_g, work) = fresh_workdir();
    let fixture = work.join("fixture.fastq.gz");
    write_fixture(&fixture, N_FIXTURE_READS, 694);
    let (in25, in20, _) = tail_collapse_rates(&fixture, 50);

    let generated = fit_and_generate_with(&work, &fixture, true);
    let (out25, out20, n_out) = tail_collapse_rates(&generated, 50);

    assert!(
        n_out > 5_000,
        "only {n_out} reads generated; too few to rate"
    );
    eprintln!(
        "TWO-POPULATION  fixture Q<25 {in25:.2}% Q<20 {in20:.2}%  |  generated Q<25 {out25:.2}% Q<20 {out20:.2}%  n={n_out}"
    );
    assert!(
        out25 > in25 / 2.0 && out25 < in25 * 2.0,
        "tail Q<25: fixture {in25:.2}%, generated {out25:.2}%. Against a single-population fit \
         at about 1.8%, landing between half and double the planted rate is the target."
    );
    assert!(
        out20 > in20 / 2.0 && out20 < in20 * 2.0,
        "tail Q<20: fixture {in20:.2}%, generated {out20:.2}%"
    );

    // The rate alone is NOT sufficient, and a mutation proved it: routing the degraded reads
    // into the healthy tensor leaves the degraded one empty, which falls back to uniform over
    // the option set. Uniform averages about Q21, so those reads collapse too and the rates
    // above still pass. It is garbage emitted at the right frequency.
    //
    // So assert the SHAPE. A fitted degraded population looks like the reads it was fitted
    // from: noisy at the head but nowhere near uniform. The fixture's collapsed reads run about
    // Q35 over their first 100 positions; uniform would be about Q21.
    let (c_mean, h_mean, c_low, h_low) = head_stats(&generated, 100, 50);
    eprintln!(
        "TWO-POPULATION  head: collapsed Q{c_mean:.2} ({c_low:.2}% low), healthy Q{h_mean:.2} ({h_low:.2}% low)"
    );
    assert!(
        c_mean > 30.0,
        "generated degraded reads average Q{c_mean:.2} over their first 100 positions. The \
         fixture's run about Q35; a uniform draw over the option set runs about Q21, which is \
         what an EMPTY degraded tensor falls back to. This reads like fallback, not a fit."
    );
    assert!(
        c_mean < h_mean,
        "degraded reads must still be worse than healthy ones at the head: collapsed \
         Q{c_mean:.2} against healthy Q{h_mean:.2}"
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
    //   Q<25  9.84% -> 1.80%   (5.5x short; real HG002 is ~27x)
    //   Q<20  4.37% -> 0.07%   (62x short; real HG002 reaches 0.00%)
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

/// FOUR POPULATIONS: R1/R2 x healthy/degraded, which #723 and #694 produce together and
/// which neither of them measured. #723 shipped a tested REFUSAL on this path and said
/// plainly that a loaded fit's fidelity was untested; this is that test.
///
/// KNOWN ANSWER, computable from the fixtures. Two libraries are planted with deliberately
/// different degraded fractions -- R2 at 2.34x R1, the ratio HG002 measures between its mates
/// -- and each generated mate must carry its OWN planted rate.
///
/// It is FALSIFIED three ways, and the last is the reason the test is shaped around the pair:
///
/// * both mates land on one rate -> one quality model is still serving both
/// * the mates' rates are swapped -> the mates are wired backwards
/// * both mates lose the degraded population -> the degradation split did not survive being
///   fitted per mate, which is exactly the interaction neither feature tested alone
#[test]
fn a_four_population_fit_keeps_the_mates_and_the_populations_apart() {
    let (_g, work) = fresh_workdir();
    let r1_fixture = work.join("fixture_r1.fastq.gz");
    let r2_fixture = work.join("fixture_r2.fastq.gz");
    // HG002 measures R2 at 2.34x R1 on collapsed-tail rate.
    write_fixture_at(&r1_fixture, N_FIXTURE_READS, 694, DEGRADED_FRACTION);
    write_fixture_at(&r2_fixture, N_FIXTURE_READS, 723, DEGRADED_FRACTION * 2.34);

    let (in1_25, in1_20, _) = tail_collapse_rates(&r1_fixture, 50);
    let (in2_25, in2_20, _) = tail_collapse_rates(&r2_fixture, 50);
    assert!(
        in2_25 > in1_25 * 1.5,
        "the fixtures must differ before anything downstream means anything: R1 {in1_25:.2}%, \
         R2 {in2_25:.2}%"
    );

    let (gen1, gen2) = fit_pair_and_generate(&work, &r1_fixture, &r2_fixture);
    let (out1_25, out1_20, n1) = tail_collapse_rates(&gen1, 50);
    let (out2_25, out2_20, n2) = tail_collapse_rates(&gen2, 50);

    eprintln!(
        "FOUR-POPULATION  R1 fixture Q<25 {in1_25:.2}% Q<20 {in1_20:.2}% -> generated \
         {out1_25:.2}% / {out1_20:.2}% (n={n1})"
    );
    eprintln!(
        "FOUR-POPULATION  R2 fixture Q<25 {in2_25:.2}% Q<20 {in2_20:.2}% -> generated \
         {out2_25:.2}% / {out2_20:.2}% (n={n2})"
    );
    assert!(
        n1 > 5_000 && n2 > 5_000,
        "too few reads to rate: {n1} / {n2}"
    );

    // Each mate carries its own planted rate, to the same factor-of-two bar the
    // two-population test uses.
    assert!(
        out1_25 > in1_25 / 2.0 && out1_25 < in1_25 * 2.0,
        "R1 tail Q<25: fixture {in1_25:.2}%, generated {out1_25:.2}%"
    );
    assert!(
        out2_25 > in2_25 / 2.0 && out2_25 < in2_25 * 2.0,
        "R2 tail Q<25: fixture {in2_25:.2}%, generated {out2_25:.2}%. Landing near R1's \
         {out1_25:.2}% instead means one quality model is serving both mates."
    );

    // THE DECISION ASSERTION. The mates must stay ordered. One model for both, or the mates
    // wired backwards, both pass every aggregate measure over the pair.
    assert!(
        out2_25 > out1_25 * 1.5,
        "R2 was planted 2.34x worse than R1 and must generate worse: R1 {out1_25:.2}%, R2 \
         {out2_25:.2}%. Equal rates mean one model serves both; reversed means the mates are \
         wired backwards."
    );

    // And the degraded population must be a FIT, not the uniform fallback an empty tensor
    // produces -- the mutation that caught this on the two-population test applies to each of
    // the four tensors independently.
    for (label, path) in [("R1", &gen1), ("R2", &gen2)] {
        let (c_mean, h_mean, c_low, h_low) = head_stats(path, 100, 50);
        eprintln!(
            "FOUR-POPULATION  {label} head: collapsed Q{c_mean:.2} ({c_low:.2}% low), healthy \
             Q{h_mean:.2} ({h_low:.2}% low)"
        );
        assert!(
            c_mean > 30.0,
            "{label} degraded reads average Q{c_mean:.2} over their first 100 positions. A \
             uniform draw over the option set runs about Q21, which is what an EMPTY degraded \
             tensor falls back to. This reads like fallback, not a fit."
        );
        assert!(
            c_mean < h_mean,
            "{label} degraded reads must still be worse at the head: collapsed Q{c_mean:.2}, \
             healthy Q{h_mean:.2}"
        );
    }
}
