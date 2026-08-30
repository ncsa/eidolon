use crate::gen_frag_length_model::{
    errors::GenFragLengthModelError, utils::config::RunConfiguration,
};
use eidolon_core::{
    file_tools::bam_reader::read_fragment_lengths, models::fragment_length::FragmentLengthModel,
};
use log::info;
use std::collections::HashMap;

// Only consider fragment lengths up to this many median absolute deviations above the median.
// Mirrors the FILTER_MEDDEV_M constant from Python NEAT.
const FILTER_MEDDEV_M: f64 = 10.0;

pub fn runner(config: &RunConfiguration) -> Result<(), GenFragLengthModelError> {
    info!("Reading fragment lengths from {:?}", config.input_file);
    let tlens = read_fragment_lengths(&config.input_file)?;
    run_from_tlens(
        tlens,
        config.min_reads,
        &config.output_file,
        config.distribution,
    )
}

/// Builds and writes a fragment length model from a pre-collected list of
/// template lengths (e.g. produced by `FragLengthObserver` during a shared
/// BAM walk in the unified `gen-bam-models` runner). Applies MAD-based
/// outlier filtering and rare-length pruning, fits a normal distribution,
/// and writes the gzipped JSON model.
/// Which distribution family the built model uses.
///
/// `Discrete` keeps the observed shape. `Normal` collapses it to mean + st_dev, which is
/// what this builder did before v3.3.0 and is kept for sparse inputs (exome, amplicon, a
/// small targeted BAM) where there is not enough data to estimate a shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DistributionKind {
    #[default]
    Discrete,
    Normal,
}

impl DistributionKind {
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "discrete" | "empirical" => Some(Self::Discrete),
            "normal" | "gaussian" => Some(Self::Normal),
            _ => None,
        }
    }
}

/// Widen the smoothing bandwidth by this much per attempt when gaps survive.
const BANDWIDTH_GROWTH: f64 = 1.6;
/// Give up after this many widenings rather than smoothing the shape into a flat line.
const MAX_SMOOTHING_PASSES: usize = 24;
/// Kernel truncation, in bandwidths. Beyond 4 sigma the Gaussian contributes < 1e-4.
const KERNEL_TRUNCATION: f64 = 4.0;
/// Trim each end while the cumulative mass there is below this.
///
/// MASS, NOT A QUANTILE. A quantile trim collapses onto the mode when the distribution is
/// concentrated -- q0.999 of "3000 reads at 400, one at 401" is 400, so the trim deletes a
/// length adjacent to the mode. What actually wants removing is the isolated stray: the
/// model this crate shipped before v3.3.0 has a bin at length 1 carrying 2.3e-9 of the
/// mass, with a 30-wide hole above it. A mass floor removes that and keeps anything with
/// real weight, whatever quantile it happens to land on.
const TRIM_MASS: f64 = 1e-5;
/// Below this many distinct observed lengths there is no shape to estimate, and smoothing
/// would only smear a handful of spikes into a plausible-looking line.
const MIN_DISTINCT_LENGTHS: usize = 10;

/// Builds and writes a fragment length model from a pre-collected list of template lengths
/// (e.g. produced by `FragLengthObserver` during a shared BAM walk in the unified
/// `gen-bam-models` runner).
///
/// EVERYTHING HERE IS COMPUTED FROM THE HISTOGRAM. The previous implementation sorted a
/// per-read Vec<usize>, built a second full-size Vec of absolute deviations, sorted that
/// too, and then rebuilt a third full-size Vec with repeat_n -- three simultaneous copies
/// and two O(n log n) sorts over one element per read pair. On a whole-genome BAM that is
/// ~800M elements, so ~19 GB before a model is fitted. Median, MAD, mean and standard
/// deviation are all functions of the counts, so the histogram is the only thing that ever
/// needs to exist: O(distinct lengths), about 800 entries.
pub fn run_from_tlens(
    tlens: Vec<usize>,
    min_reads: usize,
    output_file: &std::path::PathBuf,
    kind: DistributionKind,
) -> Result<(), GenFragLengthModelError> {
    if tlens.is_empty() {
        return Err(GenFragLengthModelError::EmptyData);
    }
    info!("Collected {} raw fragment lengths", tlens.len());

    let hist = histogram(&tlens);
    drop(tlens);

    let trimmed = trim(hist, min_reads);
    if trimmed.is_empty() {
        return Err(GenFragLengthModelError::FilteredToEmpty);
    }
    let kept: u64 = trimmed.iter().map(|&(_, c)| c).sum();
    info!(
        "Retained {} fragment lengths across {} distinct values after filtering",
        kept,
        trimmed.len()
    );

    let mean = hist_mean(&trimmed);
    let st_dev = hist_std_dev(&trimmed, mean);

    let model = match kind {
        DistributionKind::Normal => {
            info!("Fragment length model: Normal(mean={mean:.1}, std_dev={st_dev:.1})");
            FragmentLengthModel::new_normal(mean, st_dev)?
        }
        DistributionKind::Discrete => {
            // A discrete model with holes in it is the failure this guards against: it
            // builds clean, serializes clean, and then the shape it hands gen-reads is not
            // the shape that was measured. Smoothing fills them from the neighbourhood
            // rather than inventing a floor, and the result is ASSERTED gap-free below.
            let (values, weights) = smooth_to_gap_free(&trimmed, st_dev, kept)?;
            info!(
                "Fragment length model: Discrete over {}-{} ({} bins, no gaps), mean={:.1}, std_dev={:.1}",
                values[0],
                values[values.len() - 1],
                values.len(),
                mean,
                st_dev
            );
            FragmentLengthModel::new_discrete(values, weights)?
        }
    };

    model.write_file(output_file)?;
    info!("Wrote fragment length model to {output_file:?}");
    Ok(())
}

/// Counts per observed length, ascending. The only full pass over the input.
fn histogram(tlens: &[usize]) -> Vec<(usize, u64)> {
    let mut counts: HashMap<usize, u64> = HashMap::new();
    for &l in tlens {
        *counts.entry(l).or_default() += 1;
    }
    let mut out: Vec<(usize, u64)> = counts.into_iter().collect();
    out.sort_unstable_by_key(|&(l, _)| l);
    out
}

/// Drops zero-length entries, anything above the MAD outlier ceiling, and the extreme
/// quantile tails.
///
/// `min_reads` no longer deletes bins. Deleting a sparse bin is what MANUFACTURES a gap,
/// and the tail -- the part we most want to keep -- is exactly where counts are lowest. It
/// survives as a floor on the total number of observations instead.
fn trim(hist: Vec<(usize, u64)>, min_reads: usize) -> Vec<(usize, u64)> {
    let hist: Vec<(usize, u64)> = hist.into_iter().filter(|&(l, _)| l > 0).collect();
    if hist.is_empty() {
        return hist;
    }
    let total: u64 = hist.iter().map(|&(_, c)| c).sum();
    if (total as usize) < min_reads.max(1) {
        return Vec::new();
    }
    if min_reads == 0 {
        return hist;
    }

    let median = hist_quantile(&hist, 0.5) as f64;
    let mad = hist_mad(&hist, median);
    // A concentrated histogram has MAD 0 -- "3000 reads at 400, one at 401" gives a median
    // absolute deviation of exactly zero -- and `median + 10 * 0` is the median, so the
    // ceiling would delete every length above the mode. There is no outlier scale to
    // measure here, so there is no outlier rule to apply.
    let ceiling = if mad > 0.0 {
        median + FILTER_MEDDEV_M * mad
    } else {
        f64::INFINITY
    };

    let cutoff = total as f64 * TRIM_MASS;
    let mut lo = hist[0].0;
    let mut acc = 0f64;
    for &(l, c) in &hist {
        if acc + c as f64 > cutoff {
            lo = l;
            break;
        }
        acc += c as f64;
    }
    let mut hi = hist[hist.len() - 1].0;
    acc = 0f64;
    for &(l, c) in hist.iter().rev() {
        if acc + c as f64 > cutoff {
            hi = l;
            break;
        }
        acc += c as f64;
    }

    hist.into_iter()
        .filter(|&(l, _)| l >= lo && l <= hi && (l as f64) <= ceiling)
        .collect()
}

/// Densify the support and smooth until no interior bin is empty.
///
/// Bandwidth comes from Silverman's rule, which scales as n^(-1/5): sparse inputs get more
/// smoothing, which is exactly when gaps appear, and a dense whole-genome histogram gets a
/// bandwidth of a couple of bp, so the result is the empirical distribution with its holes
/// bridged rather than a reshaped one. If gaps still survive the bandwidth widens and it
/// tries again -- so "no gaps" is a checked postcondition, not a hope.
fn smooth_to_gap_free(
    hist: &[(usize, u64)],
    st_dev: f64,
    n: u64,
) -> Result<(Vec<usize>, Vec<f64>), GenFragLengthModelError> {
    let lo = hist[0].0;
    let hi = hist[hist.len() - 1].0;
    // Refuse before smoothing rather than after. Given enough bandwidth the kernel will
    // bridge ANY two spikes, and the result builds clean, serializes clean, and hands
    // gen-reads a flat smear -- a model that is broken in exactly the way that is hardest
    // to notice. Too few distinct observations is the honest answer.
    if hist.len() < MIN_DISTINCT_LENGTHS {
        return Err(GenFragLengthModelError::ConfigurationError(format!(
            "only {} distinct fragment lengths observed ({lo}-{hi}); a discrete model needs \
             at least {MIN_DISTINCT_LENGTHS} to have a shape worth estimating. Use \
             `distribution: normal` for sparse input, or supply a BAM with more pairs.",
            hist.len()
        )));
    }
    debug_assert!(
        hi > lo,
        "{MIN_DISTINCT_LENGTHS} distinct values cannot share one length"
    );
    let span = hi - lo + 1;

    let mut dense = vec![0f64; span];
    for &(l, c) in hist {
        dense[l - lo] = c as f64;
    }
    // Already contiguous: nothing to bridge, so do not touch the measured shape at all.
    if !dense.iter().any(|&w| w <= 0.0) {
        return Ok(finish(lo, dense));
    }

    let iqr = (hist_quantile(hist, 0.75) as f64) - (hist_quantile(hist, 0.25) as f64);
    let spread = if iqr > 0.0 {
        st_dev.min(iqr / 1.34)
    } else {
        st_dev
    };
    let mut h = (0.9 * spread * (n as f64).powf(-0.2)).max(0.5);

    for _ in 0..MAX_SMOOTHING_PASSES {
        let smoothed = gaussian_smooth(&dense, h);
        if !smoothed.iter().any(|&w| w <= 0.0) {
            return Ok(finish(lo, smoothed));
        }
        h *= BANDWIDTH_GROWTH;
    }
    Err(GenFragLengthModelError::ConfigurationError(format!(
        "fragment lengths {lo}-{hi} could not be smoothed into a gap-free distribution; \
         the input is too sparse for a discrete model. Use `distribution: normal`, or \
         supply a BAM with more pairs."
    )))
}

fn gaussian_smooth(dense: &[f64], h: f64) -> Vec<f64> {
    let radius = (KERNEL_TRUNCATION * h).ceil() as usize;
    let two_h2 = 2.0 * h * h;
    let kernel: Vec<f64> = (0..=radius)
        .map(|d| (-((d * d) as f64) / two_h2).exp())
        .collect();

    let mut out = vec![0f64; dense.len()];
    for (i, &c) in dense.iter().enumerate() {
        if c <= 0.0 {
            continue;
        }
        let from = i.saturating_sub(radius);
        let to = (i + radius).min(dense.len() - 1);
        for (j, slot) in out.iter_mut().enumerate().take(to + 1).skip(from) {
            *slot += c * kernel[i.abs_diff(j)];
        }
    }
    out
}

fn finish(lo: usize, weights: Vec<f64>) -> (Vec<usize>, Vec<f64>) {
    let values = (lo..lo + weights.len()).collect();
    (values, weights)
}

/// Smallest length at or below which `q` of the mass sits.
fn hist_quantile(hist: &[(usize, u64)], q: f64) -> usize {
    let total: u64 = hist.iter().map(|&(_, c)| c).sum();
    if total == 0 {
        return 0;
    }
    let target = (total as f64 * q).ceil().max(1.0) as u64;
    let mut acc = 0u64;
    for &(l, c) in hist {
        acc += c;
        if acc >= target {
            return l;
        }
    }
    hist[hist.len() - 1].0
}

/// median(|x - median|), computed from counts rather than a second full-size vector.
fn hist_mad(hist: &[(usize, u64)], median: f64) -> f64 {
    let mut devs: Vec<(usize, u64)> = hist
        .iter()
        .map(|&(l, c)| ((l as f64 - median).abs().round() as usize, c))
        .collect();
    devs.sort_unstable_by_key(|&(d, _)| d);
    // Merge equal deviations so the quantile walk sees one entry per distinct value.
    let mut merged: Vec<(usize, u64)> = Vec::with_capacity(devs.len());
    for (d, c) in devs {
        match merged.last_mut() {
            Some((pd, pc)) if *pd == d => *pc += c,
            _ => merged.push((d, c)),
        }
    }
    hist_quantile(&merged, 0.5) as f64
}

fn hist_mean(hist: &[(usize, u64)]) -> f64 {
    let total: u64 = hist.iter().map(|&(_, c)| c).sum();
    if total == 0 {
        return 0.0;
    }
    let sum: f64 = hist.iter().map(|&(l, c)| l as f64 * c as f64).sum();
    sum / total as f64
}

fn hist_std_dev(hist: &[(usize, u64)], mean: f64) -> f64 {
    let total: u64 = hist.iter().map(|&(_, c)| c).sum();
    if total == 0 {
        return 0.0;
    }
    let var: f64 = hist
        .iter()
        .map(|&(l, c)| {
            let d = l as f64 - mean;
            d * d * c as f64
        })
        .sum::<f64>()
        / total as f64;
    var.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::models::fragment_length::FragmentLengthModel;

    // ── naive references, so the histogram math is checked against something that is
    //    obviously right rather than against itself ────────────────────────────────
    fn naive_mean(d: &[usize]) -> f64 {
        d.iter().map(|&x| x as f64).sum::<f64>() / d.len() as f64
    }
    fn naive_std(d: &[usize], mean: f64) -> f64 {
        (d.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / d.len() as f64).sqrt()
    }
    fn expand(h: &[(usize, u64)]) -> Vec<usize> {
        h.iter()
            .flat_map(|&(l, c)| std::iter::repeat_n(l, c as usize))
            .collect()
    }
    /// Third standardised moment. Positive = right tail, which is the property the whole
    /// change exists to preserve.
    fn skew(d: &[usize]) -> f64 {
        let m = naive_mean(d);
        let s = naive_std(d, m);
        if s == 0.0 {
            return 0.0;
        }
        d.iter().map(|&x| ((x as f64 - m) / s).powi(3)).sum::<f64>() / d.len() as f64
    }
    fn weighted_skew(values: &[usize], weights: &[f64]) -> f64 {
        let tot: f64 = weights.iter().sum();
        let m: f64 = values
            .iter()
            .zip(weights)
            .map(|(&v, &w)| v as f64 * w)
            .sum::<f64>()
            / tot;
        let var: f64 = values
            .iter()
            .zip(weights)
            .map(|(&v, &w)| (v as f64 - m).powi(2) * w)
            .sum::<f64>()
            / tot;
        let sd = var.sqrt();
        values
            .iter()
            .zip(weights)
            .map(|(&v, &w)| ((v as f64 - m) / sd).powi(3) * w)
            .sum::<f64>()
            / tot
    }
    /// Deterministic right-skewed sample: a long thin tail above a dense mode.
    fn right_skewed() -> Vec<usize> {
        let mut v = Vec::new();
        for l in 300..=420 {
            v.extend(std::iter::repeat_n(l, 400));
        }
        for (i, l) in (421..=900).enumerate() {
            let c = (300usize).saturating_sub(i * 2).max(1);
            v.extend(std::iter::repeat_n(l, c));
        }
        v
    }

    // ── the histogram math reproduces the per-read math exactly ──────────────

    #[test]
    fn histogram_statistics_match_the_naive_per_read_computation() {
        // The refactor's whole premise: mean, sd, median and MAD are functions of the
        // counts, so dropping the three full-size vectors must change nothing.
        let data = right_skewed();
        let h = histogram(&data);
        let m = naive_mean(&data);
        assert!((hist_mean(&h) - m).abs() < 1e-9, "mean drifted");
        assert!(
            (hist_std_dev(&h, m) - naive_std(&data, m)).abs() < 1e-9,
            "sd drifted"
        );

        let mut sorted = data.clone();
        sorted.sort_unstable();
        assert_eq!(
            hist_quantile(&h, 0.5),
            sorted[(sorted.len() as f64 * 0.5).ceil().max(1.0) as usize - 1],
            "median drifted"
        );
    }

    #[test]
    fn test_compute_mean_and_std_dev_known_answer() {
        // Population std of [100, 200, 300] = sqrt(20000/3) ~ 81.65, computed by hand.
        let h = histogram(&[100usize, 200, 300]);
        let mean = hist_mean(&h);
        assert!((mean - 200.0).abs() < 1e-10);
        assert!((hist_std_dev(&h, mean) - (20000.0f64 / 3.0).sqrt()).abs() < 1e-6);
    }

    // ── trim ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_trim_zero_min_reads_passthrough() {
        let data = vec![100usize, 200, 100, 300, 50];
        assert_eq!(expand(&trim(histogram(&data), 0)), {
            let mut d = data.clone();
            d.sort_unstable();
            d
        });
    }

    #[test]
    fn test_trim_removes_outliers() {
        let mut data: Vec<usize> = (180..=220).flat_map(|v| vec![v; 5]).collect();
        data.extend(std::iter::repeat_n(100_000usize, 5));
        let kept = expand(&trim(histogram(&data), 2));
        assert!(
            kept.iter().all(|&x| x < 100_000),
            "outlier should be filtered"
        );
        assert!(!kept.is_empty());
    }

    #[test]
    fn test_trim_keeps_rare_lengths_instead_of_deleting_them() {
        // CHANGED DELIBERATELY. `min_reads` used to delete any bin with fewer than
        // min_reads observations, which is what MANUFACTURED gaps -- and it deleted them
        // hardest in the tail, the part worth keeping. Sparse bins are now handled by
        // smoothing; min_reads survives only as a floor on the total.
        let mut data: Vec<usize> = std::iter::repeat_n(100usize, 3000).collect();
        data.push(101); // a single observation, adjacent to the mode
        let kept = trim(histogram(&data), 2);
        assert!(
            kept.iter().any(|&(l, _)| l == 101),
            "a rare in-range length must survive trimming, got {kept:?}"
        );
    }

    #[test]
    fn test_trim_empty_input() {
        assert!(trim(histogram(&[]), 2).is_empty());
    }

    #[test]
    fn test_trim_refuses_when_total_is_below_min_reads() {
        assert!(trim(histogram(&[400usize]), 50).is_empty());
    }

    // ── the gap problem ──────────────────────────────────────────────────────

    #[test]
    fn discrete_model_has_no_gaps_even_from_a_gappy_histogram() {
        // Every OTHER length missing: the shape of a sparse real histogram, and the shape
        // that produced a model that built fine and behaved badly.
        let data: Vec<usize> = (200..=400)
            .step_by(2)
            .flat_map(|l| std::iter::repeat_n(l, 20))
            .collect();
        let h = trim(histogram(&data), 2);
        let mean = hist_mean(&h);
        let sd = hist_std_dev(&h, mean);
        let n: u64 = h.iter().map(|&(_, c)| c).sum();
        let (values, weights) = smooth_to_gap_free(&h, sd, n).unwrap();

        assert_eq!(
            values.len(),
            values[values.len() - 1] - values[0] + 1,
            "support must be contiguous"
        );
        assert!(
            weights.iter().all(|&w| w > 0.0),
            "no interior bin may be empty"
        );
        for (i, pair) in values.windows(2).enumerate() {
            assert_eq!(pair[1], pair[0] + 1, "gap at index {i}");
        }
    }

    #[test]
    fn a_contiguous_histogram_is_not_smoothed_at_all() {
        // The must-not-fire case. Smoothing that always fires would pass the gap test
        // above while quietly reshaping every dense whole-genome model.
        let data = right_skewed();
        let h = trim(histogram(&data), 2);
        let mean = hist_mean(&h);
        let sd = hist_std_dev(&h, mean);
        let n: u64 = h.iter().map(|&(_, c)| c).sum();
        let (values, weights) = smooth_to_gap_free(&h, sd, n).unwrap();
        for (&v, &w) in values.iter().zip(&weights) {
            let observed = h.iter().find(|&&(l, _)| l == v).map(|&(_, c)| c as f64);
            assert_eq!(
                Some(w),
                observed,
                "length {v} was altered though the input had no gaps"
            );
        }
    }

    // ── the decision under test: does the SHAPE survive? ─────────────────────

    #[test]
    fn discrete_preserves_the_right_tail_and_normal_destroys_it() {
        let data = right_skewed();
        let input_skew = skew(&data);
        assert!(
            input_skew > 0.4,
            "fixture must actually be right-skewed, got {input_skew}"
        );

        let h = trim(histogram(&data), 2);
        // Two separate claims, because trimming removes outliers BY DESIGN and skew is a
        // cubed moment, so the extreme tail dominates it. Conflating them would let a
        // tail-destroying trim hide behind a shape-preserving smoother, or vice versa.
        let trimmed_skew = skew(&expand(&h));
        assert!(
            trimmed_skew > 0.85,
            "trimming must not flatten the tail: input {input_skew:.3}, after trim {trimmed_skew:.3}"
        );

        let mean = hist_mean(&h);
        let sd = hist_std_dev(&h, mean);
        let n: u64 = h.iter().map(|&(_, c)| c).sum();
        let (values, weights) = smooth_to_gap_free(&h, sd, n).unwrap();
        let kept = weighted_skew(&values, &weights);
        assert!(
            (kept - trimmed_skew).abs() < 0.01,
            "Discrete must reproduce the shape it was given: trimmed {trimmed_skew:.3}, model {kept:.3}"
        );

        // And the escape hatch must still do the old thing -- a Normal has skew 0 by
        // construction, which is precisely why ins_skew read 0.074 against a real 0.795.
        let normal = FragmentLengthModel::new_normal(mean, sd).unwrap();
        match normal {
            FragmentLengthModel::Normal { mean: m, st_dev } => {
                assert!((m - mean).abs() < 1e-9 && (st_dev - sd).abs() < 1e-9);
            }
            _ => panic!("new_normal must produce a Normal"),
        }
    }

    #[test]
    fn a_symmetric_input_stays_symmetric() {
        // The other must-not-fire: the change must not INVENT skew.
        let data: Vec<usize> = (300..=500)
            .flat_map(|l| {
                let d = (l as i64 - 400).unsigned_abs() as usize;
                std::iter::repeat_n(l, 200usize.saturating_sub(d).max(1))
            })
            .collect();
        let h = trim(histogram(&data), 2);
        let mean = hist_mean(&h);
        let sd = hist_std_dev(&h, mean);
        let n: u64 = h.iter().map(|&(_, c)| c).sum();
        let (values, weights) = smooth_to_gap_free(&h, sd, n).unwrap();
        assert!(
            weighted_skew(&values, &weights).abs() < 0.05,
            "symmetric input must not acquire a tail"
        );
    }

    #[test]
    fn an_unsmoothable_input_is_refused_rather_than_written() {
        // Two lengths, a thousand apart: no bandwidth bridges that without flattening the
        // distribution into noise. Refusing beats writing a model that builds fine.
        let data: Vec<usize> = std::iter::repeat_n(200usize, 5)
            .chain(std::iter::repeat_n(30_000usize, 5))
            .collect();
        let h = trim(histogram(&data), 0);
        let mean = hist_mean(&h);
        let sd = hist_std_dev(&h, mean);
        let n: u64 = h.iter().map(|&(_, c)| c).sum();
        assert!(
            smooth_to_gap_free(&h, sd, n).is_err(),
            "an unsmoothable histogram must be refused"
        );
    }

    #[test]
    fn distribution_kind_parses_both_spellings_and_rejects_junk() {
        assert_eq!(
            DistributionKind::parse("discrete"),
            Some(DistributionKind::Discrete)
        );
        assert_eq!(
            DistributionKind::parse(" NORMAL "),
            Some(DistributionKind::Normal)
        );
        assert_eq!(DistributionKind::default(), DistributionKind::Discrete);
        assert_eq!(DistributionKind::parse("lognormal"), None);
    }

    // ── runner integration ────────────────────────────────────────────────────

    #[test]
    fn test_runner_with_bam() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("frags.bam");
        write_test_frag_bam(
            &bam_path,
            &[150usize, 151, 152, 150, 151, 152, 150, 151, 152, 150],
        );
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output.clone(),
            overwrite_output: true,
            min_reads: 2,
            distribution: DistributionKind::Normal,
        };
        runner(&config).unwrap();
        assert!(output.exists());
        let model = FragmentLengthModel::discrete_from_file(&output).unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert!(mean > 100.0 && mean < 200.0, "mean={mean}");
                assert!(st_dev >= 0.0, "st_dev={st_dev}");
            }
            _ => panic!("Expected Normal model"),
        }
    }

    #[test]
    fn test_runner_min_reads_zero_skips_filter() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("frags.bam");
        // Every length appears exactly once → min_reads=2 would remove all; min_reads=0 keeps all
        write_test_frag_bam(
            &bam_path,
            &[100, 200, 300, 400, 500, 150, 250, 350, 450, 160],
        );
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output.clone(),
            overwrite_output: true,
            min_reads: 0,
            distribution: DistributionKind::Normal,
        };
        runner(&config).unwrap();
        assert!(output.exists());
    }

    #[test]
    fn test_runner_zero_tlen_filtered_to_empty_data() {
        // Records with TLEN=0 must be filtered out by the BAM reader (it requires tlen > 0).
        // With every record at TLEN=0 the runner sees zero usable fragments and returns
        // EmptyData — same error variant as a literally empty BAM.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("zero_tlen.bam");
        write_test_frag_bam(&bam_path, &[0usize; 8]);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
            distribution: DistributionKind::Normal,
        };
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenFragLengthModelError::EmptyData),
            "expected EmptyData for all-zero-TLEN BAM, got {err:?}",
        );
    }

    #[test]
    fn test_runner_single_ended_bam_filtered_to_empty_data() {
        // A BAM where reads are not segmented (single-ended sequencing) yields no fragment
        // lengths — the reader's flag filter drops every record. Lock in EmptyData.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("single_ended.bam");
        write_single_ended_frag_bam(&bam_path, &[150usize, 200, 250]);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
            distribution: DistributionKind::Normal,
        };
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenFragLengthModelError::EmptyData),
            "expected EmptyData for single-ended BAM, got {err:?}",
        );
    }

    #[test]
    fn test_runner_low_mapq_filtered_to_empty_data() {
        // Reads with MAPQ <= FRAG_FILTER_MAPQUAL (10) are dropped. With every record at
        // MAPQ=5 we expect EmptyData.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("low_mapq.bam");
        write_frag_bam_with_mapq(&bam_path, &[150usize, 200, 250], 5);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
            distribution: DistributionKind::Normal,
        };
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenFragLengthModelError::EmptyData),
            "expected EmptyData for low-MAPQ BAM, got {err:?}",
        );
    }

    /// Single-ended variant of `write_test_frag_bam`: no SEGMENTED/FIRST_SEGMENT flags set,
    /// otherwise identical. Used to verify the read filter rejects single-ended BAMs.
    #[cfg(test)]
    fn write_single_ended_frag_bam(path: &std::path::PathBuf, tlens: &[usize]) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        let seq = b"ACGT";
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::empty(); // single-ended: no SEGMENTED bit
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    /// Variant of `write_test_frag_bam` that lets the test choose MAPQ. Used to verify the
    /// MAPQ filter in read_fragment_lengths drops low-confidence alignments.
    #[cfg(test)]
    fn write_frag_bam_with_mapq(path: &std::path::PathBuf, tlens: &[usize], mapq: u8) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        let seq = b"ACGT";
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(mapq).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.mate_reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    #[test]
    fn test_runner_empty_bam_errors() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("empty.bam");
        write_test_frag_bam(&bam_path, &[]);
        let output = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            input_file: bam_path,
            output_file: output,
            overwrite_output: true,
            min_reads: 2,
            distribution: DistributionKind::Normal,
        };
        assert!(runner(&config).is_err());
    }

    /// Writes a minimal BGZF BAM with paired, first-in-pair, mapq=30 records,
    /// one record per entry in `tlens`. All records are placed on reference 0.
    #[cfg(test)]
    fn write_test_frag_bam(path: &std::path::PathBuf, tlens: &[usize]) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };

        // Build a header with one reference sequence so refID=0 is valid
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();

        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();

        let seq = b"ACGT";
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.mate_reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }
}
