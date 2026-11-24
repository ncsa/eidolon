//! Sampling for a tandem duplication's extra copy (#884).
//!
//! A heterozygous duplication of `B` gives one molecule `A B₁ B₂ C`. The ordinary pass already
//! samples the whole reference at density 1.0, which covers the other chromosome and the path
//! `A B₂ C` of this one. What is left over is the population of fragments that start in `B₁`.
//! Each follows `ref[s..end)` and then continues as `ref[pos..)`, because after the junction the
//! molecule is `B₂` followed by `C`. That is always two pieces, whatever the fragment length:
//!
//! - a fragment that ends at or before `end` never meets the junction, so it is an ordinary
//!   reference fragment `[s, s + len)`;
//! - one that runs past `end` crosses it, with the junction `end - s` bases in.
//!
//! The same draws decide both kinds, so the junction reads and the in-span reads cannot disagree
//! about how many fragments there are. The earlier design raised the start density over the span
//! and emitted a separate, fixed number of junction reads. Fragments from the raised density ran on
//! into the flank, and the junction reads came on top of them: depth past the end of a
//! 500 bp duplication read 1.2 to 1.4 times a no-variant control.
//!
//! A span shorter than a fragment has a second population. A fragment that starts *before* it and
//! runs through it continues into the second copy as well, so on a carrier it too is a junction
//! read. The reference sampling draws it as an ordinary fragment that carries on into the flank,
//! so the regular pass gives it up ([`is_diverted`]) and [`sample_flank_crossers`] draws it as a
//! junction read. A span longer than every fragment has none, and nothing is diverted.
//!
//! The excess may exceed 1.0 (a copy number of 4 or more on a diploid): each additional copy
//! contributes one more population with the same geometry.

use crate::{
    eidolon_core::models::fragment_length::FragmentLengthModel,
    gen_reads::errors::GenerateReadsError,
};
use eidolon_core::rng::NeatRng;

/// One fragment from the extra copy, in reference coordinates of the duplicated span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DupFragment {
    /// 0-based start, inside `[span_start, span_end)`.
    pub start: usize,
    pub frag_len: usize,
}

impl DupFragment {
    /// Whether the fragment runs past the end of the span, into the new tandem junction.
    pub fn crosses_junction(&self, span_end: usize) -> bool {
        self.start + self.frag_len > span_end
    }

    /// Bases before the junction, for a fragment that crosses it.
    pub fn junction_offset(&self, span_end: usize) -> usize {
        span_end - self.start
    }
}

/// Extra copies' contribution to a duplication's depth: the declared multiplier less the
/// baseline 1.0. `None` when there is nothing to add.
pub fn excess_fraction(multiplier: f64) -> Option<f64> {
    let excess = multiplier - 1.0;
    (excess.is_finite() && excess > 0.0).then_some(excess)
}

/// Draws the fragments the extra copy contributes: `excess x span x coverage / denom` of them,
/// starting uniformly over `[span_start, span_end)`, each with a length from the run's model.
///
/// `denom` is `2 x read_length` paired and `read_length` single-ended, as in
/// `generate_fragments`, so the excess is in the same units as the baseline's budget.
#[allow(clippy::too_many_arguments)]
pub fn sample_extra_copy(
    span_start: usize,
    span_end: usize,
    excess: f64,
    coverage: usize,
    read_length: usize,
    paired_ended: bool,
    long_reads: bool,
    keep_short: bool,
    model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Vec<DupFragment>, GenerateReadsError> {
    let width = span_end.saturating_sub(span_start);
    if width == 0 || !excess.is_finite() || excess <= 0.0 || coverage == 0 {
        return Ok(Vec::new());
    }
    let denom = if paired_ended {
        (2 * read_length).max(1)
    } else {
        read_length.max(1)
    };
    let count = (excess * width as f64 * coverage as f64 / denom as f64).round() as usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let frag_len = if paired_ended {
            // The rule `generate_fragments` applies: redraw a fragment too short to hold two
            // reads, and fall back to the shortest that can.
            let mut len = 0usize;
            for _ in 0..100 {
                len = model.generate_fragment(rng.random()?)?;
                if long_reads || keep_short || len >= read_length + 10 {
                    break;
                }
            }
            if len < read_length && !long_reads && !keep_short {
                len = read_length + 10;
            }
            len
        } else {
            read_length
        };
        let start = span_start + ((rng.random()? * width as f64) as usize).min(width - 1);
        out.push(DupFragment { start, frag_len });
    }
    Ok(out)
}

/// The share of the molecules that carry the duplication, for the fragments that start before it.
/// Heterozygous on a diploid is 0.5, homozygous and above is 1.0.
pub fn carrier_share(excess: f64) -> f64 {
    excess.clamp(0.0, 1.0)
}

/// A duplication whose flank crossers are diverted, in reference coordinates `[start, end)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DupSpan {
    pub start: usize,
    pub end: usize,
    pub share: f64,
}

/// Fragments that start before the span and run through it. On a molecule carrying the
/// duplication they continue into the second copy, so they cross the junction exactly like the
/// ones that start inside it. Only a span shorter than a fragment has any: a fragment that starts
/// before `span_start` crosses `span_end` only if it is longer than the span.
///
/// `share` is the carrier share (see [`carrier_share`]). The count is the carriers' fragment
/// density times the window, drawn the same way as [`sample_extra_copy`] and filtered to the ones
/// that reach the end.
#[allow(clippy::too_many_arguments)]
pub fn sample_flank_crossers(
    span_start: usize,
    span_end: usize,
    share: f64,
    coverage: usize,
    read_length: usize,
    paired_ended: bool,
    long_reads: bool,
    keep_short: bool,
    model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Vec<DupFragment>, GenerateReadsError> {
    if span_end <= span_start || !share.is_finite() || share <= 0.0 || coverage == 0 {
        return Ok(Vec::new());
    }
    // How far before the span a fragment can start and still reach it: the longest fragment the
    // model draws (or one read, single-ended). The regular pass diverts every baseline fragment
    // that crosses, wherever it starts, so the window must not be shorter than that.
    let reach = model.max_length()?.max(read_length);
    let window = reach.min(span_start);
    let candidates = sample_extra_copy(
        span_start - window,
        span_start,
        share,
        coverage,
        read_length,
        paired_ended,
        long_reads,
        keep_short,
        model,
        rng,
    )?;
    Ok(candidates
        .into_iter()
        .filter(|f| f.crosses_junction(span_end))
        .collect())
}

/// Whether a baseline fragment `[start, end)` is one the duplication takes over. It starts before
/// the span and ends after it, so on a carrier it continues into the second copy instead of into
/// the flank. Diverted fragments are replaced by [`sample_flank_crossers`]; keeping both would
/// count them twice. Drawn with probability `share`; no draw when it is 1.
pub fn is_diverted(
    start: usize,
    end: usize,
    spans: &[DupSpan],
    rng: &mut NeatRng,
) -> Result<bool, GenerateReadsError> {
    for d in spans {
        let crosses = start < d.start && end > d.end;
        if crosses && (d.share >= 1.0 || rng.random()? < d.share) {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(tag: &str) -> NeatRng {
        NeatRng::new_from_seed(&vec![tag.to_string()]).unwrap()
    }

    fn model(mean: f64) -> FragmentLengthModel {
        FragmentLengthModel::new_normal(mean, 1.0).unwrap()
    }

    /// 1000 bp x 0.5 x coverage 60 over 2 x 100 bp reads is 150 fragments, whatever the lengths.
    #[test]
    fn the_extra_copy_has_excess_times_span_times_coverage_over_denominator_fragments() {
        let frags = sample_extra_copy(
            5_000,
            6_000,
            0.5,
            60,
            100,
            true,
            false,
            false,
            &model(250.0),
            &mut rng("count"),
        )
        .unwrap();
        assert_eq!(frags.len(), 150);
        assert!(frags.iter().all(|f| (5_000..6_000).contains(&f.start)));
        // Single-ended: denom is the read length, and the fragment is one read.
        let se = sample_extra_copy(
            5_000,
            6_000,
            1.0,
            60,
            100,
            false,
            false,
            false,
            &model(250.0),
            &mut rng("se"),
        )
        .unwrap();
        assert_eq!(se.len(), 600);
        assert!(se.iter().all(|f| f.frag_len == 100));
    }

    /// The decision this module exists to make. A fragment crosses the junction exactly when it
    /// starts within one fragment length of the end, so the crossing share of a 1000 bp span at
    /// 250 bp fragments is 250 / 1000, and the number crossing is the number that
    /// the old fixed junction count had to guess at: excess x coverage x length / denom.
    #[test]
    fn the_share_crossing_the_junction_is_one_fragment_length_over_the_span() {
        let frags = sample_extra_copy(
            0,
            1_000,
            2.0,
            1_000,
            100,
            true,
            false,
            false,
            &model(250.0),
            &mut rng("share"),
        )
        .unwrap();
        assert_eq!(frags.len(), 10_000);
        let crossing = frags.iter().filter(|f| f.crosses_junction(1_000)).count();
        let share = crossing as f64 / frags.len() as f64;
        assert!((share - 0.25).abs() < 0.02, "crossing share {share}");
        for f in frags.iter().filter(|f| f.crosses_junction(1_000)) {
            let off = f.junction_offset(1_000);
            assert!(
                off >= 1 && off < f.frag_len,
                "offset {off} of {}",
                f.frag_len
            );
        }
        // Does not fire: nothing that ends inside the span is called a crossing.
        assert!(
            frags
                .iter()
                .filter(|f| !f.crosses_junction(1_000))
                .all(|f| f.start + f.frag_len <= 1_000)
        );
    }

    /// The boundary, exactly: a fragment that ends AT the end of the span has not met the
    /// junction; one base more has. The offset is how many bases lie before it.
    #[test]
    fn a_fragment_ending_at_the_span_end_does_not_cross() {
        let at_end = DupFragment {
            start: 900,
            frag_len: 100,
        };
        assert!(!at_end.crosses_junction(1_000));
        let one_over = DupFragment {
            start: 900,
            frag_len: 101,
        };
        assert!(one_over.crosses_junction(1_000));
        assert_eq!(one_over.junction_offset(1_000), 100);
        let first_base = DupFragment {
            start: 999,
            frag_len: 250,
        };
        assert!(first_base.crosses_junction(1_000));
        assert_eq!(first_base.junction_offset(1_000), 1);
    }

    /// A span shorter than a fragment: every fragment crosses, and the second piece is longer than
    /// the span, which is fine because the molecule continues as reference after the junction.
    #[test]
    fn a_span_shorter_than_a_fragment_is_crossed_by_all_of_them() {
        let frags = sample_extra_copy(
            0,
            100,
            1.0,
            60,
            100,
            true,
            false,
            false,
            &model(250.0),
            &mut rng("short"),
        )
        .unwrap();
        assert_eq!(frags.len(), 30);
        assert!(frags.iter().all(|f| f.crosses_junction(100)));
    }

    #[test]
    fn nothing_is_drawn_without_an_excess_or_coverage() {
        let m = model(250.0);
        let mut r = rng("none");
        for (excess, cov, a, b) in [(0.0, 60, 0, 100), (-0.5, 60, 0, 100), (0.5, 0, 0, 100)] {
            assert!(
                sample_extra_copy(a, b, excess, cov, 100, true, false, false, &m, &mut r)
                    .unwrap()
                    .is_empty()
            );
        }
        assert!(
            sample_extra_copy(7, 7, 0.5, 60, 100, true, false, false, &m, &mut r)
                .unwrap()
                .is_empty()
        );
        assert_eq!(excess_fraction(1.0), None);
        assert_eq!(excess_fraction(0.5), None);
        assert_eq!(excess_fraction(f64::NAN), None);
        assert_eq!(excess_fraction(1.5), Some(0.5));
    }

    /// Both passes call this with the same seeded child, so they must see the same fragments.
    #[test]
    fn the_same_seed_gives_the_same_fragments() {
        let a = sample_extra_copy(
            0,
            800,
            0.5,
            30,
            100,
            true,
            false,
            false,
            &model(300.0),
            &mut rng("repeat"),
        )
        .unwrap();
        let b = sample_extra_copy(
            0,
            800,
            0.5,
            30,
            100,
            true,
            false,
            false,
            &model(300.0),
            &mut rng("repeat"),
        )
        .unwrap();
        assert_eq!(a, b);
    }

    /// A span of 100 bp under 250 bp fragments: fragments that start before it and reach past it.
    /// Each fragment length has `length - span` such starts per base of density.
    #[test]
    fn flank_crossers_start_before_the_span_and_reach_past_it() {
        let frags = sample_flank_crossers(
            5_000,
            5_100,
            1.0,
            1_000,
            100,
            true,
            false,
            false,
            &model(250.0),
            &mut rng("flank"),
        )
        .unwrap();
        // density 1000 / 200 per base, (250 - 100) starts per fragment length.
        let expected = 1_000.0 / 200.0 * 150.0;
        let n = frags.len() as f64;
        assert!((n / expected - 1.0).abs() < 0.10, "{n} vs {expected}");
        for f in &frags {
            assert!(f.start < 5_000 && f.crosses_junction(5_100), "{f:?}");
            let off = f.junction_offset(5_100);
            assert!(off > 100 && off < f.frag_len, "offset {off}");
        }
    }

    /// Does not fire: a span longer than every fragment cannot be started before and crossed.
    #[test]
    fn a_span_longer_than_a_fragment_has_no_flank_crossers() {
        let frags = sample_flank_crossers(
            5_000,
            6_000,
            1.0,
            1_000,
            100,
            true,
            false,
            false,
            &model(250.0),
            &mut rng("long"),
        )
        .unwrap();
        assert!(frags.is_empty());
    }

    #[test]
    fn the_carrier_share_is_the_excess_capped_at_one() {
        assert_eq!(carrier_share(0.5), 0.5);
        assert_eq!(carrier_share(1.0), 1.0);
        assert_eq!(carrier_share(2.0), 1.0);
        assert_eq!(carrier_share(-1.0), 0.0);
    }

    /// The decision, at every edge: only a fragment that starts before the span AND ends after it.
    #[test]
    fn only_a_fragment_that_starts_before_and_ends_after_the_span_is_diverted() {
        let spans = [DupSpan {
            start: 1_000,
            end: 1_100,
            share: 1.0,
        }];
        let mut r = rng("div");
        let mut d = |s, e| is_diverted(s, e, &spans, &mut r).unwrap();
        assert!(d(999, 1_101), "one base each side");
        assert!(!d(1_000, 1_300), "starts inside the span");
        assert!(!d(900, 1_100), "ends at the end, has not crossed");
        assert!(!d(900, 1_000), "stops before the span");
        assert!(!d(1_200, 1_500), "past the span");
        assert!(!d(0, 50));
    }

    #[test]
    fn a_diverted_fragment_is_drawn_with_the_carrier_share() {
        let spans = [DupSpan {
            start: 1_000,
            end: 1_100,
            share: 0.5,
        }];
        let mut r = rng("coin");
        let n = (0..10_000)
            .filter(|_| is_diverted(900, 1_200, &spans, &mut r).unwrap())
            .count();
        assert!((4_700..=5_300).contains(&n), "{n} of 10000");
        // A fragment that cannot be diverted draws nothing, so the stream is untouched.
        let mut a = rng("same");
        let mut b = rng("same");
        assert!(!is_diverted(0, 50, &spans, &mut a).unwrap());
        assert_eq!(a.random().unwrap(), b.random().unwrap());
    }

    /// A long-read model (every fragment 5 kb) over a 100 bp span. Fragments crossing the span
    /// start anywhere in the 4.9 kb before it, and the regular pass diverts all of them, so the
    /// replacements must come from the whole reach and not from a fixed window:
    /// share x coverage x 4900 / denom = 0.5 x 60 x 4900 / 300 = 490 (#885 review).
    #[test]
    fn flank_crossers_cover_the_whole_reach_of_a_long_fragment_model() {
        let long = FragmentLengthModel::new_discrete(vec![5_000], vec![1.0]).unwrap();
        let got = sample_flank_crossers(
            20_000,
            20_100,
            0.5,
            60,
            150,
            true,
            true,
            false,
            &long,
            &mut rng("long"),
        )
        .unwrap();
        assert!(
            (got.len() as f64 - 490.0).abs() < 490.0 * 0.05,
            "{} crossers, want about 490",
            got.len()
        );
        let far = got.iter().filter(|f| f.start < 20_000 - 2_000).count();
        assert!(far > 200, "only {far} start more than 2 kb before the span");
        // Must not fire: a short-fragment model still draws only what can reach the span.
        let short = sample_flank_crossers(
            20_000,
            20_100,
            0.5,
            60,
            150,
            true,
            false,
            false,
            &model(350.0),
            &mut rng("short"),
        )
        .unwrap();
        assert!(short.iter().all(|f| f.start + 360 > 20_000));
    }
}
