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
}
