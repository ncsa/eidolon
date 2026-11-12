//! What `coverage` means on a `target_bed` run (#578): the mean depth over the target bases.
//!
//! Placement owns a fragment by where it STARTS, so a fragment owned by a narrow target runs
//! past the target's far edge into flanking sequence, and those read bases are not on-target
//! depth. [`on_target_fraction`] estimates the share of read bases that land on target, from the
//! BED geometry and the run's fragment model; the fragment budget is divided by it.
//!
//! The estimate counts a base as on-target if ANY target covers it. That is what keeps the
//! correction local: a neighbor that catches the spill makes the fraction rise toward 1, so
//! touching targets and the segments of an SV depth modulation inside one target are barely
//! changed, and a wide window (a sharded whole-genome run) is changed by about
//! fragment length / window width.

use crate::{
    eidolon_core::models::fragment_length::FragmentLengthModel,
    gen_reads::errors::GenerateReadsError,
};
use eidolon_core::rng::NeatRng;

/// One region whose fragments the placement code owns, in contig coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct OwnedRegion {
    pub start: usize,
    pub end: usize,
    /// Real, materializable sequence past `end` before a terminus (contig end or N gap): how far
    /// a fragment's far end may run. This is `extension_budget` in `generate_fragments`.
    pub extension: usize,
}

/// A contig's targets: the regions that own fragments, and the merged union of target bases
/// that count as on-target.
#[derive(Debug, Clone, Default)]
pub struct ContigTargets {
    pub owned: Vec<OwnedRegion>,
    pub on_target: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Copy)]
pub struct SpillParams {
    pub read_len: usize,
    pub paired_ended: bool,
    /// Keep inserts shorter than `read_len + 10` (adapters or `keep_short_fragments`).
    pub keep_short: bool,
    pub long_reads: bool,
}

/// Fragments sampled for the estimate. Each contributes an on-target fraction with a standard
/// deviation near 0.4, so the standard error of the result is about 0.0004 absolute: well inside
/// the rounding of a depth people quote.
const SAMPLES: usize = 1_000_000;
/// Child stream of the run RNG, so the estimate consumes nothing from the streams that place
/// fragments and a seeded run stays reproducible.
pub const TARGET_DEPTH_STREAM: u64 = 9_200_000;
/// Redraw cap for a fragment length rejected by the paired-end minimum, as placement does.
const MAX_REDRAWS: usize = 64;

/// Fraction of read bases, over fragments the placement code would place, that fall on target.
/// `None` when no fragment can be placed, so there is nothing to correct.
///
/// This mirrors `cover_dataset`: a region owns a fragment by its start, a start is uniform over
/// `[0, min(width, width + extension - length + 1))`, and a fragment longer than
/// `width + extension` is dropped (it is in neither the numerator nor the denominator, because it
/// produces no reads).
pub fn on_target_fraction(
    contigs: &[ContigTargets],
    params: SpillParams,
    model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Option<f64>, GenerateReadsError> {
    // Cumulative owned widths, so a region is picked in proportion to the fragments it owns
    // (the budget is width x coverage).
    // A region that cannot hold a read even with the room past its far edge owns no fragments
    // (`generate_fragments` returns early), so it takes no part in the estimate.
    let min_region = if params.long_reads && params.paired_ended {
        0
    } else {
        params.read_len
    };
    let mut regions: Vec<(&ContigTargets, &OwnedRegion)> = Vec::new();
    let mut cumulative: Vec<usize> = Vec::new();
    let mut total_width = 0usize;
    for c in contigs {
        for r in &c.owned {
            if r.end > r.start && (r.end - r.start).saturating_add(r.extension) > min_region {
                total_width += r.end - r.start;
                regions.push((c, r));
                cumulative.push(total_width);
            }
        }
    }
    if total_width == 0 || params.read_len == 0 {
        return Ok(None);
    }

    let min_paired = params.read_len + 10;
    let (mut on_target, mut produced) = (0u128, 0u128);
    for _ in 0..SAMPLES {
        let pick = ((rng.random()? * total_width as f64) as usize).min(total_width - 1);
        let (contig, region) = regions[cumulative.partition_point(|&c| c <= pick)];
        let width = region.end - region.start;

        let length = if params.paired_ended {
            let mut drawn = None;
            for _ in 0..MAX_REDRAWS {
                let l = model.generate_fragment(rng.random()?)?;
                if params.long_reads || params.keep_short || l >= min_paired {
                    drawn = Some(l);
                    break;
                }
            }
            match drawn {
                Some(l) => l,
                None => continue,
            }
        } else {
            params.read_len
        };
        let far_limit = width + region.extension;
        if length == 0 || length > far_limit {
            continue;
        }
        let start_hi = width.min(far_limit - length + 1);
        let start = region.start + ((rng.random()? * start_hi as f64) as usize).min(start_hi - 1);
        let end = start + length;

        let read = params.read_len.min(length);
        let covered = if params.paired_ended {
            overlap(&contig.on_target, start, start + read)
                + overlap(&contig.on_target, end - read, end)
        } else {
            overlap(&contig.on_target, start, start + read)
        };
        on_target += covered as u128;
        produced += (read * if params.paired_ended { 2 } else { 1 }) as u128;
    }
    if produced == 0 {
        return Ok(None);
    }
    Ok(Some(on_target as f64 / produced as f64))
}

/// Bases of `[lo, hi)` inside the sorted, non-overlapping `intervals`.
fn overlap(intervals: &[(usize, usize)], lo: usize, hi: usize) -> usize {
    let first = intervals.partition_point(|&(_, e)| e <= lo);
    let mut covered = 0;
    for &(s, e) in &intervals[first..] {
        if s >= hi {
            break;
        }
        covered += e.min(hi) - s.max(lo);
    }
    covered
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng() -> NeatRng {
        NeatRng::new_from_seed(&vec!["target-depth-test".to_string()]).unwrap()
    }

    /// Every fragment is exactly `len` long.
    fn fixed(len: usize) -> FragmentLengthModel {
        FragmentLengthModel::new_discrete(vec![len], vec![1.0]).unwrap()
    }

    fn paired(read_len: usize) -> SpillParams {
        SpillParams {
            read_len,
            paired_ended: true,
            keep_short: false,
            long_reads: false,
        }
    }

    fn region(start: usize, end: usize, extension: usize) -> OwnedRegion {
        OwnedRegion {
            start,
            end,
            extension,
        }
    }

    fn f(contigs: &[ContigTargets], p: SpillParams, m: &FragmentLengthModel) -> f64 {
        on_target_fraction(contigs, p, m, &mut rng())
            .unwrap()
            .expect("a fragment can be placed")
    }

    /// One isolated 100 bp target, 50 bp reads, 200 bp fragments, plenty of sequence beyond.
    /// R2 starts at s+150 >= 150, past the target, so only R1 can score; its on-target bases are
    /// min(50, 100-s) for s uniform on 0..100. Summed: 2550 + 1225 = 3775 over 100 starts, out
    /// of 100 read bases per fragment: 0.3775. Worked by hand, independent of the code.
    #[test]
    fn isolated_target_loses_its_spill() {
        let c = ContigTargets {
            owned: vec![region(0, 100, 9_900)],
            on_target: vec![(0, 100)],
        };
        let got = f(&[c], paired(50), &fixed(200));
        assert!((got - 0.3775).abs() < 0.004, "got {got}");
    }

    /// Same geometry, but the target is two touching 100 bp regions (one merged on-target
    /// interval 0..200). Spill from the first lands on the second. Worked by hand:
    /// first region 6275 + second region 3775 = 10050 of 20000 read bases: 0.5025.
    #[test]
    fn a_neighbor_that_catches_the_spill_raises_the_fraction() {
        let c = ContigTargets {
            owned: vec![region(0, 100, 9_900), region(100, 200, 9_800)],
            on_target: vec![(0, 200)],
        };
        let got = f(&[c], paired(50), &fixed(200));
        assert!((got - 0.5025).abs() < 0.004, "got {got}");
    }

    /// A target that runs to the contig end has nothing to spill into: placement thins the
    /// starts instead (start < W + extension - L + 1), so every read is on target.
    #[test]
    fn a_target_at_the_contig_end_loses_nothing() {
        let c = ContigTargets {
            owned: vec![region(0, 100, 0)],
            on_target: vec![(0, 100)],
        };
        let got = f(&[c], paired(50), &fixed(60));
        assert!((got - 1.0).abs() < 1e-12, "got {got}");
    }

    /// Single-end: the fragment is the read. min(50, 100-s) averages 37.75 of 50: 0.755.
    #[test]
    fn single_end_uses_the_read_as_the_fragment() {
        let c = ContigTargets {
            owned: vec![region(0, 100, 9_900)],
            on_target: vec![(0, 100)],
        };
        let p = SpillParams {
            paired_ended: false,
            ..paired(50)
        };
        let got = f(&[c], p, &fixed(999));
        assert!((got - 0.755).abs() < 0.004, "got {got}");
    }

    /// A very wide window, the shape of a sharded whole-genome run: only the far edge loses
    /// anything, about fragment length / width.
    #[test]
    fn a_wide_window_is_corrected_by_almost_nothing() {
        let w = 25_000_000;
        let c = ContigTargets {
            owned: vec![region(0, w, 1_000_000)],
            on_target: vec![(0, w)],
        };
        let got = f(&[c], paired(150), &fixed(400));
        assert!(got > 1.0 - 1e-4 && got <= 1.0, "got {got}");
    }

    /// `generate_fragments` returns early for a region that cannot hold a read even with the room
    /// past its far edge, whatever the fragment length, so such a region owns nothing and must
    /// not take part. Only `keep_short` can tell: its 30 bp fragments would fit a 40 bp region.
    /// Region A alone: both reads cover the 30 bp fragment, min(30, 100-s) on target for s on
    /// 0..100, 25.65 of 30 on average, 0.855. Counting the 40 bp region (all on target) would
    /// give 0.896.
    #[test]
    fn a_region_that_owns_no_fragments_takes_no_part() {
        let c = ContigTargets {
            owned: vec![region(0, 100, 9_900), region(20_000, 20_040, 0)],
            on_target: vec![(0, 100), (20_000, 20_040)],
        };
        let p = SpillParams {
            keep_short: true,
            ..paired(50)
        };
        let got = f(&[c], p, &fixed(30));
        assert!((got - 0.855).abs() < 0.01, "got {got}");
    }

    #[test]
    fn nothing_placeable_gives_none() {
        let none = on_target_fraction(&[], paired(50), &fixed(200), &mut rng()).unwrap();
        assert!(none.is_none());
        // Fragments longer than the region plus its extension cannot be placed.
        let c = ContigTargets {
            owned: vec![region(0, 100, 0)],
            on_target: vec![(0, 100)],
        };
        let none = on_target_fraction(&[c], paired(50), &fixed(5_000), &mut rng()).unwrap();
        assert!(none.is_none());
    }

    #[test]
    fn the_estimate_is_reproducible_for_a_seed() {
        let c = ContigTargets {
            owned: vec![region(0, 178, 9_000)],
            on_target: vec![(0, 178)],
        };
        let m = FragmentLengthModel::new_normal(250.0, 30.0).unwrap();
        assert_eq!(
            f(std::slice::from_ref(&c), paired(100), &m),
            f(std::slice::from_ref(&c), paired(100), &m)
        );
    }
}
