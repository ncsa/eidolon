//! Read-level realism metrics, computed identically for a real BAM and a simulated one.
//!
//! WHY THIS EXISTS: eidolon's reads are too clean, and nothing measured how much. On matched
//! chr22 sequence at 30x, real NA12878 against eidolon 3.2.0:
//!
//! | metric | real | simulated |
//! |---|---|---|
//! | candidate breakpoints / Mb | **654.5** | **0.0** |
//! | improper pairs | 3.15% | 0.00% |
//! | soft clip >= 20 bp | 0.67% | 0.00% |
//! | depth VMR | 5.5–8.9 | 1.03–1.12 |
//!
//! Zero, not "fewer". A caller tuned on eidolon calibrates every clustering threshold and FP
//! filter against an empty background, then meets 654 candidates per Mb on real data. That is
//! how a simulator produces a tool that works on synthetic data and drowns on real data — and
//! it means eidolon-measured *precision* has never predicted real precision. Recall is a
//! different question: planted events are real signal.
//!
//! DESIGN: one implementation, both inputs. If "real" and "simulated" had separate code paths
//! they would drift, and the comparison would measure the drift. Everything here is a pure
//! function over `AlnRecord` so it is testable without a BAM, an aligner, or a reference —
//! the metrics are the part that has to be right, and they should not need 6 GB of inputs to
//! exercise.

/// The fields of an alignment this module needs. Deliberately not a `noodles` record: the
/// arithmetic below is what must be correct, and it should be testable from literals.
#[derive(Debug, Clone, PartialEq)]
pub struct AlnRecord {
    /// 0-based leftmost mapped position.
    pub pos: usize,
    pub mapq: u8,
    /// SAM FLAG.
    pub flags: u16,
    /// CIGAR as (operation, length) in order.
    pub cigar: Vec<(char, usize)>,
    /// TLEN / observed template length, signed as in SAM.
    pub tlen: i64,
}

impl AlnRecord {
    /// Reference bases consumed: `M`, `D`, `N`, `=`, `X`. Insertions and clips consume query
    /// only, which is exactly why a deletion cannot be expressed as a per-query-base op.
    pub fn reference_span(&self) -> usize {
        self.cigar
            .iter()
            .filter(|(op, _)| matches!(op, 'M' | 'D' | 'N' | '=' | 'X'))
            .map(|(_, n)| n)
            .sum()
    }

    /// Length of a soft clip at the start of the alignment, if any.
    pub fn leading_clip(&self) -> usize {
        match self.cigar.first() {
            Some(('S', n)) => *n,
            _ => 0,
        }
    }

    /// Length of a soft clip at the end of the alignment, if any.
    pub fn trailing_clip(&self) -> usize {
        match self.cigar.last() {
            Some(('S', n)) => *n,
            _ => 0,
        }
    }

    pub fn is_proper_pair(&self) -> bool {
        self.flags & 0x2 != 0
    }
}

/// What one region contributed. Reported per region and never pooled: a metric averaged over
/// loci hides the locus that disagrees, and locus-to-locus spread is what any future threshold
/// has to be derived from.
#[derive(Debug, Clone, PartialEq)]
pub struct RegionMetrics {
    /// Denominator. A rate over an unstated number of reads is not a result (rule 4).
    pub reads: usize,
    pub span_bp: usize,
    pub candidate_breakpoints: usize,
    pub improper_pairs: usize,
    pub clipped_reads: usize,
    pub mapq0_reads: usize,
    pub insert: Option<InsertStats>,
    pub depth: Option<DepthStats>,
}

impl RegionMetrics {
    /// Candidate breakpoints per megabase — the headline. This is the raw material every
    /// split-read SV caller clusters on, so it is the closest single number to "how much
    /// harder is real data than ours".
    pub fn candidates_per_mb(&self) -> f64 {
        if self.span_bp == 0 {
            return 0.0;
        }
        self.candidate_breakpoints as f64 * 1_000_000.0 / self.span_bp as f64
    }

    pub fn improper_pair_rate(&self) -> f64 {
        rate(self.improper_pairs, self.reads)
    }
    pub fn clip_rate(&self) -> f64 {
        rate(self.clipped_reads, self.reads)
    }
    pub fn mapq0_rate(&self) -> f64 {
        rate(self.mapq0_reads, self.reads)
    }
}

fn rate(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct InsertStats {
    pub n: usize,
    pub mean: f64,
    pub sd: f64,
    /// Real libraries are right-skewed; a Normal fragment model is symmetric. Discordant-pair
    /// callers threshold on insert size, so the SHAPE matters more than the mean — the mean is
    /// a config choice, the skew and tail are not.
    pub skew: f64,
    pub p99: i64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DepthStats {
    pub mean: f64,
    /// Variance-to-mean ratio. Poisson is 1.0. Real WGS sits well above it because coverage
    /// is modulated by GC, mappability and library prep; a simulator at ~1.0 is telling a CNV
    /// caller that segmentation is far easier than it is.
    pub vmr: f64,
    /// Autocorrelation of the depth track at a fixed lag. Uncorrelated noise averages away
    /// inside a bin; correlated noise does not, so this drives a binned caller's false
    /// positives independently of the VMR.
    pub autocorr_500: f64,
}

impl DepthStats {
    /// Dispersion in EXCESS of Poisson, normalized by depth: `(vmr - 1) / mean`.
    ///
    /// VMR alone is not comparable between datasets at different depths. Coverage noise is
    /// largely multiplicative — depth is the mean times a local factor from GC, mappability
    /// and library prep — so variance grows with the square of the mean and VMR grows with the
    /// mean. Measured: real NA12878 chr22 at 247x reads VMR 36.1 against eidolon's 1.04 at
    /// 30x, which looks like a 35x gap and is not one; the depths differ by 8x.
    ///
    /// Subtracting 1 removes the Poisson floor (a pure counting process has VMR 1 at any
    /// depth) and dividing by the mean removes the depth scaling, leaving the squared
    /// coefficient of variation of the underlying modulation. On the same pair that is 0.142
    /// against 0.0012 — a **118x** gap, and it says something stronger: eidolon's dispersion
    /// is almost entirely Poisson counting noise (1/30 = 0.033 of its 0.035 CV^2), meaning
    /// there is essentially no coverage modulation at all.
    ///
    /// This is the number to compare across datasets. `vmr` is kept because it is what the
    /// literature quotes and what a depth caller experiences at ITS depth.
    pub fn excess_dispersion(&self) -> f64 {
        if self.mean <= 0.0 {
            return 0.0;
        }
        (self.vmr - 1.0) / self.mean
    }
}

/// A position where at least `min_support` reads share a soft-clip boundary.
///
/// Clip boundaries, not clipped reads: a caller clusters *positions*, and scattered clipping
/// with no agreement is noise rather than a candidate junction. Leading and trailing clips are
/// counted at different coordinates — the alignment start for a leading clip, the end of the
/// reference span for a trailing one — because they mark opposite sides of a junction.
pub fn candidate_breakpoints(records: &[AlnRecord], min_clip: usize, min_support: usize) -> usize {
    candidate_sites(records, min_clip, min_support).len()
}

/// One side of a clip boundary that enough reads agree on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateSite {
    pub pos: usize,
    /// `L` for a leading clip (junction to the left), `R` for a trailing one.
    pub side: char,
    /// How many reads share this exact boundary.
    pub support: usize,
    /// Reads overlapping the position at all, clipped or not. The denominator: 3 reads
    /// agreeing out of 4 is a different thing from 3 out of 90, and only one of them looks
    /// like a real junction.
    pub depth: usize,
    /// How many of those overlapping reads are MAPQ 0. A boundary in a repeat looks
    /// different from a boundary at a structural variant, and this is what separates them.
    pub mapq0: usize,
}

/// The sites themselves, not just the count.
///
/// Split out so `--dump-candidates` reports exactly what the metric counted rather than a
/// reimplementation of it. The count is the length of this list, by construction.
pub fn candidate_sites(
    records: &[AlnRecord],
    min_clip: usize,
    min_support: usize,
) -> Vec<CandidateSite> {
    use std::collections::HashMap;
    let mut left: HashMap<usize, usize> = HashMap::new();
    let mut right: HashMap<usize, usize> = HashMap::new();
    for r in records {
        if r.leading_clip() >= min_clip {
            *left.entry(r.pos).or_insert(0) += 1;
        }
        if r.trailing_clip() >= min_clip {
            *right.entry(r.pos + r.reference_span()).or_insert(0) += 1;
        }
    }

    let mut out: Vec<CandidateSite> = Vec::new();
    for (side, map) in [('L', &left), ('R', &right)] {
        for (&pos, &support) in map.iter() {
            if support < min_support {
                continue;
            }
            let (mut depth, mut mapq0) = (0usize, 0usize);
            for r in records {
                if r.pos <= pos && pos < r.pos + r.reference_span() {
                    depth += 1;
                    if r.mapq == 0 {
                        mapq0 += 1;
                    }
                }
            }
            out.push(CandidateSite {
                pos,
                side,
                support,
                depth,
                mapq0,
            });
        }
    }
    out.sort_by_key(|c| (c.pos, c.side));
    out
}

/// Insert-size distribution over library-scale pairs only.
///
/// `max_tlen` excludes inter-chromosomal and far-apart pairs, which otherwise dominate the
/// moments: measured on real chr22 without a bound, the "insert mean" came out at 3606 bp with
/// a skew of +88, describing structural rearrangement rather than the library.
pub fn insert_stats(records: &[AlnRecord], max_tlen: i64) -> Option<InsertStats> {
    let mut v: Vec<i64> = records
        .iter()
        .map(|r| r.tlen)
        .filter(|t| *t > 0 && *t < max_tlen)
        .collect();
    if v.len() < 2 {
        return None;
    }
    v.sort_unstable();
    let n = v.len();
    let mean = v.iter().sum::<i64>() as f64 / n as f64;
    let var = v.iter().map(|x| (*x as f64 - mean).powi(2)).sum::<f64>() / n as f64;
    let sd = var.sqrt();
    let skew = if sd > 0.0 {
        v.iter().map(|x| (*x as f64 - mean).powi(3)).sum::<f64>() / n as f64 / sd.powi(3)
    } else {
        0.0
    };
    Some(InsertStats {
        n,
        mean,
        sd,
        skew,
        p99: v[(0.99 * n as f64) as usize % n],
    })
}

/// Per-base depth over `[start, start + len)` from the records' reference spans.
pub fn depth_track(records: &[AlnRecord], start: usize, len: usize) -> Vec<u32> {
    let mut d = vec![0u32; len];
    for r in records {
        let s = r.pos.max(start);
        let e = (r.pos + r.reference_span()).min(start + len);
        for slot in d.iter_mut().take(e.saturating_sub(start)).skip(s - start) {
            if e > s {
                *slot += 1;
            }
        }
    }
    d
}

pub fn depth_stats(track: &[u32], lag: usize) -> Option<DepthStats> {
    if track.len() <= lag || track.is_empty() {
        return None;
    }
    let n = track.len() as f64;
    let mean = track.iter().map(|x| *x as f64).sum::<f64>() / n;
    if mean == 0.0 {
        return None;
    }
    let var = track
        .iter()
        .map(|x| (*x as f64 - mean).powi(2))
        .sum::<f64>()
        / n;
    let denom: f64 = track.iter().map(|x| (*x as f64 - mean).powi(2)).sum();
    let num: f64 = (0..track.len() - lag)
        .map(|i| (track[i] as f64 - mean) * (track[i + lag] as f64 - mean))
        .sum();
    Some(DepthStats {
        mean,
        vmr: var / mean,
        autocorr_500: if denom > 0.0 { num / denom } else { 0.0 },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cig(s: &str) -> Vec<(char, usize)> {
        let mut out = Vec::new();
        let mut n = String::new();
        for c in s.chars() {
            if c.is_ascii_digit() {
                n.push(c);
            } else {
                out.push((c, n.parse().unwrap()));
                n.clear();
            }
        }
        out
    }

    fn rec(pos: usize, cigar: &str, mapq: u8, proper: bool, tlen: i64) -> AlnRecord {
        AlnRecord {
            pos,
            mapq,
            flags: if proper { 0x2 } else { 0 },
            cigar: cig(cigar),
            tlen,
        }
    }

    // ── reference span: the arithmetic every other metric is built on ──────────

    #[test]
    fn reference_span_counts_only_reference_consuming_ops() {
        assert_eq!(rec(0, "151M", 60, true, 0).reference_span(), 151);
        // I and S consume query, not reference.
        assert_eq!(rec(0, "50M20I81M", 60, true, 0).reference_span(), 131);
        assert_eq!(rec(0, "20S131M", 60, true, 0).reference_span(), 131);
        // D and N consume reference, not query.
        assert_eq!(rec(0, "50M500D101M", 60, true, 0).reference_span(), 651);
        // H consumes neither.
        assert_eq!(rec(0, "20H131M", 60, true, 0).reference_span(), 131);
    }

    #[test]
    fn clips_are_read_from_the_correct_end() {
        let r = rec(0, "30S91M30S", 60, true, 0);
        assert_eq!(r.leading_clip(), 30);
        assert_eq!(r.trailing_clip(), 30);
        // A hard clip is not a soft clip: its bases are absent from SEQ entirely.
        let h = rec(0, "30H121M", 60, true, 0);
        assert_eq!(h.leading_clip(), 0);
    }

    // ── candidate breakpoints: the headline metric ─────────────────────────────

    #[test]
    fn the_dumped_sites_are_exactly_what_the_count_counted() {
        // The count and the dump are two views of one thing, and a dump that disagreed with
        // the headline would send someone hunting for a cause that is not there. Asserted
        // rather than assumed, because `candidate_breakpoints` delegating to
        // `candidate_sites` is an implementation detail that a later edit could undo.
        let mut v = Vec::new();
        // three reads sharing a leading-clip boundary at 500 -> one candidate
        for _ in 0..3 {
            v.push(rec(500, "30S100M", 60, true, 300));
        }
        // two sharing one at 900 -> below min_support, not a candidate
        for _ in 0..2 {
            v.push(rec(900, "30S100M", 60, true, 300));
        }
        // four sharing a trailing boundary; 700 + 100 = 800
        for _ in 0..4 {
            v.push(rec(700, "100M30S", 60, true, 300));
        }
        let sites = candidate_sites(&v, 20, 3);
        assert_eq!(
            sites.len(),
            candidate_breakpoints(&v, 20, 3),
            "dump and count disagree"
        );
        let described: Vec<(usize, char, usize)> =
            sites.iter().map(|c| (c.pos, c.side, c.support)).collect();
        assert_eq!(described, vec![(500, 'L', 3), (800, 'R', 4)]);
    }

    #[test]
    fn a_dumped_site_carries_the_denominator_its_support_is_relative_to() {
        // 3 of 4 reads agreeing is a junction; 3 of 90 is noise. Without the local depth the
        // dump cannot tell those apart, which is the whole reason to dump rather than count.
        let mut v = Vec::new();
        for _ in 0..3 {
            v.push(rec(500, "30S100M", 60, true, 300));
        }
        // 20 unclipped reads spanning position 500, half of them unmappable
        for i in 0..20 {
            v.push(rec(450, "100M", if i < 10 { 0 } else { 60 }, true, 300));
        }
        let sites = candidate_sites(&v, 20, 3);
        assert_eq!(sites.len(), 1);
        let c = &sites[0];
        assert_eq!(c.support, 3);
        assert_eq!(c.depth, 23, "every read overlapping the position counts");
        assert_eq!(
            c.mapq0, 10,
            "a repeat-driven boundary is what mapq0 separates"
        );
    }

    #[test]
    fn candidate_breakpoints_need_agreement_not_just_clipping() {
        // Five clipped reads, all at DIFFERENT positions: scattered noise, no candidate.
        let scattered: Vec<_> = (0..5)
            .map(|i| rec(1000 + i * 37, "40S111M", 60, true, 0))
            .collect();
        assert_eq!(candidate_breakpoints(&scattered, 20, 3), 0);

        // Five clipped reads agreeing on one position: one candidate.
        let agreeing: Vec<_> = (0..5).map(|_| rec(1000, "40S111M", 60, true, 0)).collect();
        assert_eq!(candidate_breakpoints(&agreeing, 20, 3), 1);
    }

    #[test]
    fn leading_and_trailing_clips_mark_different_coordinates() {
        // Leading clips register at the alignment start; trailing at the end of the span.
        // A junction has two sides and they must not collapse into one count.
        let mut v: Vec<_> = (0..3).map(|_| rec(1000, "40S111M", 60, true, 0)).collect();
        v.extend((0..3).map(|_| rec(500, "111M40S", 60, true, 0)));
        assert_eq!(candidate_breakpoints(&v, 20, 3), 2);
        // The trailing group registers at 500 + 111 = 611, not at 500.
        let only_trailing: Vec<_> = (0..3).map(|_| rec(500, "111M40S", 60, true, 0)).collect();
        assert_eq!(candidate_breakpoints(&only_trailing, 20, 3), 1);
    }

    #[test]
    fn a_clip_below_the_threshold_is_not_a_candidate() {
        let short: Vec<_> = (0..5).map(|_| rec(1000, "5S146M", 60, true, 0)).collect();
        assert_eq!(candidate_breakpoints(&short, 20, 3), 0);
    }

    #[test]
    fn support_below_the_threshold_is_not_a_candidate() {
        let two: Vec<_> = (0..2).map(|_| rec(1000, "40S111M", 60, true, 0)).collect();
        assert_eq!(candidate_breakpoints(&two, 20, 3), 0);
    }

    // ── insert stats ───────────────────────────────────────────────────────────

    #[test]
    fn insert_stats_exclude_far_and_negative_pairs() {
        // Without a bound, one inter-chromosomal pair dominates the moments. Measured on real
        // chr22 that produced an "insert mean" of 3606 bp with skew +88 — a description of
        // rearrangement, not of the library.
        let mut v: Vec<_> = (0..100)
            .map(|i| rec(0, "151M", 60, true, 400 + (i % 7) as i64))
            .collect();
        v.push(rec(0, "151M", 60, false, 90_000_000));
        v.push(rec(0, "151M", 60, true, -400)); // the mate's negative TLEN, counted once only
        let s = insert_stats(&v, 2000).unwrap();
        assert_eq!(s.n, 100, "far and negative TLENs must be excluded");
        assert!((s.mean - 403.0).abs() < 1.0, "mean was {}", s.mean);
    }

    #[test]
    fn insert_skew_distinguishes_a_symmetric_model_from_a_real_library() {
        // This is the whole point of measuring skew: eidolon's fragment model is Normal, real
        // libraries are right-tailed, and a discordant-pair caller thresholds on the tail.
        let symmetric: Vec<_> = (0..1000)
            .map(|i| rec(0, "151M", 60, true, 400 + ((i % 21) as i64 - 10)))
            .collect();
        let sym = insert_stats(&symmetric, 5000).unwrap();
        assert!(sym.skew.abs() < 0.2, "symmetric input skewed {}", sym.skew);

        // Same centre, long right tail.
        let mut skewed: Vec<_> = (0..950).map(|_| rec(0, "151M", 60, true, 400)).collect();
        skewed.extend((0..50).map(|i| rec(0, "151M", 60, true, 900 + i as i64 * 10)));
        let sk = insert_stats(&skewed, 5000).unwrap();
        assert!(
            sk.skew > 1.0,
            "right-tailed input should skew positive, got {}",
            sk.skew
        );
        assert!(sk.p99 > sym.p99, "the tail must show up in p99");
    }

    // ── depth ──────────────────────────────────────────────────────────────────

    #[test]
    fn depth_track_counts_reference_span_not_read_length() {
        // A read with a 500 bp deletion covers 651 reference bases, not 151.
        let v = vec![rec(10, "50M500D101M", 60, true, 0)];
        let t = depth_track(&v, 0, 1000);
        assert_eq!(t[10], 1);
        assert_eq!(
            t[600], 1,
            "the deleted span is still spanned by this alignment"
        );
        assert_eq!(t[9], 0);
        assert_eq!(t[661], 0);
    }

    #[test]
    fn vmr_is_one_for_poisson_like_depth_and_higher_when_clustered() {
        // Flat depth: zero variance, VMR 0. Not realistic, and that is the point — this is
        // what "too clean" looks like at the extreme.
        let flat = vec![30u32; 5000];
        assert!(depth_stats(&flat, 500).unwrap().vmr < 0.001);

        // Alternating blocks: highly overdispersed AND correlated, like real coverage.
        // Blocks must be WIDE relative to the lag. At 700 bp blocks a 500 bp lag puts most
        // pairs across a boundary and the correlation comes out NEGATIVE — which is a true
        // statement about that signal and a useless fixture for this assertion.
        let mut blocky = Vec::new();
        for i in 0..20_000 {
            blocky.push(if (i / 4000) % 2 == 0 { 10u32 } else { 50u32 });
        }
        let b = depth_stats(&blocky, 500).unwrap();
        assert!(
            b.vmr > 10.0,
            "blocky depth should be overdispersed, got {}",
            b.vmr
        );
        assert!(
            b.autocorr_500 > 0.2,
            "and correlated at 500 bp, got {}",
            b.autocorr_500
        );
    }

    #[test]
    fn autocorrelation_separates_correlated_from_independent_noise() {
        // Independent noise: near zero at any lag. Correlated: positive. A binned CNV caller
        // averages the first away and cannot average the second, so the two must not look
        // alike no matter what their VMR is.
        let mut indep = Vec::new();
        let mut x: u64 = 12345;
        for _ in 0..5000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            indep.push(20 + ((x >> 33) % 21) as u32);
        }
        let a = depth_stats(&indep, 500).unwrap();
        assert!(
            a.autocorr_500.abs() < 0.1,
            "independent noise correlated {}",
            a.autocorr_500
        );
    }

    #[test]
    fn excess_dispersion_is_comparable_across_depths_and_vmr_is_not() {
        // The same multiplicative modulation at two depths, WITH counting noise on top.
        //
        // The counting noise is not decoration. `(vmr - 1) / mean` subtracts a Poisson floor,
        // so a fixture without one is over-corrected by exactly 1/mean — which at 30x and
        // 240x is 0.033 against 0.004 and looks like depth dependence in the statistic when it
        // is really depth dependence in the fixture. A first version of this test had no
        // counting noise and failed for that reason. The statistic is built for a process that
        // has both components, so the fixture must have both.
        let factors = [0.7f64, 1.0, 1.3, 0.85, 1.15];
        let track = |mean: f64| -> Vec<u32> {
            let mut x: u64 = 4242;
            (0..20_000)
                .map(|i| {
                    let target = mean * factors[(i / 1000) % factors.len()];
                    x = x
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    // Uniform[-1,1] scaled to unit variance, times sqrt(target): variance ~
                    // target, i.e. Poisson-like counting noise around the modulated mean.
                    let u = ((x >> 33) as f64 / (1u64 << 31) as f64) * 2.0 - 1.0;
                    (target + u * 3f64.sqrt() * target.sqrt()).round().max(0.0) as u32
                })
                .collect()
        };
        let shallow = depth_stats(&track(30.0), 500).unwrap();
        let deep = depth_stats(&track(240.0), 500).unwrap();

        // VMR is NOT comparable: 8x the depth, ~8x the VMR, same underlying modulation.
        let vmr_ratio = deep.vmr / shallow.vmr;
        assert!(
            vmr_ratio > 5.0,
            "VMR should scale with depth for multiplicative noise; ratio was {vmr_ratio}"
        );

        // Excess dispersion IS comparable: the same modulation reads the same at both depths.
        let a = shallow.excess_dispersion();
        let b = deep.excess_dispersion();
        assert!(
            (a - b).abs() / a < 0.15,
            "excess dispersion should be depth-independent: {a} vs {b}"
        );
    }

    #[test]
    fn pure_poisson_like_depth_has_near_zero_excess_dispersion() {
        // A track whose only variation is counting noise must report ~0 excess, whatever its
        // depth — that is what "no coverage modulation at all" looks like, and it is what
        // eidolon currently produces.
        let mut x: u64 = 99;
        let mut t = Vec::new();
        for _ in 0..20_000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            // mean 30, variance ~30: crude but Poisson-like in the only way that matters here
            let noise = ((x >> 33) % 21) as i64 - 10;
            t.push((30 + noise).max(0) as u32);
        }
        let d = depth_stats(&t, 500).unwrap();
        assert!(
            d.excess_dispersion().abs() < 0.02,
            "counting noise alone should leave ~0 excess, got {}",
            d.excess_dispersion()
        );
    }

    // ── the release goal itself ────────────────────────────────────────────────

    /// THE POINT OF THE PANEL. Two datasets that differ only in their artifact load must be
    /// distinguishable, and two identical ones must not be. A panel that cannot tell those
    /// apart is measuring nothing, and would have reported eidolon as realistic.
    ///
    /// The "clean" set is modelled on what eidolon 3.2.0 actually produces: every read a
    /// perfect match, every pair proper, no clipping anywhere. Measured on real chr22 the
    /// equivalents are 654.5 candidates/Mb, 3.15% improper and 0.67% clipped.
    #[test]
    fn the_panel_separates_clean_data_from_artifacted_data() {
        let clean: Vec<_> = (0..1000)
            .map(|i| rec(i * 100, "151M", 60, true, 400))
            .collect();

        let mut dirty: Vec<_> = (0..940)
            .map(|i| rec(i * 100, "151M", 60, true, 400))
            .collect();
        dirty.extend((0..30).map(|i| rec(50_000 + (i / 5) * 3000, "40S111M", 60, true, 400)));
        dirty.extend((0..20).map(|_| rec(70_000, "151M", 0, false, 400)));
        dirty.extend((0..10).map(|_| rec(80_000, "151M", 60, false, 9000)));

        let m = |v: &[AlnRecord]| RegionMetrics {
            reads: v.len(),
            span_bp: 100_000,
            candidate_breakpoints: candidate_breakpoints(v, 20, 3),
            improper_pairs: v.iter().filter(|r| !r.is_proper_pair()).count(),
            clipped_reads: v
                .iter()
                .filter(|r| r.leading_clip() >= 20 || r.trailing_clip() >= 20)
                .count(),
            mapq0_reads: v.iter().filter(|r| r.mapq == 0).count(),
            insert: insert_stats(v, 2000),
            depth: None,
        };

        let c = m(&clean);
        let d = m(&dirty);

        // Clean data must look clean — this is the eidolon side of the gap.
        assert_eq!(c.candidate_breakpoints, 0);
        assert_eq!(c.improper_pair_rate(), 0.0);
        assert_eq!(c.clip_rate(), 0.0);

        // Artifacted data must be visibly different on every axis.
        assert!(
            d.candidates_per_mb() > 0.0,
            "no candidates found in artifacted data"
        );
        assert!(
            d.improper_pair_rate() > 0.02,
            "improper rate {}",
            d.improper_pair_rate()
        );
        assert!(d.clip_rate() > 0.02, "clip rate {}", d.clip_rate());
        assert!(d.mapq0_rate() > 0.01, "mapq0 rate {}", d.mapq0_rate());
    }

    /// MUST NOT FIRE: the same input twice reports no gap. This is the real-vs-real self-test
    /// in unit form — if the metrics disagreed with themselves, every gap the panel reported
    /// would be its own noise.
    #[test]
    fn identical_inputs_report_an_identical_profile() {
        let mut v: Vec<_> = (0..500)
            .map(|i| rec(i * 200, "151M", 60, true, 400))
            .collect();
        v.extend((0..12).map(|_| rec(9_000, "40S111M", 0, false, 400)));

        let profile = |x: &[AlnRecord]| {
            (
                candidate_breakpoints(x, 20, 3),
                x.iter().filter(|r| !r.is_proper_pair()).count(),
                insert_stats(x, 2000).map(|s| (s.n, s.mean.to_bits())),
                depth_stats(&depth_track(x, 0, 100_000), 500)
                    .map(|d| (d.vmr.to_bits(), d.mean.to_bits())),
            )
        };
        assert_eq!(profile(&v), profile(&v.clone()));
    }
}
