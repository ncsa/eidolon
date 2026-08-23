//! Coordinate mapping between an altered haplotype and the reference it derives from.
//!
//! Lives in `eidolon-core` rather than beside the fragment sampler because BOTH
//! sides of read generation need it: the sampler (in the `eidolon` binary) places
//! fragments in haplotype coordinates, and the writer (`file_tools::fastq_tools`,
//! here in core) has to project those back to reference coordinates for BAM
//! records. Crate dependencies run `eidolon` -> `eidolon-core` only, so a shared
//! primitive has to live on this side (#516).
//!
//! **Why this handles MANY insertions, not one.** The first version of the rework
//! described a single insertion, and the sampler fell back to the pre-#516
//! head-only behaviour whenever a sub-region held more than one. Sub-regions are
//! large -- a whole contig, split only by coverage multipliers -- so that fallback
//! caught every multi-insertion case, including the obvious way to exercise this
//! feature: plant a range of sizes at once via `input_vcf`. Measured on three
//! insertions (200/600/1200bp) sharing a contig: head 10/14/15, middle 0/0/0,
//! tail 0/0/0. A map over a sorted set of insertions removes the fallback and the
//! restriction together.

use crate::structs::nucleotides::Nucleotide;

/// One literal insertion spliced into the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
struct InsertionEntry {
    /// Zero-based reference base represented by the VCF REF allele. Novel bases
    /// are inserted immediately AFTER it, so reference positions up to and
    /// including the anchor are unshifted.
    anchor: usize,
    bases: Vec<Nucleotide>,
    /// Haplotype offset at which this insertion's novel bases begin, i.e.
    /// `anchor + 1 + (inserted length of everything anchored before it)`.
    hap_start: usize,
}

/// Coordinate map between a reference contig and an altered haplotype carrying
/// one or more literal insertions.
///
/// Reference coordinates shift by the total inserted length anchored strictly
/// before them; haplotype positions falling inside inserted sequence have no
/// reference coordinate at all. The map owns the inserted bases so a caller
/// materializing an interval does not have to track which sequence belongs to
/// which anchor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InsertionCoordinateMap {
    reference_len: usize,
    /// Sorted by `anchor`, non-overlapping.
    insertions: Vec<InsertionEntry>,
    total_inserted: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HaplotypeSegment {
    Reference {
        hap_start: usize,
        hap_end: usize,
        ref_start: usize,
    },
    Insertion {
        hap_start: usize,
        hap_end: usize,
        /// Offset into THIS insertion's novel bases.
        insertion_start: usize,
        /// Which insertion, for callers resolving bases themselves.
        index: usize,
    },
}

impl InsertionCoordinateMap {
    /// Build a map over any number of insertions. Returns `None` if the
    /// reference is empty, if an anchor lies outside it, if any insertion is
    /// empty, or if two insertions share an anchor — all of which would make the
    /// projections ambiguous rather than merely unusual.
    pub fn new(
        reference_len: usize,
        insertions: impl IntoIterator<Item = (usize, Vec<Nucleotide>)>,
    ) -> Option<Self> {
        if reference_len == 0 {
            return None;
        }
        let mut entries: Vec<(usize, Vec<Nucleotide>)> = insertions.into_iter().collect();
        if entries.is_empty() {
            return None;
        }
        entries.sort_by_key(|(anchor, _)| *anchor);
        let mut built: Vec<InsertionEntry> = Vec::with_capacity(entries.len());
        let mut cumulative = 0usize;
        let mut previous_anchor: Option<usize> = None;
        for (anchor, bases) in entries {
            if anchor >= reference_len || bases.is_empty() {
                return None;
            }
            if previous_anchor == Some(anchor) {
                return None;
            }
            previous_anchor = Some(anchor);
            let hap_start = anchor + 1 + cumulative;
            cumulative += bases.len();
            built.push(InsertionEntry {
                anchor,
                bases,
                hap_start,
            });
        }
        Some(Self {
            reference_len,
            insertions: built,
            total_inserted: cumulative,
        })
    }

    /// Convenience for the single-insertion case.
    pub fn single(reference_len: usize, anchor: usize, bases: Vec<Nucleotide>) -> Option<Self> {
        Self::new(reference_len, [(anchor, bases)])
    }

    pub fn haplotype_len(&self) -> usize {
        self.reference_len + self.total_inserted
    }

    /// Reference anchors, ascending. The writer uses these to skip applying an
    /// insertion inline when its allele is already decided by haplotype sampling.
    pub fn anchors(&self) -> impl Iterator<Item = usize> + '_ {
        self.insertions.iter().map(|i| i.anchor)
    }

    /// Total inserted length anchored strictly before `reference_pos`.
    fn shift_at(&self, reference_pos: usize) -> usize {
        // Insertions are sorted by anchor, so this is a prefix boundary: find the
        // first anchor that is NOT before `reference_pos`. Binary search keeps
        // this O(log n) per lookup, which matters because the writer calls it per
        // read window rather than per fragment.
        let idx = self
            .insertions
            .partition_point(|i| i.anchor < reference_pos);
        self.insertions[..idx].iter().map(|i| i.bases.len()).sum()
    }

    /// Index of the insertion containing this haplotype position, if any.
    fn insertion_at(&self, haplotype_pos: usize) -> Option<usize> {
        let idx = self
            .insertions
            .partition_point(|i| i.hap_start <= haplotype_pos);
        if idx == 0 {
            return None;
        }
        let candidate = &self.insertions[idx - 1];
        (haplotype_pos < candidate.hap_start + candidate.bases.len()).then_some(idx - 1)
    }

    pub fn reference_base_to_haplotype(&self, reference_pos: usize) -> Option<usize> {
        (reference_pos < self.reference_len).then(|| reference_pos + self.shift_at(reference_pos))
    }

    /// Reference position for a haplotype position, or `None` when it falls
    /// inside inserted sequence — which genuinely has no reference coordinate.
    pub fn haplotype_base_to_reference(&self, haplotype_pos: usize) -> Option<usize> {
        if haplotype_pos >= self.haplotype_len() || self.insertion_at(haplotype_pos).is_some() {
            return None;
        }
        Some(self.reference_floor(haplotype_pos))
    }

    /// Reference position at or after `haplotype_pos`.
    ///
    /// Unlike [`Self::haplotype_base_to_reference`] this is **total**: a position
    /// inside inserted sequence yields the first reference base following that
    /// insertion. That makes it usable for projecting a half-open *window* back
    /// to reference coordinates, where a `None` would force the caller to give up
    /// and scan every variant instead of binary-searching a range.
    pub fn reference_floor(&self, haplotype_pos: usize) -> usize {
        match self.insertion_at(haplotype_pos) {
            // Inside an insertion: the next reference base is the one after its anchor.
            Some(i) => self.insertions[i].anchor + 1,
            None => {
                let consumed: usize = self
                    .insertions
                    .iter()
                    .take_while(|i| i.hap_start + i.bases.len() <= haplotype_pos)
                    .map(|i| i.bases.len())
                    .sum();
                haplotype_pos - consumed
            }
        }
    }

    /// Split a haplotype interval into alternating reference and insertion runs.
    pub fn segments_for(&self, start: usize, end: usize) -> Option<Vec<HaplotypeSegment>> {
        if start >= end || end > self.haplotype_len() {
            return None;
        }
        let mut segments = Vec::new();
        let mut cursor = start;
        while cursor < end {
            match self.insertion_at(cursor) {
                Some(i) => {
                    let entry = &self.insertions[i];
                    let seg_end = end.min(entry.hap_start + entry.bases.len());
                    segments.push(HaplotypeSegment::Insertion {
                        hap_start: cursor,
                        hap_end: seg_end,
                        insertion_start: cursor - entry.hap_start,
                        index: i,
                    });
                    cursor = seg_end;
                }
                None => {
                    // Reference run up to the next insertion, or the interval end.
                    let next = self
                        .insertions
                        .iter()
                        .map(|i| i.hap_start)
                        .find(|&s| s > cursor)
                        .unwrap_or(usize::MAX);
                    let seg_end = end.min(next);
                    segments.push(HaplotypeSegment::Reference {
                        hap_start: cursor,
                        hap_end: seg_end,
                        ref_start: self.reference_floor(cursor),
                    });
                    cursor = seg_end;
                }
            }
        }
        Some(segments)
    }

    /// Materialize a haplotype interval. The returned segments retain enough
    /// provenance for the caller to build reference-anchored BAM CIGAR operations.
    pub fn materialize_interval(
        &self,
        reference: &[Nucleotide],
        start: usize,
        end: usize,
    ) -> Option<(Vec<Nucleotide>, Vec<HaplotypeSegment>)> {
        if reference.len() != self.reference_len {
            return None;
        }
        let segments = self.segments_for(start, end)?;
        let mut sequence = Vec::with_capacity(end - start);
        for segment in &segments {
            match *segment {
                HaplotypeSegment::Reference {
                    ref_start,
                    hap_start,
                    hap_end,
                } => {
                    let len = hap_end - hap_start;
                    sequence.extend_from_slice(reference.get(ref_start..ref_start + len)?);
                }
                HaplotypeSegment::Insertion {
                    insertion_start,
                    hap_start,
                    hap_end,
                    index,
                } => {
                    let len = hap_end - hap_start;
                    let bases = &self.insertions.get(index)?.bases;
                    sequence.extend_from_slice(bases.get(insertion_start..insertion_start + len)?);
                }
            }
        }
        Some((sequence, segments))
    }

    /// Baseline CIGAR operation for every materialized base. Deliberately
    /// excludes sequencing-error operations; the read generator layers those on
    /// while keeping the insertion `I` operations supplied here.
    pub fn cigar_ops_for_segments(segments: &[HaplotypeSegment]) -> Vec<char> {
        let mut cigar = Vec::new();
        for segment in segments {
            let (len, op) = match *segment {
                HaplotypeSegment::Reference {
                    hap_start, hap_end, ..
                } => (hap_end - hap_start, 'M'),
                HaplotypeSegment::Insertion {
                    hap_start, hap_end, ..
                } => (hap_end - hap_start, 'I'),
            };
            cigar.extend(std::iter::repeat_n(op, len));
        }
        cigar
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::{A, C, G};

    /// n=1 — the original single-insertion case, preserved through the
    /// generalization so the simple geometry stays pinned.
    #[test]
    fn single_insertion_tracks_anchor_interior_and_tail() {
        let map = InsertionCoordinateMap::single(1_000, 500, vec![C; 600]).unwrap();
        assert_eq!(map.haplotype_len(), 1_600);
        // Reference through the anchor is unshifted; past it, shifted by the insert.
        assert_eq!(map.reference_base_to_haplotype(500), Some(500));
        assert_eq!(map.reference_base_to_haplotype(501), Some(1_101));
        assert_eq!(map.haplotype_base_to_reference(500), Some(500));
        // Inside the inserted sequence there is no reference coordinate at all.
        assert_eq!(map.haplotype_base_to_reference(501), None);
        assert_eq!(map.haplotype_base_to_reference(1_101), Some(501));
        // ...but the total projection still yields the next reference base.
        assert_eq!(map.reference_floor(501), 501);

        assert_eq!(
            map.segments_for(490, 520),
            Some(vec![
                HaplotypeSegment::Reference {
                    hap_start: 490,
                    hap_end: 501,
                    ref_start: 490,
                },
                HaplotypeSegment::Insertion {
                    hap_start: 501,
                    hap_end: 520,
                    insertion_start: 0,
                    index: 0,
                },
            ])
        );

        let reference = vec![A; 1_000];
        let (sequence, segments) = map.materialize_interval(&reference, 1_090, 1_120).unwrap();
        assert_eq!(
            sequence,
            vec![C; 11]
                .into_iter()
                .chain(vec![A; 19])
                .collect::<Vec<_>>()
        );
        assert_eq!(
            InsertionCoordinateMap::cigar_ops_for_segments(&segments),
            vec!['I'; 11]
                .into_iter()
                .chain(vec!['M'; 19])
                .collect::<Vec<_>>()
        );
    }

    /// THE CASE THE SINGLE-INSERTION MAP COULD NOT EXPRESS. Two insertions on one
    /// contig: coordinates past the second must shift by BOTH lengths, and an
    /// interval spanning both must materialize reference/insert/reference/insert.
    #[test]
    fn multiple_insertions_compose_their_shifts() {
        // anchors 100 (10 x C) and 200 (5 x G) on a 1000bp reference.
        let map =
            InsertionCoordinateMap::new(1_000, [(100, vec![C; 10]), (200, vec![G; 5])]).unwrap();
        assert_eq!(map.haplotype_len(), 1_015);
        assert_eq!(map.anchors().collect::<Vec<_>>(), vec![100, 200]);

        // Before either: unshifted.
        assert_eq!(map.reference_base_to_haplotype(50), Some(50));
        // Between them: shifted by the first only.
        assert_eq!(map.reference_base_to_haplotype(150), Some(160));
        // After both: shifted by both.
        assert_eq!(map.reference_base_to_haplotype(300), Some(315));

        // Round trip, and the insertion interiors report no reference coordinate.
        assert_eq!(map.haplotype_base_to_reference(160), Some(150));
        assert_eq!(map.haplotype_base_to_reference(315), Some(300));
        // Insert #0 occupies haplotype [101, 111): anchor 100, +1, no preceding shift.
        assert_eq!(map.haplotype_base_to_reference(101), None);
        assert_eq!(map.haplotype_base_to_reference(110), None);
        // Insert #1 occupies haplotype [211, 216) — anchor 200, +1, SHIFTED by the
        // 10 bases of insert #0. Getting this wrong by forgetting the shift is
        // precisely the arithmetic this map exists to centralize, so it is asserted
        // on both sides of the boundary rather than at one convenient point.
        assert_eq!(map.haplotype_base_to_reference(210), Some(200)); // the anchor itself
        assert_eq!(map.haplotype_base_to_reference(211), None); // first inserted base
        assert_eq!(map.haplotype_base_to_reference(215), None); // last inserted base
        assert_eq!(map.haplotype_base_to_reference(216), Some(201)); // back to reference
        // reference_floor stays total across both: inside an insertion it yields the
        // next reference base, outside it agrees with the partial projection.
        assert_eq!(map.reference_floor(101), 101);
        assert_eq!(map.reference_floor(211), 201);
        assert_eq!(map.reference_floor(206), 196); // between the two, plain reference

        // An interval spanning both insertions alternates correctly.
        let reference = vec![A; 1_000];
        let (sequence, segments) = map.materialize_interval(&reference, 99, 215).unwrap();
        // 2 ref (99,100) + 10 ins + 100 ref (101..201) + 4 ins = 116 bases.
        assert_eq!(sequence.len(), 116);
        assert_eq!(sequence.iter().filter(|&&b| b == C).count(), 10);
        assert_eq!(sequence.iter().filter(|&&b| b == G).count(), 4);
        assert_eq!(segments.len(), 4);
        let cigar = InsertionCoordinateMap::cigar_ops_for_segments(&segments);
        assert_eq!(cigar.len(), 116);
        assert_eq!(cigar.iter().filter(|&&c| c == 'I').count(), 14);
    }

    /// MUST NOT FIRE: ambiguous or impossible inputs are refused rather than
    /// silently producing a map whose projections disagree with each other.
    #[test]
    fn degenerate_inputs_are_refused() {
        // Two insertions at one anchor — ordering between them is undefined.
        assert!(
            InsertionCoordinateMap::new(1_000, [(100, vec![C; 5]), (100, vec![G; 5])]).is_none()
        );
        // Anchor outside the reference.
        assert!(InsertionCoordinateMap::single(100, 100, vec![C; 5]).is_none());
        // Empty insertion carries no novel sequence.
        assert!(InsertionCoordinateMap::single(1_000, 100, vec![]).is_none());
        // No insertions at all — callers should not build a haplotype for this.
        assert!(InsertionCoordinateMap::new(1_000, []).is_none());
        // Reference length must match at materialization time.
        let map = InsertionCoordinateMap::single(1_000, 100, vec![C; 5]).unwrap();
        assert!(map.materialize_interval(&vec![A; 999], 0, 10).is_none());
    }

    /// Insertions supplied out of order must behave identically to sorted input;
    /// the sampler has no reason to guarantee ordering.
    #[test]
    fn insertion_order_does_not_matter() {
        let sorted =
            InsertionCoordinateMap::new(1_000, [(100, vec![C; 10]), (200, vec![G; 5])]).unwrap();
        let shuffled =
            InsertionCoordinateMap::new(1_000, [(200, vec![G; 5]), (100, vec![C; 10])]).unwrap();
        assert_eq!(sorted, shuffled);
    }
}
