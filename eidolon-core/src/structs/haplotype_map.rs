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

/// One edit applied to the reference to produce the haplotype.
///
/// Both kinds share the VCF anchor convention: the anchor is the reference base
/// named by REF, it is itself unaffected, and the edit applies immediately AFTER
/// it. So an insertion splices novel bases in at `anchor + 1`, and a deletion
/// removes reference bases `[anchor + 1, anchor + 1 + len)`.
#[derive(Debug, Clone, PartialEq, Eq)]
enum EditKind {
    Insertion(Vec<Nucleotide>),
    /// Count of reference bases removed after the anchor.
    Deletion(usize),
}

/// One edit, as supplied by the caller.
#[derive(Debug, Clone, PartialEq, Eq)]
struct EditEntry {
    /// Zero-based reference base represented by the VCF REF allele.
    anchor: usize,
    kind: EditKind,
}

/// A contiguous run of the haplotype, or a deletion sitting between two runs.
///
/// Precomputed once at construction so every projection is a search over blocks
/// rather than arithmetic that has to special-case each edit kind. Deletions and
/// insertions are genuinely different — one consumes reference without producing
/// haplotype, the other the reverse — and two separate arithmetic paths for them
/// is exactly the "two paths that disagree" shape this file already carries a
/// warning about.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Block {
    /// Haplotype `[hap_start, hap_end)` maps to reference starting at `ref_start`.
    Reference {
        hap_start: usize,
        hap_end: usize,
        ref_start: usize,
    },
    /// Haplotype `[hap_start, hap_end)` is novel sequence from insertion `index`.
    Insertion {
        hap_start: usize,
        hap_end: usize,
        index: usize,
    },
    /// Zero width in haplotype coordinates: reference `[ref_start, ref_start + ref_len)`
    /// is absent from the haplotype. `hap_at` is the haplotype position it sits between.
    Deletion {
        hap_at: usize,
        ref_start: usize,
        ref_len: usize,
    },
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
    /// Sorted by anchor, non-overlapping: `(anchor, novel bases)`. Insertions only,
    /// in the order their `index` refers to, so `HaplotypeSegment::Insertion::index`
    /// stays stable.
    insertions: Vec<(usize, Vec<Nucleotide>)>,
    /// Every block of the haplotype, ascending, covering it exactly once.
    blocks: Vec<Block>,
    haplotype_len: usize,
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
    /// Reference bases absent from the haplotype, sitting between the segments on
    /// either side. Consumes NO haplotype (and therefore no read) sequence, which
    /// is why it cannot be represented in a per-base op vector — see
    /// [`InsertionCoordinateMap::deletion_runs_for_segments`].
    Deletion {
        hap_at: usize,
        ref_start: usize,
        ref_len: usize,
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
        Self::with_edits(
            reference_len,
            insertions
                .into_iter()
                .map(|(anchor, bases)| (anchor, EditKind::Insertion(bases))),
            std::iter::empty(),
        )
    }

    /// Build a map over insertions AND deletions.
    ///
    /// A deletion is `(anchor, len)`: reference bases `[anchor + 1, anchor + 1 + len)`
    /// are absent from the haplotype, matching the VCF convention where POS names the
    /// unaffected base before the event. This is what lets a fragment be *placed* in a
    /// coordinate space where the deleted bases do not exist — the per-base path could
    /// only ever skip them while rendering a read, which cannot remove coverage from a
    /// span fragments are still placed across (#590).
    pub fn with_deletions(
        reference_len: usize,
        insertions: impl IntoIterator<Item = (usize, Vec<Nucleotide>)>,
        deletions: impl IntoIterator<Item = (usize, usize)>,
    ) -> Option<Self> {
        Self::with_edits(
            reference_len,
            insertions
                .into_iter()
                .map(|(anchor, bases)| (anchor, EditKind::Insertion(bases))),
            deletions
                .into_iter()
                .map(|(anchor, len)| (anchor, EditKind::Deletion(len))),
        )
    }

    fn with_edits(
        reference_len: usize,
        insertions: impl IntoIterator<Item = (usize, EditKind)>,
        deletions: impl IntoIterator<Item = (usize, EditKind)>,
    ) -> Option<Self> {
        if reference_len == 0 {
            return None;
        }
        let mut entries: Vec<EditEntry> = insertions
            .into_iter()
            .chain(deletions)
            .map(|(anchor, kind)| EditEntry { anchor, kind })
            .collect();
        if entries.is_empty() {
            return None;
        }
        // Stable by anchor so an insertion and a deletion sharing one is still rejected
        // below rather than silently ordered.
        entries.sort_by_key(|e| e.anchor);

        let mut insertions_out: Vec<(usize, Vec<Nucleotide>)> = Vec::new();
        let mut blocks: Vec<Block> = Vec::new();
        let mut ref_cursor = 0usize;
        let mut hap_cursor = 0usize;
        let mut previous_anchor: Option<usize> = None;

        for entry in entries {
            // Reject anything that would make a projection ambiguous rather than merely
            // unusual: an anchor outside the reference, a degenerate edit, two edits on
            // one anchor, or an edit whose anchor was already consumed by a preceding
            // deletion (which would mean overlapping events).
            if entry.anchor >= reference_len {
                return None;
            }
            if previous_anchor == Some(entry.anchor) {
                return None;
            }
            previous_anchor = Some(entry.anchor);
            if entry.anchor + 1 < ref_cursor {
                return None;
            }

            // Emit the reference run up to and including this anchor.
            let run_end = entry.anchor + 1;
            if run_end > ref_cursor {
                blocks.push(Block::Reference {
                    hap_start: hap_cursor,
                    hap_end: hap_cursor + (run_end - ref_cursor),
                    ref_start: ref_cursor,
                });
                hap_cursor += run_end - ref_cursor;
                ref_cursor = run_end;
            }

            match entry.kind {
                EditKind::Insertion(bases) => {
                    if bases.is_empty() {
                        return None;
                    }
                    let len = bases.len();
                    blocks.push(Block::Insertion {
                        hap_start: hap_cursor,
                        hap_end: hap_cursor + len,
                        index: insertions_out.len(),
                    });
                    insertions_out.push((entry.anchor, bases));
                    hap_cursor += len;
                }
                EditKind::Deletion(len) => {
                    if len == 0 || ref_cursor + len > reference_len {
                        return None;
                    }
                    blocks.push(Block::Deletion {
                        hap_at: hap_cursor,
                        ref_start: ref_cursor,
                        ref_len: len,
                    });
                    ref_cursor += len;
                }
            }
        }

        if ref_cursor < reference_len {
            blocks.push(Block::Reference {
                hap_start: hap_cursor,
                hap_end: hap_cursor + (reference_len - ref_cursor),
                ref_start: ref_cursor,
            });
            hap_cursor += reference_len - ref_cursor;
        }
        // A haplotype with no bases left is not a coordinate space anything can be
        // placed in; refuse rather than hand back an empty one.
        if hap_cursor == 0 {
            return None;
        }

        Some(Self {
            reference_len,
            insertions: insertions_out,
            blocks,
            haplotype_len: hap_cursor,
        })
    }

    /// Convenience for the single-insertion case.
    pub fn single(reference_len: usize, anchor: usize, bases: Vec<Nucleotide>) -> Option<Self> {
        Self::new(reference_len, [(anchor, bases)])
    }

    pub fn haplotype_len(&self) -> usize {
        self.haplotype_len
    }

    /// Reference anchors, ascending. The writer uses these to skip applying an
    /// insertion inline when its allele is already decided by haplotype sampling.
    pub fn anchors(&self) -> impl Iterator<Item = usize> + '_ {
        self.insertions.iter().map(|(anchor, _)| *anchor)
    }

    /// Reference anchors of the deletions, ascending. The writer uses these the same
    /// way it uses [`Self::anchors`]: a deletion already expressed by the haplotype
    /// must not ALSO be applied inline while rendering the read, or it is applied
    /// twice.
    pub fn deletion_anchors(&self) -> impl Iterator<Item = usize> + '_ {
        self.blocks.iter().filter_map(|b| match *b {
            Block::Deletion { ref_start, .. } => Some(ref_start.saturating_sub(1)),
            _ => None,
        })
    }

    /// The block containing this haplotype position. `None` past the end.
    ///
    /// Blocks are ascending and cover the haplotype exactly once, so this is a
    /// prefix boundary — binary search keeps it O(log n), which matters because the
    /// writer calls it per read window rather than per fragment. Zero-width
    /// deletion blocks are skipped: they contain no haplotype position.
    fn block_at(&self, haplotype_pos: usize) -> Option<usize> {
        if haplotype_pos >= self.haplotype_len {
            return None;
        }
        let idx = self
            .blocks
            .partition_point(|b| Self::block_hap_start(b) <= haplotype_pos);
        // Walk back over any zero-width deletions to the block that really holds it.
        (0..idx)
            .rev()
            .find(|&i| Self::block_hap_end(&self.blocks[i]) > haplotype_pos)
    }

    fn block_hap_start(b: &Block) -> usize {
        match *b {
            Block::Reference { hap_start, .. } | Block::Insertion { hap_start, .. } => hap_start,
            Block::Deletion { hap_at, .. } => hap_at,
        }
    }

    fn block_hap_end(b: &Block) -> usize {
        match *b {
            Block::Reference { hap_end, .. } | Block::Insertion { hap_end, .. } => hap_end,
            Block::Deletion { hap_at, .. } => hap_at,
        }
    }

    /// Haplotype position of a reference base, or `None` when that base is absent
    /// from this haplotype — which is exactly what a deletion means, and the reason
    /// no fragment can be placed there.
    pub fn reference_base_to_haplotype(&self, reference_pos: usize) -> Option<usize> {
        if reference_pos >= self.reference_len {
            return None;
        }
        for block in &self.blocks {
            if let Block::Reference {
                hap_start,
                hap_end,
                ref_start,
            } = *block
            {
                let len = hap_end - hap_start;
                if reference_pos >= ref_start && reference_pos < ref_start + len {
                    return Some(hap_start + (reference_pos - ref_start));
                }
            }
        }
        None
    }

    /// Reference position for a haplotype position, or `None` when it falls
    /// inside inserted sequence — which genuinely has no reference coordinate.
    pub fn haplotype_base_to_reference(&self, haplotype_pos: usize) -> Option<usize> {
        match self.block_at(haplotype_pos).map(|i| &self.blocks[i]) {
            Some(Block::Reference {
                hap_start,
                ref_start,
                ..
            }) => Some(ref_start + (haplotype_pos - hap_start)),
            _ => None,
        }
    }

    /// Reference position at or after `haplotype_pos`.
    ///
    /// Unlike [`Self::haplotype_base_to_reference`] this is **total**: a position
    /// inside inserted sequence yields the first reference base following that
    /// insertion. That makes it usable for projecting a half-open *window* back
    /// to reference coordinates, where a `None` would force the caller to give up
    /// and scan every variant instead of binary-searching a range.
    pub fn reference_floor(&self, haplotype_pos: usize) -> usize {
        match self.block_at(haplotype_pos).map(|i| (i, &self.blocks[i])) {
            Some((
                _,
                Block::Reference {
                    hap_start,
                    ref_start,
                    ..
                },
            )) => ref_start + (haplotype_pos - hap_start),
            // Inside an insertion: the next reference base is the start of the next
            // reference block, i.e. the base after this insertion's anchor.
            Some((idx, Block::Insertion { .. })) => self.blocks[idx..]
                .iter()
                .find_map(|b| match *b {
                    Block::Reference { ref_start, .. } => Some(ref_start),
                    _ => None,
                })
                .unwrap_or(self.reference_len),
            _ => self.reference_len,
        }
    }

    /// Split a haplotype interval into alternating reference and insertion runs.
    pub fn segments_for(&self, start: usize, end: usize) -> Option<Vec<HaplotypeSegment>> {
        if start >= end || end > self.haplotype_len() {
            return None;
        }
        let mut segments = Vec::new();
        for block in &self.blocks {
            match *block {
                Block::Reference {
                    hap_start,
                    hap_end,
                    ref_start,
                } => {
                    let lo = hap_start.max(start);
                    let hi = hap_end.min(end);
                    if lo < hi {
                        segments.push(HaplotypeSegment::Reference {
                            hap_start: lo,
                            hap_end: hi,
                            ref_start: ref_start + (lo - hap_start),
                        });
                    }
                }
                Block::Insertion {
                    hap_start,
                    hap_end,
                    index,
                } => {
                    let lo = hap_start.max(start);
                    let hi = hap_end.min(end);
                    if lo < hi {
                        segments.push(HaplotypeSegment::Insertion {
                            hap_start: lo,
                            hap_end: hi,
                            insertion_start: lo - hap_start,
                            index,
                        });
                    }
                }
                Block::Deletion {
                    hap_at,
                    ref_start,
                    ref_len,
                } => {
                    // Zero width, so "intersects the window" means the window has at
                    // least one base on BOTH sides: a deletion at the very edge of an
                    // interval is not spanned by it and must not contribute a D op.
                    if hap_at > start && hap_at < end {
                        segments.push(HaplotypeSegment::Deletion {
                            hap_at,
                            ref_start,
                            ref_len,
                        });
                    }
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
                    let bases = &self.insertions.get(index)?.1;
                    sequence.extend_from_slice(bases.get(insertion_start..insertion_start + len)?);
                }
                // Contributes no bases by definition — that IS the deletion. The
                // reference bases on either side end up adjacent in the fragment,
                // which is what makes a junction-spanning read correct without the
                // per-base path having to skip anything.
                HaplotypeSegment::Deletion { .. } => {}
            }
        }
        Some((sequence, segments))
    }

    /// Deletions inside a materialized interval, as `(query offset, reference bases
    /// removed)`, ascending by offset.
    ///
    /// Deliberately separate from [`Self::cigar_ops_for_segments`]: that returns one
    /// op per materialized base, and a `D` consumes reference WITHOUT consuming
    /// query, so it cannot be expressed there at all. The read generator emits these
    /// as `D` runs after writing the base at each offset — the same shape the
    /// sequencing-error deletion path already uses.
    pub fn deletion_runs_for_segments(segments: &[HaplotypeSegment]) -> Vec<(usize, usize)> {
        let mut runs = Vec::new();
        let mut query = 0usize;
        let base = segments
            .first()
            .map(|s| match *s {
                HaplotypeSegment::Reference { hap_start, .. }
                | HaplotypeSegment::Insertion { hap_start, .. } => hap_start,
                HaplotypeSegment::Deletion { hap_at, .. } => hap_at,
            })
            .unwrap_or(0);
        for segment in segments {
            match *segment {
                HaplotypeSegment::Reference {
                    hap_start, hap_end, ..
                }
                | HaplotypeSegment::Insertion {
                    hap_start, hap_end, ..
                } => {
                    query = (hap_end - base).max(query.max(hap_start - base));
                }
                HaplotypeSegment::Deletion { ref_len, .. } => {
                    runs.push((query, ref_len));
                }
            }
        }
        runs
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
                // Consumes reference, not query: see `deletion_runs_for_segments`.
                HaplotypeSegment::Deletion { .. } => continue,
            };
            cigar.extend(std::iter::repeat_n(op, len));
        }
        cigar
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::{A, C, G, T};

    /// A reference whose base at index `i` is determined by `i % 4`, so a
    /// materialized interval can be checked against the exact reference positions it
    /// claims to come from rather than against a uniform filler that would match
    /// anywhere.
    fn patterned_reference(len: usize) -> Vec<Nucleotide> {
        (0..len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect()
    }

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

    // ── deletions (#590) ────────────────────────────────────────────────────────
    //
    // The point of putting deletions in the map is that the deleted bases DO NOT
    // EXIST in the coordinate space fragments are placed in. That is what removes
    // coverage; the per-base path could only ever skip them while rendering a read,
    // which cannot remove reads from a span fragments are still placed across.

    /// Known answer: 1000 bp reference, 200 bases deleted after anchor 499, so
    /// reference [500, 700) is gone and the haplotype is 800 long.
    #[test]
    fn a_deletion_removes_its_bases_from_the_coordinate_space() {
        let map = InsertionCoordinateMap::with_deletions(1_000, [], [(499, 200)]).unwrap();
        assert_eq!(map.haplotype_len(), 800);

        // The anchor itself is untouched — VCF POS names the base BEFORE the event.
        assert_eq!(map.reference_base_to_haplotype(499), Some(499));
        // Every deleted base has no haplotype position at all. This is the property
        // the fragment placer relies on.
        assert_eq!(map.reference_base_to_haplotype(500), None);
        assert_eq!(map.reference_base_to_haplotype(699), None);
        // The first surviving base lands immediately after the anchor.
        assert_eq!(map.reference_base_to_haplotype(700), Some(500));
        assert_eq!(map.reference_base_to_haplotype(999), Some(799));

        // ...and back again.
        assert_eq!(map.haplotype_base_to_reference(499), Some(499));
        assert_eq!(map.haplotype_base_to_reference(500), Some(700));
        assert_eq!(map.haplotype_base_to_reference(799), Some(999));
        assert_eq!(map.haplotype_base_to_reference(800), None);
    }

    /// A fragment spanning the junction must materialize the two flanks ADJACENT,
    /// and report the deletion so the read can carry a `D`.
    #[test]
    fn a_junction_spanning_interval_materializes_the_flanks_adjacent() {
        let map = InsertionCoordinateMap::with_deletions(1_000, [], [(499, 200)]).unwrap();
        let reference = patterned_reference(1_000);

        let (sequence, segments) = map.materialize_interval(&reference, 490, 510).unwrap();

        // 10 bases from reference 490..500, then 10 from 700..710 — NOT 490..510.
        let expected: Vec<Nucleotide> = reference[490..500]
            .iter()
            .chain(reference[700..710].iter())
            .copied()
            .collect();
        assert_eq!(sequence, expected, "the deleted bases must not be emitted");
        assert_eq!(sequence.len(), 20);

        assert_eq!(
            segments,
            vec![
                HaplotypeSegment::Reference {
                    hap_start: 490,
                    hap_end: 500,
                    ref_start: 490,
                },
                HaplotypeSegment::Deletion {
                    hap_at: 500,
                    ref_start: 500,
                    ref_len: 200,
                },
                HaplotypeSegment::Reference {
                    hap_start: 500,
                    hap_end: 510,
                    ref_start: 700,
                },
            ]
        );

        // A D consumes reference, not query, so it is reported separately from the
        // per-base op vector rather than being silently dropped from it.
        assert_eq!(
            InsertionCoordinateMap::cigar_ops_for_segments(&segments),
            vec!['M'; 20]
        );
        assert_eq!(
            InsertionCoordinateMap::deletion_runs_for_segments(&segments),
            vec![(10, 200)],
            "the D belongs after the 10th query base, and removes 200 reference bases"
        );
    }

    /// MUST NOT FIRE: a deletion at the very edge of an interval is not spanned by
    /// it. Emitting a `D` there would put a deletion in a read that does not cross
    /// the junction, which is a wrong alignment rather than a harmless extra op.
    #[test]
    fn a_deletion_flush_against_an_interval_edge_is_not_spanned() {
        let map = InsertionCoordinateMap::with_deletions(1_000, [], [(499, 200)]).unwrap();

        // Interval ENDS exactly at the junction: all bases precede the deletion.
        let before = map.segments_for(480, 500).unwrap();
        assert!(
            !before
                .iter()
                .any(|s| matches!(s, HaplotypeSegment::Deletion { .. })),
            "interval ending at the junction does not span it: {before:?}"
        );
        assert!(InsertionCoordinateMap::deletion_runs_for_segments(&before).is_empty());

        // Interval STARTS exactly at the junction: all bases follow the deletion.
        let after = map.segments_for(500, 520).unwrap();
        assert!(
            !after
                .iter()
                .any(|s| matches!(s, HaplotypeSegment::Deletion { .. })),
            "interval starting at the junction does not span it: {after:?}"
        );
    }

    /// Insertions and deletions compose: shifts accumulate with sign.
    #[test]
    fn insertions_and_deletions_compose_their_shifts() {
        // +10 at anchor 100, then -30 after anchor 200.
        let map = InsertionCoordinateMap::with_deletions(1_000, [(100, vec![C; 10])], [(200, 30)])
            .unwrap();
        assert_eq!(map.haplotype_len(), 1_000 + 10 - 30);

        assert_eq!(map.reference_base_to_haplotype(100), Some(100));
        // Past the insertion: +10.
        assert_eq!(map.reference_base_to_haplotype(200), Some(210));
        // Inside the deletion: gone.
        assert_eq!(map.reference_base_to_haplotype(201), None);
        assert_eq!(map.reference_base_to_haplotype(230), None);
        // Past both: +10 - 30.
        assert_eq!(map.reference_base_to_haplotype(231), Some(211));
        assert_eq!(map.haplotype_base_to_reference(211), Some(231));
        // Still inside the inserted run, which has no reference coordinate.
        assert_eq!(map.haplotype_base_to_reference(105), None);
    }

    /// Degenerate deletions are refused rather than silently reinterpreted — each of
    /// these would make a projection ambiguous.
    #[test]
    fn degenerate_deletions_are_refused() {
        // Zero-length: not an event.
        assert!(InsertionCoordinateMap::with_deletions(1_000, [], [(100, 0)]).is_none());
        // Runs past the end of the reference.
        assert!(InsertionCoordinateMap::with_deletions(1_000, [], [(900, 200)]).is_none());
        // Anchor outside the reference.
        assert!(InsertionCoordinateMap::with_deletions(1_000, [], [(1_000, 10)]).is_none());
        // Two edits on one anchor.
        assert!(
            InsertionCoordinateMap::with_deletions(1_000, [(100, vec![C; 5])], [(100, 10)])
                .is_none()
        );
        // Overlapping deletions: the second's anchor was already consumed by the first.
        assert!(
            InsertionCoordinateMap::with_deletions(1_000, [], [(100, 50), (120, 10)]).is_none()
        );
        // Deleting the entire reference leaves no coordinate space to place anything in.
        assert!(InsertionCoordinateMap::with_deletions(10, [], [(0, 9)]).is_some());
        assert!(InsertionCoordinateMap::with_deletions(1, [], [(0, 1)]).is_none());
    }
}
