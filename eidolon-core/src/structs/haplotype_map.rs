//! Coordinate mapping between an altered haplotype and the reference it derives from.
//!
//! Lives in `eidolon-core` rather than beside the fragment sampler because BOTH
//! sides of read generation need it: the sampler (in the `eidolon` binary) places
//! fragments in haplotype coordinates, and the writer (`file_tools::fastq_tools`,
//! here in core) has to project those back to reference coordinates for BAM
//! records. Crate dependencies run `eidolon` -> `eidolon-core` only, so a shared
//! primitive has to live on this side.
//!
//! Moved here from `eidolon/src/gen_reads/utils/generate_fragments.rs` unchanged
//! (#516). The previous long-insertion attempt kept a second, parallel writer
//! precisely because this type was unreachable from core; every defect that
//! attempt shipped was something the parallel copy dropped.

use crate::structs::nucleotides::Nucleotide;

/// Coordinate map for one literal insertion in an altered haplotype.
///
/// `anchor` is the zero-based reference base represented by the VCF REF allele;
/// novel bases are inserted immediately after it. Reference coordinates are
/// therefore unchanged through `anchor`, and shift by `insertion_len` after
/// that point. This map is the foundation for emitting insertion-interior reads
/// while retaining reference coordinates for BAM records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub struct InsertionCoordinateMap {
    reference_len: usize,
    anchor: usize,
    insertion_len: usize,
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
        insertion_start: usize,
    },
}

impl InsertionCoordinateMap {
    pub fn new(reference_len: usize, anchor: usize, insertion_len: usize) -> Option<Self> {
        (reference_len > 0 && anchor < reference_len && insertion_len > 0).then_some(Self {
            reference_len,
            anchor,
            insertion_len,
        })
    }

    pub fn haplotype_len(self) -> usize {
        self.reference_len + self.insertion_len
    }

    pub fn reference_base_to_haplotype(self, reference_pos: usize) -> Option<usize> {
        if reference_pos >= self.reference_len {
            return None;
        }
        Some(if reference_pos <= self.anchor {
            reference_pos
        } else {
            reference_pos + self.insertion_len
        })
    }

    pub fn haplotype_base_to_reference(self, haplotype_pos: usize) -> Option<usize> {
        if haplotype_pos >= self.haplotype_len() {
            return None;
        }
        let insertion_start = self.anchor + 1;
        let insertion_end = insertion_start + self.insertion_len;
        if haplotype_pos >= insertion_start && haplotype_pos < insertion_end {
            return None;
        }
        Some(if haplotype_pos < insertion_start {
            haplotype_pos
        } else {
            haplotype_pos - self.insertion_len
        })
    }

    pub fn segments_for(self, start: usize, end: usize) -> Option<Vec<HaplotypeSegment>> {
        if start >= end || end > self.haplotype_len() {
            return None;
        }
        let insertion_start = self.anchor + 1;
        let insertion_end = insertion_start + self.insertion_len;
        let mut segments = Vec::with_capacity(3);

        if start < insertion_start {
            let ref_end = end.min(insertion_start);
            segments.push(HaplotypeSegment::Reference {
                hap_start: start,
                hap_end: ref_end,
                ref_start: start,
            });
        }
        if end > insertion_start && start < insertion_end {
            let hap_start = start.max(insertion_start);
            let hap_end = end.min(insertion_end);
            segments.push(HaplotypeSegment::Insertion {
                hap_start,
                hap_end,
                insertion_start: hap_start - insertion_start,
            });
        }
        if end > insertion_end {
            let hap_start = start.max(insertion_end);
            segments.push(HaplotypeSegment::Reference {
                hap_start,
                hap_end: end,
                ref_start: hap_start - self.insertion_len,
            });
        }
        Some(segments)
    }

    /// Materialize a haplotype interval from reference and inserted bases.
    /// The returned segments retain enough provenance for the caller to build
    /// reference-anchored BAM CIGAR operations.
    pub fn materialize_interval(
        self,
        reference: &[Nucleotide],
        inserted: &[Nucleotide],
        start: usize,
        end: usize,
    ) -> Option<(Vec<Nucleotide>, Vec<HaplotypeSegment>)> {
        if reference.len() != self.reference_len || inserted.len() != self.insertion_len {
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
                } => sequence
                    .extend_from_slice(&reference[ref_start..ref_start + (hap_end - hap_start)]),
                HaplotypeSegment::Insertion {
                    insertion_start,
                    hap_start,
                    hap_end,
                } => sequence.extend_from_slice(
                    &inserted[insertion_start..insertion_start + (hap_end - hap_start)],
                ),
            }
        }
        Some((sequence, segments))
    }

    /// Return the baseline CIGAR operation for every materialized base. This
    /// intentionally excludes sequencing-error operations; the read generator
    /// can layer those onto the baseline while retaining the insertion `I`
    /// operations supplied by the haplotype map.
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
    use crate::structs::nucleotides::Nucleotide::{A, C};

    #[test]
    fn insertion_coordinate_map_tracks_anchor_interior_and_tail() {
        let map = InsertionCoordinateMap::new(1_000, 500, 600).unwrap();
        assert_eq!(map.haplotype_len(), 1_600);
        assert_eq!(map.reference_base_to_haplotype(500), Some(500));
        assert_eq!(map.reference_base_to_haplotype(501), Some(1_101));
        assert_eq!(map.haplotype_base_to_reference(500), Some(500));
        assert_eq!(map.haplotype_base_to_reference(501), None);
        assert_eq!(map.haplotype_base_to_reference(1_101), Some(501));

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
                },
            ])
        );
        assert_eq!(
            map.segments_for(1_090, 1_120),
            Some(vec![
                HaplotypeSegment::Insertion {
                    hap_start: 1_090,
                    hap_end: 1_101,
                    insertion_start: 589,
                },
                HaplotypeSegment::Reference {
                    hap_start: 1_101,
                    hap_end: 1_120,
                    ref_start: 501,
                },
            ])
        );

        let reference = vec![A; 1_000];
        let inserted = vec![C; 600];
        let (sequence, segments) = map
            .materialize_interval(&reference, &inserted, 1_090, 1_120)
            .unwrap();
        assert_eq!(
            sequence,
            vec![C; 11]
                .into_iter()
                .chain(vec![A; 19])
                .collect::<Vec<_>>()
        );
        assert_eq!(segments.len(), 2);
        assert_eq!(
            InsertionCoordinateMap::cigar_ops_for_segments(&segments),
            vec!['I'; 11]
                .into_iter()
                .chain(vec!['M'; 19])
                .collect::<Vec<_>>()
        );
    }
}
