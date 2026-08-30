//! Turning a BAM into the `AlnRecord`s the metrics consume.
//!
//! ONE PASS, ALL REGIONS. The panel compares a multi-gigabyte real BAM against a simulated
//! one over many loci; querying each locus separately would re-open and re-seek per region.
//! Streaming once and dispatching each record to whichever regions contain it costs a single
//! read of each file and needs no index, which also means it works on a freshly written BAM
//! before anyone has run `samtools index`.
//!
//! This is the layer where the metrics stop being provable from literals and start depending
//! on noodles decoding what I think it does — so `tests/realism_reader.rs` cross-checks it
//! against `samtools view` on the same file rather than against my expectations.

use crate::metrics::{
    AlnRecord, RegionMetrics, candidate_breakpoints, depth_stats, depth_track, insert_stats,
};
use noodles::bam;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

/// A locus to measure. Half-open, 0-based internally; the BED the harness writes is the same.
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub contig: String,
    pub start: usize,
    pub end: usize,
}

impl Region {
    pub fn span(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}

#[derive(Debug)]
pub enum RealismError {
    Io(String),
    /// A region naming a contig the BAM header does not have is a HARD error, not a zero.
    /// A silently empty region reports "no artifacts here", which is indistinguishable from
    /// clean data — the exact confusion this whole panel exists to prevent (rule 4).
    UnknownContig(String),
    /// Likewise: a region that matched no reads at all cannot be reported as 0.0 artifacts.
    EmptyRegion(String),
}

impl std::fmt::Display for RealismError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RealismError::Io(m) => write!(f, "{m}"),
            RealismError::UnknownContig(c) => write!(
                f,
                "region names contig '{c}', which is not in the BAM header — the region was \
                 not measured, and an unmeasured region must not be reported as a clean one"
            ),
            RealismError::EmptyRegion(r) => write!(
                f,
                "region {r} contained no reads. Zero artifacts and zero reads look identical \
                 in the output, so this is an error rather than a measurement"
            ),
        }
    }
}

/// Which alignments contribute to a measurement.
///
/// Split out from `to_aln` because it cannot be exercised through a golden BAM: the simulator
/// emits exactly one record per read, so a fixture built from it contains no secondary or
/// supplementary alignments and including them changes nothing. Mutating the filter inside
/// `to_aln` survived the samtools cross-check for precisely that reason. As a free function it
/// is testable from flag literals, which is where the policy actually lives.
///
/// * `0x4` unmapped — no position, so nothing to attribute to a region.
/// * `0x100` secondary — an alternative placement of a read counted elsewhere.
/// * `0x800` supplementary — the other half of a split read. Counting it would double-count
///   the very clip boundaries this panel measures, and would inflate the artifact rate of any
///   aligner that emits them relative to one that does not.
///
/// Equivalent to `samtools view -F 0x904`, which is what the cross-check test passes.
pub fn countable(flags: u16) -> bool {
    flags & 0x4 == 0 && flags & 0x100 == 0 && flags & 0x800 == 0
}

/// Convert one BAM record. Returns `None` for records that cannot contribute a position —
/// unmapped reads and secondary/supplementary alignments.
///
/// Secondary and supplementary alignments are EXCLUDED deliberately. A supplementary
/// alignment is the other half of a split read, so counting it would double-count the very
/// clip boundaries this panel measures, and inflate the artifact rate of any aligner that
/// emits them. `samtools view -F 0x900` is the equivalent, and the cross-check test uses it.
fn to_aln(record: &bam::Record, min_mapq_keep: u8) -> Option<AlnRecord> {
    let flags = record.flags().bits();
    if !countable(flags) {
        return None;
    }
    let pos = record.alignment_start()?.ok()?.get() - 1;
    let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
    if mapq < min_mapq_keep {
        return None;
    }
    let mut cigar = Vec::new();
    for op in record.cigar().iter() {
        let op = op.ok()?;
        let c = match op.kind() {
            noodles::sam::alignment::record::cigar::op::Kind::Match => 'M',
            noodles::sam::alignment::record::cigar::op::Kind::Insertion => 'I',
            noodles::sam::alignment::record::cigar::op::Kind::Deletion => 'D',
            noodles::sam::alignment::record::cigar::op::Kind::Skip => 'N',
            noodles::sam::alignment::record::cigar::op::Kind::SoftClip => 'S',
            noodles::sam::alignment::record::cigar::op::Kind::HardClip => 'H',
            noodles::sam::alignment::record::cigar::op::Kind::Pad => 'P',
            noodles::sam::alignment::record::cigar::op::Kind::SequenceMatch => '=',
            noodles::sam::alignment::record::cigar::op::Kind::SequenceMismatch => 'X',
        };
        cigar.push((c, op.len()));
    }
    Some(AlnRecord {
        pos,
        mapq,
        flags,
        cigar,
        tlen: record.template_length() as i64,
    })
}

/// Stream `path` once and measure every region.
///
/// `min_clip` and `min_support` define a candidate breakpoint; `max_tlen` bounds what counts
/// as a library insert. They are parameters rather than constants because the whole point is
/// to compare two datasets under IDENTICAL settings — a gap measured with different thresholds
/// on each side would be measuring the thresholds.
pub fn measure(
    path: &Path,
    regions: &[Region],
    min_clip: usize,
    min_support: usize,
    max_tlen: i64,
    depth_lag: usize,
) -> Result<Vec<RegionMetrics>, RealismError> {
    let file =
        File::open(path).map_err(|e| RealismError::Io(format!("{}: {e}", path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|e| RealismError::Io(format!("{}: header: {e}", path.display())))?;

    // Reference id -> name, so a region's contig can be matched without string compares per
    // record. A region naming an absent contig is caught here, before any counting.
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();
    let index_of: HashMap<&str, usize> = names
        .iter()
        .map(|n| (n.as_str(), 0))
        .enumerate()
        .map(|(i, (n, _))| (n, i))
        .collect();

    let mut want: Vec<(usize, &Region)> = Vec::new();
    for r in regions {
        match index_of.get(r.contig.as_str()) {
            Some(i) => want.push((*i, r)),
            None => return Err(RealismError::UnknownContig(r.contig.clone())),
        }
    }

    let mut buckets: Vec<Vec<AlnRecord>> = vec![Vec::new(); regions.len()];
    for result in reader.records() {
        let record =
            result.map_err(|e| RealismError::Io(format!("{}: record: {e}", path.display())))?;
        let Some(aln) = to_aln(&record, 0) else {
            continue;
        };
        let Some(Ok(rid)) = record.reference_sequence_id() else {
            continue;
        };
        for (bi, (want_rid, region)) in want.iter().enumerate() {
            if rid == *want_rid && aln.pos >= region.start && aln.pos < region.end {
                buckets[bi].push(aln.clone());
            }
        }
    }

    let mut out = Vec::with_capacity(regions.len());
    for (bi, region) in regions.iter().enumerate() {
        let v = &buckets[bi];
        if v.is_empty() {
            return Err(RealismError::EmptyRegion(format!(
                "{}:{}-{}",
                region.contig, region.start, region.end
            )));
        }
        out.push(summarize(
            v,
            region.span(),
            region.start,
            min_clip,
            min_support,
            max_tlen,
            depth_lag,
        ));
    }
    Ok(out)
}

/// Shared summarizer. Both sides of the comparison go through this, so the two can never
/// drift apart into measuring different things.
pub fn summarize(
    v: &[AlnRecord],
    span: usize,
    start: usize,
    min_clip: usize,
    min_support: usize,
    max_tlen: i64,
    depth_lag: usize,
) -> RegionMetrics {
    let track = depth_track(v, start, span);
    RegionMetrics {
        reads: v.len(),
        span_bp: span,
        candidate_breakpoints: candidate_breakpoints(v, min_clip, min_support),
        improper_pairs: v.iter().filter(|r| !r.is_proper_pair()).count(),
        clipped_reads: v
            .iter()
            .filter(|r| r.leading_clip() >= min_clip || r.trailing_clip() >= min_clip)
            .count(),
        mapq0_reads: v.iter().filter(|r| r.mapq == 0).count(),
        insert: insert_stats(v, max_tlen),
        depth: depth_stats(&track, depth_lag),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The flag policy, from literals. A golden BAM cannot exercise this — it has no secondary
    /// or supplementary records — so mutating the filter inside `to_aln` passed the samtools
    /// cross-check. This is the test that fails instead.
    #[test]
    fn countable_excludes_unmapped_secondary_and_supplementary() {
        assert!(countable(0x0), "a plain mapped read must count");
        assert!(
            countable(0x2 | 0x40),
            "proper pair, first in pair, still counts"
        );

        assert!(!countable(0x4), "unmapped has no position to attribute");
        assert!(!countable(0x100), "secondary is a duplicate placement");
        assert!(
            !countable(0x800),
            "supplementary is the other half of a split read"
        );

        // Set alongside ordinary flags, they must still exclude.
        assert!(
            !countable(0x2 | 0x800),
            "supplementary hides behind proper-pair"
        );
        assert!(!countable(0x1 | 0x40 | 0x100));
    }

    /// Matches `samtools view -F 0x904`, which is the flag set the cross-check test uses.
    /// If these two ever disagree the comparison is measuring different record sets.
    #[test]
    fn countable_matches_the_samtools_filter_the_crosscheck_uses() {
        for flags in 0u16..=0x0FFF {
            assert_eq!(
                countable(flags),
                flags & 0x904 == 0,
                "disagreed with -F 0x904 at flags {flags:#06x}"
            );
        }
    }
}
