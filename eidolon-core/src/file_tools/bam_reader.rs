use std::{
    collections::{HashMap, HashSet},
    io,
    path::PathBuf,
};

use noodles::bam;
use noodles::sam::{
    self as sam,
    alignment::record::{
        Flags,
        cigar::op::Kind as CigarKind,
        data::field::{Tag, Value},
    },
};
use thiserror::Error;

use crate::structs::nucleotides::Nucleotide;

#[derive(Debug, Error)]
pub enum BamReaderError {
    #[error("I/O error reading BAM file: {0}")]
    IoError(#[from] io::Error),
}

const SKIP_FLAGS: Flags = Flags::UNMAPPED
    .union(Flags::SECONDARY)
    .union(Flags::SUPPLEMENTARY);

/// MD tag token used during mismatch extraction.
enum MdToken {
    Matches(usize),
    Mismatch(u8),
    Deletion,
}

fn parse_md(bytes: &[u8]) -> Vec<MdToken> {
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'0'..=b'9' => {
                let start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let n: usize = std::str::from_utf8(&bytes[start..i])
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                tokens.push(MdToken::Matches(n));
            }
            b'^' => {
                i += 1;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                tokens.push(MdToken::Deletion);
            }
            b if b.is_ascii_alphabetic() => {
                tokens.push(MdToken::Mismatch(bytes[i]));
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }
    tokens
}

/// Walks the parsed MD token list in sync with CIGAR M/=/X operations.
struct MdWalker {
    tokens: Vec<MdToken>,
    idx: usize,
    match_remaining: usize,
}

impl MdWalker {
    fn new(tokens: Vec<MdToken>) -> Self {
        Self {
            tokens,
            idx: 0,
            match_remaining: 0,
        }
    }

    /// Advance one M/=/X base. Returns the reference base byte for a mismatch, None for a match.
    fn next_alignment_base(&mut self) -> Option<u8> {
        loop {
            if self.match_remaining > 0 {
                self.match_remaining -= 1;
                return None;
            }
            match self.tokens.get(self.idx) {
                None => return None,
                Some(MdToken::Matches(0)) => {
                    self.idx += 1;
                }
                Some(MdToken::Matches(n)) => {
                    self.match_remaining = n - 1;
                    self.idx += 1;
                    return None;
                }
                Some(MdToken::Mismatch(b)) => {
                    let b = *b;
                    self.idx += 1;
                    return Some(b);
                }
                Some(MdToken::Deletion) => {
                    self.idx += 1;
                }
            }
        }
    }

    /// Advance past the next Deletion token (called once per D/N CIGAR operation).
    ///
    /// Redundant by construction: `next_alignment_base` also consumes `MdToken::Deletion`, so
    /// removing this call alone changes no output (verified by mutation — see
    /// `transition_observer_matches_samtools_generated_md_including_indels`). Kept as the explicit
    /// path so a D/N op has a visible handler; be aware that a defect introduced in either site
    /// alone is masked by the other.
    fn skip_deletion(&mut self) {
        while self.idx < self.tokens.len() {
            match &self.tokens[self.idx] {
                MdToken::Matches(0) => {
                    self.idx += 1;
                }
                MdToken::Deletion => {
                    self.idx += 1;
                    return;
                }
                _ => return,
            }
        }
    }
}

/// Minimum mapping quality for a read to contribute to the fragment length model.
/// A read is kept only if its MAPQ is strictly greater than this value.
const FRAG_FILTER_MAPQUAL: u8 = 10;

// ── Walker abstraction ───────────────────────────────────────────────────────

/// Per-record filter applied by `walk_bam` before any observer sees a record.
///
/// A record is kept iff every enabled predicate passes. `min_mapq` is treated
/// as a strict lower bound (records with `mq <= min_mapq` are dropped);
/// `min_mapq = 0` disables the MAPQ check entirely, including the implicit
/// "drop records with no MAPQ" behavior.
pub struct BamWalkFilter {
    pub min_mapq: u8,
    pub skip_flags: Flags,
    pub require_paired: bool,
    pub require_first_in_pair: bool,
    pub require_mate_mapped: bool,
    pub require_same_ref_as_mate: bool,
}

impl BamWalkFilter {
    /// Matches the legacy `read_fragment_lengths_bam` filter: paired, first-in-pair,
    /// mate mapped to the same reference, MAPQ > FRAG_FILTER_MAPQUAL.
    pub fn for_frag_length() -> Self {
        Self {
            min_mapq: FRAG_FILTER_MAPQUAL,
            skip_flags: SKIP_FLAGS,
            require_paired: true,
            require_first_in_pair: true,
            require_mate_mapped: true,
            require_same_ref_as_mate: true,
        }
    }

    /// Mate-overlap fitting (#779): both mates mapped to the same contig, no duplicates or
    /// QC failures, MAPQ above `min_mapq`.
    pub fn for_overlaps(min_mapq: u8) -> Self {
        Self {
            min_mapq,
            skip_flags: SKIP_FLAGS.union(Flags::DUPLICATE).union(Flags::QC_FAIL),
            require_paired: true,
            require_first_in_pair: false,
            require_mate_mapped: true,
            require_same_ref_as_mate: true,
        }
    }

    /// Matches the legacy `read_bam_transitions` filter: skip unmapped/secondary/
    /// supplementary only — no MAPQ filter, no pairing requirements.
    pub fn for_transitions() -> Self {
        Self {
            min_mapq: 0,
            skip_flags: SKIP_FLAGS,
            require_paired: false,
            require_first_in_pair: false,
            require_mate_mapped: false,
            require_same_ref_as_mate: false,
        }
    }

    /// Matches `samtools depth` defaults for coverage accumulation:
    /// skip unmapped/secondary/supplementary, no MAPQ filter, no pairing requirements.
    pub fn for_coverage() -> Self {
        Self {
            min_mapq: 0,
            skip_flags: SKIP_FLAGS,
            require_paired: false,
            require_first_in_pair: false,
            require_mate_mapped: false,
            require_same_ref_as_mate: false,
        }
    }
}

/// Visitor invoked once per record that passes the `BamWalkFilter`.
///
/// `on_start` runs once after the header is read but before any records are
/// dispatched. Observers that need to map reference_sequence_id → contig
/// name/length (e.g. `CoverageObserver`) cache that lookup here.
pub trait RecordObserver {
    fn on_start(&mut self, _header: &sam::Header) -> Result<(), BamReaderError> {
        Ok(())
    }
    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError>;
}

/// Counts returned by `walk_bam`.
#[derive(Debug, Default, Clone, Copy)]
pub struct WalkStats {
    pub records_seen: u64,
    pub records_kept: u64,
}

/// Single-pass walk over a BGZF BAM file. Each record that survives `filter`
/// is dispatched to every observer in order. Errors from any observer abort
/// the walk.
pub fn walk_bam(
    path: &PathBuf,
    filter: &BamWalkFilter,
    observers: &mut [&mut dyn RecordObserver],
) -> Result<WalkStats, BamReaderError> {
    let mut reader = std::fs::File::open(path).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    for obs in observers.iter_mut() {
        obs.on_start(&header)?;
    }

    let mut stats = WalkStats::default();

    for result in reader.records() {
        let record = result?;
        stats.records_seen += 1;

        let flags = record.flags();
        if flags.intersects(filter.skip_flags) {
            continue;
        }
        if filter.require_paired && !flags.is_segmented() {
            continue;
        }
        if filter.require_first_in_pair && !flags.is_first_segment() {
            continue;
        }
        if filter.require_mate_mapped && flags.is_mate_unmapped() {
            continue;
        }
        if filter.min_mapq > 0 {
            let mq = match record.mapping_quality() {
                Some(mq) => u8::from(mq),
                None => continue,
            };
            if mq <= filter.min_mapq {
                continue;
            }
        }
        if filter.require_same_ref_as_mate {
            let ref_id = record.reference_sequence_id().transpose()?;
            let mate_ref_id = record.mate_reference_sequence_id().transpose()?;
            if ref_id != mate_ref_id {
                continue;
            }
        }

        stats.records_kept += 1;
        for obs in observers.iter_mut() {
            obs.observe(&record)?;
        }
    }

    Ok(stats)
}

// ── Observers ────────────────────────────────────────────────────────────────

/// Collects absolute template lengths (TLEN > 0) for paired, first-in-pair,
/// confidently-mapped reads whose mate is mapped to the same reference.
///
/// Self-filtering: rejects records that don't meet its criteria internally, so
/// it stays correct under any `BamWalkFilter`. Calling with
/// `BamWalkFilter::for_frag_length()` makes the internal checks redundant
/// no-ops; calling with a looser filter (e.g. `for_coverage()` from the
/// unified `gen-bam-models` walker) still produces the same TLEN set.
#[derive(Debug, Default)]
pub struct FragLengthObserver {
    pub tlens: Vec<usize>,
}

impl RecordObserver for FragLengthObserver {
    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError> {
        let flags = record.flags();
        if !flags.is_segmented() || !flags.is_first_segment() || flags.is_mate_unmapped() {
            return Ok(());
        }
        let mq = match record.mapping_quality() {
            Some(mq) => u8::from(mq),
            None => return Ok(()),
        };
        if mq <= FRAG_FILTER_MAPQUAL {
            return Ok(());
        }
        let ref_id = record.reference_sequence_id().transpose()?;
        let mate_ref_id = record.mate_reference_sequence_id().transpose()?;
        if ref_id != mate_ref_id {
            return Ok(());
        }
        let tlen = record.template_length().unsigned_abs() as usize;
        if tlen > 0 {
            self.tlens.push(tlen);
        }
        Ok(())
    }
}

/// Reference positions to leave out of the transition count, as contig -> 1-based positions.
/// A sample's own variants mismatch the reference in every read that carries them, and they
/// are not sequencing errors.
pub type KnownSites = HashMap<String, HashSet<usize>>;

/// Accumulates a 4×4 read-vs-reference mismatch count matrix from MD tags.
/// `counts[ref_base][read_base]` follows ALLOWED_NUCS order (A=0, C=1, G=2, T=3).
/// Mismatches at a `KnownSites` position are counted in `masked` instead.
/// Pair with `BamWalkFilter::for_transitions()`.
#[derive(Debug, Default)]
pub struct TransitionObserver {
    pub counts: [[usize; 4]; 4],
    pub masked: usize,
    mask: Option<KnownSites>,
    contigs: Vec<String>,
}

impl TransitionObserver {
    pub fn with_mask(mask: KnownSites) -> Self {
        Self {
            mask: Some(mask),
            ..Self::default()
        }
    }
}

impl RecordObserver for TransitionObserver {
    fn on_start(&mut self, header: &sam::Header) -> Result<(), BamReaderError> {
        self.contigs = header
            .reference_sequences()
            .keys()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect();
        Ok(())
    }

    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError> {
        let md_bytes: Vec<u8> = match record.data().get(&Tag::MISMATCHED_POSITIONS) {
            Some(Ok(Value::String(s))) => s.iter().copied().collect(),
            _ => return Ok(()),
        };

        let tokens = parse_md(&md_bytes);
        let sequence = record.sequence();
        let mut walker = MdWalker::new(tokens);
        let mut read_pos = 0usize;
        // The mask's positions for this record's contig, if any. `ref_pos` is 1-based, like
        // the VCF the mask comes from, and advances on M/=/X and D/N only.
        let sites = match (&self.mask, record.reference_sequence_id()) {
            (Some(mask), Some(Ok(id))) => self.contigs.get(id).and_then(|c| mask.get(c)),
            _ => None,
        };
        let mut ref_pos = match record.alignment_start() {
            Some(Ok(p)) => usize::from(p),
            _ => 0,
        };

        for op_result in record.cigar().iter() {
            let op = op_result?;
            let len = op.len();

            match op.kind() {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    for _ in 0..len {
                        if let Some(ref_b) = walker.next_alignment_base()
                            && let Some(read_b) = sequence.get(read_pos)
                        {
                            let ri: usize = Nucleotide::from(ref_b as char).into();
                            let wi: usize = Nucleotide::from(read_b as char).into();
                            if ri < 4 && wi < 4 {
                                if sites.is_some_and(|s| s.contains(&ref_pos)) {
                                    self.masked += 1;
                                } else {
                                    self.counts[ri][wi] += 1;
                                }
                            }
                        }
                        read_pos += 1;
                        ref_pos += 1;
                    }
                }
                CigarKind::Insertion | CigarKind::SoftClip => {
                    read_pos += len;
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    walker.skip_deletion();
                    ref_pos += len;
                }
                CigarKind::HardClip | CigarKind::Pad => {}
            }
        }
        Ok(())
    }
}

/// What a mate-overlap fit saw (#779). Every category is reported so an exclusion is a number,
/// not a silent drop.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct OverlapCounts {
    /// Sequencing errors, `counts[reference base][erroneous base]`, ALLOWED_NUCS order.
    pub counts: [[usize; 4]; 4],
    /// Reference positions both mates of a fragment aligned a base to.
    pub overlapped_bases: usize,
    /// Overlapped positions where the two mates agree (variants included: both mates carry them).
    pub agreements: usize,
    /// Disagreements at a known-variant site, dropped: an error there can turn alt into ref.
    pub masked: usize,
    /// Disagreements where neither mate shows the reference base, dropped.
    pub neither_reference: usize,
    /// Kept records whose mate never arrived.
    pub unpaired: usize,
    /// `counts` split by the erroneous mate's sequencing cycle, in `OVERLAP_BINS` equal
    /// fractions of read length. Overlaps sit at reads' 3' ends, so this is how to see whether
    /// the spectrum changes along the read.
    pub by_bin: [[[usize; 4]; 4]; OVERLAP_BINS],
    /// Base observations per bin: each overlapped position counts once for each mate, in the
    /// bin of that mate's cycle. The denominator for a per-bin error rate.
    pub bin_bases: [usize; OVERLAP_BINS],
}

/// Read-position bins for `OverlapCounts::by_bin`.
pub const OVERLAP_BINS: usize = 5;

/// Counts sequencing errors where the two mates of a fragment overlap and disagree (#779).
///
/// Both mates read the same molecule, so a variant, a PCR error or damage already in the
/// fragment shows identically in both reads; only sequencing errors make them differ. At a
/// disagreement the mate that matches the reference is taken as right, and the other mate's
/// base is the error. Pair with `BamWalkFilter::for_overlaps()`.
#[derive(Debug, Default)]
pub struct OverlapObserver {
    pub counts: OverlapCounts,
    mask: Option<KnownSites>,
    contigs: Vec<String>,
    /// Reads waiting for their mate, by name: their contig and aligned bases.
    pending: HashMap<Vec<u8>, (usize, Vec<AlignedBase>)>,
}

/// One aligned base: 1-based reference position, the read's base, the reference base, and the
/// read-position bin of its sequencing cycle.
type AlignedBase = (usize, u8, u8, usize);

/// A record's aligned (M/=/X) bases, with the reference base read from its MD tag. `None`
/// when the record has no MD tag or no position.
fn aligned_bases(record: &bam::Record) -> Result<Option<Vec<AlignedBase>>, BamReaderError> {
    let md_bytes: Vec<u8> = match record.data().get(&Tag::MISMATCHED_POSITIONS) {
        Some(Ok(Value::String(s))) => s.iter().copied().collect(),
        _ => return Ok(None),
    };
    let mut ref_pos = match record.alignment_start() {
        Some(Ok(p)) => usize::from(p),
        _ => return Ok(None),
    };
    let sequence = record.sequence();
    // A reverse-strand record stores its read reverse-complemented, so its first sequenced
    // base (cycle 0) is the LAST one here.
    let n = sequence.len();
    let reverse = record.flags().is_reverse_complemented();
    let bin = |read_pos: usize| {
        let cycle = if reverse { n - 1 - read_pos } else { read_pos };
        cycle * OVERLAP_BINS / n
    };
    let mut walker = MdWalker::new(parse_md(&md_bytes));
    let mut read_pos = 0usize;
    let mut out = Vec::with_capacity(n);
    for op_result in record.cigar().iter() {
        let op = op_result?;
        let len = op.len();
        match op.kind() {
            CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                for _ in 0..len {
                    let mismatch = walker.next_alignment_base();
                    if let Some(read_b) = sequence.get(read_pos) {
                        out.push((ref_pos, read_b, mismatch.unwrap_or(read_b), bin(read_pos)));
                    }
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            CigarKind::Insertion | CigarKind::SoftClip => read_pos += len,
            CigarKind::Deletion | CigarKind::Skip => {
                walker.skip_deletion();
                ref_pos += len;
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }
    Ok(Some(out))
}

impl OverlapObserver {
    /// Compares two mates' aligned bases over the positions both cover.
    fn compare(&mut self, contig: usize, a: &[AlignedBase], b: &[AlignedBase]) {
        let sites = self
            .mask
            .as_ref()
            .and_then(|m| self.contigs.get(contig).and_then(|c| m.get(c)));
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            let ((pa, ra, refa, bin_a), (pb, rb, _, bin_b)) = (a[i], b[j]);
            if pa < pb {
                i += 1;
                continue;
            }
            if pb < pa {
                j += 1;
                continue;
            }
            i += 1;
            j += 1;
            let ia: usize = Nucleotide::from(ra as char).into();
            let ib: usize = Nucleotide::from(rb as char).into();
            let ir: usize = Nucleotide::from(refa as char).into();
            if ia >= 4 || ib >= 4 || ir >= 4 {
                continue; // an N in either read or the reference says nothing
            }
            self.counts.overlapped_bases += 1;
            self.counts.bin_bases[bin_a] += 1;
            self.counts.bin_bases[bin_b] += 1;
            if ia == ib {
                self.counts.agreements += 1;
            } else if sites.is_some_and(|s| s.contains(&pa)) {
                self.counts.masked += 1;
            } else if ia == ir {
                self.counts.counts[ir][ib] += 1;
                self.counts.by_bin[bin_b][ir][ib] += 1;
            } else if ib == ir {
                self.counts.counts[ir][ia] += 1;
                self.counts.by_bin[bin_a][ir][ia] += 1;
            } else {
                self.counts.neither_reference += 1;
            }
        }
    }
}

impl OverlapObserver {
    pub fn with_mask(mask: KnownSites) -> Self {
        Self {
            mask: Some(mask),
            ..Self::default()
        }
    }
}

impl RecordObserver for OverlapObserver {
    fn on_start(&mut self, header: &sam::Header) -> Result<(), BamReaderError> {
        self.contigs = header
            .reference_sequences()
            .keys()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect();
        Ok(())
    }

    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError> {
        let (Some(name), Some(Ok(contig))) = (record.name(), record.reference_sequence_id()) else {
            return Ok(());
        };
        let Some(bases) = aligned_bases(record)? else {
            return Ok(());
        };
        let key: &[u8] = name.as_ref();
        match self.pending.remove(key) {
            Some((mate_contig, mate)) if mate_contig == contig => {
                self.compare(contig, &mate, &bases)
            }
            Some(_) => {} // mates on different contigs share no positions
            None => {
                self.pending.insert(key.to_vec(), (contig, bases));
            }
        }
        Ok(())
    }
}

/// Mate-overlap substitution counts for a BAM (#779), leaving out `mask`'s positions.
pub fn read_bam_overlap_transitions(
    path: &PathBuf,
    mask: Option<KnownSites>,
    min_mapq: u8,
) -> Result<OverlapCounts, BamReaderError> {
    let mut obs = match mask {
        Some(m) => OverlapObserver::with_mask(m),
        None => OverlapObserver::default(),
    };
    walk_bam(
        path,
        &BamWalkFilter::for_overlaps(min_mapq),
        &mut [&mut obs],
    )?;
    obs.counts.unpaired = obs.pending.len();
    Ok(obs.counts)
}

/// Accumulates per-base reference coverage from aligned records.
///
/// `on_start` snapshots the contig name / length tables from the BAM header.
/// `observe` walks each kept record's CIGAR M/=/X spans starting at
/// `alignment_start` and increments per-base depth counters. CIGAR D/N
/// advance the reference position without incrementing; I/S/H/P do not
/// consume reference positions and leave depth unchanged.
///
/// Records with no `reference_sequence_id` or no `alignment_start` are
/// silently skipped — pair with `BamWalkFilter::for_coverage()` which
/// drops unmapped records via `SKIP_FLAGS`.
///
/// Depth arrays are allocated lazily, so contigs with no observed
/// records consume no memory.
#[derive(Debug, Default)]
pub struct CoverageObserver {
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
    depths_by_id: Vec<Option<Vec<u32>>>,
}

impl CoverageObserver {
    /// Consumes the observer and returns depth arrays keyed by contig name.
    /// Only contigs that received at least one record are included.
    pub fn into_by_contig(self) -> HashMap<String, Vec<u32>> {
        self.contig_names
            .into_iter()
            .zip(self.depths_by_id)
            .filter_map(|(name, depths)| depths.map(|d| (name, d)))
            .collect()
    }
}

impl RecordObserver for CoverageObserver {
    fn on_start(&mut self, header: &sam::Header) -> Result<(), BamReaderError> {
        self.contig_names.clear();
        self.contig_lengths.clear();
        for (name, ref_seq) in header.reference_sequences() {
            let name_bytes: &[u8] = name.as_ref();
            let name_str = std::str::from_utf8(name_bytes)
                .map(|s| s.to_string())
                .unwrap_or_else(|_| String::from_utf8_lossy(name_bytes).into_owned());
            self.contig_names.push(name_str);
            self.contig_lengths.push(ref_seq.length().get());
        }
        self.depths_by_id = (0..self.contig_names.len()).map(|_| None).collect();
        Ok(())
    }

    fn observe(&mut self, record: &bam::Record) -> Result<(), BamReaderError> {
        let ref_id = match record.reference_sequence_id().transpose()? {
            Some(id) => id,
            None => return Ok(()),
        };
        if ref_id >= self.depths_by_id.len() {
            return Ok(());
        }
        let start_1based = match record.alignment_start() {
            Some(Ok(p)) => p.get(),
            _ => return Ok(()),
        };
        if start_1based == 0 {
            return Ok(());
        }
        let mut ref_pos = start_1based - 1;
        let contig_len = self.contig_lengths[ref_id];
        let depths = self.depths_by_id[ref_id].get_or_insert_with(|| vec![0u32; contig_len]);

        for op_result in record.cigar().iter() {
            let op = op_result?;
            let len = op.len();
            match op.kind() {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    let end = (ref_pos + len).min(depths.len());
                    for i in ref_pos..end {
                        depths[i] = depths[i].saturating_add(1);
                    }
                    ref_pos += len;
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    ref_pos += len;
                }
                CigarKind::Insertion
                | CigarKind::SoftClip
                | CigarKind::HardClip
                | CigarKind::Pad => {}
            }
        }
        Ok(())
    }
}

// ── Public API (wrappers over walk_bam) ──────────────────────────────────────

/// Reads a BAM or SAM file and returns the absolute template lengths (TLEN) for
/// paired, first-in-pair reads that are confidently mapped to the same reference
/// as their mate and have mapping quality > FRAG_FILTER_MAPQUAL.
///
/// Files with a `.sam` extension are read as plain-text SAM; all others are
/// treated as BGZF-compressed BAM.
pub fn read_fragment_lengths(path: &PathBuf) -> Result<Vec<usize>, BamReaderError> {
    let is_sam = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.eq_ignore_ascii_case("sam"))
        .unwrap_or(false);

    if is_sam {
        read_fragment_lengths_sam(path)
    } else {
        let mut obs = FragLengthObserver::default();
        walk_bam(path, &BamWalkFilter::for_frag_length(), &mut [&mut obs])?;
        Ok(obs.tlens)
    }
}

fn read_fragment_lengths_sam(path: &PathBuf) -> Result<Vec<usize>, BamReaderError> {
    use noodles::sam;
    use std::io::BufReader;

    let file = std::fs::File::open(path)?;
    let mut reader = sam::io::Reader::new(BufReader::new(file));
    let header = reader.read_header()?;
    let mut tlens = Vec::new();
    for result in reader.records() {
        let record = result?;
        let flags = record.flags()?;
        if !flags.is_segmented() || !flags.is_first_segment() {
            continue;
        }
        if flags.intersects(SKIP_FLAGS) || flags.is_mate_unmapped() {
            continue;
        }
        let mq: u8 = match record.mapping_quality() {
            Some(Ok(mq)) => u8::from(mq),
            _ => continue,
        };
        if mq <= FRAG_FILTER_MAPQUAL {
            continue;
        }
        let ref_id = record.reference_sequence_id(&header).transpose()?;
        let mate_ref_id = record.mate_reference_sequence_id(&header).transpose()?;
        if ref_id != mate_ref_id {
            continue;
        }
        let tlen = record.template_length()?.unsigned_abs() as usize;
        if tlen > 0 {
            tlens.push(tlen);
        }
    }
    Ok(tlens)
}

/// Reads a BAM file and accumulates a raw 4×4 SNP mismatch count matrix.
///
/// `counts[ref_base][read_base]` is incremented for each position where a read
/// base differs from the reference base. Both axes follow ALLOWED_NUCS order
/// (A=0, C=1, G=2, T=3). Bases that are not ACGT (N, IUPAC ambiguity codes,
/// etc.) are silently ignored.
///
/// Unmapped, secondary, and supplementary records are skipped. Records that
/// lack an MD tag are silently skipped; if no records with MD tags are found
/// the returned matrix will be all-zeros.
pub fn read_bam_transitions(path: &PathBuf) -> Result<[[usize; 4]; 4], BamReaderError> {
    let mut obs = TransitionObserver::default();
    walk_bam(path, &BamWalkFilter::for_transitions(), &mut [&mut obs])?;
    Ok(obs.counts)
}

/// As `read_bam_transitions`, leaving out mismatches at `mask`'s positions. Returns the counts
/// and how many mismatches were masked.
pub fn read_bam_transitions_masked(
    path: &PathBuf,
    mask: KnownSites,
) -> Result<([[usize; 4]; 4], usize), BamReaderError> {
    let mut obs = TransitionObserver::with_mask(mask);
    walk_bam(path, &BamWalkFilter::for_transitions(), &mut [&mut obs])?;
    Ok((obs.counts, obs.masked))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind as CigarOpKind;

    // (flags, mapq, tlen, reference_sequence_id, mate_reference_sequence_id)
    type TestBamRecord = (Flags, Option<u8>, i32, Option<usize>, Option<usize>);
    // (cigar_kind, cigar_op_length)
    type TestCigarOp = (CigarOpKind, usize);
    // (reference_sequence_id, alignment_start_1based, cigar_ops)
    type TestCoverageRecord<'a> = (usize, usize, &'a [TestCigarOp]);

    /// `read_bam_transitions` against a **real** aligner's output, with a known answer computed
    /// outside eidolon.
    ///
    /// Everything else that tests this path — including the samtools-calmd oracle above — feeds
    /// the reader a BAM written by the test itself: a handful of records, hand-chosen CIGARs, no
    /// surprises. Real alignments are not like that. This fixture has soft clips, indels,
    /// duplicate-marked and unmapped records, `N` bases, and MD strings emitted by bwa 0.5.9
    /// rather than by us.
    ///
    /// Fixture: 1000 Genomes HG00096 low-coverage alignment, GRCh37 `20:1000000-1050000`,
    /// unmodified. Full provenance and redistribution terms in
    /// `test_data/HG00096.chr20_1Mb.README.md`.
    ///
    /// KNOWN ANSWER. Derived two independent ways, neither of them this code:
    ///
    /// 1. An awk walk over `MD` + `CIGAR` (`samtools view -F 0x904`, matching `SKIP_FLAGS`),
    ///    tallying `counts[ref][read]` — 2548 usable records, **804** substitutions.
    /// 2. The aligner's own tags: `sum(NM) - sum(inserted + deleted bases)` = **804**.
    ///
    /// bwa computed those NM values at alignment time in 2012, knowing nothing about this
    /// codebase, so (2) is about as independent as an oracle gets. The two agreeing to the
    /// record is what licenses the per-cell literals below.
    #[test]
    fn transition_counts_match_an_independent_oracle_on_a_real_bam() {
        let path = PathBuf::from("test_data/HG00096.chr20_1Mb.bam");
        assert!(path.is_file(), "real-data fixture missing: {path:?}");

        let counts = read_bam_transitions(&path).unwrap();

        // ALLOWED_NUCS order: A=0, C=1, G=2, T=3.
        let expected: [[usize; 4]; 4] = [
            [0, 97, 100, 53],
            [62, 0, 15, 75],
            [91, 21, 0, 39],
            [46, 120, 85, 0],
        ];
        assert_eq!(
            counts, expected,
            "real-BAM substitution counts disagree with the awk + NM oracle"
        );

        // The total is the figure both oracles agree on, asserted separately: a per-cell
        // regression and a wholesale miscount are different failures, and the total is the one
        // an outside reader can re-derive with samtools alone.
        let total: usize = counts.iter().flatten().sum();
        assert_eq!(
            total, 804,
            "total substitutions must match sum(NM) - indels"
        );

        // The diagonal is not a substitution and must never be counted, whatever the MD says.
        for (i, row) in counts.iter().enumerate() {
            assert_eq!(row[i], 0, "diagonal cell [{i}][{i}] must be zero");
        }

        // Every off-diagonal cell is populated, which is what makes the equality assertion
        // above meaningful: a parser that silently dropped a whole class of substitution would
        // leave a zero here rather than failing some aggregate.
        for (i, row) in counts.iter().enumerate() {
            for (j, &c) in row.iter().enumerate() {
                if i != j {
                    assert!(
                        c > 0,
                        "cell [{i}][{j}] is zero; the fixture should populate all 12"
                    );
                }
            }
        }
    }

    #[test]
    fn test_parse_md_simple() {
        // "10A5G3": 10 matches, mismatch ref=A, 5 matches, mismatch ref=G, 3 matches
        let tokens = parse_md(b"10A5G3");
        let mut walker = MdWalker::new(tokens);
        // 10 matches
        for _ in 0..10 {
            assert!(walker.next_alignment_base().is_none());
        }
        // mismatch ref=A
        assert_eq!(walker.next_alignment_base(), Some(b'A'));
        // 5 matches
        for _ in 0..5 {
            assert!(walker.next_alignment_base().is_none());
        }
        // mismatch ref=G
        assert_eq!(walker.next_alignment_base(), Some(b'G'));
        // 3 matches
        for _ in 0..3 {
            assert!(walker.next_alignment_base().is_none());
        }
        // exhausted
        assert!(walker.next_alignment_base().is_none());
    }

    #[test]
    fn test_parse_md_with_deletion() {
        // "5^AT3": 5 matches, deletion AT, 3 matches
        let tokens = parse_md(b"5^AT3");
        let mut walker = MdWalker::new(tokens);
        for _ in 0..5 {
            assert!(walker.next_alignment_base().is_none());
        }
        walker.skip_deletion();
        for _ in 0..3 {
            assert!(walker.next_alignment_base().is_none());
        }
        assert!(walker.next_alignment_base().is_none());
    }

    #[test]
    fn test_parse_md_leading_zero_mismatch() {
        // "0T3": immediate mismatch at position 0, ref=T, 3 matches
        let tokens = parse_md(b"0T3");
        let mut walker = MdWalker::new(tokens);
        assert_eq!(walker.next_alignment_base(), Some(b'T'));
        for _ in 0..3 {
            assert!(walker.next_alignment_base().is_none());
        }
    }

    /// Builds a minimal BGZF BAM at `path` with one record per `records` entry.
    /// Fields: (flags, mapq, tlen, reference_sequence_id, mate_reference_sequence_id).
    /// Two reference sequences (chr1, chr2) are declared so mate-on-different-ref
    /// scenarios are expressible.
    fn write_test_bam(path: &std::path::PathBuf, records: &[TestBamRecord]) {
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    MappingQuality,
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
            .add_reference_sequence(
                b"chr2".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for &(flags, mapq, tlen, ref_id, mate_ref_id) in records {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = flags;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(b"ACGT".as_ref());
            if let Some(mq) = mapq {
                *record.mapping_quality_mut() = Some(MappingQuality::try_from(mq).unwrap());
            }
            *record.reference_sequence_id_mut() = ref_id;
            *record.mate_reference_sequence_id_mut() = mate_ref_id;
            *record.template_length_mut() = tlen;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    /// Build a BAM with arbitrary CIGAR + sequence + MD per record.
    ///
    /// Records are `(cigar_ops, sequence, md)`. Unlike the `write_test_bam` above — which is
    /// fixed at `4M` with no MD — this exists to exercise `TransitionObserver`'s CIGAR
    /// handling, which the transition tests otherwise never reach.
    /// The three samtools-generated reads documented on
    /// `transition_observer_matches_samtools_generated_md_including_indels`.
    fn write_samtools_md_fixture(path: &std::path::PathBuf) {
        write_md_cigar_bam(
            path,
            &[
                (&[(CigarOpKind::Match, 12)], "GCGAGTTCAAAA", "2T4A4"),
                (
                    &[
                        (CigarOpKind::Match, 5),
                        (CigarOpKind::Deletion, 2),
                        (CigarOpKind::Match, 5),
                    ],
                    "GCTAGAACAA",
                    "5^TT2A2",
                ),
                (
                    &[
                        (CigarOpKind::Match, 5),
                        (CigarOpKind::Insertion, 2),
                        (CigarOpKind::Match, 5),
                    ],
                    "GCTAGTTTGAAA",
                    "6T3",
                ),
            ],
        );
    }

    fn write_md_cigar_bam(
        path: &std::path::PathBuf,
        records: &[(&[(CigarOpKind, usize)], &str, &str)],
    ) {
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{cigar::Op, data::field::Tag},
                record_buf::{Cigar, Sequence, data::field::Value as BufValue},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"H1N1_HA".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_701).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for (ops, seq, md) in records {
            let cigar: Cigar = ops.iter().map(|&(k, n)| Op::new(k, n)).collect();
            let mut record = RecordBuf::default();
            // RecordBuf::default() is UNMAPPED, which SKIP_FLAGS drops before any observer sees it.
            *record.flags_mut() = Flags::empty();
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(seq.as_bytes());
            *record.reference_sequence_id_mut() = Some(0);
            *record.alignment_start_mut() = noodles::core::Position::new(501);
            record
                .data_mut()
                .insert(Tag::MISMATCHED_POSITIONS, BufValue::from(*md));
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    /// `TransitionObserver` against MD strings produced by **samtools**, not by us.
    ///
    /// WHY THIS EXISTS. Every other transition test builds its MD with
    /// `gen_seq_error_model::utils::runner::write_test_bam`, which *generates the MD string
    /// itself* by zipping equal-length ref/read slices. That encoder and this parser share one
    /// understanding of MD encoding, so if the understanding is wrong both are wrong identically
    /// and the tests pass regardless. It also cannot emit `^` at all, so the deletion branch of
    /// `observe` — and the insertion branch — were never reached by any test.
    ///
    /// The MD strings below were generated by **samtools 1.19.2** against the in-repo fixture
    /// `eidolon-core/test_data/H1N1.fa`, all three reads at `H1N1_HA:501`. To regenerate:
    ///
    /// ```text
    /// tr -d '\r' < eidolon-core/test_data/H1N1.fa > ref.fa   # the fixture is CRLF
    /// samtools faidx ref.fa H1N1_HA:501-530                  # GCTAGTTAAAAAAGGAAATTCATACCCAAA
    /// # a SAM with the three (CIGAR, sequence) rows below, POS 501, MAPQ 60, FLAG 0
    /// samtools view -b in.sam | samtools calmd -b - ref.fa > out.bam
    /// ```
    ///
    ///   read  CIGAR      sequence        samtools MD    NM
    ///   r1    12M        GCGAGTTCAAAA    2T4A4          2
    ///   r2    5M2D5M     GCTAGAACAA      5^TT2A2        3
    ///   r3    5M2I5M     GCTAGTTTGAAA    6T3            3
    ///
    /// Baked in as literals rather than computed at test time because CI has no samtools.
    ///
    /// KNOWN ANSWER, computed by hand from the reference window rather than from this code:
    /// r1 contributes T→G and A→C; r2 contributes A→C (its `^TT` is *deleted reference*, not a
    /// substitution); r3 contributes T→G (its two inserted bases are invisible to MD, and the
    /// mismatch sits AFTER them — a parser that ignored the `I` op would attribute the wrong
    /// read base). Total: **A→C ×2, T→G ×2, everything else zero.**
    ///
    /// NON-VACUITY, by mutation of `TransitionObserver::observe`:
    ///
    ///   mutation                                              result
    ///   `Insertion | SoftClip` arm made a no-op               CAUGHT (T row became G×1, T×1)
    ///   `Deletion` arm also advances `read_pos += len`         CAUGHT (A row became A×1, C×1)
    ///   `Deletion` arm made a no-op                           SURVIVED — see below
    ///   that no-op **plus** `MdToken::Deletion` no longer      CAUGHT (A→C fell to 1)
    ///   skipped inside `next_alignment_base`
    ///
    /// The lone survivor is an equivalent mutant, not a coverage hole: `MdToken::Deletion` is a
    /// zero-width marker that `next_alignment_base` also consumes, so no MD/CIGAR pair can
    /// distinguish `skip_deletion()` running from it not running. The last row proves the
    /// deletion *semantics* are covered — remove both copies and the count is wrong.
    #[test]
    fn transition_observer_matches_samtools_generated_md_including_indels() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("samtools_md.bam");
        write_samtools_md_fixture(&path);

        let counts = read_bam_transitions(&path).unwrap();

        // ALLOWED_NUCS order: A=0, C=1, G=2, T=3.
        let mut expected = [[0usize; 4]; 4];
        expected[0][1] = 2; // A -> C
        expected[3][2] = 2; // T -> G
        assert_eq!(
            counts, expected,
            "counts disagree with the hand-computed answer for samtools' MD.\n  got:      {counts:?}\n  expected: {expected:?}"
        );

        // MUST NOT FIRE: the deleted reference bases in r2's `^TT` are not substitutions. If the
        // Deletion arm stopped calling skip_deletion, the walker would desynchronise and T rows
        // would pick up spurious counts.
        assert_eq!(counts[3][3], 0, "T->T is not a substitution");
        assert_eq!(
            counts[3].iter().sum::<usize>(),
            2,
            "the T row must total exactly the two real T->G events; \
             extra counts mean deleted or inserted bases leaked in"
        );
    }

    /// 40 bp reference for the overlap fixtures: A x5, C x5, G x5, T x5, twice.
    const OVERLAP_REF: &[u8] = b"AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT";

    /// Paired, all-M records against OVERLAP_REF, with MD computed from the reference.
    /// Each entry is `(name, flags, 1-based start, sequence, mate start)`.
    fn write_pair_bam(path: &std::path::PathBuf, records: &[(&str, Flags, usize, &str, usize)]) {
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{MappingQuality, cigar::Op, data::field::Tag},
                record_buf::{Cigar, Sequence, data::field::Value as BufValue},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(40).unwrap()),
            )
            .build();
        let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
        writer.write_header(&header).unwrap();
        for (name, flags, start, seq, mate) in records {
            let reference = &OVERLAP_REF[start - 1..start - 1 + seq.len()];
            let (mut md, mut run) = (String::new(), 0usize);
            for (r, q) in reference.iter().zip(seq.bytes()) {
                if *r == q {
                    run += 1;
                } else {
                    md.push_str(&format!("{run}{}", *r as char));
                    run = 0;
                }
            }
            md.push_str(&run.to_string());
            let mut rec = RecordBuf::default();
            *rec.name_mut() = Some(name.as_bytes().into());
            *rec.flags_mut() = *flags;
            *rec.cigar_mut() = [Op::new(CigarOpKind::Match, seq.len())]
                .into_iter()
                .collect::<Cigar>();
            *rec.sequence_mut() = Sequence::from(seq.as_bytes());
            *rec.reference_sequence_id_mut() = Some(0);
            *rec.alignment_start_mut() = noodles::core::Position::new(*start);
            *rec.mate_reference_sequence_id_mut() = Some(0);
            *rec.mate_alignment_start_mut() = noodles::core::Position::new(*mate);
            *rec.mapping_quality_mut() = MappingQuality::new(60);
            rec.data_mut()
                .insert(Tag::MISMATCHED_POSITIONS, BufValue::from(md));
            writer.write_alignment_record(&header, &rec).unwrap();
        }
    }

    /// #779, KNOWN ANSWER built by hand. Mates are 20 bp at 1-20 and 11-30, so they overlap on
    /// reference 11-20 (GGGGGTTTTT). Each pair plants one case from the issue's table:
    ///
    ///   p1  pos 12: G / A           an error in mate 2 -> G->A. Its mate-1 mismatch at pos 3
    ///                               is outside the overlap and must NOT count.
    ///   p2  pos 17: C / C           a hom variant: the mates agree, nothing counted
    ///   p3  pos 13: C / T (ref G)   a het alt plus an error: neither is reference, dropped
    ///   p4  pos 18: A / T (ref T)   an error turned alt back into ref: masked. Without the
    ///                               mask it is misread as T->A, which is why the mask exists.
    ///   p5  duplicate, with a disagreement   skipped by the filter
    ///   p6  mate 1 only             unpaired
    ///
    /// Four pairs x 10 overlapped bases = 40; three disagree, so 37 agree.
    #[test]
    fn overlap_disagreements_count_only_sequencing_errors() {
        let first = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT;
        let second = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::LAST_SEGMENT
            | Flags::REVERSE_COMPLEMENTED;
        let m1 = "AAAAACCCCCGGGGGTTTTT";
        let m2 = "GGGGGTTTTTAAAAACCCCC";
        let with = |s: &str, at: usize, b: char| -> String {
            let mut v: Vec<char> = s.chars().collect();
            v[at] = b;
            v.into_iter().collect()
        };
        // Offsets: mate 1 covers 1-20 (index = pos - 1), mate 2 covers 11-30 (index = pos - 11).
        let p1a = with(m1, 2, 'T');
        let p1b = with(m2, 12 - 11, 'A');
        let p2a = with(m1, 17 - 1, 'C');
        let p2b = with(m2, 17 - 11, 'C');
        let p3a = with(m1, 13 - 1, 'C');
        let p3b = with(m2, 13 - 11, 'T');
        let p4a = with(m1, 18 - 1, 'A');
        let p5b = with(m2, 14 - 11, 'A');
        let dup = Flags::DUPLICATE;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("overlap.bam");
        write_pair_bam(
            &path,
            &[
                ("p1", first, 1, &p1a, 11),
                ("p2", first, 1, &p2a, 11),
                ("p3", first, 1, &p3a, 11),
                ("p4", first, 1, &p4a, 11),
                ("p5", first | dup, 1, m1, 11),
                ("p6", first, 1, m1, 11),
                ("p1", second, 11, &p1b, 1),
                ("p2", second, 11, &p2b, 1),
                ("p3", second, 11, &p3b, 1),
                ("p4", second, 11, m2, 1),
                ("p5", second | dup, 11, &p5b, 1),
            ],
        );
        let mut mask = KnownSites::new();
        mask.entry("chr1".to_string()).or_default().insert(18);

        let got = read_bam_overlap_transitions(&path, Some(mask), 0).unwrap();
        let mut counts = [[0usize; 4]; 4];
        counts[2][0] = 1; // G -> A
        // Read-position bins, by hand: over the overlap, mate 1 (forward, starts at 1) is at
        // cycles 10-19 and mate 2 (reverse, starts at 11) at cycles 19 down to 10, so each pair
        // puts 2+2 observations in bin 2 (cycles 10-11), 4+4 in bin 3 (12-15), 4+4 in bin 4.
        // p1's error is mate 2's base at reference 12, its query offset 1 = cycle 18: bin 4.
        let mut by_bin = [[[0usize; 4]; 4]; OVERLAP_BINS];
        by_bin[4][2][0] = 1;
        assert_eq!(
            got,
            OverlapCounts {
                counts,
                overlapped_bases: 40,
                agreements: 37,
                masked: 1,
                neither_reference: 1,
                unpaired: 1,
                by_bin,
                bin_bases: [0, 0, 16, 32, 32],
            }
        );

        // The same fixture unmasked: p4's alt-to-ref error is misread as T->A.
        let unmasked = read_bam_overlap_transitions(&path, None, 0).unwrap();
        counts[3][0] = 1; // T -> A
        assert_eq!((unmasked.counts, unmasked.masked), (counts, 0));
    }

    /// The samtools-MD fixture above, masked. Its four mismatches sit at known reference
    /// positions (hand-derived from the MD strings and CIGARs, 1-based):
    ///
    ///   r1  12M      T->G at 503, A->C at 508
    ///   r2  5M2D5M   A->C at 510, after the deletion of 506-507
    ///   r3  5M2I5M   T->G at 507, after two inserted bases that consume no reference
    ///
    /// Masking 503 and 510 must leave one T->G (507) and one A->C (508), and report 2 masked.
    /// That needs the reference position right across both indels: advancing it on the
    /// insertion, or not on the deletion, masks the wrong sites.
    #[test]
    fn masked_sites_are_left_out_of_the_transition_counts() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("samtools_md.bam");
        write_samtools_md_fixture(&path);
        let mask = |sites: &[(&str, usize)]| -> KnownSites {
            let mut m = KnownSites::new();
            for (c, p) in sites {
                m.entry(c.to_string()).or_default().insert(*p);
            }
            m
        };

        let (counts, masked) =
            read_bam_transitions_masked(&path, mask(&[("H1N1_HA", 503), ("H1N1_HA", 510)]))
                .unwrap();
        let mut expected = [[0usize; 4]; 4];
        expected[0][1] = 1; // A -> C at 508
        expected[3][2] = 1; // T -> G at 507
        assert_eq!(counts, expected, "masking 503 and 510 left {counts:?}");
        assert_eq!(masked, 2);

        // r3's mismatch is the one after the insertion; masking it alone tests that path.
        let (counts, masked) =
            read_bam_transitions_masked(&path, mask(&[("H1N1_HA", 507)])).unwrap();
        assert_eq!(
            (counts[3][2], counts[0][1], masked),
            (1, 2, 1),
            "masking 507 left {counts:?}"
        );

        // MUST NOT FIRE: a site with no mismatch, and the right position on the wrong contig.
        let (counts, masked) =
            read_bam_transitions_masked(&path, mask(&[("H1N1_HA", 501), ("H1N1_NA", 503)]))
                .unwrap();
        assert_eq!(counts, read_bam_transitions(&path).unwrap());
        assert_eq!(masked, 0);
    }
    #[test]
    fn test_walk_bam_dispatches_to_multiple_observers() {
        // One kept record reaches every observer in the slice.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                Some(30),
                200,
                Some(0),
                Some(0),
            )],
        );
        let mut frag = FragLengthObserver::default();
        let mut trans = TransitionObserver::default();
        let stats = walk_bam(
            &path,
            &BamWalkFilter::for_frag_length(),
            &mut [&mut frag, &mut trans],
        )
        .unwrap();
        assert_eq!(stats.records_seen, 1);
        assert_eq!(stats.records_kept, 1);
        assert_eq!(frag.tlens, vec![200]);
        // No MD tag → no mismatches counted, but the observer was reached.
        assert_eq!(trans.counts, [[0usize; 4]; 4]);
    }

    #[test]
    fn test_walk_bam_min_mapq_is_strict_lower_bound() {
        // mq == min_mapq is dropped; mq > min_mapq is kept. Pins the `mq <= min_mapq`
        // semantics inherited from the original FRAG_FILTER_MAPQUAL check.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        let flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        write_test_bam(
            &path,
            &[
                (flags, Some(10), 100, Some(0), Some(0)),
                (flags, Some(11), 200, Some(0), Some(0)),
            ],
        );
        let mut frag = FragLengthObserver::default();
        let stats = walk_bam(&path, &BamWalkFilter::for_frag_length(), &mut [&mut frag]).unwrap();
        assert_eq!(stats.records_seen, 2);
        assert_eq!(stats.records_kept, 1);
        assert_eq!(frag.tlens, vec![200]);
    }

    #[test]
    fn test_walk_bam_min_mapq_zero_keeps_records_with_no_mapq() {
        // The transitions filter sets min_mapq = 0, which must disable the MAPQ
        // check entirely — including the "no MAPQ → drop" behavior that applies
        // when min_mapq > 0.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(&path, &[(Flags::empty(), None, 0, Some(0), Some(0))]);
        let mut trans = TransitionObserver::default();
        let stats = walk_bam(&path, &BamWalkFilter::for_transitions(), &mut [&mut trans]).unwrap();
        assert_eq!(stats.records_seen, 1);
        assert_eq!(stats.records_kept, 1);
    }

    #[test]
    fn test_walk_bam_skip_flags_drops_secondary() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(Flags::SECONDARY, Some(30), 100, Some(0), Some(0))],
        );
        let mut trans = TransitionObserver::default();
        let stats = walk_bam(&path, &BamWalkFilter::for_transitions(), &mut [&mut trans]).unwrap();
        assert_eq!(stats.records_seen, 1);
        assert_eq!(stats.records_kept, 0);
    }

    #[test]
    fn test_walk_bam_require_same_ref_as_mate() {
        // Paired record with mate on a different reference is dropped by the
        // frag-length filter.
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_test_bam(
            &path,
            &[(
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                Some(30),
                100,
                Some(0),
                Some(1),
            )],
        );
        let mut frag = FragLengthObserver::default();
        let stats = walk_bam(&path, &BamWalkFilter::for_frag_length(), &mut [&mut frag]).unwrap();
        assert_eq!(stats.records_kept, 0);
    }

    /// Builds a minimal BGZF BAM at `path` with caller-specified CIGAR and
    /// alignment start per record. Each record is mapped, primary, mapq=30,
    /// flags=empty. Sequence is filled with 'A' sized to the read-consuming
    /// CIGAR ops so noodles is happy.
    fn write_coverage_test_bam(
        path: &std::path::PathBuf,
        contigs: &[(&[u8], usize)],
        records: &[TestCoverageRecord<'_>],
    ) {
        use noodles::core::Position;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let mut builder = sam::Header::builder();
        for &(name, len) in contigs {
            builder = builder.add_reference_sequence(
                name.to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(len).unwrap()),
            );
        }
        let header = builder.build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();

        for &(ref_id, start_1based, cigar_ops) in records {
            let cigar: Cigar = cigar_ops.iter().map(|&(k, l)| Op::new(k, l)).collect();
            let read_len: usize = cigar_ops
                .iter()
                .filter(|(k, _)| {
                    matches!(
                        k,
                        Kind::Match
                            | Kind::SequenceMatch
                            | Kind::SequenceMismatch
                            | Kind::Insertion
                            | Kind::SoftClip
                    )
                })
                .map(|&(_, l)| l)
                .sum();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::empty();
            *record.cigar_mut() = cigar;
            *record.reference_sequence_id_mut() = Some(ref_id);
            *record.alignment_start_mut() = Position::new(start_1based);
            *record.sequence_mut() = Sequence::from(vec![b'A'; read_len].as_slice());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }

    #[test]
    fn test_coverage_observer_single_record() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(&path, &[(b"chr1", 10)], &[(0, 1, &[(Kind::Match, 4)])]);
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let by_contig = cov.into_by_contig();
        let depths = by_contig.get("chr1").expect("chr1 should have depths");
        assert_eq!(&depths[..6], &[1, 1, 1, 1, 0, 0]);
    }

    #[test]
    fn test_coverage_observer_overlapping_reads_add() {
        // Two reads overlapping at ref positions 2..4 (0-based) should yield depth 2 there.
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(
            &path,
            &[(b"chr1", 10)],
            &[
                (0, 1, &[(Kind::Match, 4)]), // covers 0..4
                (0, 3, &[(Kind::Match, 4)]), // covers 2..6
            ],
        );
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let depths = cov.into_by_contig().remove("chr1").unwrap();
        assert_eq!(&depths[..7], &[1, 1, 2, 2, 1, 1, 0]);
    }

    #[test]
    fn test_coverage_observer_insertion_does_not_advance_reference() {
        // CIGAR 2M2I2M at pos=1: ref positions 0,1,4,5? No — insertion consumes read
        // only, so ref goes 0,1 (first 2M) → 0,1 (skip insertion) → 2,3 (second 2M).
        // Result: depth 1 at positions 0..4.
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(
            &path,
            &[(b"chr1", 10)],
            &[(
                0,
                1,
                &[(Kind::Match, 2), (Kind::Insertion, 2), (Kind::Match, 2)],
            )],
        );
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let depths = cov.into_by_contig().remove("chr1").unwrap();
        assert_eq!(&depths[..6], &[1, 1, 1, 1, 0, 0]);
    }

    #[test]
    fn test_coverage_observer_deletion_advances_reference_without_depth() {
        // CIGAR 2M2D2M at pos=1: ref goes 0,1 → skip 2,3 (deletion) → cover 4,5.
        // Result: depth 1 at positions 0,1,4,5; 0 at 2,3.
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(
            &path,
            &[(b"chr1", 10)],
            &[(
                0,
                1,
                &[(Kind::Match, 2), (Kind::Deletion, 2), (Kind::Match, 2)],
            )],
        );
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let depths = cov.into_by_contig().remove("chr1").unwrap();
        assert_eq!(&depths[..7], &[1, 1, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn test_coverage_observer_soft_clip_does_not_count() {
        // CIGAR 2S4M at pos=1: soft-clipped bases don't consume reference, so
        // the 4 matches cover ref positions 0..4.
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(
            &path,
            &[(b"chr1", 10)],
            &[(0, 1, &[(Kind::SoftClip, 2), (Kind::Match, 4)])],
        );
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let depths = cov.into_by_contig().remove("chr1").unwrap();
        assert_eq!(&depths[..6], &[1, 1, 1, 1, 0, 0]);
    }

    #[test]
    fn test_coverage_observer_multi_contig() {
        // Records on different contigs land in different depth arrays;
        // unobserved contigs are absent from the output map.
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(
            &path,
            &[(b"chr1", 10), (b"chr2", 10), (b"chr3", 10)],
            &[(0, 1, &[(Kind::Match, 3)]), (1, 5, &[(Kind::Match, 3)])],
        );
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let by_contig = cov.into_by_contig();
        assert_eq!(by_contig.len(), 2, "chr3 had no records, should be absent");
        assert_eq!(&by_contig.get("chr1").unwrap()[..4], &[1, 1, 1, 0]);
        assert_eq!(
            &by_contig.get("chr2").unwrap()[..8],
            &[0, 0, 0, 0, 1, 1, 1, 0]
        );
    }

    #[test]
    fn test_coverage_observer_clips_at_contig_end() {
        // A read whose match span runs past the contig end is silently clipped,
        // not a panic. CIGAR 6M at pos=8 on a 10-length contig: positions 7..10
        // get incremented, position 10..13 are out of bounds and dropped.
        use noodles::sam::alignment::record::cigar::op::Kind;
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("t.bam");
        write_coverage_test_bam(&path, &[(b"chr1", 10)], &[(0, 8, &[(Kind::Match, 6)])]);
        let mut cov = CoverageObserver::default();
        walk_bam(&path, &BamWalkFilter::for_coverage(), &mut [&mut cov]).unwrap();
        let depths = cov.into_by_contig().remove("chr1").unwrap();
        assert_eq!(depths.len(), 10);
        assert_eq!(&depths[..], &[0, 0, 0, 0, 0, 0, 0, 1, 1, 1]);
    }
}
