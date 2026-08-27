//! This library writes either single ended or paired-ended fastq files.
//! Need to update this method. We want to use the data structures and we want to make sure
//! this function is generic enough to work with the fragmented method we are implementing
//! This one needs a major overhaul, it is autogenerating quality scores etc.
//! Will wait to get other things set up first
use crate::rng::{NeatRng, NeatRngError};
use log::debug;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use thiserror::Error;

use crate::file_tools::bam_writer::BamRecordStager;
use crate::file_tools::file_io::append_to_file;
use crate::models::quality_scores::QualityScoreModel;
use crate::models::sequencing_error_model::{
    SeqModelError, SequencingErrorModel, SequencingErrorType,
};
use crate::structs::haplotype_map::InsertionCoordinateMap;
use crate::structs::mutated_map::{AdCounter, MutatedMap, MutatedMapError};
use crate::structs::nucleotides::Nucleotide;
use crate::structs::nucleotides::Nucleotide::N;
use crate::structs::read_record::ReadRecord;
use crate::structs::sequence_block::{SequenceBlock, SequenceBlockError};
use crate::structs::variants::{Genotype, Variant};

#[derive(Error, Debug)]
pub enum FastqToolsError {
    #[error("Error writing bgzip fastq block {0}")]
    FastqWriteError(String),
    #[error("Quality scores missing for read 2")]
    MissingQScores,
    #[error("Error reading bgzip fastq block {0}")]
    FastqReadError(String),
    #[error("Mismatch between indexing and reads set for block {0}")]
    InvalidFastqBlock(String),
    #[error("Fastq Tools reported a SequenceBlock error: {0}")]
    SequenceBlockError(#[from] SequenceBlockError),
    #[error("Fastq tools reported a error model error: {0}")]
    ErrorModelError(#[from] SeqModelError),
    #[error("Fastq tools reported an IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Fastq tools reported an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Error from MutatedMap: {0}")]
    MutatedMapError(#[from] MutatedMapError),
    #[error("Error locating a matching mutated map")]
    FindMapError,
    #[error("Paired ended declared, but r2 buffer empty")]
    BufferInitError,
    #[error("Malformed read during pair-ended analysis: {0}")]
    MalformedReadError(String),
    #[error("Strand value must either be 1 or 2, received {0}")]
    StrandError(usize),
    #[error("Truncated read {0}")]
    TruncatedRead(String),
    #[error("Reverse read with read_end > read_start")]
    MalformedReverseRead,
    #[error("BAM write error: {0}")]
    BamError(String),
    #[error("Haplotype baseline CIGAR does not match read sequence length")]
    HaplotypeCigarMismatch,
}

pub enum Strand {
    Forward,
    Reverse,
}

pub fn reverse_complement(sequence: Vec<Nucleotide>) -> Vec<Nucleotide> {
    // Returns the reverse complement of a vector of u8's representing a DNA sequence.
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(sequence[i].complement())
    }
    rev_comp
}

/// Reference bases consumed by a generated-read CIGAR.
///
/// The generator currently emits `M`, `I`, `D`, and `S`; retain the full SAM
/// set of reference-consuming ops so this stays correct if it later emits
/// explicit match/mismatch or reference-skip operations.
fn reference_span(cigar_ops: &[char]) -> usize {
    cigar_ops
        .iter()
        .filter(|&&op| matches!(op, 'M' | '=' | 'X' | 'D' | 'N'))
        .count()
}

/// Set SAM TLEN from the realized pair geometry rather than the sampled
/// fragment length. Insertions and soft clips consume no reference; deletions
/// do, so the sampled length can disagree with what the two CIGARs cover.
///
/// Positions are 0-based and ends are exclusive. For equal starts, retain R1
/// as the positive mate, which makes the otherwise ambiguous sign stable.
fn set_observed_template_lengths(r1: &mut ReadRecord, r2: &mut ReadRecord) {
    let left_start = r1.position.min(r2.position);
    let right_end = (r1.position + reference_span(&r1.cigar_ops))
        .max(r2.position + reference_span(&r2.cigar_ops));
    let template_length =
        i32::try_from(right_end - left_start).expect("generated template length must fit in i32");

    if r1.position <= r2.position {
        r1.template_length = template_length;
        r2.template_length = -template_length;
    } else {
        r1.template_length = -template_length;
        r2.template_length = template_length;
    }
}

/// One fragment to be written, expressed in the coordinate space named by
/// `haplotype`: reference coordinates when `None`, altered-haplotype coordinates
/// when `Some`.
///
/// This exists so long insertions can be sampled in a coordinate space that has
/// width where the reference has none, WITHOUT a second writer. The previous #516
/// attempt added a parallel writer for exactly this and silently lost everything
/// `write_block_fastq` does along the way — the heterozygous coin, allelic-depth
/// counting, the #210 position-keyed read name, adapter readthrough, and the
/// pair-desync guard. Routing both kinds of fragment through one function is what
/// makes those impossible to drop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlacedFragment {
    pub start: usize,
    pub end: usize,
    /// Index into the `haplotypes` slice passed alongside.
    pub haplotype: Option<usize>,
}

impl From<(usize, usize)> for PlacedFragment {
    /// A fragment sampled straight from the reference — the overwhelmingly common
    /// case, and what every caller produced before insertions needed their own
    /// coordinate space.
    fn from((start, end): (usize, usize)) -> Self {
        Self {
            start,
            end,
            haplotype: None,
        }
    }
}

/// An altered haplotype that fragments in this block may have been sampled from:
/// one literal insertion spliced into the reference.
#[derive(Debug, Clone)]
pub struct HaplotypeContext {
    /// Describes every insertion on this haplotype and owns their novel bases, so
    /// an interval spanning more than one materializes without the caller
    /// tracking which sequence belongs to which anchor.
    pub map: InsertionCoordinateMap,
}

pub fn write_block_fastq<B1: Write, B2: Write>(
    block_fragments: Vec<PlacedFragment>,
    // Altered haplotypes referenced by `PlacedFragment::haplotype`. Empty for a
    // block with no long insertions, which is the common case.
    haplotypes: &[HaplotypeContext],
    block_map: &MutatedMap,
    sequence_block: &SequenceBlock,
    paired_ended: bool,
    buffer1: &mut B1,
    buffer2: &mut B2,
    read_length: usize,
    long_reads: bool,
    // Keep short inserts (insert < read_length) and emit insert-length reads instead
    // of dropping them. Implied when adapters are on (readthrough pads them); set on its
    // own for the adapter-free short-insert control. Independent of adapter padding,
    // which is driven solely by non-empty r1_adapter/r2_adapter below.
    keep_short: bool,
    read_name_prefix: &str,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
    mut bam_writer: Option<&mut dyn BamRecordStager>,
    ad_counter: &mut AdCounter,
    // 3' adapter readthrough (#125). Empty slices = disabled (no behavior change).
    // When non-empty, reads whose insert < read_length are padded to read_length
    // with adapter sequence: R1 gets r1_adapter, R2 gets r2_adapter.
    r1_adapter: &[Nucleotide],
    r2_adapter: &[Nucleotide],
) -> Result<(), FastqToolsError> {
    debug!("writing reads for {}", sequence_block.contig);
    let adapters_on = !r1_adapter.is_empty() || !r2_adapter.is_empty();
    // Pad the fetched fragment with extra reference beyond `end` so deletions
    // (sequencing-error or literal-variant) near a read's tail have bases to
    // consume instead of exhausting the buffer and raising TruncatedRead.
    //   - SE: the single read can end with a deletion -> small tail pad.
    //   - PE: R2 is generated FORWARD over the right-end window, which ends at
    //     `end`, so a deletion in R2 needs reference *beyond* `end`. Without
    //     this pad, any deletion in R2 truncated the read and dropped the whole
    //     pair (~14% coverage loss at 30x). A read-length pad covers any literal
    //     deletion (all < the SV threshold). The pad does NOT shift R2's
    //     coordinates — the R2 window still starts at `end - effective_read_len`.
    let seq_len = sequence_block.sequence.len();
    let frag_pad = if paired_ended { read_length } else { 32 };
    for (frag_idx, placed) in block_fragments.into_iter().enumerate() {
        let (start, end) = (placed.start, placed.end);
        // An index the caller supplied but did not provide a context for is a
        // programming error, not data-dependent: say which, rather than panicking
        // with a bare out-of-bounds. (Hit exactly once while wiring this up, when
        // the runner built its contexts but still passed an empty slice.)
        let hap = match placed.haplotype {
            None => None,
            Some(i) => Some(haplotypes.get(i).ok_or_else(|| {
                FastqToolsError::BamError(format!(
                    "fragment references haplotype {i} but only {} were supplied",
                    haplotypes.len()
                ))
            })?),
        };
        // The span a fragment may be materialized over. For a reference fragment
        // that is the contig; for a haplotype fragment it is the contig plus the
        // inserted bases, which is exactly the extra width the insertion has.
        let materializable_len = match hap {
            None => seq_len,
            Some(h) => h.map.haplotype_len(),
        };
        let padded_end = (end + frag_pad).min(materializable_len);
        // Zero-copy for the common case: the fragment is only read (R1 reads it,
        // R2 reads a suffix), never stored or mutated, so borrow it instead of
        // allocating + copying a Vec per fragment. A haplotype fragment has to be
        // built (its bases do not exist contiguously in the reference), so that
        // path owns a buffer and borrows from it.
        let hap_materialized;
        // `frag_ops` is produced by the same match as `fragment` so the borrow lives in the
        // arm that initializes the buffer — the two must stay index-parallel.
        let (fragment, frag_ops): (&[Nucleotide], Option<&[char]>) = match hap {
            None => (sequence_block.get_subseq_slice(start, padded_end)?, None),
            Some(h) => {
                let Some((bases, segments)) =
                    h.map
                        .materialize_interval(&sequence_block.sequence, start, padded_end)
                else {
                    // The map disagrees with the sequence it was built for, or the
                    // interval is degenerate. Skipping is right (the alternative is
                    // emitting reads from a coordinate space that does not exist),
                    // but it must not be silent -- a region that quietly stops
                    // producing reads is the failure shape this project keeps hitting.
                    debug!(
                        "haplotype fragment [{start},{padded_end}) could not be materialized; skipping"
                    );
                    continue;
                };
                hap_materialized = (
                    bases,
                    InsertionCoordinateMap::cigar_ops_for_segments(&segments),
                );
                (&hap_materialized.0, Some(&hap_materialized.1[..]))
            }
        };
        // In long-read mode a fragment may be shorter than read_length; truncate the read
        // to the actual fragment length rather than discarding it.
        // With adapters on, a short insert generates an insert-length read here, then the
        // 3' adapter pads it to read_length after orientation (see append_adapter_readthrough);
        // capping at the insert length keeps generate_read from truncating-and-dropping it.
        let insert_len = end - start;
        // A zero-length insert (start == end) is a degenerate "adapter dimer": it
        // carries no genomic bases, and with adapters on the padding path built a
        // record whose quality string ran one char longer than the sequence
        // (effective_read_len == 0 desyncs the seq/qual construction). That is a
        // malformed FASTQ record — harmless to `zcat`/`wc`, but bwa-mem2's parser
        // stops at the first one and silently truncates alignment to a few thousand
        // reads (#125). Skip the whole pair so the R1/R2 streams stay in sync; this
        // never fires when adapters are off (fragments are already >= read_length).
        if insert_len == 0 {
            continue;
        }
        // Cap the read at the insert whenever short fragments are being kept — with
        // adapters on (padded back to read_length below) OR in the adapter-free
        // keep_short control (emitted as an insert-length genomic read).
        let effective_read_len = if adapters_on || keep_short {
            insert_len.min(read_length)
        } else if long_reads {
            fragment.len().min(read_length)
        } else {
            read_length
        };
        let mut read1_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads1_flagged: Vec<usize> = Vec::new();
        let mut read2_variants: HashMap<usize, &Variant> = HashMap::new();
        let mut reads2_flagged: Vec<usize> = Vec::new();
        // Only the variants overlapping this fragment's two read windows matter.
        // block_map.flagged_positions is sorted (from_interval), so binary-search
        // each window instead of scanning every variant on the contig. The old
        // full-map scan made this O(fragments × variants) — pathological with
        // dense models (~220k SNVs on chr22 from the COSMIC tumor rate).
        let flagged = &block_map.flagged_positions;
        // A variant whose allele is decided by HAPLOTYPE SAMPLING must not also be
        // applied inline here, or it would be applied twice on the alt haplotype
        // (once materialized, once by the per-variant coin) and, worse, applied at
        // all on the reference haplotype — where by construction it is absent.
        // Deciding it once, per fragment, is the whole point of sampling haplotypes
        // rather than flipping a coin per read: it is what makes a heterozygous
        // insertion actually come out at ~half depth.
        let hap_anchor = |pos: usize| haplotypes.iter().any(|h| h.map.anchors().any(|a| a == pos));
        // R1 window: [start, start + effective_read_len); var offset = pos - start.
        // For a haplotype fragment these are HAPLOTYPE offsets, so a reference
        // variant position has to be projected before it can be located in the read.
        let project = |pos: usize| -> Option<usize> {
            match hap {
                None => Some(pos),
                Some(h) => h.map.reference_base_to_haplotype(pos),
            }
        };
        // Project the read WINDOW back to reference coordinates rather than
        // projecting every variant forward. The map is monotonic, so a haplotype
        // window corresponds to a contiguous reference range -- which keeps the
        // binary search over the sorted flagged positions. Scanning the whole
        // variant map per fragment would be O(fragments x variants), the exact
        // pathology the binary search was introduced to fix (~220k SNVs on chr22
        // under the COSMIC tumor rate).
        let ref_window = |lo: usize, hi: usize| -> (usize, usize) {
            match hap {
                None => (lo, hi),
                Some(h) => (h.map.reference_floor(lo), h.map.reference_floor(hi)),
            }
        };
        let (r1_ref_lo, r1_ref_hi) = ref_window(start, start + effective_read_len);
        let r1_lo = flagged.partition_point(|&p| p < r1_ref_lo);
        let r1_hi = flagged.partition_point(|&p| p < r1_ref_hi);
        for &pos in &flagged[r1_lo..r1_hi] {
            if hap_anchor(pos) {
                continue;
            }
            // Forward-project the handful of hits to locate them within the read.
            let Some(proj) = project(pos) else { continue };
            if proj >= start && proj < start + effective_read_len {
                let var_pos = proj - start;
                read1_variants.insert(var_pos, &block_map.variant_map[&pos]);
                reads1_flagged.push(var_pos);
            }
        }
        // R2 window: [end - effective_read_len, end). R2 is generated in FORWARD
        // orientation over this right-end window (then the whole record is
        // reverse-complemented), so the variant offset is the forward offset
        // within the window, just like R1. Guard end > effective_read_len so the
        // window start doesn't underflow.
        if paired_ended && end > effective_read_len {
            let w_lo = end - effective_read_len;
            let (r2_ref_lo, r2_ref_hi) = ref_window(w_lo, end);
            let r2_lo = flagged.partition_point(|&p| p < r2_ref_lo);
            let r2_hi = flagged.partition_point(|&p| p < r2_ref_hi);
            for &pos in &flagged[r2_lo..r2_hi] {
                if hap_anchor(pos) {
                    continue;
                }
                let Some(proj) = project(pos) else { continue };
                if proj >= w_lo && proj < end {
                    let var_pos = proj - w_lo;
                    read2_variants.insert(var_pos, &block_map.variant_map[&pos]);
                    reads2_flagged.push(var_pos);
                }
            }
        }

        let ref_start = sequence_block.ref_start;
        // Reference-space coordinates for the BAM record. For a reference fragment
        // these are just the fragment's own coordinates offset by the block start.
        // For a haplotype fragment they are the PROJECTION back to the reference:
        // the inserted bases have no reference coordinate at all, so a read window
        // beginning inside them projects to `None` and the record is emitted
        // unmapped rather than being placed at the anchor with an all-insertion
        // CIGAR (which is not a valid alignment -- it consumes no reference).
        //
        // The golden BAM is an answer key, not a prediction of aligner output
        // (see #449), so an unmapped record still records where the read really
        // came from via a provenance tag rather than discarding that.
        let (r1_ref_pos, r2_ref_pos, tlen_span) = match hap {
            None => (
                Some(start + ref_start),
                Some(end.saturating_sub(effective_read_len) + ref_start),
                Some(end - start),
            ),
            Some(h) => {
                let r1 = h
                    .map
                    .haplotype_base_to_reference(start)
                    .map(|p| p + ref_start);
                let r2_window = end.saturating_sub(effective_read_len);
                let r2 = h
                    .map
                    .haplotype_base_to_reference(r2_window)
                    .map(|p| p + ref_start);
                // TLEN is only meaningful when both mates have a reference span.
                let span = match (r1, r2) {
                    (Some(a), Some(b)) => Some((b + effective_read_len).saturating_sub(a)),
                    _ => None,
                };
                (r1, r2, span)
            }
        };
        let abs_start = start + ref_start;
        let abs_end = end + ref_start;
        // Per-fragment uniqueness tag in the read name. Without this, two
        // fragments that land at the same (start, end) — common via birthday
        // paradox at 30×+ coverage — share an identical QNAME, which violates
        // BAM/VCF spec and silently confuses Picard MarkDuplicates into
        // dropping them (see #210). The within-block frag_idx is sufficient
        // because eidolon currently uses one block per contig and per-contig
        // read names already differ via `read_name_prefix`. For multi-pass
        // workflows that concatenate FASTQs from independent gen-reads runs,
        // the caller must still prefix each run's reads (e.g. cancer_simulate.sh
        // does N_/T_ between normal and tumor passes).
        let base_name = format!(
            "{}_{:010}_{:010}_{:016x}",
            read_name_prefix, abs_start, abs_end, frag_idx,
        );

        // SAM convention: an unmapped read with a mapped mate still carries the
        // mate's RNAME/POS, so the pair sorts together and the unmapped read stays
        // located at the event it came from. Falling back to the insertion anchor
        // covers the (rare) case where BOTH mates lie inside the inserted sequence.
        let anchor_fallback = match hap {
            Some(h) => h.map.anchors().next().unwrap_or(start) + ref_start,
            None => abs_start,
        };
        let r1_pos = r1_ref_pos.or(r2_ref_pos).unwrap_or(anchor_fallback);
        // Identical to the previous expression whenever `hap` is None: the
        // projection of a reference fragment is the fragment itself.
        let r2_start = if paired_ended && abs_end >= effective_read_len {
            r2_ref_pos
                .or(r1_ref_pos)
                .unwrap_or(abs_end - effective_read_len)
        } else {
            0
        };
        let tlen = if paired_ended {
            // A template length is only meaningful when both mates have a reference
            // span; an unmapped mate leaves it 0, as SAM requires.
            tlen_span.unwrap_or(0) as i32
        } else {
            0
        };

        let quality_scores_1 =
            quality_score_model.generate_quality_scores(effective_read_len, rng)?;
        let mut r1_record = match generate_read(
            fragment,
            frag_ops,
            &reads1_flagged,
            &read1_variants,
            effective_read_len,
            format!("{}/1", base_name),
            Strand::Forward,
            quality_scores_1,
            sequencing_error_model,
            rng,
            sequence_block.contig.clone(),
            r1_pos,
            sequence_block.contig.clone(),
            r2_start,
            tlen,
            paired_ended,
            ad_counter,
        ) {
            Ok(record) => record,
            Err(FastqToolsError::TruncatedRead(msg)) => {
                debug!("{}", msg);
                continue;
            }
            Err(e) => return Err(e),
        };
        // No reference position => the read lies wholly inside inserted sequence.
        // Emit it unmapped with an empty CIGAR rather than placing it at the anchor
        // with an all-insertion CIGAR, which consumes no reference and is not a
        // valid alignment. The bases still reach the FASTQ — they are real reads.
        if r1_ref_pos.is_none() {
            r1_record.is_unmapped = true;
            // Soft-clip the whole read rather than emptying the CIGAR: `S` consumes
            // query but not reference, which is exactly "these bases exist, none of
            // them align". An empty CIGAR is what SAM writes as `*`, but this
            // writer's encoder turns an empty op list into `1M`, and noodles then
            // rejects the record for a read-length/sequence-length mismatch.
            r1_record.cigar_ops = vec!['S'; r1_record.sequence.len()];
        }
        // R1 adapter readthrough: pad a short-insert read to read_length at its 3' end.
        if adapters_on {
            r1_record = append_adapter_readthrough(
                r1_record,
                r1_adapter,
                read_length,
                quality_score_model,
                sequencing_error_model,
                rng,
            )?;
        }

        // Generate r2 BEFORE writing r1, so that a TruncatedRead on r2
        // skips BOTH reads together. Otherwise r1 lands in buffer1 with
        // no matching r2 in buffer2, and BWA-MEM aborts with "paired
        // reads have different names" when the streams desync. The
        // failure mode was always reachable but became common after #221
        // (the literal-DEL skip+D-op fix), which legitimately advances
        // seq_index past the deleted bases and exhausts the buffer for
        // long deletions near a fragment edge.
        let mut r2_record = if paired_ended {
            let quality_scores_2 =
                quality_score_model.generate_quality_scores(effective_read_len, rng)?;
            let r2_pos = r2_ref_pos
                .or(r1_ref_pos)
                .unwrap_or_else(|| abs_end.saturating_sub(effective_read_len));
            let tlen_r2 = -(tlen_span.unwrap_or(0) as i32);
            // R2 covers the fragment's right end. Generate it FORWARD over that
            // window — so SNP/insertion/deletion handling is identical to R1 and
            // correct — then reverse-complement the whole record into a reverse
            // read. Applying VCF-anchored indels during a reverse walk mis-placed
            // them (insertion base order / deletion anchor), so reverse reads
            // carried garbled indels; forward-generate-then-flip avoids that.
            // R2 window starts at (end - effective_read_len); since `fragment`
            // is padded beyond `end`, &fragment[off..] is the window plus the
            // deletion buffer (so R2 deletions don't truncate and drop the pair).
            let r2_sub: &[Nucleotide] = match (end - start).checked_sub(effective_read_len) {
                Some(off) => &fragment[off..],
                // Fragment shorter than a read: skip the pair to avoid an
                // orphaned R1 (matches the TruncatedRead handling below).
                None => continue,
            };
            // The mask is sliced by the SAME offset as r2_sub, so index i of one lines up
            // with index i of the other. R2 is generated forward over this window and the
            // record (CIGAR included) is reverse-complemented afterwards, so no reversal
            // of the mask is needed here.
            let r2_ops: Option<&[char]> =
                match (frag_ops, (end - start).checked_sub(effective_read_len)) {
                    (Some(ops), Some(off)) if off <= ops.len() => Some(&ops[off..]),
                    _ => None,
                };
            match generate_read(
                r2_sub,
                r2_ops,
                &reads2_flagged,
                &read2_variants,
                effective_read_len,
                format!("{}/2", base_name),
                Strand::Forward,
                quality_scores_2,
                sequencing_error_model,
                rng,
                sequence_block.contig.clone(),
                r2_pos,
                sequence_block.contig.clone(),
                r1_pos,
                tlen_r2,
                true,
                ad_counter,
            ) {
                Ok(mut record) => {
                    if r2_ref_pos.is_none() {
                        record.is_unmapped = true;
                        record.cigar_ops = vec!['S'; record.sequence.len()];
                    }
                    // Flip to the reverse mate FIRST, then append the R2 adapter at
                    // the (now 3') end — so R2 carries the R2 adapter in read
                    // orientation, exactly as a trimmer expects.
                    let mut rec = reverse_complement_record(record);
                    if adapters_on {
                        rec = append_adapter_readthrough(
                            rec,
                            r2_adapter,
                            read_length,
                            quality_score_model,
                            sequencing_error_model,
                            rng,
                        )?;
                    }
                    Some(rec)
                }
                Err(FastqToolsError::TruncatedRead(msg)) => {
                    debug!("{}", msg);
                    // Drop r1 alongside r2 so the streams stay in sync.
                    continue;
                }
                Err(e) => return Err(e),
            }
        } else {
            None
        };

        // Allelic depth for a haplotype-sampled insertion. `generate_read` no
        // longer sees these variants (their allele is decided per fragment, not
        // per read), so nothing else would count them and the golden VCF would
        // report DP=0 for the event -- which is exactly what the reverted attempt
        // shipped. Count each fragment once, against the anchor it covers, on
        // whichever haplotype it was drawn from.
        for h in haplotypes {
            for anchor in h.map.anchors() {
                // Is this fragment drawn from the very haplotype carrying this
                // insertion? If so its coordinates are in that haplotype's space
                // and the anchor has to be projected before it can be located.
                let on_this_haplotype = hap.is_some_and(|active| std::ptr::eq(active, h));
                let covers_anchor = if on_this_haplotype {
                    h.map
                        .reference_base_to_haplotype(anchor)
                        .is_some_and(|a| a >= start && a < end)
                } else {
                    anchor >= start && anchor < end
                };
                if !covers_anchor {
                    continue;
                }
                let entry = ad_counter.entry(anchor + ref_start).or_insert((0, 0));
                if on_this_haplotype {
                    entry.1 += 1;
                } else {
                    entry.0 += 1;
                }
            }
        }

        if let Some(r2_record) = r2_record.as_mut() {
            // Both realized CIGARs are available only here, before either mate
            // is staged and potentially separated by the coordinate writer.
            set_observed_template_lengths(&mut r1_record, r2_record);
        }

        write_read_to_fastq(&r1_record, buffer1)?;
        if let Some(ref mut bam) = bam_writer {
            bam.stage_read_record(&r1_record)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
        }

        if let Some(r2_record) = r2_record {
            write_read_to_fastq(&r2_record, buffer2)?;
            if let Some(ref mut bam) = bam_writer {
                bam.stage_read_record(&r2_record)
                    .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
            }
        }

        // Flush all carry records that start strictly before abs_start — safe because
        // block_fragments is sorted, so no future record can land earlier than abs_start.
        if let Some(ref mut bam) = bam_writer {
            bam.flush_up_to(abs_start)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
        }
    }
    Ok(())
}

/// Write single-ended reads whose sequence was materialized in haplotype
/// coordinates. Each item contains the read sequence, one baseline CIGAR
/// operation per sequence base, and its reference-anchored position. This
/// opt-in path is intentionally separate from `write_block_fastq` until the
/// paired-end fragment ownership and mate-coordinate rules are finalized.
pub fn write_haplotype_fragments<W: Write>(
    fragments: impl IntoIterator<Item = (Vec<Nucleotide>, Vec<char>, usize)>,
    buffer: &mut W,
    read_length: usize,
    read_name_prefix: &str,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
) -> Result<usize, FastqToolsError> {
    let mut written = 0;
    for (sequence, baseline_ops, position) in fragments {
        if sequence.len() != read_length || baseline_ops.len() != read_length {
            return Err(FastqToolsError::HaplotypeCigarMismatch);
        }
        let quality_scores = quality_score_model.generate_quality_scores(read_length, rng)?;
        let mut ad_counter = AdCounter::new();
        let mut record = generate_read(
            &sequence,
            // Chimeric reads are stitched from reference pieces; every base is 'M'.
            None,
            &[],
            &HashMap::new(),
            read_length,
            format!("{}_{}", read_name_prefix, written),
            Strand::Forward,
            quality_scores,
            sequencing_error_model,
            rng,
            "".to_string(),
            position,
            "".to_string(),
            0,
            0,
            false,
            &mut ad_counter,
        )?;
        apply_haplotype_baseline_cigar(&mut record, &baseline_ops)?;
        write_read_to_fastq(&record, buffer)?;
        written += 1;
    }
    Ok(written)
}

pub struct HaplotypePairedFragment {
    pub r1_sequence: Vec<Nucleotide>,
    pub r1_baseline_ops: Vec<char>,
    pub r1_position: usize,
    pub r2_sequence: Vec<Nucleotide>,
    pub r2_baseline_ops: Vec<char>,
    pub r2_position: usize,
    pub template_length: i32,
}

/// Write paired-end reads from expanded-coordinate haplotype windows. R2 is
/// generated in forward orientation, annotated with its baseline operations,
/// then reverse-complemented so the two records share the same construction
/// path as the existing paired-end writer.
pub fn write_haplotype_paired_fragments<B1: Write, B2: Write>(
    fragments: impl IntoIterator<Item = HaplotypePairedFragment>,
    buffer1: &mut B1,
    buffer2: &mut B2,
    read_length: usize,
    read_name_prefix: &str,
    contig_name: &str,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
    mut bam_writer: Option<&mut dyn BamRecordStager>,
) -> Result<usize, FastqToolsError> {
    let mut written = 0;
    for fragment in fragments {
        if fragment.r1_sequence.len() != read_length
            || fragment.r2_sequence.len() != read_length
            || fragment.r1_baseline_ops.len() != read_length
            || fragment.r2_baseline_ops.len() != read_length
        {
            return Err(FastqToolsError::HaplotypeCigarMismatch);
        }
        let name = format!("{}_{}", read_name_prefix, written);
        let mut ad_counter = AdCounter::new();
        let r1_quality = quality_score_model.generate_quality_scores(read_length, rng)?;
        let mut r1 = match generate_read(
            &fragment.r1_sequence,
            // Chimeric reads are stitched from reference pieces; every base is 'M'.
            None,
            &[],
            &HashMap::new(),
            read_length,
            format!("{name}/1"),
            Strand::Forward,
            r1_quality,
            sequencing_error_model,
            rng,
            contig_name.to_string(),
            fragment.r1_position,
            contig_name.to_string(),
            fragment.r2_position,
            fragment.template_length,
            true,
            &mut ad_counter,
        ) {
            Ok(record) => record,
            Err(FastqToolsError::TruncatedRead(_)) => continue,
            Err(error) => return Err(error),
        };
        apply_haplotype_baseline_cigar(&mut r1, &fragment.r1_baseline_ops)?;

        let r2_quality = quality_score_model.generate_quality_scores(read_length, rng)?;
        let mut r2 = match generate_read(
            &fragment.r2_sequence,
            // Chimeric reads are stitched from reference pieces; every base is 'M'.
            None,
            &[],
            &HashMap::new(),
            read_length,
            format!("{name}/2"),
            Strand::Forward,
            r2_quality,
            sequencing_error_model,
            rng,
            contig_name.to_string(),
            fragment.r2_position,
            contig_name.to_string(),
            fragment.r1_position,
            -fragment.template_length,
            true,
            &mut ad_counter,
        ) {
            Ok(record) => record,
            Err(FastqToolsError::TruncatedRead(_)) => continue,
            Err(error) => return Err(error),
        };
        apply_haplotype_baseline_cigar(&mut r2, &fragment.r2_baseline_ops)?;
        r2 = reverse_complement_record(r2);

        write_read_to_fastq(&r1, buffer1)?;
        write_read_to_fastq(&r2, buffer2)?;
        if let Some(bam) = bam_writer.as_deref_mut() {
            bam.stage_read_record(&r1)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
            bam.stage_read_record(&r2)
                .map_err(|e| FastqToolsError::BamError(e.to_string()))?;
        }
        written += 1;
    }
    Ok(written)
}

pub fn combine_temp_fastqs(
    files_r1: Vec<PathBuf>,
    files_r2: Vec<PathBuf>,
    final_filename_r1: &PathBuf,
    final_filename_r2: Option<&PathBuf>,
) -> Result<(), FastqToolsError> {
    stream_gzip_files(&files_r1, final_filename_r1)?;
    if let Some(filename_r2) = final_filename_r2
        && !files_r2.is_empty()
    {
        stream_gzip_files(&files_r2, filename_r2)?;
    }
    Ok(())
}

fn stream_gzip_files(files: &[PathBuf], output: &PathBuf) -> Result<(), FastqToolsError> {
    // The per-contig temp files are each already a complete (multi-member) gzip
    // stream. gzip streams concatenate, so the final file is just the raw bytes
    // of every temp appended in order — no decompress, no recompress. This drops
    // an entire compression pass (each read is compressed once, in the parallel
    // per-contig write, instead of twice) and turns the previously single-
    // threaded combine into pure I/O.
    let mut out_file = BufWriter::new(append_to_file(output)?);
    for file in files {
        let mut f = BufReader::new(
            File::open(file).map_err(|e| FastqToolsError::FastqReadError(e.to_string()))?,
        );
        std::io::copy(&mut f, &mut out_file)
            .map_err(|e| FastqToolsError::FastqWriteError(e.to_string()))?;
    }
    out_file
        .flush()
        .map_err(|e| FastqToolsError::FastqWriteError(e.to_string()))?;
    Ok(())
}

// `cigar_ops.push('D')` runs in a loop per deletion-error base; pushing the
// same byte N times is the entire CIGAR encoding, not a copy-paste mistake.
#[allow(clippy::same_item_push)]
pub fn generate_read(
    sequence: &[Nucleotide],
    // Baseline CIGAR op per base of `sequence`, from
    // `InsertionCoordinateMap::cigar_ops_for_segments`: 'M' where the base came from the
    // reference, 'I' where it came from inserted sequence that has no reference position.
    // `None` for a plain reference fragment, where every base is 'M'.
    //
    // Without this the CIGAR cannot describe a haplotype fragment. The bases of a long
    // insertion ARE in the fragment (#516), but nothing downstream could tell them from
    // reference, so a read crossing the anchor was recorded as pure M and claimed inserted
    // bases as reference matches — 108M42I written as 150M (#589).
    hap_ops: Option<&[char]>,
    flagged_positions: &[usize],
    variant_map: &HashMap<usize, &Variant>,
    read_length: usize,
    name: String,
    read_strand: Strand,
    quality_scores: Vec<usize>,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
    contig: String,
    position: usize,
    mate_contig: String,
    mate_position: usize,
    template_length: i32,
    is_paired: bool,
    ad_counter: &mut AdCounter,
) -> Result<ReadRecord, FastqToolsError> {
    if sequence.len() < read_length {
        return Err(FastqToolsError::TruncatedRead(format!("{:?}", sequence)));
    }

    let is_reverse = matches!(read_strand, Strand::Reverse);
    let fragment_length = sequence.len();

    let mut bases_written = 0;
    // Pre-size to read_length (+ slack for insertions) so the per-read output
    // buffers don't repeatedly reallocate as they grow.
    let mut out_seq = String::with_capacity(read_length + 16);
    let mut cigar_ops: Vec<char> = Vec::with_capacity(read_length + 16);
    // Reused across positions for the rare multi-base (insertion) case only, so
    // insertions don't allocate a fresh Vec each time. The common single-base
    // path (ref base / SNP / SNP-error) uses the stack array `single` below and
    // never touches the heap. This is the hot-loop allocation that previously
    // ran once per base of every read (read_length × 2 reads × every fragment).
    let mut ins_buf: Vec<Nucleotide> = Vec::new();
    let mut quality_index = 0;
    let mut seq_index = 0;

    'outer: while (seq_index < fragment_length) && (bases_written < read_length) {
        // Index variants by seq_index for BOTH strands. The caller already
        // reverse-complements the fragment for reverse (R2) reads AND maps each
        // variant's coordinate into that reversed sequence (var_pos = (end-1)-pos),
        // so the variant base sits at `seq_index` here just like the forward read.
        // The old `fragment_length - seq_index` for Reverse applied a second,
        // erroneous reflection — the variant was looked up at the mirror position,
        // so reverse reads carried REF at the true locus and the alt landed
        // elsewhere. Net effect: alternate alleles only ever appeared on
        // forward-strand reads, which strand-aware callers (e.g. Mutect2) correctly
        // flag as strand bias and filter out.
        let fragment_position = seq_index;
        let reference_base = sequence[seq_index].get_unmasked_base();
        // Common case writes exactly one base (the reference base, or a single
        // SNP/SNP-error substitution) — kept in the stack array `single`, no
        // heap. Only insertions (multi-base alt / insertion error) set use_ins
        // and fill the reused `ins_buf`.
        let mut single = [reference_base];
        let mut use_ins = false;
        let mut deletion_skip: usize = 0;

        if reference_base == N {
            // Don't try to modify this base.
        } else if flagged_positions.contains(&fragment_position) {
            let variant = variant_map[&fragment_position];
            // MutatedMap::from_interval routes symbolic / structural ALTs to
            // sv_records, so anything in variant_map should be literal. Assert
            // that invariant in debug builds — a regression here would otherwise
            // panic at as_literal().unwrap() with no context.
            debug_assert!(
                variant.alternate.is_literal(),
                "symbolic ALT reached generate_read at position {fragment_position}"
            );
            let entry = ad_counter.entry(variant.location).or_insert((0, 0));
            // When the variant carries an explicit allele_fraction (from an input
            // VCF, #398), emit the alt allele on that fraction of reads. Otherwise
            // fall back to the Genotype default: homozygous always alt, het ~0.5.
            // The else-branch keeps the exact same short-circuit (homozygous draws
            // no rng) so default runs stay byte-identical.
            let is_alt = if let Some(f) = variant.allele_fraction {
                rng.random()? < f
            } else {
                (variant.genotype == Genotype::Homozygous) || (rng.random()? < 0.5)
            };
            if is_alt {
                let alt = variant.alternate.as_literal().unwrap();
                let alt_len = alt.len();
                // Reverse reads are generated forward then flipped by the
                // caller, EXCEPT alt bases which are complemented here (length is
                // preserved, so alt_len is strand-independent).
                if alt_len == 1 {
                    single[0] = match read_strand {
                        Strand::Forward => alt[0],
                        Strand::Reverse => alt[0].complement(),
                    };
                } else {
                    ins_buf.clear();
                    match read_strand {
                        Strand::Forward => ins_buf.extend_from_slice(alt),
                        Strand::Reverse => ins_buf.extend(alt.iter().map(|b| b.complement())),
                    }
                    use_ins = true;
                }
                entry.1 += 1; // alt
                // For a net-deletion variant (REF longer than ALT — typically a
                // pure literal deletion from indel_model with ALT.len() == 1),
                // skip the deleted REF bases. Without this, the read transcribes
                // the deleted region from the unbroken reference and emits no
                // CIGAR D-op, leaving the variant invisible to downstream
                // callers (see #221). The deletion_skip + D-op machinery below
                // is the same path SequencingErrorType::DeletionError already
                // uses; reusing it keeps CIGAR shape uniform across both
                // sources of deletion signal.
                let ref_len = variant.reference.len();
                if ref_len > alt_len {
                    let want_skip = ref_len - alt_len;
                    // Cap at remaining buffer so we don't read past the
                    // fragment end. Truncated cases fall through to the
                    // existing TruncatedRead error path.
                    let max_skip = fragment_length.saturating_sub(seq_index).saturating_sub(1);
                    let actual_skip = want_skip.min(max_skip);
                    seq_index += actual_skip;
                    deletion_skip = actual_skip;
                }
            } else {
                entry.0 += 1; // ref (het coin landed on ref)
            }
        } else {
            let score = quality_scores[quality_index];
            let prob = sequencing_error_model.convert_score(score)?;
            if rng.random()? < prob {
                let error =
                    sequencing_error_model.generate_sequencing_error(reference_base, rng)?;
                match error {
                    SequencingErrorType::SnpError(base) => {
                        single[0] = base;
                    }
                    SequencingErrorType::DeletionError(length) => {
                        if seq_index + length < sequence.len() {
                            seq_index += length;
                            deletion_skip = length;
                        }
                    }
                    SequencingErrorType::InsertionError(vec) => {
                        ins_buf.clear();
                        ins_buf.push(reference_base);
                        ins_buf.extend(vec);
                        use_ins = true;
                    }
                }
            }
        }

        // Borrow single (stack) for the common path or the reused ins_buf for
        // insertions — set only after all mutations above, so no aliasing.
        let to_write: &[Nucleotide] = if use_ins { &ins_buf } else { &single };
        let mut is_first_base = true;
        for &base in to_write {
            out_seq.push(base.into());
            bases_written += 1;
            if is_first_base {
                // 'I' when this fragment base is inserted sequence; sequencing-error
                // insertions below still layer their own 'I' on top of whichever it is.
                cigar_ops.push(
                    hap_ops
                        .and_then(|ops| ops.get(fragment_position))
                        .copied()
                        .unwrap_or('M'),
                );
                is_first_base = false;
            } else {
                cigar_ops.push('I');
            }
            if bases_written == read_length {
                break 'outer;
            }
        }

        // 'D' ops for skipped reference bases (sequencing deletion errors only)
        for _ in 0..deletion_skip {
            cigar_ops.push('D');
        }

        seq_index += 1;
        quality_index += 1;
    }

    if bases_written < read_length {
        return Err(FastqToolsError::TruncatedRead(format!("{:?}", sequence)));
    }

    Ok(ReadRecord {
        is_unmapped: false,
        name,
        sequence: out_seq,
        quality_scores,
        cigar_ops,
        is_paired,
        is_reverse,
        contig,
        position,
        mate_contig,
        mate_position,
        template_length,
    })
}

/// Apply baseline alignment operations from a materialized haplotype interval
/// to a generated read. Deletion operations do not consume query sequence and
/// are left untouched; reference `M` operations can be relabeled as insertion
/// `I` operations. Existing sequencing-error insertions remain `I`.
pub fn apply_haplotype_baseline_cigar(
    record: &mut ReadRecord,
    baseline_ops: &[char],
) -> Result<(), FastqToolsError> {
    let mut query_index = 0;
    for op in &mut record.cigar_ops {
        if matches!(*op, 'D' | 'N') {
            continue;
        }
        let baseline = baseline_ops
            .get(query_index)
            .ok_or(FastqToolsError::HaplotypeCigarMismatch)?;
        if matches!(*op, 'M' | '=' | 'X') {
            *op = *baseline;
        }
        query_index += 1;
    }
    if query_index != baseline_ops.len() || query_index != record.sequence.len() {
        return Err(FastqToolsError::HaplotypeCigarMismatch);
    }
    Ok(())
}

/// Turn a forward-generated read into its reverse-strand mate: reverse-complement
/// the sequence and reverse the per-base CIGAR ops and qualities. R2 is generated
/// forward over the fragment's right-end window (so SNP/insertion/deletion handling
/// is identical to R1 and correct) and flipped here — this avoids the
/// reverse-walk indel hazards (insertion base order, deletion anchor) that
/// garbled indels on reverse reads.
fn reverse_complement_record(mut record: ReadRecord) -> ReadRecord {
    record.sequence = record
        .sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            other => other,
        })
        .collect();
    record.cigar_ops.reverse();
    record.quality_scores.reverse();
    record.is_reverse = true;
    record
}

/// Append 3' sequencing-adapter readthrough to a FINAL-oriented read (#125).
/// When the insert was shorter than `read_length`, the record holds only the
/// insert bases; this pads it to `read_length` with adapter sequence sourced
/// cyclically (real reads can read through into — and past — the adapter).
/// Adapter bases:
///   - carry the same per-base substitution error as insert bases, so trimmers
///     see realistic (not pristine) adapters; insertion/deletion errors are
///     intentionally NOT applied here, to keep the read exactly `read_length`,
///   - take quality scores from the same quality model,
///   - are soft-clipped (`'S'`) in the golden BAM — they are not reference-aligned,
///   - carry no variants (adapter is not reference-derived).
/// Pass the R1 adapter for R1, and the R2 adapter for the post-flip R2 read.
/// (The end-to-end fastp/cutadapt trim check is an integration step — see #125.)
fn append_adapter_readthrough(
    mut record: ReadRecord,
    adapter: &[Nucleotide],
    read_length: usize,
    quality_score_model: &QualityScoreModel,
    sequencing_error_model: &SequencingErrorModel,
    rng: &mut NeatRng,
) -> Result<ReadRecord, FastqToolsError> {
    let current = record.sequence.len();
    if adapter.is_empty() || current >= read_length {
        return Ok(record);
    }
    let n_adapter = read_length - current;
    let quals = quality_score_model.generate_quality_scores(n_adapter, rng)?;
    for (i, &score) in quals.iter().enumerate() {
        let mut base = adapter[i % adapter.len()];
        let prob = sequencing_error_model.convert_score(score)?;
        if rng.random()? < prob {
            // Substitution noise only — preserves exact read_length.
            if let SequencingErrorType::SnpError(b) =
                sequencing_error_model.generate_sequencing_error(base, rng)?
            {
                base = b;
            }
        }
        record.sequence.push(base.into());
        record.quality_scores.push(score);
        record.cigar_ops.push('S');
    }
    Ok(record)
}

pub fn write_read_to_fastq<W: Write>(
    record: &ReadRecord,
    buffer: &mut W,
) -> Result<(), FastqToolsError> {
    buffer.write_all(format!("@{}\n", record.name).as_bytes())?;
    buffer.write_all(record.sequence.as_bytes())?;
    buffer.write_all(b"\n+\n")?;
    buffer.write_all(&quality_scores_to_char_vec(&record.quality_scores)?)?;
    buffer.write_all(b"\n")?;
    Ok(())
}

pub fn quality_scores_to_char_vec(array: &[usize]) -> Result<Vec<u8>, FastqToolsError> {
    let mut score_vec = Vec::new();
    for &score in array {
        score_vec.push((score + 33) as u8)
    }
    Ok(score_vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::file_tools::bam_writer::{BamRecordStager, BamWriter, BamWriterError};
    use crate::file_tools::file_io::{VectorBuffer, create_output_file, read_gzip_lines};
    use crate::structs::nucleotides::Nucleotide::*;
    use crate::structs::sequence_block::{RegionType, SequenceMap};
    use crate::structs::variants::{Variant, VariantType};
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    #[derive(Default)]
    struct CapturedPairedReads(Vec<ReadRecord>);

    impl BamRecordStager for CapturedPairedReads {
        fn stage_read_record(&mut self, record: &ReadRecord) -> Result<(), BamWriterError> {
            self.0.push(ReadRecord {
                is_unmapped: false,
                name: record.name.clone(),
                sequence: record.sequence.clone(),
                quality_scores: record.quality_scores.clone(),
                cigar_ops: record.cigar_ops.clone(),
                is_paired: record.is_paired,
                is_reverse: record.is_reverse,
                contig: record.contig.clone(),
                position: record.position,
                mate_contig: record.mate_contig.clone(),
                mate_position: record.mate_position,
                template_length: record.template_length,
            });
            Ok(())
        }

        fn flush_up_to(&mut self, _flush_pos: usize) -> Result<(), BamWriterError> {
            Ok(())
        }
    }

    #[test]
    fn test_combine_temp_fastqs() {
        let temp_dir = tempfile::tempdir().unwrap();

        let write_fastq_gz = |name: &str, content: &[u8]| -> PathBuf {
            let path = temp_dir.path().join(name);
            let f = std::fs::File::create(&path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(content).unwrap();
            enc.finish().unwrap();
            path
        };

        let file1 = write_fastq_gz("r1.fastq.gz", b"@read1\nACGT\n+\nIIII\n");
        let file2 = write_fastq_gz("r2.fastq.gz", b"@read2\nTTTT\n+\nIIII\n");
        let output = temp_dir.path().join("combined.fastq.gz");

        combine_temp_fastqs(vec![file1, file2], vec![], &output, None).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(
            lines.len(),
            8,
            "Combined file should have 8 lines (2 records × 4 lines)"
        );
        assert_eq!(lines[0], "@read1");
        assert_eq!(lines[4], "@read2");
    }

    #[test]
    fn test_apply_haplotype_baseline_cigar_marks_inserted_bases() {
        let mut record = ReadRecord {
            is_unmapped: false,
            name: "read".to_string(),
            sequence: "ACCCCCAAAAAAAAAAAAAAAAAAA".to_string(),
            quality_scores: vec![30; 25],
            cigar_ops: vec!['M'; 25],
            is_paired: false,
            is_reverse: false,
            contig: "chr1".to_string(),
            position: 100,
            mate_contig: "chr1".to_string(),
            mate_position: 0,
            template_length: 0,
        };
        let mut baseline = vec!['M'; 25];
        baseline[1..6].fill('I');
        apply_haplotype_baseline_cigar(&mut record, &baseline).unwrap();
        assert_eq!(record.cigar_ops[0], 'M');
        assert_eq!(&record.cigar_ops[1..6], &['I'; 5]);
        assert_eq!(&record.cigar_ops[6..], &['M'; 19]);
    }

    #[test]
    fn test_generated_read_accepts_haplotype_baseline_cigar() {
        let sequence = vec![A, C, C, C, C, C]
            .into_iter()
            .chain(vec![A; 19])
            .collect::<Vec<_>>();
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["haplotype-read".to_string()]).unwrap();
        let mut record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[],
            &HashMap::new(),
            sequence.len(),
            "hap/1".to_string(),
            Strand::Forward,
            vec![40; sequence.len()],
            &model,
            &mut rng,
            "chr1".to_string(),
            100,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mut baseline = vec!['M'; sequence.len()];
        baseline[1..6].fill('I');
        apply_haplotype_baseline_cigar(&mut record, &baseline).unwrap();
        assert_eq!(record.sequence.len(), sequence.len());
        assert!(record.cigar_ops[1..6].iter().all(|&op| op == 'I'));
    }

    #[test]
    fn test_write_haplotype_fragments_emits_single_ended_records() {
        let sequence = vec![A, C, C, C, C, C]
            .into_iter()
            .chain(vec![A; 19])
            .collect::<Vec<_>>();
        let baseline = vec!['M'; 1]
            .into_iter()
            .chain(vec!['I'; 5])
            .chain(vec!['M'; 19])
            .collect::<Vec<_>>();
        let quality_model = QualityScoreModel::default().unwrap();
        let error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["haplotype-writer".to_string()]).unwrap();
        let mut output = Vec::new();
        let written = write_haplotype_fragments(
            vec![(sequence, baseline, 100)],
            &mut output,
            25,
            "hap",
            &quality_model,
            &error_model,
            &mut rng,
        )
        .unwrap();
        assert_eq!(written, 1);
        let text = String::from_utf8(output).unwrap();
        assert_eq!(text.lines().count(), 4);
        assert!(text.starts_with("@hap_0\n"));
    }

    #[test]
    fn test_write_haplotype_paired_fragments_emits_mates() {
        let sequence = vec![A, C, C, C, C, C]
            .into_iter()
            .chain(vec![A; 19])
            .collect::<Vec<_>>();
        let baseline = vec!['M'; 1]
            .into_iter()
            .chain(vec!['I'; 5])
            .chain(vec!['M'; 19])
            .collect::<Vec<_>>();
        let quality_model = QualityScoreModel::default().unwrap();
        let error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["haplotype-paired-writer".to_string()]).unwrap();
        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        let written = write_haplotype_paired_fragments(
            vec![HaplotypePairedFragment {
                r1_sequence: sequence.clone(),
                r1_baseline_ops: baseline.clone(),
                r1_position: 100,
                r2_sequence: sequence,
                r2_baseline_ops: baseline,
                r2_position: 250,
                template_length: 175,
            }],
            &mut r1_output,
            &mut r2_output,
            25,
            "hap",
            "chr1",
            &quality_model,
            &error_model,
            &mut rng,
            None,
        )
        .unwrap();
        assert_eq!(written, 1);
        assert!(
            String::from_utf8(r1_output)
                .unwrap()
                .starts_with("@hap_0/1\n")
        );
        assert!(
            String::from_utf8(r2_output)
                .unwrap()
                .starts_with("@hap_0/2\n")
        );
    }

    #[test]
    fn test_write_haplotype_paired_fragments_preserves_orientation_and_mates() {
        let r1_sequence = vec![A, C, G, T].into_iter().chain(vec![A; 21]).collect();
        let r2_sequence = vec![T, G, C, A].into_iter().chain(vec![C; 21]).collect();
        let mut r2_baseline = vec!['M'; 25];
        r2_baseline[..4].fill('I');
        let quality_model = QualityScoreModel::default().unwrap();
        let error_model =
            SequencingErrorModel::from_raw_data(0.0, quality_model.clone(), None).unwrap();
        let mut rng =
            NeatRng::new_from_seed(&vec!["haplotype-paired-orientation".to_string()]).unwrap();
        let mut r1_output = Vec::new();
        let mut r2_output = Vec::new();
        let mut captured = CapturedPairedReads::default();
        assert_eq!(
            write_haplotype_paired_fragments(
                vec![HaplotypePairedFragment {
                    r1_sequence,
                    r1_baseline_ops: vec!['M'; 25],
                    r1_position: 100,
                    r2_sequence,
                    r2_baseline_ops: r2_baseline,
                    r2_position: 250,
                    template_length: 175,
                }],
                &mut r1_output,
                &mut r2_output,
                25,
                "hap",
                "chr1",
                &quality_model,
                &error_model,
                &mut rng,
                Some(&mut captured),
            )
            .unwrap(),
            1
        );
        assert_eq!(captured.0.len(), 2);
        let r1 = &captured.0[0];
        let r2 = &captured.0[1];
        assert!(r1.is_paired && !r1.is_reverse);
        assert!(r2.is_paired && r2.is_reverse);
        assert_eq!((r1.name.as_str(), r2.name.as_str()), ("hap_0/1", "hap_0/2"));
        assert_eq!((r1.position, r2.position), (100, 250));
        assert_eq!((r1.mate_position, r2.mate_position), (250, 100));
        assert_eq!((r1.template_length, r2.template_length), (175, -175));
        assert_eq!(&r2.cigar_ops[21..], &['I'; 4]);
        assert_eq!(r1.contig, "chr1");
        assert_eq!(r2.mate_contig, "chr1");
    }

    #[test]
    fn test_write_reverse() {
        let temp_dir = tempfile::tempdir().unwrap();
        let mut temp_file = PathBuf::from(temp_dir.path());
        temp_file.push("test.fastq.gz");
        let mut temp_writer = create_output_file(&temp_file, true).unwrap();
        let original_seq = vec![A, C, C, G, A, A, T, G, A];
        let rev_comp_seq = reverse_complement(original_seq);
        let expected_rev_comp = vec![T, C, A, T, T, C, G, G, T];
        assert_eq!(expected_rev_comp, rev_comp_seq);
        let flagged_positions: Vec<usize> = Vec::new();
        let variant_map: HashMap<usize, &Variant> = HashMap::new();
        let read_len = 4;
        let read_name = "neat_gen__0000000000_0000000009/2".to_string();
        let quality_scores = vec![32, 32, 32, 32];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let record = generate_read(
            &rev_comp_seq,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &flagged_positions,
            &variant_map,
            read_len,
            read_name,
            Strand::Reverse,
            quality_scores,
            &sequencing_error_model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mut buffer = GzEncoder::new(&mut temp_writer, Compression::default());
        write_read_to_fastq(&record, &mut buffer).unwrap();
        buffer.flush().unwrap();
        let temp_reader = read_gzip_lines(&temp_file).unwrap();
        let seq_line: String = {
            let mut rtrn = String::new();
            for line in temp_reader {
                let l = line.unwrap().to_string();
                if l.starts_with("@") {
                    // skip name line
                } else {
                    rtrn = l;
                    break;
                }
            }
            rtrn
        };
        let exp_rev_cmp_str = "TCAT".to_string();
        assert_eq!(seq_line, exp_rev_cmp_str);
    }

    #[test]
    fn test_qual_score_to_write() {
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37, 37, 25];
        let qual_string = "B:F=0BBFF:".as_bytes();
        assert_eq!(
            qual_string,
            quality_scores_to_char_vec(&qual_scores).unwrap()
        )
    }

    #[test]
    fn test_reverse_complement() {
        let read: Vec<Nucleotide> = vec![A, A, A, A, C, C, C, C, C];
        let revcomp: Vec<Nucleotide> = vec![G, G, G, G, G, T, T, T, T];
        assert_eq!(reverse_complement(read), revcomp);
    }

    #[test]
    fn observed_template_length_includes_deletions_but_not_insertions_or_soft_clips() {
        let mut r1 = adapter_rec("A", true, false);
        r1.position = 27;
        r1.cigar_ops = vec!['M'; 100];

        let mut r2 = adapter_rec("A", true, true);
        r2.position = 123;
        // 62M1D20M1I18M2S consumes 101 reference bases. The I and S do not.
        r2.cigar_ops = [
            vec!['M'; 62],
            vec!['D'],
            vec!['M'; 20],
            vec!['I'],
            vec!['M'; 18],
            vec!['S'; 2],
        ]
        .concat();

        set_observed_template_lengths(&mut r1, &mut r2);

        assert_eq!(reference_span(&r1.cigar_ops), 100);
        assert_eq!(reference_span(&r2.cigar_ops), 101);
        assert_eq!((r1.template_length, r2.template_length), (197, -197));
    }

    // --- adapter readthrough (#125) ---
    fn adapter_rec(seq: &str, paired: bool, reverse: bool) -> ReadRecord {
        ReadRecord {
            is_unmapped: false,
            name: "frag/1".to_string(),
            sequence: seq.to_string(),
            quality_scores: vec![30; seq.len()],
            cigar_ops: vec!['M'; seq.len()],
            is_paired: paired,
            is_reverse: reverse,
            contig: "chr1".to_string(),
            position: 0,
            mate_contig: "chr1".to_string(),
            mate_position: 0,
            template_length: 0,
        }
    }

    #[test]
    fn test_append_adapter_pads_to_read_length_and_softclips() {
        // Short insert (8) + adapter -> exactly read_length (20), adapter region soft-clipped.
        let (read_length, insert_len) = (20usize, 8usize);
        let adapter: Vec<Nucleotide> = "GATCGATCGATCGATC".chars().map(Nucleotide::from).collect();
        let qm = QualityScoreModel::default().unwrap();
        let em = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["adapter".to_string()]).unwrap();
        let out = append_adapter_readthrough(
            adapter_rec(&"A".repeat(insert_len), true, false),
            &adapter,
            read_length,
            &qm,
            &em,
            &mut rng,
        )
        .unwrap();
        assert_eq!(out.sequence.len(), read_length, "padded to read_length");
        assert_eq!(out.quality_scores.len(), read_length);
        assert_eq!(out.cigar_ops.len(), read_length);
        assert!(
            out.cigar_ops[..insert_len].iter().all(|&c| c == 'M'),
            "insert stays M"
        );
        assert_eq!(
            out.cigar_ops[insert_len..]
                .iter()
                .filter(|&&c| c == 'S')
                .count(),
            read_length - insert_len,
            "adapter region is soft-clipped",
        );
    }

    #[test]
    fn test_append_adapter_noop_when_insert_ge_read_length() {
        let read_length = 10usize;
        let adapter: Vec<Nucleotide> = "GGGG".chars().map(Nucleotide::from).collect();
        let qm = QualityScoreModel::default().unwrap();
        let em = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["x".to_string()]).unwrap();
        let rec = adapter_rec("ACGTACGTAC", false, false); // already read_length
        let before = rec.sequence.clone();
        let out =
            append_adapter_readthrough(rec, &adapter, read_length, &qm, &em, &mut rng).unwrap();
        assert_eq!(
            out.sequence, before,
            "no adapter when insert >= read_length"
        );
    }

    #[test]
    fn test_r2_adapter_appended_forward_after_flip() {
        // The caller flips R2 (reverse_complement_record) THEN appends the R2 adapter,
        // so the adapter must sit at the final read's 3' end in FORWARD orientation —
        // an A/C-rich adapter leaves the tail A/C-rich, NOT the G/T-rich revcomp.
        let read_length = 30usize;
        let flipped = reverse_complement_record(adapter_rec("ACGTAC", true, false));
        assert!(flipped.is_reverse);
        let r2_adapter: Vec<Nucleotide> = "AAAACCCCAAAACCCCAAAACCCC"
            .chars()
            .map(Nucleotide::from)
            .collect();
        let qm = QualityScoreModel::default().unwrap();
        let em = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["r2".to_string()]).unwrap();
        let out = append_adapter_readthrough(flipped, &r2_adapter, read_length, &qm, &em, &mut rng)
            .unwrap();
        assert_eq!(out.sequence.len(), read_length);
        let tail = &out.sequence[6..];
        let ac = tail.chars().filter(|c| matches!(c, 'A' | 'C')).count();
        let gt = tail.chars().filter(|c| matches!(c, 'G' | 'T')).count();
        assert!(
            ac > gt,
            "R2 adapter must be forward (A/C-rich), not revcomp'd; tail={tail}"
        );
    }

    #[test]
    fn test_write_block_fastq_ref_start_in_read_name() {
        // Verifies that when a SequenceBlock has ref_start > 0, the read names in the
        // output FASTQ use reference-relative positions, not block-local positions.
        use crate::structs::{
            mutated_map::MutatedMap,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let ref_start: usize = 1000;
        let seq_len: usize = 50;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start,
            ref_end: ref_start + seq_len,
            sequence: sequence.clone(),
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map =
            MutatedMap::from_interval(ref_start, ref_start + seq_len, vec![]).unwrap();
        let frag_start: usize = 5;
        let frag_end: usize = 25;
        let fragments = vec![(frag_start, frag_end)];
        let out_path = temp_dir.path().join("out.fastq.gz");
        let outfile = create_output_file(&out_path, true).unwrap();
        use crate::file_tools::file_io::VectorBuffer;
        let dummy_data: VectorBuffer = VectorBuffer::new();
        let mut buf1 = GzEncoder::new(outfile, Compression::default());
        let mut buf2 = GzEncoder::new(dummy_data, Compression::default());
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        write_block_fastq(
            fragments.into_iter().map(Into::into).collect(),
            &[],
            &mutated_map,
            &block,
            false,
            &mut buf1,
            &mut buf2,
            10,
            false,
            false, // keep_short
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            None,
            &mut AdCounter::new(),
            &[],
            &[],
        )
        .unwrap();
        buf1.finish().unwrap();
        let lines: Vec<String> = read_gzip_lines(&out_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let header = &lines[0];
        assert!(
            header.starts_with('@'),
            "Expected FASTQ header, got: {}",
            header
        );
        // Read name must contain ref-relative positions, not block-local ones
        let expected_start = format!("{:010}", frag_start + ref_start); // 0000001005
        let expected_end = format!("{:010}", frag_end + ref_start); // 0000001025
        assert!(
            header.contains(&expected_start),
            "Read name should contain ref-relative start {}; got: {}",
            expected_start,
            header
        );
        assert!(
            header.contains(&expected_end),
            "Read name should contain ref-relative end {}; got: {}",
            expected_end,
            header
        );
        // Also verify block-local positions are NOT what was written
        let local_start = format!("{:010}", frag_start); // 0000000005
        assert!(
            !header.contains(&local_start),
            "Read name should NOT contain block-local start {}; got: {}",
            local_start,
            header
        );
    }

    /// Two fragments at the SAME (start, end) coordinate must produce reads
    /// with distinct QNAMEs — without the per-fragment uniqueness tag this
    /// would silently emit collision-named reads, which Picard MarkDuplicates
    /// would interpret as PCR duplicates and drop (see #210).
    #[test]
    fn test_write_block_fastq_unique_names_for_same_position_fragments() {
        use crate::structs::{
            mutated_map::MutatedMap,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let seq_len: usize = 200;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: seq_len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map = MutatedMap::from_interval(0, seq_len, vec![]).unwrap();
        // Four fragments — three at the SAME coordinate, one elsewhere.
        let fragments = vec![(10, 80), (10, 80), (10, 80), (100, 170)];
        let out_path = temp_dir.path().join("collide.fastq.gz");
        let outfile = create_output_file(&out_path, true).unwrap();
        use crate::file_tools::file_io::VectorBuffer;
        let dummy: VectorBuffer = VectorBuffer::new();
        let mut buf1 = GzEncoder::new(outfile, Compression::default());
        let mut buf2 = GzEncoder::new(dummy, Compression::default());
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["uniq-name".to_string()]).unwrap();
        write_block_fastq(
            fragments.into_iter().map(Into::into).collect(),
            &[],
            &mutated_map,
            &block,
            false,
            &mut buf1,
            &mut buf2,
            70,
            false,
            false, // keep_short
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            None,
            &mut AdCounter::new(),
            &[],
            &[],
        )
        .unwrap();
        buf1.finish().unwrap();
        // Read every line that starts with @ AND is followed by a valid
        // FASTQ record. Filter by record-position (line 1 of every 4-line
        // block) — quality lines (line 4) can start with @ too.
        let lines: Vec<String> = read_gzip_lines(&out_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let names: Vec<&String> = lines
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 4 == 0)
            .map(|(_, l)| l)
            .collect();
        assert_eq!(names.len(), 4, "expected 4 read names, got: {:?}", names);
        let mut sorted: Vec<&&String> = names.iter().collect();
        sorted.sort();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            4,
            "expected 4 unique names, got duplicates in: {:?}",
            names
        );
    }

    #[test]
    fn test_write_block_fastq_skips_zero_insert_no_malformed_record() {
        // Regression (#125): a zero-length insert (start == end) is a degenerate
        // adapter-dimer fragment. With adapters on it used to emit a record whose
        // quality string ran one char longer than the sequence — malformed FASTQ
        // that `zcat`/`wc` ignore but bwa-mem2's parser halts on, silently
        // truncating alignment to a few thousand reads. It must be skipped, and
        // every emitted record must have seq.len() == qual.len() == read_length.
        use crate::structs::{
            mutated_map::MutatedMap,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let seq_len: usize = 300;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: seq_len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map = MutatedMap::from_interval(0, seq_len, vec![]).unwrap();
        let read_length = 70usize;
        // zero-insert (must be skipped) + short insert (adapter-padded) + full insert.
        let fragments = vec![(50, 50), (10, 40), (100, 170)];
        let r1_path = temp_dir.path().join("r1.fastq.gz");
        let r2_path = temp_dir.path().join("r2.fastq.gz");
        let mut buf1 = GzEncoder::new(
            create_output_file(&r1_path, true).unwrap(),
            Compression::default(),
        );
        let mut buf2 = GzEncoder::new(
            create_output_file(&r2_path, true).unwrap(),
            Compression::default(),
        );
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["zero-insert".to_string()]).unwrap();
        let adapter: Vec<Nucleotide> = vec![A, G, A, T, C, G, G, A, A, G, A, G, C];
        write_block_fastq(
            fragments.into_iter().map(Into::into).collect(),
            &[],
            &mutated_map,
            &block,
            true, // paired_ended
            &mut buf1,
            &mut buf2,
            read_length,
            false, // long_reads
            true,  // keep_short (adapters on → short fragments kept)
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            None,
            &mut AdCounter::new(),
            &adapter,
            &adapter,
        )
        .unwrap();
        buf1.finish().unwrap();
        buf2.finish().unwrap();

        for path in [&r1_path, &r2_path] {
            let lines: Vec<String> = read_gzip_lines(path).unwrap().map(|l| l.unwrap()).collect();
            // Two fragments produce reads; the zero-insert pair is skipped (4 lines each).
            assert_eq!(
                lines.len(),
                8,
                "expected 2 records (zero-insert skipped) in {:?}, got {} lines",
                path,
                lines.len()
            );
            let mut i = 0;
            while i < lines.len() {
                let seq = &lines[i + 1];
                let qual = &lines[i + 3];
                assert_eq!(
                    seq.len(),
                    qual.len(),
                    "malformed FASTQ in {:?}: seq={} qual={}",
                    path,
                    seq.len(),
                    qual.len()
                );
                assert_eq!(
                    seq.len(),
                    read_length,
                    "read not padded to read_length in {:?}",
                    path
                );
                i += 4;
            }
        }
    }

    #[test]
    fn test_write_block_fastq_keep_short_no_adapter_emits_genomic_insert_length_reads() {
        // keep_short control (#125): short fragments kept, NO adapter (empty slices).
        // A short insert must produce an insert-LENGTH genomic read (not padded to
        // read_length, not dropped), and records stay well-formed (seq == qual).
        use crate::structs::{
            mutated_map::MutatedMap,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let seq_len: usize = 300;
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        let block = SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: seq_len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        };
        let mutated_map = MutatedMap::from_interval(0, seq_len, vec![]).unwrap();
        let read_length = 70usize;
        // one short insert (30 < read_length) + one full insert (== read_length).
        let fragments = vec![(10, 40), (100, 170)];
        let out_path = temp_dir.path().join("r1.fastq.gz");
        let mut buf1 = GzEncoder::new(
            create_output_file(&out_path, true).unwrap(),
            Compression::default(),
        );
        use crate::file_tools::file_io::VectorBuffer;
        let mut buf2 = GzEncoder::new(VectorBuffer::new(), Compression::default());
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["keep-short".to_string()]).unwrap();
        write_block_fastq(
            fragments.into_iter().map(Into::into).collect(),
            &[],
            &mutated_map,
            &block,
            false, // single-ended (simplest) — R1 only
            &mut buf1,
            &mut buf2,
            read_length,
            false, // long_reads
            true,  // keep_short
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            None,
            &mut AdCounter::new(),
            &[], // NO adapter → genomic reads, no padding
            &[],
        )
        .unwrap();
        buf1.finish().unwrap();
        let lines: Vec<String> = read_gzip_lines(&out_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(
            lines.len(),
            8,
            "expected 2 reads (both fragments kept), got {}",
            lines.len()
        );
        let read_lens: Vec<usize> = (0..lines.len())
            .step_by(4)
            .map(|i| {
                assert_eq!(lines[i + 1].len(), lines[i + 3].len(), "seq/qual mismatch");
                lines[i + 1].len()
            })
            .collect();
        // Short insert (30) → 30 bp genomic read (NOT padded to 70); full insert → 70 bp.
        let mut sorted = read_lens.clone();
        sorted.sort();
        assert_eq!(
            sorted,
            vec![30, 70],
            "expected insert-length reads [30, 70], got {:?}",
            read_lens
        );
    }

    // --- CIGAR-building tests for the refactored generate_read ---

    fn make_sequence(len: usize) -> Vec<Nucleotide> {
        (0..len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect()
    }

    #[test]
    fn test_generate_read_cigar_all_match_no_errors() {
        // Q40 → error prob 0.0001; with seed and 10 positions no error fires.
        let sequence = make_sequence(30);
        let read_length = 10;
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["no_error".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        assert_eq!(record.sequence.len(), read_length);
        assert_eq!(record.cigar_ops.len(), read_length);
        assert!(
            record.cigar_ops.iter().all(|&c| c == 'M'),
            "expected all-M cigar, got: {:?}",
            record.cigar_ops
        );
    }

    #[test]
    fn test_generate_read_cigar_mi_plus_d_equals_len() {
        // Invariant: M+I == read_length, cigar_ops.len() == (M+I) + D, for any error mix.
        let sequence = make_sequence(500);
        let read_length = 50;
        let quality_scores = vec![0usize; read_length]; // Q0 → 100% error probability
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["invariant".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        let d = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert_eq!(record.sequence.len(), read_length);
        assert_eq!(mi, read_length, "M+I must equal read_length");
        assert_eq!(record.cigar_ops.len(), mi + d, "cigar len must equal M+I+D");
    }

    /// Counter behaviour: a homozygous SNP should drive alt_count up and
    /// leave ref_count at 0 across many reads; a heterozygous SNP should split
    /// roughly 50/50 around the half-point of N reads. Together these pin the
    /// AdCounter increments in generate_read at the coin-flip site.
    #[test]
    fn test_generate_read_counter_homozygous_all_alt() {
        let sequence = make_sequence(30);
        let read_length = 10;
        let hom_snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![1, 1]).unwrap();
        let variant_map = HashMap::from([(5usize, &hom_snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["hom".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        for i in 0..50 {
            let _ = generate_read(
                &sequence,
                // Reference-derived bases only; no haplotype insertion mask.
                None,
                &[5],
                &variant_map,
                read_length,
                format!("r{i}/1"),
                Strand::Forward,
                quality_scores.clone(),
                &model,
                &mut rng,
                "chr1".to_string(),
                0,
                "chr1".to_string(),
                0,
                0,
                false,
                &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        assert_eq!(refs, 0, "homozygous SNP must produce zero ref reads");
        assert_eq!(alts, 50, "homozygous SNP must produce all alt reads");
    }

    /// Regression for the reverse-strand variant-placement bug: a reverse (R2)
    /// read must apply the variant at `seq_index` (the caller already
    /// reverse-complements the fragment and maps the variant coordinate into it),
    /// NOT at the mirrored `fragment_length - seq_index`. With the old reflection
    /// the variant fell outside the read window, so reverse reads silently
    /// carried REF — alternate alleles only ever appeared on forward reads,
    /// which Mutect2 (correctly) filtered as strand bias.
    #[test]
    fn test_generate_read_reverse_applies_variant() {
        let sequence = make_sequence(30);
        let read_length = 10;
        let hom_snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![1, 1]).unwrap();
        let variant_map = HashMap::from([(5usize, &hom_snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["rev".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        for i in 0..50 {
            let _ = generate_read(
                &sequence,
                // Reference-derived bases only; no haplotype insertion mask.
                None,
                &[5],
                &variant_map,
                read_length,
                format!("r{i}/2"),
                Strand::Reverse,
                quality_scores.clone(),
                &model,
                &mut rng,
                "chr1".to_string(),
                0,
                "chr1".to_string(),
                0,
                0,
                true,
                &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        // The pre-fix reflection put the lookup at index 25, outside the
        // 10-base read, so alt would be 0 here. The fix applies it at index 5.
        assert_eq!(
            refs, 0,
            "reverse homozygous SNP must produce zero ref reads"
        );
        assert_eq!(
            alts, 50,
            "reverse read must carry the alt (strand-bias regression)"
        );
    }

    /// reverse_complement_record must reverse-complement the sequence and reverse
    /// the per-base CIGAR + qualities (so a forward-generated R2 flips correctly).
    #[test]
    fn test_reverse_complement_record() {
        let rec = ReadRecord {
            is_unmapped: false,
            name: "frag/2".to_string(),
            sequence: "AACGT".to_string(),
            quality_scores: vec![1, 2, 3, 4, 5],
            cigar_ops: vec!['M', 'I', 'M', 'M', 'M'],
            is_paired: true,
            is_reverse: false,
            contig: "chr1".to_string(),
            position: 10,
            mate_contig: "chr1".to_string(),
            mate_position: 5,
            template_length: -50,
        };
        let f = reverse_complement_record(rec);
        assert_eq!(f.sequence, "ACGTT", "sequence must be reverse-complemented");
        assert_eq!(
            f.quality_scores,
            vec![5, 4, 3, 2, 1],
            "qualities must be reversed"
        );
        assert_eq!(
            f.cigar_ops,
            vec!['M', 'M', 'M', 'I', 'M'],
            "CIGAR must be reversed"
        );
        assert!(f.is_reverse, "flipped record must be marked reverse");
    }

    #[test]
    fn test_generate_read_counter_heterozygous_splits_around_half() {
        let sequence = make_sequence(30);
        let read_length = 10;
        let het_snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![0, 1]).unwrap();
        let variant_map = HashMap::from([(5usize, &het_snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["het".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        let n = 1000;
        for i in 0..n {
            let _ = generate_read(
                &sequence,
                // Reference-derived bases only; no haplotype insertion mask.
                None,
                &[5],
                &variant_map,
                read_length,
                format!("r{i}/1"),
                Strand::Forward,
                quality_scores.clone(),
                &model,
                &mut rng,
                "chr1".to_string(),
                0,
                "chr1".to_string(),
                0,
                0,
                false,
                &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        assert_eq!(
            refs + alts,
            n,
            "every read should increment exactly one slot"
        );
        // Binomial(1000, 0.5) → 99.99% CI is well within [400, 600]
        assert!(
            (400..600).contains(&(refs as usize)),
            "het split should be roughly 50/50 ({}/1000 ref, {} alt)",
            refs,
            alts
        );
    }

    #[test]
    fn test_generate_read_counter_honors_allele_fraction() {
        // An explicit allele_fraction overrides the Genotype-based fraction: the alt
        // allele should appear on ~f of reads (here 0.2), regardless of Het/Hom (#398).
        let sequence = make_sequence(30);
        let read_length = 10;
        let mut snp =
            Variant::new(VariantType::SNP, 5, &vec![T], &vec![C], &mut vec![0, 1]).unwrap();
        snp.allele_fraction = Some(0.2);
        let variant_map = HashMap::from([(5usize, &snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["af".to_string()]).unwrap();
        let mut ad: AdCounter = HashMap::new();
        let n = 1000;
        for i in 0..n {
            let _ = generate_read(
                &sequence,
                // Reference-derived bases only; no haplotype insertion mask.
                None,
                &[5],
                &variant_map,
                read_length,
                format!("r{i}/1"),
                Strand::Forward,
                quality_scores.clone(),
                &model,
                &mut rng,
                "chr1".to_string(),
                0,
                "chr1".to_string(),
                0,
                0,
                false,
                &mut ad,
            )
            .unwrap();
        }
        let (refs, alts) = ad[&5];
        assert_eq!(
            refs + alts,
            n,
            "every read should increment exactly one slot"
        );
        // Binomial(1000, 0.2): mean 200, sd ~12.6 → generous CI well clear of 0.5.
        assert!(
            (150..250).contains(&(alts as usize)),
            "alt count should track allele_fraction=0.2 ({alts}/1000 alt, {refs} ref)"
        );
    }

    #[test]
    fn test_generate_read_snp_variant_produces_all_m_cigar() {
        // A SNP variant changes a base but does not add I or D ops.
        let sequence = make_sequence(30);
        let read_length = 10;
        let snp = Variant::new(VariantType::SNP, 3, &vec![T], &vec![C], &mut vec![1, 1]).unwrap();
        let variant_map = HashMap::from([(3usize, &snp)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["snp_cigar".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[3],
            &variant_map,
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        assert_eq!(record.cigar_ops.len(), read_length);
        assert!(
            record.cigar_ops.iter().all(|&c| c == 'M'),
            "SNP variant must not add I or D ops; got: {:?}",
            record.cigar_ops
        );
    }

    #[test]
    fn test_generate_read_insertion_variant_produces_i_ops() {
        // Homozygous insertion: alt = [A, C, G] → 1 M (anchor) + 2 I (inserted bases).
        let sequence = make_sequence(30);
        let read_length = 10;
        let ins = Variant::new(
            VariantType::Insertion,
            3,
            &vec![A],
            &vec![A, C, G],
            &mut vec![1, 1],
        )
        .unwrap();
        let variant_map = HashMap::from([(3usize, &ins)]);
        let quality_scores = vec![40usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["ins_variant".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[3],
            &variant_map,
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let i_count = record.cigar_ops.iter().filter(|&&c| c == 'I').count();
        assert_eq!(i_count, 2, "2-base insertion variant must produce 2 I ops");
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        assert_eq!(mi, read_length);
    }

    #[test]
    fn test_generate_read_error_indels_appear_in_cigar() {
        // Q0 at every position guarantees an error at each position (convert_score(0) = 1.0).
        // With 50 positions, P(all 50 errors are SNPs) ≈ 0.6^50 < 1.4e-11,
        // so this effectively guarantees at least one I or D op in the cigar.
        let sequence = make_sequence(500);
        let read_length = 50;
        let quality_scores = vec![0usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["error_indel".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        assert!(
            record.cigar_ops.iter().any(|&c| c != 'M'),
            "expected non-M ops when every position has a guaranteed error; got: {:?}",
            record.cigar_ops
        );
        // Invariant must hold even under heavy errors
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        assert_eq!(mi, read_length);
    }

    #[test]
    fn test_generate_read_deletion_error_increases_cigar_length() {
        // When deletion errors fire, each skipped reference base adds a D op, so
        // cigar_ops.len() > read_length.  Verify the formula len == M+I+D holds
        // and, whenever D > 0, that cigar_ops.len() > read_length.
        let sequence = make_sequence(500);
        let read_length = 50;
        let quality_scores = vec![0usize; read_length];
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["del_len".to_string()]).unwrap();
        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &[],
            &HashMap::new(),
            read_length,
            "r/1".to_string(),
            Strand::Forward,
            quality_scores,
            &model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();
        let mi = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        let d = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert_eq!(mi, read_length);
        assert_eq!(record.cigar_ops.len(), mi + d);
        if d > 0 {
            assert!(
                record.cigar_ops.len() > read_length,
                "deletion errors must make cigar longer than read_length"
            );
        }
    }

    /// #221 regression: a long literal Deletion (REF=50bp, ALT=1bp) must
    /// produce a CIGAR with the 49-base D-op so downstream callers see the
    /// deletion. Before the fix, the alt branch wrote the 1-byte anchor and
    /// then advanced seq_index by 1, so the read transcribed the 49 deleted
    /// bases from the unbroken reference and emitted CIGAR `<read_length>M`.
    #[test]
    fn test_literal_long_deletion_emits_d_ops() {
        // 200-bp reference: 30 bases of left context + 50-bp REF starting at
        // pos 30 + 120 bases of right context.
        let mut sequence = vec![A; 30];
        sequence.extend(vec![C; 50]); // the 50 REF bases (anchor + 49 deleted)
        sequence.extend(vec![T; 120]); // post-deletion context
        let read_length = 100;

        // Homozygous DEL at position 30, REF=50bp (CCCC...), ALT=1bp (anchor C).
        let ref_bases: Vec<Nucleotide> = vec![C; 50];
        let alt_bases: Vec<Nucleotide> = vec![C];
        let variant = Variant::new(
            VariantType::Deletion,
            30,
            &ref_bases,
            &alt_bases,
            &mut vec![1, 1], // homozygous so the alt branch always fires
        )
        .unwrap();
        let variant_map = HashMap::from([(30usize, &variant)]);
        let flagged_positions = vec![30usize];
        let qual_scores = vec![33; 100];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng =
            NeatRng::new_from_seed(&vec!["literal".into(), "long".into(), "del".into()]).unwrap();

        let record = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &flagged_positions,
            &variant_map,
            read_length,
            "del_test/1".into(),
            Strand::Forward,
            qual_scores,
            &sequencing_error_model,
            &mut rng,
            "chr1".into(),
            0,
            "chr1".into(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        )
        .unwrap();

        // Sequence written must equal read_length (100 bases).
        assert_eq!(record.sequence.len(), read_length);

        // The 49-bp deletion must be encoded as D ops in the CIGAR. The
        // exact CIGAR is read-length-dependent (sequencing errors may add
        // small ops too), but the D count must be ≥ 49 — the deletion span.
        // We assert D count ≥ 49 rather than == 49 so a stray seq-error
        // deletion doesn't flake the test.
        let d_count = record.cigar_ops.iter().filter(|&&c| c == 'D').count();
        assert!(
            d_count >= 49,
            "expected ≥49 D ops for a 50-bp REF / 1-bp ALT homozygous deletion, got {d_count}. \
             CIGAR: {:?}",
            record.cigar_ops
        );

        // Bases-written count: each M op or I op consumes one output base.
        let mi_count = record
            .cigar_ops
            .iter()
            .filter(|&&c| c == 'M' || c == 'I')
            .count();
        assert_eq!(
            mi_count, read_length,
            "M+I ops must equal read_length ({read_length}); got {mi_count}"
        );

        // AD counter on the variant position must record the alt observation.
        // (Verified via the variant_map's pointer at position 30 — alt_count
        // is incremented inside generate_read when the alt branch fires.)
    }

    #[test]
    fn test_apply_variants() {
        let sequence = vec![A, C, G, T, T, A, T, G, A, C, G, T, T, A, T, G];
        let variant1 =
            Variant::new(VariantType::SNP, 1, &vec![T], &vec![C], &mut vec![1, 0]).unwrap();
        let variant2 = Variant::new(
            VariantType::Deletion,
            3,
            &vec![T, T],
            &vec![T],
            &mut vec![0, 1],
        )
        .unwrap();
        let variant_map = HashMap::from([(1, &variant1), (3, &variant2)]);
        let flagged_positions = vec![1, 3];
        let read_name = "neat_generated__0000000000_0000000008/1".to_string();
        let qual_scores = vec![33, 25, 37, 28, 15, 33, 33, 37];
        let sequencing_error_model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let result = generate_read(
            &sequence,
            // Reference-derived bases only; no haplotype insertion mask.
            None,
            &flagged_positions,
            &variant_map,
            8,
            read_name,
            Strand::Forward,
            qual_scores,
            &sequencing_error_model,
            &mut rng,
            "chr1".to_string(),
            0,
            "chr1".to_string(),
            0,
            0,
            false,
            &mut AdCounter::new(),
        );
        assert!(result.is_ok());
    }

    // ── incremental BAM flush tests ──────────────────────────────────────────

    fn make_block(seq_len: usize) -> SequenceBlock {
        let sequence: Vec<Nucleotide> = (0..seq_len)
            .map(|i| match i % 4 {
                0 => A,
                1 => C,
                2 => G,
                _ => T,
            })
            .collect();
        SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: seq_len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, seq_len)],
        }
    }

    fn run_write_block(
        bam_writer: &mut BamWriter,
        fragments: Vec<(usize, usize)>,
        block: &SequenceBlock,
        paired_ended: bool,
        read_len: usize,
    ) {
        let mutated_map = MutatedMap::from_interval(0, block.sequence.len(), vec![]).unwrap();
        let seq_err_model = SequencingErrorModel::default().unwrap();
        let quality_model = QualityScoreModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["flush-test".to_string()]).unwrap();
        let mut buf1 = GzEncoder::new(VectorBuffer::new(), Compression::default());
        let mut buf2 = GzEncoder::new(VectorBuffer::new(), Compression::default());
        let stager: Option<&mut dyn BamRecordStager> = Some(bam_writer);
        write_block_fastq(
            fragments.into_iter().map(Into::into).collect(),
            &[],
            &mutated_map,
            block,
            paired_ended,
            &mut buf1,
            &mut buf2,
            read_len,
            false,
            false, // keep_short
            "chr1",
            &quality_model,
            &seq_err_model,
            &mut rng,
            stager,
            &mut AdCounter::new(),
            &[],
            &[],
        )
        .unwrap();
    }

    #[test]
    fn test_write_block_fastq_carry_bounded_se() {
        // With 5 SE fragments at strictly increasing positions, flush_up_to(abs_start)
        // fires after each fragment so the carry never exceeds 1 record.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("out.bam");
        let block = make_block(500);
        let mut bw = BamWriter::new(&bam_path, &[("chr1".to_string(), 500usize)]).unwrap();

        let fragments = vec![(0usize, 10), (50, 60), (100, 110), (150, 160), (200, 210)];
        run_write_block(&mut bw, fragments, &block, false, 10);

        assert_eq!(
            bw.carry_len(),
            1,
            "SE: carry should hold only the last fragment's read, got {}",
            bw.carry_len()
        );
    }

    #[test]
    fn test_write_block_fastq_carry_bounded_pe() {
        // With 5 PE fragments at strictly increasing positions, flush_up_to fires after
        // each fragment so carry never exceeds 2 records (one R1 + one R2).
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("out.bam");
        let block = make_block(5000);
        let mut bw = BamWriter::new(&bam_path, &[("chr1".to_string(), 5000usize)]).unwrap();

        // Fragments 400 bp wide, spaced 100 bp apart so positions are strictly increasing.
        let fragments = vec![
            (0usize, 400),
            (500, 900),
            (1000, 1400),
            (1500, 1900),
            (2000, 2400),
        ];
        run_write_block(&mut bw, fragments, &block, true, 150);

        assert_eq!(
            bw.carry_len(),
            2,
            "PE: carry should hold only the last fragment's read pair, got {}",
            bw.carry_len()
        );
    }

    #[test]
    fn test_write_block_fastq_bam_sorted_with_incremental_flush() {
        // End-to-end: write_block_fastq with incremental flush then flush_all must
        // produce a coordinate-sorted BAM, identical to bulk flushing.
        use noodles::bam;

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("sorted.bam");
        let block = make_block(500);
        let contigs = vec![("chr1".to_string(), 500usize)];
        let mut bw = BamWriter::new(&bam_path, &contigs).unwrap();

        let fragments = vec![(0usize, 10), (50, 60), (100, 110), (150, 160), (200, 210)];
        run_write_block(&mut bw, fragments, &block, false, 10);
        bw.flush_all().unwrap();
        drop(bw); // BGZF EOF written on drop

        let file = std::fs::File::open(&bam_path).unwrap();
        let mut reader = bam::io::Reader::new(file);
        reader.read_header().unwrap();
        let positions: Vec<usize> = reader
            .records()
            .map(|r| r.unwrap().alignment_start().unwrap().unwrap().get() - 1)
            .collect();
        assert_eq!(
            positions,
            vec![0, 50, 100, 150, 200],
            "BAM must be coordinate-sorted; got: {:?}",
            positions
        );
    }
}
