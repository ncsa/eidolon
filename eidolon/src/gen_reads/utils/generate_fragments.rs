// This is the core functionality of NEAT. Generate reads turns a mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to give the contig the proper
// read coverage. Generate reads uses this to create a list of coordinates to take slices from
// the mutated fasta file. These will either be read-length fragments or fragment model length
// fragments.
use crate::{
    eidolon_core::models::fragment_length::FragmentLengthModel,
    gen_reads::errors::GenerateReadsError,
};
use eidolon_core::models::gc_bias_model::GcBiasModel;
use eidolon_core::rng::NeatRng;
use eidolon_core::structs::nucleotides::Nucleotide;
use eidolon_core::structs::sequence_block::SequenceBlock;
use log::*;

pub fn generate_fragments(
    sequence_length: usize,
    read_length: usize,
    max_del_len: usize,
    start: usize,
    coverage: usize,
    paired_ended: bool,
    long_reads: bool,
    // when true, KEEP short fragments (insert < read_len) instead of discarding them:
    // adapter readthrough pads them to read_len, or the adapter-free short-insert
    // control emits them as insert-length genomic reads, in write_block_fastq (#125).
    keep_short: bool,
    fragment_model: &FragmentLengthModel,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    // Takes:
    // sequence_length: The length of the sequence to generate reads for.
    // read_length: the length of the reads for this run
    // coverage: the average depth of coverage for this run
    // paired_ended: when false, read_length is used as the fragment length so that
    //   coverage is computed in read-length steps rather than insert-size steps.
    // long_reads: when true, fragments shorter than read_length are accepted and produce
    //   truncated reads (appropriate for PacBio HiFi / ONT paired-end simulation).
    // rng: the random number generator for the run
    // Returns:
    // A vector of (start, end) pairs representing fragment coordinates.
    let mut fragment_pool: Vec<usize> = Vec::new();
    // For long-read paired-end mode the region only needs deletion padding; for all other
    // modes the full read_length must fit.
    let min_region = if long_reads && paired_ended {
        max_del_len * 2
    } else {
        read_length + max_del_len * 2
    };
    if sequence_length <= min_region {
        debug!("Sequence length was too short, maybe because of a small bed region.");
        return Ok(Vec::new());
    }

    // Each paired-end fragment produces two reads (R1 + R2), so halve the fragment count
    // so that mean per-base read depth matches the requested coverage value.
    // Use ceiling division (matching NEAT Python's ceil()) so fractional fragments are
    // rounded up rather than truncated — truncation compounds across many small regions.
    let denom = if paired_ended {
        (2 * read_length).max(1)
    } else {
        read_length
    };
    let num_frags = sequence_length.saturating_mul(coverage).div_ceil(denom);
    if paired_ended {
        // For paired-end reads the physical fragment length (insert size) determines spacing.
        // Keep sampling until num_frags valid fragments are collected so that the pool has
        // full size diversity even when the model's distribution overlaps with read_length.
        let max_attempts = num_frags.saturating_mul(10);
        let mut attempts = 0;
        fragment_pool.reserve(num_frags);
        while fragment_pool.len() < num_frags && attempts < max_attempts {
            let frag = fragment_model.generate_fragment(rng.random()?)?;
            // In long-read mode accept all fragments; in short-read paired-end mode require
            // frag >= read_length + 10 so that R1 and R2 don't fully overlap (matching
            // NEAT Python's min_frag = read_len + 10 for paired-end).
            if long_reads || keep_short || frag >= read_length + 10 {
                fragment_pool.push(frag);
            }
            attempts += 1;
        }
        if fragment_pool.len() < num_frags {
            warn!(
                "Fragment model produced only {}/{} valid fragments after {} attempts; \
                 insert-size diversity may be reduced. Check that fragment_mean >> read_length.",
                fragment_pool.len(),
                num_frags,
                attempts
            );
        }
    } else {
        // For single-end reads the fragment is exactly one read long; the fragment model
        // is irrelevant for spacing and would cause under-coverage if its mean >> read_length.
        fragment_pool = vec![read_length; 1];
    }
    if fragment_pool.is_empty() {
        debug!("Fragment pool is empty.");
        return Ok(Vec::new());
    }

    // Generate a vector of read positions
    debug!("Generating read coordinates.");
    let fragments: Vec<(usize, usize)> = cover_dataset(
        sequence_length,
        read_length,
        num_frags,
        start,
        fragment_pool,
        rng,
    )?;

    if fragments.is_empty() {
        debug!("No fragments generated!");
        Ok(Vec::new())
    } else {
        Ok(fragments)
    }
}

/// Builds an unnormalized prefix-sum CDF over per-position GC bias weights for a region.
///
/// Each entry `prefix_sum[i+1] = prefix_sum[i] + weight(window starting at region_start + i)`.
/// The model's `window_size` is clamped to `region_len` so short regions are handled gracefully.
/// Returns `None` when the region is empty; the caller maps that to zero fragments.
///
/// Uses a sliding window to maintain a running GC and N count, avoiding per-position
/// Vec allocations that would dominate cost on large chromosomes.
fn build_gc_weight_prefix_sum(
    sequence_block: &SequenceBlock,
    region_start: usize,
    region_end: usize,
    gc_bias_model: &GcBiasModel,
) -> Result<Option<Vec<f64>>, GenerateReadsError> {
    if region_end <= region_start {
        return Ok(None);
    }
    let region_len = region_end - region_start;
    let window = gc_bias_model.window_size().min(region_len);
    let n_positions = region_len - window + 1;

    let seq = &sequence_block.sequence[region_start..region_end];

    // Initialise counts for the first window.
    let mut gc_count: usize = seq[..window]
        .iter()
        .filter(|&&n| n == Nucleotide::G || n == Nucleotide::C)
        .count();
    let mut n_count: usize = seq[..window]
        .iter()
        .filter(|&&n| n == Nucleotide::N)
        .count();

    let weight_at = |gc: usize, n: usize| -> f64 {
        let called = window - n;
        if called == 0 {
            1.0 // all-N window: neutral, never filtered
        } else {
            gc_bias_model.weight_for_gc_fraction(gc as f64 / called as f64)
        }
    };

    let mut prefix_sum = vec![0.0f64; n_positions + 1];
    prefix_sum[1] = weight_at(gc_count, n_count);

    // Slide the window one base at a time, updating counts incrementally.
    for i in 1..n_positions {
        let outgoing = seq[i - 1];
        let incoming = seq[i + window - 1];

        match outgoing {
            Nucleotide::G | Nucleotide::C => gc_count -= 1,
            Nucleotide::N => n_count -= 1,
            _ => {}
        }
        match incoming {
            Nucleotide::G | Nucleotide::C => gc_count += 1,
            Nucleotide::N => n_count += 1,
            _ => {}
        }

        prefix_sum[i + 1] = prefix_sum[i] + weight_at(gc_count, n_count);
    }

    Ok(Some(prefix_sum))
}

/// Samples a 0-based position offset from the prefix-sum CDF using a single uniform draw.
fn sample_from_prefix_sum(
    prefix_sum: &[f64],
    total_weight: f64,
    rng: &mut NeatRng,
) -> Result<usize, GenerateReadsError> {
    let v = rng.random()? * total_weight;
    // Find the last index k where prefix_sum[k] <= v, i.e. the interval [k, k+1) containing v.
    let k = prefix_sum.partition_point(|&ps| ps <= v).saturating_sub(1);
    Ok(k)
}

/// Generates fragments for a region by sampling start positions directly from a per-position
/// GC bias weight CDF, with no rejection step.
///
/// Start positions are drawn with probability proportional to the model weight of the window
/// at each position, so high-weight positions are naturally over-sampled and low-weight
/// positions are under-sampled, exactly matching the sequencer's bias.
///
/// When `normalize` is true the fragment count is inflated by `1/mean_weight` (capped at
/// `100×`) so that the expected total coverage still approaches the target even in low-weight
/// regions. When `normalize` is false the fragment count equals what a uniform sequencer
/// would generate, and coverage will be proportionally lower in low-weight regions.
pub fn generate_weighted_fragments(
    sequence_block: &SequenceBlock,
    region_start: usize,
    region_end: usize,
    read_length: usize,
    max_del_len: usize,
    coverage: usize,
    gc_bias_model: &GcBiasModel,
    fragment_model: &FragmentLengthModel,
    normalize: bool,
    paired_ended: bool,
    long_reads: bool,
    // keep short fragments (readthrough or the adapter-free control); see generate_fragments.
    keep_short: bool,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    const MAX_COVERAGE_MULTIPLIER: usize = 100;

    let region_len = region_end.saturating_sub(region_start);
    let min_region = if long_reads {
        max_del_len * 2
    } else {
        read_length + max_del_len * 2
    };
    if region_len <= min_region {
        debug!("Region too short for weighted fragment generation.");
        return Ok(Vec::new());
    }

    let prefix_sum = match build_gc_weight_prefix_sum(
        sequence_block,
        region_start,
        region_end,
        gc_bias_model,
    )? {
        Some(ps) => ps,
        None => return Ok(Vec::new()),
    };

    let total_weight = *prefix_sum.last().unwrap();
    if total_weight == 0.0 {
        debug!("All positions in region have zero GC bias weight; no fragments generated.");
        return Ok(Vec::new());
    }

    let n_valid_positions = prefix_sum.len() - 1;
    let mean_weight = total_weight / n_valid_positions as f64;

    // Same paired-end correction as generate_fragments: each fragment yields two reads,
    // so halve the count so that mean per-base read depth matches requested coverage.
    // Use ceiling division so fractional fragments are rounded up, not truncated.
    let denom = if paired_ended {
        (2 * read_length).max(1)
    } else {
        read_length
    };
    let num_frags_base = region_len.saturating_mul(coverage).div_ceil(denom);
    let num_frags = if normalize {
        let inflated = (num_frags_base as f64 / mean_weight).round() as usize;
        let cap = num_frags_base.saturating_mul(MAX_COVERAGE_MULTIPLIER);
        if inflated > cap {
            warn!(
                "GC bias coverage inflation capped at {}x (mean weight {:.4}); \
                 target coverage may not be reached in this region",
                MAX_COVERAGE_MULTIPLIER, mean_weight
            );
        }
        inflated.min(cap)
    } else {
        num_frags_base
    };

    if num_frags == 0 {
        return Ok(Vec::new());
    }

    // Keep sampling until num_frags fragments are placed so that fragment model
    // rejection and edge clamping don't silently reduce coverage.
    let max_attempts = num_frags.saturating_mul(10);
    let mut attempts = 0;
    let mut fragments = Vec::with_capacity(num_frags);
    while fragments.len() < num_frags && attempts < max_attempts {
        let offset = sample_from_prefix_sum(&prefix_sum, total_weight, rng)?;
        let start = region_start + offset;
        let frag_len = fragment_model.generate_fragment(rng.random()?)?;
        attempts += 1;
        if !long_reads && !keep_short && frag_len < read_length {
            continue;
        }
        let end = (start + frag_len).min(region_end);
        if end == start {
            continue;
        }
        if !long_reads && !keep_short && end - start < read_length {
            continue;
        }
        fragments.push((start, end));
    }
    if fragments.len() < num_frags {
        warn!(
            "GC-weighted fragment generation placed {}/{} fragments after {} attempts; \
             coverage may be below target in this region.",
            fragments.len(),
            num_frags,
            attempts
        );
    }

    fragments.sort_by_key(|&(s, _)| s);
    Ok(fragments)
}

fn cover_dataset(
    span_length: usize,
    read_length: usize,
    target_count: usize,
    read_start: usize,
    mut fragment_pool: Vec<usize>,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    // Places exactly target_count fragments so that mean depth equals the requested
    // coverage. Stopping on fragment count rather than sweep count means the wildcard
    // gap can be proportional to read_length without reducing mean coverage.
    // Alternates sweep direction (3' to 5' then 5' to 3') for better realism.

    rng.shuffle_in_place(&mut fragment_pool)?;
    let pool_size = fragment_pool.len();
    let mut pool_idx = 0usize;
    let mut fragment_set: Vec<(usize, usize)> = Vec::with_capacity(target_count);
    let mut pos = (rng.rand_int()? as usize) % (read_length / 4).max(1);
    let mut forward = true;
    let mut non_placing_streak = 0usize;

    while fragment_set.len() < target_count {
        let fragment_length = fragment_pool[pool_idx];
        pool_idx += 1;
        if pool_idx == pool_size {
            pool_idx = 0;
        }

        let (start, end) = if forward {
            (pos, pos + fragment_length)
        } else {
            (pos.saturating_sub(fragment_length), pos)
        };

        if (forward && end > span_length) || (!forward && start == 0 && pos < fragment_length) {
            // Sweep complete, reverse direction and pick a random starting offset.
            forward = !forward;
            pos = if forward {
                (rng.rand_u32()? as usize) % (read_length / 4).max(1)
            } else {
                span_length - (rng.rand_u32()? as usize) % (read_length / 4).max(1)
            };

            // Only count fragments that can NEVER fit (length > span), not fragments
            // that merely overshot from a high start position and will fit after restart.
            if fragment_length > span_length {
                non_placing_streak += 1;
                // If we've seen pool_size consecutive fragments that are each individually
                // larger than the span, every fragment in the pool is too large — break.
                if non_placing_streak >= pool_size {
                    debug!(
                        "All fragments exceed span of {} bp; stopping early with {} fragments placed.",
                        span_length,
                        fragment_set.len()
                    );
                    break;
                }
            } else {
                non_placing_streak = 0;
            }
            continue;
        }
        non_placing_streak = 0;

        fragment_set.push((start + read_start, end + read_start));
        // Halved from read_length/4 to read_length/8 to tighten spacing between consecutive
        // fragment placements, reducing systematic coverage gaps between sweep passes.
        let wildcard = (rng.rand_u32()? % (read_length as u32 / 8).max(1)) as usize;

        if forward {
            pos = end + wildcard;
            if pos >= span_length {
                forward = false;
                pos = span_length - (rng.rand_u32()? as usize) % (read_length / 4).max(1);
            }
        } else {
            pos = start.saturating_sub(wildcard);
            if pos == 0 {
                forward = true;
                pos = (rng.rand_u32()? as usize) % (read_length / 4).max(1);
            }
        }
    }

    fragment_set.sort_by_key(|&(s, _)| s);
    Ok(fragment_set)
}

/// Prototype helper for #516: produce read windows whose coordinate system is
/// the inserted haplotype sequence rather than the reference contig.
///
/// This is deliberately not wired into production read generation yet. The
/// current writer still needs a mapping from these haplotype coordinates back
/// to reference/BAM coordinates. Keeping that mapping separate lets us test
/// the geometric requirement first: long insertions must have fragments whose
/// R1/R2 windows begin in the insertion interior, not only at its anchor.
#[cfg(test)]
pub(crate) fn prototype_insertion_read_windows(
    insertion_len: usize,
    read_length: usize,
    target_count: usize,
    fragment_pool: Vec<usize>,
    rng: &mut NeatRng,
) -> Result<Vec<(usize, usize)>, GenerateReadsError> {
    if insertion_len == 0 || read_length == 0 || target_count == 0 {
        return Ok(Vec::new());
    }
    let fragments = cover_dataset(
        insertion_len,
        read_length,
        target_count,
        0,
        fragment_pool,
        rng,
    )?;
    let mut windows = Vec::with_capacity(fragments.len() * 2);
    for (start, end) in fragments {
        let r1_end = (start + read_length).min(end);
        if r1_end > start {
            windows.push((start, r1_end));
        }
        let r2_start = end.saturating_sub(read_length).max(start);
        if end > r2_start {
            windows.push((r2_start, end));
        }
    }
    windows.sort_unstable();
    windows.dedup();
    Ok(windows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::models::gc_bias_model::GcBiasModel;
    use eidolon_core::rng::NeatRng;
    use eidolon_core::structs::nucleotides::Nucleotide::{A, C, G, T};
    use eidolon_core::structs::sequence_block::{RegionType, SequenceBlock, SequenceMap};

    fn make_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec!["gc-bias-test".to_string()]).unwrap()
    }

    fn make_sequence_block(
        sequence: Vec<eidolon_core::structs::nucleotides::Nucleotide>,
    ) -> SequenceBlock {
        let len = sequence.len();
        SequenceBlock {
            contig: "chr1".to_string(),
            ref_start: 0,
            ref_end: len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, len)],
        }
    }

    fn weights_with_value(value: f64) -> Vec<f64> {
        vec![value; 101]
    }

    // ── generate_fragments (uniform path) ────────────────────────────────────

    #[test]
    fn test_cover_dataset() {
        let span_length = 100;
        let read_length = 10;
        let fragment_pool = vec![10];
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let target_count = span_length * coverage / read_length;
        let cover = cover_dataset(
            span_length,
            read_length,
            target_count,
            0,
            fragment_pool,
            &mut rng,
        )
        .unwrap();
        assert_eq!(cover.len(), target_count);
        // All fragments must fit within the span.
        assert!(cover.iter().all(|&(_s, e)| e <= span_length));
    }

    #[test]
    fn test_gap_function() {
        let span_length = 100_000;
        let read_length = 100;
        let fragment_pool = vec![300];
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let target_count = span_length / read_length * coverage;
        let cover = cover_dataset(
            span_length,
            read_length,
            target_count,
            0,
            fragment_pool,
            &mut rng,
        )
        .unwrap();
        assert_eq!(cover.len(), target_count);
        // All fragments must fit within the span and have the pool fragment length.
        assert!(cover.iter().all(|&(s, e)| e <= span_length && e - s == 300));
    }

    #[test]
    fn test_generate_reads_single() {
        let read_length = 10;
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let reads = generate_fragments(
            2000,
            read_length,
            0,
            0,
            coverage,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &fragment_model,
            &mut rng,
        )
        .unwrap();
        // With paired_ended=false every fragment is exactly read_length wide.
        assert!(!reads.is_empty());
        assert!(reads.iter().all(|(s, e)| e - s == read_length));
    }

    #[test]
    fn test_seed_rng() {
        let sequence = vec![0, 0, 2, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 2, 2, 2, 4, 4, 4, 4];
        let read_length = 10;
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let run1 = generate_fragments(
            sequence.len(),
            read_length,
            0,
            0,
            coverage,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &fragment_model,
            &mut rng,
        )
        .unwrap();

        // Reset to the same seed for reproducibility check.
        let mut rng2 = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let run2 = generate_fragments(
            sequence.len(),
            read_length,
            0,
            0,
            coverage,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &fragment_model,
            &mut rng2,
        )
        .unwrap();

        assert_eq!(run1, run2)
    }

    #[test]
    fn test_generate_reads_paired() {
        let seq_len = 100_000_usize;
        let read_length = 100;
        let coverage = 1;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let reads = generate_fragments(
            seq_len,
            read_length,
            0,
            0,
            coverage,
            true,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &fragment_model,
            &mut rng,
        )
        .unwrap();
        // Paired-end produces 2 reads per fragment, so fragment count is halved to keep
        // per-base read depth equal to the requested coverage.
        // Use ceiling division to match generate_fragments behavior.
        let expected = (seq_len * coverage).div_ceil(2 * read_length);
        assert!(!reads.is_empty());
        assert_eq!(
            reads.len(),
            expected,
            "Paired reads: expected {expected} fragments, got {}",
            reads.len()
        );
        assert!(
            reads.iter().all(|&(s, e)| e <= seq_len && e > s),
            "All fragments must be within bounds and non-empty"
        );
    }

    #[test]
    fn test_paired_pool_fully_populated() {
        // Fragment model with mean 600, st_dev 30 vs read_length 100 → near-zero rejection.
        // The pool should be exactly num_frags entries so that fragment size diversity is full.
        let seq_len = 50_000_usize;
        let read_length = 100;
        let coverage = 2;
        let num_frags = (seq_len * coverage).div_ceil(2 * read_length);
        let mut rng = NeatRng::new_from_seed(&vec!["pool-fill-test".to_string()]).unwrap();
        let fragment_model = FragmentLengthModel::new_normal(600.0, 30.0).unwrap();
        let reads = generate_fragments(
            seq_len,
            read_length,
            0,
            0,
            coverage,
            true,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &fragment_model,
            &mut rng,
        )
        .unwrap();
        // When rejection is near-zero the pool is full, so coverage should match exactly.
        assert_eq!(
            reads.len(),
            num_frags,
            "Expected full pool ({num_frags} fragments), got {}",
            reads.len()
        );
    }

    // ── generate_weighted_fragments ──────────────────────────────────────────

    #[test]
    fn test_weighted_fragments_short_region_returns_empty() {
        // region_len=50, read_length=30, max_del_len=10 → guard: 50 <= 30 + 20
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([A, C, G, T], 25).flatten().collect());
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let result = generate_weighted_fragments(
            &sequence_block,
            0,
            50,
            30,
            10,
            5,
            &GcBiasModel::default(),
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        assert!(
            result.is_empty(),
            "Region too short should return no fragments"
        );
    }

    #[test]
    fn test_weighted_fragments_zero_weight_region_returns_empty() {
        // All-A sequence (0% GC), model has weight only for 50% GC → all positions weight 0
        let sequence_block = make_sequence_block(vec![A; 500]);
        let mut weights = weights_with_value(0.0);
        weights[50] = 1.0;
        let model = GcBiasModel::from_weights(weights, 10).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let result = generate_weighted_fragments(
            &sequence_block,
            0,
            500,
            10,
            0,
            5,
            &model,
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        assert!(
            result.is_empty(),
            "All-zero-weight region should produce no fragments"
        );
    }

    #[test]
    fn test_weighted_fragments_produces_fragments_for_valid_region() {
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([A, C, G, T], 500).flatten().collect());
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let fragments = generate_weighted_fragments(
            &sequence_block,
            0,
            2000,
            100,
            0,
            10,
            &GcBiasModel::default(),
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        assert!(
            !fragments.is_empty(),
            "Valid region with uniform model should produce fragments"
        );
    }

    #[test]
    fn test_weighted_fragments_all_within_bounds() {
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([A, C, G, T], 1000).flatten().collect());
        let (region_start, region_end) = (200, 3800);
        let fragment_model = FragmentLengthModel::default().unwrap();
        let mut rng = make_rng();

        let fragments = generate_weighted_fragments(
            &sequence_block,
            region_start,
            region_end,
            100,
            0,
            5,
            &GcBiasModel::default(),
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        for (start, end) in &fragments {
            assert!(
                *start >= region_start,
                "start {} < region_start {}",
                start,
                region_start
            );
            assert!(
                *end <= region_end,
                "end {} > region_end {}",
                end,
                region_end
            );
        }
    }

    #[test]
    fn test_weighted_fragments_deterministic() {
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([A, C, G, T], 500).flatten().collect());
        let fragment_model = FragmentLengthModel::default().unwrap();

        let mut rng1 = make_rng();
        let run1 = generate_weighted_fragments(
            &sequence_block,
            0,
            2000,
            100,
            0,
            5,
            &GcBiasModel::default(),
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng1,
        )
        .unwrap();

        let mut rng2 = make_rng();
        let run2 = generate_weighted_fragments(
            &sequence_block,
            0,
            2000,
            100,
            0,
            5,
            &GcBiasModel::default(),
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng2,
        )
        .unwrap();

        assert_eq!(run1, run2, "Same seed must produce identical fragment sets");
    }

    #[test]
    fn test_weighted_fragments_normalize_true_generates_more_than_false() {
        // 50% GC sequence; model mean weight ≈ 0.5 → normalize should roughly double count.
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([A, C, G, T], 2500).flatten().collect());
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[100] = 1.0;
        let model = GcBiasModel::from_weights(weights, 100).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();

        let mut rng = make_rng();
        let unnormalized = generate_weighted_fragments(
            &sequence_block,
            0,
            10000,
            100,
            0,
            10,
            &model,
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        let mut rng = make_rng();
        let normalized = generate_weighted_fragments(
            &sequence_block,
            0,
            10000,
            100,
            0,
            10,
            &model,
            &fragment_model,
            true,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        assert!(
            normalized.len() > unnormalized.len(),
            "normalize=true ({}) should produce more fragments than normalize=false ({})",
            normalized.len(),
            unnormalized.len()
        );
    }

    #[test]
    fn test_weighted_fragments_near_zero_mean_weight_is_capped() {
        // 100% GC sequence; model weights[100]=1e-6 → mean_weight ≈ 1e-6
        // → normalize would inflate by ~1e6, capped at 100x.
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([C, G], 5000).flatten().collect());
        let mut weights = weights_with_value(0.0);
        weights[100] = 1e-6;
        weights[50] = 1.0; // ensures model is valid (at least one positive weight)
        let model = GcBiasModel::from_weights(weights, 100).unwrap();
        let fragment_model = FragmentLengthModel::default().unwrap();
        let read_len = 100;
        let target_coverage = 10;

        let mut rng = make_rng();
        let result = generate_weighted_fragments(
            &sequence_block,
            0,
            sequence_block.sequence.len(),
            read_len,
            0,
            target_coverage,
            &model,
            &fragment_model,
            true,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        let base_frags = (sequence_block.sequence.len() / read_len) * target_coverage;
        let cap = base_frags * 100;
        assert!(
            result.len() <= cap,
            "Fragment count {} should not exceed cap {}",
            result.len(),
            cap
        );
    }

    #[test]
    fn test_weighted_fragments_rejection_compensation() {
        // Fragment model with mean very close to read_length causes high rejection
        // (many sampled fragments fall below read_length). The while-loop must retry
        // until num_frags fragments are placed, not stop early.
        let read_length = 200;
        let sequence_block =
            make_sequence_block(std::iter::repeat_n([A, C, G, T], 5000).flatten().collect());
        let (region_start, region_end) = (0, sequence_block.sequence.len());
        let coverage = 5;
        let region_len = region_end - region_start;
        let num_frags_expected = (region_len * coverage).div_ceil(read_length);
        // mean=220 with st_dev=40 → P(frag < 200) ≈ 31% rejection rate
        let fragment_model = FragmentLengthModel::new_normal(220.0, 40.0).unwrap();
        let mut rng = make_rng();

        let result = generate_weighted_fragments(
            &sequence_block,
            region_start,
            region_end,
            read_length,
            0,
            coverage,
            &GcBiasModel::default(),
            &fragment_model,
            false,
            false,
            false,
            false, // keep_short (off in these fragment-placement tests)
            &mut rng,
        )
        .unwrap();

        assert_eq!(
            result.len(),
            num_frags_expected,
            "Expected {num_frags_expected} fragments despite high rejection rate, got {}",
            result.len()
        );
    }

    #[test]
    fn test_cover_dataset_fragment_larger_than_span_terminates() {
        // Regression test: when every fragment in the pool is larger than the contig,
        // cover_dataset must terminate rather than loop forever.
        // Mirrors the yeast long-read scenario: fragment_mean=30000, small scaffold.
        let span_length = 5_000;
        let read_length = 1_000;
        let coverage = 5;
        // Pool of fragments all larger than the span.
        let fragment_pool = vec![30_000_usize; 10];
        let mut rng = NeatRng::new_from_seed(&vec!["termination-test".to_string()]).unwrap();

        let target_count = span_length * coverage / read_length;
        let result = cover_dataset(
            span_length,
            read_length,
            target_count,
            0,
            fragment_pool,
            &mut rng,
        )
        .unwrap();

        // No fragments should be placed (none can fit), but the call must return.
        assert!(
            result.is_empty(),
            "No fragments should be placed when all exceed span"
        );
    }

    #[test]
    fn prototype_insertion_windows_cover_interior_and_tail() {
        let mut rng = make_rng();
        let windows =
            prototype_insertion_read_windows(600, 100, 12, vec![250; 12], &mut rng).unwrap();

        assert!(!windows.is_empty());
        assert!(
            windows.iter().any(|&(start, _)| start >= 200),
            "prototype must place a read window inside a long insertion"
        );
        assert!(
            windows.iter().any(|&(_, end)| end >= 570),
            "prototype must reach the far end of a long insertion"
        );
        assert!(
            windows
                .iter()
                .all(|&(start, end)| { end > start && end - start <= 100 && end <= 600 })
        );
    }

    #[test]
    fn test_cover_dataset_depth() {
        let span_length = 10000;
        let read_length = 100;
        let fragment_pool = vec![read_length; 1000];
        let coverage = 10;
        let target_count = span_length * coverage / read_length;
        let mut rng = make_rng();

        let fragments = cover_dataset(
            span_length,
            read_length,
            target_count,
            0,
            fragment_pool,
            &mut rng,
        )
        .unwrap();

        // Calculate actual depth at each position
        let mut depth = vec![0usize; span_length];
        for (start, end) in fragments {
            let s = start.saturating_sub(0); // read_start is 0
            let e = end.saturating_sub(0);
            for i in s..e {
                if i < span_length {
                    depth[i] += 1;
                }
            }
        }

        let avg_depth = depth.iter().sum::<usize>() as f64 / span_length as f64;
        println!(
            "Requested coverage: {}, Actual average depth: {}",
            coverage, avg_depth
        );

        // We expect it to be close to the requested coverage.
        assert!(
            avg_depth > (coverage as f64 * 0.9),
            "Coverage too low: {}",
            avg_depth
        );
    }

    #[test]
    fn test_cover_dataset_depth_with_fragments() {
        let span_length = 10000;
        let read_length = 100;
        let fragment_pool = vec![300, 400, 500, 600]; // much larger than read_length
        let coverage = 10;
        // In generate_fragments, num_frags is calculated as:
        // span_length * coverage / (2 * read_length) for paired-ended
        let target_count = span_length * coverage / (2 * read_length);
        let mut rng = make_rng();

        let fragments = cover_dataset(
            span_length,
            read_length,
            target_count,
            0,
            fragment_pool,
            &mut rng,
        )
        .unwrap();

        // Calculate actual depth at each position
        let mut depth = vec![0usize; span_length];
        for (start, end) in fragments {
            for i in start..end {
                if i < span_length {
                    depth[i] += 1;
                }
            }
        }

        let avg_depth = depth.iter().sum::<usize>() as f64 / span_length as f64;
        println!(
            "Requested coverage: {}, Actual average depth: {}, target_count: {}",
            coverage, avg_depth, target_count
        );
    }

    #[test]
    fn test_generate_fragments_depth_paired() {
        let span_length = 10000;
        let read_length = 100;
        let coverage = 10;
        let mut rng = make_rng();

        // Mock fragment model with mean 450, stdev 50
        let fragment_model = FragmentLengthModel::new_normal(450.0, 50.0).unwrap();

        let fragments = generate_fragments(
            span_length,
            read_length,
            0,
            0,
            coverage,
            true,  // paired_ended
            false, // long_reads
            false, // keep_short (off in these fragment-placement tests)
            &fragment_model,
            &mut rng,
        )
        .unwrap();

        // Calculate actual depth at each position
        // In paired-ended mode, each fragment (start, end) produces TWO reads of read_length.
        // One at start, one at end-read_length.
        let mut depth = vec![0usize; span_length];
        for (start, end) in fragments {
            // R1
            for i in start..start + read_length {
                if i < span_length {
                    depth[i] += 1;
                }
            }
            // R2
            let r2_start = end.saturating_sub(read_length);
            for i in r2_start..end {
                if i < span_length {
                    depth[i] += 1;
                }
            }
        }

        let avg_depth = depth.iter().sum::<usize>() as f64 / span_length as f64;
        println!(
            "Requested coverage: {}, Actual average depth (paired): {}",
            coverage, avg_depth
        );
    }

    #[test]
    fn test_generate_weighted_fragments_depth() {
        let span_length = 10000;
        let read_length = 100;
        let coverage = 10;
        let mut rng = make_rng();

        let sequence = vec![Nucleotide::A; span_length];
        let block = make_sequence_block(sequence);
        let fragment_model = FragmentLengthModel::new_normal(450.0, 50.0).unwrap();
        let gc_bias_model = GcBiasModel::default(); // Uniform

        let fragments = generate_weighted_fragments(
            &block,
            0,
            span_length,
            read_length,
            0,
            coverage,
            &gc_bias_model,
            &fragment_model,
            true,  // normalize
            true,  // paired_ended
            false, // long_reads
            false, // keep_short
            &mut rng,
        )
        .unwrap();

        let mut depth = vec![0usize; span_length];
        for (start, end) in fragments {
            for i in start..start + read_length {
                if i < span_length {
                    depth[i] += 1;
                }
            }
            let r2_start = end.saturating_sub(read_length);
            for i in r2_start..end {
                if i < span_length {
                    depth[i] += 1;
                }
            }
        }

        let avg_depth = depth.iter().sum::<usize>() as f64 / span_length as f64;
        println!(
            "Requested coverage: {}, Actual average depth (weighted uniform): {}",
            coverage, avg_depth
        );
    }

    #[test]
    fn test_generate_weighted_fragments_depth_with_bias() {
        let span_length = 10000;
        let read_length = 100;
        let coverage = 10;
        let mut rng = make_rng();

        // High GC sequence
        let sequence = vec![Nucleotide::G; span_length];
        let block = make_sequence_block(sequence);
        let fragment_model = FragmentLengthModel::new_normal(450.0, 50.0).unwrap();

        // Model that strongly prefers high GC (e.g. 1.0 at 100% GC, 0.1 at 0% GC)
        let mut weights = vec![0.1; 101];
        weights[100] = 1.0;
        let gc_bias_model = GcBiasModel::from_weights(weights, read_length).unwrap();

        let fragments = generate_weighted_fragments(
            &block,
            0,
            span_length,
            read_length,
            0,
            coverage,
            &gc_bias_model,
            &fragment_model,
            true,  // normalize=true should bring coverage back to 10
            true,  // paired_ended
            false, // long_reads
            false, // keep_short
            &mut rng,
        )
        .unwrap();

        let mut depth = vec![0usize; span_length];
        for (start, end) in fragments {
            for i in start..start + read_length {
                if i < span_length {
                    depth[i] += 1;
                }
            }
            let r2_start = end.saturating_sub(read_length);
            for i in r2_start..end {
                if i < span_length {
                    depth[i] += 1;
                }
            }
        }

        let avg_depth = depth.iter().sum::<usize>() as f64 / span_length as f64;
        println!(
            "Requested coverage: {}, Actual average depth (weighted biased + normalized): {}",
            coverage, avg_depth
        );
        assert!(
            avg_depth > (coverage as f64 * 0.9),
            "Normalization failed to maintain coverage: {}",
            avg_depth
        );
    }
}
