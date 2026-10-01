use crate::{
    eidolon_core::{
        models::{mutation_model::MutationModel, snp_trinuc_model::TrinucFrame},
        structs::{
            sequence_block::SequenceBlock,
            variants::{Variant, VariantType},
        },
    },
    gen_reads::errors::GenerateReadsError,
};
use eidolon_core::rng::NeatRng;
use log::debug;
use std::collections::HashMap;

/// Max mutable positions held in a per-position weight array at once, bounding the
/// memory of context-weighted SNP placement (~8 bytes/position → ~32 MB at 4M). The
/// contig is processed in chunks of this size so a large contig never allocates a
/// contig-sized cumulative.
const PLACEMENT_CHUNK: usize = 4_000_000;

/// Generates random variants for a contig.
///
/// `rate_segments` is a sorted, non-overlapping list of `(start, end, rate)` triples
/// covering the mutable positions of the contig. N-regions and excluded positions
/// (e.g. positions already occupied by input variants) must be absent from the list —
/// the caller is responsible for building a compact segment representation.
pub fn generate_variants(
    sequence_block: &SequenceBlock,
    rate_segments: &[(usize, usize, f64)],
    mutation_model: &MutationModel,
    num_mutations: usize,
    ploidy: usize,
    rng: &mut NeatRng,
) -> Result<Option<Vec<Variant>>, GenerateReadsError> {
    if rate_segments.is_empty() {
        return Ok(Some(Vec::new()));
    }

    // #372: weight WHERE mutations land by the local trinucleotide's SNP mutation
    // propensity w(ctx), so context-specific signatures (e.g. APOBEC's TpC hotspots)
    // reproduce in the output. When the model's context weights are uniform (default /
    // untrained model), placement is context-neutral and we take the fast
    // segment-arithmetic path — byte-identical to the pre-#372 behavior.
    let ctx_weights = mutation_model.context_weights()?;
    if context_is_uniform(&ctx_weights) {
        return generate_variants_uniform(
            sequence_block,
            rate_segments,
            mutation_model,
            num_mutations,
            ploidy,
            rng,
        );
    }

    // w(ctx) is a SNP propensity, so it drives SNP placement only. Indels stay
    // context-neutral: their positional biology (homopolymers / repeats / microhomology)
    // is not the SNP trinucleotide signature, and applying the SNP weight to them
    // concentrates indels and degrades indel recall (regression seen on cancer runs).
    // We therefore decide each mutation's type up front and force it in generate_mutation,
    // routing SNPs to context-weighted placement and indels to uniform placement.
    let seq = &sequence_block.sequence;
    let mean_w = ctx_weights.values().sum::<f64>() / (ctx_weights.len().max(1) as f64);
    let w_at = |p: usize| -> f64 {
        // Interior positions have a trinucleotide context; contig-edge positions fall
        // back to the mean weight (no ±1 flank available).
        if p >= 1 && p + 1 < seq.len() {
            let ctx = TrinucFrame::from(&[seq[p - 1], seq[p], seq[p + 1]]);
            *ctx_weights.get(&ctx).unwrap_or(&mean_w)
        } else {
            mean_w
        }
    };

    let mut snp_count = 0usize;
    let mut indel_types: Vec<VariantType> = Vec::new();
    for _ in 0..num_mutations {
        match mutation_model.variant_dist.sample(rng.random()?)? {
            VariantType::SNP => snp_count += 1,
            other => indel_types.push(other),
        }
    }

    let mut block_variants: Vec<Variant> = Vec::with_capacity(num_mutations);

    // ── indels: context-neutral, regional-rate placement (pre-#372 math) ────────────
    if !indel_types.is_empty() {
        let (seg_cum, seg_total) = segment_cumulative(rate_segments);
        for t in indel_types {
            let location = sample_segment_position(rate_segments, &seg_cum, seg_total, rng)?;
            debug!("location (indel, neutral): {}", location);
            block_variants.push(mutation_model.generate_mutation(
                seq,
                location,
                ploidy,
                Some(t),
                rng,
            )?);
        }
    }

    // ── SNPs: context-weighted, memory-bounded ──────────────────────────────────────
    // Pass 1: chunk the mutable positions (≤ PLACEMENT_CHUNK each, never spanning a
    // segment) and record each chunk's SNP weight = Σ rate·w(ctx). O(#chunks) memory.
    if snp_count > 0 {
        let mut chunks: Vec<(usize, usize, f64)> = Vec::new(); // (start, end, rate)
        let mut chunk_cum: Vec<f64> = Vec::new();
        let mut acc = 0.0_f64;
        for &(start, end, rate) in rate_segments {
            let mut cs = start;
            while cs < end {
                let ce = (cs + PLACEMENT_CHUNK).min(end);
                let mut w = 0.0_f64;
                for p in cs..ce {
                    w += rate * w_at(p);
                }
                acc += w;
                chunks.push((cs, ce, rate));
                chunk_cum.push(acc);
                cs = ce;
            }
        }
        let chunk_total = acc;
        if chunk_total > 0.0 && !chunks.is_empty() {
            // Allocate the SNP placements across chunks by weight.
            let mut per_chunk = vec![0usize; chunks.len()];
            for _ in 0..snp_count {
                let target = rng.random()? * chunk_total;
                let ci = chunk_cum
                    .partition_point(|&c| c <= target)
                    .min(chunks.len() - 1);
                per_chunk[ci] += 1;
            }
            // Pass 2: for each needed chunk build one bounded per-position cumulative
            // and sample its allocation.
            for (ci, &k) in per_chunk.iter().enumerate() {
                if k == 0 {
                    continue;
                }
                let (cs, ce, rate) = chunks[ci];
                let mut positions: Vec<usize> = Vec::with_capacity(ce - cs);
                let mut cum: Vec<f64> = Vec::with_capacity(ce - cs);
                let mut a = 0.0_f64;
                for p in cs..ce {
                    a += rate * w_at(p);
                    positions.push(p);
                    cum.push(a);
                }
                if a <= 0.0 || positions.is_empty() {
                    continue;
                }
                for _ in 0..k {
                    let target = rng.random()? * a;
                    let idx = cum
                        .partition_point(|&c| c <= target)
                        .min(positions.len() - 1);
                    let location = positions[idx];
                    debug!("location (context-weighted SNP): {}", location);
                    block_variants.push(mutation_model.generate_mutation(
                        seq,
                        location,
                        ploidy,
                        Some(VariantType::SNP),
                        rng,
                    )?);
                }
            }
        }
    }

    Ok(Some(block_variants))
}

/// True when the per-context weights carry no signal (uniform / empty) — the default
/// context-flat case. Placement then stays context-neutral (pre-#372 behavior).
fn context_is_uniform(weights: &HashMap<TrinucFrame, f64>) -> bool {
    let mut vals = weights.values().copied();
    let Some(first) = vals.next() else {
        return true;
    };
    let (mut lo, mut hi) = (first, first);
    for v in vals {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    (hi - lo) <= 1e-9 * hi.max(1e-12)
}

/// Cumulative `length × rate` weights over the segments (for context-neutral placement).
fn segment_cumulative(rate_segments: &[(usize, usize, f64)]) -> (Vec<f64>, f64) {
    let mut cum: Vec<f64> = Vec::with_capacity(rate_segments.len());
    let mut acc = 0.0_f64;
    for &(start, end, rate) in rate_segments {
        acc += (end - start) as f64 * rate;
        cum.push(acc);
    }
    (cum, acc)
}

/// Sample one position ∝ segment rate, uniform within the chosen segment (pre-#372 math).
fn sample_segment_position(
    rate_segments: &[(usize, usize, f64)],
    cum: &[f64],
    total: f64,
    rng: &mut NeatRng,
) -> Result<usize, GenerateReadsError> {
    let target = rng.random()? * total;
    let seg_idx = cum
        .partition_point(|&c| c <= target)
        .min(rate_segments.len() - 1);
    let (seg_start, seg_end, seg_rate) = rate_segments[seg_idx];
    let weight_before = if seg_idx > 0 { cum[seg_idx - 1] } else { 0.0 };
    let weight_in_seg = (target - weight_before).max(0.0);
    let pos_in_seg = (weight_in_seg / seg_rate) as usize;
    Ok((seg_start + pos_in_seg).min(seg_end - 1))
}

/// Context-neutral placement for all types: positions weighted only by the regional
/// `rate_segments`, uniform within a segment. The original (pre-#372) algorithm, used
/// whenever the model's context weights are uniform.
fn generate_variants_uniform(
    sequence_block: &SequenceBlock,
    rate_segments: &[(usize, usize, f64)],
    mutation_model: &MutationModel,
    num_mutations: usize,
    ploidy: usize,
    rng: &mut NeatRng,
) -> Result<Option<Vec<Variant>>, GenerateReadsError> {
    let mut block_variants: Vec<Variant> = Vec::with_capacity(num_mutations);
    let (cum, total) = segment_cumulative(rate_segments);
    for _ in 0..num_mutations {
        let location = sample_segment_position(rate_segments, &cum, total, rng)?;
        debug!("location: {}", location);
        block_variants.push(mutation_model.generate_mutation(
            &sequence_block.sequence,
            location,
            ploidy,
            None,
            rng,
        )?);
    }
    Ok(Some(block_variants))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eidolon_core::{
        models::mutation_model::MutationModel,
        structs::sequence_block::{RegionType, SequenceBlock, SequenceMap},
    };
    use eidolon_core::rng::NeatRng;
    use eidolon_core::structs::nucleotides::Nucleotide;

    fn make_block(len: usize) -> SequenceBlock {
        let sequence: Vec<Nucleotide> = (0..len)
            .map(|i| match i % 4 {
                0 => Nucleotide::A,
                1 => Nucleotide::C,
                2 => Nucleotide::G,
                _ => Nucleotide::T,
            })
            .collect();
        SequenceBlock {
            contig: "test".to_string(),
            ref_start: 0,
            ref_end: len,
            sequence,
            sequence_map: vec![SequenceMap::from(RegionType::NonNRegion, 0, len)],
        }
    }

    fn make_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap()
    }

    #[test]
    fn test_generate_variants_count() {
        let block = make_block(200);
        let segments = vec![(0usize, 200usize, 1.0f64)];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let result = generate_variants(&block, &segments, &model, 5, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(result.len(), 5);
    }

    #[test]
    fn context_weighted_placement_concentrates_snps_at_upweighted_context() {
        use eidolon_core::structs::nucleotides::Nucleotide::{A, C, G};

        // make_block is ACGTACGT…, so C positions (p % 4 == 1) carry the context ACG.
        // Up-weight ACG heavily → #372 SNP placement should concentrate there;
        // context-neutral placement would put only ~1/4 at those positions. All-SNP model.
        let mut trinuc_freq: HashMap<TrinucFrame, f64> = HashMap::new();
        trinuc_freq.insert(TrinucFrame::from((A, C, G)), 1000.0);
        let model = MutationModel::from_raw_data(
            0.01,
            0.0,
            vec![1.0, 0.0, 0.0], // all SNP
            trinuc_freq,
            HashMap::new(),
            vec![],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();

        let block = make_block(400);
        let segments = vec![(0usize, 400usize, 1.0f64)];
        let mut rng = make_rng();
        let variants = generate_variants(&block, &segments, &model, 200, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(variants.len(), 200);
        let at_acg = variants.iter().filter(|v| v.location % 4 == 1).count();
        let frac = at_acg as f64 / variants.len() as f64;
        assert!(
            frac > 0.8,
            "expected >80% of SNP placements at the up-weighted ACG context, got {:.2} ({}/{})",
            frac,
            at_acg,
            variants.len()
        );
    }

    #[test]
    fn indels_stay_context_neutral_even_with_a_peaked_snp_context() {
        use eidolon_core::structs::nucleotides::Nucleotide::{A, C, G};

        // Peak the ACG SNP context but make every variant an insertion. Because w(ctx)
        // drives SNP placement only, the indels must NOT concentrate at ACG (they stay
        // roughly uniform ~1/4 at C positions), guarding the somatic_indel_recall fix.
        let mut trinuc_freq: HashMap<TrinucFrame, f64> = HashMap::new();
        trinuc_freq.insert(TrinucFrame::from((A, C, G)), 1000.0);
        let model = MutationModel::from_raw_data(
            0.01,
            0.0,
            vec![0.0, 1.0, 0.0], // all insertions
            trinuc_freq,
            HashMap::new(),
            vec![1, 2],
            vec![1.0, 1.0],
            vec![1, 2],
            vec![1.0, 1.0],
        )
        .unwrap();

        let block = make_block(400);
        let segments = vec![(0usize, 400usize, 1.0f64)];
        let mut rng = make_rng();
        let variants = generate_variants(&block, &segments, &model, 200, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert_eq!(variants.len(), 200);
        let at_acg = variants.iter().filter(|v| v.location % 4 == 1).count();
        let frac = at_acg as f64 / variants.len() as f64;
        assert!(
            frac < 0.45,
            "indels should be context-neutral (~0.25 at ACG), not concentrated; got {:.2}",
            frac
        );
    }

    #[test]
    fn test_generate_variants_positions_in_range() {
        let block = make_block(200);
        let segments = vec![(0usize, 200usize, 1.0f64)];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let variants = generate_variants(&block, &segments, &model, 10, 2, &mut rng)
            .unwrap()
            .unwrap();
        for v in variants {
            assert!(
                v.location < 200,
                "Variant location {} out of range",
                v.location
            );
        }
    }

    #[test]
    fn test_generate_variants_empty_segments_returns_empty() {
        let block = make_block(200);
        let segments: Vec<(usize, usize, f64)> = vec![];
        let model = MutationModel::default().unwrap();
        let mut rng = make_rng();
        let result = generate_variants(&block, &segments, &model, 5, 2, &mut rng)
            .unwrap()
            .unwrap();
        assert!(
            result.is_empty(),
            "Empty segments should produce no variants"
        );
    }

    #[test]
    fn test_generate_variants_deterministic() {
        let block = make_block(200);
        let segments = vec![(0usize, 200usize, 1.0f64)];
        let model = MutationModel::default().unwrap();
        let result1 = generate_variants(&block, &segments, &model, 5, 2, &mut make_rng())
            .unwrap()
            .unwrap();
        let result2 = generate_variants(&block, &segments, &model, 5, 2, &mut make_rng())
            .unwrap()
            .unwrap();
        for (v1, v2) in result1.iter().zip(result2.iter()) {
            assert_eq!(v1.location, v2.location);
            assert_eq!(v1.variant_type, v2.variant_type);
        }
    }
}
