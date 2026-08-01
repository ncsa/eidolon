//! Statistical model for generating symbolic / structural variants
//! (`<DEL>` / `<DUP>` / `<CNV>` / `<BND>` / `<INV>` / `<INS>`) de novo from a learned distribution.
//!
//! `SvModel` is fit by [`SvModel::fit_from_observations`] (called from
//! `gen-mut-model` after the SNP/indel pass completes) and lives on
//! [`crate::models::mutation_model::MutationModel`] as
//! [`Option<SvModel>`] — `None` whenever the training VCF didn't carry
//! enough usable SV observations, or whenever the model JSON was
//! written by a pre-v1.10 build (the field defaults to `None` via
//! `#[serde(default)]`).
//!
//! Supported types: `<DEL>` / `<DUP>` / `<CNV>` / `<BND>` / `<INV>` / `<INS>`.
//! `gen-reads` consumes the model via [`SvModel::sample_variants`],
//! gated behind the `sv_rate_scale` config knob (default 0.0 — opt-in).

use crate::rng::{NeatRng, NeatRngError};
use crate::structs::distributions::NormalDistribution;
use crate::structs::nucleotides::Nucleotide;
use crate::structs::variants::{
    AlternateType, Genotype, Provenance, SvData, SvType, Variant, VariantType,
};
use log::{info, warn};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Minimum SV length (in bases) treated as an SV rather than an indel.
/// Matches the conventional cutoff used by Manta, DELLY, Lumpy, GRIDSS,
/// and the rest of the SV-calling ecosystem.
///
/// [`SvModel::fit_from_observations`] drops training records below this
/// span — they're indels, not SVs.
pub const MIN_SV_LENGTH_BP: usize = 50;

/// Per-type sample count required to fit a length distribution.
/// Log-normal MLE on a single observation yields `sigma = 0` — a
/// degenerate distribution. Types with fewer observations than this are
/// dropped from the trained model rather than carrying a fake `sigma`.
const MIN_OBS_FOR_LENGTH_FIT: usize = 2;

/// Fallback copy-number distribution used by [`SvModel::sample_variants`]
/// when the trained model's `cnv_copy_number_distribution` is empty but
/// `Cnv` is still in `type_probabilities`. This happens when the
/// training corpus had `<CNV>` records but stored per-sample CN in
/// FORMAT/SAMPLE rather than `INFO/CN` (gnomAD-SV is the canonical
/// example — CN varies across the cohort, so the sites VCF leaves
/// `INFO/CN` unset). Without a fallback the sampler would emit `<CNV>`
/// records without `copy_number`, and gen-reads would warn-and-skip
/// coverage modulation for each one.
///
/// Values mirror `sv_model_defaults::default_sv_model`'s
/// `cnv_copy_number_distribution` — duplicated here as a local constant
/// to avoid a `structs::sv_model` ↔ `models::sv_model_defaults` module
/// dependency cycle. Keep the two in sync if either drifts.
const FALLBACK_CN_DISTRIBUTION: &[(u32, f64)] =
    &[(0, 0.15), (1, 0.35), (3, 0.30), (4, 0.15), (5, 0.05)];

/// Learned statistical model for symbolic-SV generation.
///
/// The five fields together specify everything a sampler needs: how
/// many SVs per base (`per_base_rate`), which type each one is
/// (`type_probabilities`), how long it is (`length_log_normal`, keyed
/// by type), what copy number a `<CNV>` carries
/// (`cnv_copy_number_distribution`), and whether it's homozygous
/// (`homozygous_frequency`).
#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
pub struct SvModel {
    /// Per-base SV rate across the training corpus
    /// (= count of fit-eligible SVs / sum of training contig lengths
    /// the trainer iterated over).
    pub per_base_rate: f64,

    /// Probability that a generated SV is of a given type. Keys are
    /// limited to `SvType::Del` / `SvType::Dup` / `SvType::Cnv` / `SvType::Bnd` / `SvType::Inv` / `SvType::Ins`;
    /// probabilities sum to ~1.0 across present keys. A type with
    /// fewer than [`MIN_OBS_FOR_LENGTH_FIT`] observations is dropped
    /// from this map (its length distribution couldn't be fit), and
    /// the remaining probabilities are renormalized.
    pub type_probabilities: HashMap<SvType, f64>,

    /// Length distribution per SV type — `(mu, sigma)` of a log-normal
    /// in log-space, fit by maximum likelihood on `ln(span)`. Same
    /// key set as `type_probabilities`.
    ///
    /// Conditioning length on type is the cheapest concession to
    /// biological realism: DELs trend smaller than DUPs / CNVs, and a
    /// single global length distribution would mis-bias whichever type
    /// the training corpus had less of.
    pub length_log_normal: HashMap<SvType, (f64, f64)>,

    /// Empirical copy-number distribution observed on `<CNV>` records
    /// (CN value → probability). Probabilities sum to ~1.0 if any
    /// CNVs in the training corpus carried `INFO/CN`; empty otherwise.
    pub cnv_copy_number_distribution: HashMap<u32, f64>,

    /// Probability that a generated SV is homozygous (vs heterozygous),
    /// estimated as the homozygous fraction across all fit-eligible
    /// observations of all types.
    pub homozygous_frequency: f64,
}

impl SvModel {
    /// Returns `true` when the model has enough learned signal to drive
    /// generation — at least one type probability and a positive
    /// per-base rate. Guards against partially-constructed or
    /// hand-rolled models reaching a sampler that would divide by zero
    /// or call a weighted pick on an empty distribution.
    pub fn is_usable(&self) -> bool {
        self.per_base_rate > 0.0 && !self.type_probabilities.is_empty()
    }

    /// Fit an `SvModel` from a set of observed SV records.
    ///
    /// Returns `None` when there isn't enough signal to build a usable
    /// model: zero in-scope observations, `total_reflen == 0`, or every
    /// type fell below [`MIN_OBS_FOR_LENGTH_FIT`]. Callers should treat
    /// `None` as "this corpus doesn't support de novo SV generation"
    /// and leave `MutationModel.sv_model` unset.
    ///
    /// Filtering applied (in order):
    ///   - records whose ALT isn't `AlternateType::Symbolic` are
    ///     ignored;
    ///   - records with no derivable span/length (no `END`, no `SVLEN`) are
    ///     dropped;
    ///   - records below [`MIN_SV_LENGTH_BP`] are dropped (they're
    ///     indels, not SVs), except for `<BND>` which uses a nominal
    ///     1bp length;
    ///   - unknown tags are dropped;
    ///   - types with fewer than [`MIN_OBS_FOR_LENGTH_FIT`] surviving
    ///     observations are dropped from the model entirely, so the
    ///     resulting `type_probabilities` and `length_log_normal` have
    ///     the same key set.
    ///
    /// Each filter stage logs a count so a "why is my SV model None?"
    /// debug is grep-able in the trainer output.
    pub fn fit_from_observations(observations: &[Variant], total_reflen: usize) -> Option<SvModel> {
        if total_reflen == 0 {
            warn!("SvModel: total_reflen is zero, cannot estimate per-base rate");
            return None;
        }

        let mut dropped_non_symbolic = 0usize;
        let mut dropped_unsupported_type = 0usize;
        let mut dropped_no_span = 0usize;
        let mut dropped_below_min_length = 0usize;

        // Per-type accumulators of (span, &Variant) so we can pull CN
        // and genotype back out without re-walking the source.
        let mut by_type: HashMap<SvType, Vec<(usize, &Variant)>> = HashMap::new();

        for v in observations {
            let sv = match v.alternate.as_symbolic() {
                Some(sv) => sv,
                None => {
                    dropped_non_symbolic += 1;
                    continue;
                }
            };
            match sv.sv_type {
                SvType::Del
                | SvType::Dup
                | SvType::Cnv
                | SvType::Bnd
                | SvType::Inv
                | SvType::Ins => {}
                _ => {
                    dropped_unsupported_type += 1;
                    continue;
                }
            }
            let length = match sv.sv_type {
                SvType::Bnd => Some(1), // BNDs have no length in the traditional sense; use 1 for stats
                _ => sv.event_length(v.location),
            };
            let length = match length {
                Some(s) if s > 0 => s,
                _ => {
                    dropped_no_span += 1;
                    continue;
                }
            };
            if sv.sv_type != SvType::Bnd && length < MIN_SV_LENGTH_BP {
                dropped_below_min_length += 1;
                continue;
            }
            by_type.entry(sv.sv_type).or_default().push((length, v));
        }

        // Drop types that don't have enough observations for a real
        // length fit. Renormalization happens below over what's left.
        let mut dropped_thin_types = 0usize;
        by_type.retain(|_, obs| {
            if obs.len() < MIN_OBS_FOR_LENGTH_FIT {
                dropped_thin_types += obs.len();
                false
            } else {
                true
            }
        });

        if dropped_non_symbolic > 0 {
            // Trainers should hand us symbolic-only records, but the
            // fitter is tolerant — a stray SNP just gets ignored.
            info!(
                "SvModel: ignored {} non-symbolic observation(s) (literal-base ALTs)",
                dropped_non_symbolic
            );
        }
        if dropped_unsupported_type > 0 {
            info!(
                "SvModel: dropped {} unknown record(s) — not yet generatable",
                dropped_unsupported_type
            );
        }
        if dropped_no_span > 0 {
            info!(
                "SvModel: dropped {} symbolic SV(s) with no INFO/END or INFO/SVLEN",
                dropped_no_span
            );
        }
        if dropped_below_min_length > 0 {
            info!(
                "SvModel: dropped {} symbolic SV(s) below {}bp (treated as indels, not SVs)",
                dropped_below_min_length, MIN_SV_LENGTH_BP
            );
        }
        if dropped_thin_types > 0 {
            info!(
                "SvModel: dropped {} observation(s) from types with <{} samples \
                 (can't fit a length distribution)",
                dropped_thin_types, MIN_OBS_FOR_LENGTH_FIT
            );
        }

        if by_type.is_empty() {
            warn!(
                "SvModel: no SV type cleared the {}-observation fit bar; \
                 returning None (MutationModel.sv_model stays unset)",
                MIN_OBS_FOR_LENGTH_FIT
            );
            return None;
        }

        let total_usable: usize = by_type.values().map(|v| v.len()).sum();
        let n = total_usable as f64;

        let type_probabilities: HashMap<SvType, f64> = by_type
            .iter()
            .map(|(k, obs)| (*k, obs.len() as f64 / n))
            .collect();

        // Log-normal MLE: μ = mean(ln x), σ² = mean((ln x − μ)²).
        let length_log_normal: HashMap<SvType, (f64, f64)> = by_type
            .iter()
            .map(|(sv_type, obs)| {
                let log_lens: Vec<f64> = obs.iter().map(|(s, _)| (*s as f64).ln()).collect();
                let mu = log_lens.iter().sum::<f64>() / log_lens.len() as f64;
                let variance =
                    log_lens.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / log_lens.len() as f64;
                (*sv_type, (mu, variance.sqrt()))
            })
            .collect();

        // CN distribution from <CNV> records that carry INFO/CN.
        let mut cn_counts: HashMap<u32, usize> = HashMap::new();
        if let Some(cnv_obs) = by_type.get(&SvType::Cnv) {
            for (_, v) in cnv_obs {
                if let Some(sv) = v.alternate.as_symbolic()
                    && let Some(cn) = sv.copy_number
                {
                    *cn_counts.entry(cn).or_default() += 1;
                }
            }
        }
        let cn_total: usize = cn_counts.values().sum();
        let cnv_copy_number_distribution: HashMap<u32, f64> = if cn_total > 0 {
            cn_counts
                .iter()
                .map(|(k, c)| (*k, *c as f64 / cn_total as f64))
                .collect()
        } else {
            HashMap::new()
        };

        let hom_count = by_type
            .values()
            .flat_map(|obs| obs.iter())
            .filter(|(_, v)| matches!(v.genotype, Genotype::Homozygous))
            .count();
        let homozygous_frequency = hom_count as f64 / n;

        let per_base_rate = n / total_reflen as f64;

        info!(
            "SvModel: fit from {} observation(s) across {} type(s); \
             per_base_rate = {:.3e}, homozygous_frequency = {:.3}",
            total_usable,
            by_type.len(),
            per_base_rate,
            homozygous_frequency
        );

        Some(SvModel {
            per_base_rate,
            type_probabilities,
            length_log_normal,
            cnv_copy_number_distribution,
            homozygous_frequency,
        })
    }

    /// Sample de novo symbolic-SV variants for a single contig, gated by
    /// the caller-supplied `rate_scale`.
    ///
    /// Algorithm:
    ///   1. Draw `n` from `Poisson(per_base_rate × contig_len × rate_scale)`.
    ///   2. For each: weighted-pick `SvType`; draw `length` from the
    ///      type's log-normal, rejecting draws outside
    ///      `[MIN_SV_LENGTH_BP, contig_len / 4]`; uniform-pick an anchor
    ///      position; sample `INFO/CN` for `<CNV>`; Bernoulli for
    ///      genotype using `homozygous_frequency`.
    ///   3. Reject candidates whose affected range overlaps any SV
    ///      already in `existing_svs` or any previously-accepted draw
    ///      in this batch — retry up to `10 × n` times then warn-and-stop.
    ///
    /// Returns an empty `Vec` when generation is disabled (any of:
    /// `rate_scale ≤ 0`, model is not [`is_usable`](Self::is_usable),
    /// or `contig_len < 2 × MIN_SV_LENGTH_BP`). Sampler failures
    /// (Poisson runaway, log-normal parameter blowups, etc.) log a
    /// `warn!` and return whatever was produced so far — gen-reads
    /// should keep generating SNP/indel reads rather than abort on a
    /// degenerate SV model.
    ///
    /// Variant positions and `SvData.end` are populated in the same
    /// 0-based-location-with-1-based-END convention used by
    /// `Variant::from_file` so the downstream `build_coverage_multipliers`
    /// path treats user-supplied and de novo SVs identically.
    ///
    /// `sequence` supplies the reference bases used to fill each
    /// generated record's `REF` field (one base, at the anchor for
    /// `<DEL>` or the first affected base for `<DUP>` / `<CNV>`). If
    /// the picked base is `N`, the record is rejected and re-drawn so
    /// we don't anchor SVs inside reference gaps.
    pub fn sample_variants(
        &self,
        contig_name: &str,
        contig_len: usize,
        existing_svs: &[Variant],
        sequence: &[Nucleotide],
        ploidy: usize,
        rate_scale: f64,
        max_length_fraction: f64,
        rng: &mut NeatRng,
    ) -> Vec<Variant> {
        if rate_scale <= 0.0 || !rate_scale.is_finite() {
            return Vec::new();
        }
        if !self.is_usable() {
            return Vec::new();
        }
        if contig_len < MIN_SV_LENGTH_BP.saturating_mul(2) {
            return Vec::new();
        }

        let lambda = self.per_base_rate * contig_len as f64 * rate_scale;
        if !lambda.is_finite() || lambda <= 0.0 {
            return Vec::new();
        }
        let n_sv = match sample_poisson(lambda, rng) {
            Ok(n) => n as usize,
            Err(e) => {
                warn!("SvModel: Poisson sampler error (λ={lambda}): {e}; no SVs emitted");
                return Vec::new();
            }
        };
        if n_sv == 0 {
            return Vec::new();
        }

        // Sorted type list + cumulative weights for repeatable picking
        // via NeatRng::random() (avoids depending on rand::WeightedIndex's
        // RNG trait shape).
        let mut types: Vec<SvType> = self.type_probabilities.keys().copied().collect();
        types.sort_by_key(|t| *t as u8);
        let type_weights: Vec<f64> = types
            .iter()
            .map(|t| self.type_probabilities[t].max(0.0))
            .collect();
        let type_cum = cumulative_normalized(&type_weights);
        if type_cum.is_empty() {
            return Vec::new();
        }

        // CN values + cumulative weights for sampling `<CNV>` copy
        // numbers. If the trained distribution is empty but `Cnv` is in
        // the model's type list, fall back to `FALLBACK_CN_DISTRIBUTION`
        // so sampled `<CNV>` records carry a sensible CN — otherwise
        // gen-reads warn-and-skips depth modulation on every one of them
        // (which is what gnomAD-SV-trained models hit in practice).
        let cnv_in_model = self.type_probabilities.contains_key(&SvType::Cnv);
        let cn_dist_is_empty = self.cnv_copy_number_distribution.is_empty();
        let using_cn_fallback = cnv_in_model && cn_dist_is_empty;
        if using_cn_fallback {
            info!(
                "SvModel: trained `cnv_copy_number_distribution` is empty but Cnv is \
                 in the type list — falling back to the bundled default CN distribution \
                 for sampled <CNV> records (this is what you want when training from a \
                 sites-only VCF that puts CN in FORMAT instead of INFO, like gnomAD-SV)."
            );
        }
        let (cn_values, cn_cum): (Vec<u32>, Vec<f64>) = if using_cn_fallback {
            let values: Vec<u32> = FALLBACK_CN_DISTRIBUTION.iter().map(|(c, _)| *c).collect();
            let weights: Vec<f64> = FALLBACK_CN_DISTRIBUTION.iter().map(|(_, p)| *p).collect();
            (values, cumulative_normalized(&weights))
        } else {
            let mut values: Vec<u32> = self.cnv_copy_number_distribution.keys().copied().collect();
            values.sort_unstable();
            let weights: Vec<f64> = values
                .iter()
                .map(|c| self.cnv_copy_number_distribution[c].max(0.0))
                .collect();
            (values, cumulative_normalized(&weights))
        };

        // Existing affected ranges seed the overlap set; new draws are
        // appended as accepted.
        let mut occupied: Vec<(usize, usize)> = existing_svs
            .iter()
            .filter_map(|v| affected_range_for_existing(v, contig_len))
            .collect();

        // #229: SV length cap = contig_len * max_length_fraction. Default 0.25
        // (contig_len/4) keeps overlap rejection tractable on chromosome-scale
        // references; raise toward 1.0 for small contigs (bacteria/viruses,
        // where 0.25 sits below the default DEL median and starves DEL/DUP/INV)
        // or native aneuploidy. A non-finite / non-positive value falls back to
        // 0.25 (the historical behavior).
        let frac = if max_length_fraction.is_finite() && max_length_fraction > 0.0 {
            max_length_fraction
        } else {
            0.25
        };
        let upper_len = ((contig_len as f64 * frac) as usize).max(MIN_SV_LENGTH_BP);
        // Bigger SVs are harder to place without overlap, so scale the retry
        // budget with the fraction relative to the 0.25 baseline.
        let retry_mult = (frac / 0.25).max(1.0);
        let max_retries = ((n_sv as f64 * 10.0 * retry_mult) as usize).max(20);
        let mut out: Vec<Variant> = Vec::with_capacity(n_sv);
        let mut retries = 0usize;
        // Count placed *events*, not emitted records: a BND is one event that
        // emits a mate PAIR (two records, see below). Gating the loop on
        // out.len() would let each BND consume two slots and halve the number of
        // SVs actually placed. For every other type events == out.len().
        let mut events = 0usize;

        while events < n_sv && retries < max_retries {
            retries += 1;

            let sv_type = match weighted_pick_with(&types, &type_cum, rng) {
                Some(t) => *t,
                None => break,
            };

            let (mu, sigma) = match self.length_log_normal.get(&sv_type) {
                Some(p) => *p,
                None => continue, // type without a length fit shouldn't occur given fit_from_observations
            };
            let length = match sample_log_normal_usize(mu, sigma, rng) {
                Ok(l) => l,
                Err(_) => continue,
            };
            if sv_type != SvType::Bnd && (length < MIN_SV_LENGTH_BP || length > upper_len) {
                continue;
            }

            let (affected_start, affected_end, anchor_0based, sv_end_1based) =
                match place_sv(sv_type, length, contig_len, rng) {
                    Some(p) => p,
                    None => continue,
                };

            if occupied
                .iter()
                .any(|&(s, e)| ranges_overlap(affected_start, affected_end, s, e))
            {
                continue;
            }

            let anchor_base = match sequence.get(anchor_0based) {
                Some(&b) if b != Nucleotide::N => b,
                _ => continue, // anchor falls in an N gap; try again
            };

            // #224: the single-base anchor check above is too loose. It
            // passes positions where the anchor is a non-N base sitting
            // right next to a long N-tract (centromere/telomere
            // boundaries). The chimeric pair generation then samples a
            // read_len-sized piece extending into the N-tract, producing
            // reads BWA can't place. Widen the check to a ±200 bp window
            // around the anchor — covers a typical PE fragment's local
            // piece. For DEL/DUP/CNV/INV the END side also needs to be
            // alignable since the right chimeric piece anchors there;
            // we check `sv_end_1based - 1` (= 0-based last affected base)
            // for those types.
            const ALIGNABLE_WINDOW: usize = 200;
            if !anchor_window_alignable(sequence, anchor_0based, ALIGNABLE_WINDOW) {
                continue;
            }
            if sv_type != SvType::Bnd
                && !anchor_window_alignable(
                    sequence,
                    sv_end_1based
                        .saturating_sub(1)
                        .min(sequence.len().saturating_sub(1)),
                    ALIGNABLE_WINDOW,
                )
            {
                continue;
            }

            let copy_number = if sv_type == SvType::Cnv && !cn_values.is_empty() {
                weighted_pick_with(&cn_values, &cn_cum, rng).copied()
            } else {
                None
            };

            let genotype = match rng.gen_bool(self.homozygous_frequency.clamp(0.0, 1.0)) {
                Ok(true) => Genotype::Homozygous,
                Ok(false) => Genotype::Heterozygous,
                Err(_) => Genotype::Heterozygous,
            };
            let genotype_str = genotype_string(&genotype, ploidy);

            // De novo <INS> records are emitted as LITERAL Insertion variants
            // (REF="A", ALT="A<novel-bases>") rather than symbolic <INS>. That
            // routes them through gen-reads' existing literal-insertion
            // machinery — fragment-loop reads spanning the locus naturally
            // pick up the inserted bases, and the output VCF carries the
            // resolved sequence (matching what real callers like Manta emit
            // when they can resolve the insertion). The novel sequence is
            // drawn from the single-base composition of a ±250 bp window
            // around the anchor; #190 calls out trinucleotide-context
            // sampling as a possible refinement.
            if sv_type == SvType::Ins {
                let novel = match sample_novel_insertion_bases(sequence, anchor_0based, length, rng)
                {
                    Ok(bases) => bases,
                    Err(_) => continue, // RNG hiccup; retry the slot
                };
                let mut alt_bases = Vec::with_capacity(length + 1);
                alt_bases.push(anchor_base);
                alt_bases.extend(novel);

                occupied.push((affected_start, affected_end));
                out.push(Variant {
                    variant_type: VariantType::Insertion,
                    location: anchor_0based,
                    reference: vec![anchor_base],
                    alternate: AlternateType::Literal(alt_bases),
                    genotype_str: genotype_str.clone(),
                    genotype,
                    allele_fraction: None,
                    id: None,
                    quality_score: None,
                    filter: None,
                    // No symbolic INFO needed — the SVLEN is implicit in
                    // the REF/ALT length difference.
                    info: None,
                    format: vec!["GT".to_string()],
                    sample: vec![genotype_str],
                    provenance: Provenance::Denovo,
                });
                continue;
            }

            let (raw_alt, m_contig, m_pos) = match sv_type {
                SvType::Del => ("<DEL>".to_string(), None, None),
                SvType::Dup => ("<DUP>".to_string(), None, None),
                SvType::Cnv => ("<CNV>".to_string(), None, None),
                SvType::Inv => ("<INV>".to_string(), None, None),
                SvType::Bnd => {
                    // For de novo BNDs, generate an intra-chromosomal translocation
                    // to a random mate position on the same contig. Retry up to
                    // 32 times if the picked mate position falls in an N-tract
                    // (so the BND truth record is recoverable by callers — see
                    // #224 for the v1.13.0 finding where ~37% of de novo BNDs
                    // landed in unalignable mate positions).
                    let mut mp_candidate: usize = 0;
                    let mut found = false;
                    for _ in 0..32 {
                        let candidate = rng
                            .range_i64(1, contig_len.saturating_sub(1) as i64)
                            .ok()
                            .unwrap_or(1) as usize;
                        if anchor_window_alignable(sequence, candidate, 200) {
                            mp_candidate = candidate;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        // Couldn't find a clean mate after 32 tries — likely
                        // the contig is mostly N (e.g. unplaced contig). Skip
                        // this BND rather than emit an unrecoverable truth.
                        continue;
                    }
                    // The BND ALT embeds the reference base at THIS breakend
                    // (VCF 4.2 §5.4). A literal "N" contradicts the REF column,
                    // which carries the real base — see #451.
                    let anchor_char: char = anchor_base.into();
                    (
                        format!("{anchor_char}]{contig_name}:{mp_candidate}]"),
                        Some(contig_name.to_string()),
                        Some(mp_candidate),
                    )
                }
                // sample_variants only emits the four supported types;
                // others were excluded from type_probabilities at fit
                // time. The match is exhaustive for the compiler.
                _ => continue,
            };
            // VCF 4.2 §5.4: a breakend is a POINT. END must equal POS and there is no
            // SVLEN, whatever length the model happened to sample. An anchor reporting
            // END = POS + length claims a multi-kb span for a 1bp event, and tools that
            // size records off END (truvari, bcftools) act on it — while the record's own
            // mate correctly carries END = POS, so one junction would ship two
            // conventions. `fit_from_observations` pins BND length to 1, so a trained
            // model never exposed this; `SvModel` is plain serde JSON and hand-edited or
            // grafted models are not covered by that.
            let sv_end_1based = if sv_type == SvType::Bnd {
                anchor_0based + 1
            } else {
                sv_end_1based
            };
            let mut sv_data = SvData::new(raw_alt, sv_type);
            sv_data.end = Some(sv_end_1based);
            sv_data.svlen = if sv_type != SvType::Bnd {
                Some(match sv_type {
                    SvType::Del => -(length as i64),
                    _ => length as i64,
                })
            } else {
                None
            };
            sv_data.copy_number = copy_number;
            // The chimeric read generator dispatches on these flags, NOT on the ALT
            // string, and they are otherwise only ever set by the input-VCF parser. A de
            // novo BND left them at their false/false defaults, which is case 4 (`]p]t`)
            // — a DIRECT join — while the ALT above says `t]p]`, case 2, head-to-head
            // with the mate piece reverse-complemented. The truth VCF therefore described
            // a different rearrangement than the reads carried, and a direct intra-contig
            // join is just a deletion or duplication: on Delta, Manta placed DEL/DUP
            // candidates within 1-2bp of all 7 planted junctions and emitted no breakend
            // at any of them, so BND recall read 0.000 for a caller that had found every
            // one. Keep these in step with the ALT emitted above.
            if sv_type == SvType::Bnd {
                sv_data.bnd_join_after = true;
                sv_data.bnd_mate_extends_right = false;
            }
            sv_data.mate_contig = m_contig;
            sv_data.mate_pos = m_pos;

            let info_field = build_info_field(sv_type, sv_end_1based, length, copy_number);

            // A BND is a two-sided event. VCF 4.2 §5.4 represents a novel
            // adjacency as a mate PAIR — one record per breakend, cross-linked by
            // MATEID — and that is what callers emit, so a one-sided truth record
            // is unmatchable by any comparison tool (#451). Build stable IDs from
            // the junction's own coordinates so they're deterministic and unique.
            let bnd_pair = if sv_type == SvType::Bnd {
                match (
                    sv_data.mate_pos,
                    sequence.get(sv_data.mate_pos.unwrap_or(0).saturating_sub(1)),
                ) {
                    (Some(mate_1based), Some(&mate_base)) if mate_base != Nucleotide::N => {
                        let anchor_1based = anchor_0based + 1;
                        Some((
                            format!("eidolon_bnd_{contig_name}_{anchor_1based}_{mate_1based}_1"),
                            format!("eidolon_bnd_{contig_name}_{anchor_1based}_{mate_1based}_2"),
                            anchor_1based,
                            mate_1based,
                            mate_base,
                        ))
                    }
                    // Mate base unusable (out of range / N). Emitting a lone
                    // breakend would recreate #451, so drop the whole event.
                    _ => continue,
                }
            } else {
                None
            };

            let (info_1, id_1) = match &bnd_pair {
                Some((id1, id2, ..)) => (format!("{info_field};MATEID={id2}"), Some(id1.clone())),
                None => (info_field, None),
            };

            occupied.push((affected_start, affected_end));
            out.push(Variant {
                variant_type: VariantType::Complex,
                location: anchor_0based,
                reference: vec![anchor_base],
                alternate: AlternateType::Symbolic(sv_data),
                genotype_str: genotype_str.clone(),
                genotype: genotype.clone(),
                allele_fraction: None,
                id: id_1,
                quality_score: None,
                filter: None,
                info: Some(info_1),
                format: vec!["GT".to_string()],
                sample: vec![genotype_str.clone()],
                provenance: Provenance::Denovo,
            });

            // The mate breakend. Reciprocal `t]p]` form, matching the spec's own
            // worked example (`G]17:198982]` <-> `A]2:321681]`) — keeping this
            // orientation means the junction sequence the chimeric read generator
            // builds is unchanged; only the record representation is fixed.
            if let Some((id1, id2, anchor_1based, mate_1based, mate_base)) = bnd_pair {
                let mate_char: char = mate_base.into();
                let mut mate_sv = SvData::new(
                    format!("{mate_char}]{contig_name}:{anchor_1based}]"),
                    SvType::Bnd,
                );
                mate_sv.end = Some(mate_1based);
                mate_sv.svlen = None;
                mate_sv.copy_number = None;
                // Same form as the anchor (`t]p]`), so the same geometry — a reciprocal
                // `t]p]` <-> `t]p]` pair is the spec's own worked example.
                mate_sv.bnd_join_after = true;
                mate_sv.bnd_mate_extends_right = false;
                mate_sv.mate_contig = Some(contig_name.to_string());
                mate_sv.mate_pos = Some(anchor_1based);
                let mate_info = format!(
                    "{};MATEID={id1}",
                    build_info_field(SvType::Bnd, mate_1based, 1, None)
                );
                out.push(Variant {
                    variant_type: VariantType::Complex,
                    location: mate_1based.saturating_sub(1),
                    reference: vec![mate_base],
                    alternate: AlternateType::Symbolic(mate_sv),
                    genotype_str: genotype_str.clone(),
                    genotype,
                    allele_fraction: None,
                    id: Some(id2),
                    quality_score: None,
                    filter: None,
                    info: Some(mate_info),
                    format: vec!["GT".to_string()],
                    sample: vec![genotype_str],
                    provenance: Provenance::Denovo,
                });
            }

            events += 1;
        }

        // Compare EVENTS, not records — a BND emits two records for one event, so
        // out.len() would mask an under-placement.
        if events < n_sv {
            warn!(
                "SvModel: requested {} de novo SV(s) but only placed {} after {} retries \
                 (overlap / N-anchor rejection saturated; cap = {}bp at \
                 max_length_fraction={:.2}). On small contigs, raise sv_max_length_fraction.",
                n_sv, events, retries, upper_len, frac
            );
        } else {
            info!(
                "SvModel: sampled {} de novo SV(s) on a {}bp contig (rate_scale={})",
                out.len(),
                contig_len,
                rate_scale
            );
        }

        out
    }
}

/// Hybrid Poisson sampler:
///   - λ < `LARGE_LAMBDA_THRESHOLD` (30): Knuth's multiplicative algorithm.
///     Exact, cheap, and preserves the RNG-draw sequence the small-λ
///     pipeline tests are seeded against.
///   - λ ≥ `LARGE_LAMBDA_THRESHOLD`: Gaussian approximation N(λ, √λ),
///     sampled via inverse-CDF on a uniform draw. At λ = 30 the
///     approximation error is below 1 part in 1000; by λ = 1000 it's
///     indistinguishable from exact Poisson in practice. This branch
///     exists because human-chromosome-scale λ (~10⁴–10⁵) overflows
///     Knuth's algorithm via `exp(-λ)` underflow around λ ≈ 745.
fn sample_poisson(lambda: f64, rng: &mut NeatRng) -> Result<u64, NeatRngError> {
    if lambda <= 0.0 || !lambda.is_finite() {
        return Ok(0);
    }
    const LARGE_LAMBDA_THRESHOLD: f64 = 30.0;
    if lambda < LARGE_LAMBDA_THRESHOLD {
        // Knuth's multiplicative algorithm. Iterates roughly λ+1 times in
        // expectation, so capped well below the underflow threshold.
        let l = (-lambda).exp();
        let mut k: u64 = 0;
        let mut p: f64 = 1.0;
        loop {
            k += 1;
            p *= rng.random()?;
            if p < l {
                return Ok(k - 1);
            }
            if k > 1_000_000 {
                return Err(NeatRngError::SamplingError(format!(
                    "Poisson sampler exceeded 1e6 iterations at λ={lambda}"
                )));
            }
        }
    } else {
        // Gaussian approximation. Truncation at 0 has probability
        // Φ(−√λ) ≈ 0 for λ ≥ 30 (Φ(−√30) ≈ 3e-7), so the rounding /
        // max(0) clamp doesn't materially bias the mean.
        let normal = NormalDistribution::new(lambda, lambda.sqrt())
            .map_err(|e| NeatRngError::SamplingError(format!("{e:?}")))?;
        let u = rng.random()?.clamp(1e-12, 1.0 - 1e-12);
        let sample = normal
            .sample(u)
            .map_err(|e| NeatRngError::SamplingError(format!("{e:?}")))?;
        if !sample.is_finite() {
            return Err(NeatRngError::SamplingError(format!(
                "Poisson Gaussian-approximation sample non-finite at λ={lambda}: {sample}"
            )));
        }
        Ok(sample.max(0.0).round() as u64)
    }
}

/// Sample `round(exp(N(mu, sigma)))` via inverse-CDF on a uniform draw.
fn sample_log_normal_usize(mu: f64, sigma: f64, rng: &mut NeatRng) -> Result<usize, NeatRngError> {
    if !mu.is_finite() || !sigma.is_finite() || sigma < 0.0 {
        return Err(NeatRngError::SamplingError(format!(
            "log-normal parameters invalid: mu={mu} sigma={sigma}"
        )));
    }
    if sigma == 0.0 {
        return Ok(mu.exp().round().max(0.0) as usize);
    }
    let normal = NormalDistribution::new(mu, sigma)
        .map_err(|e| NeatRngError::SamplingError(format!("{e:?}")))?;
    // statrs's inverse_cdf saturates at the open-interval endpoints, so
    // clamp the uniform draw away from the exact boundaries.
    let u = rng.random()?.clamp(1e-9, 1.0 - 1e-9);
    let normal_sample = normal
        .sample(u)
        .map_err(|e| NeatRngError::SamplingError(format!("{e:?}")))?;
    let v = normal_sample.exp();
    if !v.is_finite() || v <= 0.0 {
        return Err(NeatRngError::SamplingError(format!(
            "log-normal sample non-finite or non-positive: {v}"
        )));
    }
    Ok(v.round().max(0.0) as usize)
}

/// Normalize a non-negative weight vector and produce its cumulative
/// sum, suitable for inverse-CDF picking via a uniform `[0, 1)` draw.
/// Returns an empty `Vec` if all weights are zero.
fn cumulative_normalized(weights: &[f64]) -> Vec<f64> {
    let total: f64 = weights.iter().filter(|w| w.is_finite() && **w >= 0.0).sum();
    if total <= 0.0 || !total.is_finite() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(weights.len());
    let mut acc = 0.0;
    for &w in weights {
        let w = if w.is_finite() && w >= 0.0 { w } else { 0.0 };
        acc += w / total;
        out.push(acc);
    }
    out
}

/// Pick an element from a parallel `values` / `cum_weights` pair using
/// inverse-CDF sampling. `cum_weights` is the normalized cumulative
/// distribution produced by [`cumulative_normalized`].
fn weighted_pick_with<'a, T>(
    values: &'a [T],
    cum_weights: &[f64],
    rng: &mut NeatRng,
) -> Option<&'a T> {
    if values.is_empty() || cum_weights.len() != values.len() {
        return None;
    }
    let u = rng.random().ok()?;
    for (v, c) in values.iter().zip(cum_weights.iter()) {
        if u < *c {
            return Some(v);
        }
    }
    values.last()
}

/// Pick an anchor for a candidate SV of `length` bases on a contig of
/// `contig_len` bases. Returns `(affected_start, affected_end,
/// anchor_0based, sv_end_1based)` or `None` if the contig is too small
/// for the chosen length.
///
/// INVARIANT, for every symbolic SPAN type: `anchor_0based` is the base
/// immediately BEFORE the event and is NOT itself affected. The affected range
/// is always `[anchor_0based + 1, anchor_0based + 1 + length)` 0-based
/// half-open, i.e. 1-based `[POS+1, END]`, so `END - POS == |SVLEN|` holds
/// uniformly (VCF 4.2 symbolic-allele convention).
///
/// `<INS>` and `<BND>` are point events, not spans: nothing in the reference is
/// consumed, POS marks the position, and END == POS.
fn place_sv(
    sv_type: SvType,
    length: usize,
    contig_len: usize,
    rng: &mut NeatRng,
) -> Option<(usize, usize, usize, usize)> {
    if length == 0 || length >= contig_len {
        return None;
    }
    // For DEL we need an anchor BEFORE the deleted bases, so the max
    // anchor is `contig_len - length - 1`. For DUP/CNV, the anchor IS
    // the first affected base; max anchor is `contig_len - length`.
    let max_anchor_inclusive = match sv_type {
        SvType::Del => contig_len.checked_sub(length + 1)?,
        _ => contig_len.checked_sub(length)?,
    };
    if max_anchor_inclusive == 0 {
        return None;
    }
    let anchor = rng.range_i64(0, max_anchor_inclusive as i64).ok()?.max(0) as usize;
    let (start, end, vcf_anchor_0based, sv_end_1based) = match sv_type {
        SvType::Del => {
            let s = anchor + 1;
            let e = s + length;
            // 1-based END for SvData: span includes the anchor + deleted
            // bases, so END_1based = anchor_1based + length = anchor + 1 + length.
            (s, e, anchor, anchor + 1 + length)
        }
        SvType::Ins => {
            // Insertions are point events in the reference: POS carries the base
            // before the inserted sequence and END == POS.
            let s = anchor;
            let e = s + 1;
            (s, e, anchor, anchor + 1)
        }
        SvType::Bnd => {
            // A breakend is a point and POS IS the breakend base — the
            // anchor-before convention does not apply, and END is forced to POS
            // by the caller (VCF 4.2 5.4).
            let s = anchor;
            let e = s + 1;
            (s, e, anchor, anchor + 1)
        }
        _ => {
            // DUP/CNV/INV. `anchor` is the FIRST AFFECTED base, but VCF 4.2 puts
            // POS on the base BEFORE a symbolic event, so the record's anchor is
            // `anchor - 1` while the affected range stays [anchor, anchor+length).
            // Emitting `anchor` directly is what made POS mean "last unaffected"
            // for DEL and "first affected" for everything else, so END-POS equalled
            // |SVLEN| for DEL and |SVLEN|-1 for the rest in the same file. It also
            // put the first duplicated base in REF, where the preceding base belongs.
            // There is no base before position 0.
            if anchor == 0 {
                return None;
            }
            let s = anchor;
            let e = s + length;
            // 1-based END = last affected base = anchor + length.
            (s, e, anchor - 1, anchor + length)
        }
    };
    Some((start, end, vcf_anchor_0based, sv_end_1based))
}

/// Compute the 0-based half-open affected range of an existing symbolic
/// SV, for the overlap-rejection check. Returns `None` if the record
/// isn't symbolic or has no usable span.
fn affected_range_for_existing(v: &Variant, block_end: usize) -> Option<(usize, usize)> {
    let sv = v.alternate.as_symbolic()?;
    let span = match sv.sv_type {
        SvType::Ins | SvType::Bnd => Some(1),
        _ => sv.span(v.location.saturating_add(1)),
    }?;
    if span == 0 {
        return None;
    }
    let raw_end = v.location.saturating_add(span).min(block_end);
    let start = match sv.sv_type {
        // Point events consume no reference bases, so the anchor itself marks the
        // position for overlap purposes — there is no POS+1..END to speak of.
        SvType::Ins | SvType::Bnd => v.location.min(block_end),
        // Span events: POS is the base BEFORE the event (VCF 4.2), so the affected
        // range starts at POS+1. This used to be DEL-only, which is what let a DUP
        // and a DEL in the same file mean different things by POS.
        _ => v.location.saturating_add(1).min(block_end),
    };
    if start >= raw_end {
        return None;
    }
    Some((start, raw_end))
}

fn ranges_overlap(a_start: usize, a_end: usize, b_start: usize, b_end: usize) -> bool {
    a_start < b_end && b_start < a_end
}

/// True iff the ±`window` bp window around `pos` is at least 80% non-N.
/// Used by `sample_variants` to reject SV positions whose chimeric
/// junction pieces would be mostly N and therefore unalignable by BWA.
///
/// Before this gate (added for #224), gen-reads would emit truth records
/// at positions near reference N-tracts. The single-base anchor check at
/// line 473 only rejects positions where the literal anchor is N — it
/// passes positions where the anchor happens to be a clean base
/// adjacent to a long N-tract. The chimeric pair generation then samples
/// a `read_len`-sized window extending into the N-tract, producing reads
/// that BWA can't place uniquely → 0% caller recall on the resulting
/// truth records (see chr22 v1.13.0 validation: 13 of 35 somatic BND
/// truth records were unrecoverable for this reason).
///
/// The 80% threshold is loose on purpose: BWA-MEM tolerates some N
/// bases, and a few percent N is fine. The pathological case is windows
/// that are 100% N (centromere/telomere bulk) or close to it.
pub(crate) fn anchor_window_alignable(sequence: &[Nucleotide], pos: usize, window: usize) -> bool {
    let start = pos.saturating_sub(window);
    let end = pos.saturating_add(window).min(sequence.len());
    if end <= start {
        return false;
    }
    let total = end - start;
    let n_count = sequence[start..end]
        .iter()
        .filter(|&&b| b == Nucleotide::N)
        .count();
    (n_count * 100) <= (total * 20)
}

/// Draw `length` random nucleotides for a novel insertion, weighted by the
/// single-base composition of a 501 bp window centered on `anchor_0based`
/// (±250 bp, clamped to contig bounds; `N`s in the window are ignored).
/// Falls back to uniform ACGT if the window has no usable bases (e.g.
/// an all-N region). #190.
fn sample_novel_insertion_bases(
    sequence: &[Nucleotide],
    anchor_0based: usize,
    length: usize,
    rng: &mut NeatRng,
) -> Result<Vec<Nucleotide>, NeatRngError> {
    const WINDOW: usize = 250;
    let start = anchor_0based.saturating_sub(WINDOW);
    let end = anchor_0based.saturating_add(WINDOW).min(sequence.len());

    let mut counts = [0u64; 4];
    for &b in &sequence[start..end] {
        match b {
            Nucleotide::A => counts[0] += 1,
            Nucleotide::C => counts[1] += 1,
            Nucleotide::G => counts[2] += 1,
            Nucleotide::T => counts[3] += 1,
            _ => {} // N or other ambiguity codes — skip
        }
    }

    let total: u64 = counts.iter().sum();
    // Cumulative ACGT-ordered probabilities for the inverse-CDF sample.
    let cum: [f64; 4] = if total == 0 {
        [0.25, 0.5, 0.75, 1.0]
    } else {
        let t = total as f64;
        let p_a = counts[0] as f64 / t;
        let p_c = counts[1] as f64 / t;
        let p_g = counts[2] as f64 / t;
        [p_a, p_a + p_c, p_a + p_c + p_g, 1.0]
    };

    let bases = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
    let mut out = Vec::with_capacity(length);
    for _ in 0..length {
        let r = rng.random()?;
        let idx = cum.iter().position(|&c| r < c).unwrap_or(3);
        out.push(bases[idx]);
    }
    Ok(out)
}

/// Emit a single-line VCF `INFO` field carrying the structural-variant
/// fields downstream readers expect, so de novo records round-trip
/// through the writer with the same shape as user-supplied input.
fn build_info_field(
    sv_type: SvType,
    end_1based: usize,
    length: usize,
    copy_number: Option<u32>,
) -> String {
    let mut parts: Vec<String> = Vec::new();
    let svtype = match sv_type {
        SvType::Del => "DEL",
        SvType::Dup => "DUP",
        SvType::Cnv => "CNV",
        SvType::Ins => "INS",
        SvType::Inv => "INV",
        SvType::Bnd => "BND",
        SvType::Unknown => "UNKNOWN",
    };
    parts.push(format!("SVTYPE={svtype}"));
    parts.push(format!("END={end_1based}"));
    if sv_type != SvType::Bnd {
        let svlen = match sv_type {
            SvType::Del => -(length as i64),
            _ => length as i64,
        };
        parts.push(format!("SVLEN={svlen}"));
    }
    if let Some(cn) = copy_number {
        parts.push(format!("CN={cn}"));
    }
    parts.join(";")
}

/// Format a VCF genotype string for `ploidy` haplotypes. Homozygous →
/// all `1`s; heterozygous → one `1` and the rest `0`. Phased on `|`
/// would be more standard, but the existing v1.9.0 writers use `/` for
/// both unphased and phased records, so match that.
fn genotype_string(genotype: &Genotype, ploidy: usize) -> String {
    let n = ploidy.max(1);
    let parts: Vec<&str> = match genotype {
        Genotype::Homozygous => (0..n).map(|_| "1").collect(),
        Genotype::Heterozygous => {
            // One alt + (n-1) refs.
            let mut v = vec!["0"; n.saturating_sub(1)];
            v.push("1");
            v
        }
    };
    parts.join("/")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide;
    use crate::structs::variants::{AlternateType, SvData, VariantType, parse_bnd_alt};

    #[test]
    fn default_sv_model_is_not_usable() {
        let m = SvModel::default();
        assert!(!m.is_usable());
        assert_eq!(m.per_base_rate, 0.0);
        assert!(m.type_probabilities.is_empty());
        assert!(m.length_log_normal.is_empty());
        assert!(m.cnv_copy_number_distribution.is_empty());
        assert_eq!(m.homozygous_frequency, 0.0);
    }

    #[test]
    fn min_sv_length_is_the_industry_50bp_cutoff() {
        // Lock in the constant — downstream samplers and the trainer
        // both key off it; a silent change would shift the entire
        // SV/indel boundary across the pipeline.
        assert_eq!(MIN_SV_LENGTH_BP, 50);
    }

    #[test]
    fn populated_sv_model_round_trips_through_serde_json() {
        let mut m = SvModel {
            per_base_rate: 2.3e-7,
            homozygous_frequency: 0.31,
            ..Default::default()
        };
        m.type_probabilities.insert(SvType::Del, 0.6);
        m.type_probabilities.insert(SvType::Dup, 0.3);
        m.type_probabilities.insert(SvType::Cnv, 0.1);
        m.length_log_normal.insert(SvType::Del, (7.5, 1.8));
        m.length_log_normal.insert(SvType::Dup, (8.1, 1.5));
        m.cnv_copy_number_distribution.insert(0, 0.05);
        m.cnv_copy_number_distribution.insert(1, 0.45);
        m.cnv_copy_number_distribution.insert(3, 0.30);
        m.cnv_copy_number_distribution.insert(4, 0.20);

        let json = serde_json::to_string(&m).unwrap();
        let back: SvModel = serde_json::from_str(&json).unwrap();
        assert_eq!(back, m);
        assert!(back.is_usable());
    }

    #[test]
    fn is_usable_requires_both_rate_and_type_probabilities() {
        let mut only_rate = SvModel::default();
        only_rate.per_base_rate = 1e-6;
        assert!(!only_rate.is_usable());

        let mut only_types = SvModel::default();
        only_types.type_probabilities.insert(SvType::Del, 1.0);
        assert!(!only_types.is_usable());
    }

    /// Build a synthetic symbolic-SV `Variant` with a 1-based POS,
    /// a 1-based-inclusive END (or `None` to skip), the given type,
    /// genotype, and optional copy_number. Returns a value owned by
    /// the caller — tests construct a `Vec<Variant>` and pass a slice
    /// to the fitter.
    fn sv(
        pos_1based: usize,
        end_1based: Option<usize>,
        sv_type: SvType,
        genotype: Genotype,
        copy_number: Option<u32>,
    ) -> Variant {
        let raw_alt = match sv_type {
            SvType::Del => "<DEL>",
            SvType::Dup => "<DUP>",
            SvType::Cnv => "<CNV>",
            SvType::Ins => "<INS>",
            SvType::Inv => "<INV>",
            SvType::Bnd => "B[chr1:1000[",
            SvType::Unknown => "<NOVEL>",
        };
        let mut sv_data = SvData::new(raw_alt, sv_type);
        sv_data.end = end_1based;
        sv_data.copy_number = copy_number;
        let gt_str = match genotype {
            Genotype::Homozygous => "1/1",
            Genotype::Heterozygous => "0/1",
        };
        Variant {
            variant_type: VariantType::Complex,
            location: pos_1based,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(sv_data),
            genotype_str: gt_str.to_string(),
            genotype,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        }
    }

    #[test]
    fn fit_returns_none_for_empty_observations() {
        let m = SvModel::fit_from_observations(&[], 10_000);
        assert!(m.is_none());
    }

    #[test]
    fn fit_returns_none_when_total_reflen_zero() {
        // Per-base rate would be 0/0 — `None` keeps the contract that a
        // `Some(SvModel)` is always usable.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(700), SvType::Del, Genotype::Heterozygous, None),
        ];
        assert!(SvModel::fit_from_observations(&obs, 0).is_none());
    }

    #[test]
    fn fit_drops_sub_50bp_records_and_unsupported_types() {
        // 5 input records, only the two ≥50bp DELs and the two BNDs should survive into
        // the fit. The INV is dropped as unsupported; the
        // 30bp DEL is dropped as an indel-disguised-as-SV.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None), // 101bp DEL — keep
            sv(500, Some(700), SvType::Del, Genotype::Heterozygous, None), // 201bp DEL — keep
            sv(800, Some(829), SvType::Del, Genotype::Heterozygous, None), // 30bp — drop (sub-50)
            sv(1000, Some(2000), SvType::Inv, Genotype::Homozygous, None), // INV — drop
            sv(3000, None, SvType::Bnd, Genotype::Heterozygous, None),   // BND — keep
            sv(4000, None, SvType::Bnd, Genotype::Heterozygous, None), // BND — keep (need 2 for fit)
        ];
        let m = SvModel::fit_from_observations(&obs, 100_000).unwrap();
        assert!(m.is_usable());
        assert_eq!(m.type_probabilities.len(), 2);
        assert_eq!(m.type_probabilities.get(&SvType::Del), Some(&0.5));
        assert_eq!(m.type_probabilities.get(&SvType::Bnd), Some(&0.5));
        // 4 surviving SVs / 100kb = 4e-5 per base.
        assert!((m.per_base_rate - 4e-5).abs() < 1e-12);
        // 1 homozygous of 4 surviving → 0.25
        assert!((m.homozygous_frequency - 0.25).abs() < 1e-12);
    }

    #[test]
    fn fit_drops_types_with_one_observation() {
        // DEL has 2 observations → kept. DUP has 1 → dropped (can't fit σ
        // from a single sample). Probabilities renormalize over what's
        // left, so DEL becomes 1.0.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(700), SvType::Del, Genotype::Homozygous, None),
            sv(1000, Some(1500), SvType::Dup, Genotype::Heterozygous, None),
        ];
        let m = SvModel::fit_from_observations(&obs, 10_000).unwrap();
        assert_eq!(m.type_probabilities.len(), 1);
        assert!(m.type_probabilities.contains_key(&SvType::Del));
        assert!(!m.type_probabilities.contains_key(&SvType::Dup));
        assert!(m.length_log_normal.contains_key(&SvType::Del));
        assert!(!m.length_log_normal.contains_key(&SvType::Dup));
    }

    #[test]
    fn fit_recovers_log_normal_parameters_within_tolerance() {
        // For 4 deletions of length e^7, e^7, e^8, e^8 (i.e.
        // ln-lengths = 7, 7, 8, 8), MLE gives μ=7.5, σ²=0.25, σ=0.5.
        let lens = [
            (7.0_f64.exp() as usize), // ~1097
            (7.0_f64.exp() as usize),
            (8.0_f64.exp() as usize), // ~2981
            (8.0_f64.exp() as usize),
        ];
        let mut pos = 100usize;
        let mut obs = Vec::new();
        for &len in &lens {
            obs.push(sv(
                pos,
                Some(pos + len - 1),
                SvType::Del,
                Genotype::Heterozygous,
                None,
            ));
            pos += len + 1000; // separate them in coords (doesn't affect fit)
        }
        let m = SvModel::fit_from_observations(&obs, 1_000_000).unwrap();
        let (mu, sigma) = m.length_log_normal[&SvType::Del];
        // Rounding through usize introduces tiny error vs the exact
        // analytic value, but it should easily land within 5%.
        assert!((mu - 7.5).abs() < 0.05, "mu was {mu}, expected ~7.5");
        assert!(
            (sigma - 0.5).abs() < 0.05,
            "sigma was {sigma}, expected ~0.5"
        );
    }

    #[test]
    fn fit_populates_cnv_copy_number_distribution() {
        // Three CNVs: CN=1, CN=1, CN=3 → P(1)=2/3, P(3)=1/3.
        // Two CNVs without CN — they're kept but don't contribute to
        // the CN distribution (cn_counts is empty for those).
        let obs = vec![
            sv(100, Some(500), SvType::Cnv, Genotype::Heterozygous, Some(1)),
            sv(
                1000,
                Some(1500),
                SvType::Cnv,
                Genotype::Heterozygous,
                Some(1),
            ),
            sv(
                2000,
                Some(2500),
                SvType::Cnv,
                Genotype::Heterozygous,
                Some(3),
            ),
            sv(3000, Some(3500), SvType::Cnv, Genotype::Heterozygous, None),
            sv(4000, Some(4500), SvType::Cnv, Genotype::Heterozygous, None),
        ];
        let m = SvModel::fit_from_observations(&obs, 100_000).unwrap();
        let cn = &m.cnv_copy_number_distribution;
        assert_eq!(cn.len(), 2);
        assert!((cn[&1] - 2.0 / 3.0).abs() < 1e-12);
        assert!((cn[&3] - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn fit_returns_none_when_every_type_thin() {
        // Two records, two different types — neither type clears the
        // 2-observation bar, so the fit returns None.
        let obs = vec![
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(1500), SvType::Dup, Genotype::Heterozygous, None),
        ];
        assert!(SvModel::fit_from_observations(&obs, 10_000).is_none());
    }

    #[test]
    fn fit_ignores_non_symbolic_alts() {
        // A literal SNP slipped in shouldn't blow up the fitter — it's
        // not a symbolic record, so it's silently skipped.
        use crate::structs::variants::Variant as V;
        let snp = V::new(
            VariantType::SNP,
            42,
            &vec![Nucleotide::A],
            &vec![Nucleotide::C],
            &mut vec![1, 0],
        )
        .unwrap();
        let obs = vec![
            snp,
            sv(100, Some(200), SvType::Del, Genotype::Homozygous, None),
            sv(500, Some(700), SvType::Del, Genotype::Heterozygous, None),
        ];
        let m = SvModel::fit_from_observations(&obs, 10_000).unwrap();
        assert_eq!(m.type_probabilities.len(), 1);
        assert!(m.is_usable());
    }

    // ── sample_variants tests ────────────────────────────────────────

    fn deterministic_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec!["phase3".to_string(), "sample_variants".to_string()]).unwrap()
    }

    /// Build a usable `SvModel` with a single SV type so we can assert
    /// sampling produces only that type. Length distribution centered
    /// around `e^7` ≈ 1097bp.
    fn model_with_single_type(sv_type: SvType, per_base_rate: f64) -> SvModel {
        let mut m = SvModel::default();
        m.per_base_rate = per_base_rate;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.insert(sv_type, 1.0);
        let len_dist = if sv_type == SvType::Bnd {
            (0.0, 0.0)
        } else {
            (7.0, 0.4)
        };
        m.length_log_normal.insert(sv_type, len_dist);
        if sv_type == SvType::Cnv {
            m.cnv_copy_number_distribution.insert(0, 0.3);
            m.cnv_copy_number_distribution.insert(3, 0.5);
            m.cnv_copy_number_distribution.insert(4, 0.2);
        }
        m
    }

    /// An ACGT-only contig of the requested length. No `N`s, so anchor
    /// rejection never fires.
    fn acgt_sequence(len: usize) -> Vec<Nucleotide> {
        let bases = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        (0..len).map(|i| bases[i % 4]).collect()
    }

    #[test]
    fn sample_variants_returns_empty_when_rate_scale_zero() {
        let m = model_with_single_type(SvType::Del, 1e-3);
        let seq = acgt_sequence(50_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 50_000, &[], &seq, 2, 0.0, 0.25, &mut rng);
        assert!(
            out.is_empty(),
            "rate_scale=0 must disable generation; got {} SV(s)",
            out.len()
        );
    }

    #[test]
    fn sample_variants_returns_empty_when_model_unusable() {
        let m = SvModel::default();
        let seq = acgt_sequence(50_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 50_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(out.is_empty(), "unusable model must yield no SVs");
    }

    #[test]
    fn sample_variants_returns_empty_when_contig_too_small() {
        // Anything below 2 × MIN_SV_LENGTH_BP is rejected outright —
        // we can't reliably place a ≥50bp SV in a 50bp contig.
        let m = model_with_single_type(SvType::Del, 1e-3);
        let seq = acgt_sequence(MIN_SV_LENGTH_BP);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", MIN_SV_LENGTH_BP, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(out.is_empty());
    }

    /// A BND is a two-sided event: VCF 4.2 §5.4 wants a mate PAIR, cross-linked by
    /// MATEID, with each record's ALT embedding the reference base at *its own*
    /// breakend. Before #451 the sampler emitted a single record with `id: None`,
    /// no MATEID, and a hardcoded `N` in the ALT — unmatchable by any comparison
    /// tool, and `sample_variants_emits_all_supported_types` passed anyway because
    /// it only asserts the BND *type appears*.
    #[test]
    fn denovo_bnd_emits_a_cross_linked_mate_pair() {
        let mut m = SvModel::default();
        m.per_base_rate = 5e-3;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.clear();
        m.type_probabilities.insert(SvType::Bnd, 1.0);
        m.length_log_normal.insert(SvType::Bnd, (7.0, 0.3));
        let seq = acgt_sequence(50_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 50_000, &[], &seq, 2, 1.0, 0.25, &mut rng);

        let bnds: Vec<&Variant> = out
            .iter()
            .filter(|v| {
                v.alternate
                    .as_symbolic()
                    .map(|s| s.sv_type == SvType::Bnd)
                    .unwrap_or(false)
            })
            .collect();
        assert!(!bnds.is_empty(), "expected at least one de novo BND");
        assert_eq!(
            bnds.len() % 2,
            0,
            "BND records must come in pairs: {bnds:?}"
        );

        // Every record has an ID, and MATEID resolves to a record that points back.
        let by_id: std::collections::HashMap<&str, &Variant> = bnds
            .iter()
            .map(|v| (v.id.as_deref().expect("BND must carry an ID"), *v))
            .collect();
        assert_eq!(by_id.len(), bnds.len(), "BND IDs must be unique");

        for v in &bnds {
            let info = v.info.as_deref().unwrap_or("");
            let mateid = info
                .split(';')
                .find_map(|f| f.strip_prefix("MATEID="))
                .unwrap_or_else(|| panic!("BND record has no MATEID: {info}"));
            let mate = by_id
                .get(mateid)
                .unwrap_or_else(|| panic!("MATEID={mateid} does not resolve to a record"));
            let back = mate
                .info
                .as_deref()
                .unwrap_or("")
                .split(';')
                .find_map(|f| f.strip_prefix("MATEID="))
                .unwrap_or("");
            assert_eq!(
                back,
                v.id.as_deref().unwrap(),
                "MATEID linkage is not reciprocal"
            );

            // The ALT's embedded breakend base must equal this record's REF.
            let alt = v.alternate.as_symbolic().unwrap().raw_alt.clone();
            let ref_char: char = v.reference[0].into();
            assert_eq!(
                alt.chars().next().unwrap(),
                ref_char,
                "BND ALT must embed its own REF base (not a literal N): REF={ref_char} ALT={alt}"
            );

            // And the mate's position must be the mate record's own location.
            let mate_pos_1based = v.alternate.as_symbolic().unwrap().mate_pos.unwrap();
            assert_eq!(
                mate_pos_1based.saturating_sub(1),
                mate.location,
                "mate_pos must point at the mate record's location"
            );
        }
    }

    /// The chimeric read generator dispatches on `bnd_join_after` /
    /// `bnd_mate_extends_right` — NOT on the ALT string. Those flags are assigned only
    /// by the input-VCF parser (`variants.rs`), so a de novo BND leaves them at their
    /// `SvData::new` defaults of false/false, which is case 4 (`]p]t`): a DIRECT join
    /// with no reverse complement. Meanwhile the ALT written to the truth VCF is
    /// `t]p]` — case 2, head-to-head, reverse-complemented.
    ///
    /// So the truth VCF describes a different rearrangement than the reads encode. On
    /// Delta this showed up as Manta placing DEL/DUP candidates within 1-2bp of all 7
    /// planted junctions while emitting no breakend there at all, and BND recall being
    /// reported as 0.000 for a caller that had in fact found every junction.
    #[test]
    fn denovo_bnd_read_geometry_matches_the_alt_it_declares() {
        let mut m = SvModel::default();
        m.per_base_rate = 5e-3;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.clear();
        m.type_probabilities.insert(SvType::Bnd, 1.0);
        m.length_log_normal.insert(SvType::Bnd, (7.0, 0.3));
        let seq = acgt_sequence(50_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 50_000, &[], &seq, 2, 1.0, 0.25, &mut rng);

        let mut checked = 0;
        for v in &out {
            let Some(sv) = v.alternate.as_symbolic() else {
                continue;
            };
            if sv.sv_type != SvType::Bnd {
                continue;
            }
            // The same parser the input-VCF path uses, so the two cannot drift apart.
            let (_mc, _mp, j_after, m_right) = parse_bnd_alt(&sv.raw_alt);
            assert_eq!(
                sv.bnd_join_after, j_after,
                "ALT {} declares join_after={j_after}, but the read generator will use \
                 bnd_join_after={} — the reads would encode a different junction than \
                 the truth VCF describes",
                sv.raw_alt, sv.bnd_join_after
            );
            assert_eq!(
                sv.bnd_mate_extends_right, m_right,
                "ALT {} declares mate_extends_right={m_right}, but the read generator \
                 will use bnd_mate_extends_right={}",
                sv.raw_alt, sv.bnd_mate_extends_right
            );
            checked += 1;
        }
        assert!(
            checked > 0,
            "no de novo BND was produced — test proves nothing"
        );
    }

    /// The sampling loop is gated on placed *events*, not emitted records — a BND
    /// emits two records, so gating on `out.len()` would halve the SV count.
    #[test]
    fn bnd_pair_counts_as_one_event_toward_n_sv() {
        let mut m = SvModel::default();
        m.per_base_rate = 5e-3;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.clear();
        m.type_probabilities.insert(SvType::Bnd, 1.0);
        m.length_log_normal.insert(SvType::Bnd, (7.0, 0.3));
        let seq = acgt_sequence(50_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 50_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        // Junctions, not records.
        let junctions = out.len() / 2;
        assert!(
            junctions >= 2,
            "an all-BND model should place multiple junctions, got {junctions} \
             ({} records) — a loop gated on out.len() would halve this",
            out.len()
        );
    }

    #[test]
    fn sample_variants_emits_all_supported_types() {
        // The sampler must emit DEL/DUP/CNV/INS/BND/INV if they are in the model.
        let mut m = SvModel::default();
        m.per_base_rate = 5e-2;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.insert(SvType::Del, 0.16);
        m.type_probabilities.insert(SvType::Dup, 0.16);
        m.type_probabilities.insert(SvType::Cnv, 0.17);
        m.type_probabilities.insert(SvType::Ins, 0.17);
        m.type_probabilities.insert(SvType::Bnd, 0.17);
        m.type_probabilities.insert(SvType::Inv, 0.17);
        for t in [
            SvType::Del,
            SvType::Dup,
            SvType::Cnv,
            SvType::Ins,
            SvType::Bnd,
            SvType::Inv,
        ] {
            m.length_log_normal.insert(t, (7.0, 0.3));
        }
        let seq = acgt_sequence(200_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 200_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        let mut seen_types = std::collections::HashSet::new();
        let mut seen_literal_ins = false;
        for v in &out {
            // Post-#190: de novo INS records are LITERAL Insertion variants
            // (REF=anchor, ALT=anchor+novel bases) rather than symbolic
            // `<INS>`, so they live in `alternate.as_literal()` and have
            // `variant_type == Insertion`. Every other SV type is still
            // symbolic and goes via `alternate.as_symbolic()`.
            match v.alternate.as_symbolic() {
                Some(sv) => {
                    seen_types.insert(sv.sv_type);
                }
                None => {
                    assert_eq!(
                        v.variant_type,
                        VariantType::Insertion,
                        "non-symbolic ALT must be a literal INS; got {:?}",
                        v.variant_type
                    );
                    seen_literal_ins = true;
                    seen_types.insert(SvType::Ins);
                }
            }
        }
        // At this rate and probability we expect to see all six types,
        // with INS landing on the literal path.
        assert!(seen_types.contains(&SvType::Del));
        assert!(seen_types.contains(&SvType::Dup));
        assert!(seen_types.contains(&SvType::Cnv));
        assert!(seen_types.contains(&SvType::Ins));
        assert!(seen_types.contains(&SvType::Bnd));
        assert!(seen_types.contains(&SvType::Inv));
        assert!(seen_literal_ins, "expected at least one literal-ALT INS");
    }

    #[test]
    fn sample_variants_length_cap_is_configurable() {
        // #229: the SV length cap is `contig_len * max_length_fraction`.
        // Build a DEL-only model whose lengths cluster tightly around 3000bp
        // (well above the default 0.25 cap on an 8kb contig, which is 2000bp,
        // but comfortably below the 1.0 cap of 8000bp). With the default
        // fraction every DEL exceeds the cap and is rejected, so the sampler
        // starves and emits nothing. Raising the fraction to 1.0 lets the same
        // draws through. This is the small-contig / aneuploidy case from #229.
        let mut m = SvModel::default();
        m.per_base_rate = 2e-3;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.insert(SvType::Del, 1.0);
        // e^8.006 ≈ 3000bp, sigma 0.08 → essentially all draws in [2400, 3750].
        m.length_log_normal.insert(SvType::Del, (8.006, 0.08));

        let contig_len = 8000;
        let seq = acgt_sequence(contig_len);

        // Default cap (0.25 → 2000bp): every DEL is too long → starved.
        let mut rng = deterministic_rng();
        let capped = m.sample_variants("chr1", contig_len, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(
            capped.is_empty(),
            "with the default 0.25 cap all DELs exceed contig_len/4 and are rejected; \
             got {} variants",
            capped.len()
        );

        // Full-length cap (1.0 → 8000bp): the same draws now fit and are emitted.
        let mut rng = deterministic_rng();
        let uncapped = m.sample_variants("chr1", contig_len, &[], &seq, 2, 1.0, 1.0, &mut rng);
        assert!(
            !uncapped.is_empty(),
            "raising sv_max_length_fraction to 1.0 should let the ~3kb DELs through"
        );
        for v in &uncapped {
            assert_eq!(
                v.alternate.as_symbolic().expect("symbolic SV").sv_type,
                SvType::Del
            );
        }
    }

    #[test]
    fn sampled_dels_have_correct_anchor_and_end_invariant() {
        // For each generated DEL: variant.location is the 0-based anchor
        // (just before the deletion), sv.end is 1-based inclusive
        // covering the last deleted base. The span computed from those
        // two fields must satisfy:
        //   span_bases = end - (location + 1) + 1 = end - location
        // …and the affected bases ≥ MIN_SV_LENGTH_BP.
        let m = model_with_single_type(SvType::Del, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(
            !out.is_empty(),
            "expected at least one sampled DEL at this rate"
        );
        for v in &out {
            let sv = v.alternate.as_symbolic().expect("must be symbolic");
            assert_eq!(sv.sv_type, SvType::Del);
            let end = sv.end.expect("END must be populated");
            // location is 0-based; end is 1-based inclusive. Deleted
            // bases are 0-based [location+1, end), so length = end - location - 1.
            let length = end - v.location - 1;
            assert!(
                length >= MIN_SV_LENGTH_BP,
                "DEL length {} below MIN_SV_LENGTH_BP",
                length
            );
            // Fits within contig.
            assert!(end <= 100_000);
        }
    }

    #[test]
    fn sampled_dups_have_correct_pos_and_end_invariant() {
        let m = model_with_single_type(SvType::Dup, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(!out.is_empty());
        for v in &out {
            let sv = v.alternate.as_symbolic().expect("must be symbolic");
            assert_eq!(sv.sv_type, SvType::Dup);
            let end = sv.end.expect("END must be populated");
            // For DUP, the duplicated region is 0-based [location, end).
            let length = end - v.location;
            assert!(length >= MIN_SV_LENGTH_BP);
            assert!(end <= 100_000);
        }
    }

    #[test]
    fn sampled_invs_have_correct_pos_and_end_invariant() {
        let m = model_with_single_type(SvType::Inv, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(!out.is_empty());
        for v in &out {
            let sv = v.alternate.as_symbolic().expect("must be symbolic");
            assert_eq!(sv.sv_type, SvType::Inv);
            let end = sv.end.expect("END must be populated");
            // For INV, the inverted region is 0-based [location, end).
            let length = end - v.location;
            assert!(length >= MIN_SV_LENGTH_BP);
            assert!(end <= 100_000);
        }
    }

    /// Post-#190: de novo `<INS>` records are LITERAL Insertion variants,
    /// not symbolic `<INS>`. Their length lives in the REF/ALT length
    /// delta, the insertion site is `variant.location`, and the ALT
    /// content is the actual novel sequence plus the anchor base.
    /// (The richer in-memory shape test — uniqueness of fields, ACGT-only
    /// novel bases, MutatedMap routing to variant_map — lives in
    /// `de_novo_ins_emits_literal_variant_with_novel_sequence`. This one
    /// just pins the length invariant: ALT.len() - REF.len() ≥ MIN_SV.)
    #[test]
    fn sampled_ins_emit_literal_with_min_length_inserted_sequence() {
        let m = model_with_single_type(SvType::Ins, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(!out.is_empty());
        for v in &out {
            assert_eq!(v.variant_type, VariantType::Insertion);
            let alt = v
                .alternate
                .as_literal()
                .expect("de novo INS must be literal-ALT");
            assert_eq!(v.reference.len(), 1, "REF should be a single anchor base");
            let inserted_len = alt.len().saturating_sub(v.reference.len());
            assert!(
                inserted_len >= MIN_SV_LENGTH_BP,
                "inserted-sequence length {inserted_len} must be ≥ {MIN_SV_LENGTH_BP}"
            );
        }
    }

    #[test]
    fn sampled_cnvs_carry_info_cn_from_distribution() {
        let m = model_with_single_type(SvType::Cnv, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(!out.is_empty());
        let expected_cns: std::collections::HashSet<u32> =
            m.cnv_copy_number_distribution.keys().copied().collect();
        for v in &out {
            let sv = v.alternate.as_symbolic().expect("must be symbolic");
            assert_eq!(sv.sv_type, SvType::Cnv);
            let cn = sv.copy_number.expect("<CNV> must carry CN");
            assert!(
                expected_cns.contains(&cn),
                "sampled CN={cn} not in training distribution {:?}",
                expected_cns
            );
        }
    }

    #[test]
    fn sampled_cnvs_use_fallback_cn_when_trained_dist_is_empty() {
        // gnomAD-SV trains with Cnv in type_probabilities but
        // `cnv_copy_number_distribution: {}` because the per-sample CN
        // lives in FORMAT/SAMPLE rather than INFO. The sampler must
        // still emit `<CNV>` records that carry a sensible
        // `copy_number` — drawing from `FALLBACK_CN_DISTRIBUTION` —
        // otherwise gen-reads warn-and-skips every one of them.
        let mut m = SvModel::default();
        m.per_base_rate = 5e-3;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.insert(SvType::Cnv, 1.0);
        m.length_log_normal.insert(SvType::Cnv, (7.0, 0.4));
        // Intentionally leave cnv_copy_number_distribution empty.
        assert!(m.cnv_copy_number_distribution.is_empty());

        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(!out.is_empty(), "expected at least one sampled CNV");

        let fallback_cns: std::collections::HashSet<u32> =
            FALLBACK_CN_DISTRIBUTION.iter().map(|(c, _)| *c).collect();
        for v in &out {
            let sv = v.alternate.as_symbolic().expect("must be symbolic");
            assert_eq!(sv.sv_type, SvType::Cnv);
            let cn = sv
                .copy_number
                .expect("CNV must carry CN even with empty trained dist");
            assert!(
                fallback_cns.contains(&cn),
                "sampled CN={cn} not in FALLBACK_CN_DISTRIBUTION {:?}",
                fallback_cns
            );
        }
    }

    #[test]
    fn sampled_bnds_have_fixed_1bp_length() {
        // BNDs currently use a nominal fixed length of 1bp.
        let m = model_with_single_type(SvType::Bnd, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        assert!(!out.is_empty());
        for v in &out {
            let sv = v.alternate.as_symbolic().expect("must be symbolic");
            assert_eq!(sv.sv_type, SvType::Bnd);
            // Parse the ALT we emitted and compare the WHOLE tuple against the struct
            // the read generator will consume. `contains("chr1:")` and
            // `mate_pos.is_some()` are existence checks: shifting the ALT's mate
            // coordinate by 137bp — so the truth text points somewhere the reads do not
            // — satisfied both of them.
            assert_eq!(
                parse_bnd_alt(&sv.raw_alt),
                (
                    Some("chr1".to_string()),
                    sv.mate_pos,
                    sv.bnd_join_after,
                    sv.bnd_mate_extends_right
                ),
                "ALT {} disagrees with the fields the read generator uses",
                sv.raw_alt
            );
            assert_eq!(sv.mate_contig.as_deref(), Some("chr1"));
            let end = sv.end.expect("END must be populated");
            // For BND, length = 1, so end = location + 1.
            assert_eq!(end, v.location + 1);
        }
    }

    /// Sample BNDs with a NON-degenerate length distribution. `model_with_single_type`
    /// special-cases BND to `(0.0, 0.0)`, so every existing BND test runs at length 1 and
    /// cannot see a length-dependent bug.
    fn bnd_model_with_real_lengths() -> SvModel {
        let mut m = SvModel::default();
        m.per_base_rate = 5e-3;
        m.homozygous_frequency = 0.5;
        m.type_probabilities.clear();
        m.type_probabilities.insert(SvType::Bnd, 1.0);
        m.length_log_normal.insert(SvType::Bnd, (7.0, 0.3));
        m
    }

    /// VCF 4.2 §5.4: a breakend is a point — END == POS, no SVLEN — whatever length the
    /// model sampled. The anchor used to carry `END = POS + length` while its own mate
    /// carried `END = POS`, so one junction shipped two conventions and tools that size
    /// records off END (truvari, bcftools) would read a multi-kb span for a 1bp event.
    /// Invisible to every other test because `fit_from_observations` pins BND length to
    /// 1 — but `SvModel` is plain serde JSON, and grafted or hand-edited models are not.
    #[test]
    fn bnd_end_equals_pos_regardless_of_length_distribution() {
        let m = bnd_model_with_real_lengths();
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
        let mut checked = 0;
        for v in &out {
            let Some(sv) = v.alternate.as_symbolic() else {
                continue;
            };
            if sv.sv_type != SvType::Bnd {
                continue;
            }
            let pos_1based = v.location + 1;
            assert_eq!(
                sv.end,
                Some(pos_1based),
                "BND at {pos_1based} has END={:?}; a breakend is a point, END must equal POS",
                sv.end
            );
            assert_eq!(sv.svlen, None, "a breakend has no SVLEN");
            let info = v.info.as_deref().unwrap_or("");
            assert!(
                info.contains(&format!("END={pos_1based}")),
                "INFO {info} must carry END={pos_1based}"
            );
            assert!(!info.contains("SVLEN="), "INFO {info} must not carry SVLEN");
            checked += 1;
        }
        assert!(checked > 0, "no de novo BND produced — test proves nothing");
    }

    /// The two records of a junction must describe the SAME rearrangement. Checking each
    /// record's flags against its own ALT is not enough: a mate emitted as case 1
    /// (`t[p[`) against an anchor at case 2 (`t]p]`) is internally consistent on both
    /// sides and still describes two different events. That is the defect that shipped,
    /// one level up.
    #[test]
    fn denovo_bnd_mate_pair_declares_one_junction() {
        let m = bnd_model_with_real_lengths();
        let seq = acgt_sequence(100_000);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);

        let by_id: std::collections::HashMap<&str, &Variant> = out
            .iter()
            .filter(|v| {
                v.alternate
                    .as_symbolic()
                    .map(|s| s.sv_type == SvType::Bnd)
                    .unwrap_or(false)
            })
            .map(|v| (v.id.as_deref().expect("BND must carry an ID"), v))
            .collect();
        assert!(!by_id.is_empty(), "no de novo BND produced");

        let mut pairs = 0;
        for v in by_id.values() {
            let sv = v.alternate.as_symbolic().unwrap();
            let mateid = v
                .info
                .as_deref()
                .unwrap_or("")
                .split(';')
                .find_map(|f| f.strip_prefix("MATEID="))
                .expect("BND must carry MATEID");
            let mate = by_id.get(mateid).expect("MATEID must resolve");
            let mate_sv = mate.alternate.as_symbolic().unwrap();

            let (_, _, ja, mr) = parse_bnd_alt(&sv.raw_alt);
            let (_, _, ja_m, mr_m) = parse_bnd_alt(&mate_sv.raw_alt);
            assert_eq!(
                (ja, mr),
                (ja_m, mr_m),
                "junction described two ways: {} then {} — a reciprocal pair must share \
                 one geometry (VCF 4.2 §5.4's worked example is `t]p]` <-> `t]p]`)",
                sv.raw_alt,
                mate_sv.raw_alt
            );
            // ...and the flags the read generator uses must agree on both sides too,
            // or the two halves of one junction stitch different sequences.
            assert_eq!(
                (sv.bnd_join_after, sv.bnd_mate_extends_right),
                (mate_sv.bnd_join_after, mate_sv.bnd_mate_extends_right),
                "mate records carry different read-generation geometry"
            );
            pairs += 1;
        }
        assert!(pairs > 0, "no pair checked — test proves nothing");
    }

    /// VCF 4.2 puts POS on the base BEFORE a symbolic event, so `END - POS` must
    /// equal `|SVLEN|` for every span type — one rule, not one per type.
    ///
    /// `place_sv` used to return the FIRST AFFECTED base as the anchor for
    /// DUP/CNV/INV while returning the preceding base for DEL, so a single output
    /// file carried two meanings of POS: `END-POS == |SVLEN|` for DEL and
    /// `|SVLEN|-1` for the rest. It also put the first duplicated base in REF where
    /// the preceding base belongs.
    #[test]
    fn denovo_symbolic_svs_put_pos_on_the_base_before_the_event() {
        for sv_type in [SvType::Del, SvType::Dup, SvType::Cnv, SvType::Inv] {
            let m = model_with_single_type(sv_type, 5e-3);
            let seq = acgt_sequence(100_000);
            let mut rng = deterministic_rng();
            let out = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng);
            let mut checked = 0;
            for v in &out {
                let Some(sv) = v.alternate.as_symbolic() else {
                    continue;
                };
                if sv.sv_type != sv_type {
                    continue;
                }
                let pos = (v.location + 1) as i64; // 1-based POS
                let end = sv.end.expect("symbolic SV must carry END") as i64;
                let svlen = sv.svlen.expect("span SV must carry SVLEN");
                assert_eq!(
                    end - pos,
                    svlen.abs(),
                    "{sv_type:?} at POS {pos}: END-POS ({}) != |SVLEN| ({}) — POS is \
                     not on the base before the event",
                    end - pos,
                    svlen.abs()
                );
                checked += 1;
            }
            assert!(
                checked > 0,
                "{sv_type:?}: nothing sampled — test proves nothing"
            );
        }
    }

    #[test]
    fn sampled_variants_dont_overlap_each_other_or_existing_svs() {
        // Pre-load one big <DEL> at position 100 spanning ~1000bp.
        // Sample with a high rate so we exercise the overlap rejector.
        // After sampling, no two SVs (existing + sampled) share any base.
        let existing = vec![{
            let mut sd = SvData::new("<DEL>", SvType::Del);
            sd.end = Some(1100);
            Variant {
                variant_type: VariantType::Complex,
                location: 100,
                reference: vec![Nucleotide::A],
                alternate: AlternateType::Symbolic(sd),
                genotype_str: "1/1".to_string(),
                genotype: Genotype::Homozygous,
                allele_fraction: None,
                id: None,
                quality_score: None,
                filter: None,
                info: None,
                format: Vec::new(),
                sample: Vec::new(),
                provenance: Provenance::InputVcf,
            }
        }];
        let m = model_with_single_type(SvType::Dup, 1e-3);
        let seq = acgt_sequence(50_000);
        let mut rng = deterministic_rng();
        let sampled = m.sample_variants("chr1", 50_000, &existing, &seq, 2, 1.0, 0.25, &mut rng);

        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for v in existing.iter().chain(sampled.iter()) {
            if let Some(r) = affected_range_for_existing(v, 50_000) {
                ranges.push(r);
            }
        }
        for i in 0..ranges.len() {
            for j in (i + 1)..ranges.len() {
                let (s1, e1) = ranges[i];
                let (s2, e2) = ranges[j];
                assert!(
                    !(s1 < e2 && s2 < e1),
                    "SVs overlap: [{s1}, {e1}) vs [{s2}, {e2})"
                );
            }
        }
    }

    #[test]
    fn sample_variants_skips_anchors_in_n_gaps() {
        // First 10kb of the contig is all-N; sampler must place every
        // SV in the ACGT tail. With a deterministic seed and a small
        // contig + high rate, we exercise the rejection loop hard.
        let mut seq = vec![Nucleotide::N; 10_000];
        seq.extend(acgt_sequence(20_000));
        let m = model_with_single_type(SvType::Del, 2e-3);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", seq.len(), &[], &seq, 2, 1.0, 0.25, &mut rng);
        for v in &out {
            // Anchor is 0-based variant.location. Must NOT fall in the
            // all-N prefix.
            assert_ne!(seq[v.location], Nucleotide::N);
            assert!(
                v.location >= 10_000,
                "anchor {} fell in the all-N gap",
                v.location
            );
        }
    }

    #[test]
    fn anchor_window_alignable_basics() {
        // 200 bp all-ACGT → alignable everywhere except the very edges
        // where the window clips.
        let acgt = acgt_sequence(2000);
        assert!(anchor_window_alignable(&acgt, 1000, 200));
        // ±200 bp window of all-N → not alignable.
        let ns: Vec<Nucleotide> = vec![Nucleotide::N; 2000];
        assert!(!anchor_window_alignable(&ns, 1000, 200));
        // Sub-Mb N-tract pattern: anchor at a clean base sitting right
        // next to a 400-bp N-tract. The pre-v1.13.1 single-base check
        // would pass; the windowed check must reject.
        let mut mixed = acgt_sequence(1000);
        mixed.extend(vec![Nucleotide::N; 400]);
        mixed.extend(acgt_sequence(1000));
        // Anchor at 1000 (the last clean base before the N-tract).
        // ±200bp window = [800..1200]: 200bp clean + 200bp N = 50% N.
        // Threshold is 80%-non-N, so this should be REJECTED.
        assert!(!anchor_window_alignable(&mixed, 1000, 200));
    }

    #[test]
    fn sample_variants_skips_n_adjacent_clean_anchors() {
        // #224 regression: an anchor that's a clean base but sits next
        // to a long N-tract used to pass the single-base check at line
        // 473, then the chimeric pair would sample a read_len-sized
        // piece into the N-tract. Now anchor_window_alignable rejects
        // these positions too. The test contig has a 400-bp N-tract
        // sandwiched between two ACGT regions; the anchor for every
        // emitted SV must sit far enough from the N-tract that the
        // ±200bp window is mostly clean.
        let mut seq = acgt_sequence(5000);
        let n_start = 5000;
        seq.extend(vec![Nucleotide::N; 400]);
        seq.extend(acgt_sequence(10_000));
        let n_end = n_start + 400;
        let m = model_with_single_type(SvType::Del, 2e-3);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", seq.len(), &[], &seq, 2, 1.0, 0.25, &mut rng);
        for v in &out {
            let loc = v.location;
            let clipped = loc.saturating_sub(200)..loc.saturating_add(200).min(seq.len());
            let in_n_window = (clipped.start..clipped.end).any(|i| i >= n_start && i < n_end);
            assert!(
                !in_n_window || {
                    // It's fine if the window CLIPS the N-tract — as long
                    // as the majority of the window is non-N. Recheck
                    // with the same threshold the production code uses.
                    anchor_window_alignable(&seq, loc, 200)
                },
                "anchor {loc} sits in/next to the 400bp N-tract — \
                 its ±200bp window contains too many Ns"
            );
        }
    }

    #[test]
    fn sample_variants_is_deterministic_for_same_seed() {
        // Two runs with seeded RNGs must produce the exact same Variants
        // — a non-trivial property given the overlap retry loop, since
        // a wrong reuse of RNG state would diverge on the first rejection.
        let m = model_with_single_type(SvType::Del, 5e-3);
        let seq = acgt_sequence(100_000);
        let mut rng1 = deterministic_rng();
        let mut rng2 = deterministic_rng();
        let a = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng1);
        let b = m.sample_variants("chr1", 100_000, &[], &seq, 2, 1.0, 0.25, &mut rng2);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.location, y.location);
            assert_eq!(
                x.alternate.as_symbolic().map(|s| s.end),
                y.alternate.as_symbolic().map(|s| s.end)
            );
        }
    }

    #[test]
    fn build_info_field_round_trips_through_parse_sv_info() {
        use crate::structs::variants::parse_sv_info;
        // The INFO string we emit must round-trip through the existing
        // parser so the writer + downstream readers behave like
        // user-supplied input.
        let info = build_info_field(SvType::Cnv, 1500, 1000, Some(4));
        let parsed = parse_sv_info(&info);
        assert_eq!(parsed.svtype.as_deref(), Some("CNV"));
        assert_eq!(parsed.end, Some(1500));
        assert_eq!(parsed.copy_number, Some(4));
    }

    #[test]
    fn genotype_string_handles_polyploid() {
        assert_eq!(genotype_string(&Genotype::Homozygous, 2), "1/1");
        assert_eq!(genotype_string(&Genotype::Heterozygous, 2), "0/1");
        assert_eq!(genotype_string(&Genotype::Homozygous, 3), "1/1/1");
        assert_eq!(genotype_string(&Genotype::Heterozygous, 4), "0/0/0/1");
    }

    #[test]
    fn affected_range_for_existing_ins_is_point_event() {
        let mut v = sv(1000, None, SvType::Ins, Genotype::Heterozygous, None);
        if let AlternateType::Symbolic(ref mut sv) = v.alternate {
            sv.svlen = Some(500);
        }
        let range = affected_range_for_existing(&v, 10_000).expect("must have range");
        // POS=1000 (0-based in the 'sv' helper) -> location=1000
        // Insertions should only affect the anchor base in terms of reference overlap.
        assert_eq!(range, (1000, 1001));
    }

    // ── sample_poisson tests ─────────────────────────────────────────
    //
    // The hybrid sampler switches between Knuth's (λ<30) and a Gaussian
    // approximation (λ≥30). The lower branch is exercised implicitly by
    // every `sample_variants_*` test that uses small-λ models; the tests
    // below pin the upper branch, which previously errored out via
    // `exp(-λ)` underflow at λ ≈ 745 — silently disabling de-novo SV
    // generation at human-chromosome scale (regression introduced by
    // v1.10.2's full-genome `per_base_rate` refit).

    #[test]
    fn sample_poisson_large_lambda_returns_finite_count() {
        // The exact value that the v1.10.2 default per_base_rate produces
        // on chr22-scale: 4.841e-4 × 50e6 = ~24,205. Previously errored
        // via `exp(-λ)` underflow; must now succeed.
        let mut rng = deterministic_rng();
        let n = sample_poisson(24_205.0, &mut rng).expect("must not error at λ=24,205");
        // Within ±5σ of mean. σ = √24205 ≈ 156, so ±780 of 24,205.
        assert!(
            (n as f64 - 24_205.0).abs() < 780.0,
            "λ=24,205 draw {n} is more than 5σ from mean — algorithm broken?"
        );
    }

    #[test]
    fn sample_poisson_large_lambda_mean_matches() {
        // Empirical mean of many large-λ draws should approach λ.
        let mut rng = deterministic_rng();
        let lambda = 1_000.0;
        let n_draws = 200;
        let total: u64 = (0..n_draws)
            .map(|_| sample_poisson(lambda, &mut rng).unwrap())
            .sum();
        let mean = total as f64 / n_draws as f64;
        // Standard error of the mean is σ/√n = √λ/√n = √1000/√200 ≈ 2.24,
        // so allow ±10 (~4 SE).
        assert!(
            (mean - lambda).abs() < 10.0,
            "λ=1000 empirical mean over {n_draws} draws = {mean}; expected ≈ 1000"
        );
    }

    #[test]
    fn sample_poisson_extreme_lambda_does_not_overflow() {
        // chr1-scale: 4.841e-4 × 250e6 ≈ 121,000. Previously errored;
        // Gaussian approximation handles it cleanly.
        let mut rng = deterministic_rng();
        let n = sample_poisson(121_000.0, &mut rng).expect("must not error at λ=121,000");
        assert!(n > 0, "λ=121,000 must produce a positive draw");
    }

    #[test]
    fn sample_poisson_small_lambda_still_uses_knuth() {
        // Regression guard: small-λ behavior is unchanged. Knuth's
        // algorithm consumes a variable number of RNG draws (matters for
        // determinism tests downstream); the Gaussian branch consumes
        // exactly 1. If λ<30 stopped routing to Knuth, existing seeded
        // tests would silently shift output.
        let mut rng_a = deterministic_rng();
        let mut rng_b = deterministic_rng();
        let _ = sample_poisson(5.0, &mut rng_a).unwrap();
        // After the same number of Knuth iterations, both RNGs should be
        // in identical states — confirmed by drawing the same next value.
        let _ = sample_poisson(5.0, &mut rng_b).unwrap();
        let next_a = rng_a.random().unwrap();
        let next_b = rng_b.random().unwrap();
        assert_eq!(
            next_a, next_b,
            "small-λ path RNG-draw-sequence diverged between runs"
        );
    }

    #[test]
    fn sample_variants_works_at_chromosome_scale() {
        // End-to-end regression for the silent-zero bug. With the v1.10.2
        // full-genome `per_base_rate` (4.841e-4) and a chr22-sized
        // synthetic contig, the de-novo sampler must produce a non-zero
        // number of SVs within a Poisson-reasonable range.
        //
        // We don't use the bundled default model directly because its
        // length_log_normal has values up to μ=9.6 (CNVs around 15kb)
        // which would still mostly clear the contig_len/4 cap on a
        // synthetic 50Mb sequence — but to keep the test runtime
        // tractable, we hand-roll a small-length model.
        let mut m = SvModel::default();
        m.per_base_rate = 4.841e-4;
        m.homozygous_frequency = 0.2;
        m.type_probabilities.clear();
        m.type_probabilities.insert(SvType::Del, 1.0);
        m.length_log_normal.clear();
        m.length_log_normal.insert(SvType::Del, (6.5, 1.0)); // median ~665bp
        let contig_len = 50_000_000;
        let seq = acgt_sequence(contig_len);
        let mut rng = deterministic_rng();
        let out = m.sample_variants("chr1", contig_len, &[], &seq, 2, 1.0, 0.25, &mut rng);
        // Expected ≈ 24,205. Even with overlap-rejection saturation
        // capping below this (max_retries = 10 × n_sv), we should land
        // very far from zero. Anything > 1000 is a clear sign the
        // Poisson branch worked.
        assert!(
            out.len() > 1000,
            "Expected thousands of de-novo SVs at chr22 scale; got {}. \
             This used to silently emit zero due to Poisson `exp(-λ)` underflow.",
            out.len()
        );
    }

    #[test]
    fn fit_from_observations_handles_ins() {
        // Build a training set with only symbolic insertions. Fitter
        // must use SVLEN as length and produce a usable model.
        let mut obs = Vec::new();
        for i in 0..10 {
            let mut v = sv(1000 * i, None, SvType::Ins, Genotype::Heterozygous, None);
            if let AlternateType::Symbolic(ref mut sv) = v.alternate {
                sv.svlen = Some(200);
            }
            obs.push(v);
        }
        let m = SvModel::fit_from_observations(&obs, 100_000).expect("must produce a model");
        assert!(m.type_probabilities.contains_key(&SvType::Ins));
        let (mu, _) = m.length_log_normal.get(&SvType::Ins).unwrap();
        // exp(mu) should be ~200. ln(200) ≈ 5.298
        assert!((mu - 200.0f64.ln()).abs() < 1e-5);
        assert!(m.per_base_rate > 0.0);
    }

    /// De novo `<INS>` records (#190) must be emitted as LITERAL Insertion
    /// variants with a generated novel-base sequence — not symbolic
    /// `<INS>` with no allele. This pins:
    ///   - VariantType::Insertion (not Complex)
    ///   - AlternateType::Literal (not Symbolic)
    ///   - ALT bases = [anchor_base, ...novel_bases] of correct length
    ///   - Only ACGT in the novel bases (no Ns or other ambiguity)
    ///   - MutatedMap::from_interval routes the record to variant_map
    ///     (where literal insertions live) rather than sv_records
    #[test]
    fn de_novo_ins_emits_literal_variant_with_novel_sequence() {
        use crate::structs::mutated_map::MutatedMap;
        let seq = acgt_sequence(50_000);
        // INS-only model so we don't have to wade through DEL/DUP/etc.
        let m = model_with_single_type(SvType::Ins, 1e-3);
        let mut rng = deterministic_rng();
        let sampled = m.sample_variants("chr1", 50_000, &[], &seq, 2, 1.0, 0.25, &mut rng);

        assert!(!sampled.is_empty(), "expected ≥1 de novo INS");

        for v in &sampled {
            assert_eq!(
                v.variant_type,
                VariantType::Insertion,
                "INS must be a literal Insertion variant; got {:?}",
                v.variant_type
            );
            let alt = match &v.alternate {
                AlternateType::Literal(bases) => bases,
                AlternateType::Symbolic(_) => panic!("INS must be literal-ALT, not symbolic"),
            };
            assert!(
                alt.len() > 1,
                "ALT must include anchor + novel bases; got len {}",
                alt.len()
            );
            assert_eq!(
                v.reference.len(),
                1,
                "REF for an insertion is the single anchor base"
            );
            assert_eq!(
                alt[0], v.reference[0],
                "ALT[0] must equal REF (anchor base)"
            );
            for (i, &b) in alt.iter().enumerate().skip(1) {
                assert!(
                    matches!(
                        b,
                        Nucleotide::A | Nucleotide::C | Nucleotide::G | Nucleotide::T
                    ),
                    "novel-insertion base {i} must be ACGT; got {:?}",
                    b
                );
            }
        }

        // MutatedMap::from_interval routes literal Insertions to variant_map,
        // not sv_records. Confirm that's where the de novo INS lands.
        let map = MutatedMap::from_interval(0, 50_000, sampled).unwrap();
        assert!(
            map.sv_records
                .iter()
                .all(|v| !matches!(v.variant_type, VariantType::Insertion)),
            "de novo INS landed in sv_records; expected variant_map"
        );
        assert!(
            !map.variant_map.is_empty(),
            "de novo INS missing from variant_map"
        );
    }
}
