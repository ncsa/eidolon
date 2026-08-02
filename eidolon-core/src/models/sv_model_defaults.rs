//! Bundled default [`SvModel`] used by [`crate::models::mutation_model::MutationModel::default`].
//!
//! Most parameters are **refit** from a full-genome validation run
//! against the gnomAD-SV v4.1 sites VCF (1,476,325 surviving SVs across
//! all autosomes/sex chromosomes after the standard INV/CPX/CTX/BND/INS:ME
//! filter; `tools/validate_with_gnomad.sh full` produces the same numbers).
//! The two parameters that can't be sourced from a sites-only VCF —
//! `homozygous_frequency` and `cnv_copy_number_distribution` — are kept
//! at their original literature-derived approximations (gnomAD-SV
//! per-sample CN lives in FORMAT/SAMPLE, not INFO/CN; het:hom ≈ 4:1 is
//! taken from the Collins et al. paper, *Nature* 581, 444–451, 2020).
//!
//! Earlier revisions of this file sourced the same parameters from a
//! chr22-only smoke run. The chr22 fit reported `per_base_rate ≈ 8.9e-6`
//! because the trainer normalizes by the *whole* FASTA length, not by
//! the simulated chromosome — so the chr22 default was effectively
//! "fraction of all-genome SVs that happen to fall on chr22," not a
//! true per-base rate. The full-genome refit corrects that.
//!
//! Anyone doing serious benchmarking should train their own model with
//! `gen-mut-model` against the SV-rich VCF of their choice and override
//! the embedded default by passing `mutation_model: <trained>.json.gz`
//! in their gen-reads YAML.

use crate::structs::sv_model::SvModel;
use crate::structs::variants::SvType;
use std::collections::HashMap;

/// Build the bundled default `SvModel`.
///
/// Returns a fully-populated, [`SvModel::is_usable`]-passing model with
/// `<DEL>`, `<DUP>`, `<CNV>`, and `<BND>` entries. `gen-reads` consults this only
/// when `sv_rate_scale > 0.0`, so generation remains opt-in even though
/// the model itself is always loaded.
pub fn default_sv_model() -> SvModel {
    // Per-base rate: fitted from the full gnomAD-SV v4.1 sites VCF
    // (1,832,360 surviving SVs over the GRCh38 primary assembly ~3.05 Gb).
    // The trainer computes `n / total_reflen_seen`, which gen-reads then
    // multiplies by `contig_len * sv_rate_scale` to determine the
    // expected per-contig SV count.
    let per_base_rate = 6.009e-4;

    // Type breakdown: fitted from the full gnomAD-SV v4.1 run.
    // Deletions dominate; BNDs are included at ~19%.
    //
    // ⚠ UNVERIFIED PREMISE (2026-08-02). This 0.1943 is a HEURISTIC, not a gnomAD
    // measurement: `tools/validate_with_gnomad.sh:11` drops INV/CPX/CTX/BND before
    // fitting, so no BND ever reached this number. Two consequences worth checking
    // before it is trusted:
    //
    //   1. gnomAD-SV lists CTX (translocation) SEPARATELY from BND, so in that
    //      taxonomy "BND" means an unresolved junction — many of them
    //      intra-chromosomal — and NOT a translocation. eidolon now emits every BND
    //      as an inter-chromosomal translocation, a convention taken from PCAWG,
    //      where `svclass == "TRA"` is 100% inter-chromosomal
    //      (docs/pcawg_sv_measurement.md M1). That measurement is about the CANCER
    //      corpus and does not license the same reading of this germline share.
    //   2. MEASURED 2026-08-02, streaming gnomAD-SV v4.1 sites (2,154,486 records):
    //        CTX               99  (0.0046%)  — 99/99 inter-chromosomal
    //        BND          356,035  (16.5%)    — 210,397 inter (59%), 145,638 intra (41%)
    //        DEL 56.0%  INS 14.1%  DUP 12.5%  CPX 0.7%  INV 0.10%  CNV 0.03%
    //      So gnomAD's BND is genuinely mixed and is NOT a translocation class; CTX is,
    //      and it is vanishingly rare. Depending on which eidolon means:
    //        BND == resolved translocation (CTX)          -> 0.005%, this is ~4,000x high
    //        BND == any inter-chromosomal junction        -> 9.8%,   this is ~2x high
    //      Either way 0.1943 is not supportable, and 41% of germline breakends are
    //      intra-chromosomal — a class eidolon now emits none of.
    //
    // Flagged rather than silently corrected: guessing a replacement number here
    // would repeat the exact defect (see docs/claude_engineering_audit.md §5.1).
    // Inversions and Insertions are added here at nominal rates (heuristic) since they
    // were filtered in the original validation fit.
    let mut type_probabilities: HashMap<SvType, f64> = HashMap::new();
    type_probabilities.insert(SvType::Del, 0.6433);
    type_probabilities.insert(SvType::Dup, 0.1470);
    type_probabilities.insert(SvType::Cnv, 0.0004);
    type_probabilities.insert(SvType::Bnd, 0.1943);
    type_probabilities.insert(SvType::Inv, 0.0050);
    type_probabilities.insert(SvType::Ins, 0.0100);

    // Length log-normals: fitted from the full-genome run via MLE on
    // ln(span). Quick gut-check on the medians (= exp(mu)):
    //   DEL: median ≈ 704 bp
    //   DUP: median ≈ 863 bp
    //   CNV: median ≈ 15.3 kb
    //   BND: nominal 1 bp (fixed)
    //   INV: median ≈ 4.9 kb (heuristic)
    //   INS: median ≈ 300 bp (heuristic, matching Alu MEI)
    // Sigmas reflect the heavy right tail of DUPs (multi-megabase calls
    // pull sigma out to ~2.5) and the tighter clustering of DELs / CNVs.
    let mut length_log_normal: HashMap<SvType, (f64, f64)> = HashMap::new();
    length_log_normal.insert(SvType::Del, (6.557, 1.543));
    length_log_normal.insert(SvType::Dup, (6.761, 2.530));
    length_log_normal.insert(SvType::Cnv, (9.637, 1.152));
    length_log_normal.insert(SvType::Bnd, (0.0, 0.0));
    length_log_normal.insert(SvType::Inv, (8.5, 1.5));
    length_log_normal.insert(SvType::Ins, (5.7, 1.0));

    // Copy-number distribution for `<CNV>` records. gnomAD-SV's sites
    // VCF leaves `INFO/CN` unset (CN varies per-sample, in
    // FORMAT/SAMPLE), so the v1.10 chr22 validation run produced an
    // empty cnv_copy_number_distribution from training. We keep the
    // literature-derived distribution here so a user loading
    // MutationModel::default() with `sv_rate_scale > 0` still gets
    // CNVs with sensible CN values. SvModel::sample_variants also has
    // a `FALLBACK_CN_DISTRIBUTION` constant that mirrors this for the
    // case where the *user-trained* model has an empty CN dist —
    // keep both in sync.
    //
    // Rough per-CN proportions from Collins et al.: 0 → 15%, 1 → 35%,
    // 3 → 30%, 4 → 15%, 5+ → 5%. Compressing the tail into CN=5 keeps
    // the table small without changing the simulated coverage signal
    // much.
    let mut cnv_copy_number_distribution: HashMap<u32, f64> = HashMap::new();
    cnv_copy_number_distribution.insert(0, 0.15);
    cnv_copy_number_distribution.insert(1, 0.35);
    cnv_copy_number_distribution.insert(3, 0.30);
    cnv_copy_number_distribution.insert(4, 0.15);
    cnv_copy_number_distribution.insert(5, 0.05);

    // Homozygous fraction: kept at 0.20 (Collins et al.; gnomAD-SV
    // het:hom ≈ 4:1). The v1.10 chr22 validation run could not fit
    // this from data because gnomAD-SV is sites-only with no GT
    // column — the reader defaults to heterozygous and the trainer
    // reports `homozygous_frequency = 0`. A per-sample callset could
    // refit this; for now the literature value is the best estimate.
    let homozygous_frequency = 0.20;

    // Breakend geometry: uniform, and deliberately not a literature value. Unlike the
    // homozygous fraction above there is no figure here I can cite and defend, and this
    // repo has already been burned by an unexamined inherited number. Uniform states
    // "no information"; a model fitted from a corpus containing breakends replaces it.
    let bnd_geometry_weights = crate::structs::sv_model::default_bnd_geometry_weights();

    SvModel {
        per_base_rate,
        type_probabilities,
        length_log_normal,
        cnv_copy_number_distribution,
        bnd_geometry_weights,
        homozygous_frequency,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_sv_model_is_usable() {
        // Smoke check: the bundled default must pass the same usability
        // bar the sampler keys off, otherwise `MutationModel::default()`
        // returns a model whose `sv_model: Some(...)` would silently
        // produce zero SVs even with `sv_rate_scale: 1.0`.
        assert!(default_sv_model().is_usable());
    }

    #[test]
    fn default_sv_model_carries_all_six_supported_types() {
        // Sampler skips types missing from `type_probabilities`. The
        // default must include DEL, DUP, CNV, BND, INV, and INS so users see all six
        // when they enable generation.
        let m = default_sv_model();
        for t in [
            SvType::Del,
            SvType::Dup,
            SvType::Cnv,
            SvType::Bnd,
            SvType::Inv,
            SvType::Ins,
        ] {
            assert!(
                m.type_probabilities.contains_key(&t),
                "missing type probability for {:?}",
                t
            );
            assert!(
                m.length_log_normal.contains_key(&t),
                "missing length distribution for {:?}",
                t
            );
        }
    }

    #[test]
    fn default_sv_model_type_probabilities_sum_to_one() {
        // Within float tolerance — if the hand-rolled probabilities
        // ever drift, fail loudly rather than silently sample from a
        // mis-normalized distribution.
        let total: f64 = default_sv_model().type_probabilities.values().sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "type_probabilities sum to {total}, expected 1.0"
        );
    }

    #[test]
    fn default_sv_model_cn_distribution_sums_to_one() {
        let total: f64 = default_sv_model()
            .cnv_copy_number_distribution
            .values()
            .sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "cnv_copy_number_distribution sums to {total}, expected 1.0"
        );
    }

    #[test]
    fn default_sv_model_parameters_in_plausible_ranges() {
        // Cheap sanity check — catches an accidental decimal-place
        // typo that would otherwise sail through type-checking.
        let m = default_sv_model();
        assert!(m.per_base_rate > 0.0 && m.per_base_rate < 1e-3);
        assert!(m.homozygous_frequency >= 0.0 && m.homozygous_frequency <= 1.0);
        for (t, (mu, sigma)) in m.length_log_normal.iter() {
            // mu is ln(median bp); medians range from ~tens to ~tens-of-thousands.
            // BND has a fixed nominal mu=0.
            if *t == SvType::Bnd {
                assert_eq!(*mu, 0.0);
                assert_eq!(*sigma, 0.0);
            } else {
                assert!(*mu > 2.0 && *mu < 14.0, "mu={mu} out of plausible range");
                assert!(
                    *sigma > 0.0 && *sigma < 4.0,
                    "sigma={sigma} out of plausible range"
                );
            }
        }
        for cn in m.cnv_copy_number_distribution.keys() {
            // Excluding CN=2 (which is reference for diploid) is a
            // training-time invariant — the default should respect it.
            assert_ne!(*cn, 2, "default CN distribution must not include CN=2");
        }
    }
}
