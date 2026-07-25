//! Subclonal architecture for de-novo somatic variants (#405, part of realism epic #311).
//!
//! A tumor is a mixture of cell populations ("subclones"), each present at a
//! distinct **cancer-cell fraction (CCF)** — the fraction of tumor cells carrying
//! that subclone's mutations. `gen-cancer-reads` today models a tumor with a single
//! global `purity` scalar, so every somatic variant collapses to essentially one
//! effective VAF. A [`SubcloneModel`] restores that missing axis: each de-novo
//! somatic variant is assigned to a subclone (weighted draw) and takes the
//! subclone's CCF as its per-variant [`allele_fraction`](eidolon_core::structs::variants::Variant::allele_fraction).
//!
//! A subclone's CCF is a **cellular-fraction factor** that *composes* — it does not
//! replace the variant's allele dosage or tumor purity. For a variant at dosage `d`
//! (alt copies / ploidy) assigned to a subclone at CCF `f`, the caller stamps
//! `allele_fraction = d × f`, and purity — realized as the tumor/normal coverage
//! split at read-mixing time — contributes the final factor:
//!
//! ```text
//! observed VAF = purity × dosage × CCF
//! ```
//!
//! So a heterozygous somatic SNV (d = 0.5) at CCF `f` lands at `f/2` — the value a
//! subclonal-deconvolution tool inverts back to `f`. These are orthogonal axes
//! (purity = normal contamination; dosage = per-copy multiplicity; CCF = subclonal
//! cellular fraction). This model owns only the CCF factor; dosage composition lives
//! at the call site (`Variant::dosage_fraction`), so polyploid dosage (#266/#267)
//! flows in for free. No new read-mixing math is needed.
//!
//! This is the generative half of #405. The reproductive half (honor per-variant CCF
//! from an input somatic VCF) rides the existing `INFO/AF` → `allele_fraction` parse
//! path and is tracked separately. The generic AF-sampling superset is #404.

use eidolon_core::rng::NeatRng;

/// One tumor subclone: a cancer-cell fraction and its relative share of the
/// somatic variant burden.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subclone {
    /// Cancer-cell fraction in `(0.0, 1.0]`. `1.0` is a clonal (truncal) population
    /// present in every tumor cell.
    pub ccf: f64,
    /// Relative weight — the (unnormalized) share of de-novo somatic variants
    /// assigned to this subclone. Weights are normalized across the model, so only
    /// their ratios matter.
    pub weight: f64,
}

/// A tumor's subclonal architecture: a non-empty set of [`Subclone`]s over which
/// de-novo somatic variants are distributed.
#[derive(Debug, Clone, PartialEq)]
pub struct SubcloneModel {
    subclones: Vec<Subclone>,
}

/// Reasons a [`SubcloneModel`] cannot be constructed.
#[derive(Debug, Clone, PartialEq)]
pub enum SubcloneModelError {
    /// No subclones were supplied.
    Empty,
    /// A CCF fell outside `(0.0, 1.0]`.
    CcfOutOfRange(f64),
    /// A weight was non-positive (or NaN).
    NonPositiveWeight(f64),
}

impl std::fmt::Display for SubcloneModelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SubcloneModelError::Empty => write!(f, "subclones must not be empty"),
            SubcloneModelError::CcfOutOfRange(v) => {
                write!(f, "subclone ccf {v} out of range (0.0, 1.0]")
            }
            SubcloneModelError::NonPositiveWeight(v) => {
                write!(f, "subclone weight {v} must be positive")
            }
        }
    }
}

impl std::error::Error for SubcloneModelError {}

impl SubcloneModel {
    /// Build a validated model. Rejects an empty set, any CCF outside `(0.0, 1.0]`,
    /// and any non-positive/NaN weight.
    pub fn new(subclones: Vec<Subclone>) -> Result<Self, SubcloneModelError> {
        if subclones.is_empty() {
            return Err(SubcloneModelError::Empty);
        }
        for s in &subclones {
            if !(s.ccf > 0.0 && s.ccf <= 1.0) {
                return Err(SubcloneModelError::CcfOutOfRange(s.ccf));
            }
            // Reject non-finite (NaN/inf) and non-positive weights: a weight must be a
            // real, positive share. `<= 0.0` alone would let NaN slip through.
            if !s.weight.is_finite() || s.weight <= 0.0 {
                return Err(SubcloneModelError::NonPositiveWeight(s.weight));
            }
        }
        Ok(SubcloneModel { subclones })
    }

    /// The subclones, in construction order.
    pub fn subclones(&self) -> &[Subclone] {
        &self.subclones
    }

    /// Draw a CCF for one de-novo somatic variant, weighted by subclone share.
    ///
    /// Consumes exactly **one** `rng.random()` draw regardless of subclone count, so
    /// the RNG-stream perturbation is a single, predictable step per stamped variant.
    /// The final subclone catches any floating-point residue, so a valid model always
    /// returns a CCF.
    pub fn sample_ccf(&self, rng: &mut NeatRng) -> Result<f64, eidolon_core::rng::NeatRngError> {
        let total: f64 = self.subclones.iter().map(|s| s.weight).sum();
        let mut r = rng.random()? * total;
        for s in &self.subclones {
            r -= s.weight;
            if r < 0.0 {
                return Ok(s.ccf);
            }
        }
        // Unreachable for a validated (non-empty, positive-weight) model except via
        // floating-point rounding at the very top of the range; fall back to the last.
        Ok(self.subclones.last().expect("validated non-empty").ccf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_empty() {
        assert_eq!(SubcloneModel::new(vec![]), Err(SubcloneModelError::Empty));
    }

    #[test]
    fn rejects_bad_ccf() {
        let bad = vec![Subclone {
            ccf: 1.5,
            weight: 1.0,
        }];
        assert_eq!(
            SubcloneModel::new(bad),
            Err(SubcloneModelError::CcfOutOfRange(1.5))
        );
        let zero = vec![Subclone {
            ccf: 0.0,
            weight: 1.0,
        }];
        assert_eq!(
            SubcloneModel::new(zero),
            Err(SubcloneModelError::CcfOutOfRange(0.0))
        );
    }

    #[test]
    fn rejects_bad_weight() {
        let bad = vec![Subclone {
            ccf: 0.5,
            weight: 0.0,
        }];
        assert_eq!(
            SubcloneModel::new(bad),
            Err(SubcloneModelError::NonPositiveWeight(0.0))
        );
    }

    #[test]
    fn sample_ccf_respects_weights() {
        // Clonal 1.0 (weight 3), minor 0.2 (weight 1) → ~75% / ~25% split.
        let model = SubcloneModel::new(vec![
            Subclone {
                ccf: 1.0,
                weight: 3.0,
            },
            Subclone {
                ccf: 0.2,
                weight: 1.0,
            },
        ])
        .unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["ccf-test".to_string()]).unwrap();

        let n = 20_000;
        let mut clonal = 0;
        for _ in 0..n {
            let ccf = model.sample_ccf(&mut rng).unwrap();
            assert!(ccf == 1.0 || ccf == 0.2, "unexpected ccf {ccf}");
            if ccf == 1.0 {
                clonal += 1;
            }
        }
        let frac = clonal as f64 / n as f64;
        // Expected 0.75; allow generous sampling slack.
        assert!(
            (frac - 0.75).abs() < 0.03,
            "clonal fraction {frac} should be near 0.75"
        );
    }

    #[test]
    fn sample_ccf_single_subclone_is_deterministic() {
        let model = SubcloneModel::new(vec![Subclone {
            ccf: 0.4,
            weight: 1.0,
        }])
        .unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["one".to_string()]).unwrap();
        for _ in 0..100 {
            assert_eq!(model.sample_ccf(&mut rng).unwrap(), 0.4);
        }
    }
}
