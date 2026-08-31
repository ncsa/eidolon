//! In DNA sequencing, a fragment is a bit of DNA, roughly uniform in length that is sequenced
//! by the machine. Sometimes these fragments have special molecules attached to the end for ID
//! purposes. How this is done is a process called Chemistry Magic. For our purposes, we expect
//! the data to be uniform enough that a mean and standard deviation will describe the set.

use crate::models::lib::{model_reader, model_writer};
use crate::rng::NeatRngError;
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors, NormalDistribution};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use serde_json;
use std::io;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FragmentModelError {
    #[error("Fragment model error: {0}")]
    FragModelError(&'static str),
    #[error("Fragment model returned an RNG error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Fragment model returned an IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("Fragment Model attempted to load a file that it could not find: {0}")]
    FileNotFound(String),
    #[error("Fragment model reported a distribution initiation error: {0}")]
    DistributionInitError(#[from] DistributionErrors),
    #[error("Error building default model!")]
    SerdeError(#[from] serde_json::Error),
    #[error("Error creating fragments of sufficient size. Check fragment length model.")]
    FragGenerationError,
}

#[derive(Debug, Serialize, Deserialize)]
pub enum FragmentLengthModel {
    Discrete {
        distribution: DiscreteDistribution<usize>,
    },
    Normal {
        mean: f64,
        st_dev: f64,
    },
}

/// The shipped fallback, used when a config supplies neither `fragment_model` nor
/// `fragment_mean`. Built from HCC1395 normal (SEQC2 `WGS_NS_N_1`, NovaSeq, chr20/21/22,
/// 32.6M pairs) and cross-validated against chr1 of the same library.
///
/// Fragment length is set by library chemistry, so this is a real distribution rather than
/// a universal one -- see `model_data/README.md` for full provenance and for how to build
/// your own. `the_shipped_default_is_a_usable_fragment_distribution` below asserts the
/// properties this file has to have.
static DATA_FILE: &[u8] = include_bytes!("model_data/default_fragment_length_model.json.gz");

impl FragmentLengthModel {
    pub fn new_discrete(
        lengths: Vec<usize>,
        weights: Vec<f64>,
    ) -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        let fragment_dist = DiscreteDistribution::new(&weights, &lengths)?;
        Ok(FragmentLengthModel::Discrete {
            distribution: fragment_dist,
        })
    }

    pub fn new_normal(fragment_mean: f64, fragment_std: f64) -> Result<Self, FragmentModelError> {
        Ok(FragmentLengthModel::Normal {
            mean: fragment_mean,
            st_dev: fragment_std,
        })
    }

    // Returns Result because it deserializes an embedded model file; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, FragmentModelError> {
        // The parameters of the default model from the original neat
        // The lengths range from 1 to 799, though it skips from 1 to 32 before counting up.
        // The weights are just numbers between 0 and 1.
        // These are data gathered from publicly availble human data, but should reflect
        // whatever chemistry was used at the time
        let reader = GzDecoder::new(DATA_FILE);
        let data: FragmentLengthModel =
            serde_json::from_reader(reader).map_err(FragmentModelError::SerdeError)?;
        Ok(data)
    }

    pub fn default_normal() -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        Ok(FragmentLengthModel::Normal {
            mean: 300.0,
            st_dev: 30.0,
        })
    }

    pub fn discrete_from_file(filename: &PathBuf) -> Result<Self, FragmentModelError> {
        // The baseline model is really just a mathematical equation and can be reconstructed by the two input parameters
        // But reading from data, it may be better to use the discrete distribution to maintain outliers. This will load such
        // a model from a json file. The file must be of the format:
        // {
        //   "fragment_lengths": [ ... ],
        //   "fragment_weights": [ ... ]
        // }
        // Where fragment lengths are of type usize and fragment weights are of (or can be cast as) type f64.
        if !filename.exists() {
            return Err(FragmentModelError::FileNotFound(
                filename.display().to_string(),
            ));
        }

        let data: FragmentLengthModel = model_reader(filename).unwrap();

        Ok(data)
    }

    /// Mean fragment length, exact and deterministic (no RNG draw) for either
    /// variant -- `Normal`'s own parameter, or `Discrete`'s weighted average
    /// over its value/weight table. Used to scale a boundary-extension
    /// dilution correction to the model actually in use, including a trained
    /// (real-BAM-derived) discrete model, not just the synthetic Normal case.
    pub fn mean(&self) -> Result<f64, FragmentModelError> {
        match self {
            FragmentLengthModel::Normal { mean, .. } => Ok(*mean),
            FragmentLengthModel::Discrete { distribution } => {
                let values = distribution.values()?;
                let weights = distribution.weights()?;
                let total_weight: f64 = weights.iter().sum();
                if total_weight <= 0.0 {
                    return Err(FragmentModelError::FragModelError(
                        "Discrete fragment model has zero total weight",
                    ));
                }
                let weighted_sum: f64 = values
                    .iter()
                    .zip(weights.iter())
                    .map(|(&v, &w)| v as f64 * w)
                    .sum();
                Ok(weighted_sum / total_weight)
            }
        }
    }

    pub fn normal_params(&self) -> Result<(f64, f64), FragmentModelError> {
        // This returns the parameters used to initiate the model
        match self {
            FragmentLengthModel::Discrete { distribution: _ } => {
                Err(FragmentModelError::FragModelError(
                    "Called normal_params on a discrete fragment model",
                ))
            }
            FragmentLengthModel::Normal { mean, st_dev } => Ok((*mean, *st_dev)),
        }
    }

    pub fn generate_fragment(&self, rand: f64) -> Result<usize, FragmentModelError> {
        // This function generates a fragment length based on mean and standard deviation,
        // or based on a discrete distribution.
        match self {
            // The discrete one is pretty easy
            Self::Discrete { distribution } => Ok(distribution.sample(rand)?),
            // for normal we have to build the distribution, then sample. Not sure if this will be a pinch point.
            Self::Normal { mean, st_dev } => {
                let distribution = NormalDistribution::new(*mean, *st_dev)?;
                Ok(distribution.sample(rand)?.trunc() as usize)
            }
        }
    }

    pub fn write_file(&self, filename: &PathBuf) -> Result<(), FragmentModelError> {
        // serialize a model with serde and write it to file
        model_writer(self, filename)?;
        Ok(())
    }

    #[allow(dead_code)]
    fn is_discrete(&self) -> bool {
        // This is mainly for testing purposes
        if let FragmentLengthModel::Discrete { distribution } = self {
            !distribution.weights.is_empty()
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_default() {
        let model = FragmentLengthModel::default_normal().unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert_eq!(mean, 300.0);
                assert_eq!(st_dev, 30.0);
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn the_shipped_default_is_a_usable_fragment_distribution() {
        // Replaces a test that pinned all 766 values of the old default as a literal
        // array. That asserted the bytes had not changed; it said nothing about whether
        // they were any good -- and they were not. The model it pinned was LEFT-skewed
        // (-0.434) where every real size-selected library is right-skewed, truncated at
        // 799, had 33 integer lengths inside its own range with no bin at all, and carried
        // an isolated spike at fragment length 1 holding 2.3e-09 of the mass with a 30-wide
        // hole above it. Every one of those passed that test.
        //
        // These assertions are about whether the shipped default can be USED. Provenance
        // for the current file is in model_data/README.md.
        let model = FragmentLengthModel::default().unwrap();
        let FragmentLengthModel::Discrete { distribution } = model else {
            panic!(
                "the shipped default must be Discrete; a Normal cannot hold a real library's shape"
            )
        };
        let values = distribution.values().unwrap();
        let cumulative = distribution.weights().unwrap();
        assert_eq!(values.len(), cumulative.len());
        assert!(
            values.len() > 100,
            "too few bins to describe a distribution"
        );

        // NO GAPS. A hole in the support is a length the simulator can never produce
        // sitting inside a range where the real library has fragments.
        for pair in values.windows(2) {
            assert_eq!(
                pair[1],
                pair[0] + 1,
                "gap in the shipped default between {} and {}",
                pair[0],
                pair[1]
            );
        }

        // A valid inverse-CDF: non-decreasing, ending at 1.
        for w in cumulative.windows(2) {
            assert!(w[1] >= w[0], "cumulative weights must not decrease");
        }
        assert!(
            (cumulative[cumulative.len() - 1] - 1.0).abs() < 1e-9,
            "cumulative weights must end at 1.0, got {}",
            cumulative[cumulative.len() - 1]
        );

        // Per-bin mass, recovered from the cumulative array.
        let mut mass = Vec::with_capacity(cumulative.len());
        for (i, &c) in cumulative.iter().enumerate() {
            mass.push(if i == 0 { c } else { c - cumulative[i - 1] });
        }
        assert!(
            mass.iter().all(|&m| m > 0.0),
            "no bin may hold zero mass -- that is a gap wearing a bin's clothes"
        );

        // NO ISOLATED SPIKE AT THE BOTTOM. The old default's first bin was 4000x lighter
        // than the second and 31 lengths away from it: one stray read surviving the filter.
        // A real distribution's edge is a continuum.
        let ratio = mass[1] / mass[0];
        assert!(
            (0.05..=20.0).contains(&ratio),
            "the first two bins differ by {ratio:.0}x -- the support starts on an outlier, not on data"
        );

        let total: f64 = mass.iter().sum();
        let mean: f64 = values
            .iter()
            .zip(&mass)
            .map(|(&v, &m)| v as f64 * m)
            .sum::<f64>()
            / total;
        let sd = (values
            .iter()
            .zip(&mass)
            .map(|(&v, &m)| (v as f64 - mean).powi(2) * m)
            .sum::<f64>()
            / total)
            .sqrt();
        let skew = values
            .iter()
            .zip(&mass)
            .map(|(&v, &m)| ((v as f64 - mean) / sd).powi(3) * m)
            .sum::<f64>()
            / total;

        assert!(
            (150.0..=800.0).contains(&mean),
            "mean fragment length {mean:.1} is outside anything a paired-end library produces"
        );
        assert!(sd > 10.0, "a real library has spread; sd {sd:.1} does not");
        // A LEFT tail means the distribution was cut off at the top -- which is what the
        // old default looked like, and why its p50 sat 130 bp above the real library's.
        // Not asserting a specific positive value: a different chemistry can be flatter.
        assert!(
            skew > -0.1,
            "skew {skew:+.3} suggests the distribution is truncated at its upper end"
        );

        // Most of the mass must be usable at a common read length. Fragments shorter than
        // the reads get rejected in generate_fragments, and a model that is mostly
        // unusable produces silent under-coverage.
        let unusable: f64 = values
            .iter()
            .zip(&mass)
            .filter(|&(&v, _)| v < 161)
            .map(|(_, &m)| m)
            .sum::<f64>()
            / total;
        assert!(
            unusable < 0.05,
            "{:.2}% of the default's mass is below a 2x151 read pair's minimum",
            unusable * 100.0
        );
    }

    #[test]
    fn test_new_from_mean() {
        let mean = 34.33;
        let std_dev = 1.232;
        let model = FragmentLengthModel::new_normal(mean, std_dev).unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert_eq!(mean, 34.33);
                assert_eq!(st_dev, 1.232);
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_new_discrete() {
        let l_vec = vec![1, 8, 9, 10];
        let w_vec = vec![1.0, 3.0, 2.0, 1.2];
        let model = FragmentLengthModel::new_discrete(l_vec.clone(), w_vec.clone()).unwrap();
        match model {
            FragmentLengthModel::Discrete { distribution } => {
                assert_eq!(distribution.values().unwrap(), l_vec);
                assert_eq!(
                    distribution.weights().unwrap(),
                    [
                        0.1388888888888889,
                        0.5555555555555556,
                        0.8333333333333334,
                        1.0
                    ]
                );
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_new_normal() {
        let mean = 34.33;
        let std_dev = 1.232;
        let model = FragmentLengthModel::new_normal(mean, std_dev).unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert_eq!((mean, st_dev), (34.33, 1.232));
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn generate_fragment_is_the_inverse_cdf_of_the_model() {
        // Was `assert_eq!(model.generate_fragment(0.1), 295)`. That pinned one output of
        // one model, so it broke the moment the shipped default changed and said nothing
        // about whether sampling was correct in the first place. What generate_fragment
        // OWES its caller is the inverse CDF: feed it a uniform draw, get back that
        // quantile of the distribution. That is checkable against the model itself.
        let model = FragmentLengthModel::default().unwrap();
        let FragmentLengthModel::Discrete { ref distribution } = model else {
            panic!("expected the shipped default to be Discrete")
        };
        let values = distribution.values().unwrap();
        let cumulative = distribution.weights().unwrap();

        // The quantile the model itself says sits at q, computed here by walking the
        // cumulative array -- not by asking the code under test.
        let expected_at = |q: f64| -> usize {
            for (i, &c) in cumulative.iter().enumerate() {
                if c >= q {
                    return values[i];
                }
            }
            values[values.len() - 1]
        };

        for q in [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99] {
            let got = model.generate_fragment(q).unwrap();
            assert_eq!(
                got,
                expected_at(q),
                "generate_fragment({q}) returned {got}, but the model's own CDF puts that \
                 quantile at {}",
                expected_at(q)
            );
        }

        // Monotone: a larger uniform draw can never yield a shorter fragment. A sampler
        // that ignored its input entirely would still pass the equality checks above if
        // the CDF walk had the same bug, so assert the ordering independently.
        let mut prev = 0usize;
        for i in 1..100 {
            let f = model.generate_fragment(i as f64 / 100.0).unwrap();
            assert!(f >= prev, "sampling is not monotone: {prev} then {f}");
            prev = f;
        }
        assert!(prev > values[0], "sampling never left the first bin");
    }

    #[test]
    fn test_from_file() {
        let model = FragmentLengthModel::default().unwrap();
        let temp_dir = tempfile::tempdir().unwrap();
        let mut temp_file = PathBuf::from(temp_dir.path());
        let filename = "model_test.json";
        temp_file.push(filename);
        model
            .write_file(&temp_file)
            .expect("write_file should succeed");
        let model2 = FragmentLengthModel::discrete_from_file(&temp_file).unwrap();
        assert_eq!(model.is_discrete(), model2.is_discrete());
        temp_dir.close().unwrap();
    }
}
