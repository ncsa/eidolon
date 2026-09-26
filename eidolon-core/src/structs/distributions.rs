//! This DiscreteDistribution is an implementation of Zach Stephen's original neat-genReads code
//! from the py/probability.py file in tag 2.1 of github.com/ncsa/neat
//! (see also github.com/zstephens/neat-genreads). We may try the statrs Categorical distribution
//! as well, as I think it does the same thing.
use crate::rng::NeatRngError;
use log::debug;
use serde::{Deserialize, Serialize};
use statrs::{
    distribution::{ContinuousCDF, Normal},
    statistics::Distribution,
};
use std::fmt::Debug;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum DistributionErrors {
    #[error("A distribution returned a sampling error: {0}.")]
    SamplingError(&'static str),
    #[error("A distribution returned an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Requested a value vector from an index-only model")]
    VectorNotFound,
    #[error("Input values and weights vector must be of the same length")]
    InputMismatchError,
    #[error(
        "Distribution weights must be finite and non-negative with a positive total; got {0:?}"
    )]
    InvalidWeights(Vec<f64>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscreteDistribution<T> {
    pub(crate) values: Vec<T>,
    pub(crate) weights: Vec<f64>,
}

impl<T: Clone + Debug + Serialize + for<'de> Deserialize<'de>> DiscreteDistribution<T> {
    pub fn new(w_vec: &Vec<f64>, v_vec: &Vec<T>) -> Result<Self, DistributionErrors> {
        if w_vec.len() != v_vec.len() {
            return Err(DistributionErrors::InputMismatchError);
        }
        // No mass means no distribution. This used to fall back to "every draw is the first
        // value", which is how an empty transition row mutated everything to A (#760). The
        // right prior depends on what the distribution is, so the caller supplies it.
        let sum_weights: f64 = w_vec.iter().sum();
        if w_vec.iter().any(|w| !w.is_finite() || *w < 0.0) || sum_weights <= 0.0 {
            return Err(DistributionErrors::InvalidWeights(w_vec.clone()));
        }
        let mut normalized_weights: Vec<f64> = w_vec.iter().map(|w| w / sum_weights).collect();
        let cumulative_probability = cumulative_sum(&mut normalized_weights)?;

        Ok(DiscreteDistribution {
            values: v_vec.to_vec(),
            weights: cumulative_probability,
        })
    }

    pub fn values(&self) -> Result<Vec<T>, DistributionErrors> {
        // fetches the value matrix
        Ok(self.values.clone())
    }

    pub fn weights(&self) -> Result<Vec<f64>, DistributionErrors> {
        // fetches the value matrix, if it exists
        Ok(self.weights.to_vec())
    }

    pub fn sample(&self, rand: f64) -> Result<T, DistributionErrors> {
        // returns a random value for the distribution, based on cumulative probability
        // This is basically an icdf for a discrete distribution
        // We take a random number as an input to avoid copying the RNG around the program.
        let mut lo: usize = 0;
        // There was no minus 1 in the original, but I ran into an issue with this, for example,
        // lo = 47, hi = 48
        // mid = (lo + hi) / 2 = 47
        // lo = mid + 1 = 48
        // now we're out of bounds
        let mut hi: usize = self.weights.len();
        // bisect left
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.weights[mid] < rand {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo == self.values.len() {
            // All right, the edge case in question is when the quaity array has no values, because the original quality array wasn't 100% filled
            // out. In the future, when we fix the quality arrays so that this doesn't happen, we will be able to remove this safety check, I think
            debug!("Errored on {:?}", self);
            return Ok(self.values[0].clone());
        }
        Ok(self.values[lo].clone())
    }
}

/// TODO: Figure out how to implement serialize and deserialize for normal. I think it should be easy
/// enough to do, but I'm tired
#[derive(Debug, Clone)]
pub struct NormalDistribution {
    distribution: Normal,
}

impl NormalDistribution {
    // This is little more than a wrapper for the statrs normal distribution. Other distributions
    // may be more complicated, as needed. I would copy in the code for custom tailoring, but the
    // inverse function requires tons of hardcoded coefficient tables
    // and no one has time for all that
    pub fn new(mean: f64, std_dev: f64) -> Result<Self, DistributionErrors> {
        Ok(NormalDistribution {
            distribution: Normal::new(mean, std_dev).unwrap(),
        })
    }

    pub fn params(&self) -> Result<(f64, f64), DistributionErrors> {
        // returns the initiating parameters of the underlying distribution
        Ok((
            self.distribution.mean().unwrap(),
            self.distribution.std_dev().unwrap(),
        ))
    }

    pub fn inverse_cdf(&self, x: f64) -> Result<f64, DistributionErrors> {
        Ok(self.distribution.inverse_cdf(x))
    }

    pub fn sample(&self, rand: f64) -> Result<f64, DistributionErrors> {
        // Takes a statrs NormalDistribution object, uses the ICDF to start with a probability (our
        // RNG, which generates numbers between 0 and 1 with approximately normal distribution) and
        // generate the corresponding Y value from the PDF. Very handy!
        Ok(self.distribution.inverse_cdf(rand))
    }
}

fn cumulative_sum(a: &mut Vec<f64>) -> Result<Vec<f64>, DistributionErrors> {
    // This calculation allows us to sample the dataset with a simple bisect,
    // which is very fast an convenient.
    let mut acc = 0.0;
    let mut cumvec = Vec::with_capacity(a.len());
    for x in a {
        acc += *x;
        cumvec.push(acc);
    }
    Ok(cumvec)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::rng::NeatRng;

    // An all-zero vector used to become "every draw is the first value", which is how an
    // empty transition row mutated everything to A (#760). There is no correct value for a
    // distribution with no mass, so it is refused; the caller supplies its own prior.
    #[test]
    fn weights_with_no_usable_mass_are_refused() {
        let values = vec![1, 2, 3, 4];
        for weights in [
            vec![0.0, 0.0, 0.0, 0.0],
            vec![],
            vec![1.0, -1.0, 1.0, 1.0],
            vec![1.0, f64::NAN, 1.0, 1.0],
            vec![1.0, f64::INFINITY, 1.0, 1.0],
        ] {
            let v = values[..weights.len()].to_vec();
            assert!(
                matches!(
                    DiscreteDistribution::new(&weights, &v),
                    Err(DistributionErrors::InvalidWeights(_))
                ),
                "weights {weights:?} should be refused"
            );
        }
    }

    #[test]
    fn zero_weights_alongside_positive_ones_are_accepted() {
        let d = DiscreteDistribution::new(&vec![0.0, 3.0, 0.0, 1.0], &vec![1, 2, 3, 4]).unwrap();
        assert_eq!(d.weights().unwrap(), vec![0.0, 0.75, 0.75, 1.0]);
    }

    #[test]
    fn test_discrete_distribution() {
        let values = vec![1, 2, 3, 4, 5, 6];
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let d = DiscreteDistribution::new(&weights, &values).unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 6);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 4);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 4);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 2);
    }

    #[test]
    fn test_discrete_distribution_with_values() {
        let weights: Vec<f64> = vec![1.1, 2.0, 1.0, 8.0, 0.2, 2.0];
        let values: Vec<usize> = vec![2, 3, 7, 8, 9, 11];
        let d = DiscreteDistribution::new(&weights, &values).unwrap();
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 11);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 8);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 8);
        let rand = rng.random().unwrap();
        let x = d.sample(rand).unwrap();
        assert_eq!(x, 3);
    }

    #[test]
    fn test_normal_distribution_params() {
        let nd = NormalDistribution::new(100.0, 15.0).unwrap();
        let (mean, std_dev) = nd.params().unwrap();
        assert_eq!(mean, 100.0);
        assert_eq!(std_dev, 15.0);
    }

    #[test]
    fn test_normal_distribution_sample_median() {
        // At p=0.5 the inverse CDF returns the mean exactly
        let nd = NormalDistribution::new(200.0, 30.0).unwrap();
        let median = nd.sample(0.5).unwrap();
        assert!((median - 200.0).abs() < 1e-9);
    }

    #[test]
    fn test_normal_distribution_sample_range() {
        // Samples between 0 and 1 should all land within ±4 standard deviations
        let nd = NormalDistribution::new(0.0, 1.0).unwrap();
        for p in [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999] {
            let val = nd.sample(p).unwrap();
            assert!(
                val.abs() < 4.0,
                "sample({}) = {} is outside ±4 sigma",
                p,
                val
            );
        }
    }

    // ── Property tests ──────────────────────────────────────────────────────
    //
    // These exercise invariants of DiscreteDistribution that hold across the space of
    // possible inputs, not just the hand-rolled fixtures above. The CDF binary search
    // in sample() is hot-path code shared by every model, so blanketing it with
    // randomized inputs is cheap insurance against off-by-one bugs.

    use proptest::collection::vec as prop_vec;
    use proptest::prelude::*;

    /// Strategy: a (weights, values) pair with positive, well-formed weights. Length
    /// 2..=10 (small enough to keep tests fast but large enough to exercise the
    /// bisect tree depth). Weights drawn from (0.0, 10.0) — strictly positive to
    /// ensure sum > 0 and the normalization branch is exercised. Values are indices.
    fn positive_weights_strategy() -> impl Strategy<Value = (Vec<f64>, Vec<usize>)> {
        prop_vec(0.001f64..10.0, 2..=10).prop_map(|weights| {
            let values: Vec<usize> = (0..weights.len()).collect();
            (weights, values)
        })
    }

    proptest! {
        /// Whatever weights and `rand` value are thrown at it, `sample` must return one of
        /// the values from the input set — never panic, never return a value not in `values`.
        /// `rand` covers the inclusive endpoints 0.0 and 1.0 explicitly (proptest does include
        /// boundary values in its shrink space).
        #[test]
        fn proptest_sample_returns_value_from_input_set(
            (weights, values) in positive_weights_strategy(),
            rand in 0.0f64..=1.0,
        ) {
            let d = DiscreteDistribution::new(&weights, &values).unwrap();
            let sampled = d.sample(rand).unwrap();
            prop_assert!(
                values.contains(&sampled),
                "sample({rand}) returned {sampled}, which is not in input values {values:?}",
            );
        }

        /// With many samples drawn from a deterministic RNG, the empirical frequency of
        /// each value must approach its normalized weight. This stresses the bisect-left
        /// implementation in `sample()` more thoroughly than the hand-rolled fixture tests
        /// — any off-by-one in the CDF construction would show up as a systematic skew.
        #[test]
        fn proptest_empirical_frequency_approaches_weights(
            (weights, values) in positive_weights_strategy(),
        ) {
            let d = DiscreteDistribution::new(&weights, &values).unwrap();
            let total_weight: f64 = weights.iter().sum();
            let mut counts = vec![0usize; values.len()];

            // Use a deterministic seed so the property is reproducible.
            let mut rng = NeatRng::new_from_seed(&vec!["proptest empirical seed".to_string()])
                .unwrap();
            let n_samples = 20_000usize;
            for _ in 0..n_samples {
                let s = d.sample(rng.random().unwrap()).unwrap();
                counts[s] += 1;
            }

            // Tolerance: with N=20k draws and worst-case probability ~0.5, the standard
            // deviation of a count is sqrt(N*p*(1-p)) ≈ 70 (~0.35% of N). A 3% slack on
            // *N* (equivalent to ~9σ) is comfortable for any individual bin without being
            // so loose that systematic bias slips through.
            let tolerance_count = (0.03 * n_samples as f64) as usize;
            for i in 0..values.len() {
                let expected = (weights[i] / total_weight) * n_samples as f64;
                let actual = counts[i] as f64;
                prop_assert!(
                    (actual - expected).abs() < tolerance_count as f64,
                    "bin {i}: expected ~{expected:.1} samples, got {actual} \
                     (weights={weights:?}, tolerance={tolerance_count})",
                );
            }
        }

        /// All-zero weights are refused whatever their length (#760).
        #[test]
        fn proptest_all_zero_weights_are_refused(
            values in prop_vec(any::<u32>(), 2..=10),
        ) {
            let weights = vec![0.0f64; values.len()];
            prop_assert!(matches!(
                DiscreteDistribution::new(&weights, &values),
                Err(DistributionErrors::InvalidWeights(_))
            ));
        }
    }

    #[test]
    fn test_normal_distribution_inverse_cdf_symmetry() {
        // inverse_cdf(p) and inverse_cdf(1-p) should be symmetric around the mean
        let nd = NormalDistribution::new(50.0, 10.0).unwrap();
        let lo = nd.inverse_cdf(0.1).unwrap();
        let hi = nd.inverse_cdf(0.9).unwrap();
        assert!((lo + hi - 100.0).abs() < 1e-9);
    }
}
