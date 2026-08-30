//! Walking through Zac Stephens original algorithm, to try to make sure I replicate it correctly.
//!
//! * For position 1, there is a vector of weights for each score, extracted from data.
//! * For each position in the read length after that
//!   * For each possible quality score, a distribution is constructed with weights and
//!     scores, as determined by a matrix of weights
//! * For read length N and # possible quality scores Q, this creates a vector with length N
//!   * first element is a 1-D vector of weights with length Q
//!   * each subsequent element is a vector of length Q,
//!     each element of which is a vector of length N.
//!
//! To generate quality scores, they follow the following procedure:
//!
//! * Sample the first element (1D vector) for initial quality score.
//! * For next position, the previous Q score determines which N-length set of weights to use to
//!   determine the next quality score
//!
//! Advantages of this approach:
//!
//! * Does a fairly effective job of modeling shapes of the quality scores for a set read length
//!
//! Disadvantages of this approach:
//!
//! * The fact that we're working with a matrix with a different first element is
//!   extremely confusing
//! * In addition to the difficulty keeping track of indexes, I'm not sure how well this will
//!   translate to Rust. May need a custom data structure. Like seed + subsequent.
//! * Assumes a fixed read length, meaning you have to extrapolate for longer read lengths.
//! * In Python, at least, this was slow, although in retrospect it didn't eat up much memory.

use crate::rng::{NeatRng, NeatRngError};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use serde_json;
use std::fmt::{Display, Formatter};
use std::{io, path::PathBuf};
use thiserror::Error;

use crate::models::lib::{model_reader, model_writer};
use crate::models::sequencing_error_model::SeqModelError;
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};

pub const QUALITY_OFFSET: usize = 33;

#[derive(Debug, Error)]
pub enum QualityModelError {
    #[error("Quality model initiation returned a distribution error: {0}")]
    DistributionError(#[from] DistributionErrors),
    #[error("Quality score creation reported an RNG Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Quality model return an IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("Serde error building default model: {0}")]
    SerdeError(#[from] serde_json::Error),
    #[error("Invalid quality model configuration: {0}")]
    InvalidConfiguration(String),
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QualityScoreModel {
    // This is the vector of the quality scores possible in this dataset. This could be a list
    // of numbers from 1-42, for example, or bins of scores, [2, 13, 27, 33] or whatever the
    // dataset uses. This list is expected to be sorted.
    pub quality_score_options: Vec<usize>,
    // True for binned scores, false for continuous
    pub binned_scores: bool,
    // The assumed read length of this dataset. The model will assume this read length and adjust
    // on a per-run basis in a deterministic way (doubling positional weight arrays)
    pub assumed_read_length: usize,
    // Weights for the first position in the read length.
    pub seed_dist: DiscreteDistribution<usize>,
    // A matrix for each subsequent position along the read length after the first. Each row is a
    // discrete distribution, keyed the previous score. For example, for possible scores 0-41,
    // inclusive, there would be 42 vectors (one for each possible previous score), each giving the
    // distribution for the current position (one weight for each of 42 scores). This is based on
    // the original design in NEAT. Previous attempts to simplify this have not been able to
    // successfully reproduce quality scores.
    pub distros_from_one: Vec<Vec<DiscreteDistribution<usize>>>,
}

impl Display for QualityScoreModel {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // Just a basic display showing the possible quality scores and read length of the model.
        // This won't be used by generate_reads, but could be used by a generate quality score model
        write!(
            f,
            "QualityScoreModel: (rl: {})\n\
            \tscores: {:?}\n\
            \tbinned? {:?}\n",
            self.assumed_read_length, self.quality_score_options, self.binned_scores,
        )
    }
}

static DATA_FILE: &[u8] = include_bytes!("model_data/default_quality_score_model.json.gz");

impl QualityScoreModel {
    // methods for QualityScoreModel objects

    // Returns Result because it builds distributions that can fail; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, QualityModelError> {
        // This generates the default quality score model based on the original NEAT default.
        // the parameters are sort of outdated now:
        //  - quality_score_options: 0 - 41
        //  - binned_scores: false
        //  - assumed_read_length: 101
        // I won't list all the NEAT calculated weights. zcat the model file if you are interested, but
        // it's a lot of data.

        // this map_err thing I had help on. I feel like this is something rust-analyzer should have known.
        let reader = GzDecoder::new(DATA_FILE);
        let data: QualityScoreModel =
            serde_json::from_reader(reader).map_err(QualityModelError::SerdeError)?;
        Ok(data)
    }

    /// we will write subutilities that use these features, eventually
    pub fn from_file(filename: &PathBuf) -> Result<Self, QualityModelError> {
        // Uses the serde_json crate to read a quality model from file
        let data: QualityScoreModel = model_reader(filename).unwrap();
        Ok(data)
    }

    #[allow(dead_code)]
    // need a test
    pub fn display(&self) -> String {
        // Some detailed formats. I thing these will be useful for quality model generation debugging.
        format!(
            "QualityScoreModel: (rl: {})\n\
            \tscores: {:?}\n\
            \tbinned? {:?}\n\
            \tseed distribution: {:?}\n\
            \tfirst weight array: {:?}",
            self.assumed_read_length,
            self.quality_score_options,
            self.binned_scores,
            self.seed_dist,
            self.distros_from_one[1][0],
        )
    }

    #[allow(dead_code)]
    // need a test
    pub fn display_it_all(&self) -> String {
        format!(
            "QualityScoreModel: (rl: {})\n\
            \tscores: {:?}\n\
            \tbinned? {:?}\n\
            \tseed distribution: {:?}\n\
            \tscore weight array: {:?}",
            self.assumed_read_length,
            self.quality_score_options,
            self.binned_scores,
            self.seed_dist,
            self.distros_from_one,
        )
    }

    pub fn generate_quality_scores(
        &self,
        length: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<usize>, SeqModelError> {
        // Generates a list of quality scores of length run_read_length using the model. If the
        // input read length differs, we do some index magic to extrapolate the model
        // run_read_length: The desired read length for the model to generate.

        // This will be the list of scores generated
        let mut score_list: Vec<usize> = Vec::with_capacity(length);
        // sample the scores list with the seed weights applied to generate the first score.
        // Samples an index based on the weights, which then selects the quality score.
        let mut seed_score = self.seed_dist.sample(rng.random()?)?;
        // We don't want the first score to be 31. That throws off the parser, because it makes an @ symbol
        if seed_score == 31 {
            seed_score -= 1
        }
        // Adding the seed score to the list. This is safe so long as the qual score max is <= 255 (currently 40)
        score_list.push(seed_score);
        // To map from one length to another, we use the algorithm found in the original NEAT 2.0,
        // adapted to rust. See function for implementation details.
        let indexes: Vec<usize> = self.quality_index_remap(length);
        // Sort of annoying, but to account for the remap, in order to get "previous score" from
        // the score list, we need to know the current index we are filling, absolutely, in cases
        // of mismatches between model read length and run read length
        // We can skip the first one, since we already generated it above. On loop 1, we will look
        // at that seed score to get our first set of weights.
        let mut current_index = 1;
        for i in indexes {
            // The weight at this index is the score weights for position i, given a
            // previous score of score_list[i-1]

            // First get the index of the previous score from the original scores list.
            // This will match the index in the score weights table that corresponds to that score.
            let score_position = self
                .quality_score_options
                .iter()
                .position(|&x| x == score_list[current_index - 1])
                .unwrap();
            // Now we have an index (in the default case 0..<4) of a vector for the position, based
            // on the previous score. We have to subtract one from the index because the first matrix
            // (position 0) corresponds to the second quality score (position 1).
            // We need to figure out a better way if we're going to do discrete.
            let index = self.distros_from_one[i - 1][score_position].sample(rng.random()?)?;
            let score = self.quality_score_options[index];
            score_list.push(score);
            current_index += 1;
        }
        Ok(score_list)
    }

    fn quality_index_remap(&self, length: usize) -> Vec<usize> {
        // Basically, this function does integer division (truncation) to fill positions
        // in a vector the length of the desired read length.
        // for example. You are mapping from read length 6 to read length 8,
        // A: [0, 1, 2, 3, 4, 5] -> B: [0, 1, 2, 3, 4, 5, 6, 7]
        // The question is, since our model (A) only has 6 columns of values, and for each of the 8
        // items in the desire output (B), we need to know which of the 5 columns from A to use.
        // So we create a vector of length 8, which tells us, at each position, which column from
        // A to use. For example, at position 1, we take (6 * 1) // 8 (where `//` denotes integer
        // division), resulting in 6//8 = 0, so for the first quality score, we select based on
        // The 0 column from the A model. At position 2, (6 * 2) // 8 = 12//8 = 1. Repeating, this
        // for each i from 0 to read length, we end up with: C: [0, 0, 1, 2, 3, 3, 4, 5]
        // Similarly, we can remap this way from larger to smaller:
        // A: [0, 1, 2, 3, 4, 5, 6, 7] -> B: [0, 1, 2, 3, 4, 5]
        // C = [(8*0)//6 = 0, (8*1)//6 = 1, (8*2)//6 = 2, (8*3)//6 = 4, (8*4)//6 = 5]
        //   = [0, 1, 2, 4, 5]
        // Advantages: should be pretty quick. Easy calculations.
        // Disadvantages: Tends to lose info from the back of the read when downsizing. Might need
        //                to check that.
        if length == self.assumed_read_length {
            (1..length).collect()
        } else {
            let mut indexes: Vec<usize> = Vec::new();
            for i in 1..length {
                let index: usize = (self.assumed_read_length * i) / length;
                // This first value(s) will always be zero when run_read_length is longer than
                // assumed read length.
                if index < 1 {
                    indexes.push(1);
                } else if index >= (self.assumed_read_length - 1) {
                    indexes.push(self.assumed_read_length - 1)
                } else {
                    indexes.push(index);
                }
            }
            indexes
        }
    }

    pub fn from_counts(
        quality_score_options: Vec<usize>,
        read_length: usize,
        seed_weights: Vec<f64>,
        trans_weights: Vec<Vec<Vec<f64>>>,
        is_binned: bool,
    ) -> Result<Self, QualityModelError> {
        // Phred+33 maps score 31 to '@', which the FASTQ writer treats specially as a seed
        // (see generate_quality_scores). A binned model that included 31 would silently break
        // the binning invariant when the seed-`@` workaround fires, so reject it here. Config
        // validation in gen_seq_error_model also catches this earlier with a friendlier message.
        if is_binned && quality_score_options.contains(&31) {
            return Err(QualityModelError::InvalidConfiguration(
                "binned quality model cannot include score 31 ('@' under Phred+33)".to_string(),
            ));
        }
        // seed_dist values are actual quality scores (generate_quality_scores pushes the sampled
        // value directly); distros_from_one values are indices into quality_score_options
        // (generate_quality_scores indexes quality_score_options[sampled_index]).
        let seed_dist = DiscreteDistribution::new(&seed_weights, &quality_score_options)?;
        let n_scores = quality_score_options.len();
        let score_indices: Vec<usize> = (0..n_scores).collect();
        let uniform = vec![1.0f64; n_scores];
        let mut distros_from_one: Vec<Vec<DiscreteDistribution<usize>>> = Vec::new();
        for pos_weights in &trans_weights {
            let mut row: Vec<DiscreteDistribution<usize>> = Vec::new();
            for prev_weights in pos_weights {
                // A prev_score that never appeared as a transition source has all-zero weights.
                // Fall back to uniform so sample() doesn't always return the lowest score.
                let weights = if prev_weights.iter().all(|&w| w == 0.0) {
                    &uniform
                } else {
                    prev_weights
                };
                row.push(DiscreteDistribution::new(weights, &score_indices)?);
            }
            distros_from_one.push(row);
        }
        Ok(QualityScoreModel {
            quality_score_options,
            binned_scores: is_binned,
            assumed_read_length: read_length,
            seed_dist,
            distros_from_one,
        })
    }

    #[allow(unused)]
    /// we will write subutilities that use these features, eventually
    fn write_to_file(&self, filename: &PathBuf) -> std::io::Result<()> {
        // Uses the serde_json crate to write out the json form of the model. This will help us
        // create base datasets from old neat data, and give us a way to write out models that are
        // generated from user data.
        model_writer(self, filename)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_model_write_read() {
        // rewrite this test to read default model, so we aren't keeping long.json around
        let temp_dir = tempfile::tempdir().unwrap();
        let mut temp_file = PathBuf::from(temp_dir.path());
        temp_file.push("test.json.gz");
        let model: QualityScoreModel = QualityScoreModel::default().unwrap();
        assert_eq!(model.assumed_read_length, 101);
        let result = model.write_to_file(&temp_file);
        assert_eq!(result.unwrap(), ());
        temp_dir.close().unwrap();
    }

    #[test]
    fn test_model_default() {
        let model = QualityScoreModel::default().unwrap();
        assert_eq!(model.assumed_read_length, 101)
    }

    #[test]
    fn test_quality_scores_short() {
        let run_read_length = 100;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model
            .generate_quality_scores(run_read_length, &mut rng)
            .unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 100);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_quality_scores_same() {
        let run_read_length = 150;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model
            .generate_quality_scores(run_read_length, &mut rng)
            .unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 150);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_quality_scores_long() {
        let run_read_length = 200;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model
            .generate_quality_scores(run_read_length, &mut rng)
            .unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 200);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_quality_scores_vast_difference() {
        let run_read_length = 2000;
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();
        let model = QualityScoreModel::default().unwrap();
        let scores = model
            .generate_quality_scores(run_read_length, &mut rng)
            .unwrap();
        assert!(!scores.is_empty());
        assert_eq!(scores.len(), 2000);
        scores
            .iter()
            .map(|x| assert!(model.quality_score_options.contains(x)))
            .collect()
    }

    #[test]
    fn test_from_counts_basic() {
        // 2 scores (Q30, Q40), read_length 3 → 2 transition positions
        let options = vec![30usize, 40usize];
        let seed_weights = vec![3.0, 1.0]; // 75% Q30 seed, 25% Q40 seed
        let trans_weights = vec![
            // position 0→1
            vec![vec![3.0, 1.0], vec![0.0, 2.0]],
            // position 1→2
            vec![vec![1.0, 1.0], vec![1.0, 0.0]],
        ];
        let model =
            QualityScoreModel::from_counts(options.clone(), 3, seed_weights, trans_weights, false)
                .unwrap();
        assert_eq!(model.quality_score_options, options);
        assert_eq!(model.assumed_read_length, 3);
        assert!(!model.binned_scores);
        assert_eq!(model.distros_from_one.len(), 2);
        assert_eq!(model.distros_from_one[0].len(), 2);
        // seed samples actual quality scores
        let mut rng = NeatRng::new_from_seed(&vec!["seed".to_string()]).unwrap();
        let s = model.seed_dist.sample(rng.random().unwrap()).unwrap();
        assert!(options.contains(&s));
        // generate_quality_scores must produce the right length
        let scores = model.generate_quality_scores(3, &mut rng).unwrap();
        assert_eq!(scores.len(), 3);
        for &sc in &scores {
            assert!(options.contains(&sc), "unexpected score {sc}");
        }
    }

    #[test]
    fn test_from_counts_zero_transition_row_uses_uniform_fallback() {
        // Q40 appears in seed but never as a prev_score in any transition.
        // Without the fix its row would be degenerate (always returns index 0 = Q30).
        // With the fix its row is uniform: CDF should be [0.5, 1.0].
        let options = vec![30usize, 40usize];
        let seed_weights = vec![1.0, 1.0];
        let trans_weights = vec![
            // position 0→1: prev=Q30 has real data; prev=Q40 row is all zeros
            vec![vec![2.0, 1.0], vec![0.0, 0.0]],
        ];
        let model =
            QualityScoreModel::from_counts(options.clone(), 2, seed_weights, trans_weights, false)
                .unwrap();
        // distros_from_one[pos=0][prev_idx=1] is the formerly-zero row for Q40
        let zero_row = &model.distros_from_one[0][1];
        let cdf = zero_row.weights().unwrap();
        assert_eq!(cdf.len(), 2);
        // Uniform over 2 values → CDF = [0.5, 1.0]
        assert!(
            (cdf[0] - 0.5).abs() < 1e-10,
            "expected uniform CDF [0.5, 1.0], got {cdf:?}"
        );
        assert!((cdf[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_binned_model_generates_only_bins() {
        // Build a binned model with NovaSeq-style bins and verify that every sampled score
        // falls within the bin set across many draws and across read-length remapping.
        let bins = vec![2usize, 12, 23, 37];
        let n = bins.len();
        let seed_weights = vec![1.0; n];
        let uniform_row = vec![1.0; n];
        let trans_weights = vec![
            vec![uniform_row.clone(); n],
            vec![uniform_row.clone(); n],
            vec![uniform_row.clone(); n],
        ];
        let model =
            QualityScoreModel::from_counts(bins.clone(), 4, seed_weights, trans_weights, true)
                .unwrap();
        assert!(model.binned_scores);
        assert_eq!(model.quality_score_options, bins);
        let mut rng = NeatRng::new_from_seed(&vec!["binned".to_string()]).unwrap();
        for _ in 0..250 {
            let scores = model.generate_quality_scores(4, &mut rng).unwrap();
            for &s in &scores {
                assert!(bins.contains(&s), "binned model emitted non-bin score {s}");
            }
        }
        // Also check the read-length-remap path emits only bins.
        for _ in 0..250 {
            let scores = model.generate_quality_scores(11, &mut rng).unwrap();
            for &s in &scores {
                assert!(
                    bins.contains(&s),
                    "binned model emitted non-bin score {s} under remap"
                );
            }
        }
    }

    #[test]
    fn test_binned_model_serialization_round_trip() {
        // Build a binned model, write to disk, read it back, and verify the binned flag
        // and the bin set survive the JSON-GZ round trip. Catches future serde renames
        // (e.g., adding #[serde(skip)] or renaming a field).
        let bins = vec![2usize, 12, 23, 37];
        let n = bins.len();
        let seed_weights = vec![1.0, 2.0, 3.0, 4.0];
        let row = vec![1.0; n];
        let trans_weights = vec![vec![row.clone(); n]; 3];
        let model =
            QualityScoreModel::from_counts(bins.clone(), 4, seed_weights, trans_weights, true)
                .unwrap();

        let temp_dir = tempfile::tempdir().unwrap();
        let mut path = PathBuf::from(temp_dir.path());
        path.push("binned_model.json.gz");
        model.write_to_file(&path).unwrap();
        let loaded = QualityScoreModel::from_file(&path).unwrap();

        assert!(
            loaded.binned_scores,
            "binned_scores flag must survive round trip"
        );
        assert_eq!(loaded.quality_score_options, bins);
        assert_eq!(loaded.assumed_read_length, model.assumed_read_length);
        assert_eq!(
            loaded.distros_from_one.len(),
            model.distros_from_one.len(),
            "transition matrix shape must survive round trip",
        );
        temp_dir.close().unwrap();
    }

    #[test]
    fn test_from_counts_rejects_binned_with_thirty_one() {
        let options = vec![2usize, 12, 31, 37];
        let n = options.len();
        let seed_weights = vec![1.0; n];
        let uniform_row = vec![1.0; n];
        let trans_weights = vec![vec![uniform_row.clone(); n]];
        let result = QualityScoreModel::from_counts(options, 2, seed_weights, trans_weights, true);
        assert!(matches!(
            result,
            Err(QualityModelError::InvalidConfiguration(_))
        ));
    }

    #[test]
    fn test_from_counts_generates_only_observed_scores() {
        // All scores produced by generate_quality_scores must be in quality_score_options.
        let options = vec![20usize, 30usize, 40usize];
        let n = options.len();
        let seed_weights = vec![1.0; n];
        let uniform_row = vec![1.0; n];
        let trans_weights = vec![
            vec![
                uniform_row.clone(),
                uniform_row.clone(),
                uniform_row.clone(),
            ],
            vec![
                uniform_row.clone(),
                uniform_row.clone(),
                uniform_row.clone(),
            ],
        ];
        let model =
            QualityScoreModel::from_counts(options.clone(), 3, seed_weights, trans_weights, false)
                .unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["t".to_string()]).unwrap();
        for _ in 0..50 {
            let scores = model.generate_quality_scores(3, &mut rng).unwrap();
            for &s in &scores {
                assert!(options.contains(&s), "score {s} not in options");
            }
        }
    }
}
