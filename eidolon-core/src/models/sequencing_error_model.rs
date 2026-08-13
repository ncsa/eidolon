use crate::rng::{NeatRng, NeatRngError};
use crate::{
    models::{
        lib::{model_reader, model_writer},
        quality_scores::{QualityModelError, QualityScoreModel},
    },
    structs::{
        distributions::{DiscreteDistribution, DistributionErrors},
        nucleotides::{ALLOWED_NUCS, Nucleotide},
        transition_matrix::{TransitionMatrix, TransitionMatrixError},
    },
};
use serde::{Deserialize, Serialize};
use std::{io, path::PathBuf};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SeqModelError {
    #[error("Error creating sequencing error model")]
    ModelCreationError,
    #[error("Error creating transition matrix: {0}")]
    TransMatrixError(#[from] TransitionMatrixError),
    #[error("Error sampling distribution: {0}")]
    DistributionError(#[from] DistributionErrors),
    #[error("Error with rng: {0}")]
    RngError(#[from] NeatRngError),
    #[error("No RNG supplied for this model.")]
    MissingRngError,
    #[error("Sequencing Error model return an IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("Error initializing Quality Score model: {0}")]
    QualModelError(#[from] QualityModelError),
}

#[derive(Debug)]
pub enum SequencingErrorType {
    SnpError(Nucleotide),
    InsertionError(Vec<Nucleotide>),
    DeletionError(usize),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequencingErrorModel {
    // Neat only dealt with 2 types of sequencing errors: snps and small indels.
    // We will retain that idea and assume it is accurate.
    error_rate: f64,
    del_length_distribution: DiscreteDistribution<usize>,
    ins_length_distribution: DiscreteDistribution<usize>,
    indel_probability: f64,
    insertion_bias: DiscreteDistribution<Nucleotide>,
    transition_distros: TransitionMatrix,
    quality_score_model: QualityScoreModel,
}

impl SequencingErrorModel {
    // Returns Result because it builds distributions that can fail; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, SeqModelError> {
        // This is the default sequencing error model employed by NEAT2
        // Note that this was originally in a file, and we could have done it the way we did
        // the other defaults, but it was so small, I just included it in full here.
        let default_transition_distros = TransitionMatrix::from(
            [0.0, 0.4918, 0.3377, 0.1705],
            [0.5238, 0.0, 0.2661, 0.2101],
            [0.3754, 0.2355, 0.0, 0.389],
            [0.2505, 0.2552, 0.4942, 0.0],
        )?;
        let default_error_rate = 0.006638164688495656;
        let default_lengths = vec![1, 2];
        let default_ins_distr = DiscreteDistribution::new(&vec![0.999, 0.001], &default_lengths)?;
        let default_del_distr = default_ins_distr.clone();
        let default_indel_probability = 0.4;
        // default is no bias
        let default_insertion_bias =
            DiscreteDistribution::new(&vec![1.0, 1.0, 1.0, 1.0], &ALLOWED_NUCS.to_vec())?;
        let quality_score_model = QualityScoreModel::default()?;

        Ok(SequencingErrorModel {
            error_rate: default_error_rate,
            del_length_distribution: default_del_distr,
            ins_length_distribution: default_ins_distr,
            indel_probability: default_indel_probability,
            insertion_bias: default_insertion_bias,
            transition_distros: default_transition_distros,
            quality_score_model,
        })
    }

    pub fn from_file(filename: &PathBuf) -> Result<Self, SeqModelError> {
        let data: SequencingErrorModel = model_reader(filename).unwrap();
        Ok(data)
    }

    pub fn from_raw_data(
        error_rate: f64,
        quality_score_model: QualityScoreModel,
        transition_matrix: Option<TransitionMatrix>,
    ) -> Result<Self, SeqModelError> {
        let transition_distros = match transition_matrix {
            Some(tm) => tm,
            None => TransitionMatrix::from(
                [0.0, 0.4918, 0.3377, 0.1705],
                [0.5238, 0.0, 0.2661, 0.2101],
                [0.3754, 0.2355, 0.0, 0.389],
                [0.2505, 0.2552, 0.4942, 0.0],
            )?,
        };
        let default_lengths = vec![1, 2];
        let default_ins_distr = DiscreteDistribution::new(&vec![0.999, 0.001], &default_lengths)?;
        let default_del_distr = default_ins_distr.clone();
        Ok(SequencingErrorModel {
            error_rate,
            del_length_distribution: default_del_distr,
            ins_length_distribution: default_ins_distr,
            indel_probability: 0.4,
            insertion_bias: DiscreteDistribution::new(
                &vec![1.0, 1.0, 1.0, 1.0],
                &ALLOWED_NUCS.to_vec(),
            )?,
            transition_distros,
            quality_score_model,
        })
    }

    pub fn error_rate(&self) -> f64 {
        self.error_rate
    }

    pub fn write_model(&self, filename: &PathBuf) -> Result<(), SeqModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    pub fn generate_sequencing_error(
        &self,
        reference: Nucleotide,
        rng: &mut NeatRng,
    ) -> Result<SequencingErrorType, SeqModelError> {
        // This method picks an error type and determines any additional data needed
        // for the current error, based on the statistical model
        if rng.random()? < self.indel_probability {
            // Indel error
            Ok(self.generate_indel_error(rng)?)
        } else {
            // SNP error
            Ok(SequencingErrorType::SnpError(
                self.generate_snp_error(reference, rng.random()?)?,
            ))
        }
    }

    pub fn convert_score(&self, score: usize) -> Result<f64, SeqModelError> {
        // Takes a quality score, converts it to a probability of error, and returns the result
        let score = score as f64;
        Ok(10.0_f64.powf(-score / 10.0))
    }

    fn generate_snp_error(
        &self,
        reference: Nucleotide,
        rand: f64,
    ) -> Result<Nucleotide, SeqModelError> {
        // This is a basic mutation function for starting us off
        // Pick the weights list for the base that was input
        // We will use this simple model for sequence errors ultimately.
        let distro = &self.transition_distros[&reference];
        // Now we create a distribution from the weights and sample our choices.
        // We have constructed things such that this will return a valid u8. But
        // to be extra safe, we could mod by 4 and then convert
        Ok(distro.sample(rand)?)
    }

    fn generate_indel_error(
        &self,
        rng: &mut NeatRng,
    ) -> Result<SequencingErrorType, SeqModelError> {
        // Returns either an insertion (option 1) or a deletion (option 2) depending on a random selection from a list of potential
        // error lengths (-2..2). This makes an insertion of up to 2 bases as likely as a random deletion of up to 2 bases.
        if rng.random()? < 0.5 {
            // We assume fifty-fifty chance of insertion v deletion
            // insertion
            let mut sequence = Vec::new();
            let length = self.ins_length_distribution.sample(rng.random()?)?;
            for _ in 0..length {
                // We could mod this value by 4 to ensure it is a valid base. Or create a data structure.
                sequence.push(self.insertion_bias.sample(rng.random()?)?)
            }
            // Insertion of sequence
            Ok(SequencingErrorType::InsertionError(sequence))
        } else {
            // Deletion
            let length = self.del_length_distribution.sample(rng.random()?)?;
            Ok(SequencingErrorType::DeletionError(length))
        }
    }

    pub fn generate_quality_scores(
        &self,
        read_length: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<usize>, SeqModelError> {
        self.quality_score_model
            .generate_quality_scores(read_length, rng)
    }

    /// Borrow the inner quality-score model. Useful for tests and for callers that need
    /// to inspect model metadata (e.g., `binned_scores`, `quality_score_options`).
    pub fn quality_score_model(&self) -> &QualityScoreModel {
        &self.quality_score_model
    }

    /// Borrow the SNP transition matrix. Exists so a test can assert which matrix a
    /// built model actually carries — a BAM-inferred one, a TSV-supplied one, or the
    /// default. Without this the only reachable assertion was "the file exists and
    /// deserializes", which passes just as happily when the matrix is wrong.
    pub fn transition_distros(&self) -> &TransitionMatrix {
        &self.transition_distros
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::NeatRng;

    fn make_rng() -> NeatRng {
        NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap()
    }

    #[test]
    fn test_sequencing_error_model() {
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = make_rng();
        let result = model
            .generate_sequencing_error(Nucleotide::A, &mut rng)
            .unwrap();
        match result {
            SequencingErrorType::SnpError(base) => assert_ne!(base, Nucleotide::A),
            SequencingErrorType::InsertionError(seq) => assert!(!seq.is_empty()),
            SequencingErrorType::DeletionError(len) => assert!(len > 0),
        }
    }

    #[test]
    fn test_convert_score() {
        let model = SequencingErrorModel::default().unwrap();
        // Q20 → error prob 0.01
        assert!((model.convert_score(20).unwrap() - 0.01).abs() < 1e-10);
        // Q30 → error prob 0.001
        assert!((model.convert_score(30).unwrap() - 0.001).abs() < 1e-10);
        // Q0 → error prob 1.0
        assert!((model.convert_score(0).unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sequencing_error_deterministic() {
        let model = SequencingErrorModel::default().unwrap();
        let error1 = model
            .generate_sequencing_error(Nucleotide::C, &mut make_rng())
            .unwrap();
        let error2 = model
            .generate_sequencing_error(Nucleotide::C, &mut make_rng())
            .unwrap();
        let type1 = match error1 {
            SequencingErrorType::SnpError(_) => 0,
            SequencingErrorType::InsertionError(_) => 1,
            SequencingErrorType::DeletionError(_) => 2,
        };
        let type2 = match error2 {
            SequencingErrorType::SnpError(_) => 0,
            SequencingErrorType::InsertionError(_) => 1,
            SequencingErrorType::DeletionError(_) => 2,
        };
        assert_eq!(type1, type2);
    }

    #[test]
    fn test_indel_forced_path_types() {
        // With indel_probability=1.0 every call must produce an insertion or deletion, never a SNP.
        let model = SequencingErrorModel {
            error_rate: 0.1,
            del_length_distribution: DiscreteDistribution::new(&vec![1.0], &vec![1]).unwrap(),
            ins_length_distribution: DiscreteDistribution::new(&vec![1.0], &vec![1]).unwrap(),
            indel_probability: 1.0,
            insertion_bias: DiscreteDistribution::new(
                &vec![1.0, 1.0, 1.0, 1.0],
                &Vec::from(ALLOWED_NUCS),
            )
            .unwrap(),
            transition_distros: TransitionMatrix::from(
                [0.0, 0.5, 0.25, 0.25],
                [0.5, 0.0, 0.25, 0.25],
                [0.25, 0.25, 0.0, 0.5],
                [0.25, 0.25, 0.5, 0.0],
            )
            .unwrap(),
            quality_score_model: QualityScoreModel::default().unwrap(),
        };
        let mut rng = make_rng();
        let mut saw_insertion = false;
        let mut saw_deletion = false;
        for _ in 0..20 {
            match model
                .generate_sequencing_error(Nucleotide::A, &mut rng)
                .unwrap()
            {
                SequencingErrorType::SnpError(_) => {
                    panic!("should not produce SNP when indel_probability=1.0")
                }
                SequencingErrorType::InsertionError(seq) => {
                    assert!(!seq.is_empty());
                    saw_insertion = true;
                }
                SequencingErrorType::DeletionError(len) => {
                    assert!(len > 0);
                    saw_deletion = true;
                }
            }
            if saw_insertion && saw_deletion {
                break;
            }
        }
        assert!(
            saw_insertion,
            "should have seen at least one insertion in 20 calls"
        );
        assert!(
            saw_deletion,
            "should have seen at least one deletion in 20 calls"
        );
    }

    #[test]
    fn test_sequencing_error_model_file_round_trip() {
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let path = dir.path().join("seq_error_model.json.gz");
        let model = SequencingErrorModel::default().unwrap();
        model.write_model(&path).unwrap();
        let loaded = SequencingErrorModel::from_file(&path).unwrap();
        assert!((loaded.error_rate - model.error_rate).abs() < 1e-10);
        assert!((loaded.indel_probability - model.indel_probability).abs() < 1e-10);
    }

    #[test]
    fn test_convert_score_additional() {
        let model = SequencingErrorModel::default().unwrap();
        // Q40 → 10^(-40/10) = 0.0001
        assert!((model.convert_score(40).unwrap() - 0.0001).abs() < 1e-12);
        // Q10 → 10^(-10/10) = 0.1
        assert!((model.convert_score(10).unwrap() - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_from_raw_data_stores_error_rate_and_defaults() {
        use crate::models::quality_scores::QualityScoreModel;
        let quality_score_model = QualityScoreModel::default().unwrap();
        let error_rate = 0.00312;
        let model =
            SequencingErrorModel::from_raw_data(error_rate, quality_score_model, None).unwrap();
        assert!((model.error_rate() - error_rate).abs() < 1e-15);
        // indel_probability default from NEAT2 should be preserved
        assert!((model.indel_probability - 0.4).abs() < 1e-15);
        // Model must be usable
        let mut rng = NeatRng::new_from_seed(&vec!["r".to_string()]).unwrap();
        let scores = model.generate_quality_scores(100, &mut rng).unwrap();
        assert_eq!(scores.len(), 100);
    }

    #[test]
    fn test_from_raw_data_round_trips_file() {
        use crate::models::quality_scores::QualityScoreModel;
        use tempfile::tempdir;
        let dir = tempdir().unwrap();
        let path = dir.path().join("from_raw.json.gz");
        let model = SequencingErrorModel::from_raw_data(
            0.00555,
            QualityScoreModel::default().unwrap(),
            None,
        )
        .unwrap();
        model.write_model(&path).unwrap();
        let loaded = SequencingErrorModel::from_file(&path).unwrap();
        assert!((loaded.error_rate() - 0.00555).abs() < 1e-10);
    }

    #[test]
    fn test_sequencing_error_model_binned_emits_only_bins() {
        // Wrap a binned QualityScoreModel in SequencingErrorModel, round-trip through disk,
        // and sample via the wrapper (not the inner model directly). Catches regressions in
        // the wrapper's delegation path and confirms the binned flag survives the
        // SequencingErrorModel serialization layer too.
        use crate::models::quality_scores::QualityScoreModel;
        use std::collections::HashSet;
        use tempfile::tempdir;

        let bins = vec![2usize, 12, 23, 37];
        let n = bins.len();
        let row = vec![1.0; n];
        let trans_weights = vec![vec![row.clone(); n]; 3];
        let qsm =
            QualityScoreModel::from_counts(bins.clone(), 4, vec![1.0; n], trans_weights, true)
                .unwrap();
        let model = SequencingErrorModel::from_raw_data(0.001, qsm, None).unwrap();

        let dir = tempdir().unwrap();
        let path = dir.path().join("binned_seq_err.json.gz");
        model.write_model(&path).unwrap();
        let loaded = SequencingErrorModel::from_file(&path).unwrap();
        assert!(loaded.quality_score_model().binned_scores);

        let bin_set: HashSet<usize> = bins.iter().copied().collect();
        let mut rng = make_rng();
        for _ in 0..200 {
            let scores = loaded.generate_quality_scores(50, &mut rng).unwrap();
            for &s in &scores {
                assert!(
                    bin_set.contains(&s),
                    "wrapper emitted non-bin score {s}; bins={bins:?}"
                );
            }
        }
    }
}
