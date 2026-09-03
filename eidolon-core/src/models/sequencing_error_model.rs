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

fn default_insertion_fraction() -> f64 {
    0.4
}

/// Indel-error propensity by local homopolymer run length, indexed by `run - 1`, with the
/// last entry covering every run at or above its index.
///
/// Measured on HCC1395 matched normal, chr20/21/22 at 46x, over an exact background of
/// 3,999,990 reference bases (1,726 slippage events; Delta job 21674484). Each entry is a
/// normalized enrichment — the share of indel errors occurring at that run length divided
/// by the share of reference bases at that run length — so it is 1.0-centred **by
/// construction** over the human background it was measured on. Applying it therefore
/// redistributes [`SequencingErrorModel::indel_probability`] across sequence context
/// without changing the genome-wide total on human. On a reference with a different
/// homopolymer composition the realized total moves with that composition, which is the
/// intended behaviour: a genome with fewer homopolymers really does slip less.
///
/// This is a **shipped default, not a measurement of the user's data** — the same status
/// the fragment-length model carries. See `model_data/README.md`. Issue #662 makes it
/// fittable from a BAM.
///
/// Deliberately NOT the variant curve from #378: variants reach 60.44x at runs >= 10 where
/// errors reach 39.20x, and conflating the two is the mistake #378 already records.
pub(crate) const DEFAULT_INDEL_CONTEXT_CURVE: [f64; 10] = [
    0.64, 0.76, 0.82, 1.11, 1.58, 1.84, 5.64, 12.16, 24.24, 39.20,
];

fn default_indel_context_curve() -> Vec<f64> {
    DEFAULT_INDEL_CONTEXT_CURVE.to_vec()
}

/// How far a run must be measured for the SHIPPED curve before the answer stops mattering.
///
/// This describes [`DEFAULT_INDEL_CONTEXT_CURVE`] only. A caller must not use it to bound
/// its own scan — a model carrying a FITTED curve (#662) may have more buckets than the
/// default, and a scan capped here could never reach them: the file would hold 20 entries,
/// every value lookup would be correct, and the top buckets would simply never be asked
/// about. Use [`SequencingErrorModel::context_run_cap`], which reads the length off the
/// curve actually loaded.
pub const INDEL_CONTEXT_RUN_CAP: usize = DEFAULT_INDEL_CONTEXT_CURVE.len();

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequencingErrorModel {
    // Neat only dealt with 2 types of sequencing errors: snps and small indels.
    // We will retain that idea and assume it is accurate.
    error_rate: f64,
    del_length_distribution: DiscreteDistribution<usize>,
    ins_length_distribution: DiscreteDistribution<usize>,
    indel_probability: f64,
    #[serde(default = "default_insertion_fraction")]
    insertion_fraction: f64,
    /// Scales `indel_probability` by local homopolymer run length. A model file written
    /// before this field existed deserializes to the shipped curve rather than failing.
    #[serde(default = "default_indel_context_curve")]
    indel_context_curve: Vec<f64>,
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
        let default_indel_probability = 0.01;
        // default is no bias
        let default_insertion_bias =
            DiscreteDistribution::new(&vec![1.0, 1.0, 1.0, 1.0], &ALLOWED_NUCS.to_vec())?;
        let quality_score_model = QualityScoreModel::default()?;

        Ok(SequencingErrorModel {
            error_rate: default_error_rate,
            del_length_distribution: default_del_distr,
            ins_length_distribution: default_ins_distr,
            indel_probability: default_indel_probability,
            insertion_fraction: default_insertion_fraction(),
            indel_context_curve: default_indel_context_curve(),
            insertion_bias: default_insertion_bias,
            transition_distros: default_transition_distros,
            quality_score_model,
        })
    }

    pub fn from_file(filename: &PathBuf) -> Result<Self, SeqModelError> {
        Ok(model_reader(filename)?)
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
            indel_probability: 0.01,
            insertion_fraction: default_insertion_fraction(),
            indel_context_curve: default_indel_context_curve(),
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

    /// The chance that an error at this base is an indel rather than a substitution.
    ///
    /// `homopolymer_run` is the length of the maximal homopolymer run the base sits in.
    /// `None` — or `Some(0)`, which is not a meaningful run length — means "no context
    /// available" and yields the flat, context-free probability, so a caller that cannot
    /// supply context keeps its previous behaviour exactly.
    ///
    /// Runs longer than the curve saturate at its last entry: the measurement pooled every
    /// run of 10 or more into one bucket, so claiming a distinction beyond that would be
    /// inventing precision the data does not have.
    fn indel_probability_at(&self, homopolymer_run: Option<usize>) -> f64 {
        let scale = match homopolymer_run {
            Some(run) if run > 0 && !self.indel_context_curve.is_empty() => {
                let last = self.indel_context_curve.len() - 1;
                self.indel_context_curve[(run - 1).min(last)]
            }
            _ => 1.0,
        };
        // Clamped because the curve reaches 39.2x: any `indel_probability` above ~0.026
        // would otherwise exceed 1.0 at long runs. Unreachable at the shipped 0.01, but
        // #662 makes the base rate fittable and a fitted value has no such guarantee.
        (self.indel_probability * scale).clamp(0.0, 1.0)
    }

    /// How far a caller must measure a homopolymer run for THIS model's curve.
    ///
    /// Read off the loaded curve rather than a constant, so a fitted curve with more
    /// buckets than the shipped default still has its tail reached. Capping a scan at the
    /// default's length would make a longer curve's top entries unreachable without any
    /// error — the file would look complete and the strongest signal would be inert.
    ///
    /// At least 1: a degenerate empty curve must not ask for a zero-length scan.
    pub fn context_run_cap(&self) -> usize {
        self.indel_context_curve.len().max(1)
    }

    pub fn generate_sequencing_error(
        &self,
        reference: Nucleotide,
        homopolymer_run: Option<usize>,
        rng: &mut NeatRng,
    ) -> Result<SequencingErrorType, SeqModelError> {
        // This method picks an error type and determines any additional data needed
        // for the current error, based on the statistical model
        if rng.random()? < self.indel_probability_at(homopolymer_run) {
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
        if rng.random()? < self.insertion_fraction {
            // Insertion vs deletion is NEAT2's SIE_INS_FREQ (0.4), not an even split; the
            // "fifty-fifty" this comment used to claim was the hardcoded 0.5 that #660
            // replaced with `insertion_fraction`.
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
            .generate_sequencing_error(Nucleotide::A, None, &mut rng)
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
            .generate_sequencing_error(Nucleotide::C, None, &mut make_rng())
            .unwrap();
        let error2 = model
            .generate_sequencing_error(Nucleotide::C, None, &mut make_rng())
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
            insertion_fraction: 0.4,
            // The shipped curve, but the calls below pass no context, so it never
            // applies — indel_probability stays a flat 1.0 and the assertion holds.
            indel_context_curve: default_indel_context_curve(),
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
                .generate_sequencing_error(Nucleotide::A, None, &mut rng)
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
        let mut model = SequencingErrorModel::default().unwrap();
        // A present value must win over Serde's fallback for pre-field models.
        model.insertion_fraction = 0.25;
        model.write_model(&path).unwrap();
        let loaded = SequencingErrorModel::from_file(&path).unwrap();
        assert!((loaded.error_rate - model.error_rate).abs() < 1e-10);
        assert!((loaded.indel_probability - model.indel_probability).abs() < 1e-10);
        assert!((loaded.insertion_fraction - model.insertion_fraction).abs() < 1e-10);
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
        assert!((model.indel_probability - 0.01).abs() < 1e-15);
        assert!((model.insertion_fraction - 0.4).abs() < 1e-15);
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

    fn assert_distribution<T>(
        distribution: &DiscreteDistribution<T>,
        expected_values: Vec<T>,
        expected_cdf: &[f64],
        name: &str,
    ) where
        T: std::fmt::Debug
            + PartialEq
            + Clone
            + serde::Serialize
            + for<'de> serde::Deserialize<'de>,
    {
        assert_eq!(
            distribution.values().unwrap(),
            expected_values,
            "{name}: values drifted from NEAT2"
        );
        let actual_cdf = distribution.weights().unwrap();
        assert_eq!(
            actual_cdf.len(),
            expected_cdf.len(),
            "{name}: CDF width changed"
        );
        for (index, expected) in expected_cdf.iter().enumerate() {
            assert!(
                (actual_cdf[index] - expected).abs() < 1e-12,
                "{name}: CDF entry {index} was {}, expected {expected}",
                actual_cdf[index]
            );
        }
    }

    #[test]
    fn neat2_gen_seq_error_model_defaults_are_pinned_to_the_source_constants() {
        // NEAT2 neat2/utilities/genSeqErrorModel.py, in the `if PILEUP == None`
        // default branch, defines SIE_RATE = 0.01 (fraction of sequencing errors
        // that are indels) and SIE_INS_FREQ = 0.4 (fraction of those indels that
        // are insertions). Keep their meanings separate: their earlier conflation
        // made 40 times too many sequencing errors into indels.
        let default_model = SequencingErrorModel::default().unwrap();
        let fitted_model =
            SequencingErrorModel::from_raw_data(0.006, QualityScoreModel::default().unwrap(), None)
                .unwrap();

        for (name, model) in [
            ("default()", &default_model),
            ("from_raw_data()", &fitted_model),
        ] {
            assert!(
                (model.indel_probability - 0.01).abs() < f64::EPSILON,
                "{name}: NEAT2 genSeqErrorModel.py SIE_RATE must remain 0.01"
            );
            assert!(
                (model.insertion_fraction - 0.4).abs() < f64::EPSILON,
                "{name}: NEAT2 genSeqErrorModel.py SIE_INS_FREQ must remain 0.4"
            );

            // These NEAT2 defaults already translated correctly. Assert them here so
            // changing SIE_RATE cannot accidentally rewrite unrelated parameters.
            assert_distribution(
                &model.ins_length_distribution,
                vec![1, 2],
                &[0.999, 1.0],
                "NEAT2 sequencing insertion-length distribution [0.999, 0.001]",
            );
            assert_distribution(
                &model.del_length_distribution,
                vec![1, 2],
                &[0.999, 1.0],
                "NEAT2 sequencing deletion-length distribution [0.999, 0.001]",
            );
            assert_distribution(
                &model.insertion_bias,
                ALLOWED_NUCS.to_vec(),
                &[0.25, 0.5, 0.75, 1.0],
                "NEAT2 uniform sequencing insertion-base composition",
            );
        }

        let expected_matrix_cdf = [
            [0.0, 0.4918, 0.8295, 1.0],
            [0.5238, 0.5238, 0.7899, 1.0],
            [0.3754 / 0.9999, 0.6109 / 0.9999, 0.6109 / 0.9999, 1.0],
            [0.2505 / 0.9999, 0.5057 / 0.9999, 1.0, 1.0],
        ];
        for (name, model) in [
            ("default()", &default_model),
            ("from_raw_data()", &fitted_model),
        ] {
            for (base, expected) in ALLOWED_NUCS.iter().zip(expected_matrix_cdf) {
                assert_distribution(
                    &model.transition_distros[base],
                    ALLOWED_NUCS.to_vec(),
                    &expected,
                    &format!("{name}: NEAT2 sequencing-error transition matrix row {base:?}"),
                );
            }
        }
    }

    #[test]
    fn neat2_indel_and_insertion_fractions_drive_generated_error_types() {
        // Error type is selected after a caller has decided that a base has a sequencing
        // error at its fixed quality.  Sample only that conditional choice: Q does not
        // affect this split.  250k draws make both historical mutations (0.4 indels and
        // a 0.5 insertion split) unambiguously outside these intervals.
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["NEAT2 SIE regression".to_string()]).unwrap();
        let mut indels = 0usize;
        let mut insertions = 0usize;
        const DRAWS: usize = 250_000;

        for _ in 0..DRAWS {
            match model
                .generate_sequencing_error(Nucleotide::A, None, &mut rng)
                .unwrap()
            {
                SequencingErrorType::SnpError(_) => {}
                SequencingErrorType::InsertionError(_) => {
                    indels += 1;
                    insertions += 1;
                }
                SequencingErrorType::DeletionError(_) => indels += 1,
            }
        }

        let indel_fraction = indels as f64 / DRAWS as f64;
        let insertion_fraction = insertions as f64 / indels as f64;
        assert!(
            (0.008..0.012).contains(&indel_fraction),
            "SIE_RATE: observed {indel_fraction:.5}; expected about 0.01 indels per generated error"
        );
        assert!(
            (0.37..0.43).contains(&insertion_fraction),
            "SIE_INS_FREQ: observed {insertion_fraction:.5}; expected about 0.4 insertions per indel"
        );
    }

    #[test]
    fn the_shipped_curve_is_pinned_to_the_measured_values() {
        // Known answer, pinned against its source the way #660 pinned the NEAT2 constants.
        // These are enrichments from Delta job 21674484 (HCC1395 normal, chr20/21/22,
        // 1,726 slippage events over 3,999,990 reference bases). Changing one silently is
        // exactly how the NEAT2 mistranslation survived, so name them here.
        assert_eq!(
            DEFAULT_INDEL_CONTEXT_CURVE,
            [
                0.64, 0.76, 0.82, 1.11, 1.58, 1.84, 5.64, 12.16, 24.24, 39.20
            ],
            "the shipped indel-context curve drifted from job 21674484"
        );
        // Monotone, crossing 1.0 at run 4. Both properties are load-bearing: monotonicity
        // is the biological claim (longer run, more slippage), and the crossing point is
        // what makes this a redistribution rather than a rate increase.
        for window in DEFAULT_INDEL_CONTEXT_CURVE.windows(2) {
            assert!(
                window[1] > window[0],
                "curve must be monotone increasing; {:?} is not",
                window
            );
        }
        // Stated as "where does it cross" rather than two point checks: the crossing
        // point is the claim that this redistributes rather than adds.
        let crossing = DEFAULT_INDEL_CONTEXT_CURVE
            .iter()
            .position(|&value| value > 1.0)
            .expect("a curve that never exceeds 1.0 could only ever suppress");
        assert_eq!(
            crossing + 1,
            4,
            "curve must cross 1.0 at run 4; it crosses at run {}",
            crossing + 1
        );
    }

    #[test]
    fn the_scan_cap_follows_the_loaded_curve_not_the_shipped_default() {
        // The trap this guards: #662 fits a curve from a BAM, and a fitted curve need not
        // have the default's ten buckets. If the read generator bounded its run-length
        // scan by the DEFAULT's length, a longer curve's top entries could never be
        // reached — the model file would hold every value, each lookup would return the
        // right number, and the strongest buckets would simply never be asked about. That
        // failure is completely silent, which is why it gets an explicit test.
        let mut model = SequencingErrorModel::default().unwrap();
        assert_eq!(model.context_run_cap(), DEFAULT_INDEL_CONTEXT_CURVE.len());

        model.indel_context_curve = (1..=20).map(|i| i as f64).collect();
        assert_eq!(
            model.context_run_cap(),
            20,
            "a 20-bucket fitted curve must ask for a 20-deep scan"
        );
        // Every bucket must be distinguishable at the cap the model asks for, or the tail
        // is dead weight.
        assert_ne!(
            model.indel_probability_at(Some(model.context_run_cap())),
            model.indel_probability_at(Some(DEFAULT_INDEL_CONTEXT_CURVE.len())),
            "buckets past the default length are unreachable — the #662 trap"
        );

        // A shorter curve is safe in the other direction, but must still saturate rather
        // than index out of bounds.
        model.indel_context_curve = vec![0.5, 2.0];
        assert_eq!(model.context_run_cap(), 2);
        assert_eq!(
            model.indel_probability_at(Some(2)),
            model.indel_probability_at(Some(50))
        );

        // Degenerate: an empty curve must not request a zero-length scan.
        model.indel_context_curve = Vec::new();
        assert_eq!(model.context_run_cap(), 1, "cap must never be 0");
    }

    #[test]
    fn the_run_cap_and_the_curve_cannot_drift_apart() {
        // Cross-component invariant. The read generator caps its scan at
        // INDEL_CONTEXT_RUN_CAP; the model saturates at the curve's last entry. Neither
        // side is wrong on its own, and nothing else asserts they must agree — the exact
        // shape of defect CLAUDE.md requires be pinned rather than left to two literals
        // happening to match.
        assert_eq!(
            INDEL_CONTEXT_RUN_CAP,
            DEFAULT_INDEL_CONTEXT_CURVE.len(),
            "a scan capped short of the curve would never reach its top entries"
        );
        let model = SequencingErrorModel::default().unwrap();
        assert_eq!(
            model.indel_probability_at(Some(INDEL_CONTEXT_RUN_CAP)),
            model.indel_probability_at(Some(INDEL_CONTEXT_RUN_CAP + 500)),
            "runs past the cap must be indistinguishable, or capping the scan changes results"
        );
    }

    #[test]
    fn indel_probability_tracks_the_curve_across_every_run_length() {
        // The decision under test is the SCALING, so assert the whole shape. A single
        // fixture would pass just as happily for code that returns a constant.
        let model = SequencingErrorModel::default().unwrap();
        for (index, scale) in DEFAULT_INDEL_CONTEXT_CURVE.iter().enumerate() {
            let run = index + 1;
            let expected = 0.01 * scale; // computed from the table, not from the code
            let actual = model.indel_probability_at(Some(run));
            assert!(
                (actual - expected).abs() < 1e-12,
                "run {run}: expected {expected}, got {actual}"
            );
        }
    }

    #[test]
    fn absent_context_yields_exactly_the_flat_context_free_rate() {
        // Must-not-fire. `None` is the adapter path and every pre-#661 caller; it must
        // reproduce #660 behaviour bit for bit, not merely approximately.
        let model = SequencingErrorModel::default().unwrap();
        assert_eq!(model.indel_probability_at(None), model.indel_probability);
        // 0 is not a run length. A caller that computes a run over an N gets 0 back from
        // `homopolymer_run_at`, and that must mean "no context", not "index -1".
        assert_eq!(model.indel_probability_at(Some(0)), model.indel_probability);
        // An empty curve is a degenerate model file, not a panic and not an index error.
        let mut empty = SequencingErrorModel::default().unwrap();
        empty.indel_context_curve = Vec::new();
        assert_eq!(empty.indel_probability_at(Some(7)), empty.indel_probability);
    }

    #[test]
    fn a_fitted_base_rate_cannot_push_the_scaled_probability_past_one() {
        // #662 makes indel_probability fittable. 39.2x means any base rate above ~0.026
        // would otherwise produce a probability over 1.0, which `rng.random() < p` would
        // silently read as "always an indel".
        let mut model = SequencingErrorModel::default().unwrap();
        model.indel_probability = 0.5;
        let scaled = model.indel_probability_at(Some(INDEL_CONTEXT_RUN_CAP));
        assert!(
            (0.0..=1.0).contains(&scaled),
            "scaled probability {scaled} escaped [0, 1]"
        );
        assert_eq!(scaled, 1.0, "0.5 x 39.2 must clamp to exactly 1.0");
    }

    #[test]
    fn a_homopolymer_shifts_the_generated_error_mix_by_the_curves_factor() {
        // Behavioural counterpart to the arithmetic above: the curve must reach the
        // generated error TYPES, not merely the probability function. Run 10 (39.2x)
        // against no context, same seed, 250k draws each.
        fn indel_share(run: Option<usize>, draws: usize) -> f64 {
            let model = SequencingErrorModel::default().unwrap();
            let mut rng = NeatRng::new_from_seed(&vec!["indel context mix".to_string()]).unwrap();
            let mut indels = 0usize;
            for _ in 0..draws {
                match model
                    .generate_sequencing_error(Nucleotide::A, run, &mut rng)
                    .unwrap()
                {
                    SequencingErrorType::SnpError(_) => {}
                    _ => indels += 1,
                }
            }
            indels as f64 / draws as f64
        }
        const DRAWS: usize = 250_000;
        let flat = indel_share(None, DRAWS);
        let enriched = indel_share(Some(10), DRAWS);
        let suppressed = indel_share(Some(1), DRAWS);

        // Expected values come from the table (0.01 x 39.20 and 0.01 x 0.64), computed
        // independently of the code under test.
        assert!(
            (0.37..0.41).contains(&enriched),
            "run 10 should give about 0.392 indels per error; got {enriched:.5}"
        );
        assert!(
            (0.0055..0.0075).contains(&suppressed),
            "run 1 should give about 0.0064 indels per error; got {suppressed:.5}"
        );
        assert!(
            (0.008..0.012).contains(&flat),
            "no context must stay at #660's 0.01; got {flat:.5}"
        );
        assert!(
            enriched > flat && flat > suppressed,
            "ordering broke: run1 {suppressed:.5} < none {flat:.5} < run10 {enriched:.5}"
        );
    }

    #[test]
    fn models_without_insertion_fraction_deserialize_with_the_neat2_default() {
        use flate2::{Compression, write::GzEncoder};
        use std::fs::File;
        use std::io::Write;
        use tempfile::tempdir;

        // Model files written before insertion_fraction existed must continue to load.
        // Create such a file from current serialized data, instead of relying on a stale
        // fixture whose other fields could hide a deserialization failure.
        let mut old_format =
            serde_json::to_value(SequencingErrorModel::default().unwrap()).unwrap();
        old_format
            .as_object_mut()
            .unwrap()
            .remove("insertion_fraction");
        let dir = tempdir().unwrap();
        let path = dir.path().join("pre-insertion-fraction.json.gz");
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        encoder
            .write_all(&serde_json::to_vec(&old_format).unwrap())
            .unwrap();
        encoder.finish().unwrap();

        let model = SequencingErrorModel::from_file(&path).unwrap();
        assert!(
            (model.insertion_fraction - 0.4).abs() < f64::EPSILON,
            "models without insertion_fraction must default to NEAT2 SIE_INS_FREQ = 0.4"
        );
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
