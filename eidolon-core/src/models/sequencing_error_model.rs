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
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::{io, path::PathBuf};
use thiserror::Error;

/// The shipped default model, fitted from GIAB HG002 2x250. See `SequencingErrorModel::default`.
pub(crate) static DEFAULT_MODEL_FILE: &[u8] =
    include_bytes!("model_data/default_sequencing_error_model.json.gz");

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
    #[error("Error reading the shipped default sequencing error model: {0}")]
    SerdeError(serde_json::Error),
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
/// by the share of reference bases at that run length — so it is 1.0-centered **by
/// construction** over the human background it was measured on. Applying it therefore
/// redistributes [`SequencingErrorModel::indel_probability`] across sequence context
/// without changing the genome-wide total on human. On a reference with a different
/// homopolymer composition the realized total moves with that composition, which is the
/// intended behavior: a genome with fewer homopolymers really does slip less.
///
/// This is a **shipped default, not a measurement of the user's data** — the same status
/// the fragment-length model carries. See `model_data/README.md`. Issue #662 makes it
/// fittable from a BAM.
///
/// This is the sequencing-error curve. Variants have their own, measurably steeper, curve
/// (60.44x at runs >= 10 against 39.20x here); it belongs to placement — see #378.
pub(crate) const DEFAULT_INDEL_CONTEXT_CURVE: [f64; 10] = [
    0.64, 0.76, 0.82, 1.11, 1.58, 1.84, 5.64, 12.16, 24.24, 39.20,
];

fn default_indel_context_curve() -> Vec<f64> {
    DEFAULT_INDEL_CONTEXT_CURVE.to_vec()
}

/// Sequencing-error indel lengths, and the observed counts behind them.
///
/// Measured on HCC1395 matched normal, GRCh38 chr20/21/22, over the ten 400 kb loci a
/// realism-panel run placed (Delta job 21801707). Indels were separated from variants by
/// support fraction — below 10% of local depth is slippage, at or above 25% is a variant —
/// and these are the **low-support** counts. The variant length distribution is a different
/// population and belongs to placement, not here.
///
/// **The weights are raw event counts, not a normalized pmf.** `DiscreteDistribution::new`
/// divides by their sum, so what appears here is the measurement as reported rather than
/// arithmetic performed on it, and it can be checked against the job output by eye.
///
/// **Not truncated.** Every observed length is carried, out to 45 bp for deletions and 30 bp
/// for insertions. An earlier version cut this at 10 bp on the grounds that the tail bins
/// hold single observations and that "low support" means below 10% of local depth rather
/// than "proven sequencing error". Both points are true and neither justifies dropping the
/// data: sparse bins make the *shape within* the tail uncertain, not its existence, and its
/// total mass (12 of 1726 events, 0.70%) is an ordinary estimate. Cutting it also removes
/// the only part of this distribution that could ever produce a candidate breakpoint, which
/// decides #672 by construction instead of measuring it.
///
/// Consequence to be aware of: the model can now emit a 45 bp deletion as a sequencing
/// error, at p = 0.00095. If that turns out to be a mapping artifact rather than slippage,
/// the fix is a better classifier upstream, not a cutoff chosen here.
///
/// What this replaced: `[0.999, 0.001]` over lengths `[1, 2]`, inherited from NEAT2 and
/// never measured. That put 99.9% of indel errors at a single base and could not emit
/// anything above 2 bp at all, against a measured 16.4% of slippage events at 3 bp or more.
/// Deletion lengths observed, in order. Insertions have their own set: the two arms did
/// not observe the same lengths, so they cannot share one array.
pub(crate) const DEL_ERROR_LENGTHS: [usize; 26] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 27, 34, 38, 45,
];

/// Insertion lengths observed, in order.
pub(crate) const INS_ERROR_LENGTHS: [usize; 19] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 18, 19, 22, 27, 30,
];

/// Observed deletion-error counts, aligned to `DEL_ERROR_LENGTHS` (n = 1058).
pub(crate) const DEL_ERROR_LENGTH_COUNTS: [f64; 26] = [
    762.0, 135.0, 42.0, 47.0, 11.0, 8.0, 5.0, 9.0, 5.0, 6.0, 2.0, 7.0, 2.0, 2.0, 3.0, 1.0, 2.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
];

/// Observed insertion-error counts, aligned to `INS_ERROR_LENGTHS` (n = 668).
///
/// Deletions and insertions are NOT the same shape — 74.0% of deletions are 1 bp against
/// 64.8% of insertions — so these are separate distributions. They were previously one
/// distribution cloned twice.
pub(crate) const INS_ERROR_LENGTH_COUNTS: [f64; 19] = [
    426.0, 120.0, 22.0, 55.0, 13.0, 8.0, 1.0, 6.0, 1.0, 5.0, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0,
    1.0, 1.0,
];

fn indel_error_length_distributions()
-> Result<(DiscreteDistribution<usize>, DiscreteDistribution<usize>), SeqModelError> {
    Ok((
        DiscreteDistribution::new(
            &INS_ERROR_LENGTH_COUNTS.to_vec(),
            &INS_ERROR_LENGTHS.to_vec(),
        )?,
        DiscreteDistribution::new(
            &DEL_ERROR_LENGTH_COUNTS.to_vec(),
            &DEL_ERROR_LENGTHS.to_vec(),
        )?,
    ))
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
    /// A fitted SUMMARY of the quality histogram, and NOT read during generation (#724).
    ///
    /// The fitter computes it as `sum(10^(-q/10) * count[q]) / total_bases` and records it.
    /// Nothing in `gen-reads` consumes it: errors are injected per base from that base's own
    /// quality score via `convert_score`, so the quality model IS the error rate and this is a
    /// description of it. Three models differing only in this field — including 0.0 and 0.5 —
    /// generate byte-identical FASTQ.
    ///
    /// It is kept because it is good provenance: it is what makes a fitted model's rate
    /// comparable against the shipped default's, which is the whole of #695. The name reads
    /// like a knob, which is how a harness came to pin it expecting the output to move, so:
    /// **to change how many errors a run produces, change the quality model.**
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
    /// The R2 mate's quality model, when the fit was given both mates (#723).
    ///
    /// `None` means one model serves both mates, which is every model written before this
    /// field existed. `skip_serializing_if` so such a model writes the SAME BYTES it wrote
    /// before -- the same requirement, for the same reason, as `QualityScoreModel::degradation`
    /// (#694): without it every model file gains `"quality_score_model_r2": null` and the
    /// frozen baselines all move.
    ///
    /// Both mates are indexed against ONE `quality_score_options` list. A model whose two
    /// halves carried their own option sets would index the wrong scores at generation, which
    /// is why the fitter pools its counts rather than fitting the mates independently.
    #[serde(skip_serializing_if = "Option::is_none")]
    quality_score_model_r2: Option<QualityScoreModel>,
}

impl SequencingErrorModel {
    // Returns Result because it builds distributions that can fail; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    /// The shipped default, deserialized whole from a model FITTED ON REAL DATA.
    ///
    /// Until v3.4.0 this was assembled from constants plus a bare quality model whose
    /// originating sample was never recorded (NEAT2's `errorModel_toy.p`). It is now one
    /// file, fitted from GIAB HG002 `NIST_Illumina_2x250bps` — a named, public, downloadable
    /// library. Provenance, including what it is NOT, is in `model_data/README.md`.
    ///
    /// WHY THE WHOLE MODEL AND NOT JUST THE QUALITY PART. The fit carries two mates (#723)
    /// and a degraded population per mate (#694), and only the first of those fits inside a
    /// bare `QualityScoreModel`: `quality_score_model_r2` is a field out here on the outer
    /// struct. Loading the file whole is what lets the default use both.
    ///
    /// The non-quality fields in the file are the same values this function used to hardcode
    /// — the fit ran without `bam_file` or `transition_matrix_file`, so it recorded the
    /// default substitution matrix, `indel_probability` 0.01, `insertion_fraction` 0.4 and
    /// the shipped indel-context curve. Swapping to the file therefore changes the quality
    /// model and `error_rate`, and nothing else.
    ///
    /// The provenance stamp (`_eidolon`) is ignored rather than checked: this file ships
    /// inside the binary that reads it, so the two cannot disagree.
    pub fn default() -> Result<Self, SeqModelError> {
        let reader = GzDecoder::new(DEFAULT_MODEL_FILE);
        serde_json::from_reader(reader).map_err(SeqModelError::SerdeError)
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
        let (default_ins_distr, default_del_distr) = indel_error_length_distributions()?;
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
            quality_score_model_r2: None,
        })
    }

    /// Attach the R2 mate's quality model (#723).
    ///
    /// The caller must have built it against the SAME `quality_score_options` as the R1 model;
    /// that is checked here rather than trusted, because the failure it prevents is silent —
    /// generation indexes one option list, so a mismatched pair emits plausible scores drawn
    /// from the wrong distribution.
    pub fn with_mate_r2(mut self, r2: QualityScoreModel) -> Result<Self, SeqModelError> {
        if r2.quality_score_options != self.quality_score_model.quality_score_options {
            return Err(SeqModelError::QualModelError(
                QualityModelError::InvalidConfiguration(format!(
                    "the two mates carry different quality score option sets ({} values for R1, \
                     {} for R2); generation indexes one list, so they must share it",
                    self.quality_score_model.quality_score_options.len(),
                    r2.quality_score_options.len()
                )),
            ));
        }
        if r2.assumed_read_length != self.quality_score_model.assumed_read_length {
            return Err(SeqModelError::QualModelError(
                QualityModelError::InvalidConfiguration(format!(
                    "the two mates were fitted at different read lengths ({} bp for R1, {} bp \
                     for R2); both describe the same run",
                    self.quality_score_model.assumed_read_length, r2.assumed_read_length
                )),
            ));
        }
        self.quality_score_model_r2 = Some(r2);
        Ok(self)
    }

    /// The R2 quality model, when this model carries one.
    pub fn quality_score_model_r2(&self) -> Option<&QualityScoreModel> {
        self.quality_score_model_r2.as_ref()
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
    /// supply context keeps its previous behavior exactly.
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
            // Insertions are 0.4 of indel errors, not an even split (#660). Measured at
            // 0.387 on real data; see model_data/README.md.
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
            quality_score_model_r2: None,
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

    /// WHAT SHIPS, pinned.
    ///
    /// The default is now a fitted file rather than an assembly of constants, so what it
    /// carries is a fact about that file. This test exists because the thing it replaced was
    /// a model whose originating sample nobody recorded — the whole point of the change is
    /// that the default is knowable, and a silent swap would undo that without failing
    /// anything else.
    ///
    /// The numbers are from the fit: GIAB HG002 `NIST_Illumina_2x250bps`, job 22233888.
    #[test]
    fn the_shipped_default_is_the_hg002_fit() {
        let model = SequencingErrorModel::default().unwrap();
        let q = model.quality_score_model();

        assert_eq!(q.assumed_read_length, 250, "fitted from a 2x250 library");
        assert_eq!(q.quality_score_options.len(), 31);
        assert_eq!(q.quality_score_options[0], 2);
        assert_eq!(*q.quality_score_options.last().unwrap(), 40);
        assert!(
            !q.binned_scores,
            "HiSeq 2500 era: continuous scoring, not binned. A binned default is #730."
        );

        let deg = q
            .degradation
            .as_ref()
            .expect("the default carries a degraded population (#694)");
        assert!(
            (deg.read_fraction - 0.1102).abs() < 5e-4,
            "R1 degraded fraction {} is not the fitted 0.1102",
            deg.read_fraction
        );

        let r2 = model
            .quality_score_model_r2()
            .expect("the default carries a separate R2 mate (#723)");
        assert_eq!(r2.assumed_read_length, 250);
        let deg2 = r2
            .degradation
            .as_ref()
            .expect("R2 carries its own degraded population");
        assert!(
            (deg2.read_fraction - 0.2598).abs() < 5e-4,
            "R2 degraded fraction {} is not the fitted 0.2598",
            deg2.read_fraction
        );
        // The measured asymmetry is the reason per-mate exists at all. If these two ever come
        // out equal, the default has lost its R2 half and both mates are drawing from one fit.
        assert!(
            deg2.read_fraction > deg.read_fraction * 2.0,
            "R2 should be materially worse than R1: {} vs {}",
            deg2.read_fraction,
            deg.read_fraction
        );

        assert!(
            (model.error_rate() - 0.003774).abs() < 1e-5,
            "error_rate {} is not the fitted 0.003774",
            model.error_rate()
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
        // indel_probability must keep its inherited value
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
            "{name}: values drifted from their source"
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
        // Pins the two inherited constants to their source values: SIE_RATE = 0.01 (the
        // fraction of sequencing errors that are indels) and SIE_INS_FREQ = 0.4 (the
        // fraction of those that are insertions), from NEAT2's genSeqErrorModel.py. The
        // two were transposed in the Rust port; this keeps their meanings separate.
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
                "{name}: SIE_RATE must remain 0.01"
            );
            assert!(
                (model.insertion_fraction - 0.4).abs() < f64::EPSILON,
                "{name}: SIE_INS_FREQ must remain 0.4"
            );

            // These translated correctly. Asserted here so changing SIE_RATE cannot
            // rewrite unrelated parameters.
            //
            // The length distributions are asserted separately, in
            // `indel_error_lengths_are_pinned_to_the_measured_distribution`: they are a
            // measurement now (job 21801707) rather than an inherited default.
            assert_distribution(
                &model.insertion_bias,
                ALLOWED_NUCS.to_vec(),
                &[0.25, 0.5, 0.75, 1.0],
                "uniform sequencing insertion-base composition",
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
                    &format!("{name}: sequencing-error transition matrix row {base:?}"),
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
        let mut rng = NeatRng::new_from_seed(&vec!["SIE regression".to_string()]).unwrap();
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
    fn indel_error_lengths_are_pinned_to_the_measured_distribution() {
        // Known answer against the source measurement, the way #660 pinned the NEAT2
        // constants. These are the low-support (slippage) counts from Delta job 21801707,
        // truncated at 10 bp; they can be read straight off that job's length table.
        assert_eq!(
            DEL_ERROR_LENGTHS,
            [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 27,
                34, 38, 45
            ],
            "deletion lengths drifted from job 21801707"
        );
        assert_eq!(
            DEL_ERROR_LENGTH_COUNTS,
            [
                762.0, 135.0, 42.0, 47.0, 11.0, 8.0, 5.0, 9.0, 5.0, 6.0, 2.0, 7.0, 2.0, 2.0, 3.0,
                1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
            ],
            "deletion-length counts drifted from job 21801707"
        );
        assert_eq!(
            INS_ERROR_LENGTHS,
            [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 18, 19, 22, 27, 30
            ],
            "insertion lengths drifted from job 21801707"
        );
        assert_eq!(
            INS_ERROR_LENGTH_COUNTS,
            [
                426.0, 120.0, 22.0, 55.0, 13.0, 8.0, 1.0, 6.0, 1.0, 5.0, 1.0, 2.0, 1.0, 2.0, 1.0,
                1.0, 1.0, 1.0, 1.0
            ],
            "insertion-length counts drifted from job 21801707"
        );
        // Totals must match the job's reported class size, or a bin was dropped.
        assert_eq!(DEL_ERROR_LENGTH_COUNTS.iter().sum::<f64>(), 1058.0);
        assert_eq!(INS_ERROR_LENGTH_COUNTS.iter().sum::<f64>(), 668.0);
        assert_eq!(DEL_ERROR_LENGTHS.len(), DEL_ERROR_LENGTH_COUNTS.len());
        assert_eq!(INS_ERROR_LENGTHS.len(), INS_ERROR_LENGTH_COUNTS.len());

        // The counts must reach the model as a normalized distribution. Expected CDFs are
        // computed here from the counts, independently of DiscreteDistribution's own
        // arithmetic, so a normalization bug cannot pass by agreeing with itself.
        let model = SequencingErrorModel::default().unwrap();
        let check = |name: &str, counts: &[f64], lengths: Vec<usize>, distro| {
            let total: f64 = counts.iter().sum();
            let mut running = 0.0;
            let expected: Vec<f64> = counts
                .iter()
                .map(|c| {
                    running += c / total;
                    running
                })
                .collect();
            assert_distribution(
                distro,
                lengths,
                &expected,
                &format!("measured {name}-length distribution"),
            );
        };
        check(
            "deletion",
            &DEL_ERROR_LENGTH_COUNTS,
            DEL_ERROR_LENGTHS.to_vec(),
            &model.del_length_distribution,
        );
        check(
            "insertion",
            &INS_ERROR_LENGTH_COUNTS,
            INS_ERROR_LENGTHS.to_vec(),
            &model.ins_length_distribution,
        );
    }

    #[test]
    fn insertions_and_deletions_have_different_measured_shapes() {
        // They were one distribution cloned twice. The measurement says they differ:
        // 74.0% of deletions are 1 bp against 64.8% of insertions. A test that only
        // checked "both are non-empty" would pass on the cloned version.
        let model = SequencingErrorModel::default().unwrap();
        let del = model.del_length_distribution.weights().unwrap();
        let ins = model.ins_length_distribution.weights().unwrap();
        assert_ne!(
            del, ins,
            "deletion and insertion length distributions must not be the same object"
        );
        // Computed from the counts, not read off the code under test.
        let d1 = 762.0 / DEL_ERROR_LENGTH_COUNTS.iter().sum::<f64>();
        let i1 = 426.0 / INS_ERROR_LENGTH_COUNTS.iter().sum::<f64>();
        assert!(
            (del[0] - d1).abs() < 1e-12,
            "deletion P(1 bp) is {}",
            del[0]
        );
        assert!(
            (ins[0] - i1).abs() < 1e-12,
            "insertion P(1 bp) is {}",
            ins[0]
        );
        assert!(
            d1 > i1,
            "deletions are more concentrated at 1 bp than insertions: {d1} vs {i1}"
        );
    }

    #[test]
    fn the_model_can_emit_indel_errors_longer_than_two_bases() {
        // The defect this replaced. NEAT2's [0.999, 0.001] over [1, 2] could not produce a
        // 3 bp indel error at all, against a measured 16.4% of slippage events at >= 3 bp.
        // Sampling the distribution directly is the honest test: going through
        // generate_sequencing_error would need ~100 draws per indel at a 1% indel rate.
        let model = SequencingErrorModel::default().unwrap();
        let mut rng = NeatRng::new_from_seed(&vec!["indel length spread".to_string()]).unwrap();
        let mut seen_del = std::collections::BTreeSet::new();
        let mut seen_ins = std::collections::BTreeSet::new();
        for _ in 0..20_000 {
            seen_del.insert(
                model
                    .del_length_distribution
                    .sample(rng.random().unwrap())
                    .unwrap(),
            );
            seen_ins.insert(
                model
                    .ins_length_distribution
                    .sample(rng.random().unwrap())
                    .unwrap(),
            );
        }
        assert!(
            seen_del.iter().any(|&l| l >= 3) && seen_ins.iter().any(|&l| l >= 3),
            "no indel error longer than 2 bp was ever drawn: del {seen_del:?} ins {seen_ins:?}"
        );
        // The tail is the part that can produce a candidate breakpoint (#672): a >= 20 bp
        // clip needs a >= 20 bp indel. Truncating it away would decide that by
        // construction, so assert it is reachable.
        assert!(
            seen_del.iter().any(|&l| l >= 20),
            "no deletion >= 20 bp was ever drawn in 20k samples; drawn: {seen_del:?}"
        );
        // Must not fire: only lengths that were actually OBSERVED, never an interpolation.
        // DiscreteDistribution samples its value list, so a 14 bp deletion appearing would
        // mean something invented a length nobody measured.
        assert!(
            seen_del.iter().all(|l| DEL_ERROR_LENGTHS.contains(l)),
            "a deletion length outside the measured set was drawn: {seen_del:?}"
        );
        assert!(
            seen_ins.iter().all(|l| INS_ERROR_LENGTHS.contains(l)),
            "an insertion length outside the measured set was drawn: {seen_ins:?}"
        );
    }

    #[test]
    fn the_shipped_curve_is_pinned_to_the_measured_values() {
        // Known answer, pinned against its source: enrichments from Delta job 21674484
        // (HCC1395 normal, chr20/21/22, 1,726 slippage events over 3,999,990 reference
        // bases). Named here so a silent change fails.
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
        // reproduce #660 behavior bit for bit, not merely approximately.
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
        // Behavioral counterpart to the arithmetic above: the curve must reach the
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
            "models without insertion_fraction must default to 0.4"
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
    /// A quality model over `options` at `read_length` bp, whose seed puts all its weight on
    /// `seed_weights` and whose every transition row keeps the previous score. A read from it
    /// therefore carries one value end to end, so "which mate produced this" is readable off
    /// the output.
    fn mate_model(
        options: Vec<usize>,
        read_length: usize,
        seed_weights: Vec<f64>,
    ) -> QualityScoreModel {
        let n = options.len();
        let stay: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        QualityScoreModel::from_counts(
            options,
            read_length,
            seed_weights,
            vec![stay; read_length - 1],
            false,
        )
        .unwrap()
    }

    /// #723: the two mates must be built against ONE option set, because
    /// `generate_quality_scores` indexes a single list -- a mate carrying its own list would
    /// index the wrong scores and emit plausible values from the wrong distribution.
    ///
    /// `gen-seq-error-model` cannot reach this: it builds both mates from the union of
    /// `global_counts` by construction. The guard exists for any other caller, and an
    /// unreachable guard with no test is a guard nobody knows still works.
    #[test]
    fn mates_built_against_different_option_sets_are_refused() {
        let r1 = mate_model(vec![20, 40], 10, vec![0.0, 1.0]);
        let r2 = mate_model(vec![20, 30, 40], 10, vec![1.0, 0.0, 0.0]);
        let model = SequencingErrorModel::from_raw_data(0.001, r1, None).unwrap();

        let err = model
            .with_mate_r2(r2)
            .expect_err("mates carrying different option sets must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("option sets") && msg.contains('2') && msg.contains('3'),
            "the refusal must name the mismatch and both sizes: {msg}"
        );
    }

    /// The other half of the same guard: one read length describes one run, and a mate fitted
    /// at a different length would be indexed past its own tensor.
    #[test]
    fn mates_fitted_at_different_read_lengths_are_refused() {
        let r1 = mate_model(vec![20, 40], 10, vec![0.0, 1.0]);
        let r2 = mate_model(vec![20, 40], 12, vec![1.0, 0.0]);
        let model = SequencingErrorModel::from_raw_data(0.001, r1, None).unwrap();

        let err = model
            .with_mate_r2(r2)
            .expect_err("mates fitted at different read lengths must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("read lengths") && msg.contains("10") && msg.contains("12"),
            "the refusal must name both lengths: {msg}"
        );
    }

    /// MUST NOT FIRE: an agreeing pair attaches, and the attached model is R2's rather than a
    /// second copy of R1's. Asserted on what the two models EMIT -- storing R1 twice passes
    /// any check that only asks whether an R2 population is present.
    #[test]
    fn an_agreeing_pair_of_mates_attaches_and_keeps_them_apart() {
        let r1 = mate_model(vec![20, 40], 10, vec![0.0, 1.0]);
        let r2 = mate_model(vec![20, 40], 10, vec![1.0, 0.0]);
        let model = SequencingErrorModel::from_raw_data(0.001, r1, None)
            .unwrap()
            .with_mate_r2(r2)
            .expect("mates agreeing on option set and read length must be accepted");

        let mut rng = make_rng();
        let s1 = model
            .quality_score_model()
            .generate_quality_scores(10, &mut rng)
            .unwrap();
        let s2 = model
            .quality_score_model_r2()
            .expect("the pair was accepted, so an R2 population must be present")
            .generate_quality_scores(10, &mut rng)
            .unwrap();
        assert!(
            s1.iter().all(|&q| q == 40),
            "R1's seed put all its weight on Q40: {s1:?}"
        );
        assert!(
            s2.iter().all(|&q| q == 20),
            "R2's seed put all its weight on Q20; Q40 here means R1 was stored twice: {s2:?}"
        );
    }
}
