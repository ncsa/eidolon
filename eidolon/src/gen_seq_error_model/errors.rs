use eidolon_core::{
    file_tools::{bam_reader::BamReaderError, vcf_tools::VcfToolsError},
    models::{quality_scores::QualityModelError, sequencing_error_model::SeqModelError},
    structs::transition_matrix::TransitionMatrixError,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum GenSeqErrorModelError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Configuration error: {0}")]
    ConfigurationError(String),
    #[error("CLI error: {0}")]
    CliError(String),
    #[error("Malformed FASTQ: {0}")]
    MalformedFastq(String),
    #[error("Quality score model error: {0}")]
    QualModelError(#[from] QualityModelError),
    #[error("Sequencing error model error: {0}")]
    SeqModelError(#[from] SeqModelError),
    #[error("BAM reader error: {0}")]
    BamError(#[from] BamReaderError),
    #[error("Transition matrix error: {0}")]
    TransitionMatrixError(#[from] TransitionMatrixError),
    #[error("Error reading known_variants_vcf: {0}")]
    KnownVariantsError(#[from] VcfToolsError),
}
