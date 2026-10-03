use std::{
    num::{ParseFloatError, ParseIntError},
    str::ParseBoolError,
};

use eidolon_core::rng::NeatRngError;
use eidolon_core::{
    file_tools::{
        bam_writer::BamWriterError, bed_reader::BedReaderError, fasta_stream::FastaStreamError,
        fastq_tools::FastqToolsError, vcf_tools::VcfToolsError,
    },
    models::{
        fragment_length::FragmentModelError, gc_bias_model::GcBiasModelError,
        mutation_model::MutationModelError, quality_scores::QualityModelError,
        sequencing_error_model::SeqModelError,
    },
    structs::{
        distributions::DistributionErrors, mutated_map::MutatedMapError,
        sequence_block::SequenceBlockError,
    },
};
use serde_yml;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GenerateReadsError {
    #[error("Invalid CLI inputs: {0}")]
    CliError(String),
    #[error("NEAT generate-reads requires a reference to proceed.")]
    MissingReferenceError,
    #[error("Error processing Configuration for the run")]
    MainConfigurationError,
    #[error("Error while applying variants!")]
    ApplyVariantsError,
    #[error("Error reading config key: {0} accepts values of {1}")]
    ConfigReadError(String, String),
    #[error("Error generating configuration!")]
    ConfigError,
    #[error("{0}")]
    UnknownConfigKeys(String),
    #[error("Error generating fragments!")]
    GenerateFragmentsError,
    #[error("Error generating variants!")]
    GenerateVariantsError,
    #[error("Error while mutating fasta!")]
    MutateFastaError,
    #[error("Generate Reads Runner reported an error!")]
    RunnerError,
    #[error("File not found: {0}")]
    FileNotFound(String),
    #[error("IO Error called by NEAT: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Error generating fragment length model: {0}")]
    FragModelError(#[from] FragmentModelError),
    #[error("Error creating quality score model: {0}")]
    QualModelError(#[from] QualityModelError),
    #[error("Error creating mutation model: {0}")]
    MutModelError(#[from] MutationModelError),
    #[error("Error creating sequencing error model: {0}")]
    SeqModelError(#[from] SeqModelError),
    #[error("SequenceBlock error: {0}")]
    SequenceBlockError(#[from] SequenceBlockError),
    #[error("Error creating distributions: {0}")]
    GenReadsDistroError(#[from] DistributionErrors),
    #[error("Error sampling distro: {0}")]
    GenReadsRngError(#[from] NeatRngError),
    #[error("Error creating mutated map for block: {0}")]
    MutMapError(#[from] MutatedMapError),
    #[error("Error during configuration creation: {0}")]
    ConfigYamlError(#[from] serde_yml::Error),
    #[error("Configuration reader reported an error parsing an int: {0}")]
    ConfigParseIntError(#[from] ParseIntError),
    #[error("Configuration reader reported an error parsing a float: {0}")]
    ConfigParseFloatError(#[from] ParseFloatError),
    #[error("Configuration reader reported an error parsing a bool: {0}")]
    ConfigParseBoolError(#[from] ParseBoolError),
    #[error("Main threw an error: {0}")]
    MainError(#[from] Box<dyn std::error::Error + Send + Sync>),
    #[error("FastqTools error: {0}")]
    FqToolsError(#[from] FastqToolsError),
    #[error("Bed reader error: {0}")]
    BedError(#[from] BedReaderError),
    #[error("Input VCF error: {0}")]
    InputVcfError(#[from] VcfToolsError),
    #[error("BAM writer error: {0}")]
    BamWriterError(#[from] BamWriterError),
    #[error("Sequence too short to process")]
    ShortSequence,
    #[error("GC Bias model threw an error: {0}")]
    BiasModelError(#[from] GcBiasModelError),
    #[error("Error streaming FASTA: {0}")]
    FastaStreamError(#[from] FastaStreamError),
}
