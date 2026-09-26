use eidolon_core::{
    file_tools::{
        bed_reader::BedReaderError, fasta_stream::FastaStreamError, vcf_tools::VcfToolsError,
    },
    models::mutation_model::MutationModelError,
    structs::transition_matrix::TransitionMatrixError,
};
use std::num::ParseIntError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum GenMutationModelError {
    #[error("Error while filtering the vcf file {0}")]
    VcfFilterError(String),
    #[error("Error while filtering the fastq file {0}")]
    FastqFilterError(String),
    #[error("File does not exist: {0}")]
    FileNotFound(String),
    #[error("Unable to parse file name for extension {0}")]
    MalformedFileName(String),
    #[error("I/O error during filtering of files")]
    IoError(#[from] std::io::Error),
    #[error("Error parsing coordinates from read name: {0}")]
    CoordParseError(#[from] ParseIntError),
    #[error("The bed reader returned an error: {0}")]
    BedReaderErr(#[from] BedReaderError),
    #[error("Unknown char while parsing suspected read name")]
    UnknownChar,
    #[error("Cannot overwrite existing file: {0}")]
    OverwriteFileError(String),
    #[error("Invalid CLI inputs: {0}")]
    CliError(String),
    #[error("Invalid Configuration inputs: {0}")]
    ConfigurationError(String),
    #[error("Error reading the vcf file: {0}")]
    VcfReaderError(#[from] VcfToolsError),
    #[error("Error generating trinucleotide counts: {0}")]
    TrinucCountError(String),
    #[error("Training VCF has no SNP / insertion / deletion observations: {0}")]
    NoLiteralVariants(String),
    #[error("Error streaming FASTA file: {0}")]
    FastaStreamError(#[from] FastaStreamError),
    #[error("Variant location did not match reference: {0}")]
    BaseMismatch(String),
    #[error("Error building mutation model: {0}")]
    MutModelError(#[from] MutationModelError),
    #[error("Transition matrix error: {0}")]
    TransitionMatrixError(#[from] TransitionMatrixError),
    #[error(
        "The VCF has {counted} SNP(s) but none could be used: {ref_mismatch} did not match \
         the reference base and {edge} were at a contig edge. A VCF called against a \
         different reference build produces this. No model was written."
    )]
    NoUsableSnps {
        counted: usize,
        ref_mismatch: usize,
        edge: usize,
    },
}
