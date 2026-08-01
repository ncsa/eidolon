use thiserror::Error;

#[derive(Error, Debug)]
pub enum ValidateError {
    #[error("could not read {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error(
        "cannot determine the format of {0} from its extension — pass --format fastq|vcf \
         (recognised: .fq/.fastq/.vcf, optionally .gz)"
    )]
    UnknownFormat(String),
    #[error("{0} failed validation")]
    Invalid(String),
}
