use thiserror::Error;

#[derive(Error, Debug)]
pub enum CompareAfError {
    #[error("could not read {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    /// Carries the message the Python implementation passed to `sys.exit`, so callers
    /// (and the harness's line-parsing gate) see identical text.
    #[error("{0}")]
    Fatal(String),
}
