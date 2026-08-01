//! `eidolon validate` — check an emitted artifact against the formats' consumers.
//!
//! Why this exists: our own output has shipped a malformed `AF=AF=0.3000` INFO value, a
//! VCF whose header did not declare tags its records used, and unsorted records. None of
//! those were caught downstream, because the tools are permissive — bcftools silently
//! converts a type-mismatched value to `.` rather than complaining, so the data simply
//! disappeared.
//!
//! The correctness criterion, stated before it was built: for every file in
//! `test_data/validation_corpus`, this must reach the same accept/reject verdict as the
//! tool that consumes it. That is a differential test, and it is what makes this
//! falsifiable rather than merely exercised. Where we are deliberately stricter, the
//! finding is a WARNING that says no tool rejects it.

pub mod errors;
pub mod utils;

use errors::ValidateError;
use log::{error, info, warn};
use std::path::Path;
use utils::finding::{Finding, Severity};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fastq,
    Vcf,
    Bam,
}

/// Infer the format from the extension, ignoring a trailing `.gz`.
pub fn format_from_path(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    if name.ends_with(".fq") || name.ends_with(".fastq") {
        Some(Format::Fastq)
    } else if name.ends_with(".vcf") {
        Some(Format::Vcf)
    } else if name.ends_with(".bam") {
        Some(Format::Bam)
    } else {
        None
    }
}

fn read_maybe_gzip(path: &Path) -> Result<String, ValidateError> {
    let bytes = std::fs::read(path).map_err(|source| ValidateError::Io {
        path: path.display().to_string(),
        source,
    })?;
    // gzip magic; BGZF is gzip so this covers bgzipped artifacts too.
    if bytes.starts_with(&[0x1f, 0x8b]) {
        use flate2::read::MultiGzDecoder;
        use std::io::Read;
        let mut s = String::new();
        MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .map_err(|source| ValidateError::Io {
                path: path.display().to_string(),
                source,
            })?;
        Ok(s)
    } else {
        Ok(String::from_utf8_lossy(&bytes).into_owned())
    }
}

/// Validate one file, returning every finding. Callers decide what to do with warnings.
pub fn validate_file(path: &Path, format: Option<Format>) -> Result<Vec<Finding>, ValidateError> {
    let format = format
        .or_else(|| format_from_path(path))
        .ok_or_else(|| ValidateError::UnknownFormat(path.display().to_string()))?;
    // BAM is binary and is read from the path directly rather than decompressed to a
    // String — a truncated BGZF stream must be detectable, and lossy UTF-8 conversion
    // would destroy exactly the evidence that check needs.
    if format == Format::Bam {
        return Ok(utils::bam::validate_bam(path));
    }
    let text = read_maybe_gzip(path)?;
    Ok(match format {
        Format::Fastq => utils::fastq::validate_fastq(&text),
        Format::Vcf => utils::vcf::validate_vcf(&text),
        Format::Bam => unreachable!("handled above"),
    })
}

/// Run validation and report. Returns `Err` when any ERROR-severity finding was made;
/// warnings are reported but do not fail, since by definition nothing downstream
/// rejects them.
pub fn run(paths: &[std::path::PathBuf], format: Option<Format>) -> Result<(), ValidateError> {
    // Only the files that actually failed. Listing every path checked — including ones
    // that passed — would misattribute the failure, which is the opposite of what this
    // subcommand is for.
    let mut failed: Vec<String> = Vec::new();
    for path in paths {
        let findings = validate_file(path, format)?;
        let errors = findings
            .iter()
            .filter(|f| f.severity == Severity::Error)
            .count();
        let warnings = findings.len() - errors;

        if findings.is_empty() {
            info!("{}: OK", path.display());
            continue;
        }
        for f in &findings {
            match f.severity {
                Severity::Error => error!("{}\n{f}", path.display()),
                Severity::Warning => warn!("{}\n{f}", path.display()),
            }
        }
        info!(
            "{}: {errors} error(s), {warnings} warning(s)",
            path.display()
        );
        if errors > 0 {
            failed.push(path.display().to_string());
        }
    }
    if !failed.is_empty() {
        return Err(ValidateError::Invalid(failed.join(", ")));
    }
    Ok(())
}
