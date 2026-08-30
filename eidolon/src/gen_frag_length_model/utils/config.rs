use crate::gen_frag_length_model::errors::GenFragLengthModelError;
use crate::gen_frag_length_model::utils::runner::DistributionKind;
use eidolon_core::file_tools::file_io::check_overwrite;
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    /// Path to the input BAM or SAM file.
    pub input_file: PathBuf,
    /// Path to the output model file (should end in `.json.gz`).
    pub output_file: PathBuf,
    /// If false, refuse to overwrite an existing output file.
    pub overwrite_output: bool,
    /// Minimum number of reads for a fragment length to be included in the model.
    /// Default is 2 to handle smaller datasets. Set to 0 to disable filtering.
    /// For larger datasets, try 100 and adjust from there.
    pub min_reads: usize,
    /// `discrete` (default) keeps the measured shape; `normal` collapses it to mean +
    /// st_dev, which is the pre-v3.3.0 behaviour and is right for sparse inputs.
    pub distribution: DistributionKind,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenFragLengthModelError> {
        let f = fs::File::open(yml_file).map_err(|e| {
            GenFragLengthModelError::ConfigurationError(format!(
                "Cannot open config file {:?}: {}",
                yml_file, e
            ))
        })?;
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(f)
            .map_err(|e| GenFragLengthModelError::ConfigurationError(e.to_string()))?;

        let input_file = PathBuf::from(
            scrape_config
                .get("input_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    GenFragLengthModelError::ConfigurationError(
                        "Missing required field: input_file".to_string(),
                    )
                })?,
        );
        if !input_file.is_file() {
            return Err(GenFragLengthModelError::ConfigurationError(format!(
                "input_file not found: {:?}",
                input_file
            )));
        }

        let output_file = PathBuf::from(
            scrape_config
                .get("output_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    GenFragLengthModelError::ConfigurationError(
                        "Missing required field: output_file".to_string(),
                    )
                })?,
        );

        let overwrite_output = scrape_config
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        check_overwrite("output_file", &output_file, overwrite_output)
            .map_err(|e| GenFragLengthModelError::ConfigurationError(e.to_string()))?;

        let min_reads = scrape_config
            .get("min_reads")
            .and_then(|v| v.as_u64())
            .unwrap_or(2) as usize;

        let distribution = match scrape_config.get("distribution").and_then(|v| v.as_str()) {
            None => DistributionKind::default(),
            Some(raw) => DistributionKind::parse(raw).ok_or_else(|| {
                GenFragLengthModelError::ConfigurationError(format!(
                    "distribution must be `discrete` or `normal`, got `{raw}`"
                ))
            })?,
        };

        Ok(RunConfiguration {
            input_file,
            output_file,
            overwrite_output,
            min_reads,
            distribution,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_bam_placeholder(dir: &std::path::Path) -> PathBuf {
        let p = dir.join("input.bam");
        std::fs::write(&p, b"placeholder").unwrap();
        p
    }

    #[test]
    fn test_run_configuration_defaults() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let bam = format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir);
        let out_dir = tempfile::tempdir().unwrap();
        let output = out_dir.path().join("model.json.gz");
        let yaml = format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: true\n",
            bam,
            output.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.min_reads, 2);
        assert!(config.overwrite_output);
    }

    #[test]
    fn test_custom_min_reads() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let bam = format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir);
        let out_dir = tempfile::tempdir().unwrap();
        let output = out_dir.path().join("model.json.gz");
        let yaml = format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 100\n",
            bam,
            output.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.min_reads, 100);
    }

    #[test]
    fn test_missing_input_file_errors() {
        let out_dir = tempfile::tempdir().unwrap();
        let output = out_dir.path().join("model.json.gz");
        let yaml = format!(
            "input_file: /nonexistent/input.bam\noutput_file: {}\n",
            output.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let result = RunConfiguration::from(&tmp.path().to_path_buf());
        assert!(result.is_err());
    }

    #[test]
    fn test_overwrite_false_blocks_existing_output() {
        let dir = tempfile::tempdir().unwrap();
        let bam = make_bam_placeholder(dir.path());
        let output = dir.path().join("model.json.gz");
        std::fs::write(&output, b"existing").unwrap();
        let yaml = format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: false\n",
            bam.display(),
            output.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let result = RunConfiguration::from(&tmp.path().to_path_buf());
        assert!(result.is_err());
    }

    #[test]
    fn test_min_reads_zero_disables_filter() {
        let dir = tempfile::tempdir().unwrap();
        let bam = make_bam_placeholder(dir.path());
        let out = dir.path().join("model.json.gz");
        let yaml = format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 0\n",
            bam.display(),
            out.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.min_reads, 0);
    }
}
