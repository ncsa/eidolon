use crate::gen_seq_error_model::errors::GenSeqErrorModelError;
use eidolon_core::file_tools::file_io::check_overwrite;
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    pub fastq_file: PathBuf,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
    pub max_reads: usize,
    pub qual_offset: usize,
    /// Optional list of Q-score bins (e.g. `[2, 12, 23, 37]` for NovaSeq). When set,
    /// learned scores are snapped to the nearest bin and the resulting model only
    /// emits these values. Sorted/deduped at parse time. Value 31 is rejected
    /// because it encodes to `@` under Phred+33.
    pub binned_quality_bins: Option<Vec<usize>>,
    /// Optional aligned BAM file to infer the SNP transition matrix from read-vs-reference
    /// mismatches. Requires MD tags. Ignored if transition_matrix_file is also set.
    pub bam_file: Option<PathBuf>,
    /// Optional path to a 4×4 TSV specifying a custom SNP transition matrix.
    /// Rows/columns are A/C/G/T. A single header line is ignored.
    /// Takes precedence over bam_file.
    pub transition_matrix_file: Option<PathBuf>,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenSeqErrorModelError> {
        let f = fs::File::open(yml_file).map_err(|e| {
            GenSeqErrorModelError::ConfigurationError(format!(
                "Cannot open config file {:?}: {}",
                yml_file, e
            ))
        })?;
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(f)
            .map_err(|e| GenSeqErrorModelError::ConfigurationError(e.to_string()))?;

        let fastq_file = PathBuf::from(
            scrape_config
                .get("fastq_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    GenSeqErrorModelError::ConfigurationError(
                        "Missing required field: fastq_file".to_string(),
                    )
                })?,
        );
        if !fastq_file.is_file() {
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "fastq_file not found: {:?}",
                fastq_file
            )));
        }

        let output_file = PathBuf::from(
            scrape_config
                .get("output_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| {
                    GenSeqErrorModelError::ConfigurationError(
                        "Missing required field: output_file".to_string(),
                    )
                })?,
        );

        let overwrite_output = scrape_config
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        check_overwrite("output_file", &output_file, overwrite_output)
            .map_err(|e| GenSeqErrorModelError::ConfigurationError(e.to_string()))?;

        let max_reads = scrape_config
            .get("max_reads")
            .and_then(|v| v.as_u64())
            .unwrap_or(0) as usize;

        let qual_offset = scrape_config
            .get("qual_offset")
            .and_then(|v| v.as_u64())
            .unwrap_or(33) as usize;

        let binned_quality_bins = match scrape_config.get("binned_quality_bins") {
            None => None,
            Some(v) => {
                let seq = v.as_sequence().ok_or_else(|| {
                    GenSeqErrorModelError::ConfigurationError(
                        "binned_quality_bins must be a YAML list of integers".to_string(),
                    )
                })?;
                let mut bins: Vec<usize> = Vec::with_capacity(seq.len());
                for elem in seq {
                    let n = elem.as_u64().ok_or_else(|| {
                        GenSeqErrorModelError::ConfigurationError(
                            "binned_quality_bins entries must be non-negative integers".to_string(),
                        )
                    })?;
                    bins.push(n as usize);
                }
                if bins.is_empty() {
                    return Err(GenSeqErrorModelError::ConfigurationError(
                        "binned_quality_bins must contain at least one value".to_string(),
                    ));
                }
                if bins.iter().any(|&b| b >= 94) {
                    return Err(GenSeqErrorModelError::ConfigurationError(
                        "binned_quality_bins entries must be < 94".to_string(),
                    ));
                }
                if bins.contains(&31) {
                    return Err(GenSeqErrorModelError::ConfigurationError(
                        "binned_quality_bins cannot include 31 (encodes to '@' under \
                         Phred+33 and would corrupt FASTQ output)"
                            .to_string(),
                    ));
                }
                bins.sort_unstable();
                bins.dedup();
                Some(bins)
            }
        };

        let bam_file = scrape_config
            .get("bam_file")
            .and_then(|v| v.as_str())
            .map(|s| {
                let p = PathBuf::from(s);
                if !p.is_file() {
                    Err(GenSeqErrorModelError::ConfigurationError(format!(
                        "bam_file not found: {:?}",
                        p
                    )))
                } else {
                    Ok(p)
                }
            })
            .transpose()?;

        let transition_matrix_file = scrape_config
            .get("transition_matrix_file")
            .and_then(|v| v.as_str())
            .map(|s| {
                let p = PathBuf::from(s);
                if !p.is_file() {
                    Err(GenSeqErrorModelError::ConfigurationError(format!(
                        "transition_matrix_file not found: {:?}",
                        p
                    )))
                } else {
                    Ok(p)
                }
            })
            .transpose()?;

        Ok(RunConfiguration {
            fastq_file,
            output_file,
            overwrite_output,
            max_reads,
            qual_offset,
            binned_quality_bins,
            bam_file,
            transition_matrix_file,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_config(yaml: &str) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        tmp
    }

    fn make_fastq(dir: &tempfile::TempDir) -> PathBuf {
        let p = dir.path().join("test.fastq");
        std::fs::write(&p, "@read1\nACGT\n+\nIIII\n").unwrap();
        p
    }

    #[test]
    fn test_run_configuration() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 1000\nqual_offset: 33\n",
            fastq.display(),
            output.display()
        ));
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.max_reads, 1000);
        assert_eq!(config.qual_offset, 33);
        assert!(config.overwrite_output);
        assert!(config.bam_file.is_none());
        assert!(config.transition_matrix_file.is_none());
    }

    #[test]
    fn test_missing_fastq_file_errors() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().join("model.json.gz");
        let tmp = write_config(&format!(
            "fastq_file: /nonexistent/path.fastq\noutput_file: {}\n",
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_output_file_overwrite_false_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");
        std::fs::write(&output, "existing").unwrap();

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: false\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_defaults_applied() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\n",
            fastq.display(),
            output.display()
        ));
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.max_reads, 0);
        assert_eq!(config.qual_offset, 33);
        assert!(!config.overwrite_output);
    }

    #[test]
    fn test_bam_file_parsed() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let bam_path = dir.path().join("reads.bam");
        std::fs::write(&bam_path, b"dummy").unwrap();
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbam_file: {}\n",
            fastq.display(),
            output.display(),
            bam_path.display()
        ));
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert!(config.bam_file.is_some());
        assert!(config.transition_matrix_file.is_none());
    }

    #[test]
    fn test_transition_matrix_file_parsed() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let tsv_path = dir.path().join("matrix.tsv");
        std::fs::write(
            &tsv_path,
            "A\tC\tG\tT\n0\t0.5\t0.3\t0.2\n0.5\t0\t0.3\t0.2\n0.4\t0.3\t0\t0.3\n0.3\t0.3\t0.4\t0\n",
        )
        .unwrap();
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\ntransition_matrix_file: {}\n",
            fastq.display(),
            output.display(),
            tsv_path.display()
        ));
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert!(config.transition_matrix_file.is_some());
        assert!(config.bam_file.is_none());
    }

    #[test]
    fn test_missing_bam_file_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbam_file: /nonexistent/reads.bam\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_binned_quality_bins_parsed_sorted_deduped() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbinned_quality_bins: [37, 2, 23, 12, 12]\n",
            fastq.display(),
            output.display()
        ));
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(config.binned_quality_bins, Some(vec![2, 12, 23, 37]));
    }

    #[test]
    fn test_binned_quality_bins_absent_is_none() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\n",
            fastq.display(),
            output.display()
        ));
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        assert!(config.binned_quality_bins.is_none());
    }

    #[test]
    fn test_binned_quality_bins_empty_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbinned_quality_bins: []\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_binned_quality_bins_out_of_range_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbinned_quality_bins: [2, 12, 94]\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_binned_quality_bins_thirty_one_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbinned_quality_bins: [2, 12, 31, 37]\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_binned_quality_bins_non_list_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nbinned_quality_bins: \"2,12,23,37\"\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_missing_transition_matrix_file_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fastq = make_fastq(&dir);
        let output = dir.path().join("model.json.gz");

        let tmp = write_config(&format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\ntransition_matrix_file: /nonexistent/matrix.tsv\n",
            fastq.display(),
            output.display()
        ));
        assert!(RunConfiguration::from(&tmp.path().to_path_buf()).is_err());
    }
}
