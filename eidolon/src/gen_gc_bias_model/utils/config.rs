// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use crate::config_keys::check_keys;
use crate::gen_gc_bias_model::errors::GenGcBiasModelError;
use eidolon_core::{
    file_tools::{bed_reader::read_bed, file_io::check_overwrite},
    structs::bed_record::BedRecord,
};
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

/// Every top-level key the config parser reads. Anything else is rejected (#496).
pub const KNOWN_KEYS: &[&str] = &[
    "bam_file",
    "bed_file",
    "min_mapq",
    "min_windows_per_bin",
    "output_file",
    "overwrite_output",
    "reference",
    "window_size",
    "window_stride",
];

/// Parameters consumed by `run_from_coverage` — the FASTA-sweep + model-write
/// stage. Does not carry BAM/walker state because `run_from_coverage` works
/// from pre-computed per-contig coverage and never reads the BAM itself.
///
/// Both the standalone CLI path (which fills these from a YAML
/// `RunConfiguration`) and the unified `gen-bam-models` runner (which fills
/// these from its own `GcBiasSection` after one shared BAM walk) construct
/// this directly.
#[derive(Debug, Clone)]
pub struct GcBiasModelParams {
    pub reference: PathBuf,
    pub bed_table: HashMap<String, Vec<BedRecord>>,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
    pub window_size: usize,
    pub window_stride: usize,
    pub min_windows_per_bin: usize,
}

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for running GC Bias modeling
    // from a standalone YAML config: BAM walker inputs + model params.
    /// Aligned BAM. Coverage is accumulated in process via
    /// `eidolon_core::file_tools::bam_reader::CoverageObserver`.
    pub bam_file: PathBuf,
    /// Minimum mapping quality applied while accumulating coverage from `bam_file`.
    /// Default 0 matches `samtools depth` defaults.
    pub min_mapq: u8,
    pub model: GcBiasModelParams,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenGcBiasModelError> {
        let file = fs::File::open(yml_file)?;
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)?;
        check_keys(
            scrape_config.keys().map(String::as_str),
            KNOWN_KEYS,
            "gen-gc-bias-model",
        )
        .map_err(GenGcBiasModelError::ConfigError)?;

        let reference =
            PathBuf::from(scrape_config["reference"].as_str().ok_or_else(|| {
                GenGcBiasModelError::ConfigError("Missing reference".to_string())
            })?);
        if !reference.is_file() {
            return Err(GenGcBiasModelError::ConfigError(format!(
                "Invalid reference file {:?}",
                reference
            )));
        }

        let bam_file = PathBuf::from(
            scrape_config
                .get("bam_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| GenGcBiasModelError::ConfigError("Missing bam_file".to_string()))?,
        );
        if !bam_file.is_file() {
            return Err(GenGcBiasModelError::ConfigError(format!(
                "Invalid bam_file {:?}",
                bam_file
            )));
        }

        let min_mapq = scrape_config
            .get("min_mapq")
            .and_then(|v| v.as_u64())
            .unwrap_or(0) as u8;

        let bed_file_raw = scrape_config
            .get("bed_file")
            .and_then(|v| v.as_str())
            .unwrap_or(".");

        let bed_table = if bed_file_raw == "." {
            HashMap::new()
        } else {
            let bed_file = PathBuf::from(bed_file_raw);
            if !bed_file.is_file() {
                return Err(GenGcBiasModelError::ConfigError(format!(
                    "Invalid BED file {:?}",
                    bed_file
                )));
            }
            read_bed(&bed_file, false)?
        };

        let overwrite_output = scrape_config
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let output_file =
            PathBuf::from(scrape_config["output_file"].as_str().ok_or_else(|| {
                GenGcBiasModelError::ConfigError("Missing output_file".to_string())
            })?);

        check_overwrite("output_file", &output_file, overwrite_output)
            .map_err(|e| GenGcBiasModelError::ConfigError(e.to_string()))?;

        let window_size = scrape_config
            .get("window_size")
            .and_then(|v| v.as_u64())
            .unwrap_or(100) as usize;

        if window_size == 0 {
            return Err(GenGcBiasModelError::ConfigError(
                "window_size must be > 0".to_string(),
            ));
        }

        let window_stride = scrape_config
            .get("window_stride")
            .and_then(|v| v.as_u64())
            .unwrap_or(window_size as u64) as usize;

        if window_stride == 0 {
            return Err(GenGcBiasModelError::ConfigError(
                "window_stride must be > 0".to_string(),
            ));
        }

        if window_stride > window_size {
            return Err(GenGcBiasModelError::ConfigError(format!(
                "window_stride ({}) must not exceed window_size ({}); \
                 for non-overlapping windows use window_stride = window_size",
                window_stride, window_size
            )));
        }

        let min_windows_per_bin = scrape_config
            .get("min_windows_per_bin")
            .and_then(|v| v.as_u64())
            .unwrap_or(10) as usize;

        Ok(RunConfiguration {
            bam_file,
            min_mapq,
            model: GcBiasModelParams {
                reference,
                bed_table,
                output_file,
                overwrite_output,
                window_size,
                window_stride,
                min_windows_per_bin,
            },
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

    fn write_temp_file(contents: &str) -> tempfile::NamedTempFile {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", contents).unwrap();
        tmp
    }

    fn base_yaml(reference: &PathBuf, bam_file: &PathBuf, output_file: &PathBuf) -> String {
        format!(
            "\
reference: {}
bam_file: {}
bed_file: .
output_file: {}
overwrite_output: true
window_size: 100
window_stride: 100
min_windows_per_bin: 10
",
            reference.display(),
            bam_file.display(),
            output_file.display(),
        )
    }

    #[test]
    fn test_run_configuration_parses_required_fields() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = base_yaml(
            &reference.path().to_path_buf(),
            &bam.path().to_path_buf(),
            &output_path,
        );
        let config_file = write_config(&yaml);

        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert_eq!(config.model.reference, reference.path());
        assert_eq!(config.bam_file, bam.path());
        assert_eq!(config.min_mapq, 0);
        assert!(config.model.bed_table.is_empty());
        assert_eq!(config.model.output_file, output_path);
        assert!(config.model.overwrite_output);
        assert_eq!(config.model.window_size, 100);
        assert_eq!(config.model.window_stride, 100);
        assert_eq!(config.model.min_windows_per_bin, 10);
    }

    #[test]
    fn test_run_configuration_defaults_optional_fields() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n",
            reference.path().display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        assert!(config.model.bed_table.is_empty());
        assert_eq!(config.model.window_size, 100);
        assert_eq!(config.model.window_stride, 100);
        assert_eq!(config.model.min_windows_per_bin, 10);
        assert_eq!(config.min_mapq, 0);
    }

    #[test]
    fn test_run_configuration_parses_min_mapq() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\nmin_mapq: 20\n",
            reference.path().display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();
        assert_eq!(config.min_mapq, 20);
    }

    #[test]
    fn test_run_configuration_parses_bed_file() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let bed = write_temp_file("chr1\t0\t4\nchr1\t4\t8\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\nbed_file: {}\noutput_file: {}\noverwrite_output: true\n",
            reference.path().display(),
            bam.path().display(),
            bed.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();

        let chr1_records = config
            .model
            .bed_table
            .get("chr1")
            .expect("Expected chr1 BED records");
        assert_eq!(chr1_records.len(), 2);
    }

    #[test]
    fn test_run_configuration_rejects_missing_reference_file() {
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let missing_reference = PathBuf::from("/tmp/eidolon_missing_reference_for_gc_bias.fa");
        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n",
            missing_reference.display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_rejects_missing_bam_file() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let missing_bam = PathBuf::from("/tmp/eidolon_missing_bam_for_gc_bias.bam");
        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n",
            reference.path().display(),
            missing_bam.display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_rejects_bam_file_omitted() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\noutput_file: {}\noverwrite_output: true\n",
            reference.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_rejects_missing_bed_file() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let missing_bed = PathBuf::from("/tmp/eidolon_missing_bed_for_gc_bias.bed");
        let yaml = format!(
            "reference: {}\nbam_file: {}\nbed_file: {}\noutput_file: {}\noverwrite_output: true\n",
            reference.path().display(),
            bam.path().display(),
            missing_bed.display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_rejects_zero_window_size() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\nwindow_size: 0\nwindow_stride: 100\n",
            reference.path().display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_rejects_zero_window_stride() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\nwindow_size: 100\nwindow_stride: 0\n",
            reference.path().display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_rejects_window_stride_greater_than_window_size() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\nwindow_size: 100\nwindow_stride: 200\n",
            reference.path().display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_allows_window_stride_equal_to_window_size() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = tempfile::NamedTempFile::new().unwrap();
        let output_path = output.path().to_path_buf();
        drop(output);

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\nwindow_size: 100\nwindow_stride: 100\n",
            reference.path().display(),
            bam.path().display(),
            output_path.display(),
        );
        let config_file = write_config(&yaml);
        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();
        assert_eq!(config.model.window_stride, 100);
    }

    #[test]
    fn test_run_configuration_rejects_existing_output_when_overwrite_false() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = write_temp_file("already exists\n");

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: false\n",
            reference.path().display(),
            bam.path().display(),
            output.path().display(),
        );
        let config_file = write_config(&yaml);
        assert!(RunConfiguration::from(&config_file.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_run_configuration_allows_existing_output_when_overwrite_true() {
        let reference = write_temp_file(">chr1\nACGTACGT\n");
        let bam = write_temp_file("dummy");
        let output = write_temp_file("already exists\n");

        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n",
            reference.path().display(),
            bam.path().display(),
            output.path().display(),
        );
        let config_file = write_config(&yaml);
        let config = RunConfiguration::from(&config_file.path().to_path_buf()).unwrap();
        assert_eq!(config.model.output_file, output.path());
        assert!(config.model.overwrite_output);
    }
}
