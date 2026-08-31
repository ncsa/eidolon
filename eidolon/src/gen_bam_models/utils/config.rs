// Unified config for `gen-bam-models`: one BAM, optional per-model sections.
// At least one of `frag_length` or `gc_bias` must be present. The shared
// BAM walk uses `BamWalkFilter::for_coverage()` with the top-level
// `min_mapq` (default 0); each observer self-filters for its own criteria.
use crate::gen_frag_length_model::utils::runner::DistributionKind;
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use eidolon_core::{
    file_tools::{bed_reader::read_bed, file_io::check_overwrite},
    structs::bed_record::BedRecord,
};

use crate::gen_bam_models::errors::GenBamModelsError;

#[derive(Debug, Clone)]
pub struct FragLengthSection {
    pub output_file: PathBuf,
    pub overwrite_output: bool,
    pub min_reads: usize,
    pub distribution: DistributionKind,
}

#[derive(Debug, Clone)]
pub struct GcBiasSection {
    pub reference: PathBuf,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
    pub window_size: usize,
    pub window_stride: usize,
    pub min_windows_per_bin: usize,
    pub bed_table: HashMap<String, Vec<BedRecord>>,
}

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    pub bam_file: PathBuf,
    /// Walker-level MAPQ filter applied during the shared BAM pass.
    /// Default 0 (matches `samtools depth` defaults). Individual observers
    /// may apply stricter filters internally (e.g. FragLengthObserver
    /// requires MAPQ > 10 regardless of this value).
    pub min_mapq: u8,
    pub frag_length: Option<FragLengthSection>,
    pub gc_bias: Option<GcBiasSection>,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenBamModelsError> {
        let file = fs::File::open(yml_file)?;
        let scrape: HashMap<String, Value> = serde_yml::from_reader(file)?;

        let bam_file = PathBuf::from(
            scrape
                .get("bam_file")
                .and_then(|v| v.as_str())
                .ok_or_else(|| GenBamModelsError::ConfigError("Missing bam_file".to_string()))?,
        );
        if !bam_file.is_file() {
            return Err(GenBamModelsError::ConfigError(format!(
                "Invalid bam_file {:?}",
                bam_file
            )));
        }

        let min_mapq = scrape.get("min_mapq").and_then(|v| v.as_u64()).unwrap_or(0) as u8;

        let frag_length = scrape
            .get("frag_length")
            .map(parse_frag_length_section)
            .transpose()?;

        let gc_bias = scrape
            .get("gc_bias")
            .map(parse_gc_bias_section)
            .transpose()?;

        if frag_length.is_none() && gc_bias.is_none() {
            return Err(GenBamModelsError::ConfigError(
                "Must request at least one model (frag_length or gc_bias)".to_string(),
            ));
        }

        Ok(RunConfiguration {
            bam_file,
            min_mapq,
            frag_length,
            gc_bias,
        })
    }
}

fn parse_frag_length_section(v: &Value) -> Result<FragLengthSection, GenBamModelsError> {
    let map = v.as_mapping().ok_or_else(|| {
        GenBamModelsError::ConfigError("frag_length section must be a mapping".to_string())
    })?;

    let get_str = |k: &str| {
        map.get(Value::String(k.to_string()))
            .and_then(|v| v.as_str())
    };
    let get_u64 = |k: &str| {
        map.get(Value::String(k.to_string()))
            .and_then(|v| v.as_u64())
    };
    let get_bool = |k: &str| {
        map.get(Value::String(k.to_string()))
            .and_then(|v| v.as_bool())
    };

    let output_file = PathBuf::from(get_str("output_file").ok_or_else(|| {
        GenBamModelsError::ConfigError("frag_length.output_file is required".to_string())
    })?);
    let overwrite_output = get_bool("overwrite_output").unwrap_or(false);
    check_overwrite("frag_length.output_file", &output_file, overwrite_output)
        .map_err(|e| GenBamModelsError::ConfigError(e.to_string()))?;
    let min_reads = get_u64("min_reads").unwrap_or(2) as usize;
    let distribution = match get_str("distribution") {
        None => DistributionKind::default(),
        Some(raw) => DistributionKind::parse(raw).ok_or_else(|| {
            GenBamModelsError::ConfigError(format!(
                "frag_length.distribution must be `discrete` or `normal`, got `{raw}`"
            ))
        })?,
    };

    Ok(FragLengthSection {
        output_file,
        overwrite_output,
        min_reads,
        distribution,
    })
}

fn parse_gc_bias_section(v: &Value) -> Result<GcBiasSection, GenBamModelsError> {
    let map = v.as_mapping().ok_or_else(|| {
        GenBamModelsError::ConfigError("gc_bias section must be a mapping".to_string())
    })?;

    let get_str = |k: &str| {
        map.get(Value::String(k.to_string()))
            .and_then(|v| v.as_str())
    };
    let get_u64 = |k: &str| {
        map.get(Value::String(k.to_string()))
            .and_then(|v| v.as_u64())
    };
    let get_bool = |k: &str| {
        map.get(Value::String(k.to_string()))
            .and_then(|v| v.as_bool())
    };

    let reference = PathBuf::from(get_str("reference").ok_or_else(|| {
        GenBamModelsError::ConfigError("gc_bias.reference is required".to_string())
    })?);
    if !reference.is_file() {
        return Err(GenBamModelsError::ConfigError(format!(
            "Invalid gc_bias.reference {:?}",
            reference
        )));
    }

    let output_file = PathBuf::from(get_str("output_file").ok_or_else(|| {
        GenBamModelsError::ConfigError("gc_bias.output_file is required".to_string())
    })?);
    let overwrite_output = get_bool("overwrite_output").unwrap_or(false);
    check_overwrite("gc_bias.output_file", &output_file, overwrite_output)
        .map_err(|e| GenBamModelsError::ConfigError(e.to_string()))?;

    let window_size = get_u64("window_size").unwrap_or(100) as usize;
    if window_size == 0 {
        return Err(GenBamModelsError::ConfigError(
            "gc_bias.window_size must be > 0".to_string(),
        ));
    }
    let window_stride = get_u64("window_stride").unwrap_or(window_size as u64) as usize;
    if window_stride == 0 {
        return Err(GenBamModelsError::ConfigError(
            "gc_bias.window_stride must be > 0".to_string(),
        ));
    }
    if window_stride > window_size {
        return Err(GenBamModelsError::ConfigError(format!(
            "gc_bias.window_stride ({}) must not exceed window_size ({})",
            window_stride, window_size
        )));
    }
    let min_windows_per_bin = get_u64("min_windows_per_bin").unwrap_or(10) as usize;

    let bed_file_raw = get_str("bed_file").unwrap_or(".");
    let bed_table = if bed_file_raw == "." {
        HashMap::new()
    } else {
        let bed_file = PathBuf::from(bed_file_raw);
        if !bed_file.is_file() {
            return Err(GenBamModelsError::ConfigError(format!(
                "Invalid gc_bias.bed_file {:?}",
                bed_file
            )));
        }
        read_bed(&bed_file, false)
            .map_err(|e| GenBamModelsError::ConfigError(format!("BED read error: {}", e)))?
    };

    Ok(GcBiasSection {
        reference,
        output_file,
        overwrite_output,
        window_size,
        window_stride,
        min_windows_per_bin,
        bed_table,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    #[test]
    fn test_config_requires_bam_file() {
        let yaml = write_temp("frag_length:\n  output_file: /tmp/x\n");
        assert!(RunConfiguration::from(&yaml.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_config_requires_at_least_one_model() {
        let bam = write_temp("dummy");
        let yaml = write_temp(&format!("bam_file: {}\n", bam.path().display()));
        let err = RunConfiguration::from(&yaml.path().to_path_buf()).unwrap_err();
        assert!(
            format!("{err}").contains("at least one model"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_config_gc_bias_requires_reference() {
        let bam = write_temp("dummy");
        let out = write_temp("");
        // gc_bias section without reference → error
        let yaml = write_temp(&format!(
            "bam_file: {}\ngc_bias:\n  output_file: {}\n  overwrite_output: true\n",
            bam.path().display(),
            out.path().display()
        ));
        let err = RunConfiguration::from(&yaml.path().to_path_buf()).unwrap_err();
        assert!(
            format!("{err}").contains("reference"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_config_frag_length_only_parses() {
        let bam = write_temp("dummy");
        let out_dir = tempfile::tempdir().unwrap();
        let out = out_dir.path().join("frag.json.gz");
        let yaml = write_temp(&format!(
            "bam_file: {}\nfrag_length:\n  output_file: {}\n  overwrite_output: true\n  min_reads: 5\n",
            bam.path().display(),
            out.display()
        ));
        let config = RunConfiguration::from(&yaml.path().to_path_buf()).unwrap();
        assert_eq!(config.bam_file, bam.path());
        assert!(config.gc_bias.is_none());
        let fl = config.frag_length.as_ref().unwrap();
        assert_eq!(fl.output_file, out);
        assert_eq!(fl.min_reads, 5);
    }

    #[test]
    fn test_config_both_sections_parse() {
        let bam = write_temp("dummy");
        let reference = write_temp(">chr1\nACGT\n");
        let out_dir = tempfile::tempdir().unwrap();
        let frag_out = out_dir.path().join("frag.json.gz");
        let gc_out = out_dir.path().join("gc.json.gz");
        let yaml = write_temp(&format!(
            "bam_file: {}\nmin_mapq: 20\nfrag_length:\n  output_file: {}\n  overwrite_output: true\ngc_bias:\n  reference: {}\n  output_file: {}\n  overwrite_output: true\n  window_size: 50\n  window_stride: 25\n  min_windows_per_bin: 3\n",
            bam.path().display(),
            frag_out.display(),
            reference.path().display(),
            gc_out.display(),
        ));
        let config = RunConfiguration::from(&yaml.path().to_path_buf()).unwrap();
        assert_eq!(config.min_mapq, 20);
        assert!(config.frag_length.is_some());
        let gc = config.gc_bias.as_ref().unwrap();
        assert_eq!(gc.window_size, 50);
        assert_eq!(gc.window_stride, 25);
        assert_eq!(gc.min_windows_per_bin, 3);
    }

    #[test]
    fn test_config_rejects_existing_output_when_overwrite_false() {
        let bam = write_temp("dummy");
        let existing = write_temp("already here");
        let yaml = write_temp(&format!(
            "bam_file: {}\nfrag_length:\n  output_file: {}\n  overwrite_output: false\n",
            bam.path().display(),
            existing.path().display()
        ));
        assert!(RunConfiguration::from(&yaml.path().to_path_buf()).is_err());
    }
}
