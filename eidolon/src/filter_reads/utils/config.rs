// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.

use crate::config_keys::check_keys;
use crate::filter_reads::errors::FilterReadsError;
use log::*;
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::string::String;

/// Every top-level key the config parser reads. Anything else is rejected (#496).
pub const KNOWN_KEYS: &[&str] = &[
    "bed_file",
    "files_to_filter",
    "filter_key",
    "overwrite_output",
];

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for this filter-reads run. It is built from user-supplied input
    // in the form of a configuration yaml file
    //
    // bed_file: The path to the bed_file for the run.
    // files_to_filter: The list of files to filter.
    // filter_key: The key to add to the filtered file names so you know they have been filtered.
    pub bed_file: PathBuf,
    pub file_map: HashMap<PathBuf, (PathBuf, bool, bool)>,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, FilterReadsError> {
        // Reads an input configuration file from yaml using the serde package. Then sets the
        // parameters based on the inputs. A "." value means to use the default value.
        let config_err = FilterReadsError::ConfigurationError;
        let file = fs::File::open(yml_file).map_err(|error| {
            config_err(format!("cannot read the config file {yml_file:?}: {error}"))
        })?;
        // Uses serde_yml to read the file into a HashMap
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)
            .map_err(|error| config_err(format!("{yml_file:?} is not valid YAML: {error}")))?;
        check_keys(
            scrape_config.keys().map(String::as_str),
            KNOWN_KEYS,
            "filter-reads",
        )
        .map_err(config_err)?;
        let required_str = |key: &str| -> Result<&str, FilterReadsError> {
            scrape_config
                .get(key)
                .and_then(|v| v.as_str())
                .ok_or_else(|| config_err(format!("missing required key '{key}'")))
        };
        let bed_file = PathBuf::from(required_str("bed_file")?);
        if !bed_file.is_file() {
            return Err(FilterReadsError::FileNotFound(format!(
                "bed_file {bed_file:?}"
            )));
        }
        let files_to_filter_raw = scrape_config
            .get("files_to_filter")
            .and_then(|v| v.as_sequence())
            .ok_or_else(|| config_err("missing required key 'files_to_filter' (a list)".into()))?;
        let filter_key: &str = {
            let temp_key = required_str("filter_key")?;
            if temp_key == "." { "_filter" } else { temp_key }
        };
        let overwrite_output = scrape_config
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .ok_or_else(|| {
                config_err("missing required key 'overwrite_output' (true/false)".into())
            })?;
        info!("Overwrite? {overwrite_output}");
        let mut file_map: HashMap<PathBuf, (PathBuf, bool, bool)> = HashMap::new();
        for raw_value in files_to_filter_raw {
            let raw_string = raw_value.as_str().ok_or_else(|| {
                config_err(format!("files_to_filter entry {raw_value:?} is not a path"))
            })?;
            let input_path = PathBuf::from(raw_string);
            let (output_path, is_gzip, is_fastq) = filtered_output(&input_path, filter_key)
                .ok_or_else(|| {
                    FilterReadsError::MalformedFileName(format!(
                        "{raw_string:?}: expected one of .fastq, .fastq.gz, .vcf, .vcf.gz"
                    ))
                })?;
            if !input_path.is_file() {
                return Err(FilterReadsError::FileNotFound(format!("{input_path:?}")));
            }
            if !overwrite_output && output_path.is_file() {
                return Err(FilterReadsError::OverwriteFileError(format!(
                    "{output_path:?}; set overwrite_output: true to replace it"
                )));
            }
            file_map.insert(input_path, (output_path, is_gzip, is_fastq));
        }

        Ok(RunConfiguration { bed_file, file_map })
    }

    pub fn log(&mut self) {
        info!(
            "Using bed file to filter the input files: {:?}",
            self.bed_file
        );
        for key in self.file_map.keys() {
            info!("Filtering: {:?}", key);
            info!("Producing: {:?}", self.file_map[key].0);
        }
    }
}

/// The suffixes filter-reads accepts: `(suffix, is_gzip, is_fastq)`.
const INPUT_SUFFIXES: [(&str, bool, bool); 4] = [
    (".fastq.gz", true, true),
    (".fastq", false, true),
    (".vcf.gz", true, false),
    (".vcf", false, false),
];

/// Where the filtered copy of `input` goes, and whether `input` is gzipped and a FASTQ.
///
/// The name is built from the file name alone, so a dot in a directory cannot shift it: the
/// filter key goes between the stem and the suffix, in the input's directory. The writer
/// always compresses, so the output always ends in `.gz`, including for plain input:
/// `reads.fastq` becomes `reads<key>.fastq.gz`. Returns `None` for any other suffix, or for a
/// name that is nothing but a suffix (#824).
pub fn filtered_output(input: &Path, filter_key: &str) -> Option<(PathBuf, bool, bool)> {
    let name = input.file_name()?.to_str()?;
    let &(suffix, is_gzip, is_fastq) = INPUT_SUFFIXES
        .iter()
        .find(|(suffix, _, _)| name.len() > suffix.len() && name.ends_with(suffix))?;
    let stem = &name[..name.len() - suffix.len()];
    let out_suffix = if is_gzip {
        suffix.to_string()
    } else {
        format!("{suffix}.gz")
    };
    Some((
        input.with_file_name(format!("{stem}{filter_key}{out_suffix}")),
        is_gzip,
        is_fastq,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Write an empty file at `dir/name` and return its path as a String.
    fn touch(dir: &std::path::Path, name: &str) -> String {
        let path = dir.join(name);
        fs::write(&path, "").unwrap();
        path.to_str().unwrap().to_string()
    }

    #[test]
    fn filtered_output_inserts_the_key_before_the_suffix() {
        for (input, expected, is_gzip, is_fastq) in [
            ("/d/sample.fastq.gz", "/d/sample_flt.fastq.gz", true, true),
            ("/d/calls.vcf.gz", "/d/calls_flt.vcf.gz", true, false),
            // Plain input is written compressed, so its output name gains `.gz`.
            ("/d/sample.fastq", "/d/sample_flt.fastq.gz", false, true),
            ("/d/calls.vcf", "/d/calls_flt.vcf.gz", false, false),
            // A dot in a directory, or extra dots in the stem, must not move the key.
            (
                "/data/run.1/reads.fastq",
                "/data/run.1/reads_flt.fastq.gz",
                false,
                true,
            ),
            (
                "/d/sample.v2.fastq.gz",
                "/d/sample.v2_flt.fastq.gz",
                true,
                true,
            ),
            ("reads.vcf.gz", "reads_flt.vcf.gz", true, false),
        ] {
            assert_eq!(
                filtered_output(Path::new(input), "_flt"),
                Some((PathBuf::from(expected), is_gzip, is_fastq)),
                "{input}"
            );
        }
    }

    /// Must not fire: anything but the four suffixes is refused, as is a bare suffix.
    #[test]
    fn filtered_output_refuses_other_suffixes() {
        for input in [
            "/d/reads.fq.gz",
            "/d/reads.bam",
            "/d/reads.fastq.bz2",
            "/d/.fastq",
            "/d/vcf",
        ] {
            assert_eq!(filtered_output(Path::new(input), "_flt"), None, "{input}");
        }
    }

    #[test]
    fn run_configuration_maps_each_file_with_its_flags_and_default_key() {
        let temp_dir = tempfile::tempdir().unwrap();
        let dir = temp_dir.path().to_str().unwrap();
        let bed = touch(temp_dir.path(), "regions.bed");
        let fastq = touch(temp_dir.path(), "reads.fastq.gz");
        let vcf = touch(temp_dir.path(), "calls.vcf.gz");
        let yml = temp_dir.path().join("config.yml");
        // "." selects the default filter key, "_filter".
        fs::write(
            &yml,
            format!(
                "bed_file: {bed}\nfiles_to_filter:\n  - {fastq}\n  - {vcf}\n\
                 filter_key: .\noverwrite_output: false\n"
            ),
        )
        .unwrap();

        let config = RunConfiguration::from(&yml).unwrap();

        assert_eq!(config.bed_file, PathBuf::from(&bed));
        assert_eq!(config.file_map.len(), 2);
        // (output path, is_gzip, is_fastq)
        assert_eq!(
            config.file_map[&PathBuf::from(&fastq)],
            (
                PathBuf::from(format!("{dir}/reads_filter.fastq.gz")),
                true,
                true
            )
        );
        assert_eq!(
            config.file_map[&PathBuf::from(&vcf)],
            (
                PathBuf::from(format!("{dir}/calls_filter.vcf.gz")),
                true,
                false
            )
        );
    }
}
