// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use crate::gen_mut_model::errors::GenMutationModelError;
use eidolon_core::{
    file_tools::{bed_reader::read_bed, vcf_tools::read_vcf_lean},
    structs::{bed_record::BedRecord, variants::Variant},
};
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::string::String;

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for this gen_mut_model run. It is built from
    // user-supplied input, as provided by the configuration yaml file.
    pub reference: PathBuf,
    pub mutations: HashMap<String, Vec<Variant>>,
    pub bed_table: HashMap<String, Vec<BedRecord>>,
    pub output_file: PathBuf,
    pub overwrite_output: bool,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, GenMutationModelError> {
        // Reads an input configuration file from yaml using the serde package. Then sets the
        // parameters based on the inputs. A "." value means to use the default value.
        //
        // Opens file for reading
        let f = fs::File::open(yml_file);
        let file = match f {
            Ok(l) => l,
            Err(error) => panic!(
                "Problem reading the config file: {:?}. System error: {}",
                yml_file, error,
            ),
        };
        // Uses serde_yml to read the file into a HashMap
        let scrape_config: HashMap<String, Value> =
            serde_yml::from_reader(file).expect("Error reading yaml file!");
        // Fill in the bed_file first, all hinges on that
        let reference = PathBuf::from(scrape_config["reference"].as_str().unwrap());
        if !reference.is_file() {
            panic!("Invalid reference file {:?}", reference)
        }
        let vcf_file = PathBuf::from(scrape_config["vcf_file"].as_str().unwrap());
        if !vcf_file.is_file() {
            panic!("Invalid bed file {:?}", vcf_file)
        }
        let bed_file_raw = scrape_config
            .get("bed_file")
            .and_then(|v| v.as_str())
            .unwrap_or(".");
        let bed_table = if bed_file_raw == "." {
            HashMap::new()
        } else {
            let bed_file = PathBuf::from(bed_file_raw);
            if !bed_file.is_file() {
                panic!("Invalid bed file {:?}", bed_file)
            }
            read_bed(&bed_file, false).expect("Error reading bed file!")
        };
        let overwrite_output = scrape_config["overwrite_output"].as_bool().unwrap_or(false);
        let output_file = PathBuf::from(scrape_config["output_file"].as_str().unwrap());
        if !overwrite_output && output_file.is_file() {
            panic!("Attempting to overwrite an existing file {:?}", output_file)
        }
        let mutations = read_vcf_lean(vcf_file)?;

        // `transition_matrix_file` set a context-free SNP matrix that generation never
        // read, so it was removed (#758). The old template ships the key with no value;
        // that is accepted. A value means someone expects an effect it never had.
        if let Some(path) = scrape_config
            .get("transition_matrix_file")
            .and_then(|v| v.as_str())
            .filter(|s| !s.trim().is_empty() && s.trim() != ".")
        {
            panic!(
                "transition_matrix_file ({path}) is no longer supported. It never affected \
                 generated variants: SNP alt bases come from the per-trinucleotide-context \
                 model fitted from the VCF. Remove the key to continue (see #758)."
            )
        }

        Ok(RunConfiguration {
            reference,
            mutations,
            bed_table,
            output_file,
            overwrite_output,
        })
    }
}

pub fn create_map_item(
    split_string: Vec<&str>,
    raw_string: &str,
    length: usize,
    filter_key: &str,
) -> (PathBuf, PathBuf) {
    let old_element = split_string[length - 3];
    let new_element = format!("{old_element}{filter_key}");
    let mut output_name = String::new();
    // stop one short of the end
    for i in 0..length - 1 {
        if i == length - 3 {
            output_name.push_str(&new_element);
            output_name.push('.');
        } else {
            output_name.push_str(split_string[i]);
            output_name.push('.');
        }
    }
    // End with the last extension with no traling dot.
    output_name.push_str(split_string[length - 1]);
    let temp_path = PathBuf::from(raw_string);
    if !temp_path.is_file() {
        panic!("Input file not found! {:?}", temp_path)
    }
    (temp_path, PathBuf::from(output_name))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_run_configuration() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = format!("{}/test_data/references/H1N1.fa", manifest_dir);
        let vcf_file = format!("{}/test_data/vcfs/small_snps.vcf", manifest_dir);
        let output_file = format!("{}/test_data/test_run_config_output.json.gz", manifest_dir);
        let yaml = format!(
            "reference: {}\nvcf_file: {}\noutput_file: {}\nbed_file: .\noverwrite_output: true\n",
            reference, vcf_file, output_file
        );
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        let config = RunConfiguration::from(&tmp.path().to_path_buf()).unwrap();
        let h1n1_variants = config.mutations.get("H1N1_HA").expect("H1N1_HA not found");
        assert_eq!(h1n1_variants.len(), 3);
        assert!(config.bed_table.is_empty());
    }

    fn config_with(extra: &str) -> Result<RunConfiguration, GenMutationModelError> {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let out = tempfile::tempdir().unwrap();
        let yaml = format!(
            "reference: {manifest_dir}/test_data/references/H1N1.fa\n\
             vcf_file: {manifest_dir}/test_data/vcfs/small_snps.vcf\n\
             output_file: {}\n\
             overwrite_output: true\n{extra}",
            out.path().join("model.json.gz").display()
        );
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        RunConfiguration::from(&tmp.path().to_path_buf())
    }

    /// Configs copied from the old template carry `transition_matrix_file:` with no
    /// value. Those must keep working after the option's removal.
    #[test]
    fn an_empty_transition_matrix_file_key_is_accepted() {
        assert!(config_with("transition_matrix_file:\n").is_ok());
        assert!(config_with("transition_matrix_file: .\n").is_ok());
    }

    /// A value means the user expects the matrix to shape the model. It never did,
    /// so the run must stop rather than silently ignore it.
    #[test]
    #[should_panic(expected = "transition_matrix_file (/any/matrix.tsv) is no longer supported")]
    fn a_set_transition_matrix_file_is_refused() {
        let _ = config_with("transition_matrix_file: /any/matrix.tsv\n");
    }
}
