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
    /// Optional path to a 4×4 TSV specifying a custom SNP transition matrix.
    /// Rows/columns are A/C/G/T. A single header line is ignored.
    /// Overrides the transition matrix inferred from VCF data.
    pub transition_matrix_file: Option<PathBuf>,
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

        let transition_matrix_file = scrape_config
            .get("transition_matrix_file")
            .and_then(|v| v.as_str())
            .map(|s| {
                let p = PathBuf::from(s);
                if !p.is_file() {
                    panic!("transition_matrix_file not found: {:?}", p)
                }
                p
            });

        Ok(RunConfiguration {
            reference,
            mutations,
            bed_table,
            output_file,
            overwrite_output,
            transition_matrix_file,
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
}
