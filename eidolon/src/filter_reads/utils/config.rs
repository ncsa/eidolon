// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.

use log::*;
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::string::String;

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
    pub fn from(yml_file: &PathBuf) -> Self {
        // Reads an input configuration file from yaml using the serde package. Then sets the
        // parameters based on the inputs. A "." value means to use the default value.

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
        let bed_file = PathBuf::from(scrape_config["bed_file"].as_str().unwrap());
        if !bed_file.is_file() {
            panic!("Invalid bed file {:?}", bed_file)
        }
        let files_to_filter_raw = scrape_config["files_to_filter"].as_sequence().unwrap();
        let filter_key: &str = {
            let temp_key = scrape_config["filter_key"].as_str().unwrap();
            if temp_key == "." { "_filter" } else { temp_key }
        };
        let overwrite_output = scrape_config["overwrite_output"].as_bool().unwrap();
        info!("Overwrite? {overwrite_output}");
        let file_map: HashMap<PathBuf, (PathBuf, bool, bool)> = {
            let mut temp_map = HashMap::new();
            for raw_value in files_to_filter_raw {
                let mut is_gzip = false;
                let mut is_fastq = true;
                let raw_string = raw_value.as_str().unwrap();
                let split_string: Vec<&str> = raw_string.split(".").collect();
                let elem_len = split_string.len();
                let ext1 = split_string[elem_len - 1];
                if ext1 == "gz" {
                    is_gzip = true;
                    let ext2 = split_string[elem_len - 2];
                    if ext2 == "vcf" {
                        is_fastq = false;
                    } else if ext2 == "fastq" {
                        // Do nothing
                    } else {
                        panic!("Unknown File Extension! {:?}", raw_string)
                    }
                } else if ext1 == "vcf" {
                    is_fastq = false;
                } else if ext1 == "fastq" {
                    // do nothing
                } else {
                    panic!("Unknown File Extension! {:?}", raw_string)
                }
                let (input_path, output_path) =
                    create_map_item(split_string, raw_string, elem_len, filter_key);
                if !overwrite_output && output_path.is_file() {
                    panic!("Attempting to overwrite an existing file {:?}", output_path)
                }
                temp_map.insert(input_path, (output_path, is_gzip, is_fastq));
            }
            temp_map
        };

        RunConfiguration { bed_file, file_map }
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
    // use super::*;

    #[test]
    fn test_run_configuration() {
        // TODO: real assertions — stub kept so the test name is reserved.
    }
}
