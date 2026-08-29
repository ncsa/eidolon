pub mod errors;
pub mod utils;
use eidolon_core::file_tools::bed_reader::{self};
use log::*;
use std::path::PathBuf;
use utils::filter_lib::{filter_fastq, filter_vcf};

use crate::filter_reads::{errors::FilterReadsError, utils::config::RunConfiguration};

pub fn main(config: &PathBuf) -> Result<(), FilterReadsError> {
    info!("////////////// Welcome to eidolon filter reads! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    info!("Using Configuration file input: {:?}", config);
    let run_config = RunConfiguration::from(config);
    // bed_file: path to bed file to use for filtering
    // file_to_filter: fastq or vcf ready for filtering, can be gzipped or not
    // is_gzip: whether the input file is gzipped or not
    // is_fastq: wether the
    // output_file: will be gzipped.
    let bed_table = bed_reader::read_bed(&run_config.bed_file, false)?;
    info!("Running filter-reads on input file.");
    for input_file in run_config.file_map.keys() {
        let (output_file, is_gzip, is_fastq) = run_config.file_map[input_file].clone();
        if is_fastq {
            info!("Filtering input fastq file: {:?}", input_file);
            filter_fastq(&bed_table, input_file, is_gzip, &output_file)?;
        } else {
            info!("Filtering input vcf file: {:?}", input_file);
            filter_vcf(&bed_table, input_file, is_gzip, &output_file)?;
        }
        info!(
            "Successfully filtered input file {:?}, written to {:?}",
            input_file, output_file
        );
    }
    Ok(())
}
