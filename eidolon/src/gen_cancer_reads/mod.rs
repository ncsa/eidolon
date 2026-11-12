pub mod errors;
pub mod utils;

use std::path::PathBuf;

use errors::GenCancerReadsError;
use log::info;

use crate::gen_cancer_reads::utils::config::CancerConfig;
use crate::gen_cancer_reads::utils::runner::run_cancer;

/// gen-cancer-reads simulates a tumor/normal mixture: two gen-reads passes over
/// the same reference (normal at `(1-purity)·C`, tumor at `purity·C` sharing the
/// normal's germline), merged into tagged FASTQs + an origin-tagged truth VCF.
/// Native port of `tools/cancer_simulate.sh` (no bcftools/awk).
pub fn main(config: &PathBuf) -> Result<(), GenCancerReadsError> {
    info!("////////////// Welcome to eidolon cancer simulator! \\\\\\\\\\\\\\\\\\\\\\\\");
    let cfg = CancerConfig::from_yaml_file(config)?;

    if !cfg.overwrite_output {
        for name in [
            format!("{}_merged_r1.fastq.gz", cfg.output_basename),
            format!("{}_merged_truth.vcf.gz", cfg.output_basename),
        ] {
            let p = cfg.output_dir.join(&name);
            if p.is_file() {
                return Err(GenCancerReadsError::ConfigError(format!(
                    "refusing to overwrite existing {p:?}; set overwrite_output: true"
                )));
            }
        }
    }

    run_cancer(&cfg)
}
