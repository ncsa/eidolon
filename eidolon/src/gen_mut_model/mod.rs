/// This file will create a mutation model based on input data.
///
pub mod errors;
pub mod utils;

use crate::gen_mut_model::{
    errors::GenMutationModelError, utils::config::RunConfiguration, utils::runner::runner,
};
use eidolon_core::structs::{bed_record::BedRecord, variants::Variant};
use log::*;
use std::{collections::HashMap, path::PathBuf};

pub fn main(config_file: &PathBuf) -> Result<(), GenMutationModelError> {
    let run_config = RunConfiguration::from(config_file)?;
    // filter the variants with the bed_table, if applicable
    let filtered_mutations = if run_config.bed_table.is_empty() {
        run_config.mutations
    } else {
        filter_mutations_by_bed(&run_config.mutations, &run_config.bed_table)
    };
    runner(
        &run_config.reference,
        filtered_mutations,
        run_config.bed_table,
        &run_config.output_file,
    )?;

    Ok(())
}

/// Keep the variants that fall inside `bed_table`. `Variant::location` is the
/// VCF's 1-based POS (see `read_vcf_lean`), so membership goes through
/// [`BedRecord::contains_vcf_pos`]. Contigs absent from the BED are dropped.
fn filter_mutations_by_bed(
    mutations: &HashMap<String, Vec<Variant>>,
    bed_table: &HashMap<String, Vec<BedRecord>>,
) -> HashMap<String, Vec<Variant>> {
    let mut filtered: HashMap<String, Vec<Variant>> = HashMap::new();
    for (chrom, contig_mutations) in mutations {
        let Some(search_region) = bed_table.get(chrom) else {
            // Throw out all mutations in the vcf with this chrom
            debug!("Filtered out contig: {chrom}");
            continue;
        };
        for variant in contig_mutations {
            if search_region
                .iter()
                .any(|record| record.contains_vcf_pos(chrom.as_str(), variant.location))
            {
                filtered
                    .entry(chrom.clone())
                    .or_default()
                    .push(variant.clone());
            }
        }
    }
    filtered
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::file_tools::vcf_tools::read_vcf_lean;
    use std::io::Write;

    /// #770: read through the real VCF reader so the test pins the reader's
    /// coordinate convention too. BED `chr1 100 200` covers POS 101..=200.
    #[test]
    fn bed_filter_edges_use_one_based_pos() {
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(
            tmp,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG\t60\tPASS\t.\n\
             chr1\t101\t.\tA\tG\t60\tPASS\t.\n\
             chr1\t200\t.\tA\tG\t60\tPASS\t.\n\
             chr1\t201\t.\tA\tG\t60\tPASS\t.\n\
             chr2\t150\t.\tA\tG\t60\tPASS\t.\n"
        )
        .unwrap();
        let mutations = read_vcf_lean(tmp.path().to_path_buf()).unwrap();
        assert_eq!(mutations["chr1"].len(), 4, "fixture must parse fully");
        let bed = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 100, 200).unwrap()],
        )]);
        let filtered = filter_mutations_by_bed(&mutations, &bed);
        let locs: Vec<usize> = filtered["chr1"].iter().map(|v| v.location).collect();
        assert_eq!(locs, vec![101, 200]);
        assert!(
            !filtered.contains_key("chr2"),
            "contig absent from BED is dropped"
        );
    }
}
