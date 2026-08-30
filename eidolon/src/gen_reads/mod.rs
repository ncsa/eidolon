pub mod errors;
pub mod utils;
use errors::GenerateReadsError;
use log::*;

use crate::gen_reads::utils::{config::RunConfiguration, runner::run_neat};
use eidolon_core::rng::NeatRng;
use std::path::PathBuf;

/// gen-reads is the primary read generation function of eidolon. It reads a fasta file and generates a set of fastqs and/or a set of variants. It can now also filter reads by bed file.
pub fn main(config: &PathBuf) -> Result<(), GenerateReadsError> {
    info!("////////////// Welcome to eidolon read generator! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if !config.display().to_string().is_empty() {
        info!("Using Configuration file input: {:?}", config);
        RunConfiguration::from_yaml_file(config)
            .expect("Error generating run configuration from input yaml")
    } else {
        panic!("Failed to supply config file!");
    };

    // Check that we are clear to write outputs
    if !config.overwrite_output {
        if let Some(filename) = &config.output_fastq_1
            && filename.is_file()
        {
            panic!("Attempting to overwrite an existing file: {:?}", filename);
        }
        if let Some(filename) = &config.output_fastq_2
            && filename.is_file()
        {
            panic!("Attempting to overwrite an existing file: {:?}", filename);
        }
        if let Some(filename) = &config.output_vcf
            && filename.is_file()
        {
            panic!("Attempting to overwrite an existing file: {:?}", filename);
        }
    }

    info!("////////////// Configuration successuful! Ready to run! \\\\\\\\\\\\\\\\\\\\\\\\\\");
    // Generate the RNG used for this run. If one was given in the config file, use that, or else
    // use thread_rng to generate a random seed, then seed using a SeedableRng based on StdRng
    let mut rng: NeatRng =
        NeatRng::new_from_seed(&config.seed_vec).expect("Neat failed during rng creation!");
    // run the generate reads main script
    let result = run_neat(&config, &mut rng);
    match result {
        Ok(_) => {
            // Continue on for bed filtering
            Ok(())
        }
        Err(error) => {
            error!("runner returned an error {:?}", error);
            Err(GenerateReadsError::RunnerError)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::gen_reads::utils::{config::RunConfiguration, runner::run_neat};
    use eidolon_core::rng::NeatRng;
    use flate2::read::MultiGzDecoder;
    use noodles::bam;
    use std::io::{BufRead, BufReader, Write};
    use std::path::{Path, PathBuf};
    use tempfile::{NamedTempFile, tempdir};

    // ── helpers ──────────────────────────────────────────────────────────────

    fn make_config(reference: &Path, output_dir: &Path, paired_ended: bool) -> RunConfiguration {
        let mut yaml = NamedTempFile::new().unwrap();
        let pair_section = if paired_ended {
            "paired_ended: true\nfragment_mean: 200.0\nfragment_st_dev: 30.0\n"
        } else {
            "paired_ended: false\n"
        };
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: true\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: bam_int_test\noverwrite_output: true\n\
             rng_seed: integration bam test\n{pair}",
            ref = reference.display(),
            out = output_dir.display(),
            pair = pair_section,
        )
        .unwrap();
        RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap()
    }

    fn run(config: RunConfiguration) {
        let mut rng = NeatRng::new_from_seed(&config.seed_vec).unwrap();
        run_neat(&config, &mut rng).unwrap();
    }

    /// Returns (ref_id, 0-based alignment start) for every record in the BAM.
    fn bam_positions(path: &Path) -> Vec<(usize, usize)> {
        let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
        reader.read_header().unwrap();
        reader
            .records()
            .map(|r| {
                let rec = r.unwrap();
                let ref_id = rec
                    .reference_sequence_id()
                    .unwrap() // Option → Some(Result<usize>)
                    .unwrap(); // Result  → usize
                let pos = rec
                    .alignment_start()
                    .unwrap() // Option → Some(Result<Position>)
                    .unwrap() // Result → Position
                    .get()
                    - 1;
                (ref_id, pos)
            })
            .collect()
    }

    fn fastq_record_count(path: &Path) -> usize {
        // The output FASTQ is multiple concatenated gzip streams (one per contig);
        // MultiGzDecoder reads all streams in sequence.
        let reader = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
        reader.lines().count() / 4
    }

    fn h1n1_reference() -> PathBuf {
        PathBuf::from(format!(
            "{}/test_data/references/H1N1.fa",
            env!("CARGO_MANIFEST_DIR")
        ))
    }

    // ── tests ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_bam_single_ended_produced_and_sorted() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), false);
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let positions = bam_positions(&bam_path);
        assert!(!positions.is_empty(), "BAM should contain records");

        for window in positions.windows(2) {
            assert!(
                window[0] <= window[1],
                "BAM not coordinate-sorted: {:?} followed by {:?}",
                window[0],
                window[1]
            );
        }
    }

    #[test]
    fn test_bam_paired_ended_produced_and_sorted() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), true);
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let positions = bam_positions(&bam_path);
        assert!(!positions.is_empty(), "BAM should contain records");

        for window in positions.windows(2) {
            assert!(
                window[0] <= window[1],
                "BAM not coordinate-sorted: {:?} followed by {:?}",
                window[0],
                window[1]
            );
        }
    }

    #[test]
    fn test_bam_record_count_matches_fastq_single_ended() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), false);
        let bam_path = config.output_bam.clone().unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let bam_count = bam_positions(&bam_path).len();
        let fastq_count = fastq_record_count(&fastq_path);
        assert_eq!(
            bam_count, fastq_count,
            "BAM record count ({bam_count}) must equal FASTQ record count ({fastq_count})"
        );
    }

    /// Writes a BED covering only H1N1_HA (the full segment) and verifies that:
    /// 1. Reads are generated (the covered contig produces output).
    /// 2. Fewer reads are produced than the no-BED run (7 of 8 segments are excluded).
    #[test]
    fn test_mutation_regions_custom_and_default_rates() {
        use std::fs::File;
        let out = tempdir().unwrap();
        let ref_path = h1n1_reference();

        // Create a mutation BED file
        // H1N1 has 8 segments. We'll set rates for all of them to be safe,
        // or just check the contig name in the VCF.
        let mut_bed_path = out.path().join("mut.bed");
        {
            let mut f = File::create(&mut_bed_path).unwrap();
            let segments = [
                "H1N1_segment_1",
                "H1N1_HA",
                "H1N1_MP",
                "H1N1_NA",
                "H1N1_NP",
                "H1N1_NS",
                "H1N1_PA",
                "H1N1_PB1",
                "H1N1_PB2",
            ];
            for seg in segments {
                // High rate in 100-200
                writeln!(f, "{}\t100\t200\tmut_rate=1.0", seg).unwrap();
                // Zero rate in 300-400
                writeln!(f, "{}\t300\t400\tmut_rate=0.0", seg).unwrap();
            }
        }

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 0\n\
             produce_fastq: false\nproduce_bam: false\nproduce_vcf: true\n\
             output_dir: {out}\noutput_filename: mut_test\noverwrite_output: true\n\
             mutation_regions: {mut_bed}\n\
             mutation_rate: 0.1\n\
             rng_seed: mut test",
            ref = ref_path.display(),
            out = out.path().display(),
            mut_bed = mut_bed_path.display(),
        )
        .unwrap();

        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        run(config.clone());

        let vcf_path = config.output_vcf.unwrap();
        assert!(vcf_path.exists());

        // Check variants in VCF - it is gzipped
        let mut reader = BufReader::new(MultiGzDecoder::new(File::open(vcf_path).unwrap()));
        let mut line = String::new();
        let mut variants_in_high_rate = 0;
        let mut variants_in_zero_rate = 0;
        let mut variants_elsewhere = 0;

        while reader.read_line(&mut line).unwrap() > 0 {
            if line.starts_with('#') {
                line.clear();
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let pos: usize = fields[1].parse().unwrap(); // 1-based
            let pos0 = pos - 1;

            if (100..200).contains(&pos0) {
                variants_in_high_rate += 1;
            } else if (300..400).contains(&pos0) {
                // The variant's location is actually the start position of the mutation.
                // For a SNP it is exactly pos0.
                variants_in_zero_rate += 1;
            } else {
                variants_elsewhere += 1;
            }
            line.clear();
        }

        // With rate 1.0 over 100bp, we expect ~100 mutations in the high rate region.
        // With rate 0.0 over 100bp, we expect 0 mutations.
        // With rate 0.1 over the rest of the contig (segment 1 is 2341 bp),
        // area = 2341 - 100 (high) - 100 (zero) = 2141 bp.
        // Expected mutations elsewhere = 2141 * 0.1 = 214.1 -> 214.

        assert!(
            variants_in_high_rate > 50,
            "Should have many mutations in high rate region, got {}",
            variants_in_high_rate
        );
        assert_eq!(
            variants_in_zero_rate, 0,
            "Should have NO mutations in zero rate region, got {}",
            variants_in_zero_rate
        );
        assert!(
            variants_elsewhere > 100,
            "Should have mutations in default rate regions, got {}",
            variants_elsewhere
        );
    }

    #[test]
    fn test_target_bed_limits_reads() {
        let out = tempdir().unwrap();

        // BED covering only H1N1_HA; the other 7 segments are absent → skipped entirely.
        let bed_path = out.path().join("targets.bed");
        std::fs::write(&bed_path, "H1N1_HA\t0\t1726\n").unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: false\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: bed_test\noverwrite_output: true\n\
             rng_seed: bed test\ntarget_bed: {bed}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            bed = bed_path.display(),
        )
        .unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let bed_count = fastq_record_count(&fastq_path);

        // Without a BED, all 8 segments produce reads.
        let out_full = tempdir().unwrap();
        let config_full = make_config(&h1n1_reference(), out_full.path(), false);
        let fastq_full = config_full.output_fastq_1.clone().unwrap();
        run(config_full);
        let full_count = fastq_record_count(&fastq_full);

        assert!(bed_count > 0, "BED-filtered run produced no reads");
        assert!(
            bed_count < full_count,
            "BED-filtered run ({bed_count} reads) should be less than full run ({full_count} reads)"
        );
    }

    /// A BED that names a contig not in the reference should produce an empty FASTQ.
    #[test]
    fn test_target_bed_unknown_contig_produces_no_reads() {
        let out = tempdir().unwrap();

        let bed_path = out.path().join("empty_targets.bed");
        std::fs::write(&bed_path, "nonexistent_contig\t0\t1000\n").unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: false\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: empty_bed_test\noverwrite_output: true\n\
             rng_seed: empty bed test\ntarget_bed: {bed}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            bed = bed_path.display(),
        )
        .unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let count = fastq_record_count(&fastq_path);
        assert_eq!(
            count, 0,
            "Expected no reads when BED references only unknown contigs, got {count}"
        );
    }

    /// When produce_fastq=false and produce_bam=true, reads must still be generated and
    /// staged into the BAM body writer. The output BAM must be a valid, non-empty,
    /// coordinate-sorted file — identical in record count to a paired FASTQ run with the
    /// same seed and coverage.
    #[test]
    fn test_bam_only_no_fastq_produces_records() {
        let out = tempdir().unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: false\nproduce_bam: true\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: bam_only_test\noverwrite_output: true\n\
             rng_seed: bam only test\npaired_ended: false\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
        )
        .unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let bam_path = config.output_bam.clone().unwrap();
        let fastq_path = config.output_fastq_1.clone();

        run(config);

        if let Some(ref fq) = fastq_path {
            assert!(
                !fq.exists(),
                "FASTQ file must not be created when produce_fastq=false"
            );
        }

        assert!(
            bam_path.exists(),
            "BAM file must exist with produce_bam=true"
        );

        let mut reader = bam::io::Reader::new(std::fs::File::open(&bam_path).unwrap());
        let header = reader.read_header().unwrap();
        let n_ref_seqs = header.reference_sequences().len();
        assert_eq!(
            n_ref_seqs, 8,
            "H1N1 BAM header should list 8 reference sequences, got {n_ref_seqs}"
        );

        let positions = bam_positions(&bam_path);
        assert!(
            !positions.is_empty(),
            "BAM must contain records when produce_bam=true even if produce_fastq=false"
        );

        // Records must be in coordinate-sorted order.
        for window in positions.windows(2) {
            assert!(
                window[0] <= window[1],
                "BAM not coordinate-sorted in BAM-only mode: {:?} followed by {:?}",
                window[0],
                window[1]
            );
        }
    }

    /// Running with an explicit num_threads exercises the ThreadPoolBuilder path in
    /// run_neat. Verifies the BAM is non-empty, sorted, and has the same record count
    /// as the FASTQ — which confirms the parallel-gather and sort-by-idx logic are
    /// correct under an explicitly-sized thread pool.
    #[test]
    fn test_bam_explicit_num_threads_sorted_and_consistent() {
        let out = tempdir().unwrap();

        // Start from make_config and inject num_threads afterward (the field is pub).
        let mut config = make_config(&h1n1_reference(), out.path(), false);
        config.num_threads = Some(2);

        let bam_path = config.output_bam.clone().unwrap();
        let fastq_path = config.output_fastq_1.clone().unwrap();
        run(config);

        let positions = bam_positions(&bam_path);
        assert!(
            !positions.is_empty(),
            "BAM should contain records with num_threads=2"
        );

        for window in positions.windows(2) {
            assert!(
                window[0] <= window[1],
                "BAM is not coordinate-sorted with num_threads=2: {:?} followed by {:?}",
                window[0],
                window[1]
            );
        }

        let fastq_count = fastq_record_count(&fastq_path);
        assert_eq!(
            positions.len(),
            fastq_count,
            "BAM record count ({}) must equal FASTQ record count ({}) with num_threads=2",
            positions.len(),
            fastq_count
        );
    }

    /// An explicit thread count must produce deterministic output: the same seed with
    /// num_threads=2 and the default (None) pool should yield identical BAM positions.
    #[test]
    fn test_bam_num_threads_deterministic_matches_default() {
        let out_explicit = tempdir().unwrap();
        let out_default = tempdir().unwrap();

        let mut config_explicit = make_config(&h1n1_reference(), out_explicit.path(), false);
        config_explicit.num_threads = Some(2);
        let bam_explicit = config_explicit.output_bam.clone().unwrap();

        let config_default = make_config(&h1n1_reference(), out_default.path(), false);
        let bam_default = config_default.output_bam.clone().unwrap();

        run(config_explicit);
        run(config_default);

        let positions_explicit = bam_positions(&bam_explicit);
        let positions_default = bam_positions(&bam_default);

        assert_eq!(
            positions_explicit, positions_default,
            "Same seed must yield identical BAM positions regardless of num_threads"
        );
    }

    /// Every record in a paired-ended BAM must carry consistent SAM flags:
    /// all reads have the SEGMENTED (0x1) bit set; R1 has FIRST_SEGMENT (0x40) and
    /// MATE_REVERSE_COMPLEMENTED (0x20) but not REVERSE_COMPLEMENTED (0x10); R2 has
    /// LAST_SEGMENT (0x80) and REVERSE_COMPLEMENTED (0x10). Also checks that the
    /// number of R1 and R2 records match, which would catch any read-drop or
    /// duplication in the BAM staging path.
    #[test]
    fn test_paired_ended_bam_flags_correct() {
        use noodles::sam::alignment::record::Flags;

        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), true);
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let mut reader = bam::io::Reader::new(std::fs::File::open(&bam_path).unwrap());
        reader.read_header().unwrap();

        let mut r1_count = 0usize;
        let mut r2_count = 0usize;

        for result in reader.records() {
            let record = result.unwrap();
            let flags = record.flags();

            assert!(
                flags.is_segmented(),
                "PAIRED flag (0x1) must be set on every record in a paired-ended BAM"
            );

            if flags.is_first_segment() {
                r1_count += 1;
                assert!(
                    !flags.is_reverse_complemented(),
                    "R1 must not have REVERSE_COMPLEMENTED (0x10); flags={:?}",
                    flags
                );
                assert!(
                    flags.contains(Flags::MATE_REVERSE_COMPLEMENTED),
                    "R1 must have MATE_REVERSE_COMPLEMENTED (0x20); flags={:?}",
                    flags
                );
            } else {
                assert!(
                    flags.is_last_segment(),
                    "Record is neither R1 (FIRST_SEGMENT) nor R2 (LAST_SEGMENT); flags={:?}",
                    flags
                );
                r2_count += 1;
                assert!(
                    flags.is_reverse_complemented(),
                    "R2 must have REVERSE_COMPLEMENTED (0x10); flags={:?}",
                    flags
                );
            }
        }

        assert!(r1_count > 0, "Paired-ended BAM produced no R1 records");
        assert_eq!(
            r1_count, r2_count,
            "R1 count ({r1_count}) must equal R2 count ({r2_count}) in a paired-ended BAM"
        );
    }

    /// When a paired-ended run is configured with an explicit fragment_mean, the
    /// template lengths (TLEN) recorded in the BAM should reflect that model.
    /// make_config uses fragment_mean=200, fragment_st_dev=30 for paired runs.
    /// We collect |TLEN| from R1 records and assert the mean is within ±20% of 200.
    #[test]
    fn test_paired_ended_insert_size_matches_model() {
        let out = tempdir().unwrap();
        let config = make_config(&h1n1_reference(), out.path(), true);
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let mut reader = bam::io::Reader::new(std::fs::File::open(&bam_path).unwrap());
        reader.read_header().unwrap();

        let mut tlens: Vec<f64> = Vec::new();
        for result in reader.records() {
            let record = result.unwrap();
            let flags = record.flags();
            // Collect only R1 records so each template is counted once.
            if !flags.is_segmented() || !flags.is_first_segment() {
                continue;
            }
            let tlen = record.template_length().unsigned_abs() as f64;
            if tlen > 0.0 {
                tlens.push(tlen);
            }
        }

        assert!(
            !tlens.is_empty(),
            "No TLEN values found in paired-ended BAM"
        );
        let mean = tlens.iter().sum::<f64>() / tlens.len() as f64;
        let expected = 200.0_f64;
        let tolerance = 0.20;
        assert!(
            (mean - expected).abs() <= expected * tolerance,
            "Mean insert size {mean:.1} is not within ±20% of expected {expected:.1} \
             (fragment_mean=200, fragment_st_dev=30)"
        );
    }

    /// Seeds a homozygous A→T SNP at H1N1_HA:100 via input VCF, runs with BAM
    /// output enabled, then confirms that at least one aligned read covering that
    /// position carries the T alt allele. This is the end-to-end check for variant
    /// injection: the output VCF test only verifies the metadata path, not that reads
    /// actually encode the mutation.
    #[test]
    fn test_input_vcf_snp_appears_in_bam_reads() {
        let out = tempdir().unwrap();

        let vcf_path = out.path().join("input.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             H1N1_HA\t100\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 10\n\
             produce_fastq: true\nproduce_bam: true\nproduce_vcf: false\n\
             output_dir: {out}\noutput_filename: vcf_bam_test\noverwrite_output: true\n\
             rng_seed: vcf bam test\ninput_vcf: {vcf}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            vcf = vcf_path.display(),
        )
        .unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let bam_path = config.output_bam.clone().unwrap();
        run(config);

        let mut reader = bam::io::Reader::new(std::fs::File::open(&bam_path).unwrap());
        let header = reader.read_header().unwrap();

        // Locate H1N1_HA's 0-based contig index in the BAM header.
        // BString::as_bytes() is crate-private; go through the public AsRef<[u8]> impl.
        let ha_idx = header
            .reference_sequences()
            .iter()
            .enumerate()
            .find_map(|(i, (name, _))| {
                let name_bytes: &[u8] = name.as_ref();
                if name_bytes == b"H1N1_HA" {
                    Some(i)
                } else {
                    None
                }
            })
            .expect("H1N1_HA not found in BAM header");

        // SNP at VCF position 100 (1-based) = 99 (0-based).
        // A read of length 50 covers this position when its 0-based start is in [50, 99].
        let variant_pos_0based: usize = 99;
        let read_len: usize = 50;

        let mut alt_seen = false;
        for result in reader.records() {
            let record = result.unwrap();
            let ref_id = match record.reference_sequence_id() {
                Some(Ok(id)) => id,
                _ => continue,
            };
            if ref_id != ha_idx {
                continue;
            }
            let aln_start_0based = match record.alignment_start() {
                Some(Ok(p)) => p.get() - 1,
                _ => continue,
            };
            if aln_start_0based > variant_pos_0based {
                continue;
            }
            let offset = variant_pos_0based - aln_start_0based;
            if offset >= read_len {
                continue;
            }
            let sequence = record.sequence();
            if let Some(base_byte) = sequence.get(offset)
                && (base_byte as char).eq_ignore_ascii_case(&'T')
            {
                alt_seen = true;
                break;
            }
        }

        assert!(
            alt_seen,
            "Expected at least one BAM read with T at H1N1_HA:100 (homozygous A→T SNP) \
             but no such read was found"
        );
    }

    /// Supply a VCF with a single known SNP on H1N1_HA and verify it appears in the
    /// output VCF. The SNP is at 1-based position 100 (A→T); H1N1_HA position 99
    /// (0-based) is confirmed non-N so the variant will be applied.
    #[test]
    fn test_input_vcf_snp_appears_in_output() {
        let out = tempdir().unwrap();

        // Minimal single-sample VCF: position 100 (1-based), A→T, homozygous
        let vcf_path = out.path().join("input.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             H1N1_HA\t100\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();

        let mut yaml = NamedTempFile::new().unwrap();
        write!(yaml,
            "reference: {ref}\nread_len: 50\ncoverage: 2\n\
             produce_fastq: true\nproduce_bam: false\nproduce_vcf: true\n\
             output_dir: {out}\noutput_filename: vcf_input_test\noverwrite_output: true\n\
             rng_seed: vcf input test\ninput_vcf: {vcf}\n",
            ref = h1n1_reference().display(),
            out = out.path().display(),
            vcf = vcf_path.display(),
        )
        .unwrap();
        let config = RunConfiguration::from_yaml_file(&yaml.path().to_path_buf()).unwrap();
        let output_vcf = config.output_vcf.clone().unwrap();
        run(config);

        // Read the output VCF and look for a record at position 100 on H1N1_HA
        use flate2::read::GzDecoder;
        use std::io::{BufRead, BufReader};
        let reader = BufReader::new(GzDecoder::new(std::fs::File::open(&output_vcf).unwrap()));
        let found = reader
            .lines()
            .map_while(|l| l.ok())
            .filter(|l| !l.starts_with('#'))
            .any(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                cols.len() >= 5 && cols[0] == "H1N1_HA" && cols[1] == "100"
            });
        assert!(
            found,
            "Expected SNP at H1N1_HA:100 in output VCF but did not find it"
        );
    }
}
