use crate::gen_seq_error_model::{errors::GenSeqErrorModelError, utils::config::RunConfiguration};
use eidolon_core::file_tools::file_io::is_gzipped_file;
use eidolon_core::{
    file_tools::{
        bam_reader::read_bam_transitions,
        file_io::{read_gzip_lines, read_lines},
    },
    models::{quality_scores::QualityScoreModel, sequencing_error_model::SequencingErrorModel},
    structs::transition_matrix::TransitionMatrix,
};
use log::{info, warn};
#[cfg(test)]
use std::path::PathBuf;

const MAX_SCORE: usize = 94;

/// Snap a raw quality score to the nearest value in a sorted bin list.
/// Ties round toward the lower bin (deterministic).
/// `bins` must be non-empty and sorted ascending.
fn snap_to_bin(score: usize, bins: &[usize]) -> usize {
    debug_assert!(!bins.is_empty(), "snap_to_bin called with empty bin list");
    // Find the first bin >= score. The nearest bin is either that one or the previous one.
    match bins.iter().position(|&b| b >= score) {
        None => *bins.last().unwrap(),
        Some(0) => bins[0],
        Some(i) => {
            let hi = bins[i];
            let lo = bins[i - 1];
            // Tie → lower bin.
            if hi - score < score - lo { hi } else { lo }
        }
    }
}

/// Normalizes a raw 4×4 mismatch count matrix into a `TransitionMatrix`.
///
/// Each row is normalized independently. Rows with no observed mismatches get
/// equal probability distributed across the three off-diagonal positions.
fn build_transition_matrix_from_counts(
    counts: [[usize; 4]; 4],
) -> Result<TransitionMatrix, GenSeqErrorModelError> {
    let mut weights = [[0.0f64; 4]; 4];
    for i in 0..4 {
        let total: f64 = counts[i].iter().sum::<usize>() as f64;
        if total == 0.0 {
            for j in 0..4 {
                if i != j {
                    weights[i][j] = 1.0 / 3.0;
                }
            }
        } else {
            for j in 0..4 {
                weights[i][j] = counts[i][j] as f64 / total;
            }
            weights[i][i] = 0.0;
        }
    }
    Ok(TransitionMatrix::from(
        weights[0], weights[1], weights[2], weights[3],
    )?)
}

pub fn runner(config: &RunConfiguration) -> Result<(), GenSeqErrorModelError> {
    let mut iter: Box<dyn Iterator<Item = std::io::Result<String>>> =
        if is_gzipped_file(&config.fastq_file)? {
            Box::new(read_gzip_lines(&config.fastq_file)?)
        } else {
            Box::new(read_lines(&config.fastq_file)?)
        };

    // Read the first record to determine read_length before allocating transition_counts.
    for _ in 0..3 {
        match iter.next() {
            Some(Ok(_)) => {}
            Some(Err(e)) => return Err(e.into()),
            None => {
                return Err(GenSeqErrorModelError::MalformedFastq(
                    "FASTQ file has fewer than 4 lines".to_string(),
                ));
            }
        }
    }
    let first_qual = match iter.next() {
        Some(Ok(q)) => q,
        Some(Err(e)) => return Err(e.into()),
        None => {
            return Err(GenSeqErrorModelError::MalformedFastq(
                "FASTQ file has fewer than 4 lines".to_string(),
            ));
        }
    };

    let read_length = first_qual.len();
    if read_length == 0 {
        return Err(GenSeqErrorModelError::MalformedFastq(
            "Quality line in first record is empty".to_string(),
        ));
    }

    let qual_offset = config.qual_offset;
    let mut seed_counts = vec![0usize; MAX_SCORE];
    let mut transition_counts = vec![vec![vec![0usize; MAX_SCORE]; MAX_SCORE]; read_length - 1];
    let mut global_counts = vec![0usize; MAX_SCORE];
    let mut total_bases: usize = 0;
    let mut reads_processed: usize = 0;

    fn accumulate_qual(
        qual_line: &str,
        qual_offset: usize,
        read_length: usize,
        bins: Option<&[usize]>,
        seed_counts: &mut [usize],
        transition_counts: &mut [Vec<Vec<usize>>],
        global_counts: &mut [usize],
        total_bases: &mut usize,
    ) -> bool {
        let scores: Vec<usize> = qual_line
            .bytes()
            .map(|b| {
                let raw = (b as usize).saturating_sub(qual_offset).min(MAX_SCORE - 1);
                match bins {
                    Some(bins) => snap_to_bin(raw, bins),
                    None => raw,
                }
            })
            .collect();
        if scores.is_empty() {
            return false;
        }
        seed_counts[scores[0]] += 1;
        for j in 1..scores.len().min(read_length) {
            transition_counts[j - 1][scores[j - 1]][scores[j]] += 1;
        }
        for &score in &scores {
            global_counts[score] += 1;
        }
        *total_bases += scores.len();
        true
    }

    let bins_slice: Option<&[usize]> = config.binned_quality_bins.as_deref();

    if accumulate_qual(
        &first_qual,
        qual_offset,
        read_length,
        bins_slice,
        &mut seed_counts,
        &mut transition_counts,
        &mut global_counts,
        &mut total_bases,
    ) {
        reads_processed += 1;
    }

    'records: loop {
        if config.max_reads > 0 && reads_processed >= config.max_reads {
            break;
        }
        for _ in 0..3 {
            match iter.next() {
                None => break 'records,
                Some(Ok(_)) => {}
                Some(Err(e)) => return Err(e.into()),
            }
        }
        let qual = match iter.next() {
            None => break,
            Some(Ok(q)) => q,
            Some(Err(e)) => return Err(e.into()),
        };
        if accumulate_qual(
            &qual,
            qual_offset,
            read_length,
            bins_slice,
            &mut seed_counts,
            &mut transition_counts,
            &mut global_counts,
            &mut total_bases,
        ) {
            reads_processed += 1;
        }
    }

    info!(
        "Processed {} reads ({} bases)",
        reads_processed, total_bases
    );

    if total_bases == 0 {
        return Err(GenSeqErrorModelError::MalformedFastq(
            "No bases found in FASTQ file".to_string(),
        ));
    }

    // Compute error_rate
    let error_rate: f64 = (0..MAX_SCORE)
        .map(|q| 10f64.powf(-(q as f64) / 10.0) * global_counts[q] as f64)
        .sum::<f64>()
        / total_bases as f64;

    info!("Computed error_rate: {:.6}", error_rate);

    // Collect quality_score_options. When binning is enabled, use the configured bin list
    // verbatim — even bins that received zero observed counts stay in the option set so
    // that from_counts' uniform fallback can produce a sane transition row for them.
    let (quality_score_options, is_binned): (Vec<usize>, bool) = match &config.binned_quality_bins {
        Some(bins) => {
            let zero_count_bins: Vec<usize> = bins
                .iter()
                .copied()
                .filter(|&b| global_counts[b] == 0)
                .collect();
            if !zero_count_bins.is_empty() {
                warn!(
                    "Binned quality model has {} bin(s) with no observed counts: {:?}; \
                     transition rows for these will use a uniform fallback",
                    zero_count_bins.len(),
                    zero_count_bins,
                );
            }
            (bins.clone(), true)
        }
        None => {
            let opts: Vec<usize> = (0..MAX_SCORE).filter(|&q| global_counts[q] > 0).collect();
            (opts, false)
        }
    };

    let n_scores = quality_score_options.len();

    // Re-index seed counts
    let seed_weights: Vec<f64> = quality_score_options
        .iter()
        .map(|&q| seed_counts[q] as f64)
        .collect();

    // Re-index transition counts
    let trans_weights: Vec<Vec<Vec<f64>>> = (0..read_length - 1)
        .map(|pos| {
            (0..n_scores)
                .map(|p| {
                    let prev_raw = quality_score_options[p];
                    (0..n_scores)
                        .map(|c| {
                            let curr_raw = quality_score_options[c];
                            transition_counts[pos][prev_raw][curr_raw] as f64
                        })
                        .collect()
                })
                .collect()
        })
        .collect();

    let quality_score_model = QualityScoreModel::from_counts(
        quality_score_options,
        read_length,
        seed_weights,
        trans_weights,
        is_binned,
    )?;

    // Determine transition matrix: TSV > BAM inference > defaults from original Python NEAT
    // (see https://github.com/ncsa/NEAT/blob/main/neat/model_sequencing_error/runner.py)
    let snp_transition_matrix: Option<TransitionMatrix> =
        if let Some(path) = &config.transition_matrix_file {
            info!("Loading SNP transition matrix from TSV: {:?}", path);
            Some(TransitionMatrix::from_tsv(path)?)
        } else if let Some(path) = &config.bam_file {
            info!("Inferring SNP transition matrix from BAM: {:?}", path);
            let counts = read_bam_transitions(path)?;
            let total_mismatches: usize = counts.iter().flatten().sum();
            if total_mismatches == 0 {
                warn!(
                    "BAM file {:?} had no MD-tagged records; using default transition matrix",
                    path
                );
                None
            } else {
                info!(
                    "Observed {} SNP mismatches across all records",
                    total_mismatches
                );
                Some(build_transition_matrix_from_counts(counts)?)
            }
        } else {
            None
        };

    let model = SequencingErrorModel::from_raw_data(
        error_rate,
        quality_score_model,
        snp_transition_matrix,
    )?;
    model.write_model(&config.output_file)?;

    info!("Wrote sequencing error model to {:?}", config.output_file);
    Ok(())
}

/// Write a minimal BAM file with `n_records` mapped records.
/// Each record uses `ref_bases` for the MD mismatch source and `read_bases` for the SEQ field.
/// Slices must be the same length. Pass `with_md = false` to omit the MD tag entirely.
#[cfg(test)]
fn write_test_bam(
    path: &PathBuf,
    n_records: usize,
    ref_bases: &[u8],
    read_bases: &[u8],
    with_md: bool,
) {
    use noodles::bam;
    use noodles::sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags,
                cigar::{Op, op::Kind},
                data::field::Tag,
            },
            record_buf::{Cigar, Sequence, data::field::Value as BufValue},
        },
    };

    let n = read_bases.len();
    let header = sam::Header::default();

    // Build MD string: match counts interspersed with ref bases at mismatches
    let md_str: Option<String> = if with_md {
        let mut md = String::new();
        let mut match_count = 0usize;
        for (&r, &q) in ref_bases.iter().zip(read_bases.iter()) {
            if r.eq_ignore_ascii_case(&q) {
                match_count += 1;
            } else {
                md.push_str(&match_count.to_string());
                md.push(r as char);
                match_count = 0;
            }
        }
        md.push_str(&match_count.to_string());
        Some(md)
    } else {
        None
    };

    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    for _ in 0..n_records {
        let cigar: Cigar = [Op::new(Kind::Match, n)].into_iter().collect();
        let mut record = RecordBuf::default();
        *record.flags_mut() = Flags::empty();
        *record.cigar_mut() = cigar;
        *record.sequence_mut() = Sequence::from(read_bases);

        if let Some(ref md) = md_str {
            record
                .data_mut()
                .insert(Tag::MISMATCHED_POSITIONS, BufValue::from(md.as_str()));
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
}

#[cfg(test)]
fn make_test_fastq(path: &PathBuf, n_reads: usize, read_length: usize) {
    use std::io::Write;
    let seq: String = "ACGT".chars().cycle().take(read_length).collect();
    // Quality scores cycling through a range: I(40), J(41), K(42), etc.
    let qual: String = (0..read_length)
        .map(|i| char::from_u32(('!' as u32) + 33 + (i % 10) as u32).unwrap())
        .collect();
    let mut f = std::fs::File::create(path).unwrap();
    for i in 0..n_reads {
        writeln!(f, "@read{}\n{}\n+\n{}", i, seq, qual).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gen_seq_error_model::utils::config::RunConfiguration;

    fn make_config(fastq: PathBuf, output: PathBuf) -> RunConfiguration {
        RunConfiguration {
            fastq_file: fastq,
            output_file: output,
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            binned_quality_bins: None,
            bam_file: None,
            transition_matrix_file: None,
        }
    }

    #[test]
    fn test_runner_basic() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 50);
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        assert!(output_path.exists());
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let scores = model
            .generate_quality_scores(
                50,
                &mut eidolon_core::rng::NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap(),
            )
            .unwrap();
        assert_eq!(scores.len(), 50);
    }

    #[test]
    fn test_runner_max_reads() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 100, 50);
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 10;
        runner(&config).unwrap();
        assert!(output_path.exists());
    }

    #[test]
    fn test_runner_truncated_fastq_errors() {
        // FASTQ with only 3 lines (no quality line) must yield MalformedFastq, not a panic
        // or an Ok with garbage data.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("trunc.fastq");
        std::fs::write(&fastq_path, "@h\nACGT\n+\n").unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path);
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenSeqErrorModelError::MalformedFastq(_)),
            "expected MalformedFastq for truncated FASTQ, got {err:?}",
        );
    }

    #[test]
    fn test_runner_empty_first_qual_line_errors() {
        // A FASTQ where the first record's quality line is empty (length 0) must surface as
        // MalformedFastq — the code computes read_length from that line.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("empty_qual.fastq");
        std::fs::write(&fastq_path, "@h\nA\n+\n\n").unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path);
        let err = runner(&config).unwrap_err();
        assert!(
            matches!(err, GenSeqErrorModelError::MalformedFastq(_)),
            "expected MalformedFastq for empty-qual FASTQ, got {err:?}",
        );
    }

    #[test]
    fn test_runner_empty_fastq_errors() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("empty.fastq");
        std::fs::write(&fastq_path, "").unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path);
        let result = runner(&config);
        assert!(result.is_err());
    }

    #[test]
    fn test_runner_error_rate_accuracy() {
        // All bases at Q30: '?' = ASCII 63, score = 63 - 33 = 30
        // Expected error_rate = 10^(-30/10) = 0.001 exactly
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("q30.fastq");
        let seq = "A".repeat(40);
        let qual = "?".repeat(40);
        let content: String = (0..20)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-30.0 / 10.0);
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected error_rate {expected:.6}, got {:.6}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_error_rate_survives_serialization() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("q20.fastq");
        let qual = "5".repeat(30);
        let seq = "C".repeat(30);
        let content: String = (0..10)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-20.0 / 10.0);
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected {expected:.6}, got {:.6}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_zero_transition_row_handled() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("sparse.fastq");
        let content: String = (0..20).map(|i| format!("@r{i}\nAAA\n+\n?=J\n")).collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let mut rng = eidolon_core::rng::NeatRng::new_from_seed(&vec!["s".to_string()]).unwrap();
        let scores: Vec<usize> = (0..100)
            .map(|_| model.generate_quality_scores(3, &mut rng).unwrap()[2])
            .collect();
        assert!(
            scores.contains(&41),
            "expected Q41 to appear in last-position scores with uniform fallback; got {scores:?}"
        );
    }

    #[test]
    fn test_runner_qual_offset() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("phred64.fastq");
        let qual = "a".repeat(20);
        let seq = "G".repeat(20);
        let content: String = (0..15)
            .map(|i| format!("@r{i}\n{seq}\n+\n{qual}\n"))
            .collect();
        std::fs::write(&fastq_path, content).unwrap();
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.qual_offset = 64;
        runner(&config).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let expected = 10f64.powf(-33.0 / 10.0);
        assert!(
            (model.error_rate() - expected).abs() < 1e-10,
            "expected {expected:.6e}, got {:.6e}",
            model.error_rate()
        );
    }

    #[test]
    fn test_runner_gzip_h1n1_r1() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_r1_model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let error_rate = model.error_rate();
        assert!(
            error_rate > 0.0001 && error_rate < 0.05,
            "error_rate {error_rate} outside expected Illumina range [0.0001, 0.05]"
        );
        let mut rng = eidolon_core::rng::NeatRng::new_from_seed(&vec!["h1n1".to_string()]).unwrap();
        let scores = model.generate_quality_scores(151, &mut rng).unwrap();
        assert_eq!(scores.len(), 151);
        assert!(
            scores.iter().all(|&s| s <= 94),
            "all scores should be ≤ MAX_SCORE"
        );
    }

    #[test]
    fn test_runner_gzip_h1n1_r2() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read2.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_r2_model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let error_rate = model.error_rate();
        assert!(
            error_rate > 0.0001 && error_rate < 0.05,
            "error_rate {error_rate} outside expected Illumina range [0.0001, 0.05]"
        );
    }

    #[test]
    fn test_runner_max_reads_reduces_processed_count() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fastq_path = PathBuf::from(format!("{}/test_data/H1N1_read1.fq.gz", manifest_dir));
        let temp = tempfile::tempdir().unwrap();
        let output_path = temp.path().join("h1n1_maxreads_model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 50;
        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let mut rng = eidolon_core::rng::NeatRng::new_from_seed(&vec!["seed".to_string()]).unwrap();
        let scores = model.generate_quality_scores(151, &mut rng).unwrap();
        assert_eq!(scores.len(), 151);
    }

    #[test]
    fn test_transition_matrix_from_tsv() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 50);
        let output_path = temp.path().join("model.json.gz");

        // Write a valid 4×4 TSV with a header line
        let tsv_path = temp.path().join("matrix.tsv");
        std::fs::write(
            &tsv_path,
            "from\\to\tA\tC\tG\tT\n\
             0.0\t0.5\t0.3\t0.2\n\
             0.5\t0.0\t0.3\t0.2\n\
             0.4\t0.3\t0.0\t0.3\n\
             0.3\t0.3\t0.4\t0.0\n",
        )
        .unwrap();

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            binned_quality_bins: None,
            bam_file: None,
            transition_matrix_file: Some(tsv_path),
        };
        runner(&config).unwrap();
        assert!(output_path.exists());
    }

    /// The cumulative weights of one transition-matrix row, as `[A, C, G, T]`.
    ///
    /// `DiscreteDistribution` stores a CDF rather than per-value probabilities, so a row
    /// whose whole weight sits on C reads `[0.0, 1.0, 1.0, 1.0]`, and one split
    /// 0.5/0/0.25/0.25 reads `[0.5, 0.5, 0.75, 1.0]`.
    fn row_cdf(
        tm: &eidolon_core::structs::transition_matrix::TransitionMatrix,
        base: eidolon_core::structs::nucleotides::Nucleotide,
    ) -> Vec<f64> {
        tm[&base].weights().unwrap()
    }

    fn assert_row_cdf_eq(actual: &[f64], expected: &[f64; 4], what: &str) {
        assert_eq!(actual.len(), 4, "{what}: expected a 4-wide row");
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (actual[i] - exp).abs() < 1e-9,
                "{what}: position {i} was {}, expected {exp} (full row {actual:?})",
                actual[i]
            );
        }
    }

    /// The default SEQUENCING-ERROR matrix's A row as a CDF, derived rather than transcribed.
    ///
    /// A hardcoded copy would be a vacuity hazard precisely here: the assertions below use
    /// this to prove a built matrix is *not* the default. If the default changed and a
    /// transcribed constant did not, the check would compare against a value that is no
    /// longer the default — and could pass while the built matrix IS the new default.
    ///
    /// DO NOT reach for `TransitionMatrix::default()`. Two differently-valued defaults share
    /// that name: `TransitionMatrix::default()` is NEAT 2.0's *mutation* matrix (A row
    /// 0.0/0.1695/0.6878/0.1427), while the sequencing-error default is a separate literal
    /// inside `SequencingErrorModel` (A row 0.0/0.4918/0.3377/0.1705). Deriving from the
    /// wrong one makes this test fail against correct code, which is how the distinction was
    /// found.
    ///
    /// Residual, pre-existing: that literal appears TWICE in `sequencing_error_model.rs`
    /// (`default()` and `from_raw_data`'s `None` arm). This reads the former while the
    /// runner's fallback path uses the latter, so a drift between the two copies would slip
    /// past. They are identical today.
    fn default_a_row_cdf() -> [f64; 4] {
        let m = SequencingErrorModel::default()
            .expect("the default sequencing error model must be constructible");
        let row = row_cdf(
            m.transition_distros(),
            eidolon_core::structs::nucleotides::Nucleotide::A,
        );
        [row[0], row[1], row[2], row[3]]
    }

    #[test]
    fn test_runner_with_bam_md_tags_puts_all_a_weight_on_c() {
        // ref=AAAA, read=CCCC → every mismatch is A→C, so the built model's A row must
        // place its entire weight on C. Asserting the artifact's content, not that a file
        // appeared: the previous version of this test passed with any matrix at all,
        // including the default, which is precisely the failure it needed to catch.
        use eidolon_core::structs::nucleotides::Nucleotide;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("test.bam");
        write_test_bam(&bam_path, 10, b"AAAA", b"CCCC", true);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
        };
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let tm = model.transition_distros();

        // A → C with probability 1: CDF steps to 1.0 at C and stays there.
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::A),
            &[0.0, 1.0, 1.0, 1.0],
            "A row inferred from an all-A→C BAM",
        );

        // And it must not simply be the default matrix wearing a disguise.
        let a_row = row_cdf(tm, Nucleotide::A);
        assert!(
            (a_row[1] - default_a_row_cdf()[1]).abs() > 1e-6,
            "A row matched the default matrix, so the BAM was not consulted: {a_row:?}"
        );

        // Rows the BAM said nothing about fall back to uniform off-diagonal (1/3 each):
        // C row CDF = [1/3, 1/3, 2/3, 1.0].
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::C),
            &[1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 1.0],
            "C row, unobserved in the BAM",
        );
    }

    #[test]
    fn test_runner_bam_transitions_track_a_mixed_mismatch_pattern() {
        // The fixture must be ASYMMETRIC, and each row must spread over more than one
        // target base. A reciprocal pattern (ref=ACGT / read=CATG, giving A→C and C→A in
        // equal number) produces a count matrix equal to its own transpose, so swapping
        // the ref/read axes cannot be detected — the first version of this test had that
        // flaw and a transpose mutation passed it. Single-entry rows are no good either:
        // the distribution is normalized, so one nonzero cell lands on probability 1.0
        // wherever it sits.
        //
        // ref=AAAACCCC vs read=CCGTAAAG gives
        //   A row: A→C ×2, A→G ×1, A→T ×1  → 0.50 / 0.25 / 0.25
        //   C row: C→A ×3, C→G ×1          → 0.75 / 0.25
        // counts[A][C]=2 against counts[C][A]=3, so the matrix is not symmetric.
        use eidolon_core::structs::nucleotides::Nucleotide;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("mixed.bam");
        write_test_bam(&bam_path, 8, b"AAAACCCC", b"CCGTAAAG", true);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
        };
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let tm = model.transition_distros();

        // A: 0.50 C, 0.25 G, 0.25 T → CDF [0.0, 0.50, 0.75, 1.0].
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::A),
            &[0.0, 0.5, 0.75, 1.0],
            "A row (2 C, 1 G, 1 T)",
        );
        // C: 0.75 A, 0.25 G → CDF [0.75, 0.75, 1.0, 1.0].
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::C),
            &[0.75, 0.75, 1.0, 1.0],
            "C row (3 A, 1 G)",
        );
        // G and T were never the reference base here, so both fall back to uniform.
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::G),
            &[1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0],
            "G row, unobserved",
        );
    }

    #[test]
    fn test_build_transition_matrix_from_counts_is_not_transposed() {
        // Guards the ref/read axis directly on the pure function, where the fixture can be
        // made maximally asymmetric without having to express it as a BAM. Transposing
        // `counts[i][j]` changes every row here.
        use eidolon_core::structs::nucleotides::Nucleotide;
        let mut counts = [[0usize; 4]; 4];
        // A row: 6 C, 3 G, 1 T  → 0.6 / 0.3 / 0.1
        counts[0][1] = 6;
        counts[0][2] = 3;
        counts[0][3] = 1;
        // C row: 1 A, 1 G, 8 T  → 0.1 / 0.1 / 0.8
        counts[1][0] = 1;
        counts[1][2] = 1;
        counts[1][3] = 8;
        let tm = build_transition_matrix_from_counts(counts).unwrap();

        assert_row_cdf_eq(
            &row_cdf(&tm, Nucleotide::A),
            &[0.0, 0.6, 0.9, 1.0],
            "A row 0.6/0.3/0.1",
        );
        assert_row_cdf_eq(
            &row_cdf(&tm, Nucleotide::C),
            &[0.1, 0.1, 0.2, 1.0],
            "C row 0.1/0.1/0.8",
        );
    }

    #[test]
    fn test_runner_bam_no_md_tags_falls_back() {
        // BAM records without MD tags → runner succeeds and uses the default matrix.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("nomd.bam");
        write_test_bam(&bam_path, 5, b"ACGT", b"TGCA", false);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
        };
        runner(&config).unwrap();

        // The case where inference must NOT fire: no MD tags means nothing was observed,
        // so the model has to carry the default matrix. Asserting only "a file appeared"
        // could not tell that apart from silently inventing a matrix from zero evidence.
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        assert_row_cdf_eq(
            &row_cdf(
                model.transition_distros(),
                eidolon_core::structs::nucleotides::Nucleotide::A,
            ),
            &default_a_row_cdf(),
            "A row with no MD tags available",
        );
    }

    #[test]
    fn test_tsv_takes_precedence_over_bam() {
        // When both transition_matrix_file and bam_file are set, the TSV wins. The BAM is
        // all A→C; the TSV's A row is 0.5/0.3/0.2 across C/G/T. Those are distinguishable,
        // so the assertion can actually establish precedence — the earlier version only
        // checked that the runner completed, and would have passed with precedence
        // inverted, which is the single thing it existed to rule out.
        use eidolon_core::structs::nucleotides::Nucleotide;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("test.bam");
        write_test_bam(&bam_path, 10, b"AAAA", b"CCCC", true);
        let tsv_path = temp.path().join("matrix.tsv");
        std::fs::write(
            &tsv_path,
            "A\tC\tG\tT\n\
             0.0\t0.5\t0.3\t0.2\n\
             0.5\t0.0\t0.3\t0.2\n\
             0.4\t0.3\t0.0\t0.3\n\
             0.3\t0.3\t0.4\t0.0\n",
        )
        .unwrap();
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: Some(tsv_path),
        };
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let a_row = row_cdf(model.transition_distros(), Nucleotide::A);

        // TSV A row 0.0/0.5/0.3/0.2 → CDF [0.0, 0.5, 0.8, 1.0].
        assert_row_cdf_eq(&a_row, &[0.0, 0.5, 0.8, 1.0], "A row taken from the TSV");

        // The BAM would have produced [0.0, 1.0, 1.0, 1.0]. State the negative directly.
        assert!(
            (a_row[2] - 1.0).abs() > 1e-6,
            "A row looks BAM-derived, so the TSV did not take precedence: {a_row:?}"
        );
    }

    #[test]
    fn test_build_transition_matrix_from_counts_uniform_fallback() {
        use eidolon_core::structs::nucleotides::Nucleotide;
        // A row with all zeros should produce uniform off-diagonal distribution —
        // verify by sampling and checking that A is never the result and all three
        // other bases are reachable.
        let mut counts = [[0usize; 4]; 4];
        counts[1][0] = 10; // C→A
        counts[1][2] = 5; // C→G
        counts[1][3] = 5; // C→T
        let tm = build_transition_matrix_from_counts(counts).unwrap();

        // The observed row is the point of the function and went unasserted: 10/5/5 out of
        // 20 is 0.5 / 0.25 / 0.25 across A / G / T, so the C row's CDF is
        // [0.5, 0.5, 0.75, 1.0]. Without this, the counts could have been ignored entirely
        // and only the untouched A row was ever checked.
        assert_row_cdf_eq(
            &row_cdf(&tm, Nucleotide::C),
            &[0.5, 0.5, 0.75, 1.0],
            "C row fitted from 10/5/5 counts",
        );

        // Sample the A row at several evenly-spaced random values
        let a_dist = &tm[&Nucleotide::A];
        let mut seen = std::collections::HashSet::new();
        for k in 1..=99 {
            let r = k as f64 / 100.0;
            let result = a_dist.sample(r).unwrap();
            assert_ne!(
                result,
                Nucleotide::A,
                "A→A self-transition should be impossible"
            );
            seen.insert(result as usize);
        }
        assert_eq!(
            seen.len(),
            3,
            "all three off-diagonal bases should be reachable"
        );
    }

    #[test]
    fn test_snap_to_bin_basic() {
        let bins = [2usize, 12, 23, 37];
        // Below smallest bin
        assert_eq!(snap_to_bin(0, &bins), 2);
        assert_eq!(snap_to_bin(2, &bins), 2);
        // Above largest bin
        assert_eq!(snap_to_bin(40, &bins), 37);
        assert_eq!(snap_to_bin(93, &bins), 37);
        // Interior — nearest bin
        assert_eq!(snap_to_bin(11, &bins), 12);
        assert_eq!(snap_to_bin(13, &bins), 12);
        assert_eq!(snap_to_bin(20, &bins), 23);
        // Tie: midpoint between 2 and 12 is 7 — rounds to lower (2).
        assert_eq!(snap_to_bin(7, &bins), 2);
        // Midpoint between 23 and 37 is 30 — rounds to lower (23).
        assert_eq!(snap_to_bin(30, &bins), 23);
    }

    #[test]
    fn test_snap_to_bin_singleton() {
        let bins = [12usize];
        assert_eq!(snap_to_bin(0, &bins), 12);
        assert_eq!(snap_to_bin(12, &bins), 12);
        assert_eq!(snap_to_bin(50, &bins), 12);
    }

    #[test]
    fn test_runner_binned_quality_scores() {
        // Synthetic FASTQ has scores in range 33..=42 (see make_test_fastq). With bins
        // [2, 12, 23, 37], every observed score should snap to 37, so the only quality
        // option used is 37 — but the model must keep all four bins in
        // quality_score_options and flag binned_scores = true.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 30, 50);
        let output_path = temp.path().join("model.json.gz");

        let mut config = make_config(fastq_path, output_path.clone());
        config.binned_quality_bins = Some(vec![2, 12, 23, 37]);

        runner(&config).unwrap();
        assert!(output_path.exists());

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let qsm = model.quality_score_model();
        assert!(qsm.binned_scores, "model should be marked binned");
        assert_eq!(qsm.quality_score_options, vec![2, 12, 23, 37]);
    }

    #[test]
    fn test_runner_binned_emits_only_bin_values() {
        // End-to-end: train a binned model, sample many quality vectors, and verify
        // every emitted score is in the bin set.
        use eidolon_core::rng::NeatRng;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 30, 50);
        let output_path = temp.path().join("model.json.gz");

        let mut config = make_config(fastq_path, output_path.clone());
        config.binned_quality_bins = Some(vec![2, 12, 23, 37]);
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let qsm = model.quality_score_model();
        let bins: std::collections::HashSet<usize> =
            qsm.quality_score_options.iter().copied().collect();

        let mut rng = NeatRng::new_from_seed(&vec!["binned-runner".to_string()]).unwrap();
        for _ in 0..200 {
            let scores = qsm.generate_quality_scores(50, &mut rng).unwrap();
            for &s in &scores {
                assert!(bins.contains(&s), "sampled non-bin score {s}");
            }
        }
    }

    #[test]
    fn test_runner_binned_error_rate_differs_from_unbinned() {
        // Same FASTQ, same everything except the binning. Quality scores in the synthetic
        // FASTQ span 33..=42 (see make_test_fastq); binning to [2, 12, 23, 37] snaps every
        // observed score to Q37, which has a smaller per-base error probability than the
        // average of 33..=42. The serialized error_rate must reflect that — if it doesn't,
        // we've forgotten to use the snapped counts somewhere in the pipeline.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 60, 80);

        // Run 1: unbinned.
        let out_unbinned = temp.path().join("unbinned.json.gz");
        let config = make_config(fastq_path.clone(), out_unbinned.clone());
        runner(&config).unwrap();
        let er_unbinned = SequencingErrorModel::from_file(&out_unbinned)
            .unwrap()
            .error_rate();

        // Run 2: binned.
        let out_binned = temp.path().join("binned.json.gz");
        let mut binned_config = make_config(fastq_path, out_binned.clone());
        binned_config.binned_quality_bins = Some(vec![2, 12, 23, 37]);
        runner(&binned_config).unwrap();
        let er_binned = SequencingErrorModel::from_file(&out_binned)
            .unwrap()
            .error_rate();

        assert!(er_unbinned > 0.0, "unbinned error_rate should be positive");
        assert!(er_binned > 0.0, "binned error_rate should be positive");
        // Q37 → 10^-3.7 ≈ 2.00e-4. Average of Q33..=Q42 → ~2.2e-4. They must differ.
        assert!(
            (er_unbinned - er_binned).abs() > 1e-6,
            "binning must change error_rate (unbinned={er_unbinned}, binned={er_binned})",
        );
        // Expect binned < unbinned in this specific synthetic case (all snaps go up to Q37,
        // which is on the lower-error end of the 33..=42 range).
        assert!(
            er_binned < er_unbinned,
            "binning to Q37 in this fixture should reduce error_rate \
             (unbinned={er_unbinned}, binned={er_binned})",
        );
    }
}
