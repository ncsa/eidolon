// Unified runner: one BAM walk feeding multiple observers. Each requested
// model is built from its observer's accumulated state via the per-tool
// `run_from_*` helpers exposed by the standalone modules.
use log::info;
use std::path::PathBuf;

use eidolon_core::file_tools::bam_reader::{
    BamWalkFilter, CoverageObserver, FragLengthObserver, RecordObserver, walk_bam,
};

use crate::gen_bam_models::{
    errors::GenBamModelsError,
    utils::config::{GcBiasSection, RunConfiguration},
};
use crate::gen_frag_length_model::utils::runner::run_from_tlens;
use crate::gen_gc_bias_model::utils::config::GcBiasModelParams;
use crate::gen_gc_bias_model::utils::runner::run_from_coverage as gc_bias_run_from_coverage;

pub fn runner(config_path: &PathBuf) -> Result<(), GenBamModelsError> {
    let config = RunConfiguration::from(config_path)?;

    let mut frag_obs = FragLengthObserver::default();
    let mut cov_obs = CoverageObserver::default();

    // Build the observer list based on which sections are requested. Each
    // observer self-filters internally, so the walker uses the loosest
    // (coverage-style) filter modulated only by the top-level min_mapq.
    let mut observers: Vec<&mut dyn RecordObserver> = Vec::new();
    if config.frag_length.is_some() {
        observers.push(&mut frag_obs);
    }
    if config.gc_bias.is_some() {
        observers.push(&mut cov_obs);
    }

    let mut filter = BamWalkFilter::for_coverage();
    filter.min_mapq = config.min_mapq;

    info!(
        "Walking BAM {:?} once for {} observer(s)",
        config.bam_file,
        observers.len()
    );
    let stats = walk_bam(&config.bam_file, &filter, observers.as_mut_slice())?;
    info!(
        "BAM walk complete: {} records seen, {} kept",
        stats.records_seen, stats.records_kept
    );

    if let Some(fl) = &config.frag_length {
        info!("Building fragment-length model");
        run_from_tlens(
            frag_obs.tlens,
            fl.min_reads,
            &fl.output_file,
            fl.distribution,
        )?;
    }

    if let Some(gc) = &config.gc_bias {
        info!("Building GC bias model");
        let gc_params = build_gc_bias_model_params(gc);
        gc_bias_run_from_coverage(&gc_params, cov_obs.into_by_contig())?;
    }

    Ok(())
}

/// Adapts a unified `GcBiasSection` into a `GcBiasModelParams` for the
/// shared sweep + model-write stage. No walker fields are involved because
/// the unified runner has already done the BAM walk by the time we get here.
fn build_gc_bias_model_params(gc: &GcBiasSection) -> GcBiasModelParams {
    GcBiasModelParams {
        reference: gc.reference.clone(),
        bed_table: gc.bed_table.clone(),
        output_file: gc.output_file.clone(),
        overwrite_output: gc.overwrite_output,
        window_size: gc.window_size,
        window_stride: gc.window_stride,
        min_windows_per_bin: gc.min_windows_per_bin,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Writes a paired BAM with `n_pairs` proper pairs (first-in-pair at pos 1,
    /// mate at pos `tlen - read_len + 1`), plus `n_coverage_reads` unpaired
    /// reads stacked at pos 1 to inflate per-base depth. The same BAM thus
    /// exercises both FragLengthObserver (paired reads pass its filter) and
    /// CoverageObserver (all primary mapped reads count toward depth).
    fn write_test_bam(
        path: &std::path::PathBuf,
        contig: &[u8],
        contig_len: usize,
        n_pairs: usize,
        tlen: usize,
        read_len: usize,
        n_extra_coverage_reads: usize,
    ) {
        use noodles::bam;
        use noodles::core::Position;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };

        let header = sam::Header::builder()
            .add_reference_sequence(
                contig.to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(contig_len).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();

        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
        let seq = vec![b'A'; read_len];

        // First, n_pairs proper pairs.
        let paired_first = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let paired_second = Flags::SEGMENTED;
        for _ in 0..n_pairs {
            let mate_start = tlen - read_len + 1;
            // First in pair at position 1
            let mut r1 = RecordBuf::default();
            *r1.flags_mut() = paired_first;
            *r1.cigar_mut() = cigar.clone();
            *r1.reference_sequence_id_mut() = Some(0);
            *r1.alignment_start_mut() = Position::new(1);
            *r1.mate_reference_sequence_id_mut() = Some(0);
            *r1.mate_alignment_start_mut() = Position::new(mate_start);
            *r1.template_length_mut() = tlen as i32;
            *r1.sequence_mut() = Sequence::from(seq.as_slice());
            *r1.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            writer.write_alignment_record(&header, &r1).unwrap();

            // Second in pair at the mate position
            let mut r2 = RecordBuf::default();
            *r2.flags_mut() = paired_second;
            *r2.cigar_mut() = cigar.clone();
            *r2.reference_sequence_id_mut() = Some(0);
            *r2.alignment_start_mut() = Position::new(mate_start);
            *r2.mate_reference_sequence_id_mut() = Some(0);
            *r2.mate_alignment_start_mut() = Position::new(1);
            *r2.template_length_mut() = -(tlen as i32);
            *r2.sequence_mut() = Sequence::from(seq.as_slice());
            *r2.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            writer.write_alignment_record(&header, &r2).unwrap();
        }

        // Then n_extra_coverage_reads unpaired reads at pos 1 to inflate depth.
        for _ in 0..n_extra_coverage_reads {
            let mut r = RecordBuf::default();
            *r.flags_mut() = Flags::empty();
            *r.cigar_mut() = cigar.clone();
            *r.reference_sequence_id_mut() = Some(0);
            *r.alignment_start_mut() = Position::new(1);
            *r.sequence_mut() = Sequence::from(seq.as_slice());
            *r.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            writer.write_alignment_record(&header, &r).unwrap();
        }
    }

    fn write_temp(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    #[test]
    fn test_unified_runner_produces_both_models() {
        // 300bp ACGT reference, paired reads with TLEN=200, plus extra unpaired
        // coverage reads to push average depth higher. One BAM walk should
        // produce both a fragment-length model and a GC bias model.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("paired.bam");
        // Read length 100 means n=2 unique TLENs after dedup wouldn't trigger
        // min_reads=2 cleanly, so we lay 10 identical pairs to keep min_reads=2 happy.
        write_test_bam(&bam_path, b"chr1", 300, 10, 200, 100, 0);

        let seq: String = "ACGT".repeat(75); // 300bp, 50% GC
        let ref_file = write_temp(&format!(">chr1\n{}\n", seq));

        let frag_out = temp.path().join("frag.json.gz");
        let gc_out = temp.path().join("gc.json.gz");

        let yaml = write_temp(&format!(
            "bam_file: {}\nmin_mapq: 0\n\
             frag_length:\n  output_file: {}\n  overwrite_output: true\n  min_reads: 2\n  \
             distribution: normal\n\
             gc_bias:\n  reference: {}\n  output_file: {}\n  overwrite_output: true\n  \
             window_size: 100\n  window_stride: 100\n  min_windows_per_bin: 1\n",
            bam_path.display(),
            frag_out.display(),
            ref_file.path().display(),
            gc_out.display(),
        ));

        runner(&yaml.path().to_path_buf()).unwrap();
        assert!(frag_out.exists(), "fragment length model should be written");
        assert!(gc_out.exists(), "GC bias model should be written");

        // Load the GC bias model and confirm uniform coverage produces neutral weights.
        let gc_model =
            eidolon_core::models::gc_bias_model::GcBiasModel::from_file(&gc_out).unwrap();
        for gc_pct in 0..=100 {
            let w = gc_model.weight_for_gc_fraction(gc_pct as f64 / 100.0);
            assert!(
                (w - 1.0).abs() < 1e-6,
                "Expected weight ~1.0 at GC {}%, got {}",
                gc_pct,
                w
            );
        }
    }

    #[test]
    fn test_unified_runner_parity_with_per_tool_commands() {
        // The unified path must produce a fragment-length model identical to
        // what gen-frag-length-model produces on its own, since both use the
        // same FragLengthObserver against the same BAM.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("paired.bam");
        write_test_bam(&bam_path, b"chr1", 500, 10, 200, 100, 0);

        // Run unified path producing only the frag-length model.
        let unified_out = temp.path().join("frag_unified.json.gz");
        let unified_yaml = write_temp(&format!(
            "bam_file: {}\nfrag_length:\n  output_file: {}\n  overwrite_output: true\n  \
             min_reads: 2\n  distribution: normal\n",
            bam_path.display(),
            unified_out.display(),
        ));
        runner(&unified_yaml.path().to_path_buf()).unwrap();

        // Run standalone gen-frag-length-model on the same BAM.
        let standalone_out = temp.path().join("frag_standalone.json.gz");
        let standalone_yaml = write_temp(&format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 2\n\
             distribution: normal\n",
            bam_path.display(),
            standalone_out.display(),
        ));
        crate::gen_frag_length_model::main(&standalone_yaml.path().to_path_buf()).unwrap();

        // The two models should be numerically identical (same TLENs collected,
        // same filter, same normal-fit parameters).
        let m_unified =
            eidolon_core::models::fragment_length::FragmentLengthModel::discrete_from_file(
                &unified_out,
            )
            .unwrap();
        let m_standalone =
            eidolon_core::models::fragment_length::FragmentLengthModel::discrete_from_file(
                &standalone_out,
            )
            .unwrap();
        // Models serialize to JSON; the easiest equivalence check is to compare
        // the serialized representations.
        let s_unified = serde_json::to_string(&m_unified).unwrap();
        let s_standalone = serde_json::to_string(&m_standalone).unwrap();
        assert_eq!(
            s_unified, s_standalone,
            "unified and standalone frag-length models diverged"
        );
    }
}
