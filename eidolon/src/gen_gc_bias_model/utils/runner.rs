use log::*;
use std::{collections::HashMap, path::PathBuf};

use crate::gen_gc_bias_model::{
    errors::GenGcBiasModelError,
    utils::config::{GcBiasModelParams, RunConfiguration},
};
use eidolon_core::{
    file_tools::{
        bam_reader::{BamWalkFilter, CoverageObserver, walk_bam},
        fasta_stream::{FastaStream, non_n_regions},
    },
    models::gc_bias_model::GcBiasModel,
    structs::nucleotides::Nucleotide,
};

pub fn runner(path: &PathBuf) -> Result<(), GenGcBiasModelError> {
    let config = RunConfiguration::from(path)?;

    info!("Accumulating coverage from BAM {:?}", config.bam_file);
    let mut obs = CoverageObserver::default();
    let mut filter = BamWalkFilter::for_coverage();
    filter.min_mapq = config.min_mapq;
    walk_bam(&config.bam_file, &filter, &mut [&mut obs])?;
    run_from_coverage(&config.model, obs.into_by_contig())
}

/// Builds and writes a GC bias model from a pre-computed per-contig coverage
/// map (e.g. produced by `CoverageObserver` during a shared BAM walk in the
/// unified `gen-bam-models` runner). Reads the reference FASTA, sweeps
/// windows, accumulates GC% vs mean-coverage, and writes the gzipped JSON
/// model to `config.output_file`.
pub fn run_from_coverage(
    config: &GcBiasModelParams,
    mut cov_by_contig: HashMap<String, Vec<u32>>,
) -> Result<(), GenGcBiasModelError> {
    // Every window's mean coverage, by GC percent.
    let mut bins: Vec<Vec<f32>> = vec![Vec::new(); 101];

    info!("Processing reference {:?}", config.reference);
    let fasta = FastaStream::open(&config.reference)?;

    for result in fasta {
        let (contig_name, raw) = result?;
        // IUPAC codes map to N here intentionally — GC bias model training works on
        // observed coverage and doesn't require stochastic base resolution.
        let sequence: Vec<Nucleotide> = raw.chars().map(Nucleotide::from).collect();
        let contig_len = sequence.len();

        if contig_len < config.window_size {
            debug!(
                "Skipping {} (length {} < window_size {})",
                contig_name, contig_len, config.window_size
            );
            continue;
        }

        let cov = match cov_by_contig.remove(&contig_name) {
            Some(d) => d,
            None => {
                debug!("No coverage data for {}, skipping", contig_name);
                continue;
            }
        };

        let regions: Vec<(usize, usize)> =
            if let Some(bed_regions) = config.bed_table.get(&contig_name) {
                bed_regions
                    .iter()
                    .map(|r| (r.start, r.end.min(contig_len)))
                    .collect()
            } else {
                non_n_regions(&sequence)
            };

        for (region_start, region_end) in regions {
            if region_end.saturating_sub(region_start) < config.window_size {
                continue;
            }
            accumulate_region(
                &sequence,
                &cov,
                region_start,
                region_end,
                config.window_size,
                config.window_stride,
                &mut bins,
            );
        }

        debug!("Processed {}", contig_name);
    }

    let total_windows: usize = bins.iter().map(Vec::len).sum();
    if total_windows == 0 {
        return Err(GenGcBiasModelError::ConfigError(
            "No windows were processed — verify that the BAM and reference share contig names"
                .to_string(),
        ));
    }

    let any_supported = bins
        .iter()
        .any(|v| !v.is_empty() && v.len() >= config.min_windows_per_bin);
    let fitted = fit_weights(&mut bins, config.min_windows_per_bin);
    // Supported bins but nothing fitted means the overall median is zero: most windows have
    // no coverage, so no weight is defined. Refuse rather than write a neutral model.
    if any_supported && fitted.iter().all(|b| !b.fitted) {
        return Err(GenGcBiasModelError::ConfigError(format!(
            "the median window across the reference has zero coverage ({total_windows} windows), \
             so GC weights cannot be measured. The BAM covers only part of the reference, as a \
             targeted or exome library does. Set bed_file to the regions it covers."
        )));
    }
    if fitted.iter().all(|b| !b.fitted) {
        warn!(
            "No GC bins met min_windows_per_bin ({}); all weights will be neutral (1.0). \
             Try lowering min_windows_per_bin or using a larger genome/region.",
            config.min_windows_per_bin
        );
    }
    let weights: Vec<f64> = fitted.iter().map(|b| b.weight).collect();
    let mut report = config.output_file.clone().into_os_string();
    report.push(".bins.tsv");
    let report = PathBuf::from(report);
    write_bin_report(&report, &fitted)?;
    info!(
        "{} of 101 GC bins fitted from at least {} windows; the rest interpolated. \
         Per-bin support: {:?}",
        fitted.iter().filter(|b| b.fitted).count(),
        config.min_windows_per_bin,
        report
    );

    let model = GcBiasModel::from_weights(weights, config.window_size)?;
    model.write_to_file(&config.output_file)?;
    info!("GC bias model written to {:?}", config.output_file);

    Ok(())
}

/// One GC bin's support and fitted weight.
#[derive(Debug, Clone, PartialEq)]
struct GcBin {
    windows: usize,
    mean: f64,
    median: f64,
    weight: f64,
    /// True when the bin had at least `min_windows_per_bin` windows and its weight comes from
    /// its own data; false when it was filled in.
    fitted: bool,
}

/// Median of an ascending slice; the middle two averaged for an even count.
fn median_sorted(v: &[f32]) -> f64 {
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2] as f64
    } else {
        (v[n / 2 - 1] as f64 + v[n / 2] as f64) / 2.0
    }
}

/// Weights for the 101 GC bins, from each bin's window coverages. Sorts each bin in place.
///
/// A bin with at least `min_windows` windows is fitted: its median coverage over the median of
/// every window in a fitted bin. Medians, because a few pileup windows (collapsed repeats) move
/// a mean arbitrarily far, and the simulator cannot reproduce them from a single-copy reference.
/// Every other bin is interpolated linearly between the nearest fitted bins on each side, or
/// takes the nearest fitted bin's weight past the ends, so a thinly supported bin neither jumps
/// to 1.0 nor rests on its own few windows. With no fitted bin, every weight is 1.0.
fn fit_weights(bins: &mut [Vec<f32>], min_windows: usize) -> Vec<GcBin> {
    for v in bins.iter_mut() {
        v.sort_by(f32::total_cmp);
    }
    let supported: Vec<bool> = bins
        .iter()
        .map(|v| !v.is_empty() && v.len() >= min_windows)
        .collect();
    let mut pooled: Vec<f32> = bins
        .iter()
        .zip(&supported)
        .filter(|(_, s)| **s)
        .flat_map(|(v, _)| v.iter().copied())
        .collect();
    pooled.sort_by(f32::total_cmp);
    let overall = if pooled.is_empty() {
        0.0
    } else {
        median_sorted(&pooled)
    };

    let mut out: Vec<GcBin> = bins
        .iter()
        .zip(&supported)
        .map(|(v, &s)| {
            let n = v.len();
            let median = if n > 0 { median_sorted(v) } else { 0.0 };
            let fitted = s && overall > 0.0;
            GcBin {
                windows: n,
                mean: if n > 0 {
                    v.iter().map(|&x| x as f64).sum::<f64>() / n as f64
                } else {
                    0.0
                },
                median,
                weight: if fitted { median / overall } else { 1.0 },
                fitted,
            }
        })
        .collect();

    let anchors: Vec<usize> = (0..out.len()).filter(|&i| out[i].fitted).collect();
    if let (Some(&first), Some(&last)) = (anchors.first(), anchors.last()) {
        for i in 0..out.len() {
            if out[i].fitted {
                continue;
            }
            out[i].weight = if i < first {
                out[first].weight
            } else if i > last {
                out[last].weight
            } else {
                let hi = *anchors.iter().find(|&&a| a > i).unwrap();
                let lo = *anchors.iter().rev().find(|&&a| a < i).unwrap();
                let t = (i - lo) as f64 / (hi - lo) as f64;
                out[lo].weight + (out[hi].weight - out[lo].weight) * t
            };
        }
    }
    out
}

/// Writes one row per GC percent: the windows behind each weight, and whether it was fitted
/// from them or interpolated.
fn write_bin_report(path: &PathBuf, bins: &[GcBin]) -> std::io::Result<()> {
    let mut body = String::from("gc_percent\twindows\tmedian\tmean\tweight\tsource\n");
    for (gc, b) in bins.iter().enumerate() {
        body.push_str(&format!(
            "{gc}\t{}\t{:.3}\t{:.3}\t{:.6}\t{}\n",
            b.windows,
            b.median,
            b.mean,
            b.weight,
            if b.fitted { "fitted" } else { "interpolated" }
        ));
    }
    std::fs::write(path, body)
}

/// Accumulate GC% vs mean-coverage data for all windows in `[region_start, region_end)`.
///
/// Uses a sliding window for both GC counting and coverage summing — O(region_len)
/// regardless of window_size or window_stride.
fn accumulate_region(
    sequence: &[Nucleotide],
    cov: &[u32],
    region_start: usize,
    region_end: usize,
    window_size: usize,
    window_stride: usize,
    bins: &mut [Vec<f32>],
) {
    let first = &sequence[region_start..region_start + window_size];
    let mut gc_count: usize = first
        .iter()
        .filter(|n| matches!(n.get_unmasked_base(), Nucleotide::G | Nucleotide::C))
        .count();
    let mut n_count: usize = first
        .iter()
        .filter(|n| n.get_unmasked_base() == Nucleotide::N)
        .count();

    let mut cov_sum: u64 = (0..window_size)
        .map(|i| cov.get(region_start + i).copied().unwrap_or(0) as u64)
        .sum();

    let mut w = region_start;
    loop {
        let called = window_size - n_count;
        if called > 0 {
            let gc_pct = ((gc_count as f64 / called as f64) * 100.0).round() as usize;
            let gc_pct = gc_pct.min(100);
            bins[gc_pct].push((cov_sum as f64 / window_size as f64) as f32);
        }

        let next_w = w + window_stride;
        if next_w + window_size > region_end {
            break;
        }

        // Advance both sliding windows by stride steps.
        for j in 0..window_stride {
            // Soft-masked (lowercase) bases count as the base they are (#771).
            match sequence[w + j].get_unmasked_base() {
                Nucleotide::G | Nucleotide::C => gc_count -= 1,
                Nucleotide::N => n_count -= 1,
                _ => {}
            }
            match sequence[w + window_size + j].get_unmasked_base() {
                Nucleotide::G | Nucleotide::C => gc_count += 1,
                Nucleotide::N => n_count += 1,
                _ => {}
            }
            cov_sum -= cov.get(w + j).copied().unwrap_or(0) as u64;
            cov_sum += cov.get(w + window_size + j).copied().unwrap_or(0) as u64;
        }

        w = next_w;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    fn write_bam_config(
        reference: &PathBuf,
        bam_file: &PathBuf,
        output_file: &PathBuf,
        window_size: usize,
        window_stride: usize,
        min_windows: usize,
    ) -> NamedTempFile {
        let yaml = format!(
            "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n\
             window_size: {}\nwindow_stride: {}\nmin_windows_per_bin: {}\n",
            reference.display(),
            bam_file.display(),
            output_file.display(),
            window_size,
            window_stride,
            min_windows,
        );
        write_temp(&yaml)
    }

    /// Writes a BGZF BAM at `path` that produces per-base reference coverage equal
    /// to the sum of `depth` over all overlapping `(start_0based, end_0based, depth)`
    /// segments. Each segment is realized by stacking `depth` reads of length
    /// `end - start` at position `start + 1`.
    fn write_coverage_bam(
        path: &std::path::PathBuf,
        contigs: &[(&[u8], usize)],
        // (ref_id, start_0based, end_0based, depth)
        segments: &[(usize, usize, usize, usize)],
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

        let mut builder = sam::Header::builder();
        for &(name, len) in contigs {
            builder = builder.add_reference_sequence(
                name.to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(len).unwrap()),
            );
        }
        let header = builder.build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();

        for &(ref_id, start_0, end_0, depth) in segments {
            let read_len = end_0 - start_0;
            if read_len == 0 || depth == 0 {
                continue;
            }
            let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
            let seq = vec![b'A'; read_len];
            for _ in 0..depth {
                let mut record = RecordBuf::default();
                *record.flags_mut() = Flags::empty();
                *record.cigar_mut() = cigar.clone();
                *record.reference_sequence_id_mut() = Some(ref_id);
                *record.alignment_start_mut() = Position::new(start_0 + 1);
                *record.sequence_mut() = Sequence::from(seq.as_slice());
                *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
                writer.write_alignment_record(&header, &record).unwrap();
            }
        }
    }

    /// `bins[gc]` = the given window coverages; every other bin empty.
    fn bins_with(entries: &[(usize, &[f32])]) -> Vec<Vec<f32>> {
        let mut bins = vec![Vec::new(); 101];
        for (gc, v) in entries {
            bins[*gc] = v.to_vec();
        }
        bins
    }

    // CORRECTNESS CRITERION, #752. A bin's weight is its typical window, not its average: a few
    // pileup windows (collapsed repeats) swung real bins either side of 1.0 on NA12878 chr22,
    // one reading 18,138x against a median of 298. Known answer by hand: bin 40 is
    // [10,10,10,10,1000] and bin 60 is [20]*5. The overall median of those ten windows is 20,
    // so bin 40 weighs 10/20 = 0.5 and bin 60 20/20 = 1.0. The mean would put bin 40 at ~10.4.
    #[test]
    fn a_bins_weight_is_its_median_over_the_overall_median() {
        let mut bins = bins_with(&[(40, &[10.0, 10.0, 10.0, 10.0, 1000.0]), (60, &[20.0; 5])]);
        let fit = fit_weights(&mut bins, 5);
        assert_eq!(fit[40].weight, 0.5, "bin 40: {:?}", fit[40]);
        assert_eq!(fit[60].weight, 1.0, "bin 60: {:?}", fit[60]);
        assert!(fit[40].fitted && fit[60].fitted);
        assert_eq!((fit[40].windows, fit[40].median), (5, 10.0));
    }

    // A bin below min_windows_per_bin used to snap to exactly 1.0, a cliff between its
    // neighbors (0.86 then 1.00 at 90/91% GC on chr22). It now takes the value interpolated
    // between the nearest supported bins, and holds the edge value past the last one. Its own
    // windows, however extreme, do not set its weight.
    #[test]
    fn a_sparse_bin_is_interpolated_between_its_supported_neighbors() {
        // Unequal support on purpose: with equal counts the two weights average to exactly
        // 1.0, which is also what the old cliff gave, so a midpoint check could not tell them
        // apart. Twelve supported windows [10]*5 + [20]*7: overall median 20, so bin 40 weighs
        // 0.5 and bin 60 1.0. Bin 45 sits a quarter of the way across: 0.625.
        let mut bins = bins_with(&[
            (40, &[10.0; 5]),
            (45, &[5000.0]),
            (60, &[20.0; 7]),
            (95, &[1.0]),
        ]);
        let fit = fit_weights(&mut bins, 5);
        assert_eq!((fit[40].weight, fit[60].weight), (0.5, 1.0));
        assert!(
            (fit[45].weight - 0.625).abs() < 1e-12,
            "bin 45: {:?}",
            fit[45]
        );
        assert!(!fit[45].fitted);
        assert!(
            (fit[50].weight - 0.75).abs() < 1e-12,
            "empty bin 50: {:?}",
            fit[50]
        );
        // Past the supported range: the edge value, on both sides, empty bins included.
        assert_eq!(fit[0].weight, 0.5, "bin 0: {:?}", fit[0]);
        assert_eq!(fit[95].weight, 1.0, "bin 95: {:?}", fit[95]);
        assert_eq!(fit[100].weight, 1.0, "bin 100: {:?}", fit[100]);
    }

    /// With no supported bin there is nothing to measure, so every weight stays neutral.
    #[test]
    fn with_no_supported_bin_every_weight_is_neutral() {
        let mut bins = bins_with(&[(50, &[7.0, 9.0])]);
        let fit = fit_weights(&mut bins, 5);
        assert!(fit.iter().all(|b| b.weight == 1.0 && !b.fitted));
    }

    fn params(reference: &std::path::Path, out: &std::path::Path) -> GcBiasModelParams {
        GcBiasModelParams {
            reference: reference.to_path_buf(),
            bed_table: HashMap::new(),
            output_file: out.to_path_buf(),
            overwrite_output: true,
            window_size: 100,
            window_stride: 100,
            min_windows_per_bin: 1,
        }
    }

    // Every fit writes a per-bin report next to the model, so a surprising weight can be traced
    // to the windows behind it: HG002's 94% bin read 34x and nothing said how many windows made
    // it. One 0%-GC window at depth 10 and one 100%-GC window at depth 30.
    #[test]
    fn the_fit_writes_a_per_bin_report() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(
            &fasta,
            format!(">chr1\n{}{}\n", "A".repeat(100), "G".repeat(100)),
        )
        .unwrap();
        let out = dir.path().join("gc.json.gz");
        let mut cov = vec![10u32; 100];
        cov.extend(vec![30u32; 100]);
        run_from_coverage(
            &params(&fasta, &out),
            HashMap::from([("chr1".to_string(), cov)]),
        )
        .unwrap();

        let report = std::fs::read_to_string(dir.path().join("gc.json.gz.bins.tsv")).unwrap();
        let rows: Vec<&str> = report.lines().collect();
        assert_eq!(rows[0], "gc_percent\twindows\tmedian\tmean\tweight\tsource");
        assert_eq!(rows.len(), 102, "a header and one row per GC percent");
        assert_eq!(rows[1], "0\t1\t10.000\t10.000\t0.500000\tfitted");
        assert_eq!(rows[101], "100\t1\t30.000\t30.000\t1.500000\tfitted");
        assert!(rows[51].ends_with("\tinterpolated"), "{}", rows[51]);
    }

    // A BAM that covers a small part of the reference (an exome or panel fitted against a
    // whole genome) has a median window at zero coverage. Every median-based weight is then
    // undefined, and the model would silently come out neutral. That is a zero denominator,
    // so it is refused and the message says how to fix it.
    #[test]
    fn a_fit_whose_median_window_is_uncovered_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "ACGT".repeat(250))).unwrap();
        let out = dir.path().join("gc.json.gz");
        let mut cov = vec![30u32; 100];
        cov.extend(vec![0u32; 900]);
        let err = run_from_coverage(
            &params(&fasta, &out),
            HashMap::from([("chr1".to_string(), cov)]),
        )
        .expect_err("a mostly uncovered fit must be refused");
        assert!(
            format!("{err}").contains("bed_file"),
            "the error must say how to fix it: {err}"
        );
        assert!(!out.exists(), "no model may be written");
    }

    // #771: a soft-masked (lowercase) base is the same base. Counting only uppercase G/C made
    // repeat windows read as lower GC than they are, so the same reference lowercased must give
    // the same model. Must-not-fire: the uppercase fit is the reference result.
    #[test]
    fn a_soft_masked_reference_gives_the_same_model() {
        let dir = tempfile::tempdir().unwrap();
        let seq = format!(
            "{}{}{}",
            "ACGT".repeat(25),
            "GGCA".repeat(25),
            "ATTA".repeat(25)
        );
        let cov: Vec<u32> = (0..300).map(|i| 10 + (i / 100) as u32 * 7).collect();
        let fit = |name: &str, s: &str| {
            let fasta = dir.path().join(format!("{name}.fa"));
            std::fs::write(&fasta, format!(">chr1\n{s}\n")).unwrap();
            let out = dir.path().join(format!("{name}.json.gz"));
            run_from_coverage(
                &params(&fasta, &out),
                HashMap::from([("chr1".to_string(), cov.clone())]),
            )
            .unwrap();
            GcBiasModel::from_file(&out).unwrap()
        };
        let upper = fit("upper", &seq);
        let lower = fit("lower", &seq.to_lowercase());
        for gc in 0..=100 {
            let f = gc as f64 / 100.0;
            assert_eq!(
                upper.weight_for_gc_fraction(f),
                lower.weight_for_gc_fraction(f),
                "GC {gc}%"
            );
        }
    }

    #[test]
    fn test_uniform_coverage_produces_neutral_weights() {
        // Reference: 300 bp of ACGT repeating (50% GC), uniform depth 10.
        // All populated bins should have relative weight ~1.0 (no bias).
        let seq: String = "ACGT".repeat(75);
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("uniform.bam");
        write_coverage_bam(&bam_path, &[(b"chr1", 300)], &[(0, 0, 300, 10)]);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        for gc_pct in 0..=100 {
            let w = model.weight_for_gc_fraction(gc_pct as f64 / 100.0);
            assert!(
                (w - 1.0).abs() < 1e-6,
                "Expected weight ~1.0 at GC {}%, got {}",
                gc_pct,
                w
            );
        }
    }

    #[test]
    fn test_gc_correlated_coverage_produces_correct_relative_weights() {
        // Window A: 200 bp AT (0% GC) at depth 10; Window B: 200 bp GC (100% GC)
        // at depth 30. Overall mean coverage = 20, so weight[0%] ~ 0.5 and
        // weight[100%] ~ 1.5.
        let seq = format!("{}{}", "AT".repeat(100), "GC".repeat(100));
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("two_window.bam");
        write_coverage_bam(
            &bam_path,
            &[(b"chr1", 400)],
            &[(0, 0, 200, 10), (0, 200, 400, 30)],
        );

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            200,
            200,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        assert!(
            (model.weight_for_gc_fraction(0.0) - 0.5).abs() < 1e-6,
            "Expected low-GC weight ~0.5, got {}",
            model.weight_for_gc_fraction(0.0)
        );
        assert!(
            (model.weight_for_gc_fraction(1.0) - 1.5).abs() < 1e-6,
            "Expected high-GC weight ~1.5, got {}",
            model.weight_for_gc_fraction(1.0)
        );
    }

    #[test]
    fn test_sparse_bin_gets_neutral_weight() {
        // One 100bp window at GC=50%, min_windows_per_bin=10 → all weights stay 1.0.
        let seq: String = "ACGT".repeat(25);
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("sparse.bam");
        write_coverage_bam(&bam_path, &[(b"chr1", 100)], &[(0, 0, 100, 20)]);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            10,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        assert_eq!(model.weight_for_gc_fraction(0.5), 1.0);
    }

    #[test]
    fn test_all_n_contig_is_skipped() {
        // chr1 is all-N (skipped); chr2 has real coverage. Runner must not panic
        // or error on the all-N contig.
        let fasta = ">chr1\nNNNNNNNNNN\n>chr2\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let ref_file = write_temp(fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("n_skip.bam");
        write_coverage_bam(
            &bam_path,
            &[(b"chr1", 10), (b"chr2", 104)],
            &[(0, 0, 10, 10), (1, 0, 104, 10)],
        );

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();
    }

    #[test]
    fn test_sparse_bin_excluded_from_overall_mean() {
        // 2 low-GC windows at depth 10 (well-populated, min_windows=2).
        // 1 high-GC window at depth 1000 (sparse — exactly 1, below min_windows=2).
        // If the sparse bin were folded into overall_mean it would inflate the mean
        // and crush low-GC weights well below 1.0. Excluding it as the code does
        // leaves low-GC weight at exactly 1.0.
        let low_gc: String = "AT".repeat(50);
        let high_gc: String = "GC".repeat(50);
        let seq = format!("{}{}{}", low_gc, low_gc, high_gc);
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("sparse_high.bam");
        write_coverage_bam(
            &bam_path,
            &[(b"chr1", 300)],
            &[(0, 0, 200, 10), (0, 200, 300, 1000)],
        );

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            100,
            100,
            2,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        let w_low = model.weight_for_gc_fraction(0.0);
        assert!(
            (w_low - 1.0).abs() < 1e-6,
            "Expected low-GC weight ~1.0, got {} (sparse high-GC bin may have inflated mean)",
            w_low
        );
        assert_eq!(
            model.weight_for_gc_fraction(1.0),
            1.0,
            "Expected sparse high-GC bin to have neutral weight 1.0"
        );
    }

    #[test]
    fn test_overlapping_windows_correct_gc_counts() {
        // window_size=4, window_stride=2 → overlapping windows on "AAAAGGGGAAAA".
        // Windows: [0,4)=0%GC, [2,6)=50%GC, [4,8)=100%GC, [6,10)=50%GC, [8,12)=0%GC.
        // Uniform depth → all populated bins should have weight ~1.0.
        let seq = "AAAAGGGGAAAA";
        let fasta = format!(">chr1\n{}\n", seq);
        let ref_file = write_temp(&fasta);

        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("overlap.bam");
        write_coverage_bam(&bam_path, &[(b"chr1", 12)], &[(0, 0, 12, 10)]);

        let out = NamedTempFile::new().unwrap();
        let out_path = out.path().to_path_buf();
        drop(out);

        let cfg = write_bam_config(
            &ref_file.path().to_path_buf(),
            &bam_path,
            &out_path,
            4,
            2,
            1,
        );
        runner(&cfg.path().to_path_buf()).unwrap();

        let model = GcBiasModel::from_file(&out_path).unwrap();
        for &gc_frac in &[0.0f64, 0.5, 1.0] {
            let w = model.weight_for_gc_fraction(gc_frac);
            assert!(
                (w - 1.0).abs() < 1e-6,
                "Expected weight ~1.0 at GC {:.0}%, got {}",
                gc_frac * 100.0,
                w
            );
        }
    }
}
