use crate::gen_seq_error_model::{
    errors::GenSeqErrorModelError,
    utils::config::{BamMethod, RunConfiguration},
};
use eidolon_core::file_tools::file_io::is_gzipped_file;
use eidolon_core::rng::NeatRng;
use eidolon_core::{
    file_tools::{
        bam_reader::{
            OVERLAP_BINS, OverlapCounts, read_bam_overlap_transitions, read_bam_transition_report,
        },
        file_io::{read_gzip_lines, read_lines},
        vcf_tools::read_known_sites,
    },
    models::{quality_scores::QualityScoreModel, sequencing_error_model::SequencingErrorModel},
    structs::transition_matrix::TransitionMatrix,
};
use log::{info, warn};
use std::path::PathBuf;

const MAX_SCORE: usize = 94;

/// Below this many observations, a position's transition row is noise rather than measurement.
/// Not a tuned figure -- a round number chosen to make a thin tail visible, which is the point.
const THIN_POSITION_OBSERVATIONS: usize = 100;

/// How many leading positions of a transition tensor were actually trained.
///
/// A position whose whole count matrix is zero was never observed. `build_distros` turns such a
/// row into a UNIFORM draw over the option set, which is a fabricated quality profile that
/// nothing downstream can distinguish from a fitted one -- so the fitter refuses to emit one
/// rather than reporting it. Counting from the front (rather than totalling zeros) is what the
/// caller needs: it converts directly to the longest read that population carried.
fn trained_positions(counts: &[Vec<Vec<usize>>]) -> usize {
    counts
        .iter()
        .position(|row| row.iter().flatten().sum::<usize>() == 0)
        .unwrap_or(counts.len())
}

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

/// Fewest counted errors an overlap fit accepts (#779). About 2,500 per row, which puts the
/// sampling error of a cell near 0.3 at ~0.009: comparable to the round trip's tolerance.
const MIN_OVERLAP_ERRORS: usize = 10_000;

/// Substitution counts from mate-overlap disagreements (#779), reported by category and refused
/// when there is too little evidence to fit a matrix from.
fn overlap_counts(
    path: &PathBuf,
    config: &RunConfiguration,
) -> Result<[[f64; 4]; 4], GenSeqErrorModelError> {
    let mask = match &config.known_variants_vcf {
        Some(vcf) => {
            let sites = read_known_sites(vcf)?;
            let n_sites: usize = sites.values().map(|s| s.len()).sum();
            info!("Masking {n_sites} reference positions from {vcf:?}");
            Some(sites)
        }
        None => None,
    };
    let c = read_bam_overlap_transitions(path, mask, config.bam_min_mapq)?;
    let errors: usize = c.counts.iter().flatten().sum();
    info!(
        "Overlap fit: {} overlapped bases; {} agree; {errors} errors counted; {} at masked sites \
         and {} with neither mate on the reference, both dropped; {} records without a mate",
        c.overlapped_bases, c.agreements, c.masked, c.neither_reference, c.unpaired
    );
    let mut report = config.output_file.clone().into_os_string();
    report.push(".overlap_bins.tsv");
    let report = PathBuf::from(report);
    write_overlap_bin_report(&report, &c)?;
    info!("Overlap errors by read position: {report:?}");
    if errors < MIN_OVERLAP_ERRORS {
        return Err(GenSeqErrorModelError::ConfigurationError(format!(
            "bam_file {path:?} gave {errors} sequencing errors from overlapping mates, fewer than \
             the {MIN_OVERLAP_ERRORS} a substitution matrix needs ({} overlapped bases). Overlap \
             fitting needs paired reads whose fragments are shorter than the two reads combined. \
             For a library without that, set `bam_method: mismatch` to count every \
             read-vs-reference mismatch instead, or remove bam_file to use the default matrix.",
            c.overlapped_bases
        )));
    }
    Ok(per_cycle_spectrum(&c))
}

/// One row per read-position bin (#779): the observations, errors and substitution counts
/// behind the pooled matrix, so a shift in the spectrum along the read is visible.
fn write_overlap_bin_report(path: &PathBuf, c: &OverlapCounts) -> std::io::Result<()> {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut body = String::from("bin\tcycles_pct\tbases\terrors\terror_rate");
    for f in 0..4 {
        for t in 0..4 {
            if f != t {
                body.push_str(&format!("\t{}>{}", BASES[f], BASES[t]));
            }
        }
    }
    body.push('\n');
    for (b, m) in c.by_bin.iter().enumerate() {
        let errors: usize = m.iter().flatten().sum();
        let bases = c.bin_bases[b];
        let rate = if bases > 0 {
            errors as f64 / bases as f64
        } else {
            0.0
        };
        let (lo, hi) = (b * 100 / OVERLAP_BINS, (b + 1) * 100 / OVERLAP_BINS);
        body.push_str(&format!("{b}\t{lo}-{hi}\t{bases}\t{errors}\t{rate:.6}"));
        for (f, row) in m.iter().enumerate() {
            for (t, n) in row.iter().enumerate() {
                if f != t {
                    body.push_str(&format!("\t{n}"));
                }
            }
        }
        body.push('\n');
    }
    std::fs::write(path, body)
}

/// Mismatch counts as weights, logging the count first (the mismatch path's only report).
fn as_weights(counts: [[usize; 4]; 4]) -> [[f64; 4]; 4] {
    let total: usize = counts.iter().flatten().sum();
    info!("Observed {total} SNP mismatches across all records");
    counts.map(|row| row.map(|n| n as f64))
}

/// Fewest errors a read-position bin needs to be weighted on its own; sparser bins merge with
/// their neighbors first, so one noisy bin cannot carry a full bin's weight.
const MIN_BIN_ERRORS: usize = 10_000;

/// The substitution mix of a whole read, from per-bin overlap counts (#779).
///
/// Overlaps sit at reads' 3' ends, so a bin near the end contributes far more observed bases
/// than one near the start (HG002: 8.8G against 0.7G). Pooling the raw counts gives the
/// overlap region's spectrum. Dividing each bin's counts by its observations gives per-base
/// rates, and summing those counts every part of the read equally: the errors a whole read
/// gets, which is what generation applies.
fn per_cycle_spectrum(c: &OverlapCounts) -> [[f64; 4]; 4] {
    // Group adjacent bins from the read's start until each group holds MIN_BIN_ERRORS; a
    // trailing remainder joins the last group. Each group is (counts, bases, bins spanned).
    let mut groups: Vec<([[usize; 4]; 4], usize, usize)> = Vec::new();
    let mut cur = ([[0usize; 4]; 4], 0usize, 0usize);
    for (bin, m) in c.by_bin.iter().enumerate() {
        for (from, row) in m.iter().enumerate() {
            for (to, &n) in row.iter().enumerate() {
                cur.0[from][to] += n;
            }
        }
        cur.1 += c.bin_bases[bin];
        cur.2 += 1;
        if cur.0.iter().flatten().sum::<usize>() >= MIN_BIN_ERRORS {
            groups.push(std::mem::take(&mut cur));
        }
    }
    if cur.2 > 0 {
        match groups.last_mut() {
            Some(last) => {
                for from in 0..4 {
                    for to in 0..4 {
                        last.0[from][to] += cur.0[from][to];
                    }
                }
                last.1 += cur.1;
                last.2 += cur.2;
            }
            None => groups.push(cur),
        }
    }
    // Per-base rate times the number of bins the group spans: every cycle counts equally.
    let mut w = [[0.0f64; 4]; 4];
    for (counts, bases, width) in groups {
        if bases == 0 {
            continue;
        }
        for from in 0..4 {
            for to in 0..4 {
                w[from][to] += counts[from][to] as f64 / bases as f64 * width as f64;
            }
        }
    }
    w
}

/// Normalizes a raw 4×4 mismatch count matrix into a `TransitionMatrix`.
///
/// Each row is normalized independently. Rows with no observed mismatches get
/// equal probability distributed across the three off-diagonal positions.
#[cfg(test)]
fn build_transition_matrix_from_counts(
    counts: [[usize; 4]; 4],
) -> Result<TransitionMatrix, GenSeqErrorModelError> {
    build_transition_matrix_from_weights(counts.map(|row| row.map(|n| n as f64)))
}

/// As `build_transition_matrix_from_counts`, from non-negative weights rather than counts:
/// the overlap fit's per-cycle rates (#779) are not integers.
fn build_transition_matrix_from_weights(
    counts: [[f64; 4]; 4],
) -> Result<TransitionMatrix, GenSeqErrorModelError> {
    let mut weights = [[0.0f64; 4]; 4];
    for i in 0..4 {
        let total: f64 = counts[i].iter().sum::<f64>();
        if total == 0.0 {
            for j in 0..4 {
                if i != j {
                    weights[i][j] = 1.0 / 3.0;
                }
            }
        } else {
            for j in 0..4 {
                weights[i][j] = counts[i][j] / total;
            }
            weights[i][i] = 0.0;
        }
    }
    Ok(TransitionMatrix::from(
        weights[0], weights[1], weights[2], weights[3],
    )?)
}

/// Shorten a line for an error message. A corrupt "line" can be megabytes.
fn truncate_for_message(line: &str) -> String {
    const MAX: usize = 60;
    if line.chars().count() <= MAX {
        return line.to_string();
    }
    let head: String = line.chars().take(MAX).collect();
    format!("{head}...")
}

/// Read a non-header line of a record. Unlike the header, running out of input here means the
/// final record is truncated -- a corrupt file -- rather than a clean end.
fn next_record_line<I>(
    iter: &mut I,
    first_line: usize,
    this_line: usize,
) -> Result<String, GenSeqErrorModelError>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    match iter.next() {
        Some(Ok(l)) => Ok(l),
        Some(Err(e)) => Err(e.into()),
        None => Err(GenSeqErrorModelError::MalformedFastq(format!(
            "the record starting at line {first_line} is truncated: the file ends before line \
             {this_line}. A FASTQ record is four lines (header, sequence, '+', quality)."
        ))),
    }
}

/// Structural validation of one FASTQ record.
///
/// WHY THIS EXISTS (#700): the reader used to consume records in groups of four and keep only
/// the fourth line, discarding the other three unexamined. Any file whose line count is a
/// multiple of four was therefore accepted -- including a FASTA, whose SEQUENCE lines then
/// landed in the quality accumulator. Every nucleotide decodes to a plausible Phred score under
/// the offset-33 arithmetic (A->Q32, C->Q34, G->Q38, T->Q51), so the run completed at exit 0 and
/// wrote a model with an ordinary-looking error_rate. Nothing downstream could tell it apart
/// from a model fitted on reads.
///
/// These are the structural checks `eidolon validate` already applies, minus the per-base IUPAC
/// test, which is O(n) in the hot loop and buys little once '@' is enforced. That leaves one
/// known gap: a file with its sequence and quality lines transposed passes all three checks.
/// `eidolon validate` is the place for that.
fn validate_record(
    header: &str,
    seq: &str,
    plus: &str,
    qual: &str,
    first_line: usize,
) -> Result<(), GenSeqErrorModelError> {
    if !header.starts_with('@') {
        // The common mistake gets its own message rather than the generic one.
        if header.starts_with('>') {
            return Err(GenSeqErrorModelError::MalformedFastq(format!(
                "line {first_line} starts with '>', so this looks like a FASTA, not a FASTQ: \
                 {:?}. A sequencing-error model is fitted from quality scores, which a FASTA \
                 does not carry.",
                truncate_for_message(header)
            )));
        }
        return Err(GenSeqErrorModelError::MalformedFastq(format!(
            "line {first_line}: a FASTQ record header must start with '@', found {:?}",
            truncate_for_message(header)
        )));
    }
    if !plus.starts_with('+') {
        return Err(GenSeqErrorModelError::MalformedFastq(format!(
            "line {}: the third line of a FASTQ record must start with '+', found {:?}",
            first_line + 2,
            truncate_for_message(plus)
        )));
    }
    if seq.len() != qual.len() {
        return Err(GenSeqErrorModelError::MalformedFastq(format!(
            "line {}: the quality string is {} character(s) but the sequence on line {} is {}",
            first_line + 3,
            qual.len(),
            first_line + 1,
            seq.len()
        )));
    }
    Ok(())
}

/// Everything ONE mate's fit accumulates.
///
/// `global_counts` and `total_bases` are deliberately NOT here: they stay pooled across both
/// mates, so `quality_score_options` is the UNION over the pair. That is the same constraint
/// the healthy and degraded populations already live under -- `generate_quality_scores`
/// indexes a single option list, and a tensor built against its own options indexes the wrong
/// scores (#723, and NEAT2 enforced the parallel rule at `SequenceContainer.py:678`).
struct MateCounts {
    seed_counts: Vec<usize>,
    transition_counts: Vec<Vec<Vec<usize>>>,
    seed_counts_deg: Vec<usize>,
    transition_counts_deg: Vec<Vec<Vec<usize>>>,
    n_degraded: usize,
    n_healthy: usize,
    n_unclassifiable: usize,
    n_deg_low_cut: usize,
    n_deg_high_cut: usize,
    records_seen: usize,
    min_qual_len: usize,
}

impl MateCounts {
    fn new() -> Self {
        Self {
            seed_counts: vec![0usize; MAX_SCORE],
            transition_counts: Vec::new(),
            seed_counts_deg: vec![0usize; MAX_SCORE],
            transition_counts_deg: Vec::new(),
            n_degraded: 0,
            n_healthy: 0,
            n_unclassifiable: 0,
            n_deg_low_cut: 0,
            n_deg_high_cut: 0,
            records_seen: 0,
            min_qual_len: usize::MAX,
        }
    }

    /// Positions this mate's healthy tensor covers. The degraded one is checked separately.
    fn positions(&self) -> usize {
        self.transition_counts
            .len()
            .max(self.transition_counts_deg.len())
    }
}

/// Accumulate one quality line into the counts.
///
/// Every allocated position carries at least one observation BY CONSTRUCTION: a row exists
/// only because some read reached that position, and that read necessarily incremented it.
/// An untrained position is therefore unrepresentable rather than merely unlikely, which is
/// what `every_position_has_an_observation` pins.
fn accumulate_qual(
    qual_line: &str,
    qual_offset: usize,
    max_model_read_length: usize,
    bins: Option<&[usize]>,
    seed_counts: &mut [usize],
    transition_counts: &mut Vec<Vec<Vec<usize>>>,
    global_counts: &mut [usize],
    total_bases: &mut usize,
) -> Result<bool, GenSeqErrorModelError> {
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
        return Ok(false);
    }
    // Growing the tensor removed truncation, and truncation was the only bound on the
    // allocation. This restores a bound WITHOUT restoring silent data loss: a read above
    // the ceiling stops the run and says so, rather than being quietly trimmed (#697).
    if max_model_read_length > 0 && scores.len() > max_model_read_length {
        let mib =
            scores.len().saturating_sub(1) * MAX_SCORE * MAX_SCORE * std::mem::size_of::<usize>()
                / (1024 * 1024);
        return Err(GenSeqErrorModelError::ConfigurationError(format!(
            "a read is {} bp, above the {} bp limit set by `max_model_read_length`; \
                 modeling it would allocate ~{} MiB of transition counts. Raise that key (or \
                 set it to 0 for no limit) if the input really is this long. Long-read error \
                 models are tracked in #319.",
            scores.len(),
            max_model_read_length,
            mib,
        )));
    }
    seed_counts[scores[0]] += 1;
    if scores.len() > 1 && transition_counts.len() < scores.len() - 1 {
        transition_counts.resize(scores.len() - 1, vec![vec![0usize; MAX_SCORE]; MAX_SCORE]);
    }
    for j in 1..scores.len() {
        transition_counts[j - 1][scores[j - 1]][scores[j]] += 1;
    }
    for &score in &scores {
        global_counts[score] += 1;
    }
    *total_bases += scores.len();
    Ok(true)
}

/// Fixed so a capped fit is reproducible, and shared by both mates so R1 and R2 of the same
/// length keep the same record numbers.
const MAX_READS_SAMPLING_SEED: &str = "gen-seq-error-model max_reads (#721)";

/// Which records a `max_reads`-capped fit keeps (#721).
///
/// The cap used to keep the FIRST N records. Illumina FASTQs are written in flowcell order,
/// so that was one corner of one lane: the first 200,000 records of the HG002 library were
/// all lane 1, tile 1101. Each record is now kept with probability `max_reads / total`, so
/// the sample is spread across the whole file and holds about `max_reads` records.
struct ReadSampler {
    /// `None` keeps every record: no cap, or a cap at or above the file's size.
    keep_probability: Option<f64>,
    rng: NeatRng,
}

impl ReadSampler {
    fn new(max_reads: usize, total_records: usize) -> Result<Self, GenSeqErrorModelError> {
        let keep_probability = (max_reads > 0 && total_records > max_reads)
            .then(|| max_reads as f64 / total_records as f64);
        let rng = NeatRng::new_from_seed(&vec![MAX_READS_SAMPLING_SEED.to_string()])
            .map_err(|e| GenSeqErrorModelError::ConfigurationError(e.to_string()))?;
        Ok(Self {
            keep_probability,
            rng,
        })
    }

    fn keep(&mut self) -> Result<bool, GenSeqErrorModelError> {
        match self.keep_probability {
            None => Ok(true),
            Some(p) => self
                .rng
                .gen_bool(p)
                .map_err(|e| GenSeqErrorModelError::ConfigurationError(e.to_string())),
        }
    }
}

/// Records in a FASTQ, counted by lines. A first pass for `max_reads`, which needs the total
/// to sample uniformly; the fitting pass validates every record, so a malformed file is
/// still reported there.
fn count_fastq_records(path: &PathBuf) -> Result<usize, GenSeqErrorModelError> {
    let lines = if is_gzipped_file(path)? {
        read_gzip_lines(path)?.count()
    } else {
        read_lines(path)?.count()
    };
    Ok(lines / 4)
}

/// Read one mate's FASTQ into its own `MateCounts`, pooling the option-set and base counts.
///
/// `max_reads` applies PER MATE rather than across the pair: a shared budget would spend it
/// all on R1 and fit R2 from whatever was left, which is nothing at the default of 0-means-all
/// and is silently lopsided otherwise. Within a mate it is a uniform sample, not the head of
/// the file; see `ReadSampler`.
fn accumulate_one_file(
    path: &PathBuf,
    config: &RunConfiguration,
    bins_slice: Option<&[usize]>,
    m: &mut MateCounts,
    global_counts: &mut [usize],
    total_bases: &mut usize,
    reads_processed: &mut usize,
) -> Result<(), GenSeqErrorModelError> {
    let mut iter: Box<dyn Iterator<Item = std::io::Result<String>>> = if is_gzipped_file(path)? {
        Box::new(read_gzip_lines(path)?)
    } else {
        Box::new(read_lines(path)?)
    };
    let qual_offset = config.qual_offset;
    let mut sampler = if config.max_reads > 0 {
        let total = count_fastq_records(path)?;
        let sampler = ReadSampler::new(config.max_reads, total)?;
        match sampler.keep_probability {
            Some(p) => info!(
                "max_reads {}: sampling {:?} uniformly, keeping each of its {} records with \
                 probability {:.6}",
                config.max_reads, path, total, p
            ),
            None => info!(
                "max_reads {} is at least the {} records in {:?}; using every record",
                config.max_reads, total, path
            ),
        }
        sampler
    } else {
        ReadSampler::new(0, 0)?
    };
    // Every record read, kept or not: line numbers in error messages count from this, while
    // `m.records_seen` counts only the records the fit used.
    let mut records_read = 0usize;
    'records: loop {
        // Line number of this record's header, 1-based, for error messages.
        let first_line = records_read * 4 + 1;

        // The header is the ONLY line whose absence is a clean end of file. Missing any of the
        // other three means the last record is truncated, which is a corrupt file, not an end.
        let header = match iter.next() {
            None => break 'records,
            Some(Ok(l)) => l,
            Some(Err(e)) => return Err(e.into()),
        };
        let seq = next_record_line(&mut iter, first_line, first_line + 1)?;
        let plus = next_record_line(&mut iter, first_line, first_line + 2)?;
        let qual = next_record_line(&mut iter, first_line, first_line + 3)?;

        validate_record(&header, &seq, &plus, &qual, first_line)?;
        records_read += 1;
        if !sampler.keep()? {
            continue;
        }
        m.records_seen += 1;
        m.min_qual_len = m.min_qual_len.min(qual.len());

        // Classify BEFORE accumulating, so the read's counts go to one population or the other.
        // A read shorter than the window cannot be classified; it is counted as healthy and
        // reported separately rather than silently folded in.
        let window = config.degradation_tail_window;
        // Tracked separately from `degraded`: a read too short to classify is accumulated into
        // the healthy tensor, because its bases are still data, but it must NOT enter the
        // denominator of the reported fraction. Counting it as healthy understated the fraction
        // (30 of 120 rather than 30 of 100) until a fixture with short reads caught it.
        let mut classifiable = config.fit_quality_degradation;
        let degraded = if !config.fit_quality_degradation {
            false
        } else if qual.len() < window {
            m.n_unclassifiable += 1;
            classifiable = false;
            false
        } else {
            let tail: f64 = qual.as_bytes()[qual.len() - window..]
                .iter()
                .map(|&b| (b as usize).saturating_sub(qual_offset) as f64)
                .sum::<f64>()
                / window as f64;
            if tail < (config.degradation_tail_cut as f64 - 2.0) {
                m.n_deg_low_cut += 1;
            }
            if tail < (config.degradation_tail_cut as f64 + 2.0) {
                m.n_deg_high_cut += 1;
            }
            tail < config.degradation_tail_cut as f64
        };
        if classifiable {
            if degraded {
                m.n_degraded += 1;
            } else {
                m.n_healthy += 1;
            }
        }

        let MateCounts {
            ref mut seed_counts,
            ref mut transition_counts,
            ref mut seed_counts_deg,
            ref mut transition_counts_deg,
            ..
        } = *m;
        let (seed_target, trans_target) = if degraded {
            (seed_counts_deg, transition_counts_deg)
        } else {
            (seed_counts, transition_counts)
        };
        if accumulate_qual(
            &qual,
            qual_offset,
            config.max_model_read_length,
            bins_slice,
            seed_target,
            trans_target,
            global_counts,
            total_bases,
        )? {
            (*reads_processed) += 1;
        }
    }
    if sampler.keep_probability.is_some() {
        info!(
            "max_reads {}: fitted {} of {} records in {:?}",
            config.max_reads, m.records_seen, records_read, path
        );
    }
    Ok(())
}

pub fn runner(config: &RunConfiguration) -> Result<(), GenSeqErrorModelError> {
    let mut global_counts = vec![0usize; MAX_SCORE];
    let mut total_bases: usize = 0;
    let mut reads_processed: usize = 0;
    let bins_slice: Option<&[usize]> = config.binned_quality_bins.as_deref();

    // One pass per mate, sharing `global_counts` so the option set is the union over the pair.
    // Fitting them separately and stapling the results together would give each its own option
    // list, which generation cannot index.
    let mut r1 = MateCounts::new();
    accumulate_one_file(
        &config.fastq_file,
        config,
        bins_slice,
        &mut r1,
        &mut global_counts,
        &mut total_bases,
        &mut reads_processed,
    )?;
    let r2 = match &config.fastq_file_r2 {
        Some(path) => {
            let mut m = MateCounts::new();
            accumulate_one_file(
                path,
                config,
                bins_slice,
                &mut m,
                &mut global_counts,
                &mut total_bases,
                &mut reads_processed,
            )?;
            Some(m)
        }
        None => None,
    };

    // Keep the R1 names the rest of this function already uses.
    let MateCounts {
        seed_counts,
        transition_counts,
        seed_counts_deg,
        transition_counts_deg,
        n_degraded,
        n_healthy,
        n_unclassifiable,
        n_deg_low_cut,
        n_deg_high_cut,
        records_seen,
        min_qual_len,
    } = r1;

    if records_seen == 0 {
        return Err(GenSeqErrorModelError::MalformedFastq(
            "FASTQ file has fewer than 4 lines".to_string(),
        ));
    }

    // The fit's OUTPUT, not an input to it: one row per position after the first. Both
    // populations describe the same read length, because `generate_quality_scores` indexes one
    // position list for either of them.
    let positions = transition_counts
        .len()
        .max(transition_counts_deg.len())
        .max(r2.as_ref().map_or(0, |m| m.positions()));
    let read_length = positions + 1;

    // Rule 4: a rate over an unknown denominator is not a result, and an empty population is a
    // hard failure rather than a warning. A fit that was ASKED for two populations and found
    // one has not produced a two-population model; saying so beats emitting a model whose
    // degraded tensor is entirely uniform fallback.
    //
    // These run BEFORE the coverage check below. An empty population has an empty tensor, so it
    // trips that check too, and "covers 1 bp of 100" is a true but unhelpful way to say "this
    // library has no degraded reads at this cut".
    let classified = n_degraded + n_healthy;
    if config.fit_quality_degradation {
        if classified == 0 {
            return Err(GenSeqErrorModelError::MalformedFastq(
                "no reads could be classified, so no degraded population could be fitted"
                    .to_string(),
            ));
        }
        if n_degraded == 0 {
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "fit_quality_degradation is set, but none of {classified} classified reads fall \
                 below the Q{} tail cut, so there is no degraded population to fit. Either this \
                 library has none, or `degradation_tail_cut` is too low. Unset \
                 fit_quality_degradation to fit a single population.",
                config.degradation_tail_cut
            )));
        }
        if n_healthy == 0 {
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "fit_quality_degradation is set, but all {classified} classified reads fall \
                 below the Q{} tail cut, so there is no healthy population to fit. \
                 `degradation_tail_cut` is above this library's quality range. Unset \
                 fit_quality_degradation to fit a single population.",
                config.degradation_tail_cut
            )));
        }
    }

    // The shorter population is REFUSED, not padded.
    //
    // `every_modeled_position_is_trained` pins the invariant for a single population: a
    // position row exists only because some read reached it, and that read necessarily trained
    // it. A second population breaks it -- padding the shorter tensor to the common length
    // leaves all-zero rows, and `build_distros` turns an all-zero row into a UNIFORM draw over
    // the option set. 60 bp degraded reads against 100 bp healthy ones would then hand the
    // degraded population invented quality from base 62 to base 100, and no consumer of the
    // model could tell that from a fit. Rule 4: a zero denominator is a hard failure, not a
    // warning, and here it is a hard failure BEFORE the number is fabricated.
    //
    // Checked on observations rather than on tensor lengths. The two agree today -- a read of
    // length L trains every position below L -- but the length is a proxy and the observation
    // count is the thing that matters, so an untrained position cannot reappear through some
    // later change to how counts are accumulated.
    {
        // The remedy names the axis the populations came from. A run that never set
        // `fit_quality_degradation` was once told to unset it, because one message served
        // every label.
        const DEGRADED_REMEDY: &str = "trim or filter the library so the two carry the same \
                                       read lengths, or unset fit_quality_degradation to fit \
                                       a single population.";
        const MATE_REMEDY: &str = "trim or filter the two mate files so they carry the same \
                                   read lengths, or drop fastq_file_r2 to fit one quality \
                                   model for both mates.";
        const MATE_DEGRADED_REMEDY: &str = "trim or filter the two mate files so they carry \
                                            the same read lengths, or unset \
                                            fit_quality_degradation to fit one population per \
                                            mate.";

        // (label, transition counts, seed counts, remedy). The seed is position 1, which the
        // transition tensor does not hold: a population with no seed covers 0 bp. Counting
        // transitions alone could not tell that apart from 1 bp, so with 1 bp reads, where
        // there are no transitions at all, the check never fired (#767).
        let mut populations: Vec<(&str, &Vec<Vec<Vec<usize>>>, &[usize], &str)> = Vec::new();
        // R1's healthy tensor is checked only when there is a second population to be ragged
        // against. On its own it defines `positions` and cannot fall short of itself.
        if config.fit_quality_degradation {
            populations.push(("healthy", &transition_counts, &seed_counts, DEGRADED_REMEDY));
            populations.push((
                "degraded",
                &transition_counts_deg,
                &seed_counts_deg,
                DEGRADED_REMEDY,
            ));
        }
        if let Some(m) = r2.as_ref() {
            populations.push(("R1", &transition_counts, &seed_counts, MATE_REMEDY));
            populations.push(("R2", &m.transition_counts, &m.seed_counts, MATE_REMEDY));
            if config.fit_quality_degradation {
                populations.push((
                    "R2 degraded",
                    &m.transition_counts_deg,
                    &m.seed_counts_deg,
                    MATE_DEGRADED_REMEDY,
                ));
            }
        }
        for (label, counts, seeds, remedy) in populations {
            let covered = if seeds.iter().sum::<usize>() == 0 {
                0
            } else {
                trained_positions(counts) + 1
            };
            if covered < read_length {
                return Err(GenSeqErrorModelError::ConfigurationError(format!(
                    "the {label} population covers {covered} bp, but the model covers \
                     {read_length} bp: positions {} to {read_length} carry no {label} \
                     observation at all. Filling them would give that population a uniform \
                     quality distribution -- fabricated data a reader cannot distinguish \
                     from a fit. Every population must reach the model's read length: \
                     {remedy}",
                    covered + 1,
                )));
            }
        }
    }

    info!(
        "Processed {} reads ({} bases)",
        reads_processed, total_bases
    );
    info!(
        "Model read length: {} bp (maximum over {} record(s) accumulated)",
        read_length, reads_processed
    );
    if min_qual_len != usize::MAX && min_qual_len != read_length {
        warn!(
            "Read lengths vary across the {} record(s) fitted: {}-{} bp. Every position is \
             modeled; shorter reads contribute only the positions they cover.",
            reads_processed, min_qual_len, read_length
        );
    }

    // Rule 4: report the denominator, not just the metric. Every position holds at least one
    // observation -- by construction for a single population, and by the refusal above for
    // two -- but "at least one" is not "enough": a thin tail makes the late positions noise,
    // and a small `max_reads` is the usual way to get one. Reported per population, because
    // the degraded one is a minority of the library and so thins out first.
    for (label, counts) in [
        ("healthy population", &transition_counts),
        ("degraded population", &transition_counts_deg),
    ] {
        let Some(min_obs) = counts
            .iter()
            .map(|row| row.iter().flatten().sum::<usize>())
            .min()
        else {
            continue; // no degraded population was fitted
        };
        info!("Fewest transition observations at any position ({label}): {min_obs}");
        if min_obs < THIN_POSITION_OBSERVATIONS {
            warn!(
                "Some modeled position of the {label} rests on only {} observation(s) \
                 (threshold {}). Those positions are noise rather than measurement; fit more \
                 reads, or accept that the tail of the model is thin.",
                min_obs, THIN_POSITION_OBSERVATIONS
            );
        }
    }

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

    // Re-index raw MAX_SCORE-indexed counts onto the shared option list. Every population --
    // R1, its degraded half, R2 and its degraded half -- goes through THESE two closures, so
    // they cannot drift in how they are built. A pair built two different ways would index
    // different scores while looking identical.
    let reindex_seed = |counts: &[usize]| -> Vec<f64> {
        quality_score_options
            .iter()
            .map(|&q| counts[q] as f64)
            .collect()
    };
    let reindex_trans = |counts: &[Vec<Vec<usize>>]| -> Vec<Vec<Vec<f64>>> {
        (0..read_length - 1)
            .map(|pos| {
                (0..n_scores)
                    .map(|p| {
                        let prev_raw = quality_score_options[p];
                        (0..n_scores)
                            .map(|c| {
                                let curr_raw = quality_score_options[c];
                                counts[pos][prev_raw][curr_raw] as f64
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect()
    };

    let seed_weights = reindex_seed(&seed_counts);
    let trans_weights = reindex_trans(&transition_counts);

    let quality_score_model = QualityScoreModel::from_counts(
        quality_score_options.clone(),
        read_length,
        seed_weights,
        trans_weights,
        is_binned,
    )?;

    // Attach the degraded population, re-indexed against the SAME option set.
    let quality_score_model = if config.fit_quality_degradation {
        let read_fraction = n_degraded as f64 / classified as f64;

        let seed_weights_deg = reindex_seed(&seed_counts_deg);
        let trans_weights_deg = reindex_trans(&transition_counts_deg);

        info!(
            "Degraded population: {} of {} classified reads ({:.2}%), {} unclassifiable (shorter \
             than the {} base window)",
            n_degraded,
            classified,
            100.0 * read_fraction,
            n_unclassifiable,
            config.degradation_tail_window,
        );
        // The separability check #694 asks for, reported rather than assumed. If moving the cut
        // by two Q moves the fraction by a large factor, the cut is slicing one distribution
        // rather than separating two.
        let (lo, hi) = (
            100.0 * n_deg_low_cut as f64 / classified as f64,
            100.0 * n_deg_high_cut as f64 / classified as f64,
        );
        info!(
            "Separability: Q{} -> {:.2}%, Q{} -> {:.2}%, Q{} -> {:.2}%. A large swing across \
             this range means the two populations are not cleanly separable at this cut.",
            config.degradation_tail_cut - 2,
            lo,
            config.degradation_tail_cut,
            100.0 * read_fraction,
            config.degradation_tail_cut + 2,
            hi,
        );

        quality_score_model.with_degradation(read_fraction, seed_weights_deg, trans_weights_deg)?
    } else {
        quality_score_model
    };

    // Determine transition matrix: TSV > BAM inference > defaults from original Python NEAT
    // (see https://github.com/ncsa/NEAT/blob/main/neat/model_sequencing_error/runner.py)
    let snp_transition_matrix: Option<TransitionMatrix> = if let Some(path) =
        &config.transition_matrix_file
    {
        info!("Loading SNP transition matrix from TSV: {:?}", path);
        Some(TransitionMatrix::from_tsv(path)?)
    } else if let Some(path) = &config.bam_file {
        info!("Inferring SNP transition matrix from BAM: {:?}", path);
        let mut md_seq_disagreements = 0usize;
        let counts = match config.bam_method {
            BamMethod::Overlap => overlap_counts(path, config)?,
            BamMethod::Mismatch => {
                let mask = match &config.known_variants_vcf {
                    None => None,
                    Some(vcf) => {
                        let sites = read_known_sites(vcf)?;
                        let n_sites: usize = sites.values().map(|s| s.len()).sum();
                        info!("Masking {n_sites} reference positions from {vcf:?}");
                        Some(sites)
                    }
                };
                let masking = mask.is_some();
                let report = read_bam_transition_report(path, mask)?;
                if masking {
                    let kept: usize = report.counts.iter().flatten().sum();
                    let seen = kept + report.masked;
                    info!(
                        "Masked {} of {seen} mismatches ({:.2}%) at known variant sites",
                        report.masked,
                        if seen > 0 {
                            100.0 * report.masked as f64 / seen as f64
                        } else {
                            0.0
                        }
                    );
                }
                md_seq_disagreements = report.md_seq_disagreements;
                if md_seq_disagreements > 0 {
                    warn!(
                        "{md_seq_disagreements} position(s) in {path:?} are MD mismatches where \
                             the read has the reference base: the MD tags disagree with SEQ and \
                             were left out. `samtools calmd -b` rewrites MD from the reference."
                    );
                }
                as_weights(report.counts)
            }
        };
        // Every row empty means the BAM yielded no evidence at all. The overlap path has
        // already refused thin evidence and logged its counts; this is the mismatch path's
        // guard, and a backstop for both.
        let total: f64 = counts.iter().flatten().sum();
        if total == 0.0 && md_seq_disagreements > 0 {
            // Not the missing-MD case below: the tags are there but contradict the reads,
            // which used to surface as a bare InvalidWeights (#767).
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "bam_file {:?} yielded no usable mismatches: its MD tags disagree with \
                     the reads at all {md_seq_disagreements} position(s) they name as \
                     mismatches, where the read has the reference base. The tags are stale, \
                     for instance written before the reads were edited. Rewrite them with \
                     `samtools calmd -b {} reference.fa > fixed.bam`.",
                path,
                path.display()
            )));
        }
        if total == 0.0 {
            // Hard error, not a warning. `bam_file:` exists for exactly one purpose --
            // fitting the transition matrix -- so a BAM that yields no evidence cannot be
            // honoured at all. Silently substituting the default produces a model that is
            // indistinguishable from a trained one downstream. MD is a *predefined* SAM tag
            // rather than a required one, so this is easy to hit unintentionally.
            //
            // Omitting `bam_file:` is how a user asks for the default matrix, so there is
            // nothing to fall back to here: a config that names a BAM it does not use would
            // be a lie about provenance.
            return Err(GenSeqErrorModelError::ConfigurationError(format!(
                "bam_file {:?} yielded no read-vs-reference mismatches, so the SNP \
                     transition matrix cannot be inferred from it. Most likely the BAM has no \
                     MD tags, which are optional in SAM and not written by every aligner. \
                     Either add them with `samtools calmd -b {} reference.fa > with_md.bam`, \
                     or remove `bam_file:` from the config to use the built-in default matrix.",
                path,
                path.display()
            )));
        } else {
            Some(build_transition_matrix_from_weights(counts)?)
        }
    } else {
        None
    };

    let model = SequencingErrorModel::from_raw_data(
        error_rate,
        quality_score_model,
        snp_transition_matrix,
    )?;

    // Attach R2, re-indexed against the SAME option set as R1 (#723).
    let model = match r2 {
        None => model,
        Some(m) => {
            let r2_model = QualityScoreModel::from_counts(
                quality_score_options.clone(),
                read_length,
                reindex_seed(&m.seed_counts),
                reindex_trans(&m.transition_counts),
                is_binned,
            )?;
            let r2_model = if config.fit_quality_degradation {
                let classified_r2 = m.n_degraded + m.n_healthy;
                info!(
                    "R2 degraded population: {} of {} classified reads ({:.2}%), {} \
                     unclassifiable",
                    m.n_degraded,
                    classified_r2,
                    100.0 * m.n_degraded as f64 / classified_r2 as f64,
                    m.n_unclassifiable,
                );
                r2_model.with_degradation(
                    m.n_degraded as f64 / classified_r2 as f64,
                    reindex_seed(&m.seed_counts_deg),
                    reindex_trans(&m.transition_counts_deg),
                )?
            } else {
                r2_model
            };
            info!(
                "Fitted a separate R2 quality model over {} record(s)",
                m.records_seen
            );
            model.with_mate_r2(r2_model)?
        }
    };
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
    let records = vec![(read_bases, md_str.as_deref()); n_records];
    write_test_bam_records(path, &records);
}

/// A minimal BAM of `(SEQ, MD)` records, each fully aligned at position 1 of `chr1`. The MD
/// tag is written verbatim rather than derived from a reference, so a test can write one
/// that disagrees with SEQ (#767). `None` omits the tag.
#[cfg(test)]
fn write_test_bam_records(path: &PathBuf, records: &[(&[u8], Option<&str>)]) {
    use noodles::bam;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
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
    let header = sam::Header::builder()
        .add_reference_sequence(
            b"chr1".to_vec(),
            Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000).unwrap()),
        )
        .build();

    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    for &(read_bases, md_str) in records {
        let cigar: Cigar = [Op::new(Kind::Match, read_bases.len())]
            .into_iter()
            .collect();
        let mut record = RecordBuf::default();
        *record.flags_mut() = Flags::empty();
        *record.cigar_mut() = cigar;
        *record.sequence_mut() = Sequence::from(read_bases);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = noodles::core::Position::new(1);

        if let Some(md) = md_str {
            record
                .data_mut()
                .insert(Tag::MISMATCHED_POSITIONS, BufValue::from(md));
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
            fastq_file_r2: None,
            output_file: output,
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: None,
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
        }
    }

    /// A FASTQ whose FIRST record is short and whose remainder is full length.
    ///
    /// Quality is planted so the answer is known without consulting the code: the long reads
    /// are `hi` for their first half and `lo` for their second. Under the #697 bug the model's
    /// read length came from the short first record, `accumulate_qual` truncated every long
    /// read to it, and the entire `lo` half was never seen by the transition tensor — so a
    /// model built from this file could not emit a `lo` score at any position.
    fn make_short_first_fastq(
        path: &PathBuf,
        n_long: usize,
        short_len: usize,
        long_len: usize,
        hi: u8,
        lo: u8,
    ) {
        use std::io::Write;
        let mut out = String::new();
        let seq = |n: usize| -> String { "ACGT".chars().cycle().take(n).collect() };
        // The offending record: short, and first.
        out.push_str(&format!(
            "@short\n{}\n+\n{}\n",
            seq(short_len),
            (hi as char).to_string().repeat(short_len)
        ));
        let half = long_len / 2;
        let qual: String = (hi as char).to_string().repeat(half)
            + &(lo as char).to_string().repeat(long_len - half);
        for i in 0..n_long {
            out.push_str(&format!("@long{}\n{}\n+\n{}\n", i, seq(long_len), qual));
        }
        let mut f = std::fs::File::create(path).unwrap();
        f.write_all(out.as_bytes()).unwrap();
    }

    /// THE regression for #697. Fitting HG002 R2 produced a 168-position model from 250 bp
    /// reads because one short record led the file; a third of every read was discarded, and
    /// specifically the 3' end where quality degrades.
    #[test]
    fn read_length_is_not_taken_from_a_short_first_record() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("mixed.fastq");
        // Q40 = 'I', Q10 = '+' at offset 33. Long reads are Q40 for 25 bases then Q10 for 25.
        make_short_first_fastq(&fastq_path, 200, 10, 50, b'I', b'+');
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(
            q.assumed_read_length, 50,
            "read length must come from the reads, not from the short first record"
        );
        assert_eq!(
            q.distros_from_one.len(),
            49,
            "one transition row per position after the first"
        );

        // The known answer: the second half of every long read is Q10, so a model that saw
        // those positions must emit Q10 there. Under the bug it emits Q40 everywhere, because
        // positions 11-50 contributed no transitions at all.
        let mut rng = eidolon_core::rng::NeatRng::new_from_seed(&vec!["s697".to_string()]).unwrap();
        let scores = q.generate_quality_scores(50, &mut rng).unwrap();
        assert_eq!(scores.len(), 50);
        let late_low = scores[30..].iter().filter(|&&s| s == 10).count();
        assert!(
            late_low >= 15,
            "positions 31-50 were all Q10 in training; got {} of 20 at Q10: {:?}",
            late_low,
            &scores[30..]
        );
        // Must not fire: the FIRST half was Q40 and must not have been overwritten by it.
        let early_high = scores[..25].iter().filter(|&&s| s == 40).count();
        assert!(
            early_high >= 20,
            "positions 1-25 were Q40 in training; got {} of 25: {:?}",
            early_high,
            &scores[..25]
        );
    }

    /// Must not fire: a uniform-length file is unchanged by the sampling change.
    #[test]
    fn a_uniform_length_fastq_still_reports_its_own_length() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("uniform.fastq");
        make_test_fastq(&fastq_path, 50, 75);
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        assert_eq!(model.quality_score_model().assumed_read_length, 75);
    }

    /// Write `n_a` records of `len_a` with quality `qual_a`, then `n_b` of `len_b` / `qual_b`.
    /// Quality strings are supplied whole so the expected model is readable from the fixture.
    fn write_two_phase_fastq(path: &PathBuf, n_a: usize, qual_a: &str, n_b: usize, qual_b: &str) {
        use std::io::Write;
        let mut out = String::new();
        for (n, qual) in [(n_a, qual_a), (n_b, qual_b)] {
            for i in 0..n {
                out.push_str(&format!(
                    "@r{}_{}\n{}\n+\n{}\n",
                    qual.len(),
                    i,
                    "A".repeat(qual.len()),
                    qual
                ));
            }
        }
        std::fs::File::create(path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();
    }

    /// Every position of a model must emit what was planted there, across seeds. A position
    /// whose transition row went untrained falls back to a uniform draw over the option set,
    /// so with two planted values it lands on the wrong one about half the time.
    fn assert_positions_match(
        q: &eidolon_core::models::quality_scores::QualityScoreModel,
        expected: &str,
        qual_offset: usize,
    ) {
        let want: Vec<usize> = expected.bytes().map(|b| b as usize - qual_offset).collect();
        for seed in ["a", "b", "c", "d"] {
            let mut rng =
                eidolon_core::rng::NeatRng::new_from_seed(&vec![seed.to_string()]).unwrap();
            let got = q.generate_quality_scores(want.len(), &mut rng).unwrap();
            assert_eq!(got.len(), want.len());
            for (pos, (g, w)) in got.iter().zip(want.iter()).enumerate() {
                assert_eq!(
                    g,
                    w,
                    "seed {seed}, position {} emitted Q{g} but Q{w} was planted there; an \
                     untrained row falls back to a uniform draw",
                    pos + 1
                );
            }
        }
    }

    /// THE regression for #700. An ordinary 60-column FASTA used to produce a model at exit 0:
    /// its sequence lines landed in the quality accumulator, and because every nucleotide
    /// decodes to a plausible Phred score (A->Q32, C->Q34, G->Q38, T->Q51) the output was
    /// indistinguishable from a real model - `error_rate` 0.000299, quality options
    /// [32, 34, 38, 51]. Measured, on exactly this fixture.
    #[test]
    fn a_fasta_is_rejected_and_says_so() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();
        let fasta_path = temp.path().join("genome.fasta");
        let mut out = String::from(">chr1 test contig\n");
        for _ in 0..40 {
            out.push_str(&format!("{}\n", "ACGT".repeat(15)));
        }
        std::fs::File::create(&fasta_path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();

        let output_path = temp.path().join("model.json.gz");
        let err = runner(&make_config(fasta_path, output_path.clone())).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("FASTA"),
            "the common mistake deserves its own message, got: {msg}"
        );
        assert!(msg.contains("line 1"), "error must name the line: {msg}");
        assert!(
            !output_path.exists(),
            "a rejected file must not leave a model behind - writing one is the whole defect"
        );
    }

    /// A corrupt record in the MIDDLE of an otherwise good file fails at that record, naming
    /// its line. Failing at record 1 only would let damage anywhere else through.
    #[test]
    fn a_bad_header_mid_file_names_its_line() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("corrupt.fastq");
        let mut out = String::new();
        for i in 0..5 {
            out.push_str(&format!("@r{i}\nACGT\n+\nIIII\n"));
        }
        // Record 6 has lost its header. Its first line is at 6*4+1 = 21.
        out.push_str("r5_no_at\nACGT\n+\nIIII\n");
        std::fs::File::create(&fastq_path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();

        let output_path = temp.path().join("model.json.gz");
        let err = runner(&make_config(fastq_path, output_path.clone())).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("line 21"),
            "error must name line 21, got: {msg}"
        );
        assert!(
            msg.contains('@'),
            "error should say what was expected: {msg}"
        );
        assert!(
            !output_path.exists(),
            "a corrupt file must not yield a model"
        );
    }

    /// The '+' separator and the sequence/quality length agreement each get their own check,
    /// so each gets its own test.
    #[test]
    fn a_bad_separator_and_a_length_mismatch_are_both_caught() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();

        let bad_plus = temp.path().join("badplus.fastq");
        std::fs::File::create(&bad_plus)
            .unwrap()
            .write_all(b"@r0\nACGT\nNOT_A_PLUS\nIIII\n")
            .unwrap();
        let out0 = temp.path().join("m0.json.gz");
        let msg = runner(&make_config(bad_plus, out0))
            .unwrap_err()
            .to_string();
        assert!(
            msg.contains("line 3"),
            "separator error names line 3: {msg}"
        );
        assert!(
            msg.contains('+'),
            "separator error says what was expected: {msg}"
        );

        let mismatch = temp.path().join("mismatch.fastq");
        std::fs::File::create(&mismatch)
            .unwrap()
            .write_all(b"@r0\nACGTACGT\n+\nIIII\n")
            .unwrap();
        let out1 = temp.path().join("m1.json.gz");
        let msg = runner(&make_config(mismatch, out1))
            .unwrap_err()
            .to_string();
        assert!(
            msg.contains("line 4"),
            "length error names the quality line: {msg}"
        );
        assert!(
            msg.contains('8') && msg.contains('4'),
            "length error reports both lengths: {msg}"
        );
    }

    /// A truncated final record is corruption, not a clean end of file.
    #[test]
    fn a_truncated_final_record_is_an_error() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("truncated.fastq");
        // Two good records, then a header and sequence with no '+' or quality line.
        std::fs::File::create(&fastq_path)
            .unwrap()
            .write_all(b"@r0\nACGT\n+\nIIII\n@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
            .unwrap();
        let output_path = temp.path().join("model.json.gz");
        let msg = runner(&make_config(fastq_path, output_path))
            .unwrap_err()
            .to_string();
        assert!(
            msg.contains("line 9") && msg.contains("truncated"),
            "error must say the record starting at line 9 is truncated: {msg}"
        );
    }

    /// MUST NOT FIRE, and this is the case that breaks naive FASTQ parsers: a quality line may
    /// legitimately begin with '@', because '@' is Q31 under Phred+33. A parser that scanned for
    /// '@' to find record boundaries would mis-frame this file. Reading fixed four-line records
    /// does not, and this pins that it stays true.
    #[test]
    fn a_quality_line_starting_with_an_at_sign_is_not_a_header() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("at_qual.fastq");
        let mut out = String::new();
        for i in 0..50 {
            // '@' = Q31, 'I' = Q40. First base Q31, the rest Q40.
            out.push_str(&format!("@r{i}\nACGTACGTAC\n+\n@IIIIIIIII\n"));
        }
        std::fs::File::create(&fastq_path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();

        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(
            q.assumed_read_length, 10,
            "all 50 records must have been read"
        );
        // The framing is right only if the '@' line was read as QUALITY: '@' is Q31 under
        // Phred+33, so Q31 must appear in the learned option set. Had the reader mistaken that
        // line for a record header, the file would have mis-framed and Q31 would be absent.
        assert!(
            q.quality_score_options.contains(&31),
            "the '@' line must have been read as quality (Q31), got options {:?}",
            q.quality_score_options
        );
        assert!(
            q.quality_score_options.contains(&40),
            "the 'I' bases are Q40: {:?}",
            q.quality_score_options
        );
        // NOT asserted here: that the model EMITS Q31 at position 1. It deliberately does not --
        // a quality line must not begin with '@', so position 1 never carries Q31 by design.
        // The rewrite that enforces that has a bug (#705): it decrements to 30 without checking
        // 30 is in the option set, which panics on a sparse set like this one. That is a
        // generation defect, not a parsing one, so it is out of scope here.
    }

    /// Write one record whose quality line is `len` characters.
    fn write_single_read_fastq(path: &PathBuf, len: usize) {
        use std::io::Write;
        let record = format!("@long\n{}\n+\n{}\n", "A".repeat(len), "I".repeat(len));
        std::fs::File::create(path)
            .unwrap()
            .write_all(record.as_bytes())
            .unwrap();
    }

    /// A read above the ceiling STOPS the run. The limit exists to bound the allocation
    /// without reintroducing silent truncation, so the failure has to be loud.
    #[test]
    fn a_read_above_the_ceiling_is_an_error_not_a_truncation() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("toolong.fastq");
        write_single_read_fastq(&fastq_path, 1500);
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_model_read_length = 1000;

        let err = runner(&config).unwrap_err();
        let msg = err.to_string();
        // Assert the CONTENT: a reader has to learn the length, the limit, and the way out.
        assert!(
            msg.contains("1500"),
            "error must name the read length: {msg}"
        );
        assert!(msg.contains("1000"), "error must name the limit: {msg}");
        assert!(
            msg.contains("max_model_read_length"),
            "error must name the key to raise: {msg}"
        );
        assert!(
            !output_path.exists(),
            "a refused run must not leave a model file behind"
        );
    }

    /// Must not fire: a read exactly AT the ceiling is fitted. An off-by-one here would reject
    /// legitimate data, which is the failure a limit most easily introduces.
    #[test]
    fn a_read_exactly_at_the_ceiling_is_accepted() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("atlimit.fastq");
        write_single_read_fastq(&fastq_path, 300);
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_model_read_length = 300;

        runner(&config).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        assert_eq!(model.quality_score_model().assumed_read_length, 300);
    }

    /// Raising the key admits the read, and 0 disables the limit — the escape hatch the error
    /// message points at has to work, or the message is a dead end.
    #[test]
    fn the_ceiling_can_be_raised_or_disabled() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("long.fastq");
        write_single_read_fastq(&fastq_path, 1500);

        for (limit, tag) in [(2000usize, "raised"), (0usize, "disabled")] {
            let output_path = temp.path().join(format!("model_{tag}.json.gz"));
            let mut config = make_config(fastq_path.clone(), output_path.clone());
            config.max_model_read_length = limit;
            runner(&config).unwrap_or_else(|e| panic!("limit {limit} should admit 1500 bp: {e}"));
            let model = SequencingErrorModel::from_file(&output_path).unwrap();
            assert_eq!(model.quality_score_model().assumed_read_length, 1500);
        }
    }

    /// THE criterion for the #694 fitter: the fraction it reports is the fraction that is
    /// there.
    ///
    /// Known answer by construction. 30 of 100 reads carry a tail of Q10 over their last 50
    /// bases and 70 carry Q40, so the classifier's answer is 0.30 whatever the code computes.
    #[test]
    fn the_fitter_recovers_the_planted_degraded_fraction() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("bimodal.fastq");
        let healthy = "I".repeat(100);
        let degraded = "I".repeat(50) + &"+".repeat(50);
        write_two_phase_fastq(&fastq_path, 70, &healthy, 30, &degraded);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.fit_quality_degradation = true;
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        let d = q
            .degradation
            .as_ref()
            .expect("a two-population fit must produce a degraded population");
        assert!(
            (0.28..0.32).contains(&d.read_fraction),
            "planted 30 of 100 reads degraded; fitted read_fraction {:.4}",
            d.read_fraction
        );
        // The two populations share one option set, so the degraded tensor must cover the same
        // positions as the healthy one. A ragged pair indexes the wrong scores at generation.
        assert_eq!(
            d.degraded_distros.len(),
            q.distros_from_one.len(),
            "both populations describe the same read length"
        );
    }

    /// Rule 4, pinned: the reported fraction is over the reads that could be CLASSIFIED, not
    /// over every read seen.
    ///
    /// Every other fixture here has reads longer than the window, so the two denominators are
    /// identical and neither is pinned. Found by mutation: dividing by records seen instead of
    /// reads classified passed every test above. Here 20 of 120 reads are too short to classify,
    /// so the two answers are 0.30 and 0.25 and only one of them is right.
    #[test]
    fn reads_too_short_to_classify_are_out_of_the_denominator() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("ragged.fastq");
        let mut out = String::new();
        for i in 0..70 {
            out.push_str(&format!(
                "@h{i}\n{}\n+\n{}\n",
                "A".repeat(100),
                "I".repeat(100)
            ));
        }
        for i in 0..30 {
            let q = "I".repeat(50) + &"+".repeat(50);
            out.push_str(&format!("@d{i}\n{}\n+\n{q}\n", "A".repeat(100)));
        }
        // Shorter than the 50 base window, so unclassifiable in either direction.
        for i in 0..20 {
            out.push_str(&format!(
                "@s{i}\n{}\n+\n{}\n",
                "A".repeat(30),
                "I".repeat(30)
            ));
        }
        std::fs::File::create(&fastq_path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.fit_quality_degradation = true;
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let d = model.quality_score_model().degradation.clone().unwrap();
        assert!(
            (0.28..0.32).contains(&d.read_fraction),
            "30 degraded of 100 CLASSIFIABLE reads is 0.30; over all 120 records it would read \
             0.25. Got {:.4}",
            d.read_fraction
        );
    }

    /// Write a FASTQ from `(count, quality-string)` groups. The quality string sets the read
    /// length, so this is how a fixture gives the two populations DIFFERENT lengths.
    fn write_groups(path: &PathBuf, groups: &[(usize, String)]) {
        use std::io::Write;
        let mut out = String::new();
        for (gi, (n, qual)) in groups.iter().enumerate() {
            for i in 0..*n {
                out.push_str(&format!(
                    "@g{gi}_r{i}\n{}\n+\n{qual}\n",
                    "A".repeat(qual.len())
                ));
            }
        }
        std::fs::File::create(path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();
    }

    /// A population that cannot cover the model's read length must be REFUSED, not padded.
    ///
    /// `every_modeled_position_is_trained` pins the invariant for a single population: a
    /// position row exists only because a read reached it and that read trained it. Two
    /// populations break it -- the shorter one's tensor has to be padded to the common length,
    /// and `build_distros` turns an all-zero row into a UNIFORM distribution. So 60 bp degraded
    /// reads against 100 bp healthy ones would give the degraded population a fabricated,
    /// uniform quality profile from base 62 on, indistinguishable downstream from a fit.
    ///
    /// The short reads here are 60 bp, longer than the 50 base classification window, so they
    /// ARE classified. That is what separates this from
    /// `reads_too_short_to_classify_are_out_of_the_denominator`, where the short reads are
    /// unclassifiable and never reach the degraded tensor at all.
    #[test]
    fn a_degraded_population_too_short_to_cover_the_model_is_an_error() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("short_degraded.fastq");
        write_groups(
            &fastq_path,
            &[
                (70, "I".repeat(100)),
                (30, "I".repeat(10) + &"+".repeat(50)),
            ],
        );

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path);
        config.fit_quality_degradation = true;
        let err = runner(&config).expect_err(
            "a degraded population covering 60 of 100 positions must be refused, not padded \
             with uniform rows",
        );
        let msg = err.to_string();
        for needle in ["degraded", "60", "100"] {
            assert!(
                msg.contains(needle),
                "the error must name the population and both lengths; missing {needle:?} in: {msg}"
            );
        }
    }

    /// The mirror: the HEALTHY population is the short one. Symmetric by construction -- either
    /// tensor can be the one that gets padded -- and an asymmetric guard would fabricate the
    /// healthy profile instead of the degraded one.
    #[test]
    fn a_healthy_population_too_short_to_cover_the_model_is_an_error() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("short_healthy.fastq");
        write_groups(
            &fastq_path,
            &[(70, "I".repeat(50) + &"+".repeat(50)), (30, "I".repeat(60))],
        );

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path);
        config.fit_quality_degradation = true;
        let err = runner(&config)
            .expect_err("a healthy population covering 60 of 100 positions must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("healthy") && msg.contains("60") && msg.contains("100"),
            "the error must name the population and both lengths: {msg}"
        );
    }

    /// #723: the same refusal on the MATE axis. R2 covering 60 of the model's 100 positions
    /// must be refused, not padded -- `build_distros` turns a padded all-zero row into a
    /// UNIFORM draw, so R2 would carry invented quality from base 62 on, and nothing
    /// downstream could tell that from a fit.
    ///
    /// Reachable only through this path. `with_mate_r2`'s own length guard cannot fire from
    /// here, because the fitter hands both mates one `read_length` by construction; the unit
    /// tests in `sequencing_error_model.rs` cover that guard for a library caller.
    #[test]
    fn a_short_r2_is_refused_rather_than_padded() {
        let temp = tempfile::tempdir().unwrap();
        let r1_path = temp.path().join("long_r1.fastq");
        let r2_path = temp.path().join("short_r2.fastq");
        write_groups(&r1_path, &[(20, "I".repeat(100))]);
        write_groups(&r2_path, &[(20, "I".repeat(60))]);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(r1_path, output_path.clone());
        config.fastq_file_r2 = Some(r2_path);
        let err = runner(&config).expect_err(
            "an R2 covering 60 of the model's 100 positions must be refused, not padded with \
             uniform rows",
        );
        let msg = err.to_string();
        for needle in ["R2", "60", "100"] {
            assert!(
                msg.contains(needle),
                "the error must name the mate and both lengths; missing {needle:?} in: {msg}"
            );
        }
        // The remedy must name the axis this run actually used. One shared message told a run
        // that never set `fit_quality_degradation` to unset it.
        assert!(
            msg.contains("fastq_file_r2") && !msg.contains("fit_quality_degradation"),
            "a mate-length refusal must point at the mate inputs, not at a flag this run \
             never set: {msg}"
        );
        assert!(
            !output_path.exists(),
            "a refused run must not leave a model file behind"
        );
    }

    /// The mirror. Either mate can be the short one, and an asymmetric guard would fabricate
    /// R1's profile instead of R2's.
    #[test]
    fn a_short_r1_is_refused_when_the_mate_is_longer() {
        let temp = tempfile::tempdir().unwrap();
        let r1_path = temp.path().join("short_r1.fastq");
        let r2_path = temp.path().join("long_r2.fastq");
        write_groups(&r1_path, &[(20, "I".repeat(60))]);
        write_groups(&r2_path, &[(20, "I".repeat(100))]);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(r1_path, output_path);
        config.fastq_file_r2 = Some(r2_path);
        let err = runner(&config)
            .expect_err("an R1 covering 60 of the model's 100 positions must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("R1") && msg.contains("60") && msg.contains("100"),
            "the error must name the mate and both lengths: {msg}"
        );
    }

    /// MUST NOT FIRE: two mates that reach the same length are accepted, and the model carries
    /// both. The guard is about a population that cannot cover the model, not about the mates
    /// differing -- which is the entire point of fitting them separately.
    #[test]
    fn mates_that_both_reach_the_model_length_are_accepted() {
        let temp = tempfile::tempdir().unwrap();
        let r1_path = temp.path().join("r1.fastq");
        let r2_path = temp.path().join("r2.fastq");
        write_groups(&r1_path, &[(20, "I".repeat(100))]);
        write_groups(&r2_path, &[(20, "5".repeat(100))]);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(r1_path, output_path.clone());
        config.fastq_file_r2 = Some(r2_path);
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let r2 = model
            .quality_score_model_r2()
            .expect("a two-FASTQ fit must carry an R2 population");
        assert_eq!(model.quality_score_model().assumed_read_length, 100);
        assert_eq!(r2.assumed_read_length, 100);
    }

    /// The fourth population. A fully loaded fit carries R1/R2 x healthy/degraded, and
    /// `R2 degraded` is the one label that sits on both axes -- so its refusal has to name
    /// both, and a run that trimmed only its R2 mate needs to hear about the mate files as
    /// well as the flag.
    ///
    /// This drives the four-tensor path far enough to pin the refusal. It is NOT a fidelity
    /// test of a four-population model: that a loaded fit reproduces all four distributions is
    /// still unmeasured.
    #[test]
    fn a_short_r2_degraded_population_is_refused_and_names_both_axes() {
        let temp = tempfile::tempdir().unwrap();
        let r1_path = temp.path().join("both_ok_r1.fastq");
        let r2_path = temp.path().join("short_deg_r2.fastq");
        // R1: healthy and degraded both reach 100 bp.
        write_groups(
            &r1_path,
            &[
                (40, "I".repeat(100)),
                (30, "I".repeat(50) + &"+".repeat(50)),
            ],
        );
        // R2: healthy reaches 100 bp, degraded stops at 60.
        write_groups(
            &r2_path,
            &[
                (40, "I".repeat(100)),
                (30, "I".repeat(10) + &"+".repeat(50)),
            ],
        );

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(r1_path, output_path.clone());
        config.fastq_file_r2 = Some(r2_path);
        config.fit_quality_degradation = true;
        let err = runner(&config)
            .expect_err("R2's degraded population covers 60 of 100 positions and must be refused");
        let msg = err.to_string();
        for needle in ["R2 degraded", "60", "100"] {
            assert!(
                msg.contains(needle),
                "the error must name the population and both lengths; missing {needle:?} in: \
                 {msg}"
            );
        }
        assert!(
            msg.contains("mate files") && msg.contains("fit_quality_degradation"),
            "this label sits on both axes, so its remedy must name both: {msg}"
        );
        assert!(
            !output_path.exists(),
            "a refused run must not leave a model file behind"
        );
    }

    /// `max_reads` applies PER MATE (#723). A shared budget spends it all on R1 and fits R2
    /// from what is left -- nothing -- which the ragged refusal above would turn into an error
    /// rather than a silent lopsided fit, but the declared semantics are that each mate gets
    /// its own N.
    ///
    /// KNOWN ANSWER: every R1 record is Q40 and every R2 record Q20, 100 of each, capped at 10
    /// per mate. Whichever records the sampler keeps, R1 must be fitted all-Q40 and R2 must
    /// exist and be fitted all-Q20. A shared budget starves R2.
    #[test]
    fn max_reads_applies_to_each_mate() {
        let temp = tempfile::tempdir().unwrap();
        let r1_path = temp.path().join("capped_r1.fastq");
        let r2_path = temp.path().join("capped_r2.fastq");
        // 'I' = Q40, '5' = Q20 under Phred+33.
        write_groups(&r1_path, &[(100, "I".repeat(50))]);
        write_groups(&r2_path, &[(100, "5".repeat(50))]);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(r1_path, output_path.clone());
        config.fastq_file_r2 = Some(r2_path);
        config.max_reads = 10;
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q1 = model.quality_score_model();
        let q2 = model
            .quality_score_model_r2()
            .expect("R2 must be fitted from its own sample, not from R1's leftovers");
        assert_eq!(q1.quality_score_options, vec![20, 40]);

        let mut rng = eidolon_core::rng::NeatRng::new_from_seed(&vec!["723".to_string()]).unwrap();
        let s1 = q1.generate_quality_scores(50, &mut rng).unwrap();
        let s2 = q2.generate_quality_scores(50, &mut rng).unwrap();
        assert!(s1.iter().all(|&s| s == 40), "R1 is all-Q40: {s1:?}");
        assert!(s2.iter().all(|&s| s == 20), "R2 is all-Q20: {s2:?}");
    }

    /// #721 KNOWN ANSWER: the first 500 records are Q40 and the last 500 Q10, as a flowcell
    /// whose later tiles read worse would be. A cap of 100 must see both halves, so the
    /// option set is {Q10, Q40}. Taking the head of the file, as the cap did before, fits
    /// from the Q40 half only and gives {Q40}.
    #[test]
    fn a_capped_fit_samples_the_whole_file_not_its_head() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("two_halves.fastq");
        write_groups(&fastq_path, &[(500, "I".repeat(50)), (500, "+".repeat(50))]);
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 100;
        runner(&config).unwrap();

        let q = SequencingErrorModel::from_file(&output_path).unwrap();
        assert_eq!(
            q.quality_score_model().quality_score_options,
            vec![10, 40],
            "a capped fit saw only one half of the file: `max_reads` is taking the head (#721)"
        );
    }

    /// The sampler keeps about `max_reads` records, spread evenly. KNOWN ANSWER: 1,000 of
    /// 100,000 is p = 0.01, so the count is Binomial(100000, 0.01): mean 1,000, sigma ~31.5,
    /// and each tenth of the file expects 100 (sigma ~10). Bounds are +/-5 sigma. A head
    /// sampler puts all 1,000 in the first tenth.
    #[test]
    fn the_sampler_keeps_about_max_reads_spread_across_the_file() {
        let mut sampler = ReadSampler::new(1_000, 100_000).unwrap();
        let mut per_tenth = [0usize; 10];
        for i in 0..100_000 {
            if sampler.keep().unwrap() {
                per_tenth[i / 10_000] += 1;
            }
        }
        let kept: usize = per_tenth.iter().sum();
        assert!(
            (843..=1_157).contains(&kept),
            "kept {kept} of 100,000 at p = 0.01"
        );
        for (t, &n) in per_tenth.iter().enumerate() {
            assert!(
                (50..=150).contains(&n),
                "tenth {t} kept {n}, expected ~100: {per_tenth:?}"
            );
        }
    }

    /// MUST NOT FIRE: a cap at or above the file's size, or no cap, keeps every record.
    #[test]
    fn the_sampler_keeps_everything_when_the_cap_does_not_bind() {
        for (max_reads, total) in [(0, 500), (500, 500), (1_000, 500)] {
            let mut sampler = ReadSampler::new(max_reads, total).unwrap();
            assert!(
                sampler.keep_probability.is_none(),
                "max_reads {max_reads}, total {total}"
            );
            assert!((0..total).all(|_| sampler.keep().unwrap()));
        }
    }

    /// MUST NOT FIRE, end to end: a cap that does not bind writes the same model as no cap.
    #[test]
    fn a_cap_above_the_file_size_fits_the_same_model_as_no_cap() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("small.fastq");
        write_groups(&fastq_path, &[(60, "I".repeat(50)), (40, "+".repeat(50))]);
        let fit = |name: &str, max_reads: usize| {
            let out = temp.path().join(name);
            let mut config = make_config(fastq_path.clone(), out.clone());
            config.max_reads = max_reads;
            runner(&config).unwrap();
            serde_json::to_value(SequencingErrorModel::from_file(&out).unwrap()).unwrap()
        };
        assert_eq!(fit("uncapped.json.gz", 0), fit("capped.json.gz", 100));
    }

    /// MUST NOT FIRE: ragged read lengths are fine as long as BOTH populations reach the
    /// model's read length. The guard is about a population that cannot cover the model, not
    /// about variable lengths, which every real library has.
    #[test]
    fn ragged_lengths_are_accepted_when_both_populations_reach_the_model_length() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("ragged_both.fastq");
        write_groups(
            &fastq_path,
            &[
                (40, "I".repeat(100)),
                (20, "I".repeat(60)),
                (30, "I".repeat(50) + &"+".repeat(50)),
                (10, "I".repeat(10) + &"+".repeat(50)),
            ],
        );

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.fit_quality_degradation = true;
        runner(&config).expect("both populations reach 100 bp, so this must fit");

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(q.assumed_read_length, 100);
        let d = q
            .degradation
            .as_ref()
            .expect("two populations were asked for");
        assert!(
            (0.38..0.42).contains(&d.read_fraction),
            "40 of 100 classified reads are degraded; got {:.4}",
            d.read_fraction
        );
    }

    /// MUST NOT FIRE: the same bimodal input with the fit switched off produces a
    /// single-population model, unchanged from every model before #694.
    #[test]
    fn without_the_flag_the_same_input_fits_one_population() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("bimodal.fastq");
        let healthy = "I".repeat(100);
        let degraded = "I".repeat(50) + &"+".repeat(50);
        write_two_phase_fastq(&fastq_path, 70, &healthy, 30, &degraded);

        let output_path = temp.path().join("model.json.gz");
        let config = make_config(fastq_path, output_path.clone());
        assert!(!config.fit_quality_degradation, "default must be off");
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        assert!(
            model.quality_score_model().degradation.is_none(),
            "opt-in means opt-in: an unflagged fit must not attach a degraded population"
        );
    }

    /// MUST NOT FIRE, and a hard failure rather than a warning: asked for two populations on a
    /// library that has one, the fit refuses instead of emitting a degraded tensor made
    /// entirely of uniform fallback.
    #[test]
    fn a_library_with_no_degraded_reads_is_an_error_not_an_empty_population() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("healthy.fastq");
        let healthy = "I".repeat(100);
        write_two_phase_fastq(&fastq_path, 100, &healthy, 0, &healthy);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.fit_quality_degradation = true;

        let err = runner(&config).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("no degraded population"),
            "the error must say what was not found: {msg}"
        );
        assert!(
            msg.contains("100"),
            "and over what denominator, so the reader can judge it: {msg}"
        );
        assert!(
            !output_path.exists(),
            "a refused fit must not leave a model file behind"
        );
    }

    /// The mirror of the test above. `n_degraded == 0` was a hard error from the start;
    /// `n_healthy == 0` was not, so a cut above the library's range classified every read as
    /// degraded and left the HEALTHY tensor empty, which falls back to uniform. That is the
    /// same garbage the end-to-end test was written to catch, on the other population.
    #[test]
    fn a_library_with_no_healthy_reads_is_an_error_not_an_empty_population() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("all_degraded.fastq");
        let healthy = "I".repeat(100); // Q40 throughout
        write_two_phase_fastq(&fastq_path, 100, &healthy, 0, &healthy);

        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.fit_quality_degradation = true;
        config.degradation_tail_cut = 45; // above Q40, so every read classifies as degraded

        let err = runner(&config).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("no healthy population"),
            "the error must say which population was not found: {msg}"
        );
        assert!(
            msg.contains("100"),
            "and over what denominator, so the reader can judge it: {msg}"
        );
        assert!(
            !output_path.exists(),
            "a refused fit must not leave a model file behind"
        );
    }

    /// The cut is a tuning knob, so moving it must move the fitted fraction. A classifier that
    /// ignored its threshold would pass every test above.
    #[test]
    fn the_tail_cut_actually_selects_the_population() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("mid.fastq");
        // Tails at Q30: above a Q25 cut, below a Q35 cut. Nothing else distinguishes them.
        let healthy = "I".repeat(100);
        let midling = "I".repeat(50) + &"?".repeat(50);
        write_two_phase_fastq(&fastq_path, 60, &healthy, 40, &midling);

        // At Q25 the Q30 tails are healthy, so there is no degraded population at all.
        let out_low = temp.path().join("low.json.gz");
        let mut low = make_config(fastq_path.clone(), out_low);
        low.fit_quality_degradation = true;
        assert!(
            runner(&low).is_err(),
            "at a Q25 cut nothing is degraded, which is an error, not a 0% population"
        );

        // At Q35 they are degraded, and the fraction is the planted 40 of 100.
        let out_high = temp.path().join("high.json.gz");
        let mut high = make_config(fastq_path, out_high.clone());
        high.fit_quality_degradation = true;
        high.degradation_tail_cut = 35;
        runner(&high).unwrap();
        let model = SequencingErrorModel::from_file(&out_high).unwrap();
        let d = model.quality_score_model().degradation.clone().unwrap();
        assert!(
            (0.38..0.42).contains(&d.read_fraction),
            "at a Q35 cut the planted 40 of 100 are degraded; got {:.4}",
            d.read_fraction
        );
    }

    /// THE regression for the first review finding on #698: a read LONGER than everything
    /// before it, past where the old 1000-record sample stopped looking.
    ///
    /// Known answer: the long reads are Q40 for 50 positions then Q10 for 50, so a model that
    /// saw them must advertise 100 bp and emit Q10 across the tail. Under the sampling code
    /// `assumed_read_length` was 50 and positions 51-100 were dropped entirely — measured, on
    /// this exact shape, before the fix.
    #[test]
    fn a_long_record_past_the_old_sampling_boundary_is_fitted() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("late_long.fastq");
        let short_qual = "I".repeat(50);
        let long_qual = "I".repeat(50) + &"+".repeat(50);
        // 1000 short records is exactly the old sample size, so the long ones begin one record
        // past where it stopped reading.
        write_two_phase_fastq(&fastq_path, 1000, &short_qual, 200, &long_qual);
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(
            q.assumed_read_length, 100,
            "the longest read is 100 bp; the model must not stop at the 1000-record boundary"
        );
        assert_eq!(q.distros_from_one.len(), 99);
        // Positions 1-50 are Q40 in BOTH populations and 51-100 are Q10 in the only population
        // that reaches them, so the whole model is deterministic.
        assert_positions_match(q, &long_qual, 33);
    }

    /// THE regression for the second review finding on #698: `max_reads` must not let a
    /// record it excluded set the model's read length. Under the earlier sampling code the
    /// model advertised 100 bp from records the cap had skipped, and positions 51-100 were
    /// uniform noise between Q10 and Q40.
    ///
    /// Known answer: 110 records of 50 bp plus one of 100 bp, capped at 10. The sampler is
    /// seeded, so the record numbers it skips are fixed; the 100 bp record is placed on the
    /// first of them. Only 50 bp records are fitted, so the model is a 50 bp model.
    #[test]
    fn max_reads_does_not_advertise_positions_it_did_not_fit() {
        const TOTAL: usize = 111;
        let mut sampler = ReadSampler::new(10, TOTAL).unwrap();
        let kept: Vec<bool> = (0..TOTAL).map(|_| sampler.keep().unwrap()).collect();
        let skipped = kept
            .iter()
            .position(|k| !k)
            .expect("p = 10/111 skips some record");
        assert!(
            kept.iter().any(|&k| k),
            "the fixture needs at least one fitted record"
        );

        let short_qual = "I".repeat(25) + &"+".repeat(25);
        let long_qual = "I".repeat(50) + &"+".repeat(50);
        let mut out = String::new();
        for i in 0..TOTAL {
            let q = if i == skipped {
                &long_qual
            } else {
                &short_qual
            };
            out.push_str(&format!("@r{i}\n{}\n+\n{q}\n", "A".repeat(q.len())));
        }
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("capped.fastq");
        std::fs::write(&fastq_path, out).unwrap();
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 10;
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(
            q.assumed_read_length, 50,
            "the only 100 bp record was skipped by the cap, so the model describes 50 positions"
        );
        assert_eq!(q.distros_from_one.len(), 49);
        assert_positions_match(q, &short_qual, 33);
    }

    /// The invariant the growth approach buys: a position row exists only because some read
    /// reached it, and that read necessarily incremented it, so no modeled position can be
    /// untrained. Three lengths interleaved, quality alternating by position, so every one of
    /// the 90 positions has exactly one correct answer and a uniform fallback is visible.
    #[test]
    fn every_modeled_position_is_trained() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("ragged.fastq");
        let full: String = (0..90)
            .map(|i| if i % 2 == 0 { 'I' } else { '+' })
            .collect();
        use std::io::Write;
        let mut out = String::new();
        for (i, len) in [30usize, 60, 90].iter().cycle().take(60).enumerate() {
            let qual = &full[..*len];
            out.push_str(&format!("@r{}\n{}\n+\n{}\n", i, "A".repeat(*len), qual));
        }
        std::fs::File::create(&fastq_path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();

        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(q.assumed_read_length, 90);
        assert_positions_match(q, &full, 33);
    }

    /// A short record that is NOT first must not shrink the model either — the maximum over
    /// the sample is what is wanted, not the minimum or the last seen.
    #[test]
    fn a_short_record_later_in_the_file_does_not_shrink_the_model() {
        use std::io::Write;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("late_short.fastq");
        let seq = |n: usize| -> String { "ACGT".chars().cycle().take(n).collect() };
        let mut out = String::new();
        for i in 0..40 {
            let n = if i == 17 { 12 } else { 60 };
            out.push_str(&format!("@r{}\n{}\n+\n{}\n", i, seq(n), "I".repeat(n)));
        }
        std::fs::File::create(&fastq_path)
            .unwrap()
            .write_all(out.as_bytes())
            .unwrap();
        let output_path = temp.path().join("model.json.gz");
        runner(&make_config(fastq_path, output_path.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        assert_eq!(model.quality_score_model().assumed_read_length, 60);
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
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: None,
            transition_matrix_file: Some(tsv_path),
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
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
    /// 0.0/0.1695/0.6878/0.1427), while the sequencing-error default is the matrix in the
    /// shipped model file (A row 0.0/0.292/0.154/0.553, NovaSeq mean, #779). Deriving from
    /// the wrong one makes this test fail against correct code, which is how the distinction
    /// was found. `from_raw_data`'s fallback reads the same file, so there is one copy.
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
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
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
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
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

    /// #752: a sample's own variants mismatch the reference in every read carrying them, so
    /// they must not be fitted as sequencing errors. Same fixture as above (ref AAAACCCC, read
    /// CCGTAAAG at chr1:1). Masking positions 1-2 removes both A->C, leaving the A row at
    /// G 0.5 / T 0.5; the C row is untouched.
    /// #779: an overlap fit needs mates that overlap. These records are unpaired, so there is
    /// no overlap evidence at all, and the fit must refuse and say how to proceed rather than
    /// fit a matrix from nothing.
    #[test]
    fn an_overlap_fit_without_enough_evidence_is_refused() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("unpaired.bam");
        write_test_bam(&bam_path, 8, b"AAAACCCC", b"CCGTAAAG", true);
        let mut config = make_config(fastq_path, temp.path().join("model.json.gz"));
        config.bam_file = Some(bam_path);
        config.bam_method = BamMethod::Overlap;
        let err = runner(&config).expect_err("no overlap evidence must be refused");
        let msg = format!("{err}");
        assert!(
            msg.contains("bam_method: mismatch"),
            "the refusal must name the alternative: {msg}"
        );
    }

    /// KNOWN ANSWER, by hand, every bin populated as in real reads. Bins 0-1: 1,000,000 observed
    /// bases each, A->C 10,000 and A->G 30,000 each. Bins 2-4: 100,000,000 bases each, A->C
    /// 900,000 and A->G 100,000 each. Pooled counts give A->C 2.72M / 3.08M = 0.883. Per-base
    /// rates summed over the five bins give A->C 2(0.010) + 3(0.009) = 0.047 and A->G
    /// 2(0.030) + 3(0.001) = 0.063, so the whole read's A->C is 0.047 / 0.110. Every bin clears
    /// MIN_BIN_ERRORS, so none is merged.
    #[test]
    fn the_overlap_matrix_weights_every_part_of_the_read_equally() {
        let mut c = OverlapCounts::default();
        c.bin_bases = [1_000_000, 1_000_000, 100_000_000, 100_000_000, 100_000_000];
        for bin in 0..2 {
            c.by_bin[bin][0][1] = 10_000;
            c.by_bin[bin][0][2] = 30_000;
        }
        for bin in 2..5 {
            c.by_bin[bin][0][1] = 900_000;
            c.by_bin[bin][0][2] = 100_000;
        }
        let w = per_cycle_spectrum(&c);
        let share = w[0][1] / (w[0][1] + w[0][2]);
        assert!((share - 0.047 / 0.110).abs() < 1e-12, "A->C {share}");
    }

    /// A bin with too few errors is merged with its neighbor before reweighting, or its
    /// sampling noise counts as much as a full bin's (the round trip's first bin held 99 errors,
    /// about 8 per cell, and reweighting unmerged moved its worst cell from 0.007 to 0.024).
    ///
    /// KNOWN ANSWER. Bin 0: 1,000 bases, 10 errors, all A->C (a rate like its neighbor's, but a
    /// spectrum that is pure noise). Bin 1: 1,000,000 bases, 10,000 A->C and 10,000 A->G.
    /// Unmerged, the per-base rates give A->C 0.010 + 0.010 against A->G 0 + 0.010: A->C 0.667.
    /// Merged into one group (10,010 against 10,000 over 1,001,000 bases), A->C is 0.50025.
    #[test]
    fn a_sparse_bin_is_merged_before_reweighting() {
        let mut c = OverlapCounts::default();
        c.bin_bases = [1_000, 1_000_000, 0, 0, 0];
        c.by_bin[0][0][1] = 10;
        c.by_bin[1][0][1] = 10_000;
        c.by_bin[1][0][2] = 10_000;
        let w = per_cycle_spectrum(&c);
        let share = w[0][1] / (w[0][1] + w[0][2]);
        assert!((share - 10_010.0 / 20_010.0).abs() < 1e-9, "A->C {share}");
    }

    /// A merged group stands for every bin it spans. KNOWN ANSWER: bins 0 and 1 (1M bases each)
    /// hold 5,000 A->C and 5,000 A->G, each under MIN_BIN_ERRORS, so they merge into one group
    /// spanning 2 bins with per-base rates 0.0025 / 0.0025. Bins 2-4 (100M bases each) hold
    /// 900,000 A->C and 100,000 A->G. Weighted by span, A->C = 2(0.0025) + 3(0.009) = 0.032 and
    /// A->G = 2(0.0025) + 3(0.001) = 0.008: 0.8. Ignoring the span gives 0.843.
    #[test]
    fn a_merged_group_is_weighted_by_the_bins_it_spans() {
        let mut c = OverlapCounts::default();
        c.bin_bases = [1_000_000, 1_000_000, 100_000_000, 100_000_000, 100_000_000];
        c.by_bin[0][0][1] = 5_000;
        c.by_bin[1][0][2] = 5_000;
        for bin in 2..5 {
            c.by_bin[bin][0][1] = 900_000;
            c.by_bin[bin][0][2] = 100_000;
        }
        let w = per_cycle_spectrum(&c);
        let share = w[0][1] / (w[0][1] + w[0][2]);
        assert!((share - 0.8).abs() < 1e-12, "A->C {share}");
    }

    /// Must not fire: bins that share one spectrum give that spectrum, however unequal their
    /// observation counts.
    #[test]
    fn a_spectrum_that_does_not_change_along_the_read_is_left_alone() {
        let mut c = OverlapCounts::default();
        c.bin_bases = [1_000, 0, 0, 0, 100_000];
        c.by_bin[0][0][1] = 3;
        c.by_bin[0][0][2] = 1;
        c.by_bin[4][0][1] = 300;
        c.by_bin[4][0][2] = 100;
        let w = per_cycle_spectrum(&c);
        assert!((w[0][1] / (w[0][1] + w[0][2]) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn known_variant_positions_are_not_fitted_as_errors() {
        use eidolon_core::structs::nucleotides::Nucleotide;
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("mixed.bam");
        write_test_bam(&bam_path, 8, b"AAAACCCC", b"CCGTAAAG", true);
        let vcf = temp.path().join("known.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n\
             chr1\t1\t.\tA\tC\nchr1\t2\t.\tA\tC\n",
        )
        .unwrap();
        let output_path = temp.path().join("model.json.gz");
        let config = RunConfiguration {
            fastq_file: fastq_path,
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
            known_variants_vcf: Some(vcf),
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
        };
        runner(&config).unwrap();
        let tm = SequencingErrorModel::from_file(&output_path)
            .unwrap()
            .transition_distros()
            .clone();
        assert_row_cdf_eq(
            &row_cdf(&tm, Nucleotide::A),
            &[0.0, 0.0, 0.5, 1.0],
            "A row with both A->C masked (1 G, 1 T)",
        );
        assert_row_cdf_eq(
            &row_cdf(&tm, Nucleotide::C),
            &[0.75, 0.75, 1.0, 1.0],
            "C row, which the mask does not touch",
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
    fn test_runner_no_bam_file_uses_the_default_matrix() {
        // MUST-NOT-FIRE: with no `bam_file:` at all, nothing was requested and nothing errors --
        // the model carries the default matrix. This is the supported way to ask for it.
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let output_path = temp.path().join("model.json.gz");

        let config = make_config(fastq_path, output_path.clone());
        assert!(
            config.bam_file.is_none(),
            "this test's premise is no bam_file"
        );
        runner(&config).unwrap();

        // The case where inference must NOT fire: nothing was requested,
        // so the model has to carry the default matrix. Asserting only "a file appeared"
        // could not tell that apart from silently inventing a matrix from zero evidence.
        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        assert_row_cdf_eq(
            &row_cdf(
                model.transition_distros(),
                eidolon_core::structs::nucleotides::Nucleotide::A,
            ),
            &default_a_row_cdf(),
            "A row with no bam_file at all",
        );
    }

    /// MUST-FAIL: `bam_file:` that yields no mismatch evidence is an error, not a warning.
    ///
    /// The defect this pins (#529): the runner used to `warn!` and substitute the built-in
    /// default, producing a model that is byte-identical to an untrained one while the config
    /// says it was fitted from a BAM. Nothing downstream can distinguish those, which is the
    /// same shape as every other quiet failure in this repo — a value produced by a path that
    /// never ran.
    #[test]
    fn test_runner_bam_no_md_tags_is_an_error() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("nomd.bam");
        write_test_bam(&bam_path, 5, b"ACGT", b"TGCA", false);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
        };
        let err = runner(&config).expect_err("an MD-less bam_file must not silently default");
        let msg = err.to_string();

        // Assert the message names BOTH remedies. A bare "cannot infer" would leave the user
        // with no route forward, and the whole point of the flag is that it be discoverable
        // from the error rather than from the source.
        assert!(
            msg.contains("samtools calmd"),
            "error must name the MD remedy: {msg}"
        );
        assert!(
            msg.contains("remove `bam_file:`"),
            "error must name the other remedy -- dropping the key: {msg}"
        );

        // And it must fail BEFORE writing anything. A model on disk next to an error is worse
        // than either outcome alone, because a rerunning pipeline would pick it up.
        assert!(
            !output_path.exists(),
            "no model file may be written when inference fails"
        );
    }

    /// A BAM that HAS MD tags but records zero mismatches is the same failure with a different
    /// cause — a perfectly-matching alignment carries no transition evidence either. Included
    /// because the guard keys on the mismatch total, not on tag presence, and a reader could
    /// reasonably assume otherwise from the error text.
    fn write_fastq(path: &PathBuf, seqs: &[&str]) {
        let body: String = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(s.len())))
            .collect();
        std::fs::write(path, body).unwrap();
    }

    /// #767, path 2. With every read 1 bp the model has no position after the first, so the
    /// coverage check, which counted only those, could not fire. An R2 file whose records are
    /// all empty then reached the quality fit with no seed at all and failed as a bare
    /// `InvalidWeights([0.0])`. It must be refused as a population that covers 0 bp.
    #[test]
    fn an_r2_with_no_bases_is_refused_by_name_even_for_1bp_reads() {
        let temp = tempfile::tempdir().unwrap();
        let (r1, r2) = (temp.path().join("r1.fq"), temp.path().join("r2.fq"));
        write_fastq(&r1, &["A"; 20]);
        write_fastq(&r2, &[""; 20]);
        let mut config = make_config(r1, temp.path().join("m.json.gz"));
        config.fastq_file_r2 = Some(r2);
        let err = runner(&config)
            .expect_err("an R2 with no bases cannot be fitted")
            .to_string();
        assert!(
            err.contains("the R2 population covers 0 bp"),
            "the refusal must name the population and its coverage: {err}"
        );
    }

    /// MUST NOT FIRE: a pair of 1 bp files is a 1 bp model for each mate, not an error.
    #[test]
    fn a_pair_of_1bp_files_still_fits() {
        let temp = tempfile::tempdir().unwrap();
        let (r1, r2) = (temp.path().join("r1.fq"), temp.path().join("r2.fq"));
        write_fastq(&r1, &["A"; 20]);
        write_fastq(&r2, &["C"; 20]);
        let output = temp.path().join("m.json.gz");
        let mut config = make_config(r1, output.clone());
        config.fastq_file_r2 = Some(r2);
        runner(&config).unwrap();
        let model = SequencingErrorModel::from_file(&output).unwrap();
        assert_eq!(model.quality_score_model().assumed_read_length, 1);
        assert!(
            model.quality_score_model_r2().is_some(),
            "R2 gets its own model"
        );
    }

    fn mismatch_config(fastq: PathBuf, bam: PathBuf, output: PathBuf) -> RunConfiguration {
        RunConfiguration {
            fastq_file: fastq,
            fastq_file_r2: None,
            output_file: output,
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam),
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
        }
    }

    /// #767, path 1. An MD tag that disagrees with SEQ names a mismatch where the read base
    /// equals the reference base. Counted, that is a self-transition; the diagonal is then
    /// zeroed, and a row holding nothing else was left with no weight at all, which failed
    /// as a bare `InvalidWeights`. Here every MD mismatch is stale (MD says the reference at
    /// position 1 is A; SEQ has A there), so there is no evidence, and the error must say why.
    #[test]
    fn a_bam_whose_md_disagrees_with_seq_names_the_cause() {
        let temp = tempfile::tempdir().unwrap();
        let fastq = temp.path().join("t.fastq");
        make_test_fastq(&fastq, 20, 4);
        let bam = temp.path().join("stale_md.bam");
        write_test_bam_records(&bam, &vec![(&b"ACGT"[..], Some("0A3")); 5]);
        let config = mismatch_config(fastq, bam, temp.path().join("m.json.gz"));
        let err = runner(&config)
            .expect_err("stale MD tags are not evidence")
            .to_string();
        assert!(
            err.contains("MD") && err.contains("disagree"),
            "the error must name the MD/SEQ disagreement: {err}"
        );
        assert!(err.contains("5"), "and say how many positions: {err}");
    }

    /// Stale MD positions are left out, not counted: the genuine A->C mismatches fit the A
    /// row exactly, and the C row, which only stale positions named, falls back to uniform
    /// like any row with no evidence. Before the fix the C row's self-count failed the run.
    #[test]
    fn stale_md_positions_are_left_out_of_the_fit() {
        use eidolon_core::structs::nucleotides::Nucleotide;
        let temp = tempfile::tempdir().unwrap();
        let fastq = temp.path().join("t.fastq");
        make_test_fastq(&fastq, 20, 4);
        let bam = temp.path().join("mixed_md.bam");
        let mut records: Vec<(&[u8], Option<&str>)> = vec![(&b"CCGT"[..], Some("0A3")); 4];
        records.extend(vec![(&b"ACGT"[..], Some("1C2")); 3]);
        write_test_bam_records(&bam, &records);
        let output = temp.path().join("m.json.gz");
        runner(&mismatch_config(fastq, bam, output.clone())).unwrap();
        let model = SequencingErrorModel::from_file(&output).unwrap();
        let tm = model.transition_distros();
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::A),
            &[0.0, 1.0, 1.0, 1.0],
            "A row: 4 A->C",
        );
        assert_row_cdf_eq(
            &row_cdf(tm, Nucleotide::C),
            &[1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 1.0],
            "C row: only stale MD positions, so no evidence",
        );
    }

    #[test]
    fn test_runner_bam_with_md_but_no_mismatches_is_an_error() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let bam_path = temp.path().join("md_no_mismatch.bam");
        // ref == read, so the MD tag is written but encodes only matches.
        write_test_bam(&bam_path, 5, b"ACGT", b"ACGT", true);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            fastq_file_r2: None,
            output_file: output_path,
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
        };
        let err = runner(&config).expect_err("zero mismatches is no evidence, MD tags or not");
        assert!(
            err.to_string().contains("no read-vs-reference mismatches"),
            "{err}"
        );
    }

    /// End-to-end `runner` over a **real** aligned BAM — the gap #529 left open.
    ///
    /// Every other BAM in these tests is written by `write_test_bam` a few records at a time.
    /// This one is 1000 Genomes HG00096, GRCh37 `20:1000000-1050000`, aligned with bwa 0.5.9 in
    /// 2012 and left unmodified: soft clips, indels, duplicate-marked reads, unmapped mates,
    /// `N` bases, and MD strings emitted by an aligner rather than by us. Provenance in
    /// `eidolon-core/test_data/HG00096.chr20_1Mb.README.md`.
    ///
    /// KNOWN ANSWER, computed outside eidolon by an awk walk over MD+CIGAR and corroborated by
    /// `sum(NM) - sum(indel bases)` from the aligner's own tags: the A row of the raw count
    /// matrix is `[0, 97, 100, 53]`, 250 substitutions from a reference A. Normalized that is
    /// C 97/250, G 100/250, T 53/250, which as a CDF is `[0, 0.388, 0.788, 1.0]`.
    ///
    /// The assertion is on the fitted row, not merely on "a model was produced": the whole
    /// failure mode this file keeps guarding against is a model that looks trained and is not.
    #[test]
    fn test_runner_infers_from_a_real_aligner_bam() {
        let bam_path = PathBuf::from(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../eidolon-core/test_data/HG00096.chr20_1Mb.bam"
        ));
        assert!(
            bam_path.is_file(),
            "real-data fixture missing: {bam_path:?}"
        );

        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("test.fastq");
        make_test_fastq(&fastq_path, 20, 4);
        let output_path = temp.path().join("model.json.gz");

        let config = RunConfiguration {
            fastq_file: fastq_path,
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: None,
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
        };
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let a_row = row_cdf(
            model.transition_distros(),
            eidolon_core::structs::nucleotides::Nucleotide::A,
        );
        assert_row_cdf_eq(
            &a_row,
            &[0.0, 97.0 / 250.0, 197.0 / 250.0, 1.0],
            "A row fitted from the real HG00096 BAM",
        );

        // MUST DIFFER from the default. Without this the test would pass just as happily if
        // inference had silently fallen back -- which is exactly the #529 defect, and exactly
        // the kind of pass this repo has been fooled by before.
        let default_row = default_a_row_cdf();
        assert!(
            (a_row[1] - default_row[1]).abs() > 1e-6,
            "fitted A row {a_row:?} is indistinguishable from the default {default_row:?}; \
             inference did not actually fire"
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
            fastq_file_r2: None,
            output_file: output_path.clone(),
            overwrite_output: true,
            max_reads: 0,
            qual_offset: 33,
            max_model_read_length: 1000,
            fit_quality_degradation: false,
            degradation_tail_window: 50,
            degradation_tail_cut: 25,
            binned_quality_bins: None,
            bam_file: Some(bam_path),
            transition_matrix_file: Some(tsv_path),
            known_variants_vcf: None,
            bam_method: BamMethod::Mismatch,
            bam_min_mapq: 20,
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
            seen.insert(result);
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
