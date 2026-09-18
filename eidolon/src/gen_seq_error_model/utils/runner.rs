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

/// Read one mate's FASTQ into its own `MateCounts`, pooling the option-set and base counts.
///
/// `max_reads` applies PER MATE rather than across the pair: a shared budget would spend it
/// all on R1 and fit R2 from whatever was left, which is nothing at the default of 0-means-all
/// and is silently lopsided otherwise.
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
    let mate_start = *reads_processed;
    'records: loop {
        if config.max_reads > 0 && (*reads_processed) - mate_start >= config.max_reads {
            break;
        }
        // Line number of this record's header, 1-based, for error messages.
        let first_line = m.records_seen * 4 + 1;

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
        let mut populations: Vec<(&str, &Vec<Vec<Vec<usize>>>)> = Vec::new();
        // R1's healthy tensor is checked only when there is a second population to be ragged
        // against. On its own it defines `positions` and cannot fall short of itself.
        if config.fit_quality_degradation {
            populations.push(("healthy", &transition_counts));
            populations.push(("degraded", &transition_counts_deg));
        }
        if let Some(m) = r2.as_ref() {
            populations.push(("R1", &transition_counts));
            populations.push(("R2", &m.transition_counts));
            if config.fit_quality_degradation {
                populations.push(("R2 degraded", &m.transition_counts_deg));
            }
        }
        for (label, counts) in populations {
            let trained = trained_positions(counts);
            if trained < positions {
                return Err(GenSeqErrorModelError::ConfigurationError(format!(
                    "the {label} population covers {} bp, but the model covers \
                     {read_length} bp: positions {} to {read_length} carry no {label} \
                     observation at all. Filling them would give that population a uniform \
                     quality distribution -- fabricated data a reader cannot distinguish \
                     from a fit. Both populations must reach the model's read length: trim \
                     or filter the library so the two carry the same read lengths, or unset \
                     fit_quality_degradation to fit a single population.",
                    trained + 1,
                    trained + 2,
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
    let snp_transition_matrix: Option<TransitionMatrix> =
        if let Some(path) = &config.transition_matrix_file {
            info!("Loading SNP transition matrix from TSV: {:?}", path);
            Some(TransitionMatrix::from_tsv(path)?)
        } else if let Some(path) = &config.bam_file {
            info!("Inferring SNP transition matrix from BAM: {:?}", path);
            let counts = read_bam_transitions(path)?;
            let total_mismatches: usize = counts.iter().flatten().sum();
            if total_mismatches == 0 {
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

    /// THE regression for the second review finding on #698: `max_reads` smaller than the
    /// window the read length was drawn from.
    ///
    /// Known answer: only the ten 50 bp records are fitted, so the model is a 50 bp model.
    /// Under the sampling code it advertised 100 bp — the length came from records that
    /// `max_reads` excluded — and positions 51-100 were uniform noise between Q10 and Q40.
    #[test]
    fn max_reads_does_not_advertise_positions_it_did_not_fit() {
        let temp = tempfile::tempdir().unwrap();
        let fastq_path = temp.path().join("capped.fastq");
        let short_qual = "I".repeat(25) + &"+".repeat(25);
        let long_qual = "I".repeat(50) + &"+".repeat(50);
        write_two_phase_fastq(&fastq_path, 10, &short_qual, 100, &long_qual);
        let output_path = temp.path().join("model.json.gz");
        let mut config = make_config(fastq_path, output_path.clone());
        config.max_reads = 10;
        runner(&config).unwrap();

        let model = SequencingErrorModel::from_file(&output_path).unwrap();
        let q = model.quality_score_model();
        assert_eq!(
            q.assumed_read_length, 50,
            "only the ten 50 bp records were fitted, so the model describes 50 positions"
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
