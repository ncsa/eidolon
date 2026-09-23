//! #734 — both mates' quality must degrade ALONG THE READ, not just R1's.
//!
//! CORRECTNESS CRITERION. A quality score belongs to a sequencing cycle. Cycle 1 is the first
//! base the sequencer reads, and quality falls with cycle number — on R1 and on R2 alike,
//! because each mate is sequenced from its own 5' end. So for either mate, the mean quality of
//! the first cycles must EXCEED that of the last, and the observed sequencing-error rate must
//! move the other way. It is FALSIFIED if a mate's profile rises with cycle.
//!
//! This is a known answer independent of the emitter: the shipped default quality model's
//! profile declines, so reads drawn from it must decline. Nothing here reads a fitted model.
//!
//! WHY THE ERROR RATE IS ASSERTED AND NOT JUST THE Q STRING. Errors are injected from each
//! base's own quality (`fastq_tools.rs:1110`), so a fix that reorders the emitted quality
//! string while leaving the error draw mirrored would satisfy a Q-only test and still produce
//! reads whose errors sit at the wrong end. The two must move together.
//!
//! THE MEASURED DEFECT, on the shipped default model, ecoli, 150 bp paired, n=30,935 per mate:
//!
//! | | cycle 1 | cycle 150 |
//! |---|---|---|
//! | R1 | 35.49 | 28.34 |
//! | R2 | **28.50** | **35.47** |
//!
//! `reverse_complement_record` reversed `quality_scores` along with `sequence` and `cigar_ops`.
//! The first two are reference-coordinate properties and the third is not.
//!
//! Fixture: ecoli (4.6 Mb, one contig). Nothing here needs a second contig, and a real
//! chromosome's statistics are what make a per-cycle rate worth reading.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, read_gzip_fastq_lines};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const READ_LEN: usize = 150;
/// A decile is wide enough to average out per-cycle noise and narrow enough that the two ends
/// are genuinely different parts of the read.
const EDGE: usize = READ_LEN / 10;

fn ecoli() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
}

fn reference_contigs() -> HashMap<String, String> {
    let text = std::fs::read_to_string(ecoli()).unwrap();
    let mut out: HashMap<String, String> = HashMap::new();
    let mut current = String::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(header) = line.strip_prefix('>') {
            current = header.split_whitespace().next().unwrap_or("").to_string();
            out.insert(current.clone(), String::new());
        } else if !current.is_empty() {
            out.get_mut(&current).unwrap().push_str(line);
        }
    }
    out
}

struct Run {
    work: PathBuf,
    _dir: tempfile::TempDir,
}

/// Paired reads from the SHIPPED default models, with no variants of any kind, so the only
/// thing shaping quality or mismatches is the sequencing-error model.
fn generate() -> Run {
    let (dir, work) = fresh_workdir();
    let mut config = GenReadsConfig::new(ecoli(), work.clone(), "orient");
    config.read_len = READ_LEN;
    config.coverage = 2;
    config.paired_ended = true;
    // A fragment must hold both mates end to end, or this measures adapter readthrough.
    config.fragment_mean = Some(400.0);
    config.fragment_st_dev = Some(40.0);
    config.produce_fastq = true;
    config.produce_bam = true;
    config.produce_vcf = false;
    config.mutation_rate = Some(0.0);
    config.sv_rate_scale = Some(0.0);
    config.rng_seed = "734".to_string();
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    Run { work, _dir: dir }
}

/// Mean quality at each cycle of the emitted read, and the number of reads behind it. The
/// denominator is returned rather than assumed.
fn per_cycle_quality(path: &Path) -> (Vec<f64>, usize) {
    let lines = read_gzip_fastq_lines(path);
    let mut sums = vec![0f64; READ_LEN];
    let mut n = 0usize;
    for q in lines.iter().skip(3).step_by(4) {
        if q.len() != READ_LEN {
            continue;
        }
        n += 1;
        for (i, b) in q.bytes().enumerate() {
            sums[i] += (b - 33) as f64;
        }
    }
    assert!(n > 1_000, "only {n} full-length reads; too few to profile");
    (sums.iter().map(|s| s / n as f64).collect(), n)
}

fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

/// Observed sequencing-error rate at each CYCLE, per mate, read off the golden BAM.
///
/// Restricted to records that are one `M` op spanning the whole read: those carry no indel, so
/// read index maps to reference offset directly. Selecting on the CIGAR is selecting on an
/// input property, not on the outcome being measured. For a reverse-strand record the BAM
/// stores everything reference-oriented, so cycle `i` sits at index `len - 1 - i`.
///
/// Returns (r1_rate_per_cycle, n_r1, r2_rate_per_cycle, n_r2).
fn per_cycle_error_rate(bam_path: &Path) -> (Vec<f64>, usize, Vec<f64>, usize) {
    let refs = reference_contigs();
    let mut reader = bam::io::reader::Builder.build_from_path(bam_path).unwrap();
    let header = reader.read_header().unwrap();

    let mut err = [vec![0f64; READ_LEN], vec![0f64; READ_LEN]];
    let mut n = [0usize, 0usize];

    for result in reader.records() {
        let record = result.unwrap();
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        // One M op spanning the read, or we cannot map index to reference offset.
        let ops: Vec<_> = record.cigar().iter().map(|o| o.unwrap()).collect();
        if ops.len() != 1 || ops[0].kind() != Kind::Match || ops[0].len() != READ_LEN {
            continue;
        }
        let ref_id = record.reference_sequence_id().transpose().unwrap().unwrap();
        let name = header
            .reference_sequences()
            .get_index(ref_id)
            .map(|(n, _)| n.to_string())
            .unwrap();
        let Some(contig) = refs.get(&name) else {
            continue;
        };
        let begin = usize::from(start) - 1;
        if begin + READ_LEN > contig.len() {
            continue;
        }
        let refseq = &contig.as_bytes()[begin..begin + READ_LEN];
        let seq: Vec<u8> = record.sequence().iter().collect();
        if seq.len() != READ_LEN {
            continue;
        }

        let flags = record.flags();
        let mate = usize::from(flags.is_last_segment());
        let reverse = flags.is_reverse_complemented();
        n[mate] += 1;
        for i in 0..READ_LEN {
            if !seq[i].eq_ignore_ascii_case(&refseq[i]) {
                // BAM index -> cycle. Reference-oriented storage means a reverse read's
                // cycle 1 is its LAST stored base.
                let cycle = if reverse { READ_LEN - 1 - i } else { i };
                err[mate][cycle] += 1.0;
            }
        }
    }

    assert!(
        n[0] > 500 && n[1] > 500,
        "too few indel-free records to rate: R1 {} / R2 {}",
        n[0],
        n[1]
    );
    let r1: Vec<f64> = err[0].iter().map(|e| e / n[0] as f64).collect();
    let r2: Vec<f64> = err[1].iter().map(|e| e / n[1] as f64).collect();
    (r1, n[0], r2, n[1])
}

/// THE CRITERION. Both mates must get worse along the read, because both are sequenced from
/// their own 5' end.
#[test]
fn both_mates_degrade_along_the_read() {
    let run = generate();
    for (label, file) in [("R1", "orient_r1.fastq.gz"), ("R2", "orient_r2.fastq.gz")] {
        let (profile, n) = per_cycle_quality(&run.work.join(file));
        let head = mean(&profile[..EDGE]);
        let tail = mean(&profile[READ_LEN - EDGE..]);
        eprintln!(
            "{label}: n={n} cycle 1={:.2} first{EDGE}={head:.2} last{EDGE}={tail:.2} cycle {READ_LEN}={:.2}",
            profile[0],
            profile[READ_LEN - 1]
        );
        assert!(
            head > tail + 1.0,
            "{label} quality must FALL with cycle: first {EDGE} cycles average Q{head:.2}, \
             last {EDGE} average Q{tail:.2}. A mate that improves along the read has its \
             quality array in reference order rather than cycle order (#734)."
        );
    }
}

/// Errors must follow quality. A fix that reorders the emitted Q string and leaves the error
/// draw mirrored passes the test above and still puts R2's errors at the wrong end.
#[test]
fn sequencing_errors_rise_with_cycle_on_both_mates() {
    let run = generate();
    let (r1, n1, r2, n2) = per_cycle_error_rate(&run.work.join("orient.bam"));
    for (label, rate, n) in [("R1", &r1, n1), ("R2", &r2, n2)] {
        let head = mean(&rate[..EDGE]);
        let tail = mean(&rate[READ_LEN - EDGE..]);
        eprintln!(
            "{label}: n={n} indel-free records, error rate first{EDGE}={:.5} last{EDGE}={:.5}",
            head, tail
        );
        assert!(
            tail > head * 1.5,
            "{label} sequencing errors must concentrate at the END of the read, where quality \
             is worst: first {EDGE} cycles {head:.5}, last {EDGE} cycles {tail:.5}. Errors are \
             drawn from each base's own quality, so this moves with the profile (#734)."
        );
    }
}

/// MUST NOT FIRE: the golden BAM stores a reverse-strand record reference-oriented, which is
/// the reverse of the emitted FASTQ. That relationship is correct today and the fix must not
/// disturb it — it is the check that catches the fix being applied in the wrong place.
#[test]
fn the_golden_bam_stays_the_mirror_of_the_fastq_for_reverse_reads() {
    let run = generate();
    let mut fastq: HashMap<String, String> = HashMap::new();
    for file in ["orient_r1.fastq.gz", "orient_r2.fastq.gz"] {
        let lines = read_gzip_fastq_lines(&run.work.join(file));
        for chunk in lines.chunks(4) {
            if chunk.len() < 4 {
                break;
            }
            fastq.insert(
                chunk[0].trim_start_matches('@').trim().to_string(),
                chunk[3].clone(),
            );
        }
    }

    let mut reader = bam::io::reader::Builder
        .build_from_path(run.work.join("orient.bam"))
        .unwrap();
    let _ = reader.read_header().unwrap();
    let mut checked = 0usize;
    for result in reader.records() {
        let record = result.unwrap();
        if !record.flags().is_reverse_complemented() {
            continue;
        }
        let name = String::from_utf8_lossy(record.name().unwrap().as_ref()).to_string();
        let mate = if record.flags().is_last_segment() {
            2
        } else {
            1
        };
        let Some(fq_qual) = fastq.get(&format!("{name}/{mate}")) else {
            continue;
        };
        let bam_qual: String = record
            .quality_scores()
            .iter()
            .map(|q| (q + 33) as char)
            .collect();
        if bam_qual.len() != fq_qual.len() {
            continue;
        }
        let reversed: String = fq_qual.chars().rev().collect();
        assert_eq!(
            bam_qual, reversed,
            "a reverse-strand record's BAM QUAL must be the reverse of its FASTQ QUAL \
             (SAM spec: SEQ and QUAL are stored reference-oriented). Read {name}/{mate}."
        );
        checked += 1;
    }
    assert!(
        checked > 1_000,
        "only {checked} reverse-strand records compared; too few to call this checked"
    );
}
