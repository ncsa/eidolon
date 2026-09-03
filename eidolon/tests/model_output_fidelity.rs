//! Fidelity tests: a model FILE built by a builder must shape gen-reads output,
//! not merely load. Companion to `model_fragment_fidelity.rs` (fragment length);
//! this file covers the sequencing-error and mutation models.
//!
//! The Delta `model_builders.sbatch` round-trip only smoke-tests usability (reads
//! exist and are correctly-lengthed). These tests close that gap by driving a
//! measurable, distinctive signal from the built model through to the output.

mod common;

use std::fs;
use std::io::{Read as _, Write as _};
use std::path::Path;

use common::{eidolon, h1n1_reference};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
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
use serde_json::Value;

fn write_yaml(dir: &Path, tag: &str, body: &str) -> std::path::PathBuf {
    let p = dir.join(format!("cfg_{tag}.yml"));
    fs::write(&p, body).unwrap();
    p
}

/// Write a gzipped single-end FASTQ where every base has the same Phred quality.
fn write_uniform_quality_fastq_gz(path: &Path, n_reads: usize, read_len: usize, phred: u8) {
    let seq: String = "ACGT".chars().cycle().take(read_len).collect();
    let qual: String = std::iter::repeat(char::from(phred + 33))
        .take(read_len)
        .collect();
    let f = fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    for i in 0..n_reads {
        writeln!(enc, "@read{i}\n{seq}\n+\n{qual}").unwrap();
    }
    enc.finish().unwrap();
}

/// Mean Phred quality over the base-quality lines of a (single-end) FASTQ.gz.
fn mean_fastq_quality(path: &Path) -> f64 {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let mut sum = 0u64;
    let mut n = 0u64;
    for (i, line) in raw.lines().enumerate() {
        if i % 4 == 3 {
            for b in line.bytes() {
                sum += (b - 33) as u64;
                n += 1;
            }
        }
    }
    assert!(n > 0, "no quality lines in {}", path.display());
    sum as f64 / n as f64
}

/// A training FASTQ whose quality DECAYS along the read: `start_phred` at cycle 0 falling
/// linearly to `end_phred` at the last cycle. The uniform fixture cannot distinguish a
/// model that reproduces a positional profile from one that returns a constant, because
/// with one quality option every implementation looks identical.
fn write_gradient_quality_fastq_gz(
    path: &Path,
    n_reads: usize,
    read_len: usize,
    start_phred: u8,
    end_phred: u8,
) {
    let seq: String = "ACGT".chars().cycle().take(read_len).collect();
    let qual: String = (0..read_len)
        .map(|i| {
            let frac = i as f64 / (read_len - 1) as f64;
            let q = start_phred as f64 + (end_phred as f64 - start_phred as f64) * frac;
            char::from(q.round() as u8 + 33)
        })
        .collect();
    let f = fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    for i in 0..n_reads {
        writeln!(enc, "@read{i}\n{seq}\n+\n{qual}").unwrap();
    }
    enc.finish().unwrap();
}

/// Mean Phred quality at each cycle across all reads of a single-end FASTQ.gz.
fn per_cycle_mean_quality(path: &Path, read_len: usize) -> Vec<f64> {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let mut sum = vec![0u64; read_len];
    let mut n = vec![0u64; read_len];
    for (i, line) in raw.lines().enumerate() {
        if i % 4 == 3 {
            for (cycle, b) in line.bytes().enumerate().take(read_len) {
                sum[cycle] += (b - 33) as u64;
                n[cycle] += 1;
            }
        }
    }
    assert!(n[0] > 0, "no quality lines in {}", path.display());
    sum.iter()
        .zip(&n)
        .map(|(s, c)| if *c == 0 { 0.0 } else { *s as f64 / *c as f64 })
        .collect()
}

/// Build a sequencing-error model from a FASTQ whose bases are ALL Phred 35, then
/// simulate with that model file and confirm the output read qualities are ~35 —
/// i.e. the trained quality profile reaches the reads, not a built-in default.
#[test]
fn built_seq_error_model_drives_output_quality() {
    let tmp = tempfile::tempdir().unwrap();
    let read_len = 100;
    let trained_q: u8 = 35;

    // 1. training FASTQ: uniform Phred 35 everywhere.
    let fq = tmp.path().join("train.fastq.gz");
    write_uniform_quality_fastq_gz(&fq, 500, read_len, trained_q);

    // 2. build the sequencing-error model.
    let model = tmp.path().join("seq_error.json.gz");
    let build_cfg = write_yaml(
        tmp.path(),
        "seqerr_build",
        &format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\n",
            fq.display(),
            model.display()
        ),
    );
    eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&build_cfg)
        .assert()
        .success();

    // 3. simulate single-end reads using the model file.
    let sim_cfg = write_yaml(
        tmp.path(),
        "seqerr_sim",
        &format!(
            "reference: {ref}\nread_len: {rl}\ncoverage: 20\npaired_ended: false\n\
             sequence_error_model: {model}\nproduce_fastq: true\nproduce_bam: false\n\
             produce_vcf: false\noverwrite_output: true\noutput_dir: {out}\n\
             output_filename: rt\nrng_seed: seqerr fidelity\nnum_threads: 1\n",
            ref = h1n1_reference().display(),
            rl = read_len,
            model = model.display(),
            out = tmp.path().display(),
        ),
    );
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&sim_cfg)
        .assert()
        .success();

    // 4. output qualities should track the trained Phred 35.
    let out_q = mean_fastq_quality(&tmp.path().join("rt_r1.fastq.gz"));
    eprintln!(
        "[fidelity] trained quality = {trained_q}; simulated output mean quality = {out_q:.1}"
    );
    assert!(
        (out_q - trained_q as f64).abs() <= 3.0,
        "output mean quality {out_q:.1} does not track the trained Phred {trained_q} — \
         a built sequence_error_model is NOT reaching gen-reads output"
    );
}

// ── mutation model ─────────────────────────────────────────────────────────────

/// Read the top-level fitted `mutation_rate` out of a gzipped mutation model.
fn fitted_mutation_rate(path: &Path) -> f64 {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let v: Value = serde_json::from_str(&raw).unwrap();
    v["mutation_rate"]
        .as_f64()
        .expect("mutation model should carry a numeric mutation_rate")
}

/// Count variant (non-header, non-blank) lines in a gzipped VCF.
fn count_vcf_variants(path: &Path) -> usize {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    raw.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .count()
}

/// Read one contig's sequence out of a FASTA (concatenating wrapped lines).
fn read_contig(fa: &Path, want: &str) -> Vec<u8> {
    let data = fs::read_to_string(fa).unwrap();
    let mut seq = Vec::new();
    let mut in_contig = false;
    for line in data.lines() {
        if let Some(h) = line.strip_prefix('>') {
            in_contig = h.split_whitespace().next() == Some(want);
        } else if in_contig {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    assert!(
        !seq.is_empty(),
        "contig {want} not found in {}",
        fa.display()
    );
    seq
}

/// Write a VCF of SNPs at every `step`-th interior position of `contig`, with REF
/// matching the reference (so gen-mut-model doesn't skip them) and a fixed ALT.
/// Returns the number of records written.
fn write_dense_snp_vcf(path: &Path, contig: &str, seq: &[u8], step: usize) -> usize {
    let alt = |b: u8| match b.to_ascii_uppercase() {
        b'A' => Some('T'),
        b'C' => Some('G'),
        b'G' => Some('C'),
        b'T' => Some('A'),
        _ => None,
    };
    let mut f = fs::File::create(path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .unwrap();
    let mut count = 0;
    let mut pos = 100usize; // 1-based; stay clear of contig edges
    while pos < seq.len() - 100 {
        let refb = seq[pos - 1].to_ascii_uppercase();
        if let Some(a) = alt(refb) {
            writeln!(
                f,
                "{contig}\t{pos}\t.\t{}\t{a}\t60\tPASS\t.\tGT\t0/1",
                refb as char
            )
            .unwrap();
            count += 1;
        }
        pos += step;
    }
    count
}

fn build_mut_model(tmp: &Path, tag: &str, vcf: &Path) -> std::path::PathBuf {
    let model = tmp.join(format!("mut_{tag}.json.gz"));
    let cfg = write_yaml(
        tmp,
        &format!("mut_build_{tag}"),
        &format!(
            "reference: {ref}\nvcf_file: {vcf}\noutput_file: {model}\n\
             overwrite_output: true\nbed_file: .\n",
            ref = h1n1_reference().display(),
            vcf = vcf.display(),
            model = model.display(),
        ),
    );
    eidolon()
        .args(["gen-mut-model", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    model
}

fn simulate_variant_count(tmp: &Path, tag: &str, model: &Path) -> usize {
    let cfg = write_yaml(
        tmp,
        &format!("mut_sim_{tag}"),
        &format!(
            "reference: {ref}\nread_len: 100\ncoverage: 5\npaired_ended: false\nploidy: 1\n\
             mutation_model: {model}\nproduce_vcf: true\nproduce_fastq: false\n\
             produce_bam: false\noverwrite_output: true\noutput_dir: {out}\n\
             output_filename: rt_{tag}\nrng_seed: mut fidelity {tag}\nnum_threads: 1\n",
            ref = h1n1_reference().display(),
            model = model.display(),
            out = tmp.display(),
        ),
    );
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    count_vcf_variants(&tmp.join(format!("rt_{tag}.vcf.gz")))
}

/// Build two mutation models with very different fitted rates (a sparse fixture vs a
/// dense synthetic VCF) and confirm the output variant count scales with the rate —
/// i.e. the fitted `mutation_rate` actually governs how many mutations gen-reads
/// injects, not a fixed default.
#[test]
fn built_mutation_model_rate_drives_output_variant_count() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // low-rate model: a sparse VCF (3 SNPs).
    let lo_vcf = dir.join("sparse.vcf");
    let seq_pb2 = read_contig(&h1n1_reference(), "H1N1_PB2");
    write_dense_snp_vcf(&lo_vcf, "H1N1_PB2", &seq_pb2, 900); // ~2 SNPs
    let lo_model = build_mut_model(dir, "lo", &lo_vcf);

    // high-rate model: a dense VCF (a SNP every 5 bp over PB2).
    let hi_vcf = dir.join("dense.vcf");
    let hi_written = write_dense_snp_vcf(&hi_vcf, "H1N1_PB2", &seq_pb2, 5);
    let hi_model = build_mut_model(dir, "hi", &hi_vcf);

    let rate_lo = fitted_mutation_rate(&lo_model);
    let rate_hi = fitted_mutation_rate(&hi_model);
    eprintln!(
        "[fidelity] dense VCF wrote {hi_written} SNPs; fitted rate lo={rate_lo:.2e} hi={rate_hi:.2e}"
    );
    assert!(
        rate_hi > rate_lo * 10.0,
        "test setup: hi rate {rate_hi:.2e} should dwarf lo rate {rate_lo:.2e}"
    );

    let n_lo = simulate_variant_count(dir, "lo", &lo_model);
    let n_hi = simulate_variant_count(dir, "hi", &hi_model);
    eprintln!(
        "[fidelity] output variants: lo={n_lo} hi={n_hi} (rate ratio ~{:.0}x)",
        rate_hi / rate_lo.max(1e-9)
    );

    assert!(
        n_hi > n_lo * 5,
        "output variant count did not scale with the fitted mutation_rate \
         (lo model → {n_lo} variants, hi model → {n_hi}); the built mutation_model's \
         rate is NOT governing gen-reads output"
    );
}

// ── GC-bias model ────────────────────────────────────────────────────────────

/// Write a 2-contig reference: `lowgc` (~20% GC) and `highgc` (~80% GC), each
/// `reps`*10 bp on a single line, plus a matching `.fai`. A wide, clean GC split
/// (H1N1's ~40% range can't populate distinct GC bins) makes the model deterministic.
/// Returns (path, contig_len).
fn write_gc_reference(dir: &Path, reps: usize) -> (std::path::PathBuf, usize) {
    let low = "ATATATATGC".repeat(reps); //  2/10 GC = 20%
    let high = "GCGCGCGCTA".repeat(reps); // 8/10 GC = 80%
    let path = dir.join("gc_ref.fa");
    let mut buf = String::new();
    let mut fai = String::new();
    let mut offset = 0usize;
    for (name, seq) in [("lowgc", &low), ("highgc", &high)] {
        let hdr = format!(">{name}\n");
        buf.push_str(&hdr);
        offset += hdr.len();
        let seq_offset = offset;
        buf.push_str(seq);
        buf.push('\n');
        offset += seq.len() + 1;
        // name, length, seq byte offset, bases-per-line, bytes-per-line (+newline)
        fai.push_str(&format!(
            "{name}\t{}\t{seq_offset}\t{}\t{}\n",
            seq.len(),
            seq.len(),
            seq.len() + 1
        ));
    }
    fs::write(&path, &buf).unwrap();
    fs::write(dir.join("gc_ref.fa.fai"), &fai).unwrap();
    (path, low.len())
}

/// Write an unpaired-read BAM over the given contigs. `reads` are (contig_id, start).
fn write_coverage_bam(
    path: &Path,
    contigs: &[(&str, usize)],
    reads: &[(usize, usize)],
    read_len: usize,
) {
    let mut hb = sam::Header::builder();
    for (name, len) in contigs {
        hb = hb.add_reference_sequence(
            name.as_bytes().to_vec(),
            Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(*len).unwrap()),
        );
    }
    let header = hb.build();
    let mut writer = bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
    let seq = vec![b'A'; read_len];
    for &(cid, start) in reads {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::empty();
        *r.cigar_mut() = cigar.clone();
        *r.reference_sequence_id_mut() = Some(cid);
        *r.alignment_start_mut() = Position::new(start);
        *r.sequence_mut() = Sequence::from(seq.as_slice());
        *r.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
        writer.write_alignment_record(&header, &r).unwrap();
    }
}

fn gc_weights(path: &Path) -> Vec<f64> {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let v: Value = serde_json::from_str(&raw).unwrap();
    v["weights_by_percent_gc"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

/// Classify single-end FASTQ.gz reads by their own GC content: returns
/// (reads with GC < 40%, reads with GC > 60%). The 20/80% contigs sit well outside
/// this band, so a few simulated errors can't misclassify them.
fn classify_reads_by_gc(fastq_gz: &Path) -> (usize, usize) {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(fastq_gz).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let (mut low, mut high) = (0usize, 0usize);
    for (i, line) in raw.lines().enumerate() {
        if i % 4 == 1 && !line.is_empty() {
            let gc = line
                .bytes()
                .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
                .count();
            let frac = gc as f64 / line.len() as f64;
            if frac < 0.40 {
                low += 1;
            } else if frac > 0.60 {
                high += 1;
            }
        }
    }
    (low, high)
}

/// Build a GC-bias model from a BAM that covers the low-GC contig heavily and the
/// high-GC contig barely, then simulate and confirm high-GC reads are strongly
/// depleted in the output — i.e. the built weight table shapes which regions get
/// sequenced, not just whether the run succeeds.
#[test]
fn built_gc_bias_model_depletes_disfavored_gc_in_output() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let reps = 500;
    let (reference, clen) = write_gc_reference(dir, reps); // two 5 kb contigs
    let read_len = 100;

    // training coverage: dense over lowgc (id 0), almost none over highgc (id 1).
    let mut reads: Vec<(usize, usize)> = Vec::new();
    let mut start = 1;
    while start + read_len <= clen {
        reads.push((0, start)); // lowgc, stride 20 → ~5x coverage
        start += 20;
    }
    for s in [1usize, clen / 3, 2 * clen / 3] {
        reads.push((1, s)); // highgc: 3 reads → ~0 coverage
    }
    let cov_bam = dir.join("cov.bam");
    write_coverage_bam(
        &cov_bam,
        &[("lowgc", clen), ("highgc", clen)],
        &reads,
        read_len,
    );

    // build the GC-bias model.
    let model = dir.join("gc.json.gz");
    let build_cfg = write_yaml(
        dir,
        "gc_build",
        &format!(
            "reference: {ref}\nbam_file: {bam}\noutput_file: {model}\noverwrite_output: true\n\
             min_mapq: 0\nbed_file: .\nwindow_size: 100\nwindow_stride: 100\n\
             min_windows_per_bin: 1\n",
            ref = reference.display(),
            bam = cov_bam.display(),
            model = model.display(),
        ),
    );
    eidolon()
        .args(["gen-gc-bias-model", "-c"])
        .arg(&build_cfg)
        .assert()
        .success();

    // builder side: bin 20 up-weighted, bin 80 driven to ~0.
    let w = gc_weights(&model);
    eprintln!(
        "[fidelity] gc weights: w[20]={:.3} w[80]={:.3}",
        w[20], w[80]
    );
    assert!(
        w[20] > 1.2,
        "low-GC bin should be up-weighted, got w[20]={:.3}",
        w[20]
    );
    assert!(
        w[80] < 0.2,
        "high-GC bin should be ~0, got w[80]={:.3}",
        w[80]
    );

    // simulate with the model and classify output reads by GC.
    let sim_cfg = write_yaml(
        dir,
        "gc_sim",
        &format!(
            "reference: {ref}\nread_len: {rl}\ncoverage: 20\npaired_ended: false\n\
             gc_bias_model: {model}\nproduce_fastq: true\nproduce_bam: false\n\
             produce_vcf: false\noverwrite_output: true\noutput_dir: {out}\n\
             output_filename: rt\nrng_seed: gc fidelity\nnum_threads: 1\n",
            ref = reference.display(),
            rl = read_len,
            model = model.display(),
            out = dir.display(),
        ),
    );
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&sim_cfg)
        .assert()
        .success();

    let (low_reads, high_reads) = classify_reads_by_gc(&dir.join("rt_r1.fastq.gz"));
    eprintln!("[fidelity] output reads: low-GC={low_reads} high-GC={high_reads}");
    assert!(low_reads > 0, "no low-GC reads produced at all");
    assert!(
        high_reads * 5 < low_reads,
        "high-GC reads ({high_reads}) not depleted vs low-GC ({low_reads}) — the built \
         gc_bias_model's weights are NOT shaping which regions gen-reads sequences"
    );
}

/// The trained quality profile must reproduce its POSITIONAL SHAPE, not merely its mean.
///
/// `built_seq_error_model_drives_output_quality` trains on a uniform Phred-35 fixture, so
/// the model carries a single quality option and any implementation returning a constant
/// satisfies it. Two mutations confirmed that: replacing the per-position transition
/// matrix with row 0 (quality stops decaying along the read), and returning
/// `quality_score_options[0]` for every base. Both passed the whole workspace — nothing
/// checked that quality varies with cycle at all.
///
/// Real Illumina quality decays along the read, and callers weight bases by it, so a
/// simulator that emits flat quality is not modelling the thing it claims to.
#[test]
fn built_seq_error_model_reproduces_the_positional_quality_profile() {
    let tmp = tempfile::tempdir().unwrap();
    let read_len = 100;
    let (start_q, end_q) = (40u8, 20u8);

    let fq = tmp.path().join("train_gradient.fastq.gz");
    write_gradient_quality_fastq_gz(&fq, 500, read_len, start_q, end_q);

    let model = tmp.path().join("seq_error_gradient.json.gz");
    let build_cfg = write_yaml(
        tmp.path(),
        "seqerr_grad_build",
        &format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\n",
            fq.display(),
            model.display()
        ),
    );
    eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&build_cfg)
        .assert()
        .success();

    let sim_cfg = write_yaml(
        tmp.path(),
        "seqerr_grad_sim",
        &format!(
            "reference: {ref}\nread_len: {rl}\ncoverage: 20\npaired_ended: false\n\
             sequence_error_model: {model}\nproduce_fastq: true\nproduce_bam: false\n\
             produce_vcf: false\noverwrite_output: true\noutput_dir: {out}\n\
             output_filename: grad\nrng_seed: seqerr gradient\nnum_threads: 1\n",
            ref = h1n1_reference().display(),
            rl = read_len,
            model = model.display(),
            out = tmp.path().display(),
        ),
    );
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&sim_cfg)
        .assert()
        .success();

    let profile = per_cycle_mean_quality(&tmp.path().join("grad_r1.fastq.gz"), read_len);
    let head: f64 = profile[..10].iter().sum::<f64>() / 10.0;
    let tail: f64 = profile[read_len - 10..].iter().sum::<f64>() / 10.0;
    eprintln!("[fidelity] trained {start_q} -> {end_q}; simulated head {head:.1}, tail {tail:.1}");

    // The decay itself. A constant-quality implementation lands head == tail and fails
    // here whatever value it picks, which is the point: the mean alone cannot tell them
    // apart, only the shape can.
    let trained_drop = (start_q - end_q) as f64;
    assert!(
        head - tail >= trained_drop * 0.5,
        "quality falls only {:.1} Phred from cycle 0 to cycle {} but the training data \
         falls {trained_drop:.0} — the positional profile is not reaching the reads \
         (head {head:.1}, tail {tail:.1})",
        head - tail,
        read_len - 1
    );
    // ...and each end tracks the trained value, so "decays" cannot be satisfied by an
    // arbitrary downward ramp unrelated to the model.
    assert!(
        (head - start_q as f64).abs() <= 5.0,
        "cycle-0 quality {head:.1} does not track the trained {start_q}"
    );
    assert!(
        (tail - end_q as f64).abs() <= 5.0,
        "final-cycle quality {tail:.1} does not track the trained {end_q}"
    );
}

// ── sequencing-error transition matrix ─────────────────────────────────────────

/// Write a single-contig all-A reference plus its `.fai`.
///
/// A homopolymer reference makes every true base an A, so every sequencing-error
/// substitution must be drawn from the matrix's A row. That is what lets the assertion
/// below name an expected base without reconstructing each read's alignment.
fn write_homopolymer_reference(dir: &Path, len: usize) -> std::path::PathBuf {
    let seq = "A".repeat(len);
    let path = dir.join("polyA.fa");
    let hdr = ">polyA\n";
    fs::write(&path, format!("{hdr}{seq}\n")).unwrap();
    fs::write(
        dir.join("polyA.fa.fai"),
        format!("polyA\t{}\t{}\t{}\t{}\n", len, hdr.len(), len, len + 1),
    )
    .unwrap();
    path
}

/// Base composition of a gzipped FASTQ's sequence lines, as (A, C, G, T, other).
fn base_counts(path: &Path) -> (usize, usize, usize, usize, usize) {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let mut c = (0usize, 0usize, 0usize, 0usize, 0usize);
    for (i, line) in raw.lines().enumerate() {
        if i % 4 != 1 {
            continue; // sequence lines only
        }
        for b in line.bytes() {
            match b {
                b'A' => c.0 += 1,
                b'C' => c.1 += 1,
                b'G' => c.2 += 1,
                b'T' => c.3 += 1,
                _ => c.4 += 1,
            }
        }
    }
    c
}

/// Build a sequencing-error model (optionally with a forced transition matrix), simulate
/// over an all-A reference with mutations disabled, and return the output base counts.
fn simulate_over_polya(tmp: &Path, tag: &str, tsv: Option<&Path>) -> (usize, usize, usize, usize) {
    let read_len = 100;
    // Phred 8 → ~0.158 per-base error probability, so errors are plentiful.
    let fq = tmp.join(format!("{tag}_train.fastq.gz"));
    write_uniform_quality_fastq_gz(&fq, 500, read_len, 8);

    let model = tmp.join(format!("{tag}_model.json.gz"));
    let mut body = format!(
        "fastq_file: {}\noutput_file: {}\noverwrite_output: true\n",
        fq.display(),
        model.display()
    );
    if let Some(tsv) = tsv {
        body.push_str(&format!("transition_matrix_file: {}\n", tsv.display()));
    }
    let build_cfg = write_yaml(tmp, &format!("{tag}_build"), &body);
    eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&build_cfg)
        .assert()
        .success();

    let reference = write_homopolymer_reference(tmp, 20_000);
    let sim_cfg = write_yaml(
        tmp,
        &format!("{tag}_sim"),
        &format!(
            "reference: {ref}\nread_len: {rl}\ncoverage: 30\npaired_ended: false\n\
             sequence_error_model: {model}\nmutation_rate: 0.0\nproduce_fastq: true\n\
             produce_bam: false\nproduce_vcf: false\noverwrite_output: true\n\
             output_dir: {out}\noutput_filename: {tag}\nrng_seed: transition fidelity\n\
             num_threads: 1\n",
            ref = reference.display(),
            rl = read_len,
            model = model.display(),
            out = tmp.display(),
        ),
    );
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&sim_cfg)
        .assert()
        .success();

    let (a, c, g, t, other) = base_counts(&tmp.join(format!("{tag}_r1.fastq.gz")));
    assert_eq!(other, 0, "{tag}: unexpected non-ACGT bases in output");
    (a, c, g, t)
}

/// The SNP transition matrix in a built sequencing-error model must decide WHICH base a
/// sequencing error substitutes in the output reads — not merely be present in the file.
///
/// The unit tests in `gen_seq_error_model::utils::runner` establish that a BAM's MD-tagged
/// mismatches produce the right matrix, and that a TSV overrides a BAM. Neither can show
/// the matrix reaching the reads, which is the claim that matters: a model file carrying a
/// perfect matrix is worth nothing if `gen-reads` ignores it.
///
/// Setup: an all-A reference, so every true base is an A and every substitution must come
/// from the matrix's A row. Two runs at the same seed — one with the default matrix
/// (0.4918 C / 0.3377 G / 0.1705 T), one with the A row forced entirely to T.
///
/// This is differential rather than absolute because a forced run does NOT reach 100% T.
/// `insertion_bias` is uniform over ACGT, so indel errors put inserted bases into the
/// output that never consult the transition matrix. That floor is a few hundred C and G no
/// matter what the matrix says — correct behavior, and the reason an absolute `>98% T`
/// assertion would fail on working code. The control run pins where the substitution
/// spectrum sits without the override, so the comparison isolates the matrix's effect.
///
/// **Why the floor is ~39% here and not the model's 1%.** This comment used to read
/// "`indel_probability` is 0.4, so roughly 40% of errors are indels". #660 corrected that
/// constant to NEAT2's 0.01, and #661 then made the indel share depend on local
/// homopolymer run length. `simulate_over_polya` uses a 20,000-base all-A reference — a
/// homopolymer far past the curve's cap, and so its most enriched case: 0.01 x 39.20 =
/// 0.392. The old number is accidentally close to the new one for an entirely different
/// reason, which is worth stating plainly rather than leaving as a comment that happens to
/// still look right. A homopolymer-free fixture would have a ~0.6% floor instead.
#[test]
fn built_seq_error_transition_matrix_decides_the_substituted_base() {
    let tmp = tempfile::tempdir().unwrap();

    // A row: all weight on T. Diagonals are ignored by the loader. The remaining rows are
    // never consulted — an all-A reference has no other true base.
    let tsv = tmp.path().join("a_to_t.tsv");
    fs::write(
        &tsv,
        "A\tC\tG\tT\n\
         0.0\t0.0\t0.0\t1.0\n\
         1.0\t0.0\t0.0\t0.0\n\
         1.0\t0.0\t0.0\t0.0\n\
         1.0\t0.0\t0.0\t0.0\n",
    )
    .unwrap();

    let (_, c_def, g_def, t_def) = simulate_over_polya(tmp.path(), "ctrl", None);
    let (_, c_for, g_for, t_for) = simulate_over_polya(tmp.path(), "forced", Some(&tsv));

    let subs_def = c_def + g_def + t_def;
    let subs_for = c_for + g_for + t_for;
    let share_def = t_def as f64 / subs_def as f64;
    let share_for = t_for as f64 / subs_for as f64;
    eprintln!(
        "[fidelity] default matrix: C={c_def} G={g_def} T={t_def} → T share {:.3}\n\
         [fidelity] forced A→T:     C={c_for} G={g_for} T={t_for} → T share {:.3}",
        share_def, share_for
    );

    // Non-vacuity: with no substitutions at all, every ratio below is meaningless.
    assert!(
        subs_def > 500 && subs_for > 500,
        "too few substitutions to compare (default {subs_def}, forced {subs_for})"
    );

    // The control must look like the default matrix — mostly C, T in the minority.
    assert!(
        share_def < 0.35,
        "default-matrix run put {:.1}% of substitutions on T; the default A row is only \
         0.1705 T, so this run is not using the default matrix and the comparison below \
         proves nothing",
        share_def * 100.0
    );

    // Forcing the A row to T must dominate the spectrum. The residual C/G is the uniform
    // insertion floor described above, not a matrix that went unread.
    assert!(
        share_for > 0.80,
        "forcing the A row to T only moved the T share to {:.1}% (C={c_for}, G={g_for}) — \
         the model's transition matrix is NOT deciding the substituted base in output reads",
        share_for * 100.0
    );

    // State the causal claim directly: the override, not chance, moved the spectrum.
    assert!(
        share_for > share_def * 2.0,
        "T share barely moved between the default ({:.3}) and forced ({:.3}) runs",
        share_def,
        share_for
    );
}
