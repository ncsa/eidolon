//! Shared helpers for eidolon integration tests.
//!
//! Every `tests/*.rs` file is compiled as its own crate root, so this module is included
//! via `mod common;` at the top of each. Items used by only one test file may still be
//! reported as `dead_code` from the perspective of crates that don't use them — hence the
//! blanket `#[allow(dead_code)]` below.

#![allow(dead_code)]

pub mod gate2;

use assert_cmd::Command;
use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ── binary invocation ────────────────────────────────────────────────────────

/// A fresh `Command` wired to the eidolon binary built by this workspace.
///
/// Wraps `assert_cmd::Command::cargo_bin("eidolon")` so callers don't have to import
/// assert_cmd themselves for the common case.
pub fn eidolon() -> Command {
    Command::cargo_bin("eidolon").expect("eidolon binary should be built by cargo test")
}

// ── reference fixtures ──────────────────────────────────────────────────────

/// Path to the bundled H1N1 reference shipped under `eidolon/test_data/`.
/// ~14 kb across three contigs — tiny enough for full-binary integration runs.
pub fn h1n1_reference() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/H1N1.fa",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

// ── FASTQ fixture generation ─────────────────────────────────────────────────

/// Writes a small uncompressed FASTQ with `n_reads` records of length `read_length`.
/// Quality scores cycle through ASCII '!' + 33 + (i % 10), giving Phred scores 33..=42.
/// Useful for training a sequencing error model in pipeline tests.
pub fn write_synthetic_fastq(path: &Path, n_reads: usize, read_length: usize) {
    let seq: String = "ACGT".chars().cycle().take(read_length).collect();
    let qual: String = (0..read_length)
        .map(|i| char::from_u32(('!' as u32) + 33 + (i % 10) as u32).unwrap())
        .collect();
    let mut f = std::fs::File::create(path).unwrap();
    for i in 0..n_reads {
        writeln!(f, "@read{}\n{}\n+\n{}", i, seq, qual).unwrap();
    }
}

/// Same as `write_synthetic_fastq`, but the output is gzip-compressed (.fastq.gz).
pub fn write_synthetic_fastq_gz(path: &Path, n_reads: usize, read_length: usize) {
    let seq: String = "ACGT".chars().cycle().take(read_length).collect();
    let qual: String = (0..read_length)
        .map(|i| char::from_u32(('!' as u32) + 33 + (i % 10) as u32).unwrap())
        .collect();
    let f = std::fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    for i in 0..n_reads {
        writeln!(enc, "@read{}\n{}\n+\n{}", i, seq, qual).unwrap();
    }
    enc.finish().unwrap();
}

// ── config-file builders ─────────────────────────────────────────────────────

/// Configuration parameters for a `gen-reads` config YAML. Defaults mirror what the
/// existing inline runner tests use: short reads, low coverage, deterministic seed.
pub struct GenReadsConfig {
    pub reference: PathBuf,
    pub output_dir: PathBuf,
    pub output_filename: String,
    pub read_len: usize,
    pub coverage: usize,
    pub paired_ended: bool,
    pub produce_fastq: bool,
    pub produce_bam: bool,
    pub produce_vcf: bool,
    pub rng_seed: String,
    pub sequence_error_model: Option<PathBuf>,
    pub num_threads: Option<usize>,
    pub chunk_size: Option<usize>,
    pub input_vcf: Option<PathBuf>,
    pub mutation_rate: Option<f64>,
    pub mutation_model: Option<PathBuf>,
    pub sv_rate_scale: Option<f64>,
}

impl GenReadsConfig {
    pub fn new(reference: PathBuf, output_dir: PathBuf, output_filename: &str) -> Self {
        Self {
            reference,
            output_dir,
            output_filename: output_filename.to_string(),
            read_len: 50,
            coverage: 2,
            paired_ended: false,
            produce_fastq: true,
            produce_bam: false,
            produce_vcf: false,
            rng_seed: "integration phase two".to_string(),
            sequence_error_model: None,
            num_threads: None,
            chunk_size: None,
            input_vcf: None,
            mutation_rate: None,
            mutation_model: None,
            sv_rate_scale: None,
        }
    }

    /// Serialize to a temporary YAML file (returned as a `NamedTempFile` so the caller
    /// can hand its path to the binary). The file is deleted when the value is dropped.
    pub fn write_yaml(&self) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".yml").tempfile().unwrap();
        let pair_section = if self.paired_ended {
            "paired_ended: true\nfragment_mean: 200.0\nfragment_st_dev: 30.0\n"
        } else {
            "paired_ended: false\n"
        };
        let model_section = match &self.sequence_error_model {
            Some(p) => format!("sequence_error_model: {}\n", p.display()),
            None => String::new(),
        };
        let threads_section = match self.num_threads {
            Some(n) => format!("num_threads: {n}\n"),
            None => String::new(),
        };
        let chunk_size_section = match self.chunk_size {
            Some(n) => format!("chunk_size: {n}\n"),
            None => String::new(),
        };
        let input_vcf_section = match &self.input_vcf {
            Some(p) => format!("input_vcf: {}\n", p.display()),
            None => String::new(),
        };
        let mutation_rate_section = match self.mutation_rate {
            Some(r) => format!("mutation_rate: {r}\n"),
            None => String::new(),
        };
        let mutation_model_section = match &self.mutation_model {
            Some(p) => format!("mutation_model: {}\n", p.display()),
            None => String::new(),
        };
        let sv_rate_scale_section = match self.sv_rate_scale {
            Some(s) => format!("sv_rate_scale: {s}\n"),
            None => String::new(),
        };
        write!(
            f,
            "reference: {ref_}\n\
             read_len: {rl}\n\
             coverage: {cov}\n\
             produce_fastq: {pf}\n\
             produce_bam: {pb}\n\
             produce_vcf: {pv}\n\
             output_dir: {out}\n\
             output_filename: {name}\n\
             overwrite_output: true\n\
             rng_seed: {seed}\n\
             {pair}{model}{threads}{chunksz}{ivcf}{mrate}{mmodel}{svscale}",
            ref_ = self.reference.display(),
            rl = self.read_len,
            cov = self.coverage,
            pf = self.produce_fastq,
            pb = self.produce_bam,
            pv = self.produce_vcf,
            out = self.output_dir.display(),
            name = self.output_filename,
            seed = self.rng_seed,
            pair = pair_section,
            model = model_section,
            threads = threads_section,
            chunksz = chunk_size_section,
            ivcf = input_vcf_section,
            mrate = mutation_rate_section,
            mmodel = mutation_model_section,
            svscale = sv_rate_scale_section,
        )
        .unwrap();
        f
    }
}

/// Write a synthetic FASTA containing every IUPAC ambiguity code (R/Y/M/K/S/W/H/B/V/D)
/// embedded in an otherwise valid sequence. Long enough to generate 50-bp reads.
/// Used to verify that gen-reads handles IUPAC references without crashing and does not
/// leak IUPAC characters into output FASTQ sequences.
pub fn write_iupac_fasta(path: &Path) {
    use std::io::Write as _;
    // Each IUPAC code appears at least once; padded to ~120 bp so read_len=50 fits.
    let seq = concat!(
        "ACGTRACGTYMKSWACGTHBVDACGTACGTACGT",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
    );
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">chr1\n{}", seq).unwrap();
}

/// Writes a single-contig FASTA with two ACGT flanks separated by a hard N tract.
/// Layout (1-based, inclusive): left flank [1, 240], N tract [241, 480],
/// right flank [481, 720]. Used to assert the read/variant pipeline treats the
/// N tract as an assembly gap (no anchors, no variants) rather than filling it.
pub fn write_n_tract_fasta(path: &Path) {
    use std::io::Write as _;
    let flank: String = "ACGT".chars().cycle().take(240).collect();
    let n_tract: String = "N".repeat(240);
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">chr1\n{flank}{n_tract}{flank}").unwrap();
}

/// Build a `gen-seq-error-model` config YAML with a training FASTQ and (optionally)
/// a bin set. Output model is written into `output_file`.
pub fn write_gen_seq_error_model_config(
    fastq_file: &Path,
    output_file: &Path,
    binned_quality_bins: Option<&[usize]>,
) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(".yml").tempfile().unwrap();
    let bin_line = match binned_quality_bins {
        Some(bins) => format!(
            "binned_quality_bins: [{}]\n",
            bins.iter()
                .map(|b| b.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        ),
        None => String::new(),
    };
    write!(
        f,
        "fastq_file: {fq}\n\
         output_file: {out}\n\
         overwrite_output: true\n\
         {bins}",
        fq = fastq_file.display(),
        out = output_file.display(),
        bins = bin_line,
    )
    .unwrap();
    f
}

// ── output-file inspection ───────────────────────────────────────────────────

/// Decode a gzipped FASTQ (potentially multi-stream, as eidolon emits one stream per
/// contig) into a flat vector of lines.
pub fn read_gzip_fastq_lines(path: &Path) -> Vec<String> {
    use flate2::read::MultiGzDecoder;
    use std::io::{BufRead, BufReader};
    let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
    r.lines().map(|l| l.unwrap()).collect()
}

/// Load the full bytes of a file. Used by determinism tests to compare runs.
pub fn read_file_bytes(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap()
}

// ── working-directory helper ────────────────────────────────────────────────

/// Convenience wrapper: a tempdir plus a couple of paths inside it. Returned as a
/// tuple so callers can keep the `TempDir` alive for the duration of the test.
pub fn fresh_workdir() -> (TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().to_path_buf();
    (dir, out)
}
