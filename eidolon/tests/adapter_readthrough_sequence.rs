//! What a short-insert read carries after its insert (#865, the same defect as NEAT#344).
//!
//! Adapter bases used to be sourced cyclically, so a read whose tail was longer than the
//! adapter looped back into a second copy of it: a ~40 bp insert in a 150 bp read carried
//! the 33 bp TruSeq stem three times. No molecule has that sequence. Now a read runs
//! through the whole construct once (stem, index, P7/P5; 64 and 58 bp for TruSeq) and is
//! then filled per `adapters.overrun`.
//!
//! **Known answer.** With 150 bp reads and inserts of 40 +/- 5 bp, every read has at
//! least ~150 - 65 - 64 = 21 bases past the construct, so with `overrun: polyg` its last
//! 20 bases are G exactly (fill takes no substitution errors). The 13-mer stem
//! `AGATCGGAAGAGC` occurs at most once per read; sequencing error can only remove an
//! occurrence, never add one. **Must not fire:** with `overrun: random` the same reads
//! must not end in a G run.

mod common;

use common::{eidolon, h1n1_reference};
use std::path::Path;

const STEM: &str = "AGATCGGAAGAGC";

fn read_seqs(path: &Path) -> Vec<String> {
    use std::io::Read;
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s.lines().skip(1).step_by(4).map(str::to_string).collect()
}

/// R1 and R2 sequences from a short-insert TruSeq run with the given overrun mode.
fn run(dir: &Path, overrun: &str) -> Vec<String> {
    let cfg = dir.join(format!("{overrun}.yml"));
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 150\ncoverage: 30\nploidy: 1\npaired_ended: true\n\
             fragment_mean: 40\nfragment_st_dev: 5\n\
             produce_bam: false\nproduce_fastq: true\nproduce_vcf: false\n\
             overwrite_output: true\noutput_dir: {out}\noutput_filename: {overrun}\n\
             rng_seed: adapter readthrough\nnum_threads: 1\n\
             adapters:\n  enabled: true\n  preset: truseq\n  overrun: {overrun}\n",
            ref = h1n1_reference().display(),
            out = dir.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    let mut seqs = read_seqs(&dir.join(format!("{overrun}_r1.fastq.gz")));
    seqs.extend(read_seqs(&dir.join(format!("{overrun}_r2.fastq.gz"))));
    seqs
}

#[test]
fn short_insert_reads_run_through_the_adapter_once_then_poly_g() {
    let tmp = tempfile::tempdir().unwrap();
    let seqs = run(tmp.path(), "polyg");
    assert!(
        seqs.len() > 1000,
        "only {} reads; the test would be weak",
        seqs.len()
    );

    let repeated = seqs.iter().filter(|s| s.matches(STEM).count() > 1).count();
    assert_eq!(
        repeated,
        0,
        "{repeated} of {} reads carry the adapter stem more than once, e.g. {}",
        seqs.len(),
        seqs.iter()
            .find(|s| s.matches(STEM).count() > 1)
            .unwrap_or(&String::new())
    );
    // Not vacuous: the stem is there to find in nearly every read.
    let with_stem = seqs.iter().filter(|s| s.contains(STEM)).count();
    assert!(
        with_stem as f64 > 0.9 * seqs.len() as f64,
        "only {with_stem} of {} reads contain the stem",
        seqs.len()
    );

    let g20 = "G".repeat(20);
    let not_g = seqs.iter().filter(|s| !s.ends_with(&g20)).count();
    assert_eq!(
        not_g,
        0,
        "{not_g} of {} reads do not end in 20 G, e.g. {}",
        seqs.len(),
        seqs.iter()
            .find(|s| !s.ends_with(&g20))
            .unwrap_or(&String::new())
    );
}

#[test]
fn random_overrun_does_not_end_in_poly_g() {
    let tmp = tempfile::tempdir().unwrap();
    let seqs = run(tmp.path(), "random");
    assert!(seqs.len() > 1000, "only {} reads", seqs.len());
    let g20 = "G".repeat(20);
    let polyg = seqs.iter().filter(|s| s.ends_with(&g20)).count();
    assert_eq!(polyg, 0, "{polyg} reads end in 20 G with overrun: random");
    let repeated = seqs.iter().filter(|s| s.matches(STEM).count() > 1).count();
    assert_eq!(
        repeated, 0,
        "{repeated} reads repeat the stem with overrun: random"
    );
}
