//! The BAM reader must agree with an independent implementation, not with my expectations.
//!
//! `metrics.rs` is provable from literals; this layer is not. It depends on noodles decoding
//! flags, CIGARs and positions the way I believe it does, and the only honest check is to ask
//! a different tool the same question. `samtools` is that tool.
//!
//! `#[ignore]`d because it needs samtools and a generated BAM — the realigned-gates CI job
//! runs `--ignored` and has samtools installed. Run locally with:
//!
//! ```text
//! cargo test --test realism_reader -- --ignored --nocapture
//! ```

mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::path::Path;
use std::process::Command;

use eidolon_core::realism::reader::{Region, measure};

fn have(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Generate a golden BAM to measure. H1N1 is fine here: this test is about whether the reader
/// decodes a BAM correctly, which is not a question about realism or scale.
fn golden_bam(dir: &Path, name: &str) -> std::path::PathBuf {
    let mut config = GenReadsConfig::new(h1n1_reference(), dir.to_path_buf(), name);
    config.coverage = 30;
    config.produce_bam = true;
    config.produce_fastq = false;
    config.rng_seed = "realism reader".to_string();
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    dir.join(format!("{name}.bam"))
}

/// samtools' answer to the same three questions, over the same records.
///
/// `-F 0x904` excludes unmapped, secondary and supplementary alignments — the same set the
/// reader drops. A supplementary alignment is the other half of a split read, so counting it
/// would double-count the clip boundaries this panel is built to measure.
///
/// The region is filtered HERE rather than passed to `samtools view` as a region argument.
/// A region query needs a `.bai`, and a freshly written golden BAM has none — samtools then
/// returns zero records rather than an error, so the comparison would have been "0 vs 0" and
/// passed while comparing nothing. That is exactly the failure this test exists to catch, and
/// it caught it on the first run.
fn samtools_counts(bam: &Path, contig: &str, start: usize, end: usize) -> (usize, usize, usize) {
    let out = Command::new("samtools")
        .args(["view", "-F", "0x904", bam.to_str().unwrap()])
        .output()
        .expect("samtools view failed");
    let text = String::from_utf8_lossy(&out.stdout);
    let mut reads = 0;
    let mut improper = 0;
    let mut clipped = 0;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            continue;
        }
        if f[2] != contig {
            continue;
        }
        let pos0: usize = f[3].parse::<usize>().unwrap() - 1;
        if pos0 < start || pos0 >= end {
            continue;
        }
        reads += 1;
        let flag: u16 = f[1].parse().unwrap();
        if flag & 0x2 == 0 {
            improper += 1;
        }
        // Any soft clip of >= 20 bp, at either end.
        let cig = f[5];
        let mut n = String::new();
        let mut ops: Vec<(char, usize)> = Vec::new();
        for c in cig.chars() {
            if c.is_ascii_digit() {
                n.push(c);
            } else {
                ops.push((c, n.parse().unwrap_or(0)));
                n.clear();
            }
        }
        let lead = matches!(ops.first(), Some(('S', k)) if *k >= 20);
        let trail = matches!(ops.last(), Some(('S', k)) if *k >= 20);
        if lead || trail {
            clipped += 1;
        }
    }
    (reads, improper, clipped)
}

#[test]
#[ignore = "requires samtools; run with --ignored (see module docs)"]
fn the_reader_agrees_with_samtools_on_the_same_records() {
    if !have("samtools") {
        panic!(
            "samtools not on PATH — this test asserts agreement with it and cannot be skipped silently"
        );
    }
    let (_g, work) = fresh_workdir();
    let bam = golden_bam(&work, "reader_check");

    // H1N1_HA is the fixture's first contig; 1..1500 sits inside it.
    let region = Region {
        contig: "H1N1_HA".into(),
        start: 0,
        end: 1500,
    };
    let mine = measure(&bam, std::slice::from_ref(&region), 20, 3, 2000, 500)
        .expect("reader failed on a BAM it just generated");
    let m = &mine[0];

    let (reads, improper, clipped) = samtools_counts(&bam, "H1N1_HA", 0, 1500);

    assert!(
        reads > 0,
        "samtools reported no reads at all — the comparison would be vacuous"
    );
    assert_eq!(
        m.reads, reads,
        "read count disagrees with samtools: mine {} vs {reads}",
        m.reads
    );
    assert_eq!(
        m.improper_pairs, improper,
        "improper-pair count disagrees with samtools: mine {} vs {improper}",
        m.improper_pairs
    );
    assert_eq!(
        m.clipped_reads, clipped,
        "clipped-read count disagrees with samtools: mine {} vs {clipped}",
        m.clipped_reads
    );
}

/// An unmeasurable region is an ERROR, never a zero.
///
/// This is rule 4 in its most literal form: a region reported as "0 artifacts" because it
/// matched no reads is indistinguishable in the output from genuinely clean data — which is
/// the exact confusion the whole panel exists to prevent.
#[test]
#[ignore = "requires samtools; run with --ignored (see module docs)"]
fn an_unmeasurable_region_is_an_error_not_a_clean_result() {
    let (_g, work) = fresh_workdir();
    let bam = golden_bam(&work, "reader_errors");

    let missing = Region {
        contig: "NOT_A_CONTIG".into(),
        start: 0,
        end: 1000,
    };
    let err = measure(&bam, &[missing], 20, 3, 2000, 500).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("not in the BAM header"), "unexpected: {msg}");

    // A real contig, but far past its end: no reads, and that must not read as clean.
    let empty = Region {
        contig: "H1N1_HA".into(),
        start: 900_000,
        end: 901_000,
    };
    let err = measure(&bam, &[empty], 20, 3, 2000, 500).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("contained no reads"), "unexpected: {msg}");
}

/// The same BAM measured twice must give the same answer — the real-vs-real self-test at the
/// file level. If the reader disagreed with itself, every gap the panel reported would be its
/// own nondeterminism.
#[test]
#[ignore = "requires samtools; run with --ignored (see module docs)"]
fn measuring_one_file_twice_gives_one_answer() {
    let (_g, work) = fresh_workdir();
    let bam = golden_bam(&work, "reader_stable");
    let region = Region {
        contig: "H1N1_HA".into(),
        start: 0,
        end: 1700,
    };
    let a = measure(&bam, std::slice::from_ref(&region), 20, 3, 2000, 500).unwrap();
    let b = measure(&bam, std::slice::from_ref(&region), 20, 3, 2000, 500).unwrap();
    assert_eq!(a, b);
}
