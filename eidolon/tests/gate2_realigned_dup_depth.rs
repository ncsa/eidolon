//! **Gate 2 for DUP depth** — after a real aligner, does a short duplication change depth only
//! where it is (#884)?
//!
//! **Criterion.** Realigned with bwa-mem2 and compared with a no-variant run of the same
//! configuration, depth over a heterozygous duplication of 50 to 500 bp is 1.5x across the span
//! and 1.0x in the 250 bp on either side. **Falsified by** a span more than 0.12 off 1.5, or a
//! flank outside 0.95 to 1.07, on any size.
//!
//! This is the small-event half of `dup_depth_profile.rs`, which reads the golden BAM. The
//! golden BAM writes a junction read as one continuous alignment, which is wrong on paper for a
//! span shorter than a fragment, so for those sizes only a realigned read says where the bases
//! go. Measured when this was written (ecoli, 60x, 2 x 150 bp, 350 bp fragments, 40 events per
//! size), before the flank fix: 100 bp heterozygous read **1.18 to 1.26** over the span and **1.10
//! to 1.14** over the 300 bp past its end. A fragment that starts before a short span and runs
//! through it continues into the second copy, and it was drawn as a reference fragment that
//! carried on into the flank.
//!
//! Requires bwa-mem2; `#[ignore]`d so CI never reaches it, and run by `release-gates.yml`.
//!
//! ```text
//! conda activate aln
//! cargo test --test gate2_realigned_dup_depth -- --ignored --nocapture
//! ```

mod common;

use common::eidolon;
use common::gate2::bwa_mem2;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

const CONTIG: &str = "synth2";
const CONTIG_LEN: usize = 560_000;
const FLANK: usize = 250;
const REPLICATES: usize = 12;
/// (size, label). All heterozygous: the 500 bp row is the must-not-fire case, a span longer than
/// a fragment, which has nothing to divert.
const KINDS: [(usize, &str); 4] = [
    (50, "het 50 bp"),
    (100, "het 100 bp"),
    (200, "het 200 bp"),
    (500, "het 500 bp"),
];

fn synthetic_sequence(n: usize) -> String {
    let mut x: u64 = 0x1357_9BDF_2468_ACE0;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            match (x >> 33) & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        })
        .collect()
}

/// 0-based half-open `(start, end)` and kind index, spaced 10 kb apart.
fn events() -> Vec<(usize, usize, usize)> {
    let mut out = Vec::new();
    let mut at = 30_000;
    for (k, &(size, _)) in KINDS.iter().enumerate() {
        for _ in 0..REPLICATES {
            out.push((at, at + size, k));
            at += 10_000;
        }
    }
    assert!(at < CONTIG_LEN, "events run past the contig");
    out
}

/// Generates reads, aligns them, and returns per-base depth at MAPQ 20 and up (what
/// `samtools depth -Q 20` counts).
fn realigned_depth(dir: &Path, name: &str, with_events: bool, bwa: &str) -> Vec<u32> {
    let out = dir.join(name);
    std::fs::create_dir_all(&out).unwrap();
    let mut fa = format!(">{CONTIG}\n");
    for chunk in synthetic_sequence(CONTIG_LEN).as_bytes().chunks(60) {
        fa.push_str(std::str::from_utf8(chunk).unwrap());
        fa.push('\n');
    }
    let reference = out.join("synthetic.fa");
    std::fs::write(&reference, fa).unwrap();

    let mut vcf = format!(
        "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={CONTIG_LEN}>\n\
         ##ALT=<ID=DUP,Description=\"dup\">\n\
         ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n\
         ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"x\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    );
    if with_events {
        for (s, e, k) in events() {
            vcf.push_str(&format!(
                "{CONTIG}\t{s}\t.\tG\t<DUP>\t50\tPASS\tSVTYPE=DUP;END={e};SVLEN={}\tGT\t0/1\n",
                KINDS[k].0
            ));
        }
    }
    let vcf_path = out.join("in.vcf");
    std::fs::write(&vcf_path, vcf).unwrap();

    let cfg = out.join("c.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 150\ncoverage: 60\nploidy: 2\npaired_ended: true\n\
             fragment_mean: 350\nfragment_st_dev: 35\ninput_vcf: {vcf}\n\
             produce_vcf: false\nproduce_fastq: true\nproduce_bam: false\nshuffle_fastq: false\n\
             sv_rate_scale: 0.0\nmutation_rate: 0.0\noverwrite_output: true\n\
             output_dir: {out}\noutput_basename: o\nrng_seed: dup depth gate\nnum_threads: 4\n",
            ref = reference.display(),
            vcf = vcf_path.display(),
            out = out.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();

    let index = Command::new(bwa)
        .arg("index")
        .arg(&reference)
        .output()
        .expect("bwa-mem2 index failed to spawn");
    assert!(index.status.success(), "bwa-mem2 index failed");
    let r1: PathBuf = out.join("o_r1.fastq.gz");
    let r2: PathBuf = out.join("o_r2.fastq.gz");
    assert!(r1.is_file() && r2.is_file(), "expected {r1:?} and {r2:?}");
    let sam = Command::new(bwa)
        .args(["mem", "-t", "4"])
        .arg(&reference)
        .arg(&r1)
        .arg(&r2)
        .output()
        .expect("bwa-mem2 mem failed to spawn");
    assert!(sam.status.success(), "bwa-mem2 mem failed");

    let mut depth = vec![0u32; CONTIG_LEN];
    let mut placed = 0usize;
    for line in BufReader::new(sam.stdout.as_slice()).lines() {
        let line = line.unwrap();
        if line.starts_with('@') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let flag: u32 = f[1].parse().unwrap();
        let mapq: u32 = f[4].parse().unwrap();
        if flag & (0x4 | 0x100 | 0x200 | 0x400) != 0 || mapq < 20 {
            continue;
        }
        placed += 1;
        let mut at: usize = f[3].parse::<usize>().unwrap() - 1;
        let mut num = 0usize;
        for c in f[5].chars() {
            if let Some(d) = c.to_digit(10) {
                num = num * 10 + d as usize;
                continue;
            }
            match c {
                'M' | '=' | 'X' => {
                    for x in at..(at + num).min(CONTIG_LEN) {
                        depth[x] += 1;
                    }
                    at += num;
                }
                'D' | 'N' => at += num,
                _ => {}
            }
            num = 0;
        }
    }
    assert!(placed > 10_000, "{name}: only {placed} reads placed");
    depth
}

fn mean_over(depth: &[u32], regions: &[(usize, usize)]) -> f64 {
    let (sum, n) = regions.iter().fold((0u64, 0u64), |(s, n), &(a, b)| {
        (
            s + depth[a..b].iter().map(|&d| d as u64).sum::<u64>(),
            n + (b - a) as u64,
        )
    });
    sum as f64 / n as f64
}

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_a_short_duplication_changes_realigned_depth_only_where_it_is() {
    let bwa = bwa_mem2();
    let tmp = tempfile::tempdir().unwrap();
    let control = realigned_depth(tmp.path(), "control", false, &bwa);
    let with = realigned_depth(tmp.path(), "events", true, &bwa);
    let all = events();
    let mut failures = Vec::new();

    for (k, &(_, label)) in KINDS.iter().enumerate() {
        let mine: Vec<_> = all.iter().filter(|e| e.2 == k).collect();
        assert_eq!(mine.len(), REPLICATES, "{label}: events planted");
        let span: Vec<_> = mine.iter().map(|&&(s, e, _)| (s, e)).collect();
        let right: Vec<_> = mine.iter().map(|&&(_, e, _)| (e, e + FLANK)).collect();
        let left: Vec<_> = mine.iter().map(|&&(s, ..)| (s - FLANK, s)).collect();
        let ratio = |r: &[(usize, usize)]| mean_over(&with, r) / mean_over(&control, r);
        let (span_r, right_r, left_r) = (ratio(&span), ratio(&right), ratio(&left));
        println!(
            "{label:>11}: span {span_r:.2} (want 1.5)  left flank {left_r:.2}  right flank {right_r:.2}"
        );
        if (span_r - 1.5).abs() > 0.12 {
            failures.push(format!("{label}: span {span_r:.2}, want 1.5"));
        }
        if !(0.95..=1.07).contains(&right_r) {
            failures.push(format!("{label}: right flank {right_r:.2}, want 1.0"));
        }
        if !(0.95..=1.05).contains(&left_r) {
            failures.push(format!("{label}: left flank {left_r:.2}, want 1.0"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
