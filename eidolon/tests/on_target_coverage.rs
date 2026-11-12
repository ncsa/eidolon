//! `coverage` on a `target_bed` run is the mean depth over the target bases (#578).
//!
//! The depth is measured here from the golden BAM, by counting read bases inside the BED. That
//! is independent of the estimate the runner uses to scale its fragment budget, so these tests
//! can disagree with it. CI has no samtools, so the BAM is read with noodles.

mod common;

use common::eidolon;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::{Path, PathBuf};

const CONTIG: &str = "synth1";
const CONTIG_LEN: usize = 1_000_000;

/// Deterministic pseudo-random ACGT, so the fixture needs no data file.
fn synthetic_sequence(n: usize) -> String {
    let mut x: u64 = 0x1234_5678_9ABC_DEF0;
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

fn write_reference(dir: &Path) -> PathBuf {
    let path = dir.join("synthetic.fa");
    let mut out = format!(">{CONTIG}\n");
    for chunk in synthetic_sequence(CONTIG_LEN).as_bytes().chunks(60) {
        out.push_str(std::str::from_utf8(chunk).unwrap());
        out.push('\n');
    }
    std::fs::write(&path, out).unwrap();
    path
}

fn write_bed(dir: &Path, targets: &[(usize, usize)]) -> PathBuf {
    let path = dir.join("targets.bed");
    let text: String = targets
        .iter()
        .map(|(s, e)| format!("{CONTIG}\t{s}\t{e}\n"))
        .collect();
    std::fs::write(&path, text).unwrap();
    path
}

/// Runs gen-reads over the synthetic reference at `coverage: 30`, with the targets and any extra
/// VCF records given, and returns the output directory.
fn run(dir: &Path, name: &str, targets: &[(usize, usize)], vcf_records: &[String]) -> PathBuf {
    let out = dir.join(name);
    std::fs::create_dir_all(&out).unwrap();
    let reference = write_reference(&out);
    let bed = write_bed(&out, targets);
    let mut extra = String::new();
    if !vcf_records.is_empty() {
        let vcf = out.join("in.vcf");
        let mut text = format!(
            "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={CONTIG_LEN}>\n\
             ##ALT=<ID=DUP,Description=\"dup\">\n\
             ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n\
             ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"x\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        );
        for r in vcf_records {
            text.push_str(r);
            text.push('\n');
        }
        std::fs::write(&vcf, text).unwrap();
        extra = format!("input_vcf: {}\n", vcf.display());
    }
    let cfg = out.join("c.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\ntarget_bed: {bed}\nread_len: 100\ncoverage: 30\nploidy: 2\n\
             paired_ended: true\nfragment_mean: 250\nfragment_st_dev: 25\n{extra}\
             produce_vcf: false\nproduce_fastq: false\nproduce_bam: true\n\
             sv_rate_scale: 0.0\noverwrite_output: true\n\
             output_dir: {out}\noutput_basename: o\nrng_seed: on target {name}\nnum_threads: 4\n",
            ref = reference.display(),
            bed = bed.display(),
            out = out.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    out
}

/// Per-base depth over the whole contig, 0-based, from the BAM's reference spans.
fn depth_track(out: &Path) -> Vec<u32> {
    let file = std::fs::File::open(out.join("o.bam")).expect("bam produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let mut depth = vec![0u32; CONTIG_LEN];
    for rec in reader.records() {
        let rec = rec.unwrap();
        let Some(Ok(pos)) = rec.alignment_start() else {
            continue;
        };
        let span: usize = rec
            .cigar()
            .iter()
            .map(|op| op.unwrap())
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::Deletion
                        | Kind::Skip
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                )
            })
            .map(|op| op.len())
            .sum();
        let start = pos.get() - 1;
        for d in &mut depth[start..(start + span).min(CONTIG_LEN)] {
            *d += 1;
        }
    }
    depth
}

/// Mean depth over the union of `targets` (0-based half-open).
fn mean_on_target(depth: &[u32], targets: &[(usize, usize)]) -> f64 {
    let (mut sum, mut bases) = (0u64, 0u64);
    for &(s, e) in targets {
        sum += depth[s..e].iter().map(|&d| d as u64).sum::<u64>();
        bases += (e - s) as u64;
    }
    sum as f64 / bases as f64
}

fn mean_over(depth: &[u32], s: usize, e: usize) -> f64 {
    depth[s..e].iter().map(|&d| d as f64).sum::<f64>() / (e - s) as f64
}

/// 200 isolated 180 bp targets, 4 kb apart: the exome shape. 100 bp reads on ~250 bp fragments
/// put most read bases off target, so before #578 this delivered roughly 30 x 0.4 = 12x.
#[test]
fn isolated_targets_deliver_the_requested_on_target_depth() {
    let tmp = tempfile::tempdir().unwrap();
    let targets: Vec<(usize, usize)> = (0..200)
        .map(|i| (50_000 + i * 4_000, 50_000 + i * 4_000 + 180))
        .collect();
    let out = run(tmp.path(), "isolated", &targets, &[]);
    let got = mean_on_target(&depth_track(&out), &targets);
    eprintln!("[on-target] isolated 180 bp targets, coverage 30: delivered {got:.2}x");
    assert!(
        (got / 30.0 - 1.0).abs() < 0.05,
        "on-target depth {got:.2}x against a requested 30x"
    );
}

/// Targets narrower than a read: 250 isolated 80 bp targets against 100 bp reads. A real exome
/// is full of these (63% of one public IDT design's chr22 targets are 150 bp or less). Before
/// the fix they owned no fragments and were covered by spill alone, about 1x against 30x.
#[test]
fn targets_narrower_than_a_read_still_get_their_depth() {
    let tmp = tempfile::tempdir().unwrap();
    let targets: Vec<(usize, usize)> = (0..250)
        .map(|i| (50_000 + i * 3_000, 50_000 + i * 3_000 + 80))
        .collect();
    let out = run(tmp.path(), "narrow", &targets, &[]);
    let got = mean_on_target(&depth_track(&out), &targets);
    eprintln!("[on-target] 80 bp targets, 100 bp reads, coverage 30: delivered {got:.2}x");
    assert!(
        (got / 30.0 - 1.0).abs() < 0.07,
        "on-target depth {got:.2}x against a requested 30x"
    );
}

/// MUST NOT FIRE: 400 touching 1 kb intervals tile 300 kb..700 kb. Spill from each lands on the
/// next, so nothing is lost and the budget must not be inflated. The isolated-target formula
/// applied here would over-deliver by about a quarter.
#[test]
fn touching_targets_are_not_inflated() {
    let tmp = tempfile::tempdir().unwrap();
    let targets: Vec<(usize, usize)> = (0..400)
        .map(|i| (300_000 + i * 1_000, 300_000 + (i + 1) * 1_000))
        .collect();
    let out = run(tmp.path(), "touching", &targets, &[]);
    // Score away from the window's two ends, where there is no neighbor.
    let depth = depth_track(&out);
    let got = mean_over(&depth, 310_000, 690_000);
    eprintln!("[on-target] tiled window, coverage 30: delivered {got:.2}x");
    assert!(
        (got / 30.0 - 1.0).abs() < 0.05,
        "tiled-window depth {got:.2}x against a requested 30x"
    );
}

/// MUST NOT FIRE: a 100 kb het duplication (declared 1.5x) inside one wide target. The
/// correction is one factor for the whole BED, so the ratio between the event and its flank is
/// the declared multiplier, as it is without a BED.
#[test]
fn sv_depth_modulation_inside_a_target_is_unchanged() {
    let tmp = tempfile::tempdir().unwrap();
    let (anchor, svlen) = (400_000usize, 100_000usize);
    let end = anchor + svlen;
    let record = format!(
        "{CONTIG}\t{anchor}\tv\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END={end};SVLEN={svlen}\tGT\t0/1"
    );
    let targets = [(300_000, 700_000)];
    let out = run(tmp.path(), "sv", &targets, &[record]);
    let depth = depth_track(&out);
    let inside = mean_over(&depth, anchor + 1_000, end - 1_000);
    let flank = mean_over(&depth, 310_000, anchor - 1_000);
    eprintln!("[on-target] DUP inside a wide target: event {inside:.2}x, flank {flank:.2}x");
    assert!(
        (flank / 30.0 - 1.0).abs() < 0.05,
        "flank depth {flank:.2}x against a requested 30x"
    );
    assert!(
        (inside / flank / 1.5 - 1.0).abs() < 0.10,
        "event/flank ratio {:.3}, declared 1.5",
        inside / flank
    );
}
