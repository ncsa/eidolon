//! A tandem duplication changes depth only where it is (#884).
//!
//! **Criterion.** For a duplication of any size from 100 bp up, depth in reference coordinates,
//! measured against a no-variant run of the same configuration, is the declared multiplier
//! (1.5 heterozygous, 2.0 homozygous) across the duplicated span and 1.0 outside it, including
//! the first fragment length beyond each edge. **Falsified by** a flank window outside
//! 0.94 to 1.06 of the control, or a span mean off the multiplier by more than the tolerance
//! below, on any size.
//!
//! Measured before this test existed (ecoli, 60x, 500 bp heterozygous duplications): the 100 to
//! 300 bp past the end read 1.20 to 1.33 of control. Fragments started at the raised density
//! inside the span run on into the flank, and the explicit junction reads come on top of them.
//!
//! The reference is synthetic and the comparison is against a separate no-variant run, never a
//! flank in the same run. CI has no samtools, so the BAM is read with noodles.

mod common;

use common::eidolon;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::{Path, PathBuf};

const CONTIG: &str = "synth1";
const CONTIG_LEN: usize = 700_000;
const FLANK: usize = 250;

/// (size, genotype, multiplier, label)
const KINDS: [(usize, &str, f64, &str); 5] = [
    (100, "0/1", 1.5, "het 100 bp"),
    (500, "0/1", 1.5, "het 500 bp"),
    (1_000, "0/1", 1.5, "het 1 kb"),
    (5_000, "0/1", 1.5, "het 5 kb"),
    (500, "1/1", 2.0, "hom 500 bp"),
];
const REPLICATES: usize = 12;

fn synthetic_sequence(n: usize) -> String {
    let mut x: u64 = 0x2468_ACE0_1357_9BDF;
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

/// Every planted event as 0-based half-open `(start, end)` with its kind index, spaced 10 kb apart.
fn events() -> Vec<(usize, usize, usize)> {
    let mut out = Vec::new();
    let mut at = 30_000;
    for (k, &(size, ..)) in KINDS.iter().enumerate() {
        for _ in 0..REPLICATES {
            out.push((at, at + size, k));
            at += 10_000;
        }
    }
    out
}

fn run(dir: &Path, name: &str, with_events: bool) -> PathBuf {
    let out = dir.join(name);
    std::fs::create_dir_all(&out).unwrap();
    let mut fa = format!(">{CONTIG}\n");
    for chunk in synthetic_sequence(CONTIG_LEN).as_bytes().chunks(60) {
        fa.push_str(std::str::from_utf8(chunk).unwrap());
        fa.push('\n');
    }
    let reference = out.join("synthetic.fa");
    std::fs::write(&reference, fa).unwrap();

    let vcf = out.join("in.vcf");
    let mut text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={CONTIG_LEN}>\n\
         ##ALT=<ID=DUP,Description=\"dup\">\n\
         ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n\
         ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"x\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    );
    if with_events {
        for (s, e, k) in events() {
            let (size, gt, ..) = KINDS[k];
            text.push_str(&format!(
                "{CONTIG}\t{s}\t.\tG\t<DUP>\t50\tPASS\tSVTYPE=DUP;END={e};SVLEN={size}\tGT\t{gt}\n"
            ));
        }
    }
    std::fs::write(&vcf, text).unwrap();

    let cfg = out.join("c.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 100\ncoverage: 60\nploidy: 2\npaired_ended: true\n\
             fragment_mean: 250\nfragment_st_dev: 25\ninput_vcf: {vcf}\n\
             produce_vcf: false\nproduce_fastq: false\nproduce_bam: true\n\
             sv_rate_scale: 0.0\nmutation_rate: 0.0\noverwrite_output: true\n\
             output_dir: {out}\noutput_basename: o\nrng_seed: dup depth profile\nnum_threads: 4\n",
            ref = reference.display(),
            vcf = vcf.display(),
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
        // The golden BAM writes a duplication's junction read as one continuous alignment from
        // the first piece, so its bases past END lie over the flank on paper. Its sequence is the
        // end of the span followed by the start of the span, so each such base belongs at
        // POS + (x - END). Undoing that is the stitching rule, not a correction to the result.
        let wrap = rec.name().and_then(|n| junction_span(&n.to_string()));
        for x in start..(start + span).min(CONTIG_LEN) {
            let at = match wrap {
                Some((s, e)) if x >= e => s + (x - e),
                _ => x,
            };
            if at < CONTIG_LEN {
                depth[at] += 1;
            }
        }
    }
    depth
}

/// `(POS, END)` from `EIDOLON_chimeric_DUP_<contig>_<POS>_<END>_<index>/<mate>`.
fn junction_span(name: &str) -> Option<(usize, usize)> {
    if !name.starts_with("EIDOLON_chimeric_DUP_") {
        return None;
    }
    let mut fields = name.rsplit('_');
    let _index = fields.next()?;
    let end = fields.next()?.parse().ok()?;
    let pos = fields.next()?.parse().ok()?;
    Some((pos, end))
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
fn a_duplication_changes_depth_only_where_it_is() {
    let tmp = tempfile::tempdir().unwrap();
    let control = depth_track(&run(tmp.path(), "control", false));
    let events_depth = depth_track(&run(tmp.path(), "events", true));
    let all = events();
    let mut failures = Vec::new();

    for (k, &(size, _, mult, label)) in KINDS.iter().enumerate() {
        let mine: Vec<_> = all.iter().filter(|e| e.2 == k).collect();
        // The denominator: every planted event of this kind must be scored.
        assert_eq!(mine.len(), REPLICATES, "{label}: events planted");

        let span: Vec<_> = mine.iter().map(|&&(s, e, _)| (s, e)).collect();
        let right: Vec<_> = mine.iter().map(|&&(_, e, _)| (e, e + FLANK)).collect();
        let left: Vec<_> = mine.iter().map(|&&(s, ..)| (s - FLANK, s)).collect();
        let ratio = |r: &[(usize, usize)]| mean_over(&events_depth, r) / mean_over(&control, r);
        let (span_r, right_r, left_r) = (ratio(&span), ratio(&right), ratio(&left));
        eprintln!(
            "{label:>11}: span {span_r:.2} (want {mult})  left flank {left_r:.2}  right flank {right_r:.2}"
        );

        let tol = if size < 250 { 0.15 } else { 0.08 };
        if (span_r - mult).abs() > tol * mult {
            failures.push(format!("{label}: span ratio {span_r:.2}, want {mult}"));
        }
        for (name, r) in [("left", left_r), ("right", right_r)] {
            if !(0.92..=1.08).contains(&r) {
                failures.push(format!("{label}: {name} flank {r:.2} of control, want 1.0"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
