//! A tandem duplication changes depth only where it is (#884).
//!
//! **Criterion.** For a duplication of any size from 50 bp up, depth in reference coordinates,
//! measured against a no-variant run of the same configuration, is the declared multiplier
//! (1.5 heterozygous, 2.0 homozygous) across the duplicated span and 1.0 outside it, including
//! the first fragment length beyond each edge, and every carrier has the same number of
//! fragments across its junction whatever the size of the event. **Falsified by** a flank window
//! outside 0.92 to 1.08 of the control, a junction-pair count more than 15% off the carriers'
//! fragment density times the mean fragment length, or (from 250 bp up, see below) a span mean
//! off the multiplier by more than 8%, on any size.
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
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const CONTIG: &str = "synth1";
const CONTIG_LEN: usize = 800_000;
const FLANK: usize = 250;

/// (size, genotype, multiplier, label)
const KINDS: [(usize, &str, f64, &str); 6] = [
    (50, "0/1", 1.5, "het 50 bp"),
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
            "reference: {ref}\nread_len: 150\ncoverage: 60\nploidy: 2\npaired_ended: true\n\
             fragment_mean: 350\nfragment_st_dev: 35\ninput_vcf: {vcf}\n\
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

/// Depth per reference base, and how many junction pairs each event `(POS, END)` produced.
fn depth_track(out: &Path) -> (Vec<u32>, HashMap<(usize, usize), usize>) {
    let file = std::fs::File::open(out.join("o.bam")).expect("bam produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let mut depth = vec![0u32; CONTIG_LEN];
    let mut junction_pairs: HashMap<(usize, usize), usize> = HashMap::new();
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
        // the first piece, so mate 1's bases past END lie over the flank on paper. Its sequence is
        // the end of the span followed by the start of the span, so each such base belongs at
        // POS + (x - END). Undoing that is the stitching rule, not a correction to the result.
        // Mate 2 is placed over the second piece and needs nothing, until the span is shorter
        // than a fragment: then its first bases belong at the end of the span and it is not
        // placed there, which is why spans under 250 bp are not asserted here (below).
        let name = rec.name().map(|n| n.to_string());
        let first = rec.flags().is_first_segment();
        let wrap = name.as_deref().and_then(junction_span);
        if let Some(w) = wrap.filter(|_| first) {
            *junction_pairs.entry(w).or_default() += 1;
        }
        for x in start..(start + span).min(CONTIG_LEN) {
            let at = match wrap {
                Some((s, e)) if first && x >= e => s + (x - e),
                _ => x,
            };
            if at < CONTIG_LEN {
                depth[at] += 1;
            }
        }
    }
    (depth, junction_pairs)
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
    let (control, _) = depth_track(&run(tmp.path(), "control", false));
    let (events_depth, junction_pairs) = depth_track(&run(tmp.path(), "events", true));
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

        // A span shorter than a fragment is measured over realigned reads instead, where the
        // golden BAM's placement of mate 2 does not matter: `gate2_realigned_dup_depth.rs`.
        if size >= 250 && (span_r - mult).abs() > 0.08 * mult {
            failures.push(format!("{label}: span ratio {span_r:.2}, want {mult}"));
        }
        // Junction evidence does not depend on the size of the event. Every molecule that carries
        // the duplication has the same number of fragments across its junction: the carriers'
        // fragment density times the mean fragment length, whether those fragments start inside
        // the span or before it. 60x, 2 x 150 bp, 350 bp fragments: 350 x 60 / 300 = 70 for a
        // homozygous event and half that for a heterozygous one.
        let share = (mult - 1.0).min(1.0);
        let want_pairs = 350.0 * 60.0 * share / 300.0;
        let pairs: Vec<f64> = mine
            .iter()
            .map(|&&(s, e, _)| *junction_pairs.get(&(s, e)).unwrap_or(&0) as f64)
            .collect();
        let got_pairs = pairs.iter().sum::<f64>() / pairs.len() as f64;
        eprintln!("{label:>11}: junction pairs per event {got_pairs:.1} (want {want_pairs:.1})");
        if (got_pairs - want_pairs).abs() > 0.15 * want_pairs {
            failures.push(format!(
                "{label}: {got_pairs:.1} junction pairs per event, want {want_pairs:.1}"
            ));
        }
        for (name, r) in [("left", left_r), ("right", right_r)] {
            if !(0.92..=1.08).contains(&r) {
                failures.push(format!("{label}: {name} flank {r:.2} of control, want 1.0"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
