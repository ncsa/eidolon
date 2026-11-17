//! A gap narrower than a read between two SV events gets the depth its neighbors supply, no more
//! (#881).
//!
//! The narrow-region rule exists for `target_bed` targets, which own fragments even when
//! narrower than a read. Applied to the short pieces an SV event carves out of a contig it
//! over-delivers: the neighbors' reads already cover a 60 bp gap, and a fragment budget of its
//! own on top read 12.5% over a no-variant control here (8% on ecoli at 100 bp). The reference is
//! synthetic and the comparison is against the same run with no variants. The ratio is 1.04 with
//! the rule restricted to targets; the band is 0.97 to 1.07.
//!
//! CI has no samtools, so the BAM is read with noodles.

mod common;

use common::eidolon;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::{Path, PathBuf};

const CONTIG: &str = "synth1";
const CONTIG_LEN: usize = 1_000_000;
const GAP: usize = 60;
const DEL_LEN: usize = 300;

fn synthetic_sequence(n: usize) -> String {
    let mut x: u64 = 0x0FED_CBA9_8765_4321;
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

/// Gap intervals (0-based half-open): two heterozygous deletions with `GAP` bases between them,
/// every 10 kb from 50 kb.
fn gaps() -> Vec<(usize, usize)> {
    (0..90)
        .map(|i| {
            let first_end = 50_000 + i * 10_000 + DEL_LEN;
            (first_end, first_end + GAP)
        })
        .collect()
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
         ##ALT=<ID=DEL,Description=\"del\">\n\
         ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n\
         ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"x\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    );
    if with_events {
        for (first_end, second_start) in gaps() {
            let first_start = first_end - DEL_LEN;
            let second_end = second_start + DEL_LEN;
            for (s, e) in [(first_start, first_end), (second_start, second_end)] {
                text.push_str(&format!(
                    "{CONTIG}\t{s}\t.\tG\t<DEL>\t50\tPASS\tSVTYPE=DEL;END={e};SVLEN=-{DEL_LEN}\tGT\t0/1\n"
                ));
            }
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
             output_dir: {out}\noutput_basename: o\nrng_seed: sv gap depth\nnum_threads: 4\n",
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
        for d in &mut depth[start..(start + span).min(CONTIG_LEN)] {
            *d += 1;
        }
    }
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
fn a_gap_narrower_than_a_read_between_two_deletions_is_not_over_delivered() {
    let tmp = tempfile::tempdir().unwrap();
    let control = depth_track(&run(tmp.path(), "control", false));
    let events = depth_track(&run(tmp.path(), "events", true));
    let g = gaps();
    let ratio = mean_over(&events, &g) / mean_over(&control, &g);
    eprintln!(
        "gap depth {:.1} against control {:.1}, ratio {ratio:.3}",
        mean_over(&events, &g),
        mean_over(&control, &g)
    );

    // The deletions themselves are heterozygous, so the fixture must show them: about half the
    // control's depth. A run that ignored the VCF would pass the gap check below.
    let inside: Vec<(usize, usize)> = g.iter().map(|&(e, _)| (e - DEL_LEN + 20, e - 20)).collect();
    let del_ratio = mean_over(&events, &inside) / mean_over(&control, &inside);
    assert!(
        (0.40..=0.65).contains(&del_ratio),
        "deletion interior is {del_ratio:.2} of control; the events were not planted"
    );

    assert!(
        (0.97..=1.07).contains(&ratio),
        "gap depth is {ratio:.3} of the no-variant control; a budget of its own for the gap \
         reads about 1.12"
    );
}
