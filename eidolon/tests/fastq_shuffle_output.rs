//! gen-reads writes its FASTQs in random order by default (#870).
//!
//! Genome order made bwa-mem2 need 166 GB instead of 90 GB on a 21x T2T-CHM13 simulation, with
//! identical alignments either way. A real sequencer's FASTQ is unordered, so that is the default.
//!
//! **Known answer.** Read names carry each read's origin (`…_<start>_<end>_<index>/1`). In genome
//! order nearly every read's neighbor starts within a few hundred bases of it; in random order on
//! H1N1 (13.5 kb), a neighbor on the same contig within 200 bp is rare. **Must not fire:**
//! `shuffle_fastq: false` keeps genome order, and the shuffle changes nothing but order. The same
//! records come out, mates stay in register, and the truth VCF and golden BAM are byte-identical,
//! because the shuffle draws from its own RNG child.

mod common;

use common::{eidolon, h1n1_reference};
use noodles::bam;
use std::io::Read;
use std::path::Path;

fn gz_text(path: &Path) -> String {
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

fn run(dir: &Path, name: &str, shuffle: Option<bool>, threads: usize) {
    let cfg = dir.join(format!("{name}.yml"));
    let shuffle_line = shuffle
        .map(|s| format!("shuffle_fastq: {s}\n"))
        .unwrap_or_default();
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 101\ncoverage: 30\npaired_ended: true\n\
             fragment_mean: 300\nfragment_st_dev: 30\n\
             produce_fastq: true\nproduce_vcf: true\nproduce_bam: true\n\
             overwrite_output: true\noutput_dir: {out}\noutput_filename: {name}\n\
             rng_seed: fastq shuffle\nnum_threads: {threads}\n{shuffle_line}",
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
}

/// (contig, start) parsed from a read name, counting fields from the right because contig names
/// may contain underscores.
fn origin(name: &str) -> (String, i64) {
    let core = name
        .trim_start_matches('@')
        .trim_start_matches("EIDOLON_generated_");
    let core = core.rsplit_once('/').map_or(core, |(c, _)| c);
    let parts: Vec<&str> = core.rsplitn(4, '_').collect();
    (parts[3].to_string(), parts[2].parse().unwrap())
}

fn names(fastq: &str) -> Vec<String> {
    fastq.lines().step_by(4).map(str::to_string).collect()
}

/// Fraction of adjacent reads on the same contig with starts within 200 bp.
fn near_neighbors(names: &[String]) -> f64 {
    let o: Vec<(String, i64)> = names.iter().map(|n| origin(n)).collect();
    let near = o
        .windows(2)
        .filter(|w| w[0].0 == w[1].0 && (w[0].1 - w[1].1).abs() <= 200)
        .count();
    near as f64 / (o.len() - 1) as f64
}

#[test]
fn default_output_is_shuffled_and_off_keeps_genome_order() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    run(d, "on", None, 2);
    run(d, "off", Some(false), 2);

    let on1 = names(&gz_text(&d.join("on_r1.fastq.gz")));
    let on2 = names(&gz_text(&d.join("on_r2.fastq.gz")));
    let off1 = names(&gz_text(&d.join("off_r1.fastq.gz")));
    // Random order on H1N1 puts a neighbor on the same contig ~15% of the time and within 200 bp
    // of it ~23% of that (contigs ~1.7 kb): ~0.035 expected, sigma ~0.006 at 1000 reads. So 1000
    // reads keep the < 0.1 assertion below far from noise.
    assert!(
        on1.len() > 1000,
        "only {} reads; the order statistics would be weak",
        on1.len()
    );

    let (near_on, near_off) = (near_neighbors(&on1), near_neighbors(&off1));
    eprintln!("adjacent reads within 200 bp: shuffled {near_on:.4}, genome order {near_off:.4}");
    assert!(
        near_off > 0.9,
        "shuffle_fastq: false should keep genome order: {near_off:.4}"
    );
    assert!(
        near_on < 0.1,
        "the default output is not shuffled: {near_on:.4}"
    );

    // Same records, mates in register.
    for (a, b) in on1.iter().zip(&on2) {
        assert_eq!(
            a.trim_end_matches("/1"),
            b.trim_end_matches("/2"),
            "mates out of register"
        );
    }
    let (mut s_on, mut s_off) = (on1.clone(), off1.clone());
    s_on.sort();
    s_off.sort();
    assert_eq!(
        s_on, s_off,
        "the shuffle must change order only, not which reads exist"
    );

    // The shuffle draws from its own RNG child: the truth VCF and the golden BAM are unaffected.
    assert_eq!(
        gz_text(&d.join("on.vcf.gz")),
        gz_text(&d.join("off.vcf.gz")),
        "truth VCF changed"
    );
    let bam_records = |p: &Path| -> Vec<String> {
        let mut r = bam::io::Reader::new(std::fs::File::open(p).unwrap());
        r.read_header().unwrap();
        r.records().map(|x| format!("{:?}", x.unwrap())).collect()
    };
    assert_eq!(
        bam_records(&d.join("on.bam")),
        bam_records(&d.join("off.bam")),
        "golden BAM changed"
    );
}

#[test]
fn shuffled_output_is_identical_across_thread_counts() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    run(d, "t1", None, 1);
    run(d, "t4", None, 4);
    for mate in ["r1", "r2"] {
        assert_eq!(
            gz_text(&d.join(format!("t1_{mate}.fastq.gz"))),
            gz_text(&d.join(format!("t4_{mate}.fastq.gz"))),
            "{mate} differs between 1 and 4 threads"
        );
    }
}
