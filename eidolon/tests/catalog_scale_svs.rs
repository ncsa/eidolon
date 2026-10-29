//! Literal SVs at the sizes a real catalog carries (#806).
//!
//! **Why this file exists.** GIAB HG002 v5.0q, the first HG002 SV truth set on GRCh38,
//! carries 18,007 insertions (99th percentile 7.1 kb, largest 102 kb) and 10,165 deletions
//! (largest 114 kb) in its benchmark regions, all as literal REF/ALT. Replaying it through
//! `input_vcf` is #806 step 3 and #837 step 2. But the longest literal insertion ever checked
//! at read level was 600 bp (`long_insertion_rework.rs`, on H1N1, whose longest contig is
//! 2,280 bp), and literal deletions were checked at tens of bases. This pins catalog sizes.
//!
//! **Method.** One run on a 2 Mb synthetic reference (one contig, no `N`) carries three
//! insertions and three deletions, homozygous, at least 200 kb apart, with `mutation_rate: 0`
//! so probes match exactly. Insertions are checked with #516's method extended: 30-mer probes
//! at the head, the tail and nine evenly spaced interior points must all reach the reads at a
//! sensible fraction of the expected count. A realized head with a missing interior is exactly
//! the #516 failure. Each insertion is a disjoint slice of one synthetic sequence, so no probe
//! occurs in two insertions. Deletions are checked on the golden BAM (true origin, no aligner):
//! depth inside the deleted span must be near zero, and depth in the flanks must be unchanged,
//! which is the case where the check must NOT fire.
//!
//! **Known answer.** Homozygous at 30x with 151 bp reads: a 30-mer of inserted sequence is
//! wholly contained in about 30 * (151 - 29) / 151 = 24 reads; a deleted base is covered by
//! no read. Neither number depends on the code under test.
//!
//! **Why not `ecoli.fa`.** It is repeat-masked, with an `N` run about every 13 kb, so no
//! 100 kb deletion can avoid one, and a deletion spanning an `N` run is currently truncated
//! at it (#863). That defect is real and tracked there; this file measures size alone.

mod common;

use common::{eidolon, revcomp, synthetic_insert};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::Path;

const CONTIG: &str = "synth";
const CONTIG_LEN: usize = 2_000_000;
const COVERAGE: f64 = 30.0;
const READ_LEN: usize = 151;

/// A deterministic random ACGT sequence (xorshift): no `N`, no structure.
fn random_acgt(n: usize, mut state: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect()
}

fn write_reference(path: &Path, seq: &[u8]) {
    let mut text = format!(">{CONTIG}\n");
    for chunk in seq.chunks(80) {
        text.push_str(std::str::from_utf8(chunk).unwrap());
        text.push('\n');
    }
    std::fs::write(path, text).unwrap();
}

struct Planned {
    id: String,
    pos: usize, // 1-based VCF POS (the anchor base)
    kind: &'static str,
    len: usize,
    insert: Option<String>,
}

fn record(seq: &[u8], p: &Planned) -> String {
    let anchor = seq[p.pos - 1] as char;
    match p.kind {
        "INS" => format!(
            "{CONTIG}\t{}\t{}\t{anchor}\t{anchor}{}\t60\tPASS\tSVTYPE=INS\tGT\t1/1",
            p.pos,
            p.id,
            p.insert.as_ref().unwrap()
        ),
        _ => {
            let refe: String = seq[p.pos - 1..p.pos + p.len]
                .iter()
                .map(|&b| b as char)
                .collect();
            format!(
                "{CONTIG}\t{}\t{}\t{refe}\t{anchor}\t60\tPASS\tSVTYPE=DEL\tGT\t1/1",
                p.pos, p.id
            )
        }
    }
}

/// Reads whose SEQUENCE contains `probe` in either orientation (FASTQ line 2 only:
/// Phred+33 overlaps the DNA alphabet).
fn seq_hits(fastqs: &[String], probe: &str) -> usize {
    let rc = revcomp(probe);
    fastqs
        .iter()
        .map(|s| {
            s.lines()
                .skip(1)
                .step_by(4)
                .filter(|l| l.contains(probe) || l.contains(rc.as_str()))
                .count()
        })
        .sum()
}

fn read_gz(path: &Path) -> String {
    use std::io::Read;
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

/// Per-base true-origin depth from the golden BAM (reference-consuming CIGAR ops).
fn golden_depth(bam_path: &Path, len: usize) -> Vec<u16> {
    let mut depth = vec![0u16; len + 1];
    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path).unwrap());
    reader.read_header().unwrap();
    for r in reader.records() {
        let r = r.unwrap();
        let Some(Ok(start)) = r.alignment_start() else {
            continue;
        };
        let mut pos = usize::from(start); // 1-based
        for op in r.cigar().iter() {
            let op = op.unwrap();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for d in depth.iter_mut().skip(pos).take(op.len()) {
                        *d = d.saturating_add(1);
                    }
                    pos += op.len();
                }
                Kind::Deletion | Kind::Skip => pos += op.len(),
                _ => {}
            }
        }
    }
    depth
}

/// Mean depth over [from, to), 1-based positions.
fn mean(depth: &[u16], from: usize, to: usize) -> f64 {
    let s: u64 = depth[from..to].iter().map(|&d| d as u64).sum();
    s as f64 / (to - from) as f64
}

#[test]
fn catalog_scale_literal_svs_are_realized_in_the_reads() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("run");
    std::fs::create_dir_all(&out).unwrap();
    let seq = random_acgt(CONTIG_LEN, 0x5eed_0806);
    let reference = out.join("synth.fa");
    write_reference(&reference, &seq);
    // One synthetic sequence, sliced into disjoint inserts, so no probe occurs twice.
    let pool = synthetic_insert(62_000);

    // Three insertions and three deletions at catalog sizes, >= 200 kb apart.
    let mut plan: Vec<Planned> = Vec::new();
    let mut offset = 0usize;
    for (i, (len, pos)) in [
        (2_000usize, 100_000usize),
        (10_000, 350_000),
        (50_000, 600_000),
    ]
    .into_iter()
    .enumerate()
    {
        let insert = pool[offset..offset + len].to_string();
        offset += len;
        plan.push(Planned {
            id: format!("ins{i}"),
            pos,
            kind: "INS",
            len,
            insert: Some(insert),
        });
    }
    for (i, (len, pos)) in [
        (2_000usize, 900_000usize),
        (10_000, 1_150_000),
        (100_000, 1_500_000),
    ]
    .into_iter()
    .enumerate()
    {
        plan.push(Planned {
            id: format!("del{i}"),
            pos,
            kind: "DEL",
            len,
            insert: None,
        });
    }

    let vcf = out.join("in.vcf");
    let mut text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={CONTIG_LEN}>\n\
         ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    );
    for p in &plan {
        text.push_str(&record(&seq, p));
        text.push('\n');
    }
    std::fs::write(&vcf, text).unwrap();

    let cfg = out.join("c.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {}\nread_len: {READ_LEN}\ncoverage: {COVERAGE}\nploidy: 2\n\
             paired_ended: true\nfragment_mean: 400\nfragment_st_dev: 50\n\
             input_vcf: {}\nmutation_rate: 0\nsv_rate_scale: 0.0\n\
             produce_vcf: true\nproduce_fastq: true\nproduce_bam: true\n\
             overwrite_output: true\noutput_dir: {}\noutput_filename: o\n\
             rng_seed: catalog scale svs\nnum_threads: 4\n",
            reference.display(),
            vcf.display(),
            out.display()
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();

    // Every planted record is in the truth VCF.
    let truth = read_gz(&out.join("o.vcf.gz"));
    for p in &plan {
        assert!(
            truth
                .lines()
                .any(|l| l.split('\t').nth(2) == Some(p.id.as_str())),
            "{} ({} {} bp at {}) is missing from the truth VCF",
            p.id,
            p.kind,
            p.len,
            p.pos
        );
    }

    let mut failures = Vec::new();

    // Insertions: head, tail and nine interior probes, each at a sensible fraction of the
    // ~24 reads a homozygous 30x run puts over a 30-mer.
    let fastqs = vec![
        read_gz(&out.join("o_r1.fastq.gz")),
        read_gz(&out.join("o_r2.fastq.gz")),
    ];
    let expected = COVERAGE * (READ_LEN - 29) as f64 / READ_LEN as f64;
    for p in plan.iter().filter(|p| p.kind == "INS") {
        let ins = p.insert.as_ref().unwrap();
        let mut offsets: Vec<usize> = (1..10).map(|k| k * (p.len - 30) / 10).collect();
        offsets.insert(0, 0);
        offsets.push(p.len - 30);
        let hits: Vec<usize> = offsets
            .iter()
            .map(|&o| seq_hits(&fastqs, &ins[o..o + 30]))
            .collect();
        eprintln!(
            "[{} INS {} bp] probe hits head..tail: {hits:?} (expected ~{expected:.0})",
            p.id, p.len
        );
        let low = hits
            .iter()
            .filter(|&&h| (h as f64) < 0.3 * expected)
            .count();
        let high = hits
            .iter()
            .filter(|&&h| (h as f64) > 2.5 * expected)
            .count();
        if low > 0 || high > 0 {
            failures.push(format!(
                "{} INS {} bp: {low} probe(s) under 30% and {high} over 250% of the expected {expected:.0}: {hits:?}",
                p.id, p.len
            ));
        }
    }

    // Deletions: near-zero true-origin depth inside the span, unchanged flanks.
    let depth = golden_depth(&out.join("o.bam"), seq.len());
    let genome_mean = mean(&depth, 1, seq.len());
    for p in plan.iter().filter(|p| p.kind == "DEL") {
        // Deleted bases are pos+1 ..= pos+len. Skip 200 bp at each edge for read overhang.
        let inside = mean(&depth, p.pos + 1 + 200, p.pos + p.len - 200);
        let left = mean(&depth, p.pos - 6_000, p.pos - 1_000);
        let right = mean(&depth, p.pos + p.len + 1_000, p.pos + p.len + 6_000);
        eprintln!(
            "[{} DEL {} bp] depth inside {inside:.2}, flanks {left:.1} / {right:.1}, genome {genome_mean:.1}",
            p.id, p.len
        );
        if inside > 0.5 {
            failures.push(format!(
                "{} DEL {} bp: depth {inside:.2} inside the deleted span (expected ~0)",
                p.id, p.len
            ));
        }
        for (side, d) in [("left", left), ("right", right)] {
            if !(0.85 * genome_mean..=1.15 * genome_mean).contains(&d) {
                failures.push(format!(
                    "{} DEL {} bp: {side} flank depth {d:.1} vs genome {genome_mean:.1} (must be unchanged)",
                    p.id, p.len
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "catalog-scale SVs not realized:\n  {}",
        failures.join("\n  ")
    );
}
