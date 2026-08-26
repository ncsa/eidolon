//! Acceptance criteria for the #516 long-insertion rework.
//!
//! **Why a separate file from `sv_support_matrix.rs`.** That file is the support
//! *matrix* — one cell per (SV type, path) combination, answering "can eidolon
//! reproduce this at all". This file is narrower and harder: it pins the properties
//! a *correct* long-insertion implementation must have, most of which the previous
//! attempt (PRs #573/#574, reverted in #575) got wrong while still passing the
//! matrix cell.
//!
//! That distinction is the whole reason this file exists. The reverted
//! implementation was validated on Delta and genuinely worked for its primary
//! purpose — middle and tail sequence reached the reads (48/48, 49/48, 44/44 probe
//! hits on S. cerevisiae, 2026-08-22). A review then found six further defects it
//! had shipped alongside that success. Every one of them is a criterion below, so
//! "the insertion is realized" can never again be mistaken for "the insertion is
//! realized correctly".
//!
//! **Status and non-vacuity, both measured rather than asserted.** Every criterion
//! here was run against BOTH the current (unwired) `develop` and the pre-revert
//! implementation at `6c2af32`, where the wiring was present. They discriminate
//! perfectly, which is what makes the set worth having:
//!
//! | criterion | pre-revert `6c2af32` | `develop` today |
//! |---|---|---|
//! | reaches interior and tail | **passes** | **FAILS** — head 36/43/49, middle 0, tail 0 |
//! | respects zygosity | FAILS — het/hom ratio 1.02 | guard holds |
//! | reports allelic depth | FAILS — `0/1:0,0:0:.` | guard holds |
//! | neighbours keep support | FAILS — DP 63->31, 66->35 | guard holds |
//! | no zero-reference-span records | FAILS — 275 of 8232 mapped | guard holds |
//! | read names unique across chunks | FAILS — 158 names appear 4x | guard holds |
//! | junction evidence vs fragment length | FAILS — 0 hits at fragment_mean 600 | guard holds |
//!
//! Read that table as the definition of done: the rework has to reach the top-left
//! cell (realize the insertion) without giving up any of the right-hand column.
//! Exactly one criterion is red on `develop` today, and it is #516 itself; the other
//! six are regression guards, each already proven to catch the defect it names.
//!
//! Note that `long_insertion_respects_zygosity` is red on `develop` for a *derived*
//! reason — its middle probe reads 0 for both genotypes because no interior sequence
//! reaches the reads at all. It only starts measuring zygosity once criterion 1
//! passes, and its ratio assertion is what catches the missing coin at that point.

mod common;

use common::{eidolon, h1n1_reference, revcomp, synthetic_insert};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::{Path, PathBuf};

const CONTIG: &str = "H1N1_PB2";
const CONTIG_LEN: usize = 2280;
const ANCHOR: usize = 400; // 1-based VCF POS

fn ref_base_at(pos_1based: usize) -> char {
    let text = std::fs::read_to_string(h1n1_reference()).unwrap();
    let mut seq = String::new();
    let mut in_contig = false;
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            if in_contig {
                break;
            }
            in_contig = name.split_whitespace().next() == Some(CONTIG);
            continue;
        }
        if in_contig {
            seq.push_str(line.trim());
        }
    }
    seq.as_bytes()[pos_1based - 1] as char
}

struct Cell {
    dir: PathBuf,
}

/// Run gen-reads over `records` (tab-delimited VCF body lines). `extra_cfg` appends
/// raw YAML so individual criteria can vary fragment length, chunking, etc.
fn run(dir: &Path, name: &str, records: &[String], extra_cfg: &str) -> Cell {
    let out = dir.join(name);
    std::fs::create_dir_all(&out).unwrap();
    let vcf = out.join("in.vcf");
    let mut text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={CONTIG_LEN}>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    );
    for r in records {
        text.push_str(r);
        text.push('\n');
    }
    std::fs::write(&vcf, text).unwrap();

    let cfg = out.join("c.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 100\ncoverage: 60\nploidy: 2\npaired_ended: true\n\
             fragment_mean: 250\nfragment_st_dev: 30\ninput_vcf: {vcf}\n\
             produce_vcf: true\nproduce_fastq: true\nproduce_bam: true\n\
             sv_rate_scale: 0.0\noverwrite_output: true\n\
             output_dir: {out}\noutput_filename: o\nrng_seed: i516 {name}\nnum_threads: 1\n{extra_cfg}",
            ref = h1n1_reference().display(),
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
    Cell { dir: out }
}

impl Cell {
    /// Reads whose SEQUENCE contains `probe`, either orientation. FASTQ line 2 only —
    /// Phred+33 overlaps the DNA alphabet, so grepping whole records produces false hits.
    fn seq_hits(&self, probe: &str) -> usize {
        let rc = revcomp(probe);
        let mut n = 0;
        for f in ["o_r1.fastq.gz", "o_r2.fastq.gz"] {
            let path = self.dir.join(f);
            if !path.exists() {
                continue;
            }
            let file = std::fs::File::open(&path).unwrap();
            let mut s = String::new();
            {
                use std::io::Read;
                flate2::read::MultiGzDecoder::new(file)
                    .read_to_string(&mut s)
                    .unwrap();
            }
            n += s
                .lines()
                .skip(1)
                .step_by(4)
                .filter(|l| l.contains(probe) || l.contains(rc.as_str()))
                .count();
        }
        n
    }

    fn bam_records(&self) -> Vec<bam::Record> {
        let file = std::fs::File::open(self.dir.join("o.bam")).unwrap();
        let mut reader = bam::io::Reader::new(file);
        reader.read_header().unwrap();
        reader.records().map(|r| r.unwrap()).collect()
    }

    /// The FORMAT sample field of the truth record with this ID, if present.
    fn truth_format(&self, id: &str) -> Option<String> {
        let file = std::fs::File::open(self.dir.join("o.vcf.gz")).unwrap();
        let mut s = String::new();
        {
            use std::io::Read;
            flate2::read::MultiGzDecoder::new(file)
                .read_to_string(&mut s)
                .unwrap();
        }
        for line in s.lines() {
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() > 9 && f[2] == id {
                return Some(f[9].to_string());
            }
        }
        None
    }
}

fn ins_record(pos: usize, insert: &str, gt: &str, id: &str) -> String {
    let base = ref_base_at(pos);
    format!("{CONTIG}\t{pos}\t{id}\t{base}\t{base}{insert}\t60\tPASS\t.\tGT\t{gt}")
}

// ── CRITERION 1: the original defect ────────────────────────────────────────
// The insertion's interior and far end must reach the reads, not just its head.
// This is #516 itself. `sv_support_matrix.rs` currently pins the INVERSE as
// characterization (sizes >=200 must show middle==0 and tail==0); when this passes,
// that characterization must be flipped and this becomes the guard.

#[test]
fn long_insertions_reach_interior_and_tail() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for size in [200usize, 300, 600] {
        let insert = synthetic_insert(size);
        let cell = run(
            tmp.path(),
            &format!("reach_{size}"),
            &[ins_record(ANCHOR, &insert, "1/1", "ins")],
            "",
        );
        let head = cell.seq_hits(&insert[..30]);
        let mid = cell.seq_hits(&insert[size / 2 - 15..size / 2 + 15]);
        let tail = cell.seq_hits(&insert[size - 30..]);
        eprintln!("[reach {size}bp] head={head} middle={mid} tail={tail}");
        if head == 0 || mid == 0 || tail == 0 {
            failures.push(format!(
                "{size}bp: head={head} middle={mid} tail={tail} (all three must be > 0)"
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "inserted sequence does not reach the reads across its full length:\n  {}",
        failures.join("\n  ")
    );
}

// ── CRITERION 1b: SEVERAL long insertions on one contig ─────────────────────
// The natural way to test long insertions is to plant a size range at once, and
// `input_vcf` is the documented path for planting specific events -- so two or
// more long insertions sharing a contig is the common case, not an edge case.
// It is also exactly the scenario validated on Delta on 2026-08-22 (200/600/1200bp
// on one S. cerevisiae contig).
//
// This is RED for the first implementation of the rework, which sampled one
// altered haplotype per sub-region and fell back to the pre-#516 head-only
// behaviour whenever a sub-region held more than one. Sub-regions are large -- a
// whole contig, split only by coverage multipliers -- so that fallback caught
// every multi-insertion case: measured head 10/14/15, middle 0/0/0, tail 0/0/0.

#[test]
fn several_long_insertions_on_one_contig_are_all_realized() {
    let tmp = tempfile::tempdir().unwrap();
    // Distinct novel sequence per event, so a probe cannot be satisfied by the
    // wrong insertion.
    let sizes = [200usize, 600, 1200];
    let inserts: Vec<String> = sizes
        .iter()
        .enumerate()
        .map(|(i, &n)| {
            let mut s = String::with_capacity(n);
            let mut x: u64 = 0x9E37_79B9_7F4A_7C15u64.wrapping_add((i as u64).wrapping_mul(7919));
            for _ in 0..n {
                x = x
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                s.push(match (x >> 33) & 3 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                });
            }
            s
        })
        .collect();
    // Spaced far apart, but all inside the one non-N region of this contig.
    let anchors = [400usize, 900, 1500];
    let records: Vec<String> = anchors
        .iter()
        .zip(inserts.iter())
        .enumerate()
        .map(|(i, (&pos, ins))| ins_record(pos, ins, "1/1", &format!("ins{i}")))
        .collect();

    let cell = run(tmp.path(), "multi", &records, "");

    let mut failures = Vec::new();
    for (i, (&size, ins)) in sizes.iter().zip(inserts.iter()).enumerate() {
        let head = cell.seq_hits(&ins[..30]);
        let mid = cell.seq_hits(&ins[size / 2 - 15..size / 2 + 15]);
        let tail = cell.seq_hits(&ins[size - 30..]);
        eprintln!("[multi #{i} {size}bp] head={head} middle={mid} tail={tail}");
        if head == 0 || mid == 0 || tail == 0 {
            failures.push(format!(
                "#{i} ({size}bp): head={head} middle={mid} tail={tail}"
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "insertions sharing a contig are not all fully realized -- a coordinate map \
         describing only ONE insertion makes the others fall back to head-only:\n  {}",
        failures.join("\n  ")
    );
}

// ── CRITERION 2: zygosity ───────────────────────────────────────────────────
// The reverted implementation called generate_read with an EMPTY variant map, so the
// heterozygous coin never ran and 0/1 rendered identically to 1/1. Measured then:
// het head/mid/tail 47/50/42 against hom 34/47/37 — indistinguishable.

#[test]
fn long_insertion_respects_zygosity() {
    let tmp = tempfile::tempdir().unwrap();
    let size = 600usize;
    let insert = synthetic_insert(size);
    let probe = &insert[size / 2 - 15..size / 2 + 15];

    let het = run(
        tmp.path(),
        "zyg_het",
        &[ins_record(ANCHOR, &insert, "0/1", "ins")],
        "",
    )
    .seq_hits(probe);
    let hom = run(
        tmp.path(),
        "zyg_hom",
        &[ins_record(ANCHOR, &insert, "1/1", "ins")],
        "",
    )
    .seq_hits(probe);

    eprintln!(
        "[zygosity] het={het} hom={hom} ratio={:.2}",
        het as f64 / hom.max(1) as f64
    );
    assert!(
        hom > 0 && het > 0,
        "both genotypes must produce reads (het={het} hom={hom})"
    );
    let ratio = het as f64 / hom as f64;
    assert!(
        (0.3..=0.7).contains(&ratio),
        "heterozygous insertion delivered {het} probe hits against {hom} homozygous \
         (ratio {ratio:.2}); a het event sits on one haplotype so it must be ~0.5. A ratio \
         near 1.0 means the heterozygous coin is not being applied at all."
    );
}

// ── CRITERION 3: allelic depth ──────────────────────────────────────────────
// The reverted implementation used a throwaway AdCounter per fragment, so the golden
// VCF reported `0,0:0:.` — DP=0 — for every long insertion, while neighbouring SNPs
// reported normal depth.

#[test]
fn long_insertion_reports_allelic_depth() {
    let tmp = tempfile::tempdir().unwrap();
    let insert = synthetic_insert(600);
    let cell = run(
        tmp.path(),
        "ad",
        &[ins_record(ANCHOR, &insert, "0/1", "bigins")],
        "",
    );
    let fmt = cell
        .truth_format("bigins")
        .expect("the insertion must survive into the truth VCF");
    eprintln!("[allelic depth] FORMAT sample field = {fmt}");
    let dp: usize = fmt
        .split(':')
        .nth(2)
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    assert!(
        dp > 0,
        "golden VCF reports DP={dp} for a long insertion ({fmt}) — reads covering the \
         event are not being counted into the allelic-depth accumulator, so AD/DP/AF are \
         a statement about nothing."
    );
}

// ── CRITERION 4: neighbouring variants ──────────────────────────────────────
// Reads generated for the insertion carried reference bases at every other locus,
// so a nearby het SNP lost its alt support. Measured on the reverted build: a het SNP
// 40bp away went from AD 34,37 / DP 71 to AD 9,18 / DP 27.

#[test]
fn variants_near_a_long_insertion_keep_their_support() {
    let tmp = tempfile::tempdir().unwrap();
    let insert = synthetic_insert(600);
    let (l, r) = (ANCHOR - 60, ANCHOR + 60);
    let flip = |c: char| match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        _ => 'C',
    };
    let snp = |pos: usize, id: &str| {
        let b = ref_base_at(pos);
        format!(
            "{CONTIG}\t{pos}\t{id}\t{b}\t{}\t60\tPASS\t.\tGT\t0/1",
            flip(b)
        )
    };

    let alone = run(
        tmp.path(),
        "neigh_snps_only",
        &[snp(l, "snpL"), snp(r, "snpR")],
        "",
    );
    let with_ins = run(
        tmp.path(),
        "neigh_with_ins",
        &[
            snp(l, "snpL"),
            ins_record(ANCHOR, &insert, "0/1", "bigins"),
            snp(r, "snpR"),
        ],
        "",
    );

    let dp = |c: &Cell, id: &str| -> usize {
        c.truth_format(id)
            .and_then(|f| f.split(':').nth(2).map(str::to_string))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    };

    let mut failures = Vec::new();
    for id in ["snpL", "snpR"] {
        let (a, b) = (dp(&alone, id), dp(&with_ins, id));
        eprintln!("[neighbour {id}] DP alone={a}  with insertion={b}");
        if a == 0 {
            failures.push(format!("{id}: control DP is 0 — fixture broken"));
        } else if (b as f64) < 0.7 * a as f64 {
            failures.push(format!(
                "{id}: DP fell {a} -> {b} when a long insertion was added 60bp away; reads \
                 generated for the insertion are rendering reference at neighbouring loci"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n  "));
}

// ── CRITERION 5: BAM records must be well formed ────────────────────────────
// A read lying wholly inside novel sequence has no honest reference alignment. The
// reverted implementation emitted 291 records with CIGAR `100I` at a single position —
// zero reference-consuming operations, which is not a valid alignment. Per the design
// decision for this rework, such reads are emitted UNMAPPED (their mate stays mapped)
// and carry provenance in a tag, so the golden BAM keeps working as an answer key for
// #449 without asserting an alignment that does not exist.

#[test]
fn no_record_claims_a_zero_reference_span_alignment() {
    let tmp = tempfile::tempdir().unwrap();
    let insert = synthetic_insert(600);
    let cell = run(
        tmp.path(),
        "bamvalid",
        &[ins_record(ANCHOR, &insert, "1/1", "ins")],
        "",
    );

    let mut offenders = 0usize;
    let mut total_mapped = 0usize;
    for rec in cell.bam_records() {
        if rec.flags().is_unmapped() {
            continue;
        }
        total_mapped += 1;
        let consumes_reference = rec.cigar().iter().any(|op| {
            let op = op.unwrap();
            matches!(
                op.kind(),
                Kind::Match
                    | Kind::Deletion
                    | Kind::Skip
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
            ) && op.len() > 0
        });
        if !consumes_reference {
            offenders += 1;
        }
    }
    eprintln!("[bam validity] {offenders} zero-reference-span records of {total_mapped} mapped");
    assert_eq!(
        offenders, 0,
        "{offenders} of {total_mapped} MAPPED records consume no reference (e.g. an all-I \
         CIGAR). A read lying entirely inside novel inserted sequence has no reference \
         alignment and must be emitted unmapped rather than placed at the anchor."
    );
}

// ── CRITERION 6: read names ─────────────────────────────────────────────────
// The reverted writer named reads `{prefix}_{counter}` with the counter restarting per
// call, rather than the position-keyed name added for #210. With two long insertions in
// different chunks of one contig, 158 QNAMEs appeared 4 times each.

#[test]
fn long_insertion_read_names_stay_unique_across_chunks() {
    let tmp = tempfile::tempdir().unwrap();
    let a = synthetic_insert(300);
    let b = synthetic_insert(350);
    let cell = run(
        tmp.path(),
        "qname",
        &[
            ins_record(400, &a, "1/1", "insA"),
            ins_record(1200, &b, "1/1", "insB"),
        ],
        "chunk_size: 500\n",
    );

    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for rec in cell.bam_records() {
        let name = String::from_utf8_lossy(rec.name().unwrap().as_ref()).to_string();
        *counts.entry(name).or_default() += 1;
    }
    let over: Vec<_> = counts.iter().filter(|&(_, &n)| n > 2).collect();
    eprintln!("[qname] {} names appear more than twice", over.len());
    assert!(
        over.is_empty(),
        "{} read name(s) appear more than twice in the BAM (a paired read may appear at \
         most twice). Colliding QNAMEs break mate pairing and MarkDuplicates — see #210. \
         Example: {:?}",
        over.len(),
        over.iter().take(3).collect::<Vec<_>>()
    );
}

// ── CRITERION 7: robustness to fragment length ──────────────────────────────
// The reverted implementation hardcoded its sampling window to 2*read_len regardless of
// fragment length. At fragment_mean=600 it produced 482 reads carrying NO junction
// evidence at all (max insertion run = 1bp, i.e. sequencing error), and a real aligner
// found zero soft clips >=30bp — strictly worse than before the change, where the head
// at least reached the reads.

#[test]
fn junction_evidence_survives_longer_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let size = 600usize;
    let insert = synthetic_insert(size);
    let head = &insert[..30];

    let mut failures = Vec::new();
    for frag in [250usize, 400, 600] {
        let cell = run(
            tmp.path(),
            &format!("frag_{frag}"),
            &[ins_record(ANCHOR, &insert, "1/1", "ins")],
            &format!("fragment_mean: {frag}\n"),
        );
        // A junction read carries reference AND >=30 contiguous novel bases; the head
        // probe is exactly that test, and it is what a caller assembles from.
        let hits = cell.seq_hits(head);
        eprintln!("[fragment_mean {frag}] head-probe hits = {hits}");
        if hits == 0 {
            failures.push(format!(
                "fragment_mean={frag}: no read carries the insertion's first 30 novel bases, \
                 so there is no assemblable junction evidence at all"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n  "));
}

// ── CRITERION 9: the golden BAM must DESCRIBE the insertion, not just carry it ──────
//
// #516 got the inserted bases into the reads. #589 was the other half: nothing told the
// CIGAR builder which fragment bases were inserted, so a read crossing the anchor was
// written as pure `M` and claimed inserted bases as reference matches. Measured before
// the fix on a 500 bp insertion at read_len 150: the largest `I` operation anywhere in
// the BAM was **2 bp** — sequencing-error indels only — and one read at POS 693 carried
// 108 reference bases followed by 42 insertion bases (zero mismatches against the novel
// sequence) under the CIGAR `150M`.
//
// The general property, and the one asserted here, is stronger than "an I op exists":
// **every mapped read must reconstruct against the reference from its own POS and
// CIGAR.** That is what a golden BAM being alignment truth means, and it is the check
// that caught this (it also cleared eidolon of the sibling defect in ncsa/neat#326,
// where every gapped reverse read fails to reconstruct).

/// The full base sequence of CONTIG, uppercased — the reference a golden-BAM alignment
/// is supposed to describe.
fn contig_bases() -> Vec<u8> {
    let text = std::fs::read_to_string(h1n1_reference()).unwrap();
    let mut seq = String::new();
    let mut in_contig = false;
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            if in_contig {
                break;
            }
            in_contig = name.split_whitespace().next() == Some(CONTIG);
            continue;
        }
        if in_contig {
            seq.push_str(line.trim());
        }
    }
    seq.to_ascii_uppercase().into_bytes()
}

/// Walk a read's CIGAR against the reference and count mismatches over the M/=/X blocks.
/// Returns `(mismatches, compared)`.
fn reconstruct(
    reference: &[u8],
    pos_1based: usize,
    cigar: &[(char, usize)],
    seq: &[u8],
) -> (usize, usize) {
    let (mut ri, mut qi) = (pos_1based - 1, 0usize);
    let (mut mism, mut total) = (0usize, 0usize);
    for &(op, n) in cigar {
        match op {
            'M' | '=' | 'X' => {
                for k in 0..n {
                    if ri + k < reference.len() && qi + k < seq.len() {
                        total += 1;
                        if reference[ri + k].to_ascii_uppercase()
                            != seq[qi + k].to_ascii_uppercase()
                        {
                            mism += 1;
                        }
                    }
                }
                ri += n;
                qi += n;
            }
            'I' | 'S' => qi += n,
            'D' | 'N' => ri += n,
            _ => {}
        }
    }
    (mism, total)
}

#[test]
fn the_golden_bam_cigar_describes_a_long_insertion() {
    let dir = tempfile::tempdir().unwrap();
    let insert = synthetic_insert(500);
    let cell = run(
        dir.path(),
        "bam_cigar_ins",
        &[ins_record(ANCHOR, &insert, "1/1", "ins500")],
        "read_len: 150\nfragment_mean: 300\n",
    );

    let contig_seq = contig_bases();
    let mut largest_i = 0usize;
    let mut broken = 0usize;
    let mut checked = 0usize;

    // H1N1 has EIGHT contigs and the BAM carries reads from all of them. Reconstructing a
    // read against the wrong contig's bases reports ~4340 of 5212 "broken" and means
    // nothing — so resolve each record's reference name and keep only CONTIG's.
    let ref_names: Vec<String> = {
        let file = std::fs::File::open(cell.dir.join("o.bam")).unwrap();
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header().unwrap();
        header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k.as_ref()).to_string())
            .collect()
    };

    for record in cell.bam_records() {
        let on_target = record
            .reference_sequence_id()
            .and_then(|r| r.ok())
            .and_then(|id| ref_names.get(id))
            .is_some_and(|n| n == CONTIG);
        if !on_target {
            continue;
        }
        let flags = record.flags();
        if flags.is_unmapped() {
            continue; // reads wholly inside the insertion are unmapped by design
        }
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        let ops: Vec<(char, usize)> = record
            .cigar()
            .iter()
            .map(|o| {
                let o = o.unwrap();
                let c = match o.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::SoftClip => 'S',
                    Kind::Skip => 'N',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                    _ => '?',
                };
                (c, o.len())
            })
            .collect();
        for &(c, n) in &ops {
            if c == 'I' && n > largest_i {
                largest_i = n;
            }
        }
        let seq: Vec<u8> = record.sequence().iter().collect();
        let (mism, total) = reconstruct(&contig_seq, usize::from(start), &ops, &seq);
        if total > 0 {
            checked += 1;
            if (mism as f64) / (total as f64) > 0.20 {
                broken += 1;
            }
        }
    }

    assert!(
        checked > 100,
        "only {checked} mapped reads examined — too few to conclude"
    );

    // The insertion has to appear AS an insertion. A 150 bp read can carry at most 149
    // inserted bases and still have one reference-anchored base, so that is the ceiling;
    // anything in the tens proves the op is real rather than a stray error indel (which
    // top out at a couple of bases).
    assert!(
        largest_i >= 50,
        "the largest insertion CIGAR op in the golden BAM is {largest_i} bp for a planted \
         500 bp insertion — the bases are in the reads but no CIGAR says so, which is #589. \
         Sequencing-error indels alone produce 1-2 bp ops."
    );

    // The property that generalises: the recorded alignment must describe the read.
    assert_eq!(
        broken, 0,
        "{broken} of {checked} mapped reads do not reconstruct against the reference from \
         their own POS+CIGAR. A read crossing the insertion anchor recorded as pure M \
         claims inserted bases are reference matches."
    );
}
