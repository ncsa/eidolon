//! CI-enforced SV support matrix: what eidolon can reproduce from `input_vcf`.
//!
//! Companion to `docs/sv_support_matrix.md`. Every cell there is a test here.
//!
//! **Why this file exists.** Four input-fidelity defects (#474, #498, #499, #500) were
//! found by hand, and **not one of them is visible to a VCF-only round-trip check** —
//! every broken record survives the truth VCF perfectly. `bnd_roundtrip.rs` has four
//! assertions and none on content, which is why it caught none of them. The contract that
//! matters is *"the variant I supplied reached the reads"*, and only read-level evidence
//! can check it.
//!
//! **Method, and it is load-bearing.** Two earlier measurements were wrong because the
//! method was wrong, so it is fixed here:
//!
//! 1. **Compare against a separate no-variant control run**, never an in-run control
//!    span. Coverage varies enough along a contig that an in-run span misleads.
//! 2. **Use events much larger than the fragment length**, and measure the deep interior.
//!    A 300 bp inversion with 250 bp fragments sits entirely inside the junction dip zone
//!    and reads as a 37% coverage defect that does not exist.
//! 3. **Never grep a FASTQ without isolating the sequence line.** Phred+33 overlaps the
//!    DNA alphabet (`G` is Q38), so a homopolymer probe matches quality strings — that
//!    produced 71 false hits where the true answer was 0.
//!
//! Cells that are currently BROKEN are pinned as *characterization* tests: they assert
//! the wrong-but-current behaviour and say so. If one starts failing, the underlying bug
//! was probably fixed — flip the assertion and update the doc.

mod common;

use common::{eidolon, h1n1_reference};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const VCF_HEADER_TYPES: &[&str] = &["DEL", "DUP", "INV", "CNV", "INS", "BND"];

/// Contig name -> length, read from the H1N1 fixture.
fn contigs() -> &'static Vec<(String, usize)> {
    static C: OnceLock<Vec<(String, usize)>> = OnceLock::new();
    C.get_or_init(|| {
        let text = std::fs::read_to_string(h1n1_reference()).expect("H1N1 fixture");
        let mut out: Vec<(String, usize)> = Vec::new();
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix('>') {
                out.push((rest.split_whitespace().next().unwrap().to_string(), 0));
            } else if let Some(last) = out.last_mut() {
                last.1 += line.trim().len();
            }
        }
        out
    })
}

fn vcf_header() -> String {
    let mut h = String::from("##fileformat=VCFv4.2\n");
    for (name, len) in contigs() {
        h.push_str(&format!("##contig=<ID={name},length={len}>\n"));
    }
    for t in VCF_HEADER_TYPES {
        h.push_str(&format!("##ALT=<ID={t},Description=\"{t}\">\n"));
    }
    for (id, num, ty) in [
        ("SVTYPE", "1", "String"),
        ("END", "1", "Integer"),
        ("SVLEN", ".", "Integer"),
        ("CN", "1", "Integer"),
        ("MATEID", ".", "String"),
    ] {
        h.push_str(&format!(
            "##INFO=<ID={id},Number={num},Type={ty},Description=\"x\">\n"
        ));
    }
    h.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n");
    h
}

/// Run gen-reads over `records` (already tab-delimited VCF body lines) and return the
/// output directory. `records` empty => the no-variant control.
fn run_cell(dir: &Path, name: &str, records: &[&str]) -> PathBuf {
    let out = dir.join(name);
    std::fs::create_dir_all(&out).unwrap();
    let vcf = out.join("in.vcf");
    let mut text = vcf_header();
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
             output_dir: {out}\noutput_filename: o\nrng_seed: matrix {name}\nnum_threads: 1\n",
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
    out
}

/// Mean per-base depth over `[start, end]` (1-based inclusive), read from the BAM.
/// Deliberately reimplemented over noodles rather than shelling out: CI has no samtools.
fn mean_depth(out: &Path, contig: &str, start: usize, end: usize) -> f64 {
    let idx = contigs()
        .iter()
        .position(|(n, _)| n == contig)
        .expect("contig in fixture");
    let file = std::fs::File::open(out.join("o.bam")).expect("bam produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let mut depth = vec![0u32; end - start + 1];
    for rec in reader.records() {
        let rec = rec.unwrap();
        let Some(Ok(rid)) = rec.reference_sequence_id() else {
            continue;
        };
        if rid != idx {
            continue;
        }
        let Some(Ok(pos)) = rec.alignment_start() else {
            continue;
        };
        // Reference span: only ops that consume reference.
        let mut span = 0usize;
        for op in rec.cigar().iter() {
            let op = op.unwrap();
            if matches!(
                op.kind(),
                Kind::Match
                    | Kind::Deletion
                    | Kind::Skip
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
            ) {
                span += op.len();
            }
        }
        let (s, e) = (pos.get(), pos.get() + span.saturating_sub(1));
        for p in s.max(start)..=e.min(end) {
            depth[p - start] += 1;
        }
    }
    depth.iter().map(|&d| d as f64).sum::<f64>() / depth.len() as f64
}

/// Reads whose CIGAR contains the given op, within a reference window.
fn reads_with_op(out: &Path, contig: &str, start: usize, end: usize, kind: Kind) -> usize {
    let idx = contigs().iter().position(|(n, _)| n == contig).unwrap();
    let file = std::fs::File::open(out.join("o.bam")).unwrap();
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let mut n = 0;
    for rec in reader.records() {
        let rec = rec.unwrap();
        let (Some(Ok(rid)), Some(Ok(pos))) = (rec.reference_sequence_id(), rec.alignment_start())
        else {
            continue;
        };
        if rid != idx || pos.get() < start || pos.get() > end {
            continue;
        }
        if rec.cigar().iter().any(|o| o.unwrap().kind() == kind) {
            n += 1;
        }
    }
    n
}

/// Chimeric reads, counted from the FASTQ *name* lines only.
fn chimeric_reads(out: &Path) -> usize {
    let path = out.join("o_r1.fastq.gz");
    let f = std::fs::File::open(path).unwrap();
    let mut s = String::new();
    {
        use std::io::Read;
        flate2::read::MultiGzDecoder::new(f)
            .read_to_string(&mut s)
            .unwrap();
    }
    s.lines()
        .step_by(4)
        .filter(|l| l.contains("EIDOLON_chimeric"))
        .count()
}

/// Reads whose SEQUENCE contains `probe`. Counts FASTQ line 2 of each record only —
/// never the quality line. Phred+33 overlaps the DNA alphabet (`G` is Q38), so probing
/// a whole FASTQ matches quality strings: a G-run probe once reported 71 false hits
/// where the true answer was 0.
fn seq_lines_containing(out: &Path, probe: &str) -> usize {
    let mut n = 0;
    for f in ["o_r1.fastq.gz", "o_r2.fastq.gz"] {
        let Ok(file) = std::fs::File::open(out.join(f)) else {
            continue;
        };
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
            .filter(|l| l.contains(probe))
            .count();
    }
    n
}

/// Truth-VCF records carrying a given ID.
fn truth_has_id(out: &Path, id: &str) -> bool {
    let f = std::fs::File::open(out.join("o.vcf.gz")).unwrap();
    let mut s = String::new();
    {
        use std::io::Read;
        flate2::read::MultiGzDecoder::new(f)
            .read_to_string(&mut s)
            .unwrap();
    }
    s.lines()
        .filter(|l| !l.starts_with('#'))
        .any(|l| l.split('\t').nth(2) == Some(id))
}

// ── The matrix ──────────────────────────────────────────────────────────────────
// A 1200bp event on PB2 (2280bp), depth sampled over 900-1300 — the deep interior,
// far enough from both breakpoints that junction effects do not reach it.
const EV: (&str, usize, usize) = ("H1N1_PB2", 500, 1700);
const INTERIOR: (usize, usize) = (900, 1300);

fn rec(sv: &str, info: &str, gt: &str, id: &str) -> String {
    format!(
        "{}\t{}\t{id}\tA\t<{sv}>\t60\tPASS\t{info}\tGT\t{gt}",
        EV.0, EV.1
    )
}

/// Depth of a cell relative to the no-variant control, over the deep interior.
fn ratio_vs_control(dir: &Path, name: &str, records: &[&str], control: &Path) -> f64 {
    let out = run_cell(dir, name, records);
    let a = mean_depth(&out, EV.0, INTERIOR.0, INTERIOR.1);
    let b = mean_depth(control, EV.0, INTERIOR.0, INTERIOR.1);
    assert!(
        b > 10.0,
        "control depth {b} implausibly low — fixture broken"
    );
    a / b
}

#[test]
fn copy_number_events_scale_depth_as_declared() {
    let tmp = tempfile::tempdir().unwrap();
    let control = run_cell(tmp.path(), "control", &[]);

    // (label, record, expected multiplier)
    let del_het = rec("DEL", "SVTYPE=DEL;END=1700;SVLEN=-1200", "0/1", "v");
    let del_hom = rec("DEL", "SVTYPE=DEL;END=1700;SVLEN=-1200", "1/1", "v");
    let dup_het = rec("DUP", "SVTYPE=DUP;END=1700;SVLEN=1200", "0/1", "v");
    let cnv0 = rec("CNV", "SVTYPE=CNV;END=1700;SVLEN=1200;CN=0", "0/1", "v");
    let cnv4 = rec("CNV", "SVTYPE=CNV;END=1700;SVLEN=1200;CN=4", "0/1", "v");
    let inv = rec("INV", "SVTYPE=INV;END=1700;SVLEN=1200", "0/1", "v");

    let cases: Vec<(&str, &str, f64)> = vec![
        ("del_het", del_het.as_str(), 0.50),
        ("del_hom", del_hom.as_str(), 0.00),
        ("dup_het", dup_het.as_str(), 1.50),
        ("cnv_cn0", cnv0.as_str(), 0.00),
        ("cnv_cn4", cnv4.as_str(), 2.00),
        // A balanced inversion changes no copy number: the inverted sequence is still
        // present, so the deep interior must be at full depth.
        ("inv_het", inv.as_str(), 1.00),
    ];

    let mut report = String::new();
    let mut failures = Vec::new();
    for (label, record, expected) in cases {
        let got = ratio_vs_control(tmp.path(), label, &[record], &control);
        report.push_str(&format!(
            "  {label:9} expected {expected:.2}  got {got:.2}\n"
        ));
        // Tolerance covers Poisson noise AND the known ~8% over-delivery of #499. The
        // point of this test is that the multiplier is *applied and roughly right* —
        // #499 pins the 8% separately, so tightening here would duplicate it.
        let ok = if expected == 0.0 {
            got < 0.02
        } else {
            (got / expected - 1.0).abs() < 0.20
        };
        if !ok {
            failures.push(format!("{label}: expected ~{expected:.2}, got {got:.2}"));
        }
    }
    eprintln!("[copy-number matrix]\n{report}");
    assert!(
        failures.is_empty(),
        "depth does not match the declared copy number:\n  {}",
        failures.join("\n  ")
    );
}

#[test]
fn coverage_modulation_over_delivers_by_about_eight_percent() {
    // CHARACTERIZATION of #499. Only the degenerate multipliers (0 = emit nothing,
    // 1 = change nothing) are exact; everything requiring real scaling runs ~8% high.
    // If this test starts FAILING, #499 was probably fixed — verify, then tighten
    // `copy_number_events_scale_depth_as_declared` and delete this test.
    let tmp = tempfile::tempdir().unwrap();
    let control = run_cell(tmp.path(), "control", &[]);
    let cnv4 = rec("CNV", "SVTYPE=CNV;END=1700;SVLEN=1200;CN=4", "0/1", "v");
    let got = ratio_vs_control(tmp.path(), "cnv4_bias", &[cnv4.as_str()], &control);
    let bias = got / 2.00;
    eprintln!("[#499] CN=4 delivered {got:.3} against an expected 2.00 (bias {bias:.3})");
    assert!(
        bias > 1.03,
        "CN=4 delivered {got:.3} (bias {bias:.3}) — if this is now ~1.00, #499 is FIXED: \
         tighten copy_number_events_scale_depth_as_declared and remove this test"
    );
}

/// Deterministic, varied ACGT filler. A 30-mer of it has a ~4^-30 chance of occurring in the
/// reference by accident, so finding one in a read is proof the NOVEL sequence reached the
/// output rather than something spliced from the fixture.
fn synthetic_insert(n: usize) -> String {
    let mut s = String::with_capacity(n);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
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
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            o => o,
        })
        .collect()
}

/// Insertions LARGER than a read must still reach the reads — the size regime the rest of
/// this file never covered.
///
/// WHY THIS EXISTS: campaign 20925151 planted 22 de novo insertions of 61–2155 bp. Manta
/// reported nothing at any of the 21 above 61 bp — not even the `IMPRECISE <INS>` its own
/// user guide says it emits from a breakend signature for insertions it cannot fully
/// assemble, and its documented full-assembly ceiling (~2x fragment size) sits above most of
/// the planted set. That points at the reads rather than the caller, and the question was
/// being chased through 8-hour whole-genome runs when it is answerable here in seconds.
///
/// `literal_indels_reach_the_reads` above covers an ELEVEN base insert. A de novo INS is
/// emitted as the same `VariantType::Insertion` + `AlternateType::Literal` variant that an
/// `input_vcf` literal insertion parses to, so both share this machinery downstream — which
/// makes a size sweep here a valid probe of the de novo path.
///
/// KNOWN ANSWER: at coverage 60, a heterozygous insertion gets ~30x over the inserted
/// sequence, so reads falling wholly inside it must exist for any size above a read length.
/// The probe is a 30-mer from the MIDDLE of the insert, checked in both orientations.
#[test]
fn large_literal_insertions_reach_the_reads() {
    let tmp = tempfile::tempdir().unwrap();
    // MECHANISM (located, not inferred — see #516). Fragments and read windows are chosen
    // purely in REFERENCE offsets (`cover_dataset`, generate_fragments.rs), and an insertion
    // has zero reference width, so no read window can ever BEGIN inside one. Reads are then
    // assembled per-read by walking a reference slice and expanding variants inline, capped by
    // `bases_written < read_length` (fastq_tools.rs:429) with `break 'outer` at :555 dropping
    // the remainder of `ins_buf` silently. Hence the invariant:
    //
    //     novel bases visible in one read = min(L, read_length - anchor_offset - 1)
    //
    // so AT MOST read_length-1 = 99 inserted bases can be realized at ANY declared SVLEN.
    // Probe visibility follows: head needs anchor_offset <= 69 (no L term), middle needs
    // anchor_offset <= 84 - L/2 (impossible for L >= 170), tail needs <= 99 - L (impossible
    // for L >= 100).
    //
    // MEASURED, het (this test):     head 16/21/25/26/27   middle 13/5/0/0/0   tail 8/0/0/0/0
    // MEASURED, hom control (1/1):   head 42/41/50/50/50   middle 36/6/0/0/0   tail 29/0/0/0/0
    //
    // The head count is SIZE-INDEPENDENT above saturation — 50/50/50 for 200/300/600 with the
    // het coin removed. An earlier writeup of this bug claimed the head count "rises with
    // size" and read meaning into 16 -> 27; that was RNG drift (a larger insert fires
    // `break 'outer` sooner, consuming fewer sequencing-error draws and shifting the stream).
    // Corrected by rerunning with GT 1/1, which short-circuits the coin at fastq_tools.rs:471.
    let mut carried: Vec<(usize, usize, usize, usize)> = Vec::new();
    for size in [50usize, 150, 200, 300, 600] {
        let insert = synthetic_insert(size);
        let at = size / 2 - 15;
        let probe = &insert[at..at + 30];

        let rec = format!(
            "{}\t500\tinsbig{size}\tA\tA{insert}\t60\tPASS\t.\tGT\t0/1",
            EV.0
        );
        let out = run_cell(tmp.path(), &format!("ins_big_{size}"), &[rec.as_str()]);

        let hits = seq_lines_containing(&out, probe) + seq_lines_containing(&out, &revcomp(probe));
        // Probe the START and END as well as the middle. With 100bp reads from 250bp
        // fragments there is a ~50bp unsequenced gap inside every fragment, so a middle-only
        // probe could read zero for GEOMETRIC reasons. If the start appears and the middle
        // does not, the insert is spliced but only partly sequenced; if none appear, the
        // novel sequence never reached the output at all. Different bugs, different fixes.
        let head = &insert[..30];
        let tail = &insert[size - 30..];
        let h = seq_lines_containing(&out, head) + seq_lines_containing(&out, &revcomp(head));
        let t = seq_lines_containing(&out, tail) + seq_lines_containing(&out, &revcomp(tail));
        eprintln!("[INS {size}bp] 30-mer hits — head={h} middle={hits} tail={t}");

        // The truth VCF preserves every size — which is the problem, not the reassurance.
        assert!(
            truth_has_id(&out, &format!("insbig{size}")),
            "{size}bp insertion lost from the truth VCF"
        );
        carried.push((size, h, hits, t));
    }

    // MEASURED, read_len=100 / fragment_mean=250:
    //   size  head middle tail
    //     50    16     13    8   fully realized
    //    150    21      5    0   middle thinning, tail already gone
    //    200    25      0    0   only the head survives
    //    300    26      0    0
    //    600    27      0    0
    //
    // REGRESSION GUARD. The head must ALWAYS reach the reads — that is what proves the
    // insertion is spliced into the haplotype at all, and it is the part that makes this bug
    // look like a working feature.
    for &(size, h, _, _) in &carried {
        assert!(
            h > 0,
            "{size}bp insertion: not even its first 30 bases reach the reads. This is worse \
             than the known partial-realization gap — the splice itself is broken."
        );
    }
    // Small insertions are fully realized, and must stay that way.
    for &(size, _, m, _) in &carried {
        if size <= 150 {
            assert!(
                m > 0,
                "{size}bp insertion no longer reaches the reads mid-sequence — a regression, \
                 not the known #516 gap"
            );
        }
    }

    // Production runner coverage: long literal insertions must reach both their
    // interior and their far end. The paired haplotype writer now replaces
    // regular reads touching these anchors and carries insertion-aware baseline
    // CIGAR operations into the BAM.
    for &(size, _, m, t) in &carried {
        if size >= 200 {
            assert!(
                m > 0 && t > 0,
                "{size}bp insertion does not reach both the interior and tail \
                 (middle={m}, tail={t})"
            );
        }
    }
}

#[test]
fn literal_indels_reach_the_reads() {
    let tmp = tempfile::tempdir().unwrap();
    let control = run_cell(tmp.path(), "control", &[]);
    let base_i = reads_with_op(&control, EV.0, 450, 560, Kind::Insertion);
    let base_d = reads_with_op(&control, EV.0, 450, 560, Kind::Deletion);

    let ins = format!(
        "{}\t500\tinslit\tA\tACTGACTGCATG\t60\tPASS\t.\tGT\t0/1",
        EV.0
    );
    let out_i = run_cell(tmp.path(), "ins_literal", &[ins.as_str()]);

    // CONTENT is the real assertion: the inserted bases must appear in read sequence.
    // A CIGAR-I count alone is a proxy — sequencing-error indels raise it too — so the
    // probe is the decisive check and the count is corroboration.
    let hits = seq_lines_containing(&out_i, "CTGACTGCATG");
    assert!(
        hits > 0,
        "no read carries the inserted bases CTGACTGCATG — a literal insertion supplied \
         via input_vcf is not reaching the reads"
    );
    let got_i = reads_with_op(&out_i, EV.0, 450, 560, Kind::Insertion);
    assert!(
        got_i > base_i,
        "literal insertion produced {got_i} reads with a CIGAR I against a background of \
         {base_i} — expected a clear rise"
    );
    assert!(
        truth_has_id(&out_i, "inslit"),
        "record lost from the truth VCF"
    );

    // Literal deletion: REF spans the deleted bases, ALT is the anchor alone. The REF
    // segment is read from the fixture so the record is consistent with the reference —
    // a mismatched REF is a different failure and would confound this cell.
    let refseg = {
        let text = std::fs::read_to_string(h1n1_reference()).unwrap();
        let mut seq = String::new();
        let mut on = false;
        for line in text.lines() {
            if let Some(r) = line.strip_prefix('>') {
                on = r.split_whitespace().next() == Some(EV.0);
            } else if on {
                seq.push_str(line.trim());
            }
        }
        seq[499..511].to_string()
    };
    let del = format!("{}\t500\tdellit\t{refseg}\tA\t60\tPASS\t.\tGT\t0/1", EV.0);
    let out_d = run_cell(tmp.path(), "del_literal", &[del.as_str()]);
    let got_d = reads_with_op(&out_d, EV.0, 450, 560, Kind::Deletion);
    assert!(
        got_d > base_d,
        "literal deletion produced {got_d} reads with a CIGAR D against a background of \
         {base_d} — the deletion is not reaching the reads"
    );
    assert!(
        truth_has_id(&out_d, "dellit"),
        "record lost from the truth VCF"
    );
}

#[test]
fn symbolic_ins_is_currently_a_silent_no_op() {
    // CHARACTERIZATION of #500. `<INS>` carries no sequence, so nothing is spliced —
    // but the truth VCF still declares a 60bp insertion, so a benchmark built from it
    // scores a caller as having missed something that was never in the data.
    // If this FAILS, symbolic <INS> is now realized: flip it to a positive assertion.
    let tmp = tempfile::tempdir().unwrap();
    let control = run_cell(tmp.path(), "control", &[]);
    let base = reads_with_op(&control, EV.0, 450, 560, Kind::Insertion);

    let r = rec("INS", "SVTYPE=INS;SVLEN=60", "0/1", "inssym");
    let out = run_cell(tmp.path(), "ins_symbolic", &[r.as_str()]);
    assert!(
        truth_has_id(&out, "inssym"),
        "record lost from the truth VCF"
    );
    let got = reads_with_op(&out, EV.0, 450, 560, Kind::Insertion);
    eprintln!("[#500] symbolic <INS>: {got} reads with CIGAR I vs {base} in the control");
    assert!(
        got <= base + 2,
        "symbolic <INS> produced {got} insertions against a background of {base} — it may \
         now be realized. If so, #500 is FIXED: flip this to assert the insertion appears"
    );
}

#[test]
fn a_single_breakend_currently_destroys_coverage() {
    // CHARACTERIZATION of #500. An unresolved-partner breakend is valid VCF 4.2. It
    // produces no junction reads AND silently removes ~40% of local depth, so a depth
    // caller sees a partial deletion that no record describes.
    let tmp = tempfile::tempdir().unwrap();
    let control = run_cell(tmp.path(), "control", &[]);
    let r = format!("{}\t500\tsbnd\tA\tA.\t60\tPASS\tSVTYPE=BND\tGT\t0/1", EV.0);
    let out = run_cell(tmp.path(), "bnd_single", &[r.as_str()]);
    let got = mean_depth(&out, EV.0, 495, 515);
    let base = mean_depth(&control, EV.0, 495, 515);
    let ratio = got / base;
    eprintln!(
        "[#500] single breakend: depth ratio {ratio:.2}, chimeric reads {}",
        chimeric_reads(&out)
    );
    assert_eq!(
        chimeric_reads(&out),
        0,
        "a single breakend now makes junction reads — #500 may be fixed"
    );
    assert!(
        ratio < 0.85,
        "single breakend depth ratio {ratio:.2} is no longer depressed — #500 may be \
         FIXED: verify and replace this with the intended behaviour"
    );
}

#[test]
fn inter_chromosomal_breakends_produce_junction_reads() {
    let tmp = tempfile::tempdir().unwrap();
    let a =
        format!("H1N1_HA\t600\tb1\tA\tA[H1N1_PB2:900[\t60\tPASS\tSVTYPE=BND;MATEID=b2\tGT\t0/1");
    let b =
        format!("H1N1_PB2\t900\tb2\tC\t]H1N1_HA:600]C\t60\tPASS\tSVTYPE=BND;MATEID=b1\tGT\t0/1");
    let out = run_cell(tmp.path(), "bnd_inter", &[a.as_str(), b.as_str()]);
    assert!(
        truth_has_id(&out, "b1") && truth_has_id(&out, "b2"),
        "BND pair lost from truth"
    );
    let chim = chimeric_reads(&out);
    assert!(
        chim >= 5,
        "inter-chromosomal breakend produced only {chim} chimeric read(s)"
    );
}

#[test]
fn bnd_inserted_sequence_is_currently_dropped_from_reads() {
    // CHARACTERIZATION of #498. The ALT carries novel bases at the junction; the truth
    // VCF keeps them and the reads do not. Probe is deliberately NOT a homopolymer, and
    // is matched against sequence lines only — a G-run probe matched quality strings and
    // reported 71 false hits when the true answer was 0.
    let tmp = tempfile::tempdir().unwrap();
    let a =
        "H1N1_HA\t600\td1\tA\tACTGACTGCATG[H1N1_PB2:900[\t60\tPASS\tSVTYPE=BND;MATEID=d2\tGT\t0/1";
    let b =
        "H1N1_PB2\t900\td2\tC\t]H1N1_HA:600]CTGACTGCATGC\t60\tPASS\tSVTYPE=BND;MATEID=d1\tGT\t0/1";
    let out = run_cell(tmp.path(), "bnd_insert", &[a, b]);

    let f = std::fs::File::open(out.join("o_r1.fastq.gz")).unwrap();
    let mut s = String::new();
    {
        use std::io::Read;
        flate2::read::MultiGzDecoder::new(f)
            .read_to_string(&mut s)
            .unwrap();
    }
    // Sequence lines only: lines 2, 6, 10, ... of the FASTQ.
    let hits = s
        .lines()
        .skip(1)
        .step_by(4)
        .filter(|l| l.contains("CTGACTGCATG"))
        .count();
    eprintln!("[#498] reads carrying the junction insert: {hits}");
    assert_eq!(
        hits, 0,
        "the junction insert now appears in {hits} read(s) — #498 may be FIXED: flip \
         this to assert the insert IS present"
    );
}
