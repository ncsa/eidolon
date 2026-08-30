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

/// Same as `ratio_vs_control`, averaged over `n` independent seeds. A single seed at
/// this fixture's scale (H1N1, 1200bp event, 400bp margin) is noisy enough to flip a
/// borderline case: a 3-seed spread of [1.259, 1.347, 1.377] measured for dup_het
/// during the release/fragment-placement investigation (2026-08-22) means any ONE of
/// those draws alone could read anywhere from "fine" to "clearly off" depending on
/// which one the test happens to land on. Averaging is what makes the tolerance below
/// mean what it says instead of depending on luck.
fn ratio_vs_control_averaged(dir: &Path, label: &str, record: &str, replicates: usize) -> f64 {
    let mut ratios = Vec::with_capacity(replicates);
    for i in 0..replicates {
        let control = run_cell(dir, &format!("{label}_ctl_{i}"), &[]);
        let r = ratio_vs_control(dir, &format!("{label}_var_{i}"), &[record], &control);
        ratios.push(r);
    }
    ratios.iter().sum::<f64>() / ratios.len() as f64
}

#[test]
fn copy_number_events_scale_depth_as_declared() {
    let tmp = tempfile::tempdir().unwrap();

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
        let got = if expected == 0.0 {
            // Degenerate case: a single replicate is enough, the multiplier is 0 either
            // way and there's no boundary-scale effect to average out.
            let control = run_cell(tmp.path(), &format!("{label}_ctl"), &[]);
            ratio_vs_control(tmp.path(), label, &[record], &control)
        } else {
            ratio_vs_control_averaged(tmp.path(), label, record, 3)
        };
        report.push_str(&format!(
            "  {label:9} expected {expected:.2}  got {got:.2}\n"
        ));
        // Tolerance covers Poisson noise plus a real, understood, SCALE-DEPENDENT
        // boundary effect -- not one mechanism but (at least) two, both vanishing at
        // realistic event/genome size, both irrelevant to genome-scale campaigns:
        //
        //   1. Chimeric junction reads land outside the coverage-multiplied budget
        //      (#499): excess scales ~1/event_length, ~8% at this fixture's 1200bp,
        //      confirmed negligible (~0.3%) at 100kb on a real chr22 window.
        //   2. Fragment placement (release/fragment-placement, 2026-08-22) now lets a
        //      fragment's end extend past its owning coverage-multiplier segment into
        //      a differently-multiplied neighbor -- correct and necessary (it is what
        //      removes an artificial dead zone at every such boundary), but it costs
        //      some of the segment's own declared depth to redistribution. The size
        //      the ~20% tolerance needs to cover for THIS fixture specifically: at
        //      1200bp on H1N1 (narrow ~400bp flanks), correction needed ~1.13; the
        //      SAME 1200bp event on a real chr22 window with megabase-scale flanks
        //      needs only ~1.07; at a realistic 100kb event size, ~1.003 (see
        //      `depth_modulation_is_accurate_at_realistic_scale` below, and
        //      docs/claude_engineering_audit.md SS5.6's 2026-08-22 addendum).
        //
        // H1N1 cannot host an event "much larger than the fragment length" AND leave
        // wide flanks at the same time -- the contig is 2280bp. This test's job is
        // "the multiplier is applied and roughly right", fast; it is not, and cannot
        // be, a precision check at this scale.
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

/// Same as `mean_depth`, for a BAM whose reference has exactly one contig -- avoids
/// depending on the H1N1 fixture's `contigs()` index, which a standalone synthetic
/// reference (see `depth_modulation_is_accurate_at_realistic_scale`) is not part of.
fn mean_depth_single_contig(out: &Path, start: usize, end: usize) -> f64 {
    let file = std::fs::File::open(out.join("o.bam")).expect("bam produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let mut depth = vec![0u32; end - start + 1];
    for rec in reader.records() {
        let rec = rec.unwrap();
        if rec.reference_sequence_id().is_none() {
            continue;
        }
        let Some(Ok(pos)) = rec.alignment_start() else {
            continue;
        };
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
        if e < start || s > end {
            continue;
        }
        for p in s.max(start)..=e.min(end) {
            depth[p - start] += 1;
        }
    }
    depth.iter().map(|&d| d as f64).sum::<f64>() / depth.len() as f64
}

/// Deterministic pseudo-random ACGT sequence, same generator as `synthetic_insert`
/// below, used here to build a whole standalone reference rather than an insert.
fn synthetic_sequence(n: usize) -> String {
    let mut s = String::with_capacity(n);
    let mut x: u64 = 0x1234_5678_9ABC_DEF0;
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

fn write_fasta(path: &Path, contig: &str, seq: &str) {
    let mut out = String::with_capacity(seq.len() + seq.len() / 60 + 16);
    out.push('>');
    out.push_str(contig);
    out.push('\n');
    for chunk in seq.as_bytes().chunks(60) {
        out.push_str(std::str::from_utf8(chunk).unwrap());
        out.push('\n');
    }
    std::fs::write(path, out).unwrap();
}

/// REGRESSION GUARD for the assembly-gap boundary. `regions_of_interest` is built
/// from `get_non_n_regions()`, so N bases are deliberately excluded from read
/// generation -- but fragment placement is allowed to extend a fragment's END past
/// its owning region's edge (that is what removes the artificial dead zone at chunk
/// and coverage-multiplier boundaries). An assembly gap is the one boundary that
/// extension must treat as a true terminus, or reads carry fabricated gap sequence.
///
/// Caught in review of the fragment-placement branch, then measured: 103 of 6000
/// reads contained `N` on a reference with a 2 kb gap, against 0 both before that
/// branch and at its step 1 (before extension was wired into runner.rs) -- so it
/// was introduced by the extension wiring specifically, and is not pre-existing.
/// Nothing in the suite covered N-gaps at all, which is why it reached a PR;
/// `docs/sv_polish_roadmap.md`'s Phase 1 item 3 had already flagged that blind spot.
#[test]
fn fragments_do_not_extend_across_an_assembly_gap() {
    let tmp = tempfile::tempdir().unwrap();
    let (left, gap, right) = (5_000usize, 2_000usize, 5_000usize);
    let contig_len = left + gap + right;

    // [0, 5000) real | [5000, 7000) N | [7000, 12000) real
    let mut seq = synthetic_sequence(left);
    seq.push_str(&"N".repeat(gap));
    seq.push_str(&synthetic_sequence(right));
    let reference = tmp.path().join("ngap.fa");
    write_fasta(&reference, "ngap1", &seq);

    let out = tmp.path().join("run");
    std::fs::create_dir_all(&out).unwrap();
    let cfg = out.join("c.yml");
    std::fs::write(
        &cfg,
        format!(
            "reference: {ref}\nread_len: 100\ncoverage: 60\nploidy: 2\npaired_ended: true\n\
             fragment_mean: 250\nfragment_st_dev: 30\n\
             produce_vcf: false\nproduce_fastq: false\nproduce_bam: true\n\
             sv_rate_scale: 0.0\nmutation_rate: 0.0\noverwrite_output: true\n\
             output_dir: {out}\noutput_filename: o\nrng_seed: ngap guard\nnum_threads: 1\n",
            ref = reference.display(),
            out = out.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();

    // Read the BAM directly rather than the FASTQ: a placed read's reference span is
    // the thing that must not enter the gap, and grepping sequence for 'N' would also
    // be satisfied by a sequencing-error N, which is a different mechanism.
    let file = std::fs::File::open(out.join("o.bam")).expect("bam produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    let (mut total, mut into_gap) = (0usize, 0usize);
    for rec in reader.records() {
        let rec = rec.unwrap();
        let Some(Ok(pos)) = rec.alignment_start() else {
            continue;
        };
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
        total += 1;
        // 1-based inclusive read span vs the 1-based gap [left+1, left+gap].
        let (rs, re) = (pos.get(), pos.get() + span.saturating_sub(1));
        if rs <= left + gap && re > left {
            into_gap += 1;
        }
    }

    assert!(
        total > 100,
        "fixture produced too few reads ({total}) to be meaningful"
    );
    assert_eq!(
        into_gap,
        0,
        "{into_gap} of {total} reads overlap the assembly gap at [{}, {}] -- fragment \
         extension crossed an N-region boundary, so those reads carry fabricated gap \
         sequence. Extension must stop at a gap even though it correctly crosses chunk \
         and coverage-multiplier boundaries.",
        left + 1,
        left + gap
    );
    // Must-not-fire half: the gap must not have suppressed generation elsewhere.
    assert!(
        contig_len > 0 && total >= 4_000,
        "expected roughly full coverage of the two real intervals, got {total} reads"
    );
}

/// CONFIRMS the size-dependence recorded above and in
/// docs/claude_engineering_audit.md §5.6 (2026-08-22 addendum) at the scale that
/// actually matters: H1N1 cannot host a "much larger than fragment length" event with
/// wide flanks at the same time (it is 2280bp total), so `copy_number_events_scale_
/// depth_as_declared` can only ever be a fast "roughly right" mechanism check, never a
/// precision one. This test is the precision check, on a standalone 1Mb synthetic
/// reference (generated here, not checked in) with a 100kb event and wide flanks on
/// both sides -- matching the `docs/sv_polish_roadmap.md` item 1 size sweep
/// (1kb/100kb/1Mb) this investigation finally ran. Measured during that investigation:
/// correction needed was ~1.13 at H1N1 scale (1200bp event, ~400bp flanks), ~1.07 for
/// the SAME 1200bp event with megabase flanks (a real chr22 window), and ~1.003 at
/// 100kb with megabase flanks -- confirming the effect is real, understood (chimeric
/// junction reads outside the coverage-multiplied budget, plus fragments legitimately
/// extending across a coverage-multiplier boundary), and irrelevant at genome scale.
#[test]
fn depth_modulation_is_accurate_at_realistic_scale() {
    let tmp = tempfile::tempdir().unwrap();
    let contig_len = 1_000_000usize;
    let anchor = 400_000usize; // 1-based POS
    let svlen = 100_000usize;
    let end = anchor + svlen;
    let margin = 1_000usize;

    let reference = tmp.path().join("synthetic.fa");
    write_fasta(&reference, "synth1", &synthetic_sequence(contig_len));

    let run = |name: &str, records: &[String]| -> PathBuf {
        let out = tmp.path().join(name);
        std::fs::create_dir_all(&out).unwrap();
        let vcf = out.join("in.vcf");
        let mut text = format!(
            "##fileformat=VCFv4.2\n##contig=<ID=synth1,length={contig_len}>\n\
             ##ALT=<ID=DUP,Description=\"dup\">\n\
             ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n\
             ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"x\">\n\
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
                 produce_vcf: false\nproduce_fastq: false\nproduce_bam: true\n\
                 sv_rate_scale: 0.0\noverwrite_output: true\n\
                 output_dir: {out}\noutput_filename: o\nrng_seed: realistic {name}\n\
                 num_threads: 4\n",
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
    };

    let dup_record = format!(
        "synth1\t{anchor}\tv\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END={end};SVLEN={svlen}\tGT\t0/1"
    );
    let control = run("control", &[]);
    let dup = run("dup", &[dup_record]);

    let a = mean_depth_single_contig(&dup, anchor + margin, end - margin);
    let b = mean_depth_single_contig(&control, anchor + margin, end - margin);
    assert!(
        b > 10.0,
        "control depth {b} implausibly low — fixture broken"
    );
    let got = a / b;
    eprintln!("[realistic-scale] 100kb het DUP delivered {got:.3} against declared 1.50");
    assert!(
        (got / 1.5 - 1.0).abs() < 0.10,
        "100kb het DUP with wide flanks delivered {got:.3}, expected ~1.50 within 10% -- \
         if this is failing, the boundary effect documented above no longer vanishes at \
         realistic scale, which would be a real regression worth investigating"
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

    // FIXED, and now guarded. This block used to be a CHARACTERIZATION of #516,
    // asserting the inverse — that beyond ~a read length an insertion's interior and
    // far end were never sequenced while the truth VCF kept declaring the full SVLEN.
    // It is why campaign 20925151 measured Manta INS recall at 1/22: only the 61bp
    // event had evidence for its declared length, and Manta called only that one.
    //
    // Long insertions are now sampled in altered-haplotype coordinates, where the
    // inserted sequence HAS width and a fragment can begin inside it, so the whole
    // event reaches the reads. Flipped per this block's own former instructions
    // ("promote the size into the guards above"). The properties that make the fix
    // CORRECT rather than merely present — zygosity, allelic depth, neighbouring
    // variants, BAM validity, read-name uniqueness, fragment-length robustness —
    // are pinned separately in `eidolon/tests/long_insertion_rework.rs`, because the
    // previous attempt satisfied this cell while breaking all six of them.
    for &(size, _, m, t) in &carried {
        if size >= 200 {
            assert!(
                m > 0 && t > 0,
                "{size}bp insertion no longer reaches the reads mid-sequence or at its \
                 tail (middle={m}, tail={t}) — this is a REGRESSION of #516, not the old \
                 known gap."
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
fn symbolic_ins_is_realized_with_novel_sequence() {
    // CHARACTERIZATION of #500. `<INS>` carries no sequence, so nothing is spliced —
    // but the truth VCF still declares a 60bp insertion, so a benchmark built from it
    // Was a CHARACTERIZATION of #500, now the regression guard for its fix: a symbolic
    // <INS> is realized with synthesised novel sequence, so it reaches the reads instead of
    // being preserved in the truth VCF while the reads match a no-variant control.
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
        got > base + 2,
        "symbolic <INS> produced {got} reads with an I op against a background of {base} — \
         indistinguishable from the no-variant control, so the record is a silent no-op and \
         a benchmark built from this truth scores a caller as having missed an insertion \
         that was never in the data (#500)"
    );
}

#[test]
fn a_single_breakend_is_rejected_and_leaves_coverage_intact() {
    // Was a CHARACTERIZATION of #500, now the regression guard for its fix. An
    // unresolved-partner breakend is valid VCF 4.2, and eidolon cannot build a junction read
    // without a partner. It used to be accepted anyway: no junction reads AND ~40% of local
    // depth silently removed, so a depth caller saw a partial deletion no record described.
    // It is now rejected at input_vcf filtering, which both keeps it out of the truth VCF
    // and leaves coverage alone.
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
    assert!(
        !truth_has_id(&out, "sbnd"),
        "a single breakend reached the truth VCF — it declares a junction the reads cannot \
         contain, because there is no mate to build one from (#500)"
    );
    assert_eq!(
        chimeric_reads(&out),
        0,
        "a single breakend produced junction reads without a mate to join to"
    );
    // The point of rejecting it: coverage must be untouched. 0.85 was the old FAILING
    // threshold — the defect measured 0.59 — so this asserts the opposite of what it did.
    assert!(
        ratio > 0.85,
        "single breakend depth ratio {ratio:.2} is still depressed — the record is being \
         dropped from the truth but is still eating local coverage (#500)"
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
fn bnd_inserted_sequence_reaches_the_reads() {
    // Was a CHARACTERIZATION of #498, now the regression guard for its fix: the ALT carries
    // novel bases at the junction and the reads must contain them. Probe is deliberately NOT
    // a homopolymer, and is matched against sequence lines only — a G-run probe matched
    // quality strings and reported 71 false hits when the true answer was 0.
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
    assert!(
        hits > 0,
        "the junction insert appears in NO read. The truth VCF keeps the full ALT, so a \
         benchmark built from this data asserts an insertion the reads do not contain — \
         #498 has regressed."
    );
}
