//! **Gate 2b** — does the golden BAM agree with the FASTQ it was generated alongside, and with the
//! reference it claims to be aligned to?
//!
//! See `docs/sv_polish_roadmap.md`. Gate 2a asks whether a *realigned* BAM carries the evidence a
//! caller needs, and requires an aligner. This gate asks a different and cheaper question: are
//! eidolon's own two outputs mutually consistent? They are produced by one generation pass but
//! written by different code (`fastq_tools.rs`, `bam_writer.rs`), and nothing asserted they agree
//! — the shape that let `sv_model.rs` and `runner.rs` disagree about BND geometry for eight
//! releases.
//!
//! **No aligner, no variants, no SVs required, so this runs in CI on every push.** That makes it a
//! stronger guarantee than Gate 2a, which is `#[ignore]`d and only runs when someone remembers.
//!
//! ## What it pins
//!
//! | assertion | defect it would have caught |
//! |---|---|
//! | `SEQ` matches the reference at `POS` on **both** strands | [#550](https://github.com/ncsa/eidolon/issues/550) |
//! | reverse-flagged `SEQ`/`QUAL` are the revcomp/reverse of the FASTQ's | #550 |
//! | no QNAME ends in `/1` or `/2`, and mates share a name | [#551](https://github.com/ncsa/eidolon/issues/551) |
//!
//! The reference check is the independent one: it needs neither the FASTQ nor any assumption about
//! orientation conventions. `SEQ` either matches the sequence it is aligned to or it does not, and
//! #550 was 0.14 identity against 1.00 after reverse-complementing.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;
use std::path::Path;

fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            other => other,
        })
        .collect()
}

/// Whole-contig sequences from the H1N1 fixture. The fixture is CRLF, so `\r` is stripped —
/// a detail that has already cost one debugging session via `samtools faidx`.
fn reference_contigs() -> HashMap<String, String> {
    let text = std::fs::read_to_string(h1n1_reference()).unwrap();
    let mut out: HashMap<String, String> = HashMap::new();
    let mut current = String::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(header) = line.strip_prefix('>') {
            current = header.split_whitespace().next().unwrap_or("").to_string();
            out.insert(current.clone(), String::new());
        } else if !current.is_empty() {
            out.get_mut(&current).unwrap().push_str(line);
        }
    }
    out
}

/// FASTQ name (including its `/1` or `/2`) -> (sequence, quality).
fn read_fastq(path: &Path) -> HashMap<String, (String, String)> {
    let lines = read_gzip_fastq_lines(path);
    let mut out = HashMap::new();
    for chunk in lines.chunks(4) {
        if chunk.len() < 4 {
            break;
        }
        let name = chunk[0].trim_start_matches('@').trim().to_string();
        out.insert(name, (chunk[1].clone(), chunk[3].clone()));
    }
    out
}

struct Run {
    bam: std::path::PathBuf,
    r1: HashMap<String, (String, String)>,
    r2: HashMap<String, (String, String)>,
}

fn generate(paired: bool) -> (tempfile::TempDir, Run) {
    let (dir, work) = fresh_workdir();
    let tag = if paired { "paired" } else { "single" };
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), tag);
    config.coverage = 5;
    config.read_len = 100;
    config.paired_ended = paired;
    config.produce_fastq = true;
    config.produce_bam = true;
    config.produce_vcf = false;
    // No variants at all. This gate is about the two writers agreeing, so the reads should be
    // plain reference sequence plus the sequencing-error model — nothing else to explain a
    // mismatch away with.
    config.mutation_rate = Some(0.0);
    config.sv_rate_scale = Some(0.0);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let r1 = read_fastq(&work.join(format!("{tag}_r1.fastq.gz")));
    let r2 = if paired {
        read_fastq(&work.join(format!("{tag}_r2.fastq.gz")))
    } else {
        HashMap::new()
    };
    let bam = work.join(format!("{tag}.bam"));
    assert!(bam.is_file(), "golden BAM not produced at {bam:?}");
    (dir, Run { bam, r1, r2 })
}

/// `SEQ` must match the reference at `POS` — on BOTH strands.
///
/// This is the assertion that needs no FASTQ and no convention: a record either agrees with the
/// sequence it declares itself aligned to, or it does not. #550 failed it at 0.14 identity for
/// every reverse-strand record while forward records sat at ~0.99.
#[test]
fn golden_bam_sequences_match_the_reference_on_both_strands() {
    let (_dir, run) = generate(true);
    let contigs = reference_contigs();

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&run.bam)
        .unwrap();
    let header = reader.read_header().unwrap();

    let mut checked = (0usize, 0usize); // (forward, reverse)
    let mut worst = (1.0f64, String::new());
    for result in reader.records() {
        let record = result.unwrap();
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        let ref_id = record.reference_sequence_id().transpose().unwrap().unwrap();
        let ref_name = header
            .reference_sequences()
            .get_index(ref_id)
            .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).to_string())
            .unwrap();
        let seq: String = record
            .sequence()
            .iter()
            .map(|b| char::from(b).to_ascii_uppercase())
            .collect();
        let contig = &contigs[&ref_name];

        // CIGAR-aware comparison. The sequencing-error model emits indels, so a positional
        // comparison frameshifts on the first `I`/`D` and reports a correct read as ~0.27
        // identity — which is how the first version of this assertion failed on good data.
        // Only M/=/X runs are compared; I/S consume read bases, D/N consume reference.
        let mut read_pos = 0usize;
        let mut ref_pos = usize::from(start) - 1;
        let mut matched = 0usize;
        let mut compared_bases = 0usize;
        for op in record.cigar().iter() {
            let op = op.unwrap();
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for k in 0..len {
                        if let (Some(r), Some(q)) = (
                            contig.as_bytes().get(ref_pos + k),
                            seq.as_bytes().get(read_pos + k),
                        ) {
                            compared_bases += 1;
                            if r.eq_ignore_ascii_case(q) {
                                matched += 1;
                            }
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                }
                Kind::Insertion | Kind::SoftClip => read_pos += len,
                Kind::Deletion | Kind::Skip => ref_pos += len,
                _ => {}
            }
        }
        if compared_bases == 0 {
            continue;
        }
        let expected = contig
            [(usize::from(start) - 1)..(usize::from(start) - 1 + seq.len()).min(contig.len())]
            .to_uppercase();

        let id = matched as f64 / compared_bases as f64;
        let reverse = record.flags().is_reverse_complemented();
        if reverse {
            checked.1 += 1;
        } else {
            checked.0 += 1;
        }
        if id < worst.0 {
            worst = (
                id,
                format!(
                    "{ref_name}:{} reverse={reverse}\n    SEQ {}\n    ref {}\n    rc  {}",
                    usize::from(start),
                    &seq[..seq.len().min(50)],
                    &expected[..expected.len().min(50)],
                    &revcomp(&seq)[..seq.len().min(50)],
                ),
            );
        }
    }

    // Both strands must be represented, or the assertion could pass by never testing one.
    assert!(
        checked.0 > 20 && checked.1 > 20,
        "expected both strands to be well represented; got {} forward and {} reverse",
        checked.0,
        checked.1
    );
    // 0.90 leaves room for the sequencing-error model (~0.99 in practice) while sitting far above
    // a reverse-complement mismatch, which lands near 0.14 — the two are not close.
    assert!(
        worst.0 > 0.90,
        "a golden BAM record disagrees with the reference it is aligned to (identity {:.2}).\n  \
         {}\n  If `rc` matches `ref` and `SEQ` does not, SEQ is in READ orientation and SAM 1.4 \
         requires REFERENCE orientation when 0x10 is set (#550).",
        worst.0,
        worst.1
    );
}

/// The golden BAM and the FASTQ describe the same reads, so they must agree exactly — modulo the
/// orientation convention that differs between the two formats.
#[test]
fn golden_bam_and_fastq_describe_the_same_reads() {
    let (_dir, run) = generate(true);
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&run.bam)
        .unwrap();
    let _ = reader.read_header().unwrap();

    let mut compared = 0usize;
    for result in reader.records() {
        let record = result.unwrap();
        let flags = record.flags();
        let qname = String::from_utf8_lossy(record.name().unwrap().as_ref()).to_string();
        let mate_suffix = if flags.is_last_segment() { "/2" } else { "/1" };
        let source = if flags.is_last_segment() {
            &run.r2
        } else {
            &run.r1
        };
        let key = format!("{qname}{mate_suffix}");
        let Some((fq_seq, fq_qual)) = source.get(&key) else {
            panic!(
                "golden BAM record {qname:?} (mate {mate_suffix}) has no FASTQ counterpart at \
                 key {key:?}. Either the two writers disagree about read names, or QNAME still \
                 carries its own /1 or /2 suffix (#551)."
            );
        };
        let seq: String = record
            .sequence()
            .iter()
            .map(|b| char::from(b).to_ascii_uppercase())
            .collect();
        let qual: String = record
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&q| char::from(q + 33))
            .collect();

        let (want_seq, want_qual) = if flags.is_reverse_complemented() {
            (revcomp(fq_seq), fq_qual.chars().rev().collect::<String>())
        } else {
            (fq_seq.clone(), fq_qual.clone())
        };
        assert_eq!(
            seq,
            want_seq,
            "SEQ mismatch for {qname} (reverse={}): the BAM and FASTQ disagree about the bases \
             of the same read",
            flags.is_reverse_complemented()
        );
        assert_eq!(
            qual,
            want_qual,
            "QUAL mismatch for {qname} (reverse={}): quality must be reversed exactly when the \
             sequence is reverse-complemented",
            flags.is_reverse_complemented()
        );
        compared += 1;
    }
    assert!(
        compared > 100,
        "only {compared} record(s) compared — too few for this to mean anything"
    );
}

/// SAM 1.4 field 1: mates share a QNAME. The `/1` and `/2` suffixes are a FASTQ convention and
/// must not survive into the BAM, or nothing can pair the segments (#551).
#[test]
fn golden_bam_qnames_are_shared_between_mates() {
    let (_dir, run) = generate(true);
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&run.bam)
        .unwrap();
    let _ = reader.read_header().unwrap();

    let mut seen: HashMap<String, (usize, usize)> = HashMap::new(); // name -> (first, last)
    let mut records = 0usize;
    for result in reader.records() {
        let record = result.unwrap();
        let flags = record.flags();
        let qname = String::from_utf8_lossy(record.name().unwrap().as_ref()).to_string();
        assert!(
            !qname.ends_with("/1") && !qname.ends_with("/2"),
            "QNAME {qname:?} carries a /1 or /2 suffix. SAM 1.4 field 1: segments with the same \
             QNAME are one template, and mate identity is carried by FLAG 0x40/0x80. With the \
             suffix, mates have DIFFERENT names and nothing can pair them (#551)"
        );
        let e = seen.entry(qname).or_insert((0, 0));
        if flags.is_last_segment() {
            e.1 += 1;
        } else {
            e.0 += 1;
        }
        records += 1;
    }

    assert!(records > 100, "only {records} record(s); too few to judge");
    assert_eq!(
        seen.len() * 2,
        records,
        "a paired run must have exactly half as many distinct QNAMEs as records; got {} names \
         for {records} records",
        seen.len()
    );
    for (name, (first, last)) in &seen {
        assert_eq!(
            (*first, *last),
            (1, 1),
            "QNAME {name:?} should appear once as FIRST_SEGMENT and once as LAST_SEGMENT; got \
             first={first} last={last}"
        );
    }
}

/// MUST-NOT-FIRE: a single-ended run has no mates and no reverse-strand convention to apply, so
/// its records must be byte-identical to the FASTQ. Establishes that the comparison itself is
/// sound before any orientation logic is involved.
#[test]
fn single_ended_golden_bam_is_identical_to_its_fastq() {
    let (_dir, run) = generate(false);
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&run.bam)
        .unwrap();
    let _ = reader.read_header().unwrap();

    let mut compared = 0usize;
    for result in reader.records() {
        let record = result.unwrap();
        let qname = String::from_utf8_lossy(record.name().unwrap().as_ref()).to_string();
        let seq: String = record
            .sequence()
            .iter()
            .map(|b| char::from(b).to_ascii_uppercase())
            .collect();
        // Single-ended output names have no mate suffix, so look the name up as-is first.
        let entry = run
            .r1
            .get(&qname)
            .or_else(|| run.r1.get(&format!("{qname}/1")));
        let (fq_seq, _) = entry.unwrap_or_else(|| {
            panic!("single-ended golden BAM record {qname:?} has no FASTQ counterpart")
        });
        assert_eq!(
            &seq, fq_seq,
            "single-ended records carry no reverse-strand convention and must match the FASTQ \
             exactly"
        );
        compared += 1;
    }
    assert!(compared > 50, "only {compared} record(s) compared");
}
