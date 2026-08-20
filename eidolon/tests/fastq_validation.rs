//! Validate that generated FASTQ files are well-formed under a parser that is
//! *external* to eidolon's own writer. The point isn't to re-test the parser — it's to
//! catch any output-side regression (off-by-one in name encoding, wrong qual offset,
//! truncated trailing newline, etc.) by running the bytes through a strict reader.
//!
//! We use a small hand-rolled noodles-free parser here to keep dev-deps minimal. The
//! parser is intentionally strict: any deviation from the 4-line record format yields
//! a clear error rather than a silent skip.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};

/// One parsed FASTQ record. Owned strings — these tests are small enough.
#[derive(Debug)]
struct FastqRecord {
    name: String,
    sequence: String,
    quality: String,
}

/// Strict 4-line FASTQ parser. Returns an error string on any deviation rather than
/// panicking, so the caller can assert it didn't fire.
fn parse_fastq_strict(lines: &[String]) -> Result<Vec<FastqRecord>, String> {
    if !lines.len().is_multiple_of(4) {
        return Err(format!(
            "FASTQ has {} lines, not a multiple of 4 (truncated record at end?)",
            lines.len(),
        ));
    }
    let mut out = Vec::with_capacity(lines.len() / 4);
    for (i, chunk) in lines.chunks(4).enumerate() {
        if !chunk[0].starts_with('@') {
            return Err(format!(
                "record {i}: name line must start with '@', got {:?}",
                chunk[0]
            ));
        }
        if chunk[2] != "+" {
            return Err(format!(
                "record {i}: separator line must be '+', got {:?}",
                chunk[2]
            ));
        }
        if chunk[1].is_empty() {
            return Err(format!("record {i}: sequence is empty"));
        }
        if chunk[1].len() != chunk[3].len() {
            return Err(format!(
                "record {i}: seq length {} != qual length {}",
                chunk[1].len(),
                chunk[3].len(),
            ));
        }
        // Sequence must be ACGTN-only (eidolon doesn't emit lowercase or masked bases on output).
        if !chunk[1]
            .chars()
            .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'N'))
        {
            return Err(format!(
                "record {i}: sequence contains non-ACGTN characters: {:?}",
                chunk[1],
            ));
        }
        // Phred+33 qualities are printable ASCII 33..=126.
        if !chunk[3].bytes().all(|b| (33..=126).contains(&b)) {
            return Err(format!(
                "record {i}: quality contains non-printable bytes: {:?}",
                chunk[3],
            ));
        }
        out.push(FastqRecord {
            name: chunk[0].clone(),
            sequence: chunk[1].clone(),
            quality: chunk[3].clone(),
        });
    }
    Ok(out)
}

#[test]
fn single_ended_fastq_is_well_formed() {
    let (_dir, out) = fresh_workdir();
    let config = GenReadsConfig::new(h1n1_reference(), out.clone(), "validation_se");
    let yaml = config.write_yaml();

    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let path = out.join("validation_se_r1.fastq.gz");
    let lines = read_gzip_fastq_lines(&path);
    let records = parse_fastq_strict(&lines).expect("FASTQ failed strict parse");
    assert!(!records.is_empty(), "single-ended FASTQ has no records");

    // Every record's sequence length must equal the configured read_len.
    for r in &records {
        assert_eq!(
            r.sequence.len(),
            config.read_len,
            "read {:?} has sequence length {} but read_len is {}",
            r.name,
            r.sequence.len(),
            config.read_len,
        );
    }
}

#[test]
fn paired_ended_fastq_pair_invariants() {
    // Paired-end output must:
    //   * produce both R1 and R2 files,
    //   * have equal record counts in R1 and R2,
    //   * have read names that suffix /1 in R1 and /2 in R2,
    //   * have matching name *stems* between paired records at the same index.
    let (_dir, out) = fresh_workdir();
    let mut config = GenReadsConfig::new(h1n1_reference(), out.clone(), "validation_pe");
    config.paired_ended = true;
    let yaml = config.write_yaml();

    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let r1_path = out.join("validation_pe_r1.fastq.gz");
    let r2_path = out.join("validation_pe_r2.fastq.gz");
    assert!(
        r1_path.exists() && r2_path.exists(),
        "paired-end mode must produce both R1 and R2"
    );

    let r1 = parse_fastq_strict(&read_gzip_fastq_lines(&r1_path)).expect("R1 failed strict parse");
    let r2 = parse_fastq_strict(&read_gzip_fastq_lines(&r2_path)).expect("R2 failed strict parse");

    assert!(!r1.is_empty(), "R1 has no records");
    assert_eq!(
        r1.len(),
        r2.len(),
        "R1/R2 record counts differ: r1={} r2={}",
        r1.len(),
        r2.len(),
    );

    // Suffix and stem invariants. The /1 vs /2 suffix distinguishes mates; the stem
    // (everything before the suffix) is shared between paired reads.
    for (i, (a, b)) in r1.iter().zip(r2.iter()).enumerate() {
        assert!(
            a.name.ends_with("/1"),
            "R1[{i}] name {:?} should end with /1",
            a.name
        );
        assert!(
            b.name.ends_with("/2"),
            "R2[{i}] name {:?} should end with /2",
            b.name
        );
        let stem_a = a.name.trim_end_matches("/1");
        let stem_b = b.name.trim_end_matches("/2");
        assert_eq!(
            stem_a, stem_b,
            "paired records at index {i} have different name stems: {stem_a:?} vs {stem_b:?}",
        );
    }
}

#[test]
fn all_qualities_are_valid_phred33() {
    // The strict parser already asserts qualities are printable ASCII; this test makes
    // the Phred+33 mapping explicit: every quality character decodes to a Q-score in
    // [0, 93] (the legal range for the Phred+33 encoding).
    let (_dir, out) = fresh_workdir();
    let config = GenReadsConfig::new(h1n1_reference(), out.clone(), "validation_phred");
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    let records = parse_fastq_strict(&read_gzip_fastq_lines(
        &out.join("validation_phred_r1.fastq.gz"),
    ))
    .expect("FASTQ failed strict parse");

    // `q <= 93` cannot fail: parse_fastq_strict already rejects any byte outside
    // 33..=126, so `b - 33` is <= 93 BY CONSTRUCTION. The assertion was unfalsifiable
    // given its own precondition — switching the encoder to Phred+64, the classic
    // offset bug, passed it. Assert on the DISTRIBUTION instead, which is what actually
    // pins the offset: the same bytes decode to a mean near 65 under Phred+64.
    let mut sum = 0u64;
    let mut n = 0u64;
    let mut above_45 = 0u64;
    for r in &records {
        for &b in r.quality.as_bytes() {
            let q = b.saturating_sub(33) as u64;
            sum += q;
            n += 1;
            if q > 45 {
                above_45 += 1;
            }
        }
    }
    assert!(n > 0, "no quality bytes to check");
    let mean = sum as f64 / n as f64;
    assert!(
        (25.0..=42.0).contains(&mean),
        "mean decoded quality {mean:.1} is outside the plausible Phred+33 range \
         [25, 42] — check the encoding offset (Phred+64 would land near 65)"
    );
    // Phred+33 tops out at 41-42 for real instruments; a systematic excess above 45
    // means the offset is wrong even if the mean happens to land in range.
    let pct_above = 100.0 * above_45 as f64 / n as f64;
    assert!(
        pct_above <= 1.0,
        "{pct_above:.1}% of bases decode above Phred 45 — implausible for Phred+33"
    );
}
