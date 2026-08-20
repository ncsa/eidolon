//! FASTQ structural validation, matched against observed `samtools import` behaviour.
//!
//! Verdict parity is the correctness criterion: for every corpus file, this must reach
//! the same accept/reject decision samtools does. Where we are STRICTER, the finding is
//! a warning that says so — being stricter than the tool is defensible, but silently
//! calling a file invalid that samtools happily reads is not.
//!
//! Observed (samtools 1.18, 1.19.2, 1.22.1 — identical verdicts):
//!   reject: quality/sequence length mismatch, missing `+`, truncated record,
//!           header not starting `@`, non-IUPAC base (X, Z)
//!   accept: CRLF line endings, quality bytes outside the printable range,
//!           full IUPAC ambiguity codes, lowercase bases

use super::finding::{
    Finding, NOTHING_REJECTS_TOLERATED, SAMTOOLS_IMPORT_TRUNCATED, SAMTOOLS_TOLERATES,
};

/// IUPAC nucleotide codes. `X` and `Z` are amino-acid codes and are NOT valid here —
/// which is exactly where samtools draws the line, verified rather than assumed.
fn is_iupac_nucleotide(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
    )
}

pub fn validate_fastq(text: &str) -> Vec<Finding> {
    let mut findings = Vec::new();
    let had_crlf = text.contains("\r\n");
    let lines: Vec<&str> = text.lines().map(|l| l.trim_end_matches('\r')).collect();

    // Trailing blank lines are not a record; drop them before the multiple-of-4 check
    // so a file ending in a newline is not reported as truncated.
    let mut end = lines.len();
    while end > 0 && lines[end - 1].is_empty() {
        end -= 1;
    }
    let lines = &lines[..end];

    if had_crlf {
        findings.push(Finding::tolerated(
            1,
            "file uses CRLF line endings",
            SAMTOOLS_TOLERATES,
        ));
    }

    if lines.is_empty() {
        findings.push(Finding::error(
            1,
            "file contains no records",
            SAMTOOLS_IMPORT_TRUNCATED,
        ));
        return findings;
    }

    if lines.len() % 4 != 0 {
        findings.push(Finding::error(
            lines.len(),
            format!(
                "file has {} lines, not a multiple of 4 — the final record is truncated",
                lines.len()
            ),
            SAMTOOLS_IMPORT_TRUNCATED,
        ));
    }

    for (rec, chunk) in lines.chunks(4).enumerate() {
        let base = rec * 4 + 1; // 1-based line number of this record's header
        if chunk.len() < 4 {
            break; // already reported by the multiple-of-4 check above
        }
        let (header, seq, plus, qual) = (chunk[0], chunk[1], chunk[2], chunk[3]);

        if !header.starts_with('@') {
            findings.push(Finding::error(
                base,
                format!("header does not start with '@': {:?}", truncate(header)),
                SAMTOOLS_IMPORT_TRUNCATED,
            ));
        }
        if !plus.starts_with('+') {
            findings.push(Finding::error(
                base + 2,
                format!(
                    "separator line does not start with '+': {:?}",
                    truncate(plus)
                ),
                SAMTOOLS_IMPORT_TRUNCATED,
            ));
        }
        if seq.len() != qual.len() {
            findings.push(Finding::error(
                base + 3,
                format!(
                    "quality string is {} character(s) but the sequence is {}",
                    qual.len(),
                    seq.len()
                ),
                SAMTOOLS_IMPORT_TRUNCATED,
            ));
        }
        if let Some(pos) = seq.bytes().position(|b| !is_iupac_nucleotide(b)) {
            findings.push(Finding::error(
                base + 1,
                format!(
                    "sequence contains {:?} at offset {pos}, which is not an IUPAC \
                     nucleotide code",
                    seq.as_bytes()[pos] as char
                ),
                SAMTOOLS_IMPORT_TRUNCATED,
            ));
        }
        // Stricter than samtools on purpose, and labelled as such: a quality byte
        // outside the printable range decodes to a nonsensical Phred score, but
        // samtools reads the file without complaint on every version tested.
        if let Some(pos) = qual.bytes().position(|b| !(33..=126).contains(&b)) {
            findings.push(Finding::warning(
                base + 3,
                format!(
                    "quality byte 0x{:02x} at offset {pos} is outside the printable \
                     Phred+33 range (33..=126)",
                    qual.as_bytes()[pos]
                ),
                NOTHING_REJECTS_TOLERATED,
            ));
        }
    }
    findings
}

fn truncate(s: &str) -> String {
    if s.len() <= 40 {
        s.to_string()
    } else {
        format!("{}...", &s[..40])
    }
}
