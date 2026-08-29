//! BAM structural validation, matched against observed samtools behaviour.
//!
//! Observed on samtools 1.22.1 (the version Delta runs), by exit status:
//!
//! | defect              | quickcheck | view   | index  |
//! |---------------------|------------|--------|--------|
//! | well-formed         | accept     | accept | accept |
//! | missing EOF block   | REJECT     | accept | accept |
//! | truncated mid-file  | REJECT     | REJECT | REJECT |
//! | bad magic bytes     | REJECT     | REJECT | REJECT |
//! | unsorted            | accept     | accept | REJECT |
//! | SEQ/QUAL mismatch   | —          | REJECT | —      |
//! | CIGAR/SEQ mismatch  | —          | REJECT | —      |
//! | POS beyond contig   | accept     | accept | accept |  <-- NOTHING catches this
//!
//! The last row is why this is worth writing. A record positioned past the end of its
//! own contig is nonsense — it aligns a read off the end of the reference — and every
//! samtools operation tested waves it through. Same category as a type-mismatched INFO
//! value: no downstream failure will ever surface it.

use super::finding::{
    Finding, NOTHING_REJECTS_BAM_POS, NOTHING_REJECTS_TOLERATED, SAMTOOLS_INDEX_UNSORTED,
    SAMTOOLS_QUICKCHECK_EOF, SAMTOOLS_VIEW_CIGAR_LEN, SAMTOOLS_VIEW_QUAL_LEN,
    SAMTOOLS_VIEW_UNREADABLE,
};
use noodles::bam;
use std::path::Path;

/// The canonical 28-byte BGZF EOF marker every complete BGZF file ends with.
/// `samtools quickcheck` exists largely to check for this — a truncated transfer is the
/// single most common way a BAM goes bad, and it is invisible to `samtools view`.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

pub fn validate_bam(path: &Path) -> Vec<Finding> {
    let mut findings = Vec::new();

    let bytes = match std::fs::read(path) {
        Ok(b) => b,
        Err(e) => {
            findings.push(Finding::error(
                1,
                format!("cannot read the file: {e}"),
                SAMTOOLS_VIEW_UNREADABLE,
            ));
            return findings;
        }
    };
    if !bytes.starts_with(&[0x1f, 0x8b]) {
        findings.push(Finding::error(
            1,
            "file does not begin with the gzip/BGZF magic bytes 0x1f 0x8b",
            SAMTOOLS_VIEW_UNREADABLE,
        ));
        return findings;
    }
    // Checked before parsing: a file can read back perfectly and still be truncated,
    // which is exactly the case `samtools view` accepts and `quickcheck` rejects.
    if bytes.len() < BGZF_EOF.len() || !bytes.ends_with(&BGZF_EOF) {
        findings.push(Finding::error(
            1,
            "missing the BGZF EOF marker — the file is truncated or was not closed",
            SAMTOOLS_QUICKCHECK_EOF,
        ));
    }

    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(e) => {
            findings.push(Finding::error(
                1,
                format!("cannot open the file: {e}"),
                SAMTOOLS_VIEW_UNREADABLE,
            ));
            return findings;
        }
    };
    let mut reader = bam::io::Reader::new(file);
    let header = match reader.read_header() {
        Ok(h) => h,
        Err(e) => {
            findings.push(Finding::error(
                1,
                format!("header could not be read: {e}"),
                SAMTOOLS_VIEW_UNREADABLE,
            ));
            return findings;
        }
    };

    // Contig lengths, indexed the same way records reference them.
    let contigs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| {
            (
                String::from_utf8_lossy(name).into_owned(),
                seq.length().get(),
            )
        })
        .collect();
    if contigs.is_empty() {
        findings.push(Finding::warning(
            1,
            "header declares no reference sequences (@SQ)",
            NOTHING_REJECTS_TOLERATED,
        ));
    }

    // SO lives in @HD's other-fields map in this noodles version; there is no typed
    // accessor. Checked case-insensitively because the spec's value is lowercase
    // "coordinate" but writers vary.
    let coordinate_sorted = header
        .header()
        .map(|hd| {
            hd.other_fields().iter().any(|(tag, value)| {
                tag.as_ref() == b"SO"
                    && String::from_utf8_lossy(value).eq_ignore_ascii_case("coordinate")
            })
        })
        .unwrap_or(false);

    let mut prev: Option<(usize, usize)> = None; // (ref index, start)
    let mut n = 0usize;

    for (idx, result) in reader.records().enumerate() {
        let recno = idx + 1;
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                findings.push(Finding::error(
                    recno,
                    format!("record {recno} could not be parsed: {e}"),
                    SAMTOOLS_VIEW_UNREADABLE,
                ));
                break; // the stream is unusable from here
            }
        };
        n += 1;

        let seq_len = record.sequence().len();
        let qual_len = record.quality_scores().as_ref().len();
        // QUAL may legitimately be absent ("*"), which BAM stores as length 0.
        if qual_len != 0 && qual_len != seq_len {
            findings.push(Finding::error(
                recno,
                format!("record {recno}: SEQ is {seq_len} base(s) but QUAL is {qual_len}"),
                SAMTOOLS_VIEW_QUAL_LEN,
            ));
        }

        // CIGAR must consume exactly as many query bases as SEQ has.
        let mut query_consumed = 0usize;
        let mut ref_span = 0usize;
        let mut cigar_ok = true;
        for op in record.cigar().iter() {
            match op {
                Ok(op) => {
                    if op.kind().consumes_read() {
                        query_consumed += op.len();
                    }
                    if op.kind().consumes_reference() {
                        ref_span += op.len();
                    }
                }
                Err(_) => {
                    cigar_ok = false;
                    break;
                }
            }
        }
        if cigar_ok && query_consumed != 0 && seq_len != 0 && query_consumed != seq_len {
            findings.push(Finding::error(
                recno,
                format!(
                    "record {recno}: CIGAR consumes {query_consumed} query base(s) but SEQ is {seq_len}"
                ),
                SAMTOOLS_VIEW_CIGAR_LEN,
            ));
        }

        // Placement checks only apply to mapped records.
        let ref_id = record.reference_sequence_id().and_then(|r| r.ok());
        let start = record
            .alignment_start()
            .and_then(|r| r.ok())
            .map(|p| p.get());
        if let (Some(ref_id), Some(start)) = (ref_id, start) {
            if let Some((name, len)) = contigs.get(ref_id) {
                let end = start + ref_span.saturating_sub(1).max(0);
                if start > *len || end > *len {
                    findings.push(Finding::warning(
                        recno,
                        format!(
                            "record {recno}: aligned at {name}:{start} spanning {ref_span} base(s), \
                             which ends past the contig's declared length of {len}"
                        ),
                        NOTHING_REJECTS_BAM_POS,
                    ));
                }
            } else {
                findings.push(Finding::error(
                    recno,
                    format!(
                        "record {recno}: reference id {ref_id} has no matching @SQ entry \
                         (header declares {})",
                        contigs.len()
                    ),
                    SAMTOOLS_VIEW_UNREADABLE,
                ));
            }
            if coordinate_sorted
                && let Some((pref, pstart)) = prev
                && (ref_id < pref || (ref_id == pref && start < pstart))
            {
                findings.push(Finding::error(
                    recno,
                    format!(
                        "record {recno}: position {start} follows {pstart} on reference \
                             {ref_id} — the header declares SO:coordinate but records are unsorted"
                    ),
                    SAMTOOLS_INDEX_UNSORTED,
                ));
            }
            prev = Some((ref_id, start));
        }
    }

    if n == 0 {
        findings.push(Finding::warning(
            1,
            "file contains no alignment records",
            NOTHING_REJECTS_TOLERATED,
        ));
    }
    findings
}
