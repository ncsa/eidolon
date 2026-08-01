//! A single validation finding, and the downstream failure it predicts.
//!
//! The design rule for this module: **every finding either names a tool and operation
//! that will reject the file, or says explicitly that none will.** A validator that
//! reports "invalid" without saying what breaks is just an opinion, and the second case
//! is not a weaker finding than the first — it is usually a stronger one. `bcftools`
//! silently converts a type-mismatched INFO value to `.`, so nothing downstream ever
//! complains and the data is simply gone. That is precisely how `AF=AF=0.3000` survived
//! in this project's own output.
//!
//! Citations are recorded from observed behaviour, not from documentation — see
//! `eidolon/test_data/validation_corpus/observed_tool_behaviour.tsv`.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// A named tool+operation rejects this file. Fixing it is not optional.
    Error,
    /// A spec violation that every tool tested tolerates — or worse, silently
    /// swallows. No downstream failure will surface it.
    Warning,
}

impl Severity {
    pub fn label(self) -> &'static str {
        match self {
            Severity::Error => "ERROR",
            Severity::Warning => "WARNING",
        }
    }
}

/// The concrete downstream failure a finding predicts.
///
/// `operation` is deliberately the full invocation, not just the tool name: `bcftools
/// view` accepts unsorted records while `tabix -p vcf` rejects them, and `bcftools view`
/// accepts a non-ACGTN REF while `bcftools norm -f` rejects it. A citation naming only
/// the binary would be wrong about half the time.
#[derive(Debug, Clone, Copy)]
pub struct Citation {
    pub operation: &'static str,
    pub message: &'static str,
}

#[derive(Debug, Clone)]
pub struct Finding {
    pub severity: Severity,
    /// 1-based line number in the source file, where the check can localise it.
    pub line: Option<usize>,
    pub what: String,
    pub citation: Option<Citation>,
    /// Used when nothing rejects the file, to say so rather than leave it implied.
    pub note: Option<&'static str>,
}

impl Finding {
    pub fn error(line: usize, what: impl Into<String>, citation: Citation) -> Self {
        Self {
            severity: Severity::Error,
            line: Some(line),
            what: what.into(),
            citation: Some(citation),
            note: None,
        }
    }

    pub fn warning(line: usize, what: impl Into<String>, note: &'static str) -> Self {
        Self {
            severity: Severity::Warning,
            line: Some(line),
            what: what.into(),
            citation: None,
            note: Some(note),
        }
    }

    /// A warning that a specific operation tolerates, recorded so the reader knows the
    /// tolerance was measured rather than assumed.
    pub fn tolerated(line: usize, what: impl Into<String>, citation: Citation) -> Self {
        Self {
            severity: Severity::Warning,
            line: Some(line),
            what: what.into(),
            citation: Some(citation),
            note: None,
        }
    }
}

impl fmt::Display for Finding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line {
            Some(n) => write!(f, "{:<7} line {n}: {}", self.severity.label(), self.what)?,
            None => write!(f, "{:<7} {}", self.severity.label(), self.what)?,
        }
        if let Some(c) = self.citation {
            match self.severity {
                Severity::Error => write!(f, "\n          `{}` rejects this:", c.operation)?,
                Severity::Warning => write!(f, "\n          `{}` tolerates this:", c.operation)?,
            }
            write!(f, "\n            {}", c.message)?;
        }
        if let Some(n) = self.note {
            write!(f, "\n          {n}")?;
        }
        Ok(())
    }
}

// ── Citations, all observed rather than quoted from documentation ────────────────
// samtools 1.18, 1.19.2 and 1.22.1 produced identical verdicts on the whole corpus;
// bcftools 1.17, 1.19 and 1.22 likewise. 1.22 is what Delta runs. Messages can drift
// between versions, which is why the differential test asserts verdict agreement and
// treats these strings as documentation to re-check, not as contract.

/// samtools reports FIVE distinct FASTQ defects with this one message. Saying which
/// line and which defect is the entire reason this validator exists.
pub const SAMTOOLS_IMPORT_TRUNCATED: Citation = Citation {
    operation: "samtools import -0",
    message: "samtools import: truncated file. Aborting",
};

pub const SAMTOOLS_TOLERATES: Citation = Citation {
    operation: "samtools import -0",
    message: "accepted without complaint (verified on 1.18, 1.19.2, 1.22.1)",
};

pub const BCFTOOLS_BROKEN_COLUMNS: Citation = Citation {
    operation: "bcftools view",
    message: "[E::bcf_write] Broken VCF record, the number of columns does not match",
};

pub const BCFTOOLS_CONCAT_UNDECLARED_TAG: Citation = Citation {
    operation: "bcftools concat",
    message: "[E::bcf_translate] Unchecked error (2 Tag not defined in header), exiting",
};

pub const BCFTOOLS_FILTER_UNDECLARED_TAG: Citation = Citation {
    operation: "bcftools view -i 'INFO/<TAG>...'",
    message: "[filter.c] Error: the tag \"<TAG>\" is not defined in the VCF header",
};

pub const TABIX_UNSORTED: Citation = Citation {
    operation: "tabix -p vcf",
    message: "[E::hts_idx_push] Unsorted positions on sequence #N",
};

pub const BCFTOOLS_NORM_BAD_REF: Citation = Citation {
    operation: "bcftools norm -f <reference>",
    message: "Non-ACGTN reference allele at <pos> .. REF_SEQ:'<x>' vs VCF:'<y>'",
};

/// The most dangerous category: nothing fails, and the data quietly disappears.
pub const NOTHING_REJECTS_SILENT_LOSS: &str = "No tool rejects this. bcftools silently converts the value to `.` on every path \
     tested (view, view -O b, query), so it is lost WITHOUT a warning.";

pub const SAMTOOLS_QUICKCHECK_EOF: Citation = Citation {
    operation: "samtools quickcheck",
    message: "<file> was missing EOF block when one should be present.",
};

pub const SAMTOOLS_VIEW_UNREADABLE: Citation = Citation {
    operation: "samtools view",
    message: "[E::hts_hopen] Failed to open file / truncated stream",
};

/// Observed on the SAM path; samtools refuses to WRITE such a BAM, so the BAM-path
/// wording could not be captured. eidolon writes BAM directly via noodles and is not
/// bound by that refusal, which is why the check exists at all.
pub const SAMTOOLS_VIEW_QUAL_LEN: Citation = Citation {
    operation: "samtools view",
    message: "[E::sam_parse1] SEQ and QUAL are of different length",
};

/// The SAM and BAM readers word this differently; the BAM one is quoted because this
/// citation is only used for BAM. Verified on 1.22.1 against a hand-patched record.
pub const SAMTOOLS_VIEW_CIGAR_LEN: Citation = Citation {
    operation: "samtools view",
    message: "[E::bam_read1] CIGAR and query sequence lengths differ for <read>",
};

pub const SAMTOOLS_INDEX_UNSORTED: Citation = Citation {
    operation: "samtools index",
    message: "[E::hts_idx_push] Unsorted positions on sequence #N",
};

/// A read aligned past the end of its own contig. quickcheck, view AND index all
/// accept it — verified on 1.22.1 — so nothing downstream will ever surface it.
pub const NOTHING_REJECTS_BAM_POS: &str = "No tool rejects this: samtools quickcheck, view and index all accept it. The \
     record aligns past the end of the reference, so any consumer computing a span \
     from it reads off the end of the contig.";

pub const NOTHING_REJECTS_TOLERATED: &str =
    "No tool tested rejects this; it is a spec deviation only.";
