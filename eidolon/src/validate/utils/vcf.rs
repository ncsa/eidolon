//! VCF structural validation, matched against observed bcftools/tabix behaviour.
//!
//! Observed (bcftools 1.17, 1.19, 1.22 — identical verdicts; 1.22 is the Delta target):
//!   `bcftools view`      rejects a wrong column count; only WARNS on an undeclared
//!                        INFO tag or contig; accepts unsorted records, a non-ACGTN
//!                        REF, and a type-mismatched INFO value
//!   `bcftools concat`    HARD FAILS on an undeclared INFO tag (bcf_translate)
//!   `tabix -p vcf`       rejects unsorted records
//!   `bcftools norm -f`   rejects a non-ACGTN REF
//!
//! Severity follows the strictest operation this project actually runs, because a
//! finding that only matters to a command nobody invokes is noise. `bcftools concat |
//! bcftools sort` is exactly what broke in #444, so an undeclared INFO tag is an ERROR
//! here even though `bcftools view` merely warns.

use super::finding::{
    BCFTOOLS_BROKEN_COLUMNS, BCFTOOLS_CONCAT_UNDECLARED_TAG, BCFTOOLS_NORM_BAD_REF, Finding,
    NOTHING_REJECTS_SILENT_LOSS, NOTHING_REJECTS_TOLERATED, TABIX_UNSORTED,
};
use std::collections::{HashMap, HashSet};

/// Minimum VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO.
const MIN_COLUMNS: usize = 8;

fn declared_id(line: &str, prefix: &str) -> Option<String> {
    let rest = line.strip_prefix(prefix)?;
    let rest = rest.strip_prefix("<")?;
    rest.split(',')
        .find_map(|f| f.trim().strip_prefix("ID="))
        .map(|v| v.trim_matches('"').to_string())
}

/// `Type=` of an `##INFO` declaration, for the value-conformance check.
fn declared_type(line: &str) -> Option<String> {
    line.split(',')
        .find_map(|f| f.trim().strip_prefix("Type="))
        .map(|v| v.trim_matches(|c| c == '"' || c == '>').to_string())
}

pub fn validate_vcf(text: &str) -> Vec<Finding> {
    let mut findings = Vec::new();
    let mut info_types: HashMap<String, String> = HashMap::new();
    let mut contigs: HashSet<String> = HashSet::new();
    let mut n_header_cols = 0usize;
    let mut seen_chrom_line = false;

    for line in text.lines() {
        if let Some(id) = declared_id(line, "##INFO=") {
            let ty = declared_type(line).unwrap_or_else(|| "String".to_string());
            info_types.insert(id, ty);
        } else if let Some(id) = declared_id(line, "##contig=") {
            contigs.insert(id);
        } else if line.starts_with("#CHROM") {
            seen_chrom_line = true;
            n_header_cols = line.split('\t').count();
        }
    }
    if !seen_chrom_line {
        findings.push(Finding::error(
            1,
            "no #CHROM header line",
            BCFTOOLS_BROKEN_COLUMNS,
        ));
        return findings;
    }

    // (contig, pos) of the previous record, for the sortedness check.
    let mut prev: Option<(String, u64)> = None;
    let mut contig_order: Vec<String> = Vec::new();

    for (idx, line) in text.lines().enumerate() {
        let lineno = idx + 1;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < MIN_COLUMNS || (n_header_cols > MIN_COLUMNS && cols.len() != n_header_cols)
        {
            findings.push(Finding::error(
                lineno,
                format!(
                    "record has {} column(s); the #CHROM header declares {n_header_cols}",
                    cols.len()
                ),
                BCFTOOLS_BROKEN_COLUMNS,
            ));
            continue;
        }
        let (chrom, pos_s, refb, info) = (cols[0], cols[1], cols[3], cols[7]);
        let Ok(pos) = pos_s.parse::<u64>() else {
            findings.push(Finding::error(
                lineno,
                format!("POS {pos_s:?} is not an integer"),
                BCFTOOLS_BROKEN_COLUMNS,
            ));
            continue;
        };

        if !contigs.is_empty() && !contigs.contains(chrom) {
            findings.push(Finding::warning(
                lineno,
                format!("contig {chrom:?} is not declared in the header"),
                NOTHING_REJECTS_TOLERATED,
            ));
        }

        // REF must be ACGTN. bcftools view accepts anything; norm -f is what catches it.
        if let Some(p) = refb
            .bytes()
            .position(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
        {
            findings.push(Finding::error(
                lineno,
                format!(
                    "REF contains {:?}, which is not an ACGTN base",
                    refb.as_bytes()[p] as char
                ),
                BCFTOOLS_NORM_BAD_REF,
            ));
        }

        // Sortedness: positions must ascend within a contig, and a contig must not
        // reappear after another one has started.
        if let Some((pchrom, ppos)) = &prev {
            if pchrom == chrom && pos < *ppos {
                findings.push(Finding::error(
                    lineno,
                    format!("POS {pos} follows {ppos} on {chrom} — records are unsorted"),
                    TABIX_UNSORTED,
                ));
            } else if pchrom != chrom && contig_order.iter().any(|c| c == chrom) {
                findings.push(Finding::error(
                    lineno,
                    format!(
                        "contig {chrom:?} reappears after another contig — records are unsorted"
                    ),
                    TABIX_UNSORTED,
                ));
            }
        }
        if !contig_order.iter().any(|c| c == chrom) {
            contig_order.push(chrom.to_string());
        }
        prev = Some((chrom.to_string(), pos));

        // INFO: undeclared tags, and values that do not match their declared type.
        if info != "." {
            for field in info.split(';') {
                if field.is_empty() {
                    continue;
                }
                let (key, val) = match field.split_once('=') {
                    Some((k, v)) => (k, Some(v)),
                    None => (field, None),
                };
                match info_types.get(key) {
                    None => findings.push(Finding::error(
                        lineno,
                        format!("INFO tag {key:?} is not declared in the header"),
                        BCFTOOLS_CONCAT_UNDECLARED_TAG,
                    )),
                    Some(ty) => {
                        if let Some(v) = val
                            && !value_matches_type(v, ty)
                        {
                            findings.push(Finding::warning(
                                lineno,
                                format!("INFO/{key} is declared Type={ty} but the value is {v:?}"),
                                NOTHING_REJECTS_SILENT_LOSS,
                            ));
                        }
                    }
                }
            }
        }
    }
    findings
}

/// Does `value` conform to its declared INFO type? Comma-separated lists are checked
/// element-wise; `.` is the missing value and always conforms.
fn value_matches_type(value: &str, ty: &str) -> bool {
    value.split(',').all(|v| {
        if v == "." {
            return true;
        }
        match ty {
            "Integer" => v.parse::<i64>().is_ok(),
            "Float" => v.parse::<f64>().is_ok(),
            // Flags carry no value; String/Character accept anything printable.
            _ => true,
        }
    })
}
