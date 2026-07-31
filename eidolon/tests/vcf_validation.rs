//! Well-formedness validation for the golden VCF `gen-reads` emits.
//!
//! The VCF counterpart to `fastq_validation.rs`. It exists because a whole class of
//! defect shipped with a green suite: every prior VCF test asserted that *specific
//! expected content was present* (`lines.iter().any(...)`), and none asserted the
//! file as a whole was **well formed**. That let through
//!
//!   * `MATEID` emitted with no `##INFO=<ID=MATEID,…>` declaration,
//!   * the record ID column hardcoded to `.`, so `MATEID` pointed at nothing,
//!   * de-novo BNDs emitted as a lone breakend with no mate record (#451),
//!   * a BND ALT embedding a literal `N` that contradicted its own REF column.
//!
//! An undeclared INFO tag is the sharp edge: bcftools tolerates it streaming
//! VCF->VCF but hard-fails on BCF translation, and `bcftools concat | sort` does
//! that internally — so it surfaces far downstream inside a harness rather than here.
//!
//! Dependency-free (no bcftools) so it runs in CI. Checks collect every violation
//! instead of panicking on the first.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

struct Vcf {
    info_declared: HashSet<String>,
    alt_declared: HashSet<String>,
    contig_declared: HashSet<String>,
    column_header: Vec<String>,
    records: Vec<Record>,
}

struct Record {
    line_no: usize,
    chrom: String,
    pos: usize,
    id: String,
    reference: String,
    alt: String,
    info: String,
    format: String,
    samples: Vec<String>,
    raw: String,
}

impl Record {
    fn info_has(&self, field: &str) -> bool {
        self.info.split(';').any(|f| f == field)
    }
    fn mateid(&self) -> Option<&str> {
        self.info.split(';').find_map(|f| f.strip_prefix("MATEID="))
    }
    fn is_bnd(&self) -> bool {
        self.info_has("SVTYPE=BND") || self.alt.contains('[') || self.alt.contains(']')
    }
}

fn read_lines(path: &Path) -> Vec<String> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        BufReader::new(MultiGzDecoder::new(f))
            .lines()
            .map(|l| l.unwrap())
            .collect()
    } else {
        BufReader::new(f).lines().map(|l| l.unwrap()).collect()
    }
}

/// Pull `ID=` out of a `##INFO=<ID=FOO,...>` style header line.
fn header_id(line: &str) -> Option<String> {
    let open = line.find('<')?;
    line[open + 1..]
        .trim_end_matches('>')
        .split(',')
        .find_map(|f| f.strip_prefix("ID=").map(|v| v.to_string()))
}

fn parse_vcf(path: &Path) -> Vcf {
    let mut v = Vcf {
        info_declared: HashSet::new(),
        alt_declared: HashSet::new(),
        contig_declared: HashSet::new(),
        column_header: Vec::new(),
        records: Vec::new(),
    };
    for (i, line) in read_lines(path).into_iter().enumerate() {
        let line_no = i + 1;
        if line.starts_with("##INFO=") {
            v.info_declared.extend(header_id(&line));
        } else if line.starts_with("##ALT=") {
            v.alt_declared.extend(header_id(&line));
        } else if line.starts_with("##contig=") {
            v.contig_declared.extend(header_id(&line));
        } else if line.starts_with("#CHROM") {
            v.column_header = line.split('\t').map(str::to_string).collect();
        } else if !line.starts_with('#') && !line.trim().is_empty() {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 8, "line {line_no}: <8 columns: {line:?}");
            v.records.push(Record {
                line_no,
                chrom: f[0].to_string(),
                pos: f[1]
                    .parse()
                    .unwrap_or_else(|_| panic!("line {line_no}: POS not numeric: {:?}", f[1])),
                id: f[2].to_string(),
                reference: f[3].to_string(),
                alt: f[4].to_string(),
                info: f[7].to_string(),
                format: f.get(8).map(|s| s.to_string()).unwrap_or_default(),
                samples: f.iter().skip(9).map(|s| s.to_string()).collect(),
                raw: line.clone(),
            });
        }
    }
    v
}

/// Every INFO key a record uses must be declared. This is the BCF-translation class.
fn check_info_declared(v: &Vcf) -> Vec<String> {
    let mut p = Vec::new();
    for r in &v.records {
        if r.info == "." || r.info.is_empty() {
            continue;
        }
        for field in r.info.split(';').filter(|f| !f.is_empty()) {
            let key = field.split('=').next().unwrap_or(field);
            if !v.info_declared.contains(key) {
                p.push(format!(
                    "line {}: INFO key `{key}` used but never declared (bcftools rejects \
                     this on BCF translation)",
                    r.line_no
                ));
            }
        }
    }
    p
}

/// Symbolic ALTs (`<DEL>`, `<CNV>`, …) need a matching `##ALT=<ID=…>`.
fn check_alt_declared(v: &Vcf) -> Vec<String> {
    let mut p = Vec::new();
    for r in &v.records {
        if r.alt.starts_with('<') && r.alt.ends_with('>') {
            let id = &r.alt[1..r.alt.len() - 1];
            if !v.alt_declared.contains(id) {
                p.push(format!(
                    "line {}: symbolic ALT `{}` has no ##ALT=<ID={id}> declaration",
                    r.line_no, r.alt
                ));
            }
        }
    }
    p
}

/// REF must equal the reference FASTA at POS, and a BND's ALT embeds that same base.
fn check_ref(v: &Vcf, reference: &HashMap<String, String>) -> Vec<String> {
    let mut p = Vec::new();
    for r in &v.records {
        let Some(seq) = reference.get(&r.chrom) else {
            p.push(format!(
                "line {}: CHROM `{}` absent from the reference",
                r.line_no, r.chrom
            ));
            continue;
        };
        let (start, end) = (r.pos - 1, r.pos - 1 + r.reference.len());
        if end > seq.len() {
            p.push(format!(
                "line {}: REF runs past end of {}",
                r.line_no, r.chrom
            ));
            continue;
        }
        if !seq[start..end].eq_ignore_ascii_case(&r.reference) {
            p.push(format!(
                "line {}: REF `{}` != reference `{}` at {}:{}",
                r.line_no,
                r.reference,
                &seq[start..end],
                r.chrom,
                r.pos
            ));
        }
        if r.is_bnd() {
            let embedded: String = r
                .alt
                .chars()
                .take_while(char::is_ascii_alphabetic)
                .collect();
            if !embedded.is_empty() && !embedded.eq_ignore_ascii_case(&r.reference) {
                p.push(format!(
                    "line {}: BND ALT embeds `{embedded}` but REF is `{}` ({})",
                    r.line_no, r.reference, r.alt
                ));
            }
        }
    }
    p
}

/// A BND must have a MATEID that resolves to a record pointing back at it. A lone
/// breakend is unmatchable by any comparison tool.
fn check_mateid(v: &Vcf) -> Vec<String> {
    let mut p = Vec::new();
    let by_id: HashMap<&str, &Record> = v
        .records
        .iter()
        .filter(|r| r.id != ".")
        .map(|r| (r.id.as_str(), r))
        .collect();
    for r in &v.records {
        if r.is_bnd() && r.mateid().is_none() {
            p.push(format!(
                "line {}: BND with no MATEID — a lone breakend cannot be matched",
                r.line_no
            ));
            continue;
        }
        if let Some(mid) = r.mateid() {
            match by_id.get(mid) {
                None => p.push(format!(
                    "line {}: MATEID=`{mid}` resolves to no record ID in this file",
                    r.line_no
                )),
                Some(mate) if mate.mateid() != Some(r.id.as_str()) => p.push(format!(
                    "line {}: MATEID linkage not reciprocal (-> `{mid}` -> `{:?}`)",
                    r.line_no,
                    mate.mateid()
                )),
                _ => {}
            }
        }
    }
    p
}

/// Position-sorted within each contig — tabix requires it, and unsorted output
/// broke the cancer merge once already (#185).
fn check_sorted(v: &Vcf) -> Vec<String> {
    let mut p = Vec::new();
    let mut last: HashMap<&str, usize> = HashMap::new();
    for r in &v.records {
        if let Some(&prev) = last.get(r.chrom.as_str()) {
            if r.pos < prev {
                p.push(format!(
                    "line {}: POS {} < previous {} on {} (tabix needs sorted input)",
                    r.line_no, r.pos, prev, r.chrom
                ));
            }
        }
        last.insert(r.chrom.as_str(), r.pos);
    }
    p
}

/// Column count matches #CHROM, and FORMAT key count matches each SAMPLE's values.
fn check_columns(v: &Vcf) -> Vec<String> {
    let mut p = Vec::new();
    let expected = v.column_header.len();
    for r in &v.records {
        let got = 9 + r.samples.len();
        if expected > 0 && got != expected {
            p.push(format!(
                "line {}: {got} columns but #CHROM declares {expected}",
                r.line_no
            ));
        }
        if r.format.is_empty() || r.format == "." {
            continue;
        }
        let nkeys = r.format.split(':').count();
        for (i, s) in r.samples.iter().enumerate() {
            if s.split(':').count() != nkeys {
                p.push(format!(
                    "line {}: FORMAT has {nkeys} keys but sample {i} has {} values ({} vs {s})",
                    r.line_no,
                    s.split(':').count(),
                    r.format
                ));
            }
        }
    }
    p
}

fn check_contigs(v: &Vcf) -> Vec<String> {
    if v.contig_declared.is_empty() {
        return vec![
            "header declares no ##contig lines — BCF translation and any tool building a \
             sequence dictionary from the VCF need them"
                .to_string(),
        ];
    }
    v.records
        .iter()
        .filter(|r| !v.contig_declared.contains(&r.chrom))
        .map(|r| {
            format!(
                "line {}: CHROM `{}` has no ##contig declaration",
                r.line_no, r.chrom
            )
        })
        .collect()
}

fn load_reference(path: &Path) -> HashMap<String, String> {
    let mut out = HashMap::new();
    let (mut name, mut seq) = (String::new(), String::new());
    for line in read_lines(path) {
        if let Some(h) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.insert(name.clone(), std::mem::take(&mut seq));
            }
            name = h.split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() {
        out.insert(name, seq);
    }
    out
}

fn report(label: &str, problems: Vec<String>) -> Vec<String> {
    if !problems.is_empty() {
        eprintln!("--- {label}: {} problem(s) ---", problems.len());
        for x in problems.iter().take(20) {
            eprintln!("  {x}");
        }
    }
    problems
}

fn golden_vcf(test_name: &str, sv_rate_scale: Option<f64>) -> (Vcf, HashMap<String, String>) {
    let (_dir, work) = fresh_workdir();
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), test_name);
    config.coverage = 20;
    config.read_len = 101;
    config.produce_fastq = false;
    config.produce_vcf = true;
    config.sv_rate_scale = sv_rate_scale;
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let path = work.join(format!("{test_name}.vcf.gz"));
    assert!(path.exists(), "golden VCF not produced at {path:?}");
    let parsed = parse_vcf(&path);
    assert!(
        !parsed.records.is_empty(),
        "golden VCF has no records — validating an empty file proves nothing"
    );
    // Keep the tempdir alive by leaking it: the Vcf is fully parsed into memory, so
    // the directory can go, but the reference must still be readable.
    let reference = load_reference(&h1n1_reference());
    std::mem::forget(_dir);
    (parsed, reference)
}

fn validate_all(v: &Vcf, reference: &HashMap<String, String>) -> Vec<String> {
    let mut p = Vec::new();
    p.extend(report("INFO declared", check_info_declared(v)));
    p.extend(report("ALT declared", check_alt_declared(v)));
    p.extend(report("REF matches reference", check_ref(v, reference)));
    p.extend(report("MATEID resolves", check_mateid(v)));
    p.extend(report("sorted", check_sorted(v)));
    p.extend(report("columns", check_columns(v)));
    p
}

#[test]
fn golden_vcf_without_svs_is_well_formed() {
    let (v, r) = golden_vcf("vcfval_nosv", None);
    let p = validate_all(&v, &r);
    assert!(
        p.is_empty(),
        "golden VCF not well formed ({} problem(s)); first: {}",
        p.len(),
        p[0]
    );
}

/// With SVs the file carries symbolic ALTs, SV INFO fields and BND mate pairs — the
/// surface where every previous defect lived.
#[test]
fn golden_vcf_with_svs_is_well_formed() {
    let (v, r) = golden_vcf("vcfval_sv", Some(8.0));
    let svs = v
        .records
        .iter()
        .filter(|x| x.alt.starts_with('<') || x.is_bnd())
        .count();
    assert!(svs > 0, "no SV records emitted — test would vacuously pass");
    let p = validate_all(&v, &r);
    assert!(
        p.is_empty(),
        "golden VCF with SVs not well formed ({} problem(s)); first: {}",
        p.len(),
        p[0]
    );
}

/// Every BND must be a reciprocal mate pair whose two ALTs name each other's
/// position. Before #451 the de-novo path emitted a single ID-less breakend.
#[test]
fn every_bnd_is_a_reciprocal_mate_pair() {
    let (v, _r) = golden_vcf("vcfval_bnd", Some(8.0));
    let bnds: Vec<&Record> = v.records.iter().filter(|r| r.is_bnd()).collect();
    assert!(
        !bnds.is_empty(),
        "no BND records emitted — test would vacuously pass"
    );
    assert_eq!(
        bnds.len() % 2,
        0,
        "odd BND count ({}) — a breakend has no mate",
        bnds.len()
    );
    let by_id: HashMap<&str, &Record> = bnds.iter().map(|r| (r.id.as_str(), *r)).collect();
    for r in &bnds {
        assert_ne!(r.id, ".", "BND record has no ID: {}", r.raw);
        let mid = r
            .mateid()
            .unwrap_or_else(|| panic!("BND without MATEID: {}", r.raw));
        let mate = by_id
            .get(mid)
            .unwrap_or_else(|| panic!("MATEID={mid} unresolved: {}", r.raw));
        assert!(
            r.alt.contains(&format!("{}:{}", mate.chrom, mate.pos)),
            "BND ALT `{}` does not name its mate at {}:{}",
            r.alt,
            mate.chrom,
            mate.pos
        );
    }
}

/// Missing `##contig` lines are a real interoperability gap — assert on it rather
/// than leaving it an unstated assumption.
#[test]
fn golden_vcf_declares_contigs() {
    let (v, _r) = golden_vcf("vcfval_contig", None);
    let p = report("contigs", check_contigs(&v));
    assert!(p.is_empty(), "contig declarations missing: {}", p[0]);
}
