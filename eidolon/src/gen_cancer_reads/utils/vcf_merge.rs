//! Merge the two per-pass golden VCFs into one origin-tagged truth VCF (replaces
//! the `bcftools isec` + `awk` annotate step in `tools/cancer_simulate.sh`).
//!
//! Both inputs are eidolon-golden VCFs with identical record representation (the
//! tumor pass carries normal's germline verbatim via `input_vcf`, tagged
//! `INFO/EIDOLON_PROVENANCE=input`; somatic de-novo records are `=denovo`, and
//! reproductive somatic-VCF records (#405) are `=somatic_input`). So an exact
//! `(contig,pos,ref,alt)` key is sufficient — no positional/normalized matching
//! needed. Origin rules:
//!   - tumor `denovo` / `somatic_input` -> `somatic`  (tumor-only)
//!   - tumor `input` (or other)         -> `shared`   (germline carried into tumor)
//!   - normal key absent in tumor       -> `germline` (germline-only; lost in tumor)

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

use crate::gen_cancer_reads::errors::GenCancerReadsError;

const EIDOLON_ORIGIN_DECL: &str = "##INFO=<ID=EIDOLON_ORIGIN,Number=1,Type=String,\
Description=\"Origin in tumor/normal mix: germline | somatic | shared\">";

/// SV INFO declarations injected if the per-pass headers don't already carry
/// them, so symbolic-SV records in the merged truth validate.
const SV_INFO_DECLS: &[(&str, &str)] = &[
    (
        "SVTYPE",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    ),
    (
        "END",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">",
    ),
    (
        "SVLEN",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT\">",
    ),
    (
        "CN",
        "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of segment\">",
    ),
    (
        "MATEID",
        "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">",
    ),
];

type Key = (String, String, String, String);

fn key_of(cols: &[&str]) -> Option<Key> {
    if cols.len() < 5 {
        return None;
    }
    Some((
        cols[0].into(),
        cols[1].into(),
        cols[3].into(),
        cols[4].into(),
    ))
}

/// Origin for a tumor-pass record, from its INFO/EIDOLON_PROVENANCE.
///
/// De-novo-sampled variants and reproductive somatic-VCF variants (#405) are both
/// tumor-only → `somatic`; anything else in the tumor pass is germline carried
/// through `input_vcf` → `shared`.
pub fn tumor_origin(info: &str) -> &'static str {
    let is_somatic = info
        .split(';')
        .any(|f| f == "EIDOLON_PROVENANCE=denovo" || f == "EIDOLON_PROVENANCE=somatic_input");
    if is_somatic { "somatic" } else { "shared" }
}

/// Append `EIDOLON_ORIGIN=<origin>` to an INFO column.
pub fn append_origin(info: &str, origin: &str) -> String {
    if info == "." || info.is_empty() {
        format!("EIDOLON_ORIGIN={origin}")
    } else {
        format!("{info};EIDOLON_ORIGIN={origin}")
    }
}

struct Rec {
    contig: String,
    pos: i64,
    line: String,
}

fn make_rec(
    cols: &[&str],
    origin: &str,
    order: &mut HashMap<String, usize>,
    next: &mut usize,
) -> Rec {
    let contig = cols[0].to_string();
    let pos = cols[1].parse::<i64>().unwrap_or(0);
    if !order.contains_key(&contig) {
        order.insert(contig.clone(), *next);
        *next += 1;
    }
    let mut new_cols: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
    if new_cols.len() > 7 {
        new_cols[7] = append_origin(&new_cols[7], origin);
    }
    Rec {
        contig,
        pos,
        line: new_cols.join("\t"),
    }
}

fn parse_info_id(line: &str) -> Option<String> {
    let rest = line.strip_prefix("##INFO=<ID=")?;
    Some(rest.split([',', '>']).next()?.to_string())
}

/// (meta `##` lines, `#CHROM` line, set of already-declared INFO IDs).
fn read_header(path: &Path) -> Result<(Vec<String>, String, HashSet<String>), GenCancerReadsError> {
    let r = BufReader::new(MultiGzDecoder::new(File::open(path)?));
    let mut meta = Vec::new();
    let mut chrom = String::new();
    let mut info_ids = HashSet::new();
    for line in r.lines() {
        let line = line?;
        if line.starts_with("##") {
            if let Some(id) = parse_info_id(&line) {
                info_ids.insert(id);
            }
            meta.push(line);
        } else if line.starts_with("#CHROM") {
            chrom = line;
            break;
        }
    }
    if chrom.is_empty() {
        chrom = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE".to_string();
    }
    Ok((meta, chrom, info_ids))
}

fn data_records(path: &Path) -> Result<Vec<String>, GenCancerReadsError> {
    let r = BufReader::new(MultiGzDecoder::new(File::open(path)?));
    let mut v = Vec::new();
    for line in r.lines() {
        let line = line?;
        if !line.starts_with('#') && !line.is_empty() {
            v.push(line);
        }
    }
    Ok(v)
}

pub fn merge_goldens(
    normal_vcf: &Path,
    tumor_vcf: &Path,
    out: &Path,
) -> Result<(), GenCancerReadsError> {
    let (meta, chrom_line, declared) = read_header(normal_vcf)?;

    let mut recs: Vec<Rec> = Vec::new();
    let mut tumor_keys: HashSet<Key> = HashSet::new();
    let mut order: HashMap<String, usize> = HashMap::new();
    let mut next = 0usize;

    for line in data_records(tumor_vcf)? {
        let cols: Vec<&str> = line.split('\t').collect();
        let key = match key_of(&cols) {
            Some(k) => k,
            None => continue,
        };
        tumor_keys.insert(key);
        let origin = tumor_origin(cols.get(7).copied().unwrap_or("."));
        recs.push(make_rec(&cols, origin, &mut order, &mut next));
    }
    for line in data_records(normal_vcf)? {
        let cols: Vec<&str> = line.split('\t').collect();
        let key = match key_of(&cols) {
            Some(k) => k,
            None => continue,
        };
        if tumor_keys.contains(&key) {
            continue; // germline carried into tumor → already emitted as `shared`
        }
        recs.push(make_rec(&cols, "germline", &mut order, &mut next));
    }

    recs.sort_by(|a, b| {
        order[&a.contig]
            .cmp(&order[&b.contig])
            .then(a.pos.cmp(&b.pos))
    });

    let f = File::create(out)?;
    let mut w = GzEncoder::new(BufWriter::new(f), Compression::default());
    for l in &meta {
        writeln!(w, "{l}")?;
    }
    if !declared.contains("EIDOLON_ORIGIN") {
        writeln!(w, "{EIDOLON_ORIGIN_DECL}")?;
    }
    for (id, decl) in SV_INFO_DECLS {
        if !declared.contains(*id) {
            writeln!(w, "{decl}")?;
        }
    }
    writeln!(w, "{chrom_line}")?;
    for r in &recs {
        writeln!(w, "{}", r.line)?;
    }
    w.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn tumor_origin_classifies_provenance() {
        assert_eq!(tumor_origin("EIDOLON_PROVENANCE=denovo"), "somatic");
        assert_eq!(
            tumor_origin("SVTYPE=DEL;EIDOLON_PROVENANCE=denovo"),
            "somatic"
        );
        // #405 reproductive somatic-VCF variants classify as somatic, not shared,
        // even though (like germline input) they came from a file.
        assert_eq!(tumor_origin("EIDOLON_PROVENANCE=somatic_input"), "somatic");
        assert_eq!(
            tumor_origin("AF=0.3;EIDOLON_PROVENANCE=somatic_input"),
            "somatic"
        );
        assert_eq!(tumor_origin("EIDOLON_PROVENANCE=input"), "shared");
        assert_eq!(tumor_origin("."), "shared");
    }

    #[test]
    fn append_origin_handles_empty_and_existing() {
        assert_eq!(append_origin(".", "somatic"), "EIDOLON_ORIGIN=somatic");
        assert_eq!(append_origin("", "germline"), "EIDOLON_ORIGIN=germline");
        assert_eq!(
            append_origin("SVTYPE=DEL", "shared"),
            "SVTYPE=DEL;EIDOLON_ORIGIN=shared"
        );
    }

    fn write_gz(path: &Path, body: &str) {
        let f = File::create(path).unwrap();
        let mut w = GzEncoder::new(f, Compression::default());
        w.write_all(body.as_bytes()).unwrap();
        w.finish().unwrap();
    }

    fn read_gz(path: &Path) -> String {
        let mut s = String::new();
        MultiGzDecoder::new(File::open(path).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        s
    }

    #[test]
    fn merge_classifies_germline_somatic_shared() {
        let dir = tempfile::tempdir().unwrap();
        let normal = dir.path().join("n.vcf.gz");
        let tumor = dir.path().join("t.vcf.gz");
        let out = dir.path().join("merged.vcf.gz");

        // normal: two germline SNPs (denovo in the normal pass).
        write_gz(
            &normal,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr1\t100\t.\tA\tG\t.\tPASS\tEIDOLON_PROVENANCE=denovo\tGT\t0/1\n\
             chr1\t200\t.\tC\tT\t.\tPASS\tEIDOLON_PROVENANCE=denovo\tGT\t0/1\n",
        );
        // tumor: pos100 carried (input→shared), pos200 dropped (germline-only),
        // pos300 somatic (denovo).
        write_gz(
            &tumor,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr1\t100\t.\tA\tG\t.\tPASS\tEIDOLON_PROVENANCE=input\tGT\t0/1\n\
             chr1\t300\t.\tT\tA\t.\tPASS\tEIDOLON_PROVENANCE=denovo\tGT\t1/1\n",
        );

        merge_goldens(&normal, &tumor, &out).unwrap();
        let body = read_gz(&out);

        // origin tags
        let line = |pos: &str| {
            body.lines()
                .find(|l| l.starts_with(&format!("chr1\t{pos}\t")))
                .unwrap()
        };
        assert!(
            line("100").contains("EIDOLON_ORIGIN=shared"),
            "pos100 = {}",
            line("100")
        );
        assert!(
            line("200").contains("EIDOLON_ORIGIN=germline"),
            "pos200 = {}",
            line("200")
        );
        assert!(
            line("300").contains("EIDOLON_ORIGIN=somatic"),
            "pos300 = {}",
            line("300")
        );
        // header declares EIDOLON_ORIGIN and is position-sorted (100<200<300)
        assert!(body.contains("##INFO=<ID=EIDOLON_ORIGIN"));
        let p100 = body.find("chr1\t100").unwrap();
        let p200 = body.find("chr1\t200").unwrap();
        let p300 = body.find("chr1\t300").unwrap();
        assert!(
            p100 < p200 && p200 < p300,
            "records must be position-sorted"
        );
    }
}
