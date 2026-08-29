//! Per-allele AF extraction and the two correlation statistics.
//!
//! Ported from `scripts/delta/scn_af_compare.py`. The comments that explain WHY each
//! rule exists are carried over deliberately — every one of them records a defect:
//! per-allele arity (#450's multi-allelic mismatch), the `<*>` non-ref placeholder that
//! must not disqualify the real base beside it, and zero-fill for a covered site with no
//! observed alt.

use std::collections::HashMap;

/// One VCF's alleles, keyed by (chrom, pos, ref, alt), in file order.
///
/// Insertion order is preserved to mirror the Python `OrderedDict` exactly. None of the
/// statistics depend on it, but keeping it means a diff of the two implementations'
/// output is a diff of the arithmetic, not of iteration order.
pub struct Alleles {
    pub keys: Vec<(String, String, String, String)>,
    pub map: HashMap<(String, String, String, String), (f64, Option<f64>)>,
}

impl Alleles {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            map: HashMap::new(),
        }
    }
    fn insert(&mut self, k: (String, String, String, String), v: (f64, Option<f64>)) {
        if self.map.insert(k.clone(), v).is_none() {
            self.keys.push(k);
        }
    }
    pub fn len(&self) -> usize {
        self.keys.len()
    }
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }
}

fn field_index(fmt: &str, key: &str) -> Option<usize> {
    fmt.split(':').position(|p| p == key)
}

/// Select one allele's value from a per-allele VCF field.
///
/// `kind` is `A` (one value per ALT, e.g. AF) or `R` (ref first then one per ALT, e.g.
/// AD). Returns `None` unless the field's arity actually matches the ALT count: a field
/// of the wrong length cannot be indexed safely, and guessing would silently attribute
/// one allele's number to another. The pre-#450 code took `[0]` unconditionally, which
/// reported the FIRST allele's fraction for every allele of a multi-allelic record.
fn pick_per_allele(raw: &str, alt_idx: usize, n_alts: usize, kind: char) -> Option<f64> {
    let parts: Vec<&str> = raw.split(',').collect();
    let want = if kind == 'A' { n_alts } else { n_alts + 1 };
    if parts.len() != want {
        return None;
    }
    let idx = if kind == 'A' { alt_idx } else { alt_idx + 1 };
    parts.get(idx)?.parse::<f64>().ok()
}

fn depth_from_ad(fmt: &str, sample: &str) -> Option<f64> {
    if fmt.is_empty() || sample.is_empty() {
        return None;
    }
    let ad_i = field_index(fmt, "AD")?;
    let vals: Vec<&str> = sample.split(':').collect();
    let raw = vals.get(ad_i)?;
    let mut total = 0.0;
    for p in raw.split(',') {
        total += p.parse::<f64>().ok()?;
    }
    Some(total)
}

/// Alt fraction for the ALT at 0-based `alt_idx`, and the site's total depth.
/// Order: FORMAT/AF, INFO/AF, then FORMAT/AD.
fn af_for_allele(
    info: &str,
    fmt: &str,
    sample: &str,
    alt_idx: usize,
    n_alts: usize,
) -> (Option<f64>, Option<f64>) {
    let depth = depth_from_ad(fmt, sample);
    if !fmt.is_empty()
        && !sample.is_empty()
        && let Some(af_i) = field_index(fmt, "AF")
    {
        let vals: Vec<&str> = sample.split(':').collect();
        if let Some(raw) = vals.get(af_i)
            && let Some(v) = pick_per_allele(raw, alt_idx, n_alts, 'A')
        {
            return (Some(v), depth);
        }
    }
    for kv in info.split(';') {
        if let Some(rest) = kv.strip_prefix("AF=")
            && let Some(v) = pick_per_allele(rest, alt_idx, n_alts, 'A')
        {
            return (Some(v), depth);
        }
    }
    // FORMAT/AD (Number=R). This is the path #450 depends on: reading the observed side
    // straight from `bcftools mpileup` means AD carries one count per allele, so this
    // allele's own count over the site total is its observed fraction.
    if !fmt.is_empty()
        && !sample.is_empty()
        && let Some(ad_i) = field_index(fmt, "AD")
    {
        let vals: Vec<&str> = sample.split(':').collect();
        if let Some(raw) = vals.get(ad_i) {
            let parts: Vec<&str> = raw.split(',').collect();
            if parts.len() == n_alts + 1 {
                let mut counts = Vec::with_capacity(parts.len());
                for p in &parts {
                    match p.parse::<f64>() {
                        Ok(v) => counts.push(v),
                        Err(_) => return (None, depth),
                    }
                }
                let total: f64 = counts.iter().sum();
                if total > 0.0 {
                    return (Some(counts[alt_idx + 1] / total), Some(total));
                }
            }
        }
    }
    (None, depth)
}

/// `(alleles, sites)` from a VCF. Multi-allelic records are expanded one entry per
/// LITERAL ALT — that expansion is what makes the mpileup-derived observed side
/// joinable (#450): mpileup reports `ALT=T,<*>` with `AD=93,7,0`, so the truth's single
/// ALT `T` must be matched against the ALT *list* and its own AD element selected.
/// `sites` maps (chrom, pos, ref) -> total depth, so a truth ALT the observed side never
/// listed can be scored as 0.0 observed rather than dropped.
pub fn load_af(text: &str) -> (Alleles, HashMap<(String, String, String), f64>) {
    let mut alleles = Alleles::new();
    let mut sites = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let (chrom, pos, refb, alt_field) = (f[0], f[1], f[3], f[4]);
        let info = f.get(7).copied().unwrap_or(".");
        let fmt = f.get(8).copied().unwrap_or("");
        let sample = f.get(9).copied().unwrap_or("");
        let alts: Vec<&str> = alt_field.split(',').collect();
        if let Some(d) = depth_from_ad(fmt, sample) {
            sites.insert((chrom.to_string(), pos.to_string(), refb.to_string()), d);
        }
        for (i, alt) in alts.iter().enumerate() {
            // Skip symbolic / breakend ALTERNATIVES but keep literal ones from the same
            // record: mpileup's `<*>` non-ref placeholder sits beside a real base, and
            // dropping the whole record over it was the bug.
            if alt.starts_with('<') || alt.contains('[') || alt.contains(']') || *alt == "." {
                continue;
            }
            let (af, dp) = af_for_allele(info, fmt, sample, i, alts.len());
            if let Some(af) = af {
                alleles.insert(
                    (
                        chrom.to_string(),
                        pos.to_string(),
                        refb.to_string(),
                        alt.to_string(),
                    ),
                    (af, dp),
                );
            }
        }
    }
    (alleles, sites)
}

pub fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    if n < 2 {
        return f64::NAN;
    }
    let (mx, my) = (
        xs.iter().sum::<f64>() / n as f64,
        ys.iter().sum::<f64>() / n as f64,
    );
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let syy: f64 = ys.iter().map(|y| (y - my).powi(2)).sum();
    if sxx == 0.0 || syy == 0.0 {
        return f64::NAN;
    }
    sxy / (sxx * syy).sqrt()
}

/// Average ranks, ties shared — matching the Python implementation exactly.
fn ranks(vals: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..vals.len()).collect();
    order.sort_by(|&a, &b| {
        vals[a]
            .partial_cmp(&vals[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut out = vec![0.0; vals.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && vals[order[j + 1]] == vals[order[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0;
        for k in i..=j {
            out[order[k]] = avg;
        }
        i = j + 1;
    }
    out
}

pub fn spearman(xs: &[f64], ys: &[f64]) -> f64 {
    if xs.len() < 2 {
        return f64::NAN;
    }
    pearson(&ranks(xs), &ranks(ys))
}

/// Python's `%g` for the one place the output uses it (`min-depth {:g}`).
pub fn fmt_g(v: f64) -> String {
    if v == v.trunc() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        let s = format!("{v}");
        s
    }
}
