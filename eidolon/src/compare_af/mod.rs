//! `eidolon compare-af` — per-allele AF correlation between a truth and a simulated VCF.
//!
//! A port of `scripts/delta/scn_af_compare.py`, which determined the numbers §3.12 cites
//! while being exercised by nothing automatic (#466). Two reasons to move it into Rust
//! rather than bolt tests onto the Python:
//!
//!   * it gets CI coverage for free, which was the whole point of #466; and
//!   * it needs no external tool — it parses VCFs and does arithmetic — so under this
//!     repo's standing rule it should not have been Python in the first place.
//!     `sbs96_compare.py` parses SigProfiler's output and legitimately stays.
//!
//! Output is byte-identical to the Python by design: `run_subclonal_vaf_validation.sh`
//! parses specific lines out of it to build its PASS/FAIL verdict, and a differential
//! test pins that equality rather than trusting it.

pub mod errors;
pub mod utils;

use errors::CompareAfError;
use std::collections::HashMap;
use std::path::Path;
use utils::{fmt_g, load_af, pearson, spearman};

fn read_maybe_gzip(path: &Path) -> Result<String, CompareAfError> {
    let bytes = std::fs::read(path).map_err(|source| CompareAfError::Io {
        path: path.display().to_string(),
        source,
    })?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        use flate2::read::MultiGzDecoder;
        use std::io::Read;
        let mut s = String::new();
        MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .map_err(|source| CompareAfError::Io {
                path: path.display().to_string(),
                source,
            })?;
        Ok(s)
    } else {
        Ok(String::from_utf8_lossy(&bytes).into_owned())
    }
}

/// Python's `f"{x:.1%}"`.
fn pct1(x: f64) -> String {
    format!("{:.1}%", x * 100.0)
}

fn bin_of(af: f64) -> usize {
    ((af * 10.0) as usize).min(9)
}

pub fn run(
    truth: &Path,
    sim: &Path,
    min_depth: f64,
    max_uncovered_frac: f64,
) -> Result<(), CompareAfError> {
    let (a, _a_sites) = load_af(&read_maybe_gzip(truth)?);
    let (mut b, b_sites) = load_af(&read_maybe_gzip(sim)?);

    // A truth ALT the observed side never listed, at a position it DID cover, means zero
    // reads carried that allele — an observed fraction of 0.0, which is a measurement.
    // Dropping those would recreate exactly the VAF-dependent exclusion #450 fixed, just
    // at a lower threshold: the sites most likely to have no alt reads at all are the
    // lowest-VAF ones, so excluding them biases the result optimistically.
    let mut shared: Vec<(String, String, String, String)> = Vec::new();
    let mut zero_filled = 0usize;
    for k in &a.keys {
        if b.map.contains_key(k) {
            shared.push(k.clone());
            continue;
        }
        let site = (k.0.clone(), k.1.clone(), k.2.clone());
        if let Some(&dp) = b_sites.get(&site)
            && dp > 0.0
        {
            b.map.insert(k.clone(), (0.0, Some(dp)));
            b.keys.push(k.clone());
            shared.push(k.clone());
            zero_filled += 1;
        }
    }
    let only_a = a.len() - shared.len();
    let only_b = b.len() - shared.len();

    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    let mut gated = 0usize;
    for k in &shared {
        let (af_a, dp_a) = a.map[k];
        let (af_b, dp_b) = b.map[k];
        if min_depth > 0.0
            && (dp_a.is_some_and(|d| d < min_depth) || dp_b.is_some_and(|d| d < min_depth))
        {
            gated += 1;
            continue;
        }
        xs.push(af_a);
        ys.push(af_b);
    }

    println!(
        "truth sites={}  sim sites={}  shared={}  (only-truth={only_a}, only-sim={only_b})",
        a.len(),
        b.len(),
        shared.len()
    );
    if zero_filled > 0 {
        println!(
            "  {zero_filled} truth allele(s) had coverage but zero observed reads \
             -> scored as observed AF 0.0 (not dropped)"
        );
    }
    if only_a > 0 {
        eprintln!(
            "  {only_a} truth allele(s) had NO coverage on the observed side and are \
             excluded ({} of the planted set).",
            pct1(only_a as f64 / a.len() as f64)
        );
    }
    if min_depth > 0.0 {
        println!(
            "min-depth {}: {gated} shared sites gated, {} compared",
            fmt_g(min_depth),
            xs.len()
        );
    }
    if xs.len() < 2 {
        return Err(CompareAfError::Fatal(
            "fewer than 2 comparable sites — check inputs / --min-depth".to_string(),
        ));
    }

    let diffs: Vec<f64> = xs.iter().zip(&ys).map(|(x, y)| y - x).collect();
    let n = diffs.len() as f64;
    let mae = diffs.iter().map(|d| d.abs()).sum::<f64>() / n;
    let rmse = (diffs.iter().map(|d| d * d).sum::<f64>() / n).sqrt();
    println!(
        "n={}  Pearson r={:.4}  Spearman rho={:.4}",
        xs.len(),
        pearson(&xs, &ys),
        spearman(&xs, &ys)
    );
    println!(
        "MAE={mae:.4}  RMSE={rmse:.4}  mean(sim-truth)={:+.4}",
        diffs.iter().sum::<f64>() / n
    );
    println!("  target: r>=0.95 and per-bin MAE within AF-estimation noise at that coverage");

    let mut planted: HashMap<usize, usize> = HashMap::new();
    for k in &a.keys {
        *planted.entry(bin_of(a.map[k].0)).or_insert(0) += 1;
    }

    println!("per-AF-decile coverage and MAE (truth bin):");
    println!(
        "  {:<12} {:>8} {:>8} {:>9}   MAE",
        "bin", "planted", "scored", "unscored"
    );
    let mut shortfall = 0usize;
    let total_planted: usize = planted.values().sum();
    for i in 0..10usize {
        let n_p = *planted.get(&i).unwrap_or(&0);
        if n_p == 0 {
            continue;
        }
        let (lo, hi) = (i as f64 / 10.0, i as f64 / 10.0 + 0.1);
        let bd: Vec<f64> = xs
            .iter()
            .zip(&ys)
            .filter(|(x, _)| bin_of(**x) == i)
            .map(|(x, y)| (y - x).abs())
            .collect();
        let n_s = bd.len();
        shortfall += n_p - n_s;
        let mae_s = if bd.is_empty() {
            "n/a — NOTHING SCORED".to_string()
        } else {
            format!("{:.4}", bd.iter().sum::<f64>() / bd.len() as f64)
        };
        let flag = if n_p == n_s { "" } else { "  <-- incomplete" };
        println!(
            "  [{lo:.1},{hi:.1})   {n_p:>8} {n_s:>8} {:>9}   {mae_s}{flag}",
            n_p - n_s
        );
    }
    println!(
        "  total planted={total_planted}  scored={}  unscored={shortfall} ({})",
        xs.len(),
        pct1(shortfall as f64 / total_planted as f64)
    );

    // Enforced, not advised: a bias/MAE computed over a subset that silently omits the
    // lowest-VAF stratum is exactly the #450 failure, and it read as a clean PASS.
    let uncovered = (a.len() - xs.len()) as f64 / a.len() as f64;
    if uncovered > max_uncovered_frac {
        return Err(CompareAfError::Fatal(format!(
            "FAIL: {} of planted truth alleles went unscored (limit {}). The reported \
             bias/MAE cover a subset of the planted set and must not be quoted. \
             Attribute the shortfall per stratum above before believing any of it.",
            pct1(uncovered),
            pct1(max_uncovered_frac)
        )));
    }
    Ok(())
}
