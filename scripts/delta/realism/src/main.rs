//! `realism-panel` — measure how far eidolon's reads sit from real sequencing data.
//!
//! A Delta validation helper, not part of the shipped tool. See `Cargo.toml` for why.
//!
//! WHY IT EXISTS. Measured on matched chr22 sequence at 30x, real NA12878 against eidolon
//! 3.2.0:
//!
//! ```text
//!                                  real      simulated
//!   candidate breakpoints / Mb     654.5         0.0
//!   improper pairs                 3.15%        0.00%
//!   soft clip >= 20 bp             0.67%        0.00%
//!   depth VMR                    5.5-8.9      1.03-1.12
//! ```
//!
//! Zero, not "fewer". A caller tuned on eidolon calibrates its false-positive filters against
//! an empty background and then meets 654 candidates per megabase on real data. Which means
//! eidolon-measured PRECISION has never predicted real precision — recall is a separate
//! question, since planted events are real signal.
//!
//! SCOPE. Somatic simulation is human (COSMIC/PCAWG models, human tissue types), so realism
//! is measured against real human data. Germline is species-agnostic and its correctness is
//! checked a different way — across deliberately varied contig shapes, which is where #625 and
//! #607's H1N1 artifacts were caught. Those are correctness questions, not realism ones.
//!
//! USAGE
//!
//! ```text
//!   realism-panel --bam <file> --regions <bed> [--label NAME]
//!       [--min-clip 20] [--min-support 3] [--max-tlen 2000] [--depth-lag 500]
//! ```
//!
//! Emits one TSV row per region on stdout. Both sides of a comparison must be run with
//! IDENTICAL thresholds — a gap measured with different settings on each side is measuring
//! the settings. The wrapper script passes one set to both.

mod metrics;
mod reader;

use reader::{Region, measure};
use std::path::PathBuf;
use std::process::ExitCode;

struct Args {
    bam: PathBuf,
    regions: PathBuf,
    label: String,
    min_clip: usize,
    min_support: usize,
    max_tlen: i64,
    depth_lag: usize,
}

fn usage() -> &'static str {
    "usage: realism-panel --bam <file> --regions <bed> [--label NAME] \
     [--min-clip N] [--min-support N] [--max-tlen N] [--depth-lag N]"
}

fn parse_args() -> Result<Args, String> {
    let mut bam = None;
    let mut regions = None;
    let mut label = "unlabelled".to_string();
    let mut min_clip = 20usize;
    let mut min_support = 3usize;
    let mut max_tlen = 2000i64;
    let mut depth_lag = 500usize;

    let mut it = std::env::args().skip(1);
    while let Some(a) = it.next() {
        let mut val = || it.next().ok_or_else(|| format!("{a} needs a value"));
        match a.as_str() {
            "--bam" => bam = Some(PathBuf::from(val()?)),
            "--regions" => regions = Some(PathBuf::from(val()?)),
            "--label" => label = val()?,
            "--min-clip" => min_clip = val()?.parse().map_err(|e| format!("--min-clip: {e}"))?,
            "--min-support" => {
                min_support = val()?.parse().map_err(|e| format!("--min-support: {e}"))?
            }
            "--max-tlen" => max_tlen = val()?.parse().map_err(|e| format!("--max-tlen: {e}"))?,
            "--depth-lag" => depth_lag = val()?.parse().map_err(|e| format!("--depth-lag: {e}"))?,
            "-h" | "--help" => return Err(usage().to_string()),
            other => return Err(format!("unknown argument {other}\n{}", usage())),
        }
    }
    Ok(Args {
        bam: bam.ok_or_else(|| format!("--bam is required\n{}", usage()))?,
        regions: regions.ok_or_else(|| format!("--regions is required\n{}", usage()))?,
        label,
        min_clip,
        min_support,
        max_tlen,
        depth_lag,
    })
}

/// BED: `contig<TAB>start<TAB>end`, 0-based half-open. Blank lines and `#` comments skipped.
fn read_regions(path: &PathBuf) -> Result<Vec<Region>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            return Err(format!(
                "{}:{}: expected 3 tab-separated fields",
                path.display(),
                i + 1
            ));
        }
        let start: usize = f[1]
            .parse()
            .map_err(|e| format!("{}:{}: start: {e}", path.display(), i + 1))?;
        let end: usize = f[2]
            .parse()
            .map_err(|e| format!("{}:{}: end: {e}", path.display(), i + 1))?;
        if end <= start {
            return Err(format!(
                "{}:{}: end must exceed start",
                path.display(),
                i + 1
            ));
        }
        out.push(Region {
            contig: f[0].to_string(),
            start,
            end,
        });
    }
    if out.is_empty() {
        // An empty region list would emit a header and no rows, which reads as "measured,
        // nothing found" rather than "measured nothing" (rule 4).
        return Err(format!("{}: no regions", path.display()));
    }
    Ok(out)
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{e}");
            return ExitCode::FAILURE;
        }
    };
    let regions = match read_regions(&args.regions) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{e}");
            return ExitCode::FAILURE;
        }
    };
    let measured = match measure(
        &args.bam,
        &regions,
        args.min_clip,
        args.min_support,
        args.max_tlen,
        args.depth_lag,
    ) {
        Ok(m) => m,
        Err(e) => {
            // Every error from `measure` is a refusal to report an unmeasured region as a
            // clean one. Exiting non-zero keeps that from becoming a missing TSV row nobody
            // notices.
            eprintln!("realism-panel: {e}");
            return ExitCode::FAILURE;
        }
    };

    println!(
        "label\tcontig\tstart\tend\treads\tspan_bp\tcand_bp\tcand_per_mb\timproper_pct\
         \tclip_pct\tmapq0_pct\tdepth_mean\tdepth_vmr\tdepth_excess\tdepth_acf\tins_n\tins_mean\tins_sd\
         \tins_skew\tins_p99"
    );
    for (r, m) in regions.iter().zip(measured.iter()) {
        // depth_excess is the cross-dataset number: VMR is not comparable between BAMs at
        // different depths. See DepthStats::excess_dispersion.
        let (dm, dv, dx, da) = match &m.depth {
            Some(d) => (d.mean, d.vmr, d.excess_dispersion(), d.autocorr_500),
            None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
        };
        let (n, im, isd, isk, ip) = match &m.insert {
            Some(s) => (s.n as f64, s.mean, s.sd, s.skew, s.p99 as f64),
            None => (0.0, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
        };
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.4}\t{:.4}\t{:.4}\t{:.2}\t{:.3}\t{:.5}\
             \t{:.3}\t{:.0}\t{:.1}\t{:.1}\t{:+.3}\t{:.0}",
            args.label,
            r.contig,
            r.start,
            r.end,
            m.reads,
            m.span_bp,
            m.candidate_breakpoints,
            m.candidates_per_mb(),
            m.improper_pair_rate(),
            m.clip_rate(),
            m.mapq0_rate(),
            dm,
            dv,
            dx,
            da,
            n,
            im,
            isd,
            isk,
            ip
        );
    }
    ExitCode::SUCCESS
}
