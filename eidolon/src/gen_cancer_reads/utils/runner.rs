//! Orchestrates the tumor/normal two-pass simulation: run the normal pass, feed
//! its golden VCF as the tumor pass's germline, run the tumor pass, then merge
//! the FASTQs (read-tagged) and the golden VCFs (origin-tagged). Mirrors
//! `tools/cancer_simulate.sh` but native — no shell-out, no bcftools/awk.

use std::path::{Path, PathBuf};

use eidolon_core::rng::NeatRng;
use log::{info, warn};

use crate::gen_cancer_reads::errors::GenCancerReadsError;
use crate::gen_cancer_reads::utils::config::CancerConfig;
use crate::gen_cancer_reads::utils::{fastq_merge, vcf_merge};
use crate::gen_reads::utils::config::RunConfiguration;
use crate::gen_reads::utils::runner::run_neat;

pub fn run_cancer(cfg: &CancerConfig) -> Result<(), GenCancerReadsError> {
    let (normal_cov, tumor_cov) = cfg.per_pass_coverage();

    // ── Pass 1: normal ──────────────────────────────────────────────────
    let normal = cfg.normal_pass()?;
    info!(
        ">> Pass 1 (normal): {normal_cov}x coverage → {}",
        normal.output_filename
    );
    let mut rng_n = NeatRng::new_from_seed(&normal.seed_vec)
        .map_err(|e| GenCancerReadsError::Rng(format!("{e:?}")))?;
    run_neat(&normal, &mut rng_n).map_err(|e| GenCancerReadsError::GenReadsPass {
        pass: "normal",
        source: e,
    })?;

    // Shared germline: user-supplied, else the normal pass's golden VCF. This is
    // the crux — the tumor cells must carry the same germline as the normal.
    let germline_vcf = match &cfg.germline_vcf {
        Some(p) => p.clone(),
        None => normal.output_vcf.clone().ok_or_else(|| {
            GenCancerReadsError::MergeError("normal pass produced no golden VCF".into())
        })?,
    };
    if !germline_vcf.is_file() {
        return Err(GenCancerReadsError::MergeError(format!(
            "germline VCF not found after normal pass: {germline_vcf:?}"
        )));
    }

    // ── Pass 2: tumor ───────────────────────────────────────────────────
    let tumor = cfg.tumor_pass(germline_vcf)?;
    match cfg.tumor_mutation_rate {
        Some(r) => info!(
            ">> Pass 2 (tumor): {tumor_cov}x coverage, somatic rate {r} → {}",
            tumor.output_filename
        ),
        None => info!(
            ">> Pass 2 (tumor): {tumor_cov}x coverage, model-fitted rate → {}",
            tumor.output_filename
        ),
    }
    let mut rng_t = NeatRng::new_from_seed(&tumor.seed_vec)
        .map_err(|e| GenCancerReadsError::Rng(format!("{e:?}")))?;
    run_neat(&tumor, &mut rng_t).map_err(|e| GenCancerReadsError::GenReadsPass {
        pass: "tumor",
        source: e,
    })?;

    // ── Merge: tagged FASTQs + origin-tagged truth VCF ──────────────────
    merge_fastqs(cfg, &normal, &tumor)?;
    let truth = cfg
        .output_dir
        .join(format!("{}_merged_truth.vcf.gz", cfg.output_prefix));
    let normal_vcf = normal.output_vcf.clone().expect("normal golden vcf path");
    let tumor_vcf = tumor.output_vcf.clone().expect("tumor golden vcf path");
    vcf_merge::merge_goldens(&normal_vcf, &tumor_vcf, &truth)?;
    info!(">> Wrote origin-tagged truth VCF: {truth:?}");

    if !cfg.keep_per_pass {
        cleanup_per_pass(&normal, &tumor);
    }

    print_summary(cfg, &truth);
    Ok(())
}

fn merge_fastqs(
    cfg: &CancerConfig,
    normal: &RunConfiguration,
    tumor: &RunConfiguration,
) -> Result<(), GenCancerReadsError> {
    let need = |p: &Option<PathBuf>, what: &str| -> Result<PathBuf, GenCancerReadsError> {
        p.clone()
            .ok_or_else(|| GenCancerReadsError::MergeError(format!("missing {what} path")))
    };

    let n_r1 = need(&normal.output_fastq_1, "normal R1")?;
    let t_r1 = need(&tumor.output_fastq_1, "tumor R1")?;
    let merged_r1 = cfg
        .output_dir
        .join(format!("{}_merged_r1.fastq.gz", cfg.output_prefix));
    fastq_merge::tag_and_concat(&[(&n_r1, "N"), (&t_r1, "T")], &merged_r1)?;
    info!(">> Wrote merged R1: {merged_r1:?}");

    if cfg.paired_ended {
        let n_r2 = need(&normal.output_fastq_2, "normal R2")?;
        let t_r2 = need(&tumor.output_fastq_2, "tumor R2")?;
        let merged_r2 = cfg
            .output_dir
            .join(format!("{}_merged_r2.fastq.gz", cfg.output_prefix));
        fastq_merge::tag_and_concat(&[(&n_r2, "N"), (&t_r2, "T")], &merged_r2)?;
        info!(">> Wrote merged R2: {merged_r2:?}");
    }
    Ok(())
}

/// Remove the per-pass FASTQs (keep the per-pass golden VCFs — cheap truth).
fn cleanup_per_pass(normal: &RunConfiguration, tumor: &RunConfiguration) {
    for p in [
        &normal.output_fastq_1,
        &normal.output_fastq_2,
        &tumor.output_fastq_1,
        &tumor.output_fastq_2,
    ]
    .into_iter()
    .flatten()
    {
        if let Err(e) = std::fs::remove_file(p) {
            warn!("could not remove per-pass FASTQ {p:?}: {e}");
        }
    }
}

fn print_summary(cfg: &CancerConfig, truth: &Path) {
    let (n, t) = cfg.per_pass_coverage();
    info!("──────────────────────────────────────────────");
    info!("Cancer simulation complete");
    info!(
        "  total coverage {}x, purity {} (normal {n}x / tumor {t}x)",
        cfg.total_coverage, cfg.purity
    );
    info!(
        "  merged FASTQ:  {}_merged_r1.fastq.gz{}",
        cfg.output_prefix,
        if cfg.paired_ended { " (+ _r2)" } else { "" }
    );
    info!("  merged truth:  {truth:?} (INFO/EIDOLON_ORIGIN = germline | somatic | shared)");
    info!("──────────────────────────────────────────────");
}
