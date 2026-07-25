//! Configuration for the native cancer (tumor/normal) read simulator.
//!
//! A `CancerConfig` is parsed from a single YAML file and derives the two
//! `RunConfiguration`s that drive the normal and tumor passes — no temp YAML
//! round-trip (unlike `tools/cancer_simulate.sh`, which shells out twice).

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use log::warn;
use serde_yml::Value;

use crate::gen_cancer_reads::errors::GenCancerReadsError;
use crate::gen_reads::utils::config::RunConfiguration;
use crate::gen_reads::utils::subclone::{Subclone, SubcloneModel};

/// Default tumor-pass somatic SNP/indel rate (typical solid tumor; see #235).
/// The de-novo mutations added in the tumor pass are somatic, so this — not the
/// model's (often corpus-aggregated) fitted rate — is the sensible default.
const DEFAULT_TUMOR_MUTATION_RATE: f64 = 1e-5;

#[derive(Debug, Clone)]
pub struct CancerConfig {
    pub reference: PathBuf,
    pub output_dir: PathBuf,
    pub output_prefix: String,
    pub total_coverage: usize,
    /// Tumor cell fraction in (0, 1). Tumor pass gets `purity * total`, normal
    /// pass the rest.
    pub purity: f64,
    pub read_len: usize,
    pub paired_ended: bool,
    pub fragment_mean: Option<f64>,
    pub fragment_st_dev: Option<f64>,
    /// A built fragment-length model. Like gen-reads (#355), it satisfies the
    /// paired-end fragment-source requirement on its own and takes precedence over
    /// fragment_mean/st_dev at runtime.
    pub fragment_model: Option<PathBuf>,
    pub rng_seed_root: String,
    pub normal_model: Option<PathBuf>,
    pub tumor_model: Option<PathBuf>,
    pub normal_mutation_rate: Option<f64>,
    /// `Some(r)` → use `r`; `None` → defer to the tumor model's fitted rate
    /// (YAML `tumor_mutation_rate: model`). Defaults to `Some(1e-5)`.
    pub tumor_mutation_rate: Option<f64>,
    /// Shared germline VCF. If absent, pass 1 generates one de novo and pass 2
    /// consumes it (guaranteeing tumor cells carry the same germline as normal).
    pub germline_vcf: Option<PathBuf>,
    pub sv_rate_scale: f64,
    pub keep_per_pass: bool,
    pub overwrite_output: bool,
    /// Optional subclonal architecture (#405). When set, the tumor pass distributes
    /// its de-novo somatic variants across these subclones; each subclone's
    /// cancer-cell fraction (CCF) composes with the variant's dosage, so the
    /// observed VAF in the merged output is `purity × dosage × CCF`. `None` → the
    /// pre-#405 behavior (somatic VAF driven by dosage and purity alone).
    pub subclones: Option<SubcloneModel>,
    /// Reproductive somatic input (#405): a VCF of somatic variants to replay in the
    /// tumor pass, honored at their observed VAF (`INFO/AF` or `FORMAT/AD`). Each
    /// variant's fraction is divided by `purity` so it reproduces after tumor/normal
    /// mixing, and it is tagged `somatic` in the merged truth. Composes with de-novo
    /// somatic generation (set `tumor_mutation_rate: 0` for pure replay).
    pub somatic_vcf: Option<PathBuf>,
}

impl Default for CancerConfig {
    fn default() -> Self {
        CancerConfig {
            reference: PathBuf::new(),
            output_dir: PathBuf::from("."),
            output_prefix: "neat_cancer".to_string(),
            total_coverage: 30,
            purity: 0.5,
            read_len: 151,
            paired_ended: false,
            fragment_mean: None,
            fragment_st_dev: None,
            fragment_model: None,
            rng_seed_root: "cancer-simulate".to_string(),
            normal_model: None,
            tumor_model: None,
            normal_mutation_rate: None,
            tumor_mutation_rate: Some(DEFAULT_TUMOR_MUTATION_RATE),
            germline_vcf: None,
            sv_rate_scale: 0.0,
            keep_per_pass: true,
            overwrite_output: false,
            subclones: None,
            somatic_vcf: None,
        }
    }
}

impl CancerConfig {
    pub fn from_yaml_file(yaml_file: &PathBuf) -> Result<CancerConfig, GenCancerReadsError> {
        let file = fs::File::open(yaml_file).map_err(GenCancerReadsError::Io)?;
        let scrape: HashMap<String, Value> = serde_yml::from_reader(file)
            .map_err(|e| GenCancerReadsError::ConfigError(format!("YAML parse error: {e}")))?;
        Self::from_scrape(scrape)
    }

    fn from_scrape(scrape: HashMap<String, Value>) -> Result<CancerConfig, GenCancerReadsError> {
        let mut cfg = CancerConfig::default();

        // A subclonal architecture may be given inline (`subclones:`) OR loaded from a
        // deconvolution-tool cluster table (`subclones_file:`), but not both.
        if scrape.contains_key("subclones") && scrape.contains_key("subclones_file") {
            return Err(GenCancerReadsError::ConfigError(
                "specify either `subclones` (inline) or `subclones_file` (path), not both".into(),
            ));
        }

        let req_path = |v: &Value, key: &str| -> Result<PathBuf, GenCancerReadsError> {
            let s = v.as_str().ok_or_else(|| {
                GenCancerReadsError::ConfigError(format!("{key} must be a path string"))
            })?;
            Ok(PathBuf::from(s))
        };

        for (key, value) in &scrape {
            if value.as_str() == Some(".") {
                continue; // dot = keep default (matches gen-reads convention)
            }
            match key.as_str() {
                "reference" => cfg.reference = req_path(value, "reference")?,
                "output_dir" => cfg.output_dir = req_path(value, "output_dir")?,
                "output_prefix" | "output_filename" => {
                    cfg.output_prefix = value
                        .as_str()
                        .ok_or_else(|| {
                            GenCancerReadsError::ConfigError(
                                "output_prefix must be a string".into(),
                            )
                        })?
                        .to_string();
                }
                "total_coverage" | "coverage" => {
                    cfg.total_coverage = as_usize(value, "total_coverage")?;
                }
                "purity" => cfg.purity = as_f64(value, "purity")?,
                "read_len" => cfg.read_len = as_usize(value, "read_len")?,
                "paired_ended" => cfg.paired_ended = as_bool(value, "paired_ended")?,
                "fragment_mean" => cfg.fragment_mean = Some(as_f64(value, "fragment_mean")?),
                "fragment_st_dev" => cfg.fragment_st_dev = Some(as_f64(value, "fragment_st_dev")?),
                "fragment_model" => cfg.fragment_model = Some(req_path(value, "fragment_model")?),
                "rng_seed" | "rng_seed_root" => {
                    cfg.rng_seed_root = value
                        .as_str()
                        .map(str::to_string)
                        .or_else(|| value.as_i64().map(|n| n.to_string()))
                        .ok_or_else(|| {
                            GenCancerReadsError::ConfigError("rng_seed must be a string/int".into())
                        })?;
                }
                "normal_model" => cfg.normal_model = Some(req_path(value, "normal_model")?),
                "tumor_model" => cfg.tumor_model = Some(req_path(value, "tumor_model")?),
                "normal_mutation_rate" => {
                    cfg.normal_mutation_rate = Some(as_f64(value, "normal_mutation_rate")?);
                }
                "tumor_mutation_rate" => {
                    // `model` sentinel → defer to the tumor model's fitted rate.
                    cfg.tumor_mutation_rate = if value.as_str() == Some("model") {
                        None
                    } else {
                        Some(as_f64(value, "tumor_mutation_rate")?)
                    };
                }
                "germline_vcf" => cfg.germline_vcf = Some(req_path(value, "germline_vcf")?),
                "somatic_vcf" => cfg.somatic_vcf = Some(req_path(value, "somatic_vcf")?),
                "subclones" => cfg.subclones = Some(parse_subclones(value)?),
                "subclones_file" => {
                    cfg.subclones =
                        Some(parse_subclones_file(&req_path(value, "subclones_file")?)?);
                }
                "sv_rate_scale" => cfg.sv_rate_scale = as_f64(value, "sv_rate_scale")?,
                "keep_per_pass" => cfg.keep_per_pass = as_bool(value, "keep_per_pass")?,
                "overwrite_output" => cfg.overwrite_output = as_bool(value, "overwrite_output")?,
                _ => continue,
            }
        }

        cfg.validate()?;
        // Create the output directory up front (unlike gen-reads, the cancer path
        // builds its per-pass RunConfigurations directly and bypasses the YAML
        // parser's dir-creation arm). create_dir_all is a no-op if it exists.
        std::fs::create_dir_all(&cfg.output_dir).map_err(GenCancerReadsError::Io)?;
        Ok(cfg)
    }

    fn validate(&self) -> Result<(), GenCancerReadsError> {
        if !self.reference.is_file() {
            return Err(GenCancerReadsError::ConfigError(format!(
                "reference not found: {:?}",
                self.reference
            )));
        }
        if !(self.purity > 0.0 && self.purity < 1.0) {
            return Err(GenCancerReadsError::PurityOutOfRange(self.purity));
        }
        // A paired-end run needs a fragment-length source: a fragment_model, OR both
        // fragment_mean and fragment_st_dev. Mirrors gen-reads check_and_log_config
        // (#355) — the model takes precedence at runtime, so requiring mean/st_dev when
        // a model is present would reject a config the runtime handles fine.
        if self.paired_ended
            && self.fragment_model.is_none()
            && (self.fragment_mean.is_none() || self.fragment_st_dev.is_none())
        {
            return Err(GenCancerReadsError::ConfigError(
                "paired_ended requires a fragment_model, or both fragment_mean and fragment_st_dev"
                    .into(),
            ));
        }
        let (n, t) = self.per_pass_coverage();
        if n < 1 || t < 1 {
            return Err(GenCancerReadsError::PerPassCoverageZero {
                normal: n,
                tumor: t,
            });
        }
        Ok(())
    }

    /// (normal, tumor) per-pass integer coverage, rounded to nearest.
    pub fn per_pass_coverage(&self) -> (usize, usize) {
        let total = self.total_coverage as f64;
        let normal = (total * (1.0 - self.purity)).round() as usize;
        let tumor = (total * self.purity).round() as usize;
        (normal, tumor)
    }

    /// Fields common to both passes.
    fn shared_run_config(&self) -> RunConfiguration {
        RunConfiguration {
            reference: self.reference.clone(),
            read_len: self.read_len,
            ploidy: 2,
            paired_ended: self.paired_ended,
            fragment_mean: self.fragment_mean,
            fragment_st_dev: self.fragment_st_dev,
            fragment_model: self.fragment_model.clone(),
            produce_fastq: true,
            produce_vcf: true,
            overwrite_output: self.overwrite_output,
            output_dir: self.output_dir.clone(),
            ..Default::default()
        }
    }

    /// Normal pass: `(1-purity)*total` coverage, no de-novo SVs, germline-only
    /// golden VCF. Output stem `<prefix>_normal`.
    pub fn normal_pass(&self) -> Result<RunConfiguration, GenCancerReadsError> {
        let (normal_cov, _) = self.per_pass_coverage();
        let mut c = RunConfiguration {
            coverage: normal_cov,
            mutation_rate: self.normal_mutation_rate,
            mutation_model: self.normal_model.clone(),
            input_vcf: self.germline_vcf.clone(),
            sv_rate_scale: 0.0,
            output_filename: format!("{}_normal", self.output_prefix),
            rng_seed: Some(format!("{}-normal", self.rng_seed_root)),
            ..self.shared_run_config()
        };
        finalize(&mut c)?;
        Ok(c)
    }

    /// Tumor pass: `purity*total` coverage, somatic SNP/indel + SV on top of the
    /// shared germline (`input_vcf`). Output stem `<prefix>_tumor`.
    pub fn tumor_pass(
        &self,
        germline_vcf: PathBuf,
    ) -> Result<RunConfiguration, GenCancerReadsError> {
        let (_, tumor_cov) = self.per_pass_coverage();
        let mut c = RunConfiguration {
            coverage: tumor_cov,
            mutation_rate: self.tumor_mutation_rate,
            mutation_model: self.tumor_model.clone(),
            input_vcf: Some(germline_vcf),
            sv_rate_scale: self.sv_rate_scale,
            // #405: only the tumor pass carries the subclonal architecture — the
            // normal pass has no somatic variants to stamp.
            subclone_model: self.subclones.clone(),
            // #405 reproductive: replay the somatic VCF only in the tumor pass. Its
            // observed VAFs are divided by purity so they reproduce after mixing.
            somatic_vcf: self.somatic_vcf.clone(),
            somatic_af_scale: 1.0 / self.purity,
            // #405: record each somatic variant's intended observed VAF
            // (INFO/NEAT_VAF = purity × allele_fraction) for direct comparison
            // against a caller's VAF on the merged reads.
            merged_vaf_purity: Some(self.purity),
            output_filename: format!("{}_tumor", self.output_prefix),
            rng_seed: Some(format!("{}-tumor", self.rng_seed_root)),
            ..self.shared_run_config()
        };
        finalize(&mut c)?;
        Ok(c)
    }
}

/// Run `check_and_log_config` to populate derived output paths + seed_vec.
fn finalize(c: &mut RunConfiguration) -> Result<(), GenCancerReadsError> {
    RunConfiguration::check_and_log_config(c)
        .map_err(|e| GenCancerReadsError::ConfigError(format!("invalid derived pass config: {e}")))
}

fn as_usize(v: &Value, key: &str) -> Result<usize, GenCancerReadsError> {
    v.as_u64().map(|n| n as usize).ok_or_else(|| {
        GenCancerReadsError::ConfigError(format!("{key} must be a non-negative integer"))
    })
}

fn as_f64(v: &Value, key: &str) -> Result<f64, GenCancerReadsError> {
    v.as_f64()
        .or_else(|| v.as_i64().map(|n| n as f64))
        .ok_or_else(|| GenCancerReadsError::ConfigError(format!("{key} must be a number")))
}

fn as_bool(v: &Value, key: &str) -> Result<bool, GenCancerReadsError> {
    v.as_bool()
        .ok_or_else(|| GenCancerReadsError::ConfigError(format!("{key} must be a boolean")))
}

/// Parse the optional `subclones:` YAML list into a validated [`SubcloneModel`] (#405).
///
/// Expected shape — a non-empty sequence of `{ccf, weight}` mappings:
/// ```yaml
/// subclones:
///   - {ccf: 1.0, weight: 0.6}   # clonal / truncal
///   - {ccf: 0.4, weight: 0.3}   # major subclone
///   - {ccf: 0.15, weight: 0.1}  # minor subclone
/// ```
/// `weight` is optional and defaults to `1.0` (equal share). CCF/weight validity is
/// enforced by `SubcloneModel::new`.
fn parse_subclones(v: &Value) -> Result<SubcloneModel, GenCancerReadsError> {
    let seq = v.as_sequence().ok_or_else(|| {
        GenCancerReadsError::ConfigError("subclones must be a list of {ccf, weight} entries".into())
    })?;
    let mut subclones = Vec::with_capacity(seq.len());
    for (i, entry) in seq.iter().enumerate() {
        let map = entry.as_mapping().ok_or_else(|| {
            GenCancerReadsError::ConfigError(format!(
                "subclones[{i}] must be a mapping with a ccf (and optional weight)"
            ))
        })?;
        let ccf_val = map.get(Value::String("ccf".into())).ok_or_else(|| {
            GenCancerReadsError::ConfigError(format!("subclones[{i}] is missing required 'ccf'"))
        })?;
        let ccf = as_f64(ccf_val, &format!("subclones[{i}].ccf"))?;
        let weight = match map.get(Value::String("weight".into())) {
            Some(w) => as_f64(w, &format!("subclones[{i}].weight"))?,
            None => 1.0,
        };
        subclones.push(Subclone { ccf, weight });
    }
    SubcloneModel::new(subclones)
        .map_err(|e| GenCancerReadsError::ConfigError(format!("invalid subclones: {e}")))
}

/// Column-name synonyms accepted in a `subclones_file` cluster table (matched
/// case-insensitively). Kept tool-agnostic: PyClone/PyClone-VI use
/// `cellular_prevalence`; PCAWG-11 / DPClust cluster tables use `ccf` + `n_ssms`.
const CLUSTER_COLS: &[&str] = &["cluster_id", "cluster", "clusterid"];
const CCF_COLS: &[&str] = &["cellular_prevalence", "ccf", "cancer_cell_fraction", "cp"];
const WEIGHT_COLS: &[&str] = &[
    "n_ssms",
    "size",
    "weight",
    "n_mutations",
    "n_snvs",
    "n_variants",
];

/// Build a [`SubcloneModel`] from a deconvolution-tool cluster table (#405, B1).
///
/// Accepts a tab-separated file with a header and folds the two shapes real tools
/// emit into `{ccf, weight}` clusters:
///
/// - **Cluster table** (PCAWG-11 / CSR / DPClust): one row per cluster, a `ccf`
///   column and a size column (`n_ssms` / `size` / `weight`) → used directly.
/// - **Per-mutation table** (PyClone / PyClone-VI): one row per mutation with
///   `cluster_id` + `cellular_prevalence` and no size column → grouped by cluster,
///   `ccf` = mean prevalence, `weight` = mutation count.
///
/// Robustness for real files: CCF > 1.0 is clamped to 1.0 (noisy clonal clusters);
/// rows with non-finite or ≤ 0 CCF are skipped; both are warned once with a count.
/// Extra columns (std, assignment probability, …) are ignored.
fn parse_subclones_file(path: &Path) -> Result<SubcloneModel, GenCancerReadsError> {
    let text = fs::read_to_string(path)
        .map_err(|e| GenCancerReadsError::ConfigError(format!("subclones_file {path:?}: {e}")))?;
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));

    let header = lines.next().ok_or_else(|| {
        GenCancerReadsError::ConfigError(format!("subclones_file {path:?} is empty"))
    })?;
    let cols: Vec<String> = header
        .split('\t')
        .map(|c| c.trim().to_lowercase())
        .collect();
    let find = |names: &[&str]| cols.iter().position(|c| names.contains(&c.as_str()));

    let cluster_idx = find(CLUSTER_COLS).ok_or_else(|| {
        GenCancerReadsError::ConfigError(format!(
            "subclones_file {path:?}: no cluster column (one of {CLUSTER_COLS:?})"
        ))
    })?;
    let ccf_idx = find(CCF_COLS).ok_or_else(|| {
        GenCancerReadsError::ConfigError(format!(
            "subclones_file {path:?}: no CCF column (one of {CCF_COLS:?})"
        ))
    })?;
    let weight_idx = find(WEIGHT_COLS);

    // Aggregate by cluster id, preserving first-seen order for determinism.
    struct Acc {
        ccf_sum: f64,
        rows: usize,
        weight_sum: f64,
    }
    let mut order: Vec<String> = Vec::new();
    let mut accs: HashMap<String, Acc> = HashMap::new();
    let (mut clamped, mut skipped) = (0usize, 0usize);

    for (n, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).map(|s| s.trim());
        let parse_f = |name: &str, s: Option<&str>| -> Result<f64, GenCancerReadsError> {
            s.and_then(|s| s.parse::<f64>().ok()).ok_or_else(|| {
                GenCancerReadsError::ConfigError(format!(
                    "subclones_file {path:?}: row {} has an unparseable {name}",
                    n + 2 // +1 header, +1 to 1-index
                ))
            })
        };
        let cluster = match get(cluster_idx) {
            Some(c) if !c.is_empty() => c.to_string(),
            _ => continue,
        };
        let mut ccf = parse_f("ccf", get(ccf_idx))?;
        if !ccf.is_finite() || ccf <= 0.0 {
            skipped += 1;
            continue;
        }
        if ccf > 1.0 {
            ccf = 1.0;
            clamped += 1;
        }
        let weight = match weight_idx {
            Some(i) => parse_f("weight", get(i))?,
            None => 1.0,
        };
        let acc = accs.entry(cluster.clone()).or_insert_with(|| {
            order.push(cluster.clone());
            Acc {
                ccf_sum: 0.0,
                rows: 0,
                weight_sum: 0.0,
            }
        });
        acc.ccf_sum += ccf;
        acc.rows += 1;
        acc.weight_sum += weight;
    }

    if clamped > 0 {
        warn!("subclones_file {path:?}: clamped {clamped} CCF value(s) > 1.0 to 1.0");
    }
    if skipped > 0 {
        warn!("subclones_file {path:?}: skipped {skipped} row(s) with non-positive/invalid CCF");
    }

    let subclones: Vec<Subclone> = order
        .iter()
        .map(|k| {
            let a = &accs[k];
            Subclone {
                ccf: a.ccf_sum / a.rows as f64,
                // No size column → this is a per-mutation table; weight = row count.
                weight: if weight_idx.is_some() {
                    a.weight_sum
                } else {
                    a.rows as f64
                },
            }
        })
        .collect();

    SubcloneModel::new(subclones)
        .map_err(|e| GenCancerReadsError::ConfigError(format!("subclones_file {path:?}: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h1n1() -> PathBuf {
        PathBuf::from(format!(
            "{}/test_data/references/H1N1.fa",
            env!("CARGO_MANIFEST_DIR")
        ))
    }

    fn base_scrape() -> HashMap<String, Value> {
        let mut s = HashMap::new();
        s.insert(
            "reference".into(),
            Value::String(h1n1().to_string_lossy().into()),
        );
        s.insert("output_dir".into(), Value::String("/tmp".into()));
        s
    }

    #[test]
    fn default_tumor_rate_is_1e_5() {
        let cfg = CancerConfig::from_scrape(base_scrape()).unwrap();
        assert_eq!(cfg.tumor_mutation_rate, Some(1e-5));
    }

    #[test]
    fn tumor_rate_model_sentinel_defers_to_model() {
        let mut s = base_scrape();
        s.insert("tumor_mutation_rate".into(), Value::String("model".into()));
        let cfg = CancerConfig::from_scrape(s).unwrap();
        assert_eq!(cfg.tumor_mutation_rate, None);
    }

    #[test]
    fn per_pass_coverage_splits_by_purity() {
        let mut s = base_scrape();
        s.insert("total_coverage".into(), Value::Number(30.into()));
        s.insert("purity".into(), Value::from(0.7));
        let cfg = CancerConfig::from_scrape(s).unwrap();
        assert_eq!(cfg.per_pass_coverage(), (9, 21)); // (1-0.7)*30=9, 0.7*30=21
    }

    #[test]
    fn rejects_purity_endpoints() {
        let mut s = base_scrape();
        s.insert("purity".into(), Value::Number(1.into()));
        assert!(matches!(
            CancerConfig::from_scrape(s),
            Err(GenCancerReadsError::PurityOutOfRange(_))
        ));
    }

    #[test]
    fn rejects_zero_per_pass_coverage() {
        let mut s = base_scrape();
        s.insert("total_coverage".into(), Value::Number(1.into()));
        s.insert("purity".into(), Value::from(0.5));
        // 1*0.5 rounds to 0 for the... 0.5 rounds to 1 actually; use extreme purity
        s.insert("purity".into(), Value::from(0.1));
        // (1-0.1)*1=0.9->1 normal, 0.1*1=0.1->0 tumor
        assert!(matches!(
            CancerConfig::from_scrape(s),
            Err(GenCancerReadsError::PerPassCoverageZero { .. })
        ));
    }

    #[test]
    fn from_scrape_creates_missing_output_dir() {
        let tmp = tempfile::tempdir().unwrap();
        let nested = tmp.path().join("a/b/cancer_out");
        assert!(!nested.exists());
        let mut s = base_scrape();
        s.insert(
            "output_dir".into(),
            Value::String(nested.to_string_lossy().into()),
        );
        let cfg = CancerConfig::from_scrape(s).unwrap();
        assert!(
            nested.is_dir(),
            "output_dir should be created (create_dir_all)"
        );
        assert_eq!(cfg.output_dir, nested);
    }

    #[test]
    fn paired_ended_accepts_fragment_model_alone() {
        // #355 gave gen-reads this; gen-cancer-reads now matches. A fragment_model
        // satisfies the paired-end fragment-source requirement with no mean/st_dev,
        // and must flow through to BOTH the normal and tumor passes.
        let mut s = base_scrape();
        s.insert("paired_ended".into(), Value::Bool(true));
        s.insert(
            "fragment_model".into(),
            Value::String("/tmp/frag.json.gz".into()),
        );
        let cfg =
            CancerConfig::from_scrape(s).expect("fragment_model alone should satisfy paired_ended");
        let want = Some(PathBuf::from("/tmp/frag.json.gz"));
        assert_eq!(cfg.fragment_model, want);
        assert_eq!(cfg.normal_pass().unwrap().fragment_model, want);
        assert_eq!(
            cfg.tumor_pass(PathBuf::from("/tmp/g.vcf.gz"))
                .unwrap()
                .fragment_model,
            want
        );
    }

    #[test]
    fn paired_ended_requires_a_fragment_source() {
        // Neither a fragment_model nor mean/st_dev → still rejected.
        let mut s = base_scrape();
        s.insert("paired_ended".into(), Value::Bool(true));
        assert!(matches!(
            CancerConfig::from_scrape(s),
            Err(GenCancerReadsError::ConfigError(_))
        ));
    }

    #[test]
    fn normal_pass_has_no_svs_tumor_has_scale() {
        let mut s = base_scrape();
        s.insert("sv_rate_scale".into(), Value::from(5.0));
        let cfg = CancerConfig::from_scrape(s).unwrap();
        let normal = cfg.normal_pass().unwrap();
        let tumor = cfg.tumor_pass(PathBuf::from("/tmp/g.vcf.gz")).unwrap();
        assert_eq!(normal.sv_rate_scale, 0.0);
        assert_eq!(tumor.sv_rate_scale, 5.0);
        assert!(normal.output_filename.ends_with("_normal"));
        assert!(tumor.output_filename.ends_with("_tumor"));
        assert_eq!(tumor.input_vcf, Some(PathBuf::from("/tmp/g.vcf.gz")));
    }

    #[test]
    fn no_subclones_by_default() {
        let cfg = CancerConfig::from_scrape(base_scrape()).unwrap();
        assert!(cfg.subclones.is_none());
        // And it must not leak into either pass.
        assert!(cfg.normal_pass().unwrap().subclone_model.is_none());
        assert!(
            cfg.tumor_pass(PathBuf::from("/tmp/g.vcf.gz"))
                .unwrap()
                .subclone_model
                .is_none()
        );
    }

    #[test]
    fn subclones_parse_and_land_on_tumor_pass_only() {
        let mut s = base_scrape();
        let sub: Value =
            serde_yml::from_str("- {ccf: 1.0, weight: 3.0}\n- {ccf: 0.2, weight: 1.0}").unwrap();
        s.insert("subclones".into(), sub);
        let cfg = CancerConfig::from_scrape(s).unwrap();

        let model = cfg.subclones.as_ref().expect("subclones parsed");
        assert_eq!(model.subclones().len(), 2);
        assert_eq!(
            model.subclones()[0],
            Subclone {
                ccf: 1.0,
                weight: 3.0
            }
        );
        assert_eq!(
            model.subclones()[1],
            Subclone {
                ccf: 0.2,
                weight: 1.0
            }
        );

        // The normal pass has no somatic variants → no model; the tumor pass carries it.
        assert!(cfg.normal_pass().unwrap().subclone_model.is_none());
        assert_eq!(
            cfg.tumor_pass(PathBuf::from("/tmp/g.vcf.gz"))
                .unwrap()
                .subclone_model,
            Some(model.clone())
        );
    }

    #[test]
    fn subclones_weight_defaults_to_one() {
        let mut s = base_scrape();
        let sub: Value = serde_yml::from_str("- {ccf: 0.5}\n- {ccf: 0.3}").unwrap();
        s.insert("subclones".into(), sub);
        let cfg = CancerConfig::from_scrape(s).unwrap();
        let model = cfg.subclones.unwrap();
        assert_eq!(
            model.subclones()[0],
            Subclone {
                ccf: 0.5,
                weight: 1.0
            }
        );
        assert_eq!(
            model.subclones()[1],
            Subclone {
                ccf: 0.3,
                weight: 1.0
            }
        );
    }

    #[test]
    fn subclones_reject_bad_ccf() {
        let mut s = base_scrape();
        let sub: Value = serde_yml::from_str("- {ccf: 1.5, weight: 1.0}").unwrap();
        s.insert("subclones".into(), sub);
        assert!(matches!(
            CancerConfig::from_scrape(s),
            Err(GenCancerReadsError::ConfigError(_))
        ));
    }

    #[test]
    fn subclones_reject_empty_list() {
        let mut s = base_scrape();
        s.insert("subclones".into(), Value::Sequence(vec![]));
        assert!(matches!(
            CancerConfig::from_scrape(s),
            Err(GenCancerReadsError::ConfigError(_))
        ));
    }

    // ── subclones_file (B1: ingest a deconvolution cluster table) ────────────

    fn write_tmp(name: &str, body: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("eidolon_sctest_{name}.tsv"));
        fs::write(&p, body).unwrap();
        p
    }

    #[test]
    fn subclones_file_pyclone_per_mutation_groups_by_cluster() {
        // PyClone-VI shape: one row per mutation, no size column → weight = count,
        // ccf = mean cellular_prevalence. Extra columns are ignored.
        let tsv = write_tmp(
            "pyclone",
            "mutation_id\tsample_id\tcluster_id\tcellular_prevalence\tcluster_assignment_prob\n\
             m1\tS\t0\t1.0\t0.99\n\
             m2\tS\t0\t1.0\t0.98\n\
             m3\tS\t0\t1.0\t0.97\n\
             m4\tS\t1\t0.4\t0.95\n\
             m5\tS\t1\t0.4\t0.90\n",
        );
        let m = parse_subclones_file(&tsv).unwrap();
        assert_eq!(m.subclones().len(), 2);
        assert_eq!(
            m.subclones()[0],
            Subclone {
                ccf: 1.0,
                weight: 3.0
            }
        );
        assert_eq!(
            m.subclones()[1],
            Subclone {
                ccf: 0.4,
                weight: 2.0
            }
        );
    }

    #[test]
    fn subclones_file_pcawg_cluster_table_uses_size_as_weight() {
        // PCAWG-11 / DPClust shape: one row per cluster with an n_ssms size column.
        let tsv = write_tmp(
            "pcawg",
            "cluster\tn_ssms\tccf\n\
             1\t1200\t1.0\n\
             2\t300\t0.55\n\
             3\t80\t0.2\n",
        );
        let m = parse_subclones_file(&tsv).unwrap();
        assert_eq!(m.subclones().len(), 3);
        assert_eq!(
            m.subclones()[0],
            Subclone {
                ccf: 1.0,
                weight: 1200.0
            }
        );
        assert_eq!(
            m.subclones()[2],
            Subclone {
                ccf: 0.2,
                weight: 80.0
            }
        );
    }

    #[test]
    fn subclones_file_clamps_ccf_above_one_and_skips_nonpositive() {
        let tsv = write_tmp(
            "clamp",
            "cluster_id\tcellular_prevalence\n\
             0\t1.03\n\
             0\t0.97\n\
             1\t0.0\n\
             2\t0.5\n",
        );
        let m = parse_subclones_file(&tsv).unwrap();
        // cluster 0: 1.03 clamped to 1.0, mean(1.0, 0.97) = 0.985; cluster 1 skipped
        // (ccf 0.0); cluster 2 kept.
        assert_eq!(m.subclones().len(), 2);
        assert!((m.subclones()[0].ccf - 0.985).abs() < 1e-9);
        assert_eq!(m.subclones()[0].weight, 2.0);
        assert_eq!(
            m.subclones()[1],
            Subclone {
                ccf: 0.5,
                weight: 1.0
            }
        );
    }

    #[test]
    fn subclones_file_missing_ccf_column_errors() {
        let tsv = write_tmp("nocol", "cluster_id\tfoo\n0\t1.0\n");
        assert!(matches!(
            parse_subclones_file(&tsv),
            Err(GenCancerReadsError::ConfigError(_))
        ));
    }

    #[test]
    fn subclones_file_lands_on_tumor_pass() {
        let tsv = write_tmp("wire", "cluster\tn_ssms\tccf\n1\t10\t1.0\n2\t5\t0.3\n");
        let mut s = base_scrape();
        s.insert(
            "subclones_file".into(),
            Value::String(tsv.to_string_lossy().into()),
        );
        let cfg = CancerConfig::from_scrape(s).unwrap();
        assert_eq!(cfg.subclones.as_ref().unwrap().subclones().len(), 2);
        assert!(cfg.normal_pass().unwrap().subclone_model.is_none());
        assert!(
            cfg.tumor_pass(PathBuf::from("/tmp/g.vcf.gz"))
                .unwrap()
                .subclone_model
                .is_some()
        );
    }

    #[test]
    fn somatic_vcf_lands_on_tumor_pass_with_purity_scale() {
        let mut s = base_scrape();
        s.insert("purity".into(), Value::from(0.8));
        s.insert("somatic_vcf".into(), Value::String("/tmp/som.vcf".into()));
        let cfg = CancerConfig::from_scrape(s).unwrap();
        assert_eq!(cfg.somatic_vcf, Some(PathBuf::from("/tmp/som.vcf")));

        // Only the tumor pass replays it, at 1/purity.
        let normal = cfg.normal_pass().unwrap();
        let tumor = cfg.tumor_pass(PathBuf::from("/tmp/g.vcf.gz")).unwrap();
        assert!(normal.somatic_vcf.is_none());
        assert_eq!(tumor.somatic_vcf, Some(PathBuf::from("/tmp/som.vcf")));
        assert!((tumor.somatic_af_scale - 1.25).abs() < 1e-9); // 1/0.8
    }

    #[test]
    fn subclones_inline_and_file_are_mutually_exclusive() {
        let tsv = write_tmp("mutex", "cluster\tccf\n1\t1.0\n");
        let mut s = base_scrape();
        s.insert(
            "subclones".into(),
            serde_yml::from_str("- {ccf: 1.0}").unwrap(),
        );
        s.insert(
            "subclones_file".into(),
            Value::String(tsv.to_string_lossy().into()),
        );
        assert!(matches!(
            CancerConfig::from_scrape(s),
            Err(GenCancerReadsError::ConfigError(_))
        ));
    }
}
