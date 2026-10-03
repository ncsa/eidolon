// Configuration for the `compare-vcfs` subcommand. Stores validated paths and
// flag values only; the VCFs themselves are read inside `runner` so this
// struct stays cheap to construct/test.
use crate::compare_vcfs::errors::CompareVcfsError;
use crate::config_keys::check_keys;
use serde_yml::Value;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

/// Every top-level key the config parser reads. Anything else is rejected (#496).
pub const KNOWN_KEYS: &[&str] = &[
    "called_vcf",
    "chrom_aliases",
    "contigs_simulated",
    "equivalence_window",
    "fast",
    "golden_vcf",
    "include_filtered",
    "include_homs",
    "mutation_bed",
    "output_dir",
    "overwrite_output",
    "reference",
    "target_bed",
    "write_fp_vcf",
];

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    pub golden_vcf: PathBuf,
    pub called_vcf: PathBuf,
    pub reference: PathBuf,
    pub output_dir: PathBuf,
    pub overwrite_output: bool,
    pub target_bed: Option<PathBuf>,
    /// `None` means "treat every contig in the golden VCF as simulated".
    pub contigs_simulated: Option<Vec<String>>,
    pub include_homs: bool,
    pub include_filtered: bool,
    /// Half-width of the ±N bp window used for equivalence detection.
    /// Default 50 (matches NEAT 2.1's `EV_BPRANGE`).
    pub equivalence_window: usize,
    /// Skip the equivalence sweep entirely. Matches NEAT 2.1's `--fast`.
    pub fast: bool,
    /// Mutation regions BED used by the simulator. Used only for FN
    /// attribution — FNs whose position falls outside these regions are
    /// tagged `outside_mutation_bed`. None when not configured.
    pub mutation_bed: Option<PathBuf>,
    /// Two-column TSV remapping BED chrom names to the reference's
    /// canonical names (e.g. `1\tchr1`). Applied to both `mutation_bed`
    /// and `target_bed` at load time. None when not configured.
    pub chrom_aliases: Option<PathBuf>,
    /// Also write `<output_dir>/FP.vcf` listing every false-positive call.
    /// Off by default — FP sets can be large on noisy callers and the
    /// `comparison_summary.json`/`.txt` already report counts.
    pub write_fp_vcf: bool,
}

impl RunConfiguration {
    pub fn from(yml_file: &PathBuf) -> Result<Self, CompareVcfsError> {
        let file = fs::File::open(yml_file).map_err(|e| {
            CompareVcfsError::ConfigurationError(format!(
                "could not open config {}: {e}",
                yml_file.display()
            ))
        })?;
        let scrape: HashMap<String, Value> = serde_yml::from_reader(file)?;
        check_keys(
            scrape.keys().map(String::as_str),
            KNOWN_KEYS,
            "compare-vcfs",
        )
        .map_err(CompareVcfsError::ConfigurationError)?;

        let golden_vcf = require_path(&scrape, "golden_vcf")?;
        reject_bcf(&golden_vcf)?;
        let called_vcf = require_path(&scrape, "called_vcf")?;
        reject_bcf(&called_vcf)?;
        let reference = require_path(&scrape, "reference")?;

        let output_dir = require_str(&scrape, "output_dir").map(PathBuf::from)?;
        if output_dir.is_file() {
            return Err(CompareVcfsError::ConfigurationError(format!(
                "output_dir {} exists as a file, not a directory",
                output_dir.display()
            )));
        }

        let overwrite_output = scrape
            .get("overwrite_output")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let target_bed = optional_path(&scrape, "target_bed")?;
        let contigs_simulated = optional_string_list(&scrape, "contigs_simulated");
        let include_homs = scrape
            .get("include_homs")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);
        let include_filtered = scrape
            .get("include_filtered")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);
        let equivalence_window = scrape
            .get("equivalence_window")
            .and_then(|v| v.as_u64())
            .unwrap_or(50) as usize;
        let fast = scrape
            .get("fast")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);
        let mutation_bed = optional_path(&scrape, "mutation_bed")?;
        let chrom_aliases = optional_path(&scrape, "chrom_aliases")?;
        let write_fp_vcf = scrape
            .get("write_fp_vcf")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        Ok(RunConfiguration {
            golden_vcf,
            called_vcf,
            reference,
            output_dir,
            overwrite_output,
            target_bed,
            contigs_simulated,
            include_homs,
            include_filtered,
            equivalence_window,
            fast,
            mutation_bed,
            chrom_aliases,
            write_fp_vcf,
        })
    }
}

fn require_str(scrape: &HashMap<String, Value>, key: &str) -> Result<String, CompareVcfsError> {
    scrape
        .get(key)
        .and_then(|v| v.as_str())
        .map(str::to_string)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| CompareVcfsError::ConfigurationError(format!("{key} is required in config")))
}

fn require_path(scrape: &HashMap<String, Value>, key: &str) -> Result<PathBuf, CompareVcfsError> {
    let raw = require_str(scrape, key)?;
    let p = PathBuf::from(raw);
    if !p.is_file() {
        return Err(CompareVcfsError::ConfigurationError(format!(
            "{key} not found: {}",
            p.display()
        )));
    }
    Ok(p)
}

fn optional_path(
    scrape: &HashMap<String, Value>,
    key: &str,
) -> Result<Option<PathBuf>, CompareVcfsError> {
    let raw = scrape.get(key).and_then(|v| v.as_str()).unwrap_or(".");
    if raw == "." || raw.is_empty() {
        return Ok(None);
    }
    let p = PathBuf::from(raw);
    if !p.is_file() {
        return Err(CompareVcfsError::ConfigurationError(format!(
            "{key} not found: {}",
            p.display()
        )));
    }
    Ok(Some(p))
}

fn optional_string_list(scrape: &HashMap<String, Value>, key: &str) -> Option<Vec<String>> {
    let raw = scrape.get(key).and_then(|v| v.as_str()).unwrap_or(".");
    if raw == "." || raw.is_empty() {
        return None;
    }
    let items: Vec<String> = raw
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if items.is_empty() { None } else { Some(items) }
}

fn reject_bcf(p: &Path) -> Result<(), CompareVcfsError> {
    if p.extension().and_then(|s| s.to_str()) == Some("bcf") {
        return Err(CompareVcfsError::ConfigurationError(format!(
            "{} is a BCF; only .vcf and .vcf.gz are supported",
            p.display()
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn manifest_path(rel: &str) -> String {
        format!("{}/{rel}", env!("CARGO_MANIFEST_DIR"))
    }

    fn write_yaml(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".yml").tempfile().unwrap();
        write!(f, "{content}").unwrap();
        f
    }

    #[test]
    fn test_minimal_config_parses() {
        let h1n1 = manifest_path("test_data/references/H1N1.fa");
        let vcf = manifest_path("test_data/vcfs/small_snps.vcf");
        let tmp_dir = tempfile::tempdir().unwrap();
        let yml = write_yaml(&format!(
            "golden_vcf: {vcf}\ncalled_vcf: {vcf}\nreference: {h1n1}\noutput_dir: {}\n",
            tmp_dir.path().display()
        ));
        let cfg = RunConfiguration::from(&yml.path().to_path_buf()).unwrap();
        assert!(cfg.target_bed.is_none());
        assert!(cfg.contigs_simulated.is_none());
        assert!(!cfg.include_homs);
        assert!(!cfg.include_filtered);
        assert!(!cfg.overwrite_output);
    }

    #[test]
    fn test_missing_required_field_errors() {
        let yml = write_yaml("golden_vcf: /tmp/nope.vcf\n");
        let err = RunConfiguration::from(&yml.path().to_path_buf()).unwrap_err();
        assert!(matches!(err, CompareVcfsError::ConfigurationError(_)));
    }

    #[test]
    fn test_contigs_simulated_parses_csv() {
        let h1n1 = manifest_path("test_data/references/H1N1.fa");
        let vcf = manifest_path("test_data/vcfs/small_snps.vcf");
        let tmp_dir = tempfile::tempdir().unwrap();
        let yml = write_yaml(&format!(
            "golden_vcf: {vcf}\ncalled_vcf: {vcf}\nreference: {h1n1}\noutput_dir: {}\n\
             contigs_simulated: chr1, chr2,  chrX\n",
            tmp_dir.path().display()
        ));
        let cfg = RunConfiguration::from(&yml.path().to_path_buf()).unwrap();
        assert_eq!(
            cfg.contigs_simulated.unwrap(),
            vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()]
        );
    }

    #[test]
    fn test_bcf_rejected() {
        let h1n1 = manifest_path("test_data/references/H1N1.fa");
        let bcf = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
        let tmp_dir = tempfile::tempdir().unwrap();
        let yml = write_yaml(&format!(
            "golden_vcf: {}\ncalled_vcf: {}\nreference: {h1n1}\noutput_dir: {}\n",
            bcf.path().display(),
            bcf.path().display(),
            tmp_dir.path().display()
        ));
        let err = RunConfiguration::from(&yml.path().to_path_buf()).unwrap_err();
        match err {
            CompareVcfsError::ConfigurationError(msg) => assert!(msg.contains("BCF")),
            _ => panic!("expected ConfigurationError, got {err:?}"),
        }
    }
}
