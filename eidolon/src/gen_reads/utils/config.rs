// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.
use chrono::Utc;

use crate::eidolon_core::file_tools::folder_tools::check_create_dir;
use crate::gen_reads::errors::GenerateReadsError;
use crate::gen_reads::utils::subclone;
use log::{error, info, warn};
use serde_yml::Value;
use std::collections::HashMap;
use std::path::PathBuf;
use std::string::String;
use std::{env, fs};

/// Optional 3' sequencing-adapter readthrough config (#125). When `enabled`,
/// reads whose insert is shorter than `read_len` are padded at the 3' end with
/// adapter sequence (R1 gets `r1`, R2 gets `r2`), matching real Illumina data so
/// downstream adapter-trim QC stages are actually exercised. Default disabled.
/// `r1`/`r2` are resolved 5'->3' uppercase-ACGT sequences (from a preset or custom).
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    pub enabled: bool,
    pub r1: String,
    pub r2: String,
}

impl Default for AdapterConfig {
    fn default() -> Self {
        AdapterConfig {
            enabled: false,
            r1: String::new(),
            r2: String::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RunConfiguration {
    // This struct holds all the parameters for this particular run. It is derived from input either
    // from a configuration file or from command line inputs. This is not built directly in the code,
    // but is constructed by a builder to enable default values.
    //
    // reference: The path to the reference for the run.
    // read_len: The length of reads in the output fastq.
    // coverage: The average depth of coverage for the output fastq file.
    // mutation_rate: The rate of mutation for the file.
    // ploidy: The number of copies of each chromosome in the target organism. Mutation process will
    // be replicated this number of times.
    // paired_ended: If the run will be in paired-ended mode.
    // fragment_mean: Mean size of the fragments in paired-ended mode.
    // fragment_st_dev: Standard deviation of the fragment mean describing the sample set to sample
    // fragments from.
    // produce_fastq: True or false on whether to produce an output fastq file.
    // produce_vcf: True or false on whether to produce an output VCF file, with genotyped variants.
    // produce_bam: True or false on whether to produce an output BAM file, which will be aligned to
    // the reference.
    // overwrite_output: if true, will overwrite output. If false will error and exit you attempt to
    // overwrite files with the same name.
    // output_dir: The directory, relative or absolute, path to the directory to place output.
    // output_prefix: The name to use for the output files.
    pub reference: PathBuf,
    pub read_len: usize,
    pub coverage: usize,
    pub mutation_rate: Option<f64>,
    pub ploidy: usize,
    pub paired_ended: bool,
    pub fragment_mean: Option<f64>,
    pub fragment_st_dev: Option<f64>,
    pub produce_fastq: bool,
    pub produce_vcf: bool,
    pub produce_bam: bool,
    pub rng_seed: Option<String>,
    pub seed_vec: Vec<String>,
    pub overwrite_output: bool,
    pub minimum_mutations: usize,
    pub output_dir: PathBuf,
    pub output_filename: String,
    pub output_fastq_1: Option<PathBuf>,
    pub output_fastq_2: Option<PathBuf>,
    pub output_vcf: Option<PathBuf>,
    pub output_bam: Option<PathBuf>,
    // model input
    pub quality_score_model: Option<PathBuf>,
    pub mutation_model: Option<PathBuf>,
    pub fragment_model: Option<PathBuf>,
    pub sequence_error_model: Option<PathBuf>,
    // optional BED filter applied during generation (not post-processing)
    pub target_bed: Option<PathBuf>,
    // Option BED-style file with columns: contig, start, end, other
    // where "other" contains the key "mut_rate=X.XXXX" where X.XXXX represents any arbitrary
    // float. That will become the mutation rate for that particular region
    pub mutation_regions: Option<PathBuf>,
    // optional VCF of variants to force into the simulation
    pub input_vcf: Option<PathBuf>,
    // optional gc-bias model
    pub(crate) gc_bias_model: Option<PathBuf>,
    // normalize coverage true or false
    pub(crate) gc_bias_normalize_coverage: bool,
    // maximum threads for parallel contig processing (None = rayon default = all cores)
    pub num_threads: Option<usize>,
    // Target size (bp) of a parallel sub-contig chunk. None (default) or Some(0)
    // = chunking disabled — one chunk spans the whole contig (see
    // resolve_chunk_size, which maps both to usize::MAX); Some(n) = fixed n bp.
    // Independent of num_threads, so output is identical across thread counts.
    pub chunk_size: Option<usize>,
    // when true, fragments shorter than read_len are kept and produce truncated reads
    // (long-read platforms); when false, such fragments are discarded (short-read default)
    pub long_reads: bool,
    // when true, short-read PE keeps fragments whose insert < read_len and emits
    // insert-length GENOMIC reads (no adapter). This is the short-insert library WITHOUT
    // modeling 3' adapter readthrough — the isolation control for adapter validation
    // (#125). adapters.enabled implies this (readthrough needs short fragments kept);
    // set it on its own to get the adapter-free short-insert arm. Default off →
    // short fragments discarded, output byte-identical to before.
    pub keep_short_fragments: bool,
    // Scale applied to the SvModel's per-base SV rate at generation time.
    // 0.0 (default) disables de novo SV generation entirely — only SVs
    // from `input_vcf` flow through. Values > 0 enable sampling from
    // `MutationModel.sv_model` (if present); 1.0 reproduces the trained
    // corpus rate, larger values stress-test SV handling.
    pub sv_rate_scale: f64,
    // Max de novo SV length as a fraction of contig length (#229). Default 0.25
    // (contig_len/4) keeps overlap rejection tractable on chromosome-scale
    // references; raise toward 1.0 for small contigs (bacteria/viruses) or
    // native aneuploidy. Only affects de novo SV sampling, not input_vcf SVs.
    pub sv_max_length_fraction: f64,
    // Optional 3' sequencing-adapter readthrough (#125). Default disabled; when
    // disabled, output is byte-identical to pre-adapter behavior.
    pub adapters: AdapterConfig,
    // Optional subclonal architecture for de-novo variants (#405). When `Some`,
    // each de-novo variant is assigned a subclone whose cancer-cell fraction (CCF)
    // *composes* with the variant's dosage into `allele_fraction = dosage × CCF`,
    // yielding a realistic somatic VAF spectrum. Set only on the cancer *tumor*
    // pass; `None` everywhere else keeps output byte-identical to pre-#405 behavior.
    pub subclone_model: Option<subclone::SubcloneModel>,
    // Reproductive somatic input (#405): a VCF of somatic variants to replay in the
    // tumor pass. Loaded as `Provenance::SomaticVcf` (→ merged truth origin
    // `somatic`, not `shared`) with their `allele_fraction` honored. Set only on the
    // cancer tumor pass; `None` elsewhere.
    pub somatic_vcf: Option<PathBuf>,
    // Multiplier applied to `somatic_vcf` variants' `allele_fraction` on load. The
    // cancer tumor pass sets `1/purity` so a supplied *observed* VAF reproduces after
    // tumor/normal mixing (merged VAF = purity × af = the input VAF). Clamped to 1.0.
    // Default 1.0 (no scaling).
    pub somatic_af_scale: f64,
}

impl Default for RunConfiguration {
    fn default() -> Self {
        RunConfiguration {
            reference: PathBuf::new(),
            read_len: 151,
            coverage: 10,
            mutation_rate: None,
            ploidy: 2,
            paired_ended: false,
            fragment_mean: None,
            fragment_st_dev: None,
            produce_fastq: true,
            produce_vcf: false,
            produce_bam: false,
            rng_seed: None,
            seed_vec: Vec::new(),
            overwrite_output: false,
            minimum_mutations: 0,
            output_dir: env::current_dir().unwrap_or_default(),
            output_filename: "neat_out".to_string(),
            output_fastq_1: None,
            output_fastq_2: None,
            output_vcf: None,
            output_bam: None,
            quality_score_model: None,
            mutation_model: None,
            fragment_model: None,
            sequence_error_model: None,
            target_bed: None,
            mutation_regions: None,
            input_vcf: None,
            gc_bias_model: None,
            gc_bias_normalize_coverage: true,
            num_threads: None,
            chunk_size: None,
            long_reads: false,
            keep_short_fragments: false,
            sv_rate_scale: 0.0,
            sv_max_length_fraction: 0.25,
            adapters: AdapterConfig::default(),
            subclone_model: None,
            somatic_vcf: None,
            somatic_af_scale: 1.0,
        }
    }
}

impl RunConfiguration {
    pub fn from_yaml_file(yaml_file: &PathBuf) -> Result<RunConfiguration, GenerateReadsError> {
        // Reads an input configuration file from yaml using the serde package. Then sets the
        // parameters based on the inputs. A "." value means to use the default value.

        // Opens file for reading
        let f = fs::File::open(yaml_file);
        let file = match f {
            Ok(l) => l,
            Err(error) => panic!(
                "Problem reading the config file: {:?}. System error: {}",
                &yaml_file, error,
            ),
        };
        // Uses serde_yml to read the file into a HashMap
        let scrape_config: HashMap<String, Value> = serde_yml::from_reader(file)?;
        Self::from_scrape_config(scrape_config)
    }

    pub fn from_scrape_config(
        scrape_config: HashMap<String, Value>,
    ) -> Result<RunConfiguration, GenerateReadsError> {
        // Fill in the reference first, all hinges on that
        let reference = scrape_config
            .get("reference")
            .ok_or(GenerateReadsError::MissingReferenceError)?;

        let ref_str = reference.as_str();
        let ref_buf = match ref_str {
            Some(str) => {
                let path = PathBuf::from(&str);
                if !path.is_file() {
                    return Err(GenerateReadsError::FileNotFound(str.to_owned()));
                }
                path
            }
            None => return Err(GenerateReadsError::MissingReferenceError),
        };

        let mut config = RunConfiguration {
            reference: ref_buf,
            ..Default::default()
        };

        info!(
            "Running eidolon to generate reads on {:?} with...",
            &config.reference
        );

        for (key, value) in scrape_config {
            match &value.as_str() {
                // Any item with a . for a value we keep the default
                Some(".") => continue,
                _ => match key.as_str() {
                    "read_len" => {
                        config.read_len = value.as_u64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "read_len".to_string(),
                                "integer".to_string(),
                            )
                        })? as usize;
                    }
                    "coverage" => {
                        config.coverage = value.as_u64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "coverage".to_string(),
                                "integer".to_string(),
                            )
                        })? as usize;
                    }
                    "mutation_rate" => {
                        config.mutation_rate = Some(value.as_f64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "mutation_rate".to_string(),
                                "float".to_string(),
                            )
                        })?);
                    }
                    "ploidy" => {
                        config.ploidy = value.as_u64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "ploidy".to_string(),
                                "integer".to_string(),
                            )
                        })? as usize;
                    }
                    "paired_ended" => {
                        config.paired_ended = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "paired_ended".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "fragment_mean" => {
                        config.fragment_mean = Some(value.as_f64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "fragment_mean".to_string(),
                                "float".to_string(),
                            )
                        })?);
                    }
                    "fragment_st_dev" => {
                        config.fragment_st_dev = Some(value.as_f64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "frament_st_dev".to_string(),
                                "float".to_string(),
                            )
                        })?);
                    }
                    "produce_fastq" => {
                        config.produce_fastq = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "produce_fastq".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "produce_vcf" => {
                        config.produce_vcf = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "produce_vcf".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "produce_bam" => {
                        config.produce_bam = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "produce_bam".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "rng_seed" => {
                        // Accept either a quoted string (whitespace-separated multi-term seeds)
                        // or a bare integer literal; integers are coerced to their string form
                        // since the downstream seeder consumes whitespace-split string tokens.
                        let seed_string = if let Some(s) = value.as_str() {
                            s.to_string()
                        } else if let Some(n) = value.as_u64() {
                            n.to_string()
                        } else if let Some(n) = value.as_i64() {
                            n.to_string()
                        } else {
                            return Err(GenerateReadsError::ConfigReadError(
                                "rng_seed".to_string(),
                                "String or integer".to_string(),
                            ));
                        };
                        config.rng_seed = Some(seed_string);
                    }
                    "overwrite_output" => {
                        config.overwrite_output = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "overwrite_output".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "minimum_mutations" => {
                        config.minimum_mutations = value.as_u64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "minimum_mutations".to_string(),
                                "Integer".to_string(),
                            )
                        })? as usize;
                    }
                    "output_dir" => {
                        let output_path = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "output_dir".to_string(),
                                "String".to_string(),
                            )
                        })?;
                        config.output_dir = PathBuf::from(output_path);
                        if !config.output_dir.is_dir() {
                            info!("Creating output directory: {:?}", config.output_dir);
                            check_create_dir(&config.output_dir);
                        }
                    }
                    // output_prefix is for backward compatability
                    "output_filename" | "output_prefix" => {
                        config.output_filename = value
                            .as_str()
                            .ok_or_else(|| {
                                GenerateReadsError::ConfigReadError(
                                    "output_filename".to_string(),
                                    "String".to_string(),
                                )
                            })?
                            .to_string();
                    }
                    "quality_score_model" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "quality_score_model".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.quality_score_model = Some(filename);
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "mutation_model" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "mutation_model".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.mutation_model = Some(filename);
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "fragment_model" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "fragment_model".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.fragment_model = Some(filename);
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "sequence_error_model" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "sequence_error_model".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.sequence_error_model = Some(filename);
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "target_bed" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "target_bed".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.target_bed = Some(filename);
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "mutation_regions" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "mutation_regions".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.mutation_regions = Some(filename)
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "input_vcf" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "input_vcf".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.input_vcf = Some(filename);
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "gc_bias_model" => {
                        let name = value.as_str().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "gc_bias_model".to_string(),
                                "path to file".to_string(),
                            )
                        })?;
                        let filename = PathBuf::from(name);
                        if filename.is_file() {
                            config.gc_bias_model = Some(filename)
                        } else {
                            return Err(GenerateReadsError::FileNotFound(name.to_string()));
                        }
                    }
                    "gc_bias_normalize_coverage" => {
                        config.gc_bias_normalize_coverage = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "gc_bias_normalize_coverage".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "num_threads" => {
                        config.num_threads = Some(value.as_u64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "num_threads".to_string(),
                                "integer".to_string(),
                            )
                        })? as usize);
                    }
                    "chunk_size" => {
                        config.chunk_size = Some(value.as_u64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "chunk_size".to_string(),
                                "integer".to_string(),
                            )
                        })? as usize);
                    }
                    "long_reads" => {
                        config.long_reads = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "long_reads".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "keep_short_fragments" => {
                        config.keep_short_fragments = value.as_bool().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "keep_short_fragments".to_string(),
                                "boolean".to_string(),
                            )
                        })?;
                    }
                    "sv_rate_scale" => {
                        config.sv_rate_scale = value.as_f64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "sv_rate_scale".to_string(),
                                "float".to_string(),
                            )
                        })?;
                    }
                    "sv_max_length_fraction" => {
                        config.sv_max_length_fraction = value.as_f64().ok_or_else(|| {
                            GenerateReadsError::ConfigReadError(
                                "sv_max_length_fraction".to_string(),
                                "float".to_string(),
                            )
                        })?;
                    }
                    "adapters" => {
                        // Nested map: { enabled, preset: truseq|nextera|custom, r1, r2 }.
                        let enabled = value
                            .get("enabled")
                            .and_then(|v| v.as_bool())
                            .unwrap_or(false);
                        config.adapters.enabled = enabled;
                        if enabled {
                            let preset = value
                                .get("preset")
                                .and_then(|v| v.as_str())
                                .unwrap_or("truseq")
                                .to_lowercase();
                            match preset.as_str() {
                                "truseq" => {
                                    config.adapters.r1 =
                                        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_string();
                                    config.adapters.r2 =
                                        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string();
                                }
                                "nextera" => {
                                    config.adapters.r1 = "CTGTCTCTTATACACATCT".to_string();
                                    config.adapters.r2 = "CTGTCTCTTATACACATCT".to_string();
                                }
                                "custom" => {
                                    config.adapters.r1 = value
                                        .get("r1")
                                        .and_then(|v| v.as_str())
                                        .unwrap_or("")
                                        .to_uppercase();
                                    config.adapters.r2 = value
                                        .get("r2")
                                        .and_then(|v| v.as_str())
                                        .unwrap_or("")
                                        .to_uppercase();
                                }
                                other => {
                                    error!(
                                        "Unknown adapter preset '{other}'; expected truseq|nextera|custom"
                                    );
                                    return Err(GenerateReadsError::ConfigError);
                                }
                            }
                        }
                    }
                    _ => continue,
                },
            }
        }

        RunConfiguration::check_and_log_config(&mut config)?;
        Ok(config)
    }

    pub fn check_and_log_config(config: &mut RunConfiguration) -> Result<(), GenerateReadsError> {
        // This does a final check of the configuration for valid items. It will print info
        // messages of the items, to work as a record and to assist in debugging any issues that
        // come up.
        info!("\t>read length: {}", config.read_len);
        info!("\t>coverage: {}", config.coverage);
        info!("\t>ploidy: {}", config.ploidy);
        if let Some(rate) = &config.mutation_rate {
            info!("\t>mutation rate: {}", rate)
        }
        info!("\t>paired ended: {}", config.paired_ended);
        if config.long_reads {
            info!("\t>long reads mode: enabled (short fragments produce truncated reads)");
        }
        if config.keep_short_fragments {
            if config.long_reads {
                error!(
                    "Invalid config: keep_short_fragments=true is redundant/incompatible with \
                     long_reads=true (long-read mode already keeps short fragments). Set one."
                );
                return Err(GenerateReadsError::ConfigError);
            }
            info!(
                "\t>keep short fragments: enabled (short inserts kept as genomic reads, no adapter)"
            );
        }
        if config.adapters.enabled {
            if config.long_reads {
                error!(
                    "Invalid config: adapters.enabled=true is incompatible with long_reads=true. \
                     Illumina-style 3' adapter readthrough does not occur on long-read platforms. \
                     Set long_reads=false (short-read Illumina) or adapters.enabled=false."
                );
                return Err(GenerateReadsError::ConfigError);
            }
            let acgt =
                |s: &str| !s.is_empty() && s.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'));
            if config.adapters.r1.is_empty() || config.adapters.r2.is_empty() {
                error!(
                    "Invalid config: adapters.enabled=true but an adapter sequence is empty \
                     (preset 'custom' requires non-empty r1 and r2)."
                );
                return Err(GenerateReadsError::ConfigError);
            }
            if !acgt(&config.adapters.r1) || !acgt(&config.adapters.r2) {
                error!(
                    "Invalid config: adapter sequences must contain only A/C/G/T (uppercase); got r1='{}', r2='{}'.",
                    config.adapters.r1, config.adapters.r2
                );
                return Err(GenerateReadsError::ConfigError);
            }
            info!(
                "\t>adapters: enabled (R1={}, R2={})",
                config.adapters.r1, config.adapters.r2
            );
        }
        if config.overwrite_output {
            warn!("\t>Overwriting any existing files.")
        }
        if config.minimum_mutations > 0 {
            info!(
                "\t>minimum mutations per contig: {}",
                config.minimum_mutations
            )
        }

        // No point in running if we aren't producing files
        if !(config.produce_fastq || config.produce_vcf || config.produce_bam) {
            error!("All file types set to false, no files would be produced.");
            return Err(GenerateReadsError::ConfigError);
        }

        if config.paired_ended {
            // A paired-end run needs a fragment-length source. That can be an
            // explicit fragment_mean + fragment_st_dev, OR a fragment_model — the
            // runner uses the model in preference to mean/st_dev (see runner.rs),
            // so requiring mean/st_dev when a model is present would reject a config
            // the runtime handles fine. Accept whichever is provided; only error if
            // neither is.
            let have_mean_sd = config.fragment_mean.is_some() && config.fragment_st_dev.is_some();
            let have_model = config.fragment_model.is_some();
            if !have_mean_sd && !have_model {
                error!(
                    "Paired ended is set to true, but no fragment-length source was \
                    provided: set both fragment_mean and fragment_st_dev, or supply a \
                    fragment_model."
                );
                return Err(GenerateReadsError::ConfigError);
            }
            if have_model && have_mean_sd {
                warn!(
                    "\t>Both a fragment_model and fragment_mean/fragment_st_dev are set; \
                    the fragment_model takes precedence and mean/st_dev are ignored."
                );
            }
            match (config.fragment_mean, config.fragment_st_dev) {
                (Some(mean), Some(st_dev)) => {
                    info!("\t>fragment mean: {}", mean);
                    info!("\t>fragment standard deviation: {}", st_dev);
                }
                _ => info!("\t>fragment length: from fragment_model"),
            }
            if config.produce_fastq {
                info!("\t>Producing fastq files:");
                let mut fastq_1 = config.output_dir.clone();
                fastq_1.push(format!("{}_r1.fastq.gz", config.output_filename));
                info!("\t\t> {:?}", &fastq_1);
                config.output_fastq_1 = Some(fastq_1);
                let mut fastq_2 = config.output_dir.clone();
                fastq_2.push(format!("{}_r2.fastq.gz", config.output_filename));
                info!("\t\t> {:?}", fastq_2);
                config.output_fastq_2 = Some(fastq_2);
            }
        } else {
            // single-ended
            let mut fastq_1 = config.output_dir.clone();
            fastq_1.push(format!("{}_r1.fastq.gz", config.output_filename));
            info!("\t>Producing fastq file: {:?}", &fastq_1);
            config.output_fastq_1 = Some(fastq_1);
        }
        if config.produce_vcf {
            let mut vcf = config.output_dir.clone();
            vcf.push(format!("{}.vcf.gz", config.output_filename));
            info!("\t>Producing VCF file: {:?}", &vcf);
            config.output_vcf = Some(vcf);
        }
        if config.produce_bam {
            let mut bam = config.output_dir.clone();
            bam.push(format!("{}.bam", config.output_filename));
            info!("\t>Producing BAM file: {:?}", &bam);
            config.output_bam = Some(bam);
        }

        if let Some(bed) = &config.target_bed {
            info!("\t>Target BED (generation-time filter): {:?}", bed);
        }
        if let Some(bed) = &config.mutation_regions {
            info!("\t>Custom mutation regions: {:?}", bed);
        }
        if let Some(vcf) = &config.input_vcf {
            info!("\t>Input VCF (forced variants): {:?}", vcf);
        }

        if let Some(gc_model) = &config.gc_bias_model {
            info!("\t>GC bias model: {:?}", gc_model);
            info!(
                "\t>GC bias coverage normalization: {}",
                config.gc_bias_normalize_coverage
            );
        }
        if let Some(qual_score_model) = &config.quality_score_model {
            info!("\t>Quality score model: {:?}", qual_score_model)
        }
        if let Some(mutation_model) = &config.mutation_model {
            info!("\t>Mutation model: {:?}", mutation_model)
        }
        if let Some(fragment_model) = &config.fragment_model {
            info!("\t>Fragment model: {:?}", fragment_model)
        }
        if let Some(seq_err_model) = &config.sequence_error_model {
            info!("\t>Sequence error model: {:?}", seq_err_model)
        }

        if !config.produce_bam {
            match config.num_threads {
                Some(n) => {
                    info!("\t>Maximum number of threads to use: {}", n)
                }
                None => info!("\t>parallel threads: auto (all available cores)"),
            }
        }

        match &config.rng_seed {
            Some(seed) => {
                // User supplied seed
                // Convert the string to a string vec split at the white space.
                // Seeds can be any whitespace-separated series of strings.
                for seed_term in seed.split_whitespace() {
                    config.seed_vec.push(seed_term.to_string());
                }
                info!("\t>Reproducibility seed string: {}", &seed);
                Ok(())
            }

            None => {
                // Since no seed was provided, we'll use a datetime stamp with nanoseconds
                // The seed can be any space-separated or tab-separated series of strings
                // e.g., "Every good boy does Fine"
                // seeds are case-sensitive
                let raw_string = Utc::now().format("%Y %m %d %H %M %S %f").to_string();
                for item in raw_string.split_whitespace() {
                    config.seed_vec.push(item.to_string());
                }
                info!("\t>Reproducibility seed string: {}", &raw_string);
                Ok(())
            }
        }
    }
}

#[cfg(test)]
// RunConfiguration has many fields; test setups touch only 2–4 of them, so
// `let mut c = Default::default(); c.foo = ...;` reads better than the
// struct-init form clippy suggests. Allow the lint for the whole test mod.
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn test_run_configuration() {
        let test_configuration = RunConfiguration {
            reference: PathBuf::from("Hello.world"),
            read_len: 100,
            coverage: 22,
            mutation_rate: Some(0.09),
            ploidy: 3,
            paired_ended: true,
            fragment_mean: Option::from(333.0),
            fragment_st_dev: Option::from(33.0),
            produce_fastq: false,
            produce_bam: true,
            produce_vcf: true,
            rng_seed: None,
            seed_vec: Vec::new(),
            overwrite_output: true,
            minimum_mutations: 0,
            output_dir: PathBuf::from("/my/my"),
            output_filename: String::from("Hey.hey"),
            output_fastq_1: None,
            output_fastq_2: None,
            output_vcf: None,
            output_bam: None,
            quality_score_model: None,
            mutation_model: None,
            fragment_model: None,
            sequence_error_model: None,
            target_bed: None,
            mutation_regions: None,
            input_vcf: None,
            gc_bias_model: None,
            gc_bias_normalize_coverage: true,
            num_threads: None,
            chunk_size: None,
            long_reads: false,
            keep_short_fragments: false,
            sv_rate_scale: 0.0,
            sv_max_length_fraction: 0.25,
            adapters: AdapterConfig::default(),
            subclone_model: None,
            somatic_vcf: None,
            somatic_af_scale: 1.0,
        };

        println!("{:?}", test_configuration);
        assert_eq!(test_configuration.reference, PathBuf::from("Hello.world"));
        assert_eq!(test_configuration.read_len, 100);
        assert_eq!(test_configuration.coverage, 22);
        assert_eq!(test_configuration.mutation_rate, Some(0.09));
        assert_eq!(test_configuration.ploidy, 3);
        assert!(test_configuration.paired_ended);
        assert_eq!(test_configuration.fragment_mean.unwrap(), 333.0);
        assert_eq!(test_configuration.fragment_st_dev.unwrap(), 33.0);
        assert!(!test_configuration.produce_fastq);
        assert!(test_configuration.produce_vcf);
        assert!(test_configuration.produce_bam);
        assert_eq!(test_configuration.rng_seed, None);
        assert!(test_configuration.overwrite_output);
        assert_eq!(test_configuration.output_dir, PathBuf::from("/my/my"));
        assert_eq!(test_configuration.output_filename, "Hey.hey".to_string());
    }

    #[test]
    fn test_overwrite_warn() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.overwrite_output = true;
        RunConfiguration::check_and_log_config(&mut config).unwrap();
    }

    #[test]
    fn test_paired_ended_fastq_paths() {
        // Verifies that check_and_log_config sets both fastq output paths when paired_ended=true
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        config.fragment_st_dev = Some(10.0);
        RunConfiguration::check_and_log_config(&mut config).unwrap();
        assert!(config.output_fastq_1.is_some());
        assert!(config.output_fastq_2.is_some());
    }

    #[test]
    fn test_paired_ended_fragment_model_only() {
        // A fragment_model alone satisfies the paired-end requirement — no
        // fragment_mean/st_dev needed (the runner uses the model in preference to
        // mean/st_dev). The validator must accept it and still set both fastq paths.
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.paired_ended = true;
        config.fragment_model = Some(PathBuf::from(
            "test_data/baseline_models/frag_length.canonical.json.gz",
        ));
        RunConfiguration::check_and_log_config(&mut config).unwrap();
        assert!(config.output_fastq_1.is_some());
        assert!(config.output_fastq_2.is_some());
    }

    #[test]
    fn test_no_files() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.produce_fastq = false;
        assert!(matches!(
            RunConfiguration::check_and_log_config(&mut config),
            Err(GenerateReadsError::ConfigError)
        ));
    }

    #[test]
    fn test_no_frag_mean_or_stdev() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.paired_ended = true;
        assert!(matches!(
            RunConfiguration::check_and_log_config(&mut config),
            Err(GenerateReadsError::ConfigError)
        ));
    }

    #[test]
    fn test_no_frag_mean() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.paired_ended = true;
        config.fragment_st_dev = Some(10.0);
        assert!(matches!(
            RunConfiguration::check_and_log_config(&mut config),
            Err(GenerateReadsError::ConfigError)
        ));
    }

    #[test]
    fn test_no_stdev() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        assert!(matches!(
            RunConfiguration::check_and_log_config(&mut config),
            Err(GenerateReadsError::ConfigError)
        ));
    }

    #[test]
    fn test_read_config_yaml() {
        let yaml = PathBuf::from("test_data/configs/neat_test.yml");
        let test_config = RunConfiguration::from_yaml_file(&yaml).unwrap();
        assert_eq!(
            test_config.reference,
            PathBuf::from("test_data/references/ecoli.fa")
        );
        assert_eq!(test_config.coverage, 3);
    }

    #[test]
    #[should_panic]
    fn test_bad_yaml() {
        let yaml = PathBuf::from("test_data/fake_file.yml");
        RunConfiguration::from_yaml_file(&yaml).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_missing_ref() {
        let yaml = PathBuf::from("test_data/configs/simple_template.yml");
        RunConfiguration::from_yaml_file(&yaml).unwrap();
    }

    #[test]
    fn test_creates_out_dir() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_dir = temp_dir.path().join("output");
        // Write a minimal valid yaml that points output_dir at a path that doesn't exist yet
        let yaml_content = format!(
            "reference: test_data/references/H1N1.fa\n\
             paired_ended: true\n\
             fragment_mean: 10.0\n\
             fragment_st_dev: 11.0\n\
             produce_vcf: true\n\
             overwrite_output: true\n\
             output_dir: {}\n",
            output_dir.display()
        );
        let yaml_file = temp_dir.path().join("test.yml");
        fs::write(&yaml_file, &yaml_content).unwrap();

        assert!(
            !output_dir.is_dir(),
            "output dir should not exist before the call"
        );
        RunConfiguration::from_yaml_file(&yaml_file).unwrap();
        assert!(
            output_dir.is_dir(),
            "from_yaml_file should have created the output dir"
        );
    }

    #[test]
    fn test_default_config() {
        let config = RunConfiguration::default();
        assert_eq!(config.read_len, 151);
        assert_eq!(config.coverage, 10);
        assert!(!config.paired_ended);
        assert!(config.produce_fastq);
        // sv_rate_scale defaults to 0.0 — de novo SV generation is opt-in.
        // A v1.9 upgrader with no YAML change must see identical
        // behavior to v1.9 (no SVs emitted).
        assert_eq!(config.sv_rate_scale, 0.0);
    }

    #[test]
    fn test_yaml_parses_sv_rate_scale() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = dir.path().join("c.yml");
        std::fs::write(
            &yaml,
            "reference: test_data/references/H1N1.fa\n\
             output_dir: ./\n\
             sv_rate_scale: 0.75\n",
        )
        .unwrap();
        let cfg = RunConfiguration::from_yaml_file(&yaml).unwrap();
        assert!((cfg.sv_rate_scale - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_sv_max_length_fraction_default() {
        // Default 0.25 reproduces the historical contig_len/4 cap bit-for-bit
        // (#229) so existing configs see no behavior change.
        let config = RunConfiguration::default();
        assert_eq!(config.sv_max_length_fraction, 0.25);
    }

    #[test]
    fn test_yaml_parses_sv_max_length_fraction() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = dir.path().join("c.yml");
        std::fs::write(
            &yaml,
            "reference: test_data/references/H1N1.fa\n\
             output_dir: ./\n\
             sv_max_length_fraction: 1.0\n",
        )
        .unwrap();
        let cfg = RunConfiguration::from_yaml_file(&yaml).unwrap();
        assert!((cfg.sv_max_length_fraction - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_single_ended_fastq_path() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.paired_ended = false;
        RunConfiguration::check_and_log_config(&mut config).unwrap();
        assert!(config.output_fastq_1.is_some());
        assert!(config.output_fastq_2.is_none());
    }

    #[test]
    fn test_vcf_and_bam_paths() {
        let mut config = RunConfiguration::default();
        config.reference = PathBuf::from("test_data/references/H1N1.fa");
        config.produce_vcf = true;
        config.produce_bam = true;
        RunConfiguration::check_and_log_config(&mut config).unwrap();
        assert!(config.output_vcf.is_some());
        assert!(config.output_bam.is_some());
    }

    #[test]
    fn test_from_scrape_config() {
        // We need a real file for the reference check in from_scrape_config
        let ref_path = "test_data/references/H1N1.fa";

        let mut scrape_config = HashMap::new();
        scrape_config.insert("reference".to_string(), Value::String(ref_path.to_string()));
        scrape_config.insert("read_len".to_string(), Value::Number(200.into()));
        scrape_config.insert("coverage".to_string(), Value::Number(30.into()));
        scrape_config.insert("paired_ended".to_string(), Value::Bool(true));
        scrape_config.insert("fragment_mean".to_string(), Value::Number(400.into()));
        scrape_config.insert("fragment_st_dev".to_string(), Value::Number(40.into()));

        let config = RunConfiguration::from_scrape_config(scrape_config).unwrap();

        assert_eq!(config.reference, PathBuf::from(ref_path));
        assert_eq!(config.read_len, 200);
        assert_eq!(config.coverage, 30);
        assert!(config.paired_ended);
        assert_eq!(config.fragment_mean, Some(400.0));
        assert_eq!(config.fragment_st_dev, Some(40.0));
    }

    #[test]
    fn test_from_scrape_config_dot_notation() {
        let ref_path = "test_data/references/H1N1.fa";

        let mut scrape_config = HashMap::new();
        scrape_config.insert("reference".to_string(), Value::String(ref_path.to_string()));
        // Use default for read_len
        scrape_config.insert("read_len".to_string(), Value::String(".".to_string()));

        let config = RunConfiguration::from_scrape_config(scrape_config).unwrap();

        assert_eq!(config.reference, PathBuf::from(ref_path));
        assert_eq!(config.read_len, RunConfiguration::default().read_len);
    }

    #[test]
    fn test_config_invalid_types() {
        let ref_path = "test_data/references/H1N1.fa";
        let mut scrape_config = HashMap::new();
        scrape_config.insert("reference".to_string(), Value::String(ref_path.to_string()));

        // Invalid integer
        let mut sc1 = scrape_config.clone();
        sc1.insert(
            "read_len".to_string(),
            Value::String("not_an_int".to_string()),
        );
        assert!(RunConfiguration::from_scrape_config(sc1).is_err());

        // Invalid boolean
        let mut sc2 = scrape_config.clone();
        sc2.insert(
            "paired_ended".to_string(),
            Value::String("not_a_bool".to_string()),
        );
        assert!(RunConfiguration::from_scrape_config(sc2).is_err());

        // Invalid float
        let mut sc3 = scrape_config.clone();
        sc3.insert(
            "mutation_rate".to_string(),
            Value::String("not_a_float".to_string()),
        );
        assert!(RunConfiguration::from_scrape_config(sc3).is_err());
    }

    #[test]
    fn test_config_missing_reference() {
        let scrape_config = HashMap::new();
        assert!(RunConfiguration::from_scrape_config(scrape_config).is_err());
    }

    #[test]
    fn test_config_file_not_found() {
        let mut scrape_config = HashMap::new();
        scrape_config.insert(
            "reference".to_string(),
            Value::String("non_existent_file.fasta".to_string()),
        );
        assert!(RunConfiguration::from_scrape_config(scrape_config).is_err());
    }

    #[test]
    fn test_rng_seed_accepts_integer() {
        // YAML `rng_seed: 42` (bare integer) used to panic with
        // ConfigReadError("rng_seed", "String"). It now coerces to "42" so the
        // downstream whitespace-tokenized seeder sees a single seed term.
        let ref_path = "test_data/references/H1N1.fa";
        let mut scrape_config = HashMap::new();
        scrape_config.insert("reference".to_string(), Value::String(ref_path.to_string()));
        scrape_config.insert("rng_seed".to_string(), Value::Number(42.into()));
        let config = RunConfiguration::from_scrape_config(scrape_config).unwrap();
        assert_eq!(config.rng_seed.as_deref(), Some("42"));
    }

    fn write_cfg(dir: &std::path::Path, body: &str) -> PathBuf {
        let yaml = dir.join("c.yml");
        std::fs::write(
            &yaml,
            format!("reference: test_data/references/H1N1.fa\noutput_dir: ./\n{body}"),
        )
        .unwrap();
        yaml
    }

    #[test]
    fn test_adapters_default_disabled() {
        // Default preserves pre-adapter behavior: disabled.
        assert!(!RunConfiguration::default().adapters.enabled);
    }

    #[test]
    fn test_adapters_truseq_preset() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = write_cfg(dir.path(), "adapters:\n  enabled: true\n  preset: truseq\n");
        let cfg = RunConfiguration::from_yaml_file(&yaml).unwrap();
        assert!(cfg.adapters.enabled);
        assert_eq!(cfg.adapters.r1, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        assert_eq!(cfg.adapters.r2, "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
    }

    #[test]
    fn test_adapters_nextera_preset() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = write_cfg(
            dir.path(),
            "adapters:\n  enabled: true\n  preset: nextera\n",
        );
        let cfg = RunConfiguration::from_yaml_file(&yaml).unwrap();
        assert_eq!(cfg.adapters.r1, "CTGTCTCTTATACACATCT");
        assert_eq!(cfg.adapters.r2, "CTGTCTCTTATACACATCT");
    }

    #[test]
    fn test_adapters_custom_uppercased() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = write_cfg(
            dir.path(),
            "adapters:\n  enabled: true\n  preset: custom\n  r1: acgtacgt\n  r2: tgcatgca\n",
        );
        let cfg = RunConfiguration::from_yaml_file(&yaml).unwrap();
        assert_eq!(cfg.adapters.r1, "ACGTACGT");
        assert_eq!(cfg.adapters.r2, "TGCATGCA");
    }

    #[test]
    fn test_adapters_reject_long_reads_combo() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = write_cfg(
            dir.path(),
            "long_reads: true\nadapters:\n  enabled: true\n  preset: truseq\n",
        );
        assert!(RunConfiguration::from_yaml_file(&yaml).is_err());
    }

    #[test]
    fn test_adapters_reject_custom_missing_seq() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = write_cfg(dir.path(), "adapters:\n  enabled: true\n  preset: custom\n");
        assert!(RunConfiguration::from_yaml_file(&yaml).is_err());
    }

    #[test]
    fn test_adapters_reject_non_acgt() {
        let dir = tempfile::tempdir().unwrap();
        let yaml = write_cfg(
            dir.path(),
            "adapters:\n  enabled: true\n  preset: custom\n  r1: ACGTX\n  r2: ACGT\n",
        );
        assert!(RunConfiguration::from_yaml_file(&yaml).is_err());
    }

    #[test]
    fn test_rng_seed_accepts_quoted_string() {
        // Multi-term seeds remain available via a quoted string.
        let ref_path = "test_data/references/H1N1.fa";
        let mut scrape_config = HashMap::new();
        scrape_config.insert("reference".to_string(), Value::String(ref_path.to_string()));
        scrape_config.insert(
            "rng_seed".to_string(),
            Value::String("42 hello world".to_string()),
        );
        let config = RunConfiguration::from_scrape_config(scrape_config).unwrap();
        assert_eq!(config.rng_seed.as_deref(), Some("42 hello world"));
    }
}
