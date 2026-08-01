extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate simplelog;

pub mod compare_vcfs;
pub mod filter_reads;
pub mod gen_bam_models;
pub mod gen_cancer_reads;
pub mod gen_frag_length_model;
pub mod gen_gc_bias_model;
pub mod gen_mut_model;
pub mod gen_reads;
pub mod gen_seq_error_model;
pub mod validate;

use eidolon_core::{self, file_tools::file_io::create_output_file};
use std::env;
use std::time;

use clap::{Arg, ArgAction, Command, value_parser};
use eidolon_core::file_tools::folder_tools::check_parent;
use log::*;
use simplelog::{ColorChoice, CombinedLogger, Config, TermLogger, TerminalMode, WriteLogger};
use std::path::PathBuf;
use thiserror::Error;

use crate::{
    compare_vcfs::errors::CompareVcfsError, filter_reads::errors::FilterReadsError,
    gen_bam_models::errors::GenBamModelsError, gen_cancer_reads::errors::GenCancerReadsError,
    gen_frag_length_model::errors::GenFragLengthModelError,
    gen_gc_bias_model::errors::GenGcBiasModelError, gen_mut_model::errors::GenMutationModelError,
    gen_reads::errors::GenerateReadsError, gen_seq_error_model::errors::GenSeqErrorModelError,
    validate::errors::ValidateError,
};
/// This script parses arguments and checks them before submitting to the submodules, which currently
/// include `gen-reads` and `filter-files`. As more are added, this can be expanded or refactored
/// argument parsing, timing, and other concerns will occur here.
#[derive(Error, Debug)]
pub enum NeatErrors {
    // Errors specific to each submodule
    #[error("Error while generating read dataset {0}")]
    GenerateReadsError(#[from] GenerateReadsError),
    #[error("Error while filtering reads with bed file {0}")]
    FilterReadsError(#[from] FilterReadsError),
    #[error("Error while generating mutation model {0}")]
    GenMutModelError(#[from] GenMutationModelError),
    #[error("Error while generating sequencing error model {0}")]
    GenSeqErrorModel(#[from] GenSeqErrorModelError),
    #[error("Error while generating fragment length model {0}")]
    GenFragLengthModel(#[from] GenFragLengthModelError),
    #[error("Error while generating GC bias model {0}")]
    GenGcBiasModel(#[from] GenGcBiasModelError),
    #[error("Error while generating BAM-derived models {0}")]
    GenBamModels(#[from] GenBamModelsError),
    #[error("Error while comparing VCFs {0}")]
    CompareVcfs(#[from] CompareVcfsError),
    #[error("Error while generating cancer read dataset {0}")]
    GenCancerReads(#[from] GenCancerReadsError),
    #[error("Validation failed: {0}")]
    Validate(#[from] ValidateError),
}

fn neat_commands() -> [Command; 10] {
    // These are the submodule commands. Any new commands added should go here.
    let configuration_arg = Arg::new("configuration_yaml")
        .long("configuration-yaml")
        .short('c')
        .help("Path to configuration file.")
        .action(ArgAction::Set)
        .exclusive(true)
        .default_missing_value("")
        .value_parser(value_parser!(PathBuf))
        .required(true);
    [
        Command::new("gen-reads")
            .about("Generates reads for an input dataset")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("filter-reads")
            .about("Filters the output of gen-reads")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("gen-mut-model")
            .about("Generate a mutation model based on an input mutations file")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("gen-seq-error-model")
            .about("Generate a sequencing error model from a FASTQ file")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("gen-frag-length-model")
            .about("Generate a fragment length model from a BAM or SAM file")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("gen-gc-bias-model")
            .about("Generate a GC bias model from a reference FASTA and aligned BAM")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("gen-bam-models")
            .about("Build multiple models (frag-length, GC bias) from one BAM in a single pass")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("compare-vcfs")
            .about("Compare a NEAT-simulated golden VCF against a downstream variant-caller VCF")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        Command::new("gen-cancer-reads")
            .about("Simulate a tumor/normal mixture (two gen-reads passes + merge)")
            .arg_required_else_help(true)
            .arg(&configuration_arg),
        // Takes file paths directly rather than a config YAML: it is a diagnostic run
        // ad hoc and from harness scripts, where writing a YAML per check is friction.
        Command::new("validate")
            .about("Check an emitted FASTQ/VCF against what its consumers accept")
            .arg_required_else_help(true)
            .arg(
                Arg::new("files")
                    .help("Artifacts to validate (.fq/.fastq/.vcf/.bam, optionally .gz)")
                    .action(ArgAction::Append)
                    .required(true)
                    .value_parser(value_parser!(PathBuf)),
            )
            .arg(
                Arg::new("format")
                    .long("format")
                    .help("Override format detection when the extension is absent or wrong")
                    .action(ArgAction::Set)
                    .value_parser(["fastq", "vcf", "bam"]),
            ),
    ]
}

fn main() -> Result<(), NeatErrors> {
    // parse the arguments from the command line
    let cmd = Command::new(env!("CARGO_CRATE_NAME"))
        .multicall(true)
        .subcommand(
            Command::new("eidolon")
                .version(env!("CARGO_PKG_VERSION"))
                .arg_required_else_help(true)
                .subcommand_value_name("SUB-COMMAND")
                .subcommand_help_heading("SUB-COMMANDS")
                .arg(
					// Verbosity of the WRITTEN .neat.log (WriteLogger). Default info: the
					// per-base debug/trace events in gen-reads fire ~coverage×read_len×ref_bp
					// times, so a default-trace run wrote a multi-GB log and burned most of the
					// wall-clock on log I/O. Pass --log-level debug/trace to opt into verbosity.
					// The on-screen log is always info. (default_value must agree with the info
					// fallback below — a default_value of "trace" here silently overrode it.)
                    Arg::new("log_level")
                        .long("log-level")
                        .help("Verbosity of the written .neat.log (trace|debug|info|warn|error|off; default info). The on-screen log is always info.")
                        .action(ArgAction::Set)
                        .num_args(0..=1)
                        .value_parser(["trace", "debug", "info", "warn", "error", "off"])
                        .default_value("info")
                        .default_missing_value("debug")
                )
                .arg(
					// This arg controls where the log is written. The default is <current_working_dir>/.neat.log
                    Arg::new("log_dest")
                        .long("log-dest")
                        .help("Sets the log destination (full path with full filename) for the written log")
                        .action(ArgAction::Set)
                        .value_parser(value_parser!(PathBuf))
                        .default_missing_value("")
                )
                .subcommands(neat_commands())
        )
        .subcommands(neat_commands());

    // This parses the command arguments
    let matches = cmd.get_matches();
    let mut subcommand = matches.subcommand();
    // Default is `info`, not `trace`. The per-base debug events in
    // gen-reads (sequencing-error generation, read-position logs, etc.)
    // fire on the order of `coverage × read_len × reference_bp`, so on a
    // human-sized reference a default-trace run produces a multi-GB
    // `.neat.log` and burns most of the wall-clock time on log writes.
    // Users who actually want verbose output pass `--log-level debug` or
    // `trace` explicitly.
    let mut level_filter = "info".to_string();
    let mut log_dest = env::current_dir().unwrap();
    log_dest.push(".neat.log");
    if let Some(("eidolon", cmd)) = subcommand {
        if cmd.contains_id("log_level") {
            level_filter = cmd
                .get_one::<String>("log_level")
                .expect("Must enter value for log-level")
                .to_string();
        }
        if cmd.contains_id("log_dest") {
            log_dest = cmd
                .get_one::<PathBuf>("log_dest")
                .expect("Must provide a path with log-dest")
                .to_path_buf();
        }
        subcommand = cmd.subcommand();
    }

    // log filter
    let lev_filter_choice = {
        match level_filter.to_lowercase().as_str() {
            "trace" => LevelFilter::Trace,
            "debug" => LevelFilter::Debug,
            "info" => LevelFilter::Info,
            "warn" => LevelFilter::Warn,
            "error" => LevelFilter::Error,
            "off" => LevelFilter::Off,
            _ => unreachable!("Parser should only allow one of the above values."),
        }
    };

    // Check that the parent dir exists
    let log_destination =
        check_parent(&log_dest, true).expect("Log destination parent path doesn't exist");
    let log_file = create_output_file(log_destination, true)
        .unwrap_or_else(|_| panic!("Error creating log file: {:?}", log_destination));
    // Set up the logger for the run
    CombinedLogger::init(vec![
        TermLogger::new(
            LevelFilter::Info,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        ),
        WriteLogger::new(lev_filter_choice, Config::default(), log_file),
    ])
    .expect("Error creating log file!");
    info!("////////////// Welcome to eidolon! \\\\\\\\\\\\\\\\\\\\\\\\\\\\");
    let start = time::Instant::now();
    match subcommand {
        Some(("gen-reads", _)) => {
            if let Some(("gen-reads", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::FilterReadsError(FilterReadsError::CliError(
                        "Must supply either a configuration file or a reference file to run NEAT!"
                            .to_string(),
                    )));
                }
                info!("Running gen-reads to generate read data.");
                let result = gen_reads::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenerateReadsError(error)),
                    Ok(()) => info!("eidolon gen-reads completed successfully"),
                }
            }
        }
        Some(("filter-reads", _)) => {
            if let Some(("filter-reads", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon filter-reads");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::FilterReadsError(FilterReadsError::CliError(
                        "Must supply either a configuration file or a reference file to run NEAT!"
                            .to_string(),
                    )));
                }
                info!("Running filter-reads to filter eidolon-generated read data.");
                let result = filter_reads::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::FilterReadsError(error)),
                    Ok(()) => info!("eidolon filter-reads completed succesfully"),
                }
            }
        }
        Some(("gen-mut-model", _)) => {
            if let Some(("gen-mut-model", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon gen-mut-model");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::FilterReadsError(FilterReadsError::CliError(
                        "Must supply either a configuration file or a reference file to run NEAT!"
                            .to_string(),
                    )));
                }
                info!("Running gen-mut-model to generate a model for use in eidolon.");
                let result = gen_mut_model::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenMutModelError(error)),
                    Ok(()) => info!("eidolon gen-mut-model completed succesfully"),
                }
            }
        }
        Some(("gen-seq-error-model", _)) => {
            if let Some(("gen-seq-error-model", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon gen-seq-error-model");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::FilterReadsError(FilterReadsError::CliError(
                        "Must supply a configuration file to run gen-seq-error-model!".to_string(),
                    )));
                }
                info!("Running gen-seq-error-model to generate a sequencing error model.");
                let result = gen_seq_error_model::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenSeqErrorModel(error)),
                    Ok(()) => info!("eidolon gen-seq-error-model completed successfully"),
                }
            }
        }
        Some(("gen-frag-length-model", _)) => {
            if let Some(("gen-frag-length-model", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon gen-frag-length-model");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::FilterReadsError(FilterReadsError::CliError(
                        "Must supply a configuration file to run gen-frag-length-model!"
                            .to_string(),
                    )));
                }
                info!("Running gen-frag-length-model to generate a fragment length model.");
                let result = gen_frag_length_model::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenFragLengthModel(error)),
                    Ok(()) => info!("eidolon gen-frag-length-model completed successfully"),
                }
            }
        }
        Some(("gen-gc-bias-model", _)) => {
            if let Some(("gen-gc-bias-model", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon gen-gc-bias-model");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::GenGcBiasModel(
                        GenGcBiasModelError::ConfigError(
                            "Must supply a valid configuration file to run gen-gc-bias-model!"
                                .to_string(),
                        ),
                    ));
                }
                info!("Running gen-gc-bias-model to generate a GC bias model.");
                let result = gen_gc_bias_model::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenGcBiasModel(error)),
                    Ok(()) => info!("eidolon gen-gc-bias-model completed successfully"),
                }
            }
        }
        Some(("gen-bam-models", _)) => {
            if let Some(("gen-bam-models", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon gen-bam-models");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::GenBamModels(GenBamModelsError::ConfigError(
                        "Must supply a valid configuration file to run gen-bam-models!".to_string(),
                    )));
                }
                info!("Running gen-bam-models to build multiple BAM-derived models in one pass.");
                let result = gen_bam_models::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenBamModels(error)),
                    Ok(()) => info!("eidolon gen-bam-models completed successfully"),
                }
            }
        }
        Some(("compare-vcfs", _)) => {
            if let Some(("compare-vcfs", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                info!("Running eidolon compare-vcfs");
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::CompareVcfs(
                        CompareVcfsError::ConfigurationError(
                            "Must supply a valid configuration file to run compare-vcfs!"
                                .to_string(),
                        ),
                    ));
                }
                info!("Running compare-vcfs to classify golden vs called variants.");
                let result = compare_vcfs::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::CompareVcfs(error)),
                    Ok(()) => info!("eidolon compare-vcfs completed successfully"),
                }
            }
        }
        Some(("gen-cancer-reads", _)) => {
            if let Some(("gen-cancer-reads", cmd)) = subcommand
                && cmd.contains_id("configuration_yaml")
            {
                let file = cmd
                    .get_one::<PathBuf>("configuration_yaml")
                    .expect("Must provide a path with configuration-yaml")
                    .to_path_buf();

                if !file.is_file() {
                    return Err(NeatErrors::GenCancerReads(
                        GenCancerReadsError::ConfigError(
                            "Must supply a valid configuration file to run gen-cancer-reads!"
                                .to_string(),
                        ),
                    ));
                }
                info!("Running gen-cancer-reads to simulate a tumor/normal mixture.");
                let result = gen_cancer_reads::main(&file);
                match result {
                    Err(error) => return Err(NeatErrors::GenCancerReads(error)),
                    Ok(()) => info!("eidolon gen-cancer-reads completed successfully"),
                }
            }
        }
        Some(("validate", _)) => {
            if let Some(("validate", cmd)) = subcommand {
                let files: Vec<PathBuf> = cmd
                    .get_many::<PathBuf>("files")
                    .expect("clap marks files as required")
                    .cloned()
                    .collect();
                let format = cmd.get_one::<String>("format").map(|f| match f.as_str() {
                    "fastq" => validate::Format::Fastq,
                    "bam" => validate::Format::Bam,
                    _ => validate::Format::Vcf,
                });
                info!("Running eidolon validate on {} file(s)", files.len());
                validate::run(&files, format)?;
                info!("eidolon validate completed successfully");
            }
        }
        _ => unreachable!("Parser should ensure no other subcommand choice is made."),
    }

    let elapsed_time = time::Instant::now() - start;
    if elapsed_time.as_millis() < 1000 {
        info!(
            "Processing finished in {} milliseconds",
            elapsed_time.as_millis()
        );
    } else if elapsed_time.as_secs() > 300 {
        info!(
            "Processing finished in {} minutes",
            ((elapsed_time.as_secs() as f64) / 60.0)
        )
    } else {
        info!("Processing finished in {} seconds", elapsed_time.as_secs());
    }
    Ok(())
}
