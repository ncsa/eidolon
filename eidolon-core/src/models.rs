//! Here are various models from the original NEAT project, written by Zach Stephens in Python2
//! These are rust implementations, and contain a lot of statistical data extracted from the original
//! compressed model data. Some of this may need to be converted to Json ultimately.

pub mod fragment_length;
pub mod gc_bias_model;
pub mod indel_model;
mod lib;
// The module stays private -- it is internal plumbing for the other models. Only the pieces a
// caller legitimately needs are re-exported: the provenance stamp (#708), which model_parity
// must strip before comparing, and the format version it is checked against.
pub use lib::{MODEL_FORMAT_VERSION, ModelProvenance, PROVENANCE_KEY};
pub mod mutation_model;
pub mod quality_scores;
pub mod sequencing_error_model;
pub mod snp_trinuc_model;
pub mod sv_model_defaults;
