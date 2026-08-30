//! Measuring how far eidolon's output sits from real sequencing data.
//!
//! See `metrics` for the numbers and why each one is here. The short version: a caller tuned
//! on data with no artifacts calibrates its false-positive filters against nothing.
pub mod metrics;
pub mod reader;
