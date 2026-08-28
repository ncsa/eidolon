pub mod file_tools;
pub mod models;
pub mod rng;
pub mod structs;
extern crate log;
extern crate serde;
extern crate serde_json;

// TEMPORARY — proves the clippy CI gate actually fails. Removed in the next commit.
#[allow(dead_code)]
fn clippy_gate_probe(x: usize) -> usize {
    return x + 1;
}
