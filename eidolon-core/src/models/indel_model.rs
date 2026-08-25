//! The following section are the models for each type of variant. In order to create the variant,
//! we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
//! Note that indels are actually two types, insertions and deletions, but they are usually classed
//! together in the literature and are usually short. Insertions are often slips that cause sections
//! to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
//! deletions, since they are treated similarly by variant calling software.
//!
//! We have the variants separated out in the variants struct, into Insertion and Deletion, but they
//! will share this model to avoid duplicating code, since code-wise they are closely linked.
use crate::rng::{NeatRng, NeatRngError};
use crate::{
    models::lib::{model_reader, model_writer},
    structs::{
        distributions::{DiscreteDistribution, DistributionErrors},
        nucleotides::{Nucleotide, allowed_vec},
    },
};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::{io, path::PathBuf};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IndelModelError {
    #[error("Indel Model returned an Rng Error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Indel model returnd an error while initializng distributions.")]
    DistributionInitializationError(#[from] DistributionErrors),
    #[error("Error during input/output: {0}")]
    IoError(#[from] io::Error),
    #[error("Serde error building default model: {0}")]
    SerdeError(#[from] serde_json::Error),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndelModel {
    // NEAT 2.1 stored a single `INDEL_FREQ` dict keyed by *signed* indel length, so
    // the sign carried insertion-vs-deletion and the magnitude carried length.
    // Splitting that into two unsigned distributions dropped the sign, which is why
    // an explicit `insertion_probability` was added here. It was redundant: the
    // ins:del balance is already the second and third weights of the mutation
    // model's `variant_dist`, which is what actually drives generation. The field
    // was computed from `variant_dist`, serialized, and never read.
    ins_dist: DiscreteDistribution<usize>,
    del_dist: DiscreteDistribution<usize>,
}

static DATA_FILE: &[u8] = include_bytes!("model_data/default_indel_model.json.gz");

impl IndelModel {
    pub fn from_raw_data(
        ins_lens: Vec<usize>,
        ins_weights: Vec<f64>,
        del_lens: Vec<usize>,
        del_weights: Vec<f64>,
    ) -> Result<Self, IndelModelError> {
        let ins_dist = DiscreteDistribution::new(&ins_weights, &ins_lens)?;
        let del_dist = DiscreteDistribution::new(&del_weights, &del_lens)?;
        Ok(IndelModel { ins_dist, del_dist })
    }
    // Returns Result because it deserializes an embedded model file; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, IndelModelError> {
        // Default Indel model from the original NEAT, scaled by the lowest value and rounded.
        let reader = GzDecoder::new(DATA_FILE);
        let data: IndelModel =
            serde_json::from_reader(reader).map_err(IndelModelError::SerdeError)?;
        Ok(data)
    }

    pub fn from(filename: &PathBuf) -> Result<Self, IndelModelError> {
        let data: IndelModel = model_reader(filename)?;
        Ok(data)
    }

    pub fn write_to_file(&self, filename: &PathBuf) -> Result<(), IndelModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    pub fn new_insert_length(&self, rand: f64) -> Result<usize, IndelModelError> {
        // This function uses the insertion weights to choose an insertion length from the list.
        Ok(self.ins_dist.sample(rand)?)
    }

    pub fn new_delete_length(&self, rand: f64) -> Result<usize, IndelModelError> {
        // This function is the same as above, but uses the deletion lengths instead.
        Ok(self.del_dist.sample(rand)?)
    }

    pub fn generate_random_insertion(
        &self,
        length: usize,
        rng: &mut NeatRng,
    ) -> Result<Vec<Nucleotide>, IndelModelError> {
        // We could refine this with a nucleotide bias matrix. Maybe it would make a difference,
        // but probably not, since the presence of the insertion is more important than it's content,
        // for this use. If there was a call for it, maybe.
        let mut insertion_vec = Vec::new();
        let allowed_nucs = allowed_vec();
        for _ in 0..length {
            // since our range is restricted, we are covered
            let nuc = rng.choose(&allowed_nucs).unwrap();
            insertion_vec.push(nuc)
        }
        Ok(insertion_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::*;
    use tempfile::tempdir;

    fn make_model() -> IndelModel {
        IndelModel {
            ins_dist: DiscreteDistribution::new(&vec![100.0, 1.0], &vec![1, 2]).unwrap(),
            del_dist: DiscreteDistribution::new(&vec![1.0, 1000.0], &vec![3, 4]).unwrap(),
        }
    }

    #[test]
    fn test_indel_model() {
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();

        let indel_model = IndelModel {
            ins_dist: DiscreteDistribution::new(&vec![100.0, 1.0], &vec![1, 2]).unwrap(),
            del_dist: DiscreteDistribution::new(&vec![1.0, 1000.0], &vec![3, 4]).unwrap(),
        };

        let length = indel_model
            .new_insert_length(rng.random().unwrap())
            .unwrap();
        if length > 0 {
            assert!((length == 1) || (length == 2));
        }
    }

    #[test]
    fn test_generate_insertion() {
        let mut rng = NeatRng::new_from_seed(&vec![
            "Hello".to_string(),
            "Cruel".to_string(),
            "World".to_string(),
        ])
        .unwrap();

        let indel_model = IndelModel {
            ins_dist: DiscreteDistribution::new(&vec![100.0, 1.0], &vec![1, 2]).unwrap(),
            del_dist: DiscreteDistribution::new(&vec![1.0, 1000.0], &vec![3, 4]).unwrap(),
        };

        let insertion = indel_model.generate_random_insertion(10, &mut rng).unwrap();
        assert_eq!(insertion.len(), 10);
        assert!(!insertion.contains(&Nucleotide::A));

        let insertion2 = indel_model.generate_random_insertion(12, &mut rng).unwrap();
        assert_eq!(insertion2.len(), 12);
        assert_eq!(insertion2, Vec::from([G, T, T, G, G, G, T, C, T, C, T, T]));
    }

    #[test]
    fn test_new_delete_length() {
        let model = make_model();
        // del_dist heavily weights length 4 (1000:1)
        let length = model.new_delete_length(0.5).unwrap();
        assert!(length == 3 || length == 4);
        let length_high = model.new_delete_length(0.999).unwrap();
        assert_eq!(length_high, 4);
    }

    #[test]
    fn test_indel_model_file_round_trip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("indel_model.json.gz");
        let model = IndelModel::default().unwrap();
        model.write_to_file(&path).unwrap();
        let loaded = IndelModel::from(&path).unwrap();
        // Assert the surviving content: both length distributions must sample
        // identically after a write/read cycle.
        for r in [0.0, 0.25, 0.5, 0.75, 0.999] {
            assert_eq!(
                loaded.new_insert_length(r).unwrap(),
                model.new_insert_length(r).unwrap(),
                "insertion length distribution changed across a file round trip at r={r}"
            );
            assert_eq!(
                loaded.new_delete_length(r).unwrap(),
                model.new_delete_length(r).unwrap(),
                "deletion length distribution changed across a file round trip at r={r}"
            );
        }
    }

    /// Every model file built before `insertion_probability` was removed still
    /// carries the key — the bundled `tools/cosmic_*.json.gz` included, and so do
    /// `model_data/default_mutation_model.json.gz` and its `_bkup` sibling. That is
    /// deliberate: the shipped files are NOT regenerated to drop the key, because a
    /// user's own models cannot be, and the two must stay readable by the same code
    /// path. `IndelModel` therefore carries no `deny_unknown_fields`, and this test
    /// pins that serde keeps ignoring the key rather than erroring.
    ///
    /// The key is inert wherever it appears: nothing reads it, and the ins:del
    /// balance comes from the mutation model's `variant_dist`. Do not "tidy" it out
    /// of the shipped data files — that would buy nothing and would make a
    /// freshly-built model differ from a legacy one for no functional reason.
    #[test]
    fn legacy_model_with_insertion_probability_still_loads() {
        let legacy = r#"{
            "insertion_probability": 0.4538979885714955,
            "ins_dist": {"values": [1, 2], "weights": [0.5, 1.0]},
            "del_dist": {"values": [3, 4], "weights": [0.5, 1.0]}
        }"#;
        let model: IndelModel =
            serde_json::from_str(legacy).expect("a pre-removal model file must still deserialize");
        // and the distributions it does carry must be usable
        assert_eq!(model.new_insert_length(0.0).unwrap(), 1);
        assert_eq!(model.new_insert_length(0.9).unwrap(), 2);
        assert_eq!(model.new_delete_length(0.0).unwrap(), 3);
        assert_eq!(model.new_delete_length(0.9).unwrap(), 4);
    }

    #[test]
    fn test_indel_model_from_raw_data() {
        let ins_lens = vec![1, 2, 3];
        let ins_weights = vec![0.5, 0.3, 0.2];
        let del_lens = vec![1, 5];
        let del_weights = vec![0.9, 0.1];
        let model =
            IndelModel::from_raw_data(ins_lens, ins_weights, del_lens, del_weights).unwrap();
        assert_eq!(model.new_insert_length(0.0).unwrap(), 1);
        assert_eq!(model.new_delete_length(0.0).unwrap(), 1);
        assert_eq!(model.new_delete_length(0.95).unwrap(), 5);
    }

    #[test]
    fn test_generate_random_insertion_empty() {
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        let model = IndelModel::default().unwrap();
        let insertion = model.generate_random_insertion(0, &mut rng).unwrap();
        assert!(insertion.is_empty());
    }

    #[test]
    fn test_generate_random_insertion_content() {
        let mut rng = NeatRng::new_from_seed(&vec!["test".to_string()]).unwrap();
        let model = IndelModel::default().unwrap();
        let insertion = model.generate_random_insertion(100, &mut rng).unwrap();
        // Check that only valid nucleotides are generated
        let allowed = allowed_vec();
        for nuc in insertion {
            assert!(allowed.contains(&nuc));
        }
    }
}
