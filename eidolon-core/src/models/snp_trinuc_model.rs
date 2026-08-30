//! The following section are the models for each type of variant. In order to create the variant,
//! we need to model its statistical property. NEAT included two types of variants: SNPs and Indels.
//! Note that indels are actually two types, insertions and deletions, but they are usually classed
//! together in the literature and are usually short. Insertions are often slips that cause sections
//! to be duplicated, but NEAT made no attempt to distinguish between the types of insertions or
//! deletions, since they are treated similarly by variant calling software.
use crate::models::lib::{model_reader, model_writer};
use crate::rng::NeatRngError;
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors};
use crate::structs::nucleotides::Nucleotide;
use crate::structs::nucleotides::Nucleotide::{A, C, G, N, T};
use crate::structs::transition_matrix::{TransitionMatrix, TransitionMatrixError};
use flate2::read::GzDecoder;
use itertools::Itertools;
use lazy_static::lazy_static;
use log::debug;
use serde;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::Index;
use std::{fmt, io, path::PathBuf};
use thiserror::Error;
use vectorize;

#[derive(Error, Debug)]
pub enum SnpTrinucError {
    #[error("SNP variant model reported an RNG error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("SNP variant model reported an error from Distributions: {0}")]
    DistributionError(#[from] DistributionErrors),
    #[error("SNP variant model reported an error with transition matrices: {0}")]
    TransMatrixError(#[from] TransitionMatrixError),
    #[error("SNP variant model reported an error generating a SNP.")]
    GenerateSnpError,
    #[error("SNP Trinuc model return an IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("Error Indexing the trinucleotide: {0}")]
    IndexError(usize),
    #[error("Serde error building default model: {0}")]
    SerdeError(#[from] serde_json::Error),
}

/// These are all the trinucleotide combinations.
const SNP_TRINUCS: [(Nucleotide, Nucleotide, Nucleotide); 64] = [
    (A, A, A), // AAA
    (A, A, C), // AAC
    (A, A, G), // AAG
    (A, A, T), // AAT
    (A, C, A), // ACA
    (A, C, C), // ACC
    (A, C, G), // ACG
    (A, C, T), // ACT
    (A, G, A), // AGA
    (A, G, C), // AGC
    (A, G, G), // AGG
    (A, G, T), // AGT
    (A, T, A), // ATA
    (A, T, C), // ATC
    (A, T, G), // ATG
    (A, T, T), // ATT
    (C, A, A), // CAA
    (C, A, C), // CAC
    (C, A, G), // CAG
    (C, A, T), // CAT
    (C, C, A), // CCA
    (C, C, C), // CCC
    (C, C, G), // CCG
    (C, C, T), // CCT
    (C, G, A), // CGA
    (C, G, C), // CGC
    (C, G, G), // CGG
    (C, G, T), // CGT
    (C, T, A), // CTA
    (C, T, C), // CTC
    (C, T, G), // CTG
    (C, T, T), // CTT
    (G, A, A), // GAA
    (G, A, C), // GAC
    (G, A, G), // GAG
    (G, A, T), // GAT
    (G, C, A), // GCA
    (G, C, C), // GCC
    (G, C, G), // GCG
    (G, C, T), // GCT
    (G, G, A), // GGA
    (G, G, C), // GGC
    (G, G, G), // GGG
    (G, G, T), // GGT
    (G, T, A), // GTA
    (G, T, C), // GTC
    (G, T, G), // GTG
    (G, T, T), // GTT
    (T, A, A), // TAA
    (T, A, C), // TAC
    (T, A, G), // TAG
    (T, A, T), // TAT
    (T, C, A), // TCA
    (T, C, C), // TCC
    (T, C, G), // TCG
    (T, C, T), // TCT
    (T, G, A), // TGA
    (T, G, C), // TGC
    (T, G, G), // TGG
    (T, G, T), // TGT
    (T, T, A), // TTA
    (T, T, C), // TTC
    (T, T, G), // TTG
    (T, T, T), // TTT
];

lazy_static! {
    static ref ALL_FRAMES: Vec<TrinucFrame> = {
        let mut all_frames = Vec::new();
        for trinuc in SNP_TRINUCS {
            all_frames.push(TrinucFrame::from(trinuc))
        }
        all_frames
    };

    static ref ALIAS_MAP: HashMap<TrinucFrame, TrinucFrame> = {
        // This builds an alias map. For each frame in all frames,
        // We assign a value of the context featureing N in the middle position.
        //      AAA => ANA ([0, 0, 0]: [0, 5, 0])
        //      ATA => ANA ([0, 3, 0]: [0, 5, 0])
        //      GTC => GNC ([2, 3, 1]: [2, 5, 1])
        // and so on. Building this to save time
        let all_frames = ALL_FRAMES.clone();
        let mut alias_map = HashMap::new();
        for frame in all_frames {
            alias_map.insert(
                frame,
                TrinucFrame::from((frame[0], N, frame[2]))
            );
        }
        alias_map
    };

    static ref ALL_CONTEXTS: Vec<TrinucFrame> = {
        let alias_map = ALIAS_MAP.clone();
        let frames: Vec<TrinucFrame> = alias_map
            .values()
            .map(|s| s.convert())
            .collect();

        frames.iter().unique().cloned().collect::<Vec<_>>()
    };

    static ref CONTEXT_FRAME_MAP: HashMap<TrinucFrame, Vec<TrinucFrame>> = {
        // This builds the contexts for the trinucleotides. The 5 in the middle
        // represents any one of 4 trinculeotides (unknown).
        //    ANA => [AAA, ACA, AGA, ATA]
        //    ANT => [AAT, ACT, AGT, ATT]
        // and so on;
        let mut frame_map: HashMap<TrinucFrame, Vec<TrinucFrame>> = HashMap::new();
        for i in 0..4 {
            for k in 0..4 {
                let mut frame_list = Vec::new();
                for frame in ALL_FRAMES.clone().into_iter() {
                    let frame_0: usize = frame[0].into();
                    let frame_2: usize = frame[2].into();
                    if frame_0 == i && frame_2 == k {
                        frame_list.push(frame)
                    }
                }
                frame_map.insert(
                    TrinucFrame::from((i.into(), N, k.into())),
                    frame_list
                );
            }
        }
        frame_map
    };

}

// We need all these derivations to write these out to file
#[derive(Debug, Eq, PartialEq, Hash, Clone, Copy, Serialize, Deserialize)]
pub enum TrinucFrame {
    Frame(Nucleotide, Nucleotide, Nucleotide),
}

impl fmt::Display for TrinucFrame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut display: String = String::new();
        for i in 0..3 {
            display.push(self[i].into())
        }
        write!(f, "{}", display)
    }
}

impl From<(Nucleotide, Nucleotide, Nucleotide)> for TrinucFrame {
    fn from((n1, n2, n3): (Nucleotide, Nucleotide, Nucleotide)) -> Self {
        Self::Frame(n1, n2, n3)
    }
}

impl TrinucFrame {
    pub fn convert(self) -> TrinucFrame {
        self
    }
}

impl From<[&Nucleotide; 3]> for TrinucFrame {
    fn from(trinuc: [&Nucleotide; 3]) -> Self {
        Self::Frame(*trinuc[0], *trinuc[1], *trinuc[2])
    }
}

impl From<&[Nucleotide; 3]> for TrinucFrame {
    fn from(trinuc: &[Nucleotide; 3]) -> Self {
        Self::Frame(trinuc[0], trinuc[1], trinuc[2])
    }
}

impl Index<usize> for TrinucFrame {
    type Output = Nucleotide;
    fn index(&self, i: usize) -> &Nucleotide {
        let TrinucFrame::Frame(b0, b1, b2) = self;
        match i {
            0 => b0,
            1 => b1,
            2 => b2,
            _ => panic!("Index out of range: {}", i),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpTrinucModel {
    // Relative weights given to each SNP frame. Ultimately this will be imputed from data.
    snp_distro: DiscreteDistribution<usize>,
    // The transition matrix is the chance of mutating the middle base from A, C, T, or G to a
    // different base (4x4 matrix with 0s on the diagonal).
    #[serde(with = "vectorize")]
    trinuc_distros: HashMap<TrinucFrame, TransitionMatrix>,
}

static DATA_FILE: &[u8] = include_bytes!("model_data/default_trinuc_model.json.gz");

impl SnpTrinucModel {
    pub fn from_raw_data(
        trinuc_frequency: HashMap<TrinucFrame, f64>,
        trinuc_transition_frequency: HashMap<(TrinucFrame, TrinucFrame), f64>,
    ) -> Result<Self, SnpTrinucError> {
        let all_contexts = ALL_CONTEXTS.clone();
        let alias_map = ALIAS_MAP.clone();
        let all_frames = ALL_FRAMES.clone();

        // Per-context 4x4 accumulators: [from_base_idx][to_base_idx] += prob
        let mut accum: HashMap<TrinucFrame, [[f64; 4]; 4]> = HashMap::new();
        for context in &all_contexts {
            accum.insert(*context, [[0.0f64; 4]; 4]);
        }

        for ((from_trinuc, to_trinuc), &prob) in &trinuc_transition_frequency {
            if let Some(&context) = alias_map.get(from_trinuc) {
                let row: usize = from_trinuc[1].into();
                let col: usize = to_trinuc[1].into();
                if row < 4
                    && col < 4
                    && let Some(matrix) = accum.get_mut(&context)
                {
                    matrix[row][col] += prob;
                }
            }
        }

        let mut trinuc_distros: HashMap<TrinucFrame, TransitionMatrix> = HashMap::new();
        for context in &all_contexts {
            let w = accum[context];
            let has_data = w.iter().any(|row| row.iter().sum::<f64>() > 0.0);
            let matrix = if has_data {
                TransitionMatrix::from(w[0], w[1], w[2], w[3])?
            } else {
                TransitionMatrix::default()?
            };
            trinuc_distros.insert(*context, matrix);
        }

        // Weight each of the 64 frames by its observed mutation frequency
        let snp_weights: Vec<f64> = all_frames
            .iter()
            .map(|frame| *trinuc_frequency.get(frame).unwrap_or(&1e-6f64))
            .collect();
        let snp_distro = DiscreteDistribution::new(&snp_weights, &(0..all_frames.len()).collect())?;

        Ok(SnpTrinucModel {
            snp_distro,
            trinuc_distros,
        })
    }

    /// Relative per-context mutation propensity `w(ctx)`, reconstructed from the
    /// 64-frame `snp_distro` (built at fit time from each trinucleotide's observed
    /// mutation frequency = mutations/occurrences). Returns normalized weights keyed
    /// by frame. gen-reads uses this to weight WHERE mutations land so trinucleotide
    /// signatures reproduce in the output (#372). A uniform `snp_distro` (default /
    /// untrained model) yields all-equal weights → context-neutral placement (the
    /// pre-#372 behavior).
    pub fn context_weights(&self) -> Result<HashMap<TrinucFrame, f64>, SnpTrinucError> {
        // snp_distro holds cumulative-normalized weights over the frames in ALL_FRAMES
        // order (its values are 0..len). Difference consecutive cumulatives to recover
        // each frame's normalized weight.
        let cumulative = self.snp_distro.weights()?;
        let mut weights = HashMap::with_capacity(ALL_FRAMES.len());
        let mut prev = 0.0f64;
        for (idx, frame) in ALL_FRAMES.iter().enumerate() {
            let cum = cumulative.get(idx).copied().unwrap_or(prev);
            weights.insert(*frame, (cum - prev).max(0.0));
            prev = cum;
        }
        Ok(weights)
    }

    pub fn default_minimal() -> Result<Self, SnpTrinucError> {
        // Creating the default trinuc bias model for snps. In this model, all trinucleotides
        // mutate with equal probability and middle base mutates with the same probability no matter
        // the context (the default transition matrix).

        // One weight for each snp frame, making them all equally
        let mut snp_weights: Vec<f64> = Vec::with_capacity(16);
        let mut trinuc_distros: HashMap<TrinucFrame, TransitionMatrix> = HashMap::new();
        let all_contexts = ALL_CONTEXTS.clone();
        debug!("All contexts: {:?}", all_contexts);
        for frame in &all_contexts {
            trinuc_distros.insert(*frame, TransitionMatrix::default()?);
            snp_weights.push(1.0);
        }
        let snp_distr =
            DiscreteDistribution::new(&snp_weights, &(0..all_contexts.len()).collect())?;
        Ok(SnpTrinucModel {
            snp_distro: snp_distr,
            trinuc_distros,
        })
    }

    // Returns Result because it deserializes an embedded model file; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, SnpTrinucError> {
        // This default will read in data from the original NEAT trinuc model to create a more realistic trinuc
        // bias.
        let reader = GzDecoder::new(DATA_FILE);
        let data: SnpTrinucModel =
            serde_json::from_reader(reader).map_err(SnpTrinucError::SerdeError)?;
        Ok(data)
    }

    #[allow(unused)]
    /// we will write utilities to use this, eventually
    fn from_file(filename: &PathBuf) -> Result<Self, SnpTrinucError> {
        // Uses the serde_json crate to read a quality model from file
        let data: SnpTrinucModel = model_reader(filename).unwrap();
        Ok(data)
    }

    #[allow(unused)]
    /// we will write utilities to use this, eventually
    fn write_out(&self, filename: &PathBuf) -> std::io::Result<()> {
        // Uses the serde_json crate to write out the json form of the model. This will help us
        // create base datasets from old neat data, and give us a way to write out models that are
        // generated from user data.
        model_writer(self, filename).unwrap();
        Ok(())
    }

    pub fn generate_snp(
        &self,
        rand: f64,
        trinuc_reference: &[Nucleotide; 3],
    ) -> Result<Nucleotide, SnpTrinucError> {
        // We shouldn't have N's here. Basically, this matches the correct trinuc from the enum,
        // then uses that as the index for the trinuc matrix of interest.
        // eliminating the need to copy the rng or a pointer to it, if possible
        let trinuc = TrinucFrame::from(trinuc_reference);
        let alias_map = &ALIAS_MAP.clone();
        let alias_trinuc = alias_map[&trinuc];
        let matrix = self.trinuc_distros[&alias_trinuc].clone();

        let nuc: Nucleotide = matrix[&trinuc[1]].sample(rand)?;
        Ok(nuc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_model_write_read() {
        let output_file: PathBuf = PathBuf::from("test.json.gz");
        let model = SnpTrinucModel::default().unwrap();
        let frame = TrinucFrame::from((A, N, G));
        assert!(model.trinuc_distros.contains_key(&frame));
        let result = model.write_out(&output_file);
        assert_eq!(result.unwrap(), ());
        fs::remove_file(output_file).unwrap();
    }

    #[test]
    fn context_weights_reflect_trinuc_frequency() {
        // A model fit with one context (TCT) 10x more frequent than another (AAA)
        // should expose w(TCT) ≈ 10 × w(AAA), normalized, so gen-reads can weight
        // placement toward the mutable contexts (#372).
        let hot = TrinucFrame::from((T, C, T));
        let cold = TrinucFrame::from((A, A, A));
        let mut freq: HashMap<TrinucFrame, f64> = HashMap::new();
        freq.insert(hot, 10.0);
        freq.insert(cold, 1.0);
        let model = SnpTrinucModel::from_raw_data(freq, HashMap::new()).unwrap();
        let w = model.context_weights().unwrap();
        assert!(
            w[&hot] > w[&cold],
            "w(hot)={} !> w(cold)={}",
            w[&hot],
            w[&cold]
        );
        assert!(
            w[&hot] / w[&cold] > 5.0,
            "ratio {} should be ~10",
            w[&hot] / w[&cold]
        );
        let total: f64 = w.values().sum();
        assert!(
            (total - 1.0).abs() < 1e-6,
            "weights should be normalized, sum={total}"
        );
    }

    #[test]
    fn test_generate_snp_changes_base() {
        let model = SnpTrinucModel::default().unwrap();
        let cases = [
            ([A, A, G], A),
            ([G, C, T], C),
            ([C, G, A], G),
            ([T, T, C], T),
        ];
        for (trinuc, middle) in &cases {
            let result = model.generate_snp(0.5, trinuc).unwrap();
            assert_ne!(result, *middle, "SNP must change the middle base");
            assert!(
                matches!(result, A | C | G | T),
                "result must be a standard base"
            );
        }
    }

    #[test]
    fn test_generate_snp_deterministic() {
        let model = SnpTrinucModel::default().unwrap();
        let trinuc = [G, A, T];
        let r1 = model.generate_snp(0.3, &trinuc).unwrap();
        let r2 = model.generate_snp(0.3, &trinuc).unwrap();
        assert_eq!(r1, r2);
    }

    #[test]
    fn test_default_minimal_has_16_contexts() {
        let model = SnpTrinucModel::default_minimal().unwrap();
        // 4 flanking bases × 4 flanking bases = 16 ANA-style contexts
        assert_eq!(model.trinuc_distros.len(), 16);
    }

    #[test]
    fn test_default_minimal_generate_snp() {
        let model = SnpTrinucModel::default_minimal().unwrap();
        let trinuc = [A, C, G];
        let result = model.generate_snp(0.5, &trinuc).unwrap();
        assert_ne!(result, C, "minimal model SNP must still change the base");
        assert!(matches!(result, A | C | G | T));
    }

    #[test]
    fn test_model_write_read_round_trip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("trinuc.json.gz");
        let model = SnpTrinucModel::default().unwrap();
        model.write_out(&path).unwrap();
        let loaded: SnpTrinucModel = model_reader(&path).unwrap();
        // Both models should have the same 16 ANG-style context keys
        assert_eq!(loaded.trinuc_distros.len(), model.trinuc_distros.len());
        let frame = TrinucFrame::from((A, N, G));
        assert!(loaded.trinuc_distros.contains_key(&frame));
    }
}
