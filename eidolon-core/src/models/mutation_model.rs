//! The mutation model
use crate::rng::{NeatRng, NeatRngError};
use flate2::read::GzDecoder;
use log::error;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::{self, collections::HashMap};
use thiserror::Error;

use crate::models::{
    indel_model::{IndelModel, IndelModelError},
    lib::{model_reader, model_writer},
    quality_scores::QualityModelError,
    sequencing_error_model::SeqModelError,
    snp_trinuc_model::{SnpTrinucError, SnpTrinucModel},
    sv_model_defaults::default_sv_model,
};
use crate::{
    models::snp_trinuc_model::TrinucFrame,
    structs::{
        distributions::{DiscreteDistribution, DistributionErrors},
        nucleotides::{Nucleotide, allowed_vec},
        sv_model::SvModel,
        transition_matrix::{TransitionMatrix, TransitionMatrixError},
        variants::{Variant, VariantError, VariantType},
    },
};

#[derive(Error, Debug)]
pub enum MutationModelError {
    #[error("Error in inputs to MutationModel")]
    InputError,
    #[error("Mutation model returned RNG error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Error retrieving reference sequence block")]
    ReferenceRetrievalError,
    #[error("Mutation model IO Error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Error Generating mutations")]
    GenerateMutationError,
    #[error("Mutation model returned an error during variant creation: {0}")]
    VariantInitationError(#[from] VariantError),
    #[error("Mutation model returned an error during SNP initation: {0}")]
    SnpeGnerationError(#[from] SnpTrinucError),
    #[error("Mutation model returned error during transition matrix initation: {0}")]
    TransMatrixInitationError(#[from] TransitionMatrixError),
    #[error("Mutation model returned error during indel model initation: {0}")]
    IndelModelInitationError(#[from] IndelModelError),
    #[error("Mutation model returned error during quality model initiation: {0}")]
    QualityModelInitationError(#[from] QualityModelError),
    #[error("Mutation model returned error during sequence error model initation: {0}")]
    SeqErrModelError(#[from] SeqModelError),
    #[error("Mutation model returned an error with a distribution model: {0}")]
    DistributionError(#[from] DistributionErrors),
    #[error("Mutation model return and error: {0}")]
    ModelInitiationError(String),
    #[error("Error reading Mutation Model from file {0}")]
    SerdeError(#[from] serde_json::Error),
}

fn allowed_variant_types() -> Vec<VariantType> {
    vec![
        VariantType::SNP,
        VariantType::Insertion,
        VariantType::Deletion,
    ]
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MutationModel {
    // This is the model for mutations, the same construct used by the python version, basically.
    // The idea is to store all the data we need to create mutations in one construct.
    //
    // mutation rate is the average rate of mutation for the dataset. It is used to calculate how
    // many variants to add to the dataset.
    pub mutation_rate: f64,
    // Homozygous frequency is the fraction of mutations that are homozygous, meaning the alternate
    // was inherited from both parents (in the case of humans). The definition of "homozygous"
    // is ambiguous in polyploid organisms. We'll take it to mean "on all ploids"
    pub homozygous_frequency: f64,
    // If a variant occurs, this is how frequently each type occurs. There are only 3 types at the
    // moment, but this will get expanded out in time. This data item will need to be expanded to
    // account for it. The canonical order is:
    //    1. SNP
    //    2. Insertion
    //    3. Deletion
    pub variant_dist: DiscreteDistribution<VariantType>,
    // These hold all the statistical model data we need to apply the mutations with this model
    statistical_models: StatisticalModels,
    // Optional learned model for symbolic / structural variants
    // (`<DEL>` / `<DUP>` / `<CNV>`). Populated by `gen-mut-model` when
    // the training VCF carries enough symbolic-SV observations to fit;
    // `None` otherwise, and `None` whenever this model was deserialized
    // from a JSON file that predates the field (the `#[serde(default)]`
    // attribute handles the missing-field case).
    #[serde(default)]
    pub sv_model: Option<SvModel>,
    // Store a reference to the NeatRng for this run
}

static DATA_FILE: &[u8] = include_bytes!("model_data/default_mutation_model.json.gz");

impl MutationModel {
    pub fn from_raw_data(
        average_mutation_rate: f64,
        homozygous_frequency: f64,
        variant_probs: Vec<f64>,
        snp_transition_frequency: HashMap<(Nucleotide, Nucleotide), f64>,
        trinuc_frequency: HashMap<TrinucFrame, f64>,
        trinuc_transition_frequency: HashMap<(TrinucFrame, TrinucFrame), f64>,
        ins_lengths: Vec<usize>,
        ins_weights: Vec<f64>,
        del_lengths: Vec<usize>,
        del_weights: Vec<f64>,
        transition_matrix_override: Option<TransitionMatrix>,
    ) -> Result<Self, MutationModelError> {
        let transition_matrix = if let Some(tm) = transition_matrix_override {
            tm
        } else {
            let mut temp_trans_matrix: [[f64; 4]; 4] = [[0.0; 4]; 4];
            for (key, value) in snp_transition_frequency {
                temp_trans_matrix[key.0 as usize][key.1 as usize] = value;
            }
            TransitionMatrix::from(
                temp_trans_matrix[Nucleotide::A as usize],
                temp_trans_matrix[Nucleotide::C as usize],
                temp_trans_matrix[Nucleotide::G as usize],
                temp_trans_matrix[Nucleotide::T as usize],
            )?
        };
        // build transition matrices from data for snps and trinucs
        let snp_trinuc_model =
            SnpTrinucModel::from_raw_data(trinuc_frequency, trinuc_transition_frequency)?;
        // Fall back to the default indel model when no indel data was observed.
        // The insertion-vs-deletion balance is NOT stored here: it is already
        // variant_probs[1] and variant_probs[2] (the `variant_dist` weights), which
        // is what drives generation. This model carries the length distributions only.
        let indel_denom = variant_probs[1] + variant_probs[2];
        let indel_model = if indel_denom == 0.0 || ins_lengths.is_empty() || del_lengths.is_empty()
        {
            IndelModel::default()?
        } else {
            IndelModel::from_raw_data(ins_lengths, ins_weights, del_lengths, del_weights)?
        };
        let statistical_models = StatisticalModels {
            transition_matrix,
            indel_model,
            snp_trinuc_model,
        };
        // SNP, Insertion, Deletion, respectively.
        if variant_probs.len() != 3 {
            error!(
                "Input to mutation model from raw incorrect. Variant probs should have a length of 3: SNP, INS, DEL"
            );
            return Err(MutationModelError::InputError);
        }
        let variant_distribution =
            DiscreteDistribution::new(&variant_probs, &(allowed_variant_types()))?;
        Ok(MutationModel {
            mutation_rate: average_mutation_rate,
            variant_dist: variant_distribution,
            homozygous_frequency,
            statistical_models,
            // SV training fills this in via a separate code path on the
            // gen-mut-model side; constructors that don't see SV data
            // start it as None (callers can mutate the field after).
            sv_model: None,
        })
    }

    // Returns Result because it deserializes an embedded model file; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, MutationModelError> {
        // Creating the default model based on the default for the original NEAT.
        let reader = GzDecoder::new(DATA_FILE);
        let mut data: MutationModel =
            serde_json::from_reader(reader).map_err(MutationModelError::SerdeError)?;
        // Attach the bundled default SV model in code rather than baking it
        // into `default_mutation_model.json.gz` — keeps the SV defaults
        // auditable in `sv_model_defaults.rs` instead of buried in a
        // gzipped blob, and means a parameter tweak doesn't require
        // regenerating the embedded file. `sv_rate_scale` defaults to
        // 0.0 in `RunConfiguration`, so generation stays opt-in.
        if data.sv_model.is_none() {
            data.sv_model = Some(default_sv_model());
        }
        Ok(data)
    }

    pub fn from_file(filename: &PathBuf) -> Result<Self, MutationModelError> {
        let data: MutationModel = model_reader(filename)?;
        Ok(data)
    }

    pub fn write_to_file(&self, filename: &PathBuf) -> Result<(), MutationModelError> {
        model_writer(self, filename)?;
        Ok(())
    }

    fn generate_genotype(
        &self,
        ploidy: usize,
        is_homozygous: bool,
    ) -> Result<Vec<usize>, MutationModelError> {
        // "Homozygous" is ambiguous for polyploid organisms, so we'll just take "heterozygous" to
        // mean roughly half the reads will have the variant, to keep it simple
        // The is_homozygous flag is expected to be randomly determined in practice
        if is_homozygous {
            Ok(vec![1; ploidy])
        } else {
            let num_ploids = ploidy / 2;
            Ok([vec![1; num_ploids], vec![0; ploidy - num_ploids]].concat())
        }
    }

    /// Relative per-trinucleotide-context mutation propensity `w(ctx)`, for
    /// context-weighted mutation placement (#372). Delegates to the SNP trinucleotide
    /// model's `snp_distro`. A uniform distribution → context-neutral placement.
    pub fn context_weights(&self) -> Result<HashMap<TrinucFrame, f64>, MutationModelError> {
        Ok(self.statistical_models.snp_trinuc_model.context_weights()?)
    }

    pub fn generate_mutation(
        &self,
        // A pointer to the reference sequence to be mutated
        reference_sequence: &Vec<Nucleotide>,
        // This is the start position of the variant, the POS field of a VCF
        variant_location: usize,
        // The ploidy of this organism. We'll make determinations about heterzygous
        // v homozygous based on this
        ploidy: usize,
        // Force this variant's type instead of sampling it. Used by context-weighted
        // placement (#372), which routes SNPs and indels to different position
        // distributions and so must know the type before this call. `None` samples
        // the type from the model's variant_dist as before.
        forced_type: Option<VariantType>,
        // The run's rng. Needed to determine type then generate the mutation
        rng: &mut NeatRng,
    ) -> Result<Variant, MutationModelError> {
        // Select a genotype for the variant
        let mut genotype =
            self.generate_genotype(ploidy, rng.gen_bool(self.homozygous_frequency)?)?;
        // Select a type of mutation (or use the caller-forced type).
        let mut variant_type = match forced_type {
            Some(t) => t,
            None => self.variant_dist.sample(rng.random()?)?,
        };
        // Unmask before use: soft-masked (a/c/g/t) and post-N-substitution bases
        // must not appear as-is in VCF REF/ALT or FASTQ reads.
        let ref_base = check_base(reference_sequence[variant_location]);
        // For insertions and SNPs, the reference will be the ref_base
        let mut reference: Vec<Nucleotide> = vec![ref_base];
        // For deletions, the alternate will be the ref_base
        let mut alternate: Vec<Nucleotide> = vec![ref_base];
        match variant_type {
            VariantType::SNP => {
                if (variant_location < 1) || ((variant_location + 1) >= reference_sequence.len()) {
                    // We don't want to accidentally go out of bounds, so we'll
                    // just pick a random base.
                    alternate = pick_random_snp(ref_base, rng)?;
                } else {
                    // In this case we use the trinculeotide model
                    let trinuc_reference = [
                        check_base(reference_sequence[variant_location - 1]),
                        check_base(reference_sequence[variant_location]),
                        check_base(reference_sequence[variant_location + 1]),
                    ];
                    let alternate_base = self
                        .statistical_models
                        .snp_trinuc_model
                        .generate_snp(rng.random()?, &trinuc_reference)?;

                    // Degenerate model context (sparse data, all-zero row): fall back to
                    // a uniform random pick rather than producing a same-base "mutation".
                    if alternate_base == ref_base {
                        alternate = pick_random_snp(ref_base, rng)?;
                    } else {
                        alternate = vec![alternate_base];
                    }
                }
            }
            VariantType::Insertion => {
                let length = self
                    .statistical_models
                    .indel_model
                    .new_insert_length(rng.random()?)?;
                let insertion_vec = self
                    .statistical_models
                    .indel_model
                    .generate_random_insertion(length, rng)?;
                alternate = [reference.clone(), insertion_vec.clone()].concat();
            }
            VariantType::Deletion => {
                // todo Deletions are a bitch. I need to think about them.
                let length = self
                    .statistical_models
                    .indel_model
                    .new_delete_length(rng.random()?)?;
                // +1 is so that we grab a base for the reference in the VCF. This is similar
                // to how we appended bases to the reference in the insertion model.
                if (variant_location + length + 1) > reference_sequence.len() {
                    // Too close to the end, so let's skip this
                    let ref_base = check_base(reference_sequence[variant_location]);
                    alternate = pick_random_snp(ref_base, rng)?;
                    variant_type = VariantType::SNP;
                } else {
                    reference = reference_sequence
                        .get(variant_location..(variant_location + length + 1))
                        .unwrap()
                        .iter()
                        .map(|&b| check_base(b))
                        .collect();
                }
            }
            _ => {
                panic!("Unsupported option for generating mutation")
            }
        };

        if reference.is_empty() || alternate.is_empty() || reference == alternate {
            error!(
                "Malformed mutation: {:?} {:?} {:?}",
                variant_type, reference, alternate
            );
            return Err(MutationModelError::GenerateMutationError);
        }
        Ok(Variant::new(
            variant_type,
            variant_location,
            &reference,
            &alternate,
            &mut genotype,
        )?)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct StatisticalModels {
    // This struct links the Mutation model to the other statistical models for modeling different
    // variant types. If new variant types are added, then we will need to expand this struct to
    // include them.
    transition_matrix: TransitionMatrix,
    indel_model: IndelModel,
    snp_trinuc_model: SnpTrinucModel,
}

fn pick_random_snp(
    ref_base: Nucleotide,
    rng: &mut NeatRng,
) -> Result<Vec<Nucleotide>, MutationModelError> {
    // Grab allowed nucs
    let mut allowed: Vec<Nucleotide> = allowed_vec();
    // Throw out the reference base
    allowed.retain(|&x| x != ref_base);
    // Pick one of the remaining
    let alt = rng.choose(&allowed)?;
    // sanity check
    assert_ne!(alt, ref_base);
    Ok(vec![alt])
}

fn check_base(nuc: Nucleotide) -> Nucleotide {
    // For now we'll replace occasional N's with A's.
    if nuc.is_masked() {
        return nuc.get_unmasked_base();
    }
    match nuc {
        Nucleotide::N => Nucleotide::A,
        _ => nuc,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::snp_trinuc_model::TrinucFrame;
    use crate::structs::nucleotides::Nucleotide::{A, C, G, T};

    #[test]
    fn test_model_read_write() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_file = temp_dir.path().join("test.json.gz");
        let model: MutationModel = MutationModel::default().unwrap();
        assert_eq!(model.mutation_rate, 0.0010987132390211135);
        model.write_to_file(&output_file).unwrap();
        let loaded = MutationModel::from_file(&output_file).unwrap();
        assert_eq!(loaded.mutation_rate, model.mutation_rate);
    }

    /// The embedded default model must carry real trinucleotide biology, not a
    /// context-neutral or corrupted weight vector.
    ///
    /// #372 wired `context_weights()` into placement, but that only helps if the
    /// model FILE carries meaningful weights. `default_mutation_model_bkup.json.gz`
    /// holds the older context-neutral default (CpG share exactly 6.25% = uniform);
    /// the replacement must actually elevate CpG, which is where SBS1 lives and the
    /// single most mutable context in any real genome.
    ///
    /// This is a data assertion, not a code one: it fails if the shipped default is
    /// regenerated wrongly, which is invisible to Ts/Tv and to every caller-level
    /// metric.
    #[test]
    fn default_model_context_weights_elevate_cpg() {
        let model = MutationModel::default().unwrap();
        let weights = model.context_weights().unwrap();
        assert_eq!(weights.len(), 64, "expected 64 trinucleotide frames");

        let mean: f64 = weights.values().sum::<f64>() / weights.len() as f64;
        assert!(mean > 0.0, "context weights are all zero");

        let cpg = [(A, C, G), (C, C, G), (G, C, G), (T, C, G)];
        let mut ratios = Vec::new();
        for ctx in cpg {
            let w = *weights
                .get(&TrinucFrame::from(ctx))
                .expect("CpG frame missing from context weights");
            ratios.push(w / mean);
        }
        let group = ratios.iter().sum::<f64>() / ratios.len() as f64;
        eprintln!(
            "[default-model] CpG weight/mean ratios ACG={:.2} CCG={:.2} GCG={:.2} TCG={:.2}; \
             group mean {:.2}x",
            ratios[0], ratios[1], ratios[2], ratios[3], group
        );

        assert!(
            group >= 2.0,
            "the embedded default model does not elevate CpG: the four CpG contexts \
             average {group:.2}x the mean context weight (expected >=2x; a correctly \
             built model runs ~4.4x, and 1.0x means context-neutral placement). \
             A monotonically-increasing-with-index weight vector here is the \
             signature of a snp_distro that was cumulated twice — de-cumulating it \
             a second time recovers default_trinuc_model.json.gz exactly."
        );
    }

    #[test]
    fn test_default_homozygous_frequency_is_realistic() {
        // Guards against the NEAT-lineage default that silently shipped for years:
        // NEAT2 used 0.010 and NEAT4 0.001, both faithfully ported here, which
        // produce a het/hom ratio of ~99-999 — unrealistic for any diploid, where
        // a resequenced sample is ~1/3 homozygous-alt. A near-zero value here means
        // the embedded default model regressed; keep it in a sane diploid range.
        let model = MutationModel::default().unwrap();
        assert!(
            (0.2..=0.5).contains(&model.homozygous_frequency),
            "default homozygous_frequency {} is outside the realistic diploid range \
             (0.2..=0.5); a low value inflates the het/hom ratio",
            model.homozygous_frequency
        );
    }

    #[test]
    fn test_generate_genotype_homozygous() {
        let model = MutationModel::default().unwrap();
        let geno = model.generate_genotype(4, true).unwrap();
        assert_eq!(geno, vec![1, 1, 1, 1]);
    }

    #[test]
    fn test_generate_genotype_heterozygous() {
        let model = MutationModel::default().unwrap();
        let geno = model.generate_genotype(4, false).unwrap();
        // half 1s, half 0s — order may vary but count must hold
        assert_eq!(geno.iter().filter(|&&x| x == 1).count(), 2);
        assert_eq!(geno.iter().filter(|&&x| x == 0).count(), 2);
    }

    #[test]
    fn test_check_base_passes_through_acgt() {
        assert_eq!(check_base(A), A);
        assert_eq!(check_base(C), C);
        assert_eq!(check_base(G), G);
        assert_eq!(check_base(T), T);
    }

    #[test]
    fn test_check_base_replaces_n() {
        assert_eq!(check_base(Nucleotide::N), Nucleotide::A);
    }

    #[test]
    fn test_check_base_unmasks_masked() {
        // Masked bases should be replaced with their unmasked equivalent
        assert_eq!(check_base(Nucleotide::Maskeda), Nucleotide::A);
        assert_eq!(check_base(Nucleotide::Maskedc), Nucleotide::C);
    }

    #[test]
    fn test_generate_mutation_produces_variant() {
        let model = MutationModel::default().unwrap();
        let mut rng =
            NeatRng::new_from_seed(&vec!["test".to_string(), "seed".to_string()]).unwrap();
        // 20-base sequence with room on both sides of position 5
        let seq = vec![A, C, G, T, A, C, G, T, A, C, G, T, A, C, G, T, A, C, G, T];
        let variant = model.generate_mutation(&seq, 5, 2, None, &mut rng).unwrap();
        // alternate must differ from reference
        assert_ne!(
            variant.reference.as_slice(),
            variant.alternate.as_literal().unwrap()
        );
        assert_eq!(variant.location, 5);
    }

    #[test]
    fn test_generate_mutation_deterministic() {
        let model = MutationModel::default().unwrap();
        let seq = vec![A, C, G, T, A, C, G, T, A, C, G, T, A, C, G, T, A, C, G, T];
        let seed = vec!["det".to_string(), "seed".to_string()];
        let mut rng1 = NeatRng::new_from_seed(&seed).unwrap();
        let mut rng2 = NeatRng::new_from_seed(&seed).unwrap();
        let v1 = model
            .generate_mutation(&seq, 5, 2, None, &mut rng1)
            .unwrap();
        let v2 = model
            .generate_mutation(&seq, 5, 2, None, &mut rng2)
            .unwrap();
        assert_eq!(v1, v2);
    }

    #[test]
    fn test_from_raw_data_no_indels() {
        // Exercises the indel_denom == 0 || empty vecs fallback to IndelModel::default()
        let snp_trans: HashMap<(Nucleotide, Nucleotide), f64> = HashMap::from([
            ((Nucleotide::A, Nucleotide::C), 0.25),
            ((Nucleotide::A, Nucleotide::G), 0.5),
            ((Nucleotide::A, Nucleotide::T), 0.25),
        ]);
        let trinuc_freq: HashMap<TrinucFrame, f64> = HashMap::new();
        let trinuc_trans: HashMap<(TrinucFrame, TrinucFrame), f64> = HashMap::new();
        let result = MutationModel::from_raw_data(
            0.001,
            0.5,
            vec![1.0, 0.0, 0.0], // SNP only, no indels
            snp_trans,
            trinuc_freq,
            trinuc_trans,
            vec![], // no insertion data
            vec![],
            vec![], // no deletion data
            vec![],
            None,
        );
        assert!(
            result.is_ok(),
            "Expected Ok with no indels, got {:?}",
            result
        );
        let model = result.unwrap();
        assert_eq!(model.mutation_rate, 0.001);
    }

    #[test]
    fn from_raw_data_leaves_sv_model_none() {
        // `from_raw_data` is the SNP/indel-only training path. SV stats
        // come through a separate channel on `gen-mut-model::runner`, so
        // the constructor must start with `sv_model: None`. Callers that
        // observed SVs assign the field after construction.
        let snp_trans = HashMap::from([
            ((Nucleotide::A, Nucleotide::C), 0.25),
            ((Nucleotide::A, Nucleotide::G), 0.5),
            ((Nucleotide::A, Nucleotide::T), 0.25),
        ]);
        let model = MutationModel::from_raw_data(
            0.001,
            0.5,
            vec![1.0, 0.0, 0.0],
            snp_trans,
            HashMap::new(),
            HashMap::new(),
            vec![],
            vec![],
            vec![],
            vec![],
            None,
        )
        .unwrap();
        assert!(model.sv_model.is_none());
    }

    #[test]
    fn pre_v1_10_model_json_without_sv_field_loads_with_sv_model_none() {
        // Back-compat contract: a MutationModel JSON written by v1.9 has
        // no `sv_model` key. The `#[serde(default)]` attribute must
        // populate `None` so existing model files keep loading after the
        // schema bump. Round-trip through serde to lock this in.
        let original = MutationModel::default().unwrap();
        let mut as_value = serde_json::to_value(&original).unwrap();
        // Strip the new field to simulate a v1.9-era file.
        as_value.as_object_mut().unwrap().remove("sv_model");
        assert!(
            !as_value.as_object().unwrap().contains_key("sv_model"),
            "test setup: sv_model should be absent after strip"
        );
        let reloaded: MutationModel = serde_json::from_value(as_value).unwrap();
        assert!(reloaded.sv_model.is_none());
        // Other fields must still be intact — i.e., the strip didn't
        // mangle anything else and the load wasn't lenient about extra
        // missing fields.
        assert_eq!(reloaded.mutation_rate, original.mutation_rate);
        assert_eq!(reloaded.homozygous_frequency, original.homozygous_frequency);
    }

    #[test]
    fn populated_sv_model_round_trips_through_mutation_model_json() {
        // End-to-end JSON round-trip with a non-empty SvModel attached.
        // Catches accidental `#[serde(skip)]` on the field or any
        // misconfiguration that would silently drop SV data on write.
        use crate::structs::sv_model::SvModel;
        use crate::structs::variants::SvType;

        let mut model = MutationModel::default().unwrap();
        let mut sv = SvModel::default();
        sv.per_base_rate = 1.5e-7;
        sv.homozygous_frequency = 0.4;
        sv.type_probabilities.insert(SvType::Del, 0.7);
        sv.type_probabilities.insert(SvType::Dup, 0.3);
        sv.length_log_normal.insert(SvType::Del, (7.0, 1.5));
        model.sv_model = Some(sv);

        let json = serde_json::to_string(&model).unwrap();
        let back: MutationModel = serde_json::from_str(&json).unwrap();
        let back_sv = back.sv_model.expect("sv_model dropped on round-trip");
        assert!((back_sv.per_base_rate - 1.5e-7).abs() < 1e-12);
        assert_eq!(back_sv.type_probabilities.get(&SvType::Del), Some(&0.7));
        assert_eq!(
            back_sv.length_log_normal.get(&SvType::Del),
            Some(&(7.0, 1.5))
        );
    }
}
