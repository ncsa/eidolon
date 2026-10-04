use crate::models::lib::{model_reader, model_writer};
use crate::structs::nucleotides::Nucleotide;
use serde_derive::{Deserialize, Serialize};
use std::path::PathBuf;
use thiserror::Error;

// Used as the serde default for legacy model files that pre-date the window_size field.
// The Default impl also uses this, but the default model is always uniform (is_uniform=true),
// so the value is never read during simulation.
fn default_window_size() -> usize {
    100
}

#[derive(Error, Debug)]
pub enum GcBiasModelError {
    #[error("GC bias model must contain exactly 101 weights")]
    InvalidBinCount,
    #[error("GC bias weights must be finite and non-negative")]
    InvalidWeight,
    #[error("GC bias model must contain at least one positive weight")]
    EmptyModel,
    #[error("GC bias window size must be at least 1")]
    InvalidWindowSize,
    #[error("GC bias model IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Per-GC-content weight table used to bias fragment retention during read simulation.
///
/// Weights are stored as 101 `f64` values indexed by integer GC percentage (0 = 0%, 100 = 100%).
/// The scale is relative, not absolute: only the ratio of a fragment's weight to `max_weight`
/// matters. A weight of 0.0 means that GC% is never sequenced; equal weights mean no bias.
///
/// `window_size` records the number of bases used when computing per-window GC fractions during
/// training (and later during coverage estimation). Storing it in the model ensures the same
/// window is used at application time regardless of read length.
///
/// The `is_uniform` flag is computed at construction time and lets callers short-circuit
/// sampling loops when no bias is present.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcBiasModel {
    weights_by_percent_gc: Vec<f64>,
    is_uniform: bool,
    #[serde(default = "default_window_size")]
    window_size: usize,
}

impl Default for GcBiasModel {
    fn default() -> Self {
        Self {
            weights_by_percent_gc: vec![1.0; 101],
            is_uniform: true,
            window_size: default_window_size(),
        }
    }
}

impl GcBiasModel {
    pub fn from_file(path: &PathBuf) -> Result<Self, GcBiasModelError> {
        let data: GcBiasModel = model_reader(path)?;
        // Re-validate rather than trusting stored is_uniform; catches hand-edited files
        Self::from_weights(data.weights_by_percent_gc, data.window_size)
    }

    pub fn write_to_file(&self, path: &PathBuf) -> Result<(), GcBiasModelError> {
        model_writer(self, path)?;
        Ok(())
    }

    pub fn from_weights(
        weights_by_percent_gc: Vec<f64>,
        window_size: usize,
    ) -> Result<Self, GcBiasModelError> {
        if window_size == 0 {
            return Err(GcBiasModelError::InvalidWindowSize);
        }
        if weights_by_percent_gc.len() != 101 {
            return Err(GcBiasModelError::InvalidBinCount);
        }
        let first_item = weights_by_percent_gc[0];
        if first_item.is_infinite() || first_item.is_nan() || first_item < 0.0 {
            return Err(GcBiasModelError::InvalidWeight);
        }
        let mut is_uniform = true;
        for item in &weights_by_percent_gc[1..] {
            if item.is_nan() || (*item < 0.0) || item.is_infinite() {
                return Err(GcBiasModelError::InvalidWeight);
            } else if *item != first_item && is_uniform {
                is_uniform = false;
            }
        }
        if is_uniform && first_item == 0.0 {
            return Err(GcBiasModelError::EmptyModel);
        }
        Ok(Self {
            weights_by_percent_gc,
            is_uniform,
            window_size,
        })
    }

    pub fn is_uniform(&self) -> bool {
        self.is_uniform
    }

    pub fn window_size(&self) -> usize {
        self.window_size
    }

    pub fn weight_for_gc_fraction(&self, gc_fraction: f64) -> f64 {
        let index: usize = if gc_fraction < 0.0 {
            0
        } else if gc_fraction > 1.0 {
            100
        } else {
            (gc_fraction * 100.0).round() as usize
        };
        self.weights_by_percent_gc[index]
    }

    /// Returns the model weight for the GC fraction of `sequence`, ignoring N bases.
    /// An empty slice or an all-N slice returns 1.0 (neutral) so N-masked regions are
    /// never filtered out regardless of the model's other weights. Soft-masked bases count
    /// as the base they mask, as `gen-gc-bias-model` counts them (#771).
    pub fn weight_for_sequence(&self, sequence: &[Nucleotide]) -> f64 {
        let gc_count = sequence
            .iter()
            .filter(|nuc| matches!(nuc, Nucleotide::G | Nucleotide::C))
            .count() as f64;
        let sequence_count = sequence
            .iter()
            .filter(|base| **base != Nucleotide::N)
            .count() as f64;
        if sequence_count == 0.0 {
            return 1.0;
        }
        let gc_fraction = gc_count / sequence_count;
        self.weight_for_gc_fraction(gc_fraction)
    }

    pub fn max_weight(&self) -> f64 {
        *self
            .weights_by_percent_gc
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap()
    }

    pub fn mean_weight(&self) -> f64 {
        self.weights_by_percent_gc.iter().sum::<f64>() / self.weights_by_percent_gc.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::nucleotides::Nucleotide::{A, C, G, N, T};

    fn weights_with_value(value: f64) -> Vec<f64> {
        vec![value; 101]
    }

    #[test]
    fn test_default_model_has_101_unit_weights() {
        let model = GcBiasModel::default();

        assert_eq!(model.max_weight(), 1.0);
        assert!((model.mean_weight() - 1.0).abs() < 1e-12);

        for percent_gc in 0..=100 {
            let gc_fraction = percent_gc as f64 / 100.0;
            assert!(
                (model.weight_for_gc_fraction(gc_fraction) - 1.0).abs() < 1e-12,
                "Expected uniform weight 1.0 for {}% GC",
                percent_gc
            );
        }
    }

    #[test]
    fn test_from_weights_rejects_too_few_bins() {
        let result = GcBiasModel::from_weights(vec![1.0; 100], 150);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidBinCount)),
            "Expected InvalidBinCount for 100 weights; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_too_many_bins() {
        let result = GcBiasModel::from_weights(vec![1.0; 102], 150);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidBinCount)),
            "Expected InvalidBinCount for 102 weights; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_negative_weight() {
        let mut weights = weights_with_value(1.0);
        weights[30] = -0.1;

        let result = GcBiasModel::from_weights(weights, 150);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWeight)),
            "Expected InvalidWeight for negative weight; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_nan_weight() {
        let mut weights = weights_with_value(1.0);
        weights[50] = f64::NAN;

        let result = GcBiasModel::from_weights(weights, 150);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWeight)),
            "Expected InvalidWeight for NaN weight; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_infinite_weight() {
        let mut weights = weights_with_value(1.0);
        weights[50] = f64::INFINITY;

        let result = GcBiasModel::from_weights(weights, 150);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWeight)),
            "Expected InvalidWeight for infinite weight; got {:?}",
            result
        );
    }

    #[test]
    fn test_from_weights_rejects_all_zero_weights() {
        let result = GcBiasModel::from_weights(weights_with_value(0.0), 150);

        assert!(
            matches!(result, Err(GcBiasModelError::EmptyModel)),
            "Expected EmptyModel for all-zero weights; got {:?}",
            result
        );
    }

    #[test]
    fn test_weight_for_gc_fraction_uses_percent_bins() {
        let mut weights = weights_with_value(0.0);
        weights[0] = 0.1;
        weights[30] = 0.3;
        weights[50] = 1.0;
        weights[100] = 0.2;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        assert!((model.weight_for_gc_fraction(0.00) - 0.1).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(0.30) - 0.3).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(0.50) - 1.0).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(1.00) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_gc_fraction_rounds_to_nearest_percent() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.5;
        weights[51] = 0.51;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        assert!((model.weight_for_gc_fraction(0.504) - 0.5).abs() < 1e-12);
        assert!((model.weight_for_gc_fraction(0.505) - 0.51).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_gc_fraction_clamps_below_zero() {
        let mut weights = weights_with_value(0.0);
        weights[0] = 0.25;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        assert!((model.weight_for_gc_fraction(-0.10) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_gc_fraction_clamps_above_one() {
        let mut weights = weights_with_value(0.0);
        weights[100] = 0.75;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        assert!((model.weight_for_gc_fraction(1.10) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_zero_percent_gc() {
        let mut weights = weights_with_value(0.0);
        weights[0] = 0.2;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();
        let sequence = vec![A, A, T, T, A, T];

        assert!((model.weight_for_sequence(&sequence) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_fifty_percent_gc() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 1.0;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();
        let sequence = vec![A, C, G, T];

        assert!((model.weight_for_sequence(&sequence) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_one_hundred_percent_gc() {
        let mut weights = weights_with_value(0.0);
        weights[100] = 0.4;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();
        let sequence = vec![C, G, G, C, C, G];

        assert!((model.weight_for_sequence(&sequence) - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_ignores_n_bases_in_gc_fraction() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.8;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        // Ignoring N bases leaves A, C, G, T => 2 GC / 4 called bases = 50% GC.
        let sequence = vec![A, C, G, T, N, N];

        assert!((model.weight_for_sequence(&sequence) - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_empty_slice_returns_neutral() {
        let model = GcBiasModel::default();
        // No bases at all: sequence_count == 0, returns the neutral fallback 1.0
        assert!((model.weight_for_sequence(&[]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_weight_for_sequence_all_n_returns_neutral() {
        let mut weights = weights_with_value(0.0);
        weights[50] = 0.75;
        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        // All-N: sequence_count == 0, returns 1.0 regardless of the model's other weights.
        // N-masked regions are treated as neutral so they are never filtered out.
        let sequence = vec![N, N, N, N];
        assert!((model.weight_for_sequence(&sequence) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_max_weight() {
        let mut weights = weights_with_value(0.1);
        weights[30] = 0.5;
        weights[50] = 1.25;
        weights[70] = 0.75;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        assert!((model.max_weight() - 1.25).abs() < 1e-12);
    }

    #[test]
    fn test_mean_weight() {
        let mut weights = weights_with_value(1.0);
        weights[0] = 0.0;
        weights[100] = 2.0;

        let model = GcBiasModel::from_weights(weights, 150).unwrap();

        // Sum is still 101.0 because one bin changed 1 -> 0 and one changed 1 -> 2.
        assert!((model.mean_weight() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_roundtrip_write_and_read() {
        let mut weights = weights_with_value(0.5);
        weights[30] = 0.3;
        weights[50] = 1.0;
        weights[70] = 0.8;

        let original = GcBiasModel::from_weights(weights, 200).unwrap();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = PathBuf::from(tmp.path());

        original.write_to_file(&path).unwrap();
        let loaded = GcBiasModel::from_file(&path).unwrap();

        assert_eq!(original.is_uniform(), loaded.is_uniform());
        assert_eq!(original.window_size(), loaded.window_size());
        assert!((original.max_weight() - loaded.max_weight()).abs() < 1e-12);
        for pct in 0..=100 {
            let gc = pct as f64 / 100.0;
            assert!(
                (original.weight_for_gc_fraction(gc) - loaded.weight_for_gc_fraction(gc)).abs()
                    < 1e-12,
                "Mismatch at {}% GC after roundtrip",
                pct
            );
        }
    }

    #[test]
    fn test_from_weights_rejects_zero_window_size() {
        let result = GcBiasModel::from_weights(weights_with_value(1.0), 0);

        assert!(
            matches!(result, Err(GcBiasModelError::InvalidWindowSize)),
            "Expected InvalidWindowSize for window_size=0; got {:?}",
            result
        );
    }

    #[test]
    fn test_window_size_accessor_returns_stored_value() {
        let model = GcBiasModel::from_weights(weights_with_value(1.0), 75).unwrap();
        assert_eq!(model.window_size(), 75);
    }

    #[test]
    fn test_default_window_size_is_100() {
        assert_eq!(GcBiasModel::default().window_size(), 100);
    }

    #[test]
    fn test_legacy_file_without_window_size_gets_default() {
        // A model file written before window_size was added should deserialize
        // with window_size == default_window_size() == 150.
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = PathBuf::from(tmp.path());

        // JSON without the window_size field, gzip-compressed to match model_writer format.
        let legacy_json = r#"{"weights_by_percent_gc":[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0],"is_uniform":true}"#;
        let file = std::fs::File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(legacy_json.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let loaded = GcBiasModel::from_file(&path).unwrap();
        assert_eq!(
            loaded.window_size(),
            100,
            "Legacy file missing window_size should default to 100"
        );
    }
}
