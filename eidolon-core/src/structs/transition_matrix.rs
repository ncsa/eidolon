//! The transition matrix extends the notion of a distribution. Specifically, it denotes a transition
//! probability from one base to another. The trinucleotide SNP model requires 4 layers of these transition
//! matrices. Indexing is now implemented for the Tranisition matrix, which should make the code more straightforward.
use crate::structs::{
    distributions::{DiscreteDistribution, DistributionErrors},
    nucleotides::{ALLOWED_NUCS, Nucleotide},
};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::ops::Index;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum TransitionMatrixError {
    #[error("Transition matrix reported a distribution error: {0}")]
    DistributionError(DistributionErrors),
    #[error("Input weights and lengths were of unequal length")]
    UnequalWeightsError,
    #[error("I/O error reading transition matrix file: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Transition matrix file has {0} data rows (expected 4)")]
    InvalidRowCount(usize),
    #[error("Transition matrix file, line {line}: {reason}")]
    InvalidRow { line: usize, reason: String },
}

impl From<DistributionErrors> for TransitionMatrixError {
    fn from(error: DistributionErrors) -> Self {
        TransitionMatrixError::DistributionError(error)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransitionMatrix {
    // Nucleotide transition matrix. Rows represent the base we are mutating and the weights are
    // in the standard nucleotide order (in the same a, c, g, t order). This structure is
    // fundamental to the others and to the mutation model in general.
    //
    // This defines a transition from one Nucleotide to another.
    pub a: DiscreteDistribution<Nucleotide>,
    pub c: DiscreteDistribution<Nucleotide>,
    pub g: DiscreteDistribution<Nucleotide>,
    pub t: DiscreteDistribution<Nucleotide>,
}

impl Index<usize> for TransitionMatrix {
    type Output = DiscreteDistribution<Nucleotide>;
    fn index(&self, i: usize) -> &DiscreteDistribution<Nucleotide> {
        match i {
            0 => &self.a,
            1 => &self.c,
            2 => &self.g,
            3 => &self.t,
            _ => panic!("index out of range: {} with length 4", i),
        }
    }
}

impl Index<&Nucleotide> for TransitionMatrix {
    type Output = DiscreteDistribution<Nucleotide>;
    fn index(&self, n: &Nucleotide) -> &DiscreteDistribution<Nucleotide> {
        match n {
            Nucleotide::A => &self.a,
            Nucleotide::C => &self.c,
            Nucleotide::G => &self.g,
            Nucleotide::T => &self.t,
            _ => panic!("Nucleotide not found: {:?}", n),
        }
    }
}

impl TransitionMatrix {
    // Returns Result because it builds distributions that can fail; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, TransitionMatrixError> {
        // Default transition matrix for mutations from the original NEAT 2.0
        Ok(Self {
            a: DiscreteDistribution::new(
                &vec![
                    0.0,
                    0.16952785544157142,
                    0.6878218371885525,
                    0.1426503073698761,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
            c: DiscreteDistribution::new(
                &vec![
                    0.1615595177675533,
                    0.0,
                    0.1664600510853558,
                    0.6719804311470908,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
            g: DiscreteDistribution::new(
                &vec![
                    0.6704692217289652,
                    0.16706325641014994,
                    0.0,
                    0.1624675218608849,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
            t: DiscreteDistribution::new(
                &vec![
                    0.14281397280584493,
                    0.6885459433549382,
                    0.16864008383921686,
                    0.0,
                ],
                &Vec::from(ALLOWED_NUCS),
            )?,
        })
    }

    /// Load a 4×4 SNP transition matrix from a whitespace-delimited TSV file.
    ///
    /// Rows and columns correspond to A/C/G/T (from-base and to-base).
    /// The first non-blank line is skipped as a header if its first token is non-numeric.
    /// Diagonal values are zeroed so self-transitions are impossible.
    ///
    /// Every data row must hold exactly four finite, non-negative numbers with some
    /// off-diagonal mass, and there must be exactly four rows. Anything else is refused,
    /// naming the line: a row with no off-diagonal mass would otherwise send every draw to
    /// the first base (#760).
    pub fn from_tsv(path: &PathBuf) -> Result<Self, TransitionMatrixError> {
        let content = std::fs::read_to_string(path)?;
        let mut rows: Vec<[f64; 4]> = Vec::new();
        let mut seen_first_line = false;

        for (i, line) in content.lines().enumerate() {
            let line_no = i + 1;
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }
            let is_first = !seen_first_line;
            seen_first_line = true;
            if is_first && tokens[0].parse::<f64>().is_err() {
                continue; // header
            }
            let invalid = |reason: String| TransitionMatrixError::InvalidRow {
                line: line_no,
                reason,
            };
            if tokens.len() != 4 {
                return Err(invalid(format!(
                    "expected 4 values (A C G T), found {}: {line:?}",
                    tokens.len()
                )));
            }
            let mut row = [0.0f64; 4];
            for (slot, token) in row.iter_mut().zip(&tokens) {
                *slot = token
                    .parse::<f64>()
                    .ok()
                    .filter(|v| v.is_finite() && *v >= 0.0)
                    .ok_or_else(|| {
                        invalid(format!(
                            "{token:?} is not a finite, non-negative number: {line:?}"
                        ))
                    })?;
            }
            let r = rows.len();
            if r < 4 {
                row[r] = 0.0; // zero diagonal so self-transitions are impossible
                if row.iter().sum::<f64>() <= 0.0 {
                    return Err(invalid(format!(
                        "the {:?} row has no weight off the diagonal, so it cannot say what \
                         {:?} mutates to: {line:?}",
                        ALLOWED_NUCS[r], ALLOWED_NUCS[r]
                    )));
                }
            }
            rows.push(row);
        }

        if rows.len() != 4 {
            return Err(TransitionMatrixError::InvalidRowCount(rows.len()));
        }

        Self::from(rows[0], rows[1], rows[2], rows[3])
    }

    pub fn from(
        a_weights: [f64; 4],
        c_weights: [f64; 4],
        g_weights: [f64; 4],
        t_weights: [f64; 4],
    ) -> Result<Self, TransitionMatrixError> {
        let weights_test = [a_weights, c_weights, g_weights, t_weights];
        for vector in weights_test.iter() {
            if vector.len() != 4 {
                return Err(TransitionMatrixError::UnequalWeightsError);
            }
        }
        Ok(Self {
            a: DiscreteDistribution::new(&a_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
            c: DiscreteDistribution::new(&c_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
            g: DiscreteDistribution::new(&g_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
            t: DiscreteDistribution::new(&t_weights.to_vec(), &Vec::from(ALLOWED_NUCS))?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::NeatRng;

    #[test]
    fn test_transition_matrix_build() {
        let a_weights = [0.0, 20.0, 1.0, 20.0];
        let c_weights = [20.0, 0.0, 1.0, 1.0];
        let g_weights = [1.0, 1.0, 0.0, 20.0];
        let t_weights = [20.0, 1.0, 20.0, 0.0];
        let model = TransitionMatrix::from(a_weights, c_weights, g_weights, t_weights).unwrap();
        // Index by usize and Nucleotide must reference the same distribution
        assert_eq!(
            model[0].values().unwrap(),
            model[&Nucleotide::A].values().unwrap()
        );
        assert_eq!(
            model[1].values().unwrap(),
            model[&Nucleotide::C].values().unwrap()
        );
        assert_eq!(
            model[2].values().unwrap(),
            model[&Nucleotide::G].values().unwrap()
        );
        assert_eq!(
            model[3].values().unwrap(),
            model[&Nucleotide::T].values().unwrap()
        );
        // Each row's values should be the four ACGT nucleotides
        assert_eq!(
            model[&Nucleotide::A].values().unwrap(),
            Vec::from(ALLOWED_NUCS)
        );
    }

    #[test]
    fn test_from_tsv_valid() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("matrix.tsv");
        std::fs::write(
            &tsv,
            "A\tC\tG\tT\n\
             0.0\t0.5\t0.3\t0.2\n\
             0.5\t0.0\t0.3\t0.2\n\
             0.4\t0.3\t0.0\t0.3\n\
             0.3\t0.3\t0.4\t0.0\n",
        )
        .unwrap();
        let tm = TransitionMatrix::from_tsv(&tsv).unwrap();
        // Check each row can sample one of the four ACGT nucleotides
        let mut rng = NeatRng::new_from_seed(&vec!["t".to_string()]).unwrap();
        for nuc in ALLOWED_NUCS {
            let sample = tm[&nuc].sample(rng.random().unwrap()).unwrap();
            assert!(ALLOWED_NUCS.contains(&sample));
            // Self-transition must be impossible: diagonal was zeroed
            assert_ne!(sample, nuc);
        }
    }

    #[test]
    fn test_from_tsv_too_few_rows_errors() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("short.tsv");
        std::fs::write(&tsv, "0.0\t0.5\t0.3\t0.2\n0.5\t0.0\t0.3\t0.2\n").unwrap();
        assert!(matches!(
            TransitionMatrix::from_tsv(&tsv),
            Err(TransitionMatrixError::InvalidRowCount(2))
        ));
    }

    const TSV_HEADER: &str = "A\tC\tG\tT\n";
    const GOOD_ROWS: [&str; 4] = [
        "0.0\t0.5\t0.3\t0.2",
        "0.5\t0.0\t0.3\t0.2",
        "0.4\t0.3\t0.0\t0.3",
        "0.3\t0.3\t0.4\t0.0",
    ];

    /// Writes the header and four good rows, with `row` (0 = A) replaced by `bad`.
    fn tsv_with_row(dir: &tempfile::TempDir, row: usize, bad: &str) -> PathBuf {
        let mut rows = GOOD_ROWS.map(String::from);
        rows[row] = bad.to_string();
        let path = dir.path().join(format!("row{row}.tsv"));
        std::fs::write(&path, format!("{TSV_HEADER}{}\n", rows.join("\n"))).unwrap();
        path
    }

    // A user-supplied row with no off-diagonal mass used to become "every draw is A" (#760):
    // errors at A vanished and errors at C, G or T all became A. Each malformed row must be
    // refused, and the error must name the line so the user can find it. Line 1 is the
    // header, so row A is line 2 and row C is line 3.
    #[test]
    fn a_malformed_transition_row_is_refused_with_its_line() {
        let dir = tempfile::tempdir().unwrap();
        for (row, bad, line) in [
            (1, "0\t0\t0\t0", 3),            // all zero
            (1, "0\t5\t0\t0", 3),            // only the diagonal, which is ignored
            (0, "0\t0.5\t-0.3\t0.2", 2),     // negative
            (2, "0.4\tnan\t0\t0.3", 4),      // non-finite
            (3, "0.3\t0.3\t0.4", 5),         // three values
            (3, "0.3\t0.3\t0.4\t0\t0.1", 5), // five values
            (0, "0\t0.5\tx\t0.2", 2),        // non-numeric
        ] {
            let path = tsv_with_row(&dir, row, bad);
            match TransitionMatrix::from_tsv(&path) {
                Err(TransitionMatrixError::InvalidRow { line: got, .. }) => {
                    assert_eq!(got, line, "row {bad:?} reported at the wrong line")
                }
                other => panic!("row {bad:?} should be refused at line {line}, got {other:?}"),
            }
        }
    }

    #[test]
    fn a_transition_file_with_extra_rows_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("five.tsv");
        std::fs::write(
            &path,
            format!("{TSV_HEADER}{}\n{}\n", GOOD_ROWS.join("\n"), GOOD_ROWS[0]),
        )
        .unwrap();
        assert!(matches!(
            TransitionMatrix::from_tsv(&path),
            Err(TransitionMatrixError::InvalidRowCount(5))
        ));
    }

    // Must not fire: a good file loads with each row's off-diagonal weights exactly as given.
    #[test]
    fn a_good_transition_file_loads_as_given() {
        let dir = tempfile::tempdir().unwrap();
        let path = tsv_with_row(&dir, 0, GOOD_ROWS[0]);
        let tm = TransitionMatrix::from_tsv(&path).unwrap();
        assert_eq!(tm.c.weights().unwrap(), vec![0.5, 0.5, 0.8, 1.0]);
    }

    #[test]
    fn test_transition_matrix_default() {
        let model = TransitionMatrix::default().unwrap();
        // Spot-check: sampling from A row should return one of [A, C, G, T]
        let mut rng = NeatRng::new_from_seed(&vec!["seed".to_string()]).unwrap();
        let sample = model[&Nucleotide::A].sample(rng.random().unwrap()).unwrap();
        assert!(ALLOWED_NUCS.contains(&sample));
    }
}
