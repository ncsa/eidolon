//! Representing bases by a simple u8 num
//!     A = 0
//!     C = 1
//!     G = 2
//!     T = 3
//!     N (or other unknown chars) = 4
//! This is intended to make it easier to store them. Note that the bases are always in alphabetical order
//! with N tacked on at the end.
use core::fmt;

use crate::structs::distributions::DiscreteDistribution;
use serde::{Deserialize, Serialize};

// The following are equivalent
pub const ALLOWED_NUCS: [Nucleotide; 4] =
    [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
pub const ALLOWED_USIZE: [usize; 4] = [0, 1, 2, 3];

pub fn allowed_vec() -> Vec<Nucleotide> {
    // This gives us the vec form for iteration and copying into models
    Vec::from(ALLOWED_NUCS)
}

pub fn allowed_usize() -> Vec<usize> {
    // This makes the distribution part work easier, so that can stay more generic
    Vec::from(ALLOWED_USIZE)
}

/// True when every byte of `s` is one of A/C/G/T/N (any case). Useful for
/// distinguishing literal-base VCF REF/ALT fields from symbolic / structural
/// alleles (e.g. `<DUP:TANDEM>`, `<DEL>`, breakends like `G]17:198982]`).
pub fn is_acgtn(s: &str) -> bool {
    s.bytes().all(|b| {
        matches!(
            b,
            b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
        )
    })
}

pub struct NucleotideSelector {
    // Struct for selecting a random nucleotide
    distribution: DiscreteDistribution<Nucleotide>,
}

impl Default for NucleotideSelector {
    fn default() -> Self {
        Self::new()
    }
}

impl NucleotideSelector {
    pub fn new() -> Self {
        let allowed_nucs: Vec<Nucleotide> = allowed_vec();
        let weights: Vec<f64> = vec![0.25, 0.25, 0.25, 0.25];
        NucleotideSelector {
            distribution: DiscreteDistribution::new(&weights, &allowed_nucs)
                .expect("Error creating distribution"),
        }
    }

    pub fn sample_bases(&self, rand: f64) -> Nucleotide {
        self.distribution
            .sample(rand)
            .expect("Error sampling bases")
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Serialize, Deserialize)]
pub enum Nucleotide {
    // A, C, G, T is the canonical order inherited from NEAT2, and every table indexed
    // by base assumes it. The discriminants match `ALLOWED_NUCS` and the `usize`
    // conversions below, so `as usize` and `usize::from` agree for every base. X is the
    // one exception: it is buffer filler, not a base, and `usize::from(X)` is 4.
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
    Maskeda = 5,
    Maskedc = 6,
    Maskedg = 7,
    Maskedt = 8,
    X = 9, // This is purely used to fill out buffers when writing files.
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let name: char = self.to_owned().into();
        write!(f, "{}", name)
    }
}

impl From<char> for Nucleotide {
    fn from(c: char) -> Self {
        match c {
            'A' => Self::A,
            'C' => Self::C,
            'G' => Self::G,
            'T' => Self::T,
            'a' => Self::Maskeda,
            'c' => Self::Maskedc,
            'g' => Self::Maskedg,
            't' => Self::Maskedt,
            _ => Self::N,
        }
    }
}

impl From<usize> for Nucleotide {
    fn from(value: usize) -> Self {
        match value {
            0 => Self::A,
            1 => Self::C,
            2 => Self::G,
            3 => Self::T,
            5 => Self::Maskeda,
            6 => Self::Maskedc,
            7 => Self::Maskedg,
            8 => Self::Maskedt,
            _ => Self::N,
        }
    }
}

impl From<Nucleotide> for usize {
    fn from(val: Nucleotide) -> Self {
        match val {
            Nucleotide::A => 0,
            Nucleotide::C => 1,
            Nucleotide::G => 2,
            Nucleotide::T => 3,
            Nucleotide::Maskeda => 5,
            Nucleotide::Maskedc => 6,
            Nucleotide::Maskedg => 7,
            Nucleotide::Maskedt => 8,
            _ => 4,
        }
    }
}

impl From<Nucleotide> for char {
    fn from(val: Nucleotide) -> Self {
        match val {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'T',
            Nucleotide::Maskeda => 'a',
            Nucleotide::Maskedc => 'c',
            Nucleotide::Maskedg => 'g',
            Nucleotide::Maskedt => 't',
            _ => 'N',
        }
    }
}

impl Nucleotide {
    pub fn complement(&self) -> Self {
        // matches with the complement of each nucleotide.
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::Maskeda => Self::Maskedt,
            Self::Maskedc => Self::Maskedg,
            Self::Maskedg => Self::Maskedc,
            Self::Maskedt => Self::Maskeda,
            _ => *self,
        }
    }

    pub fn get_unmasked_base(&self) -> Self {
        match self {
            Self::Maskeda => Nucleotide::A,
            Self::Maskedc => Nucleotide::C,
            Self::Maskedg => Nucleotide::G,
            Self::Maskedt => Nucleotide::T,
            Self::X => Self::N,
            _ => *self,
        }
    }

    pub fn is_masked(&self) -> bool {
        matches!(
            self,
            Self::Maskeda | Self::Maskedc | Self::Maskedg | Self::Maskedt | Self::X
        )
    }

    pub fn get_masked(&self) -> Nucleotide {
        match self {
            Self::A => Self::Maskeda,
            Self::C => Self::Maskedc,
            Self::G => Self::Maskedg,
            Self::T => Self::Maskedt,
            Self::N => Self::X,
            _ => *self,
        }
    }
}

#[allow(unused)]
pub enum AminoAcid {
    A = 255,
    R = 254,
    N = 253,
    D = 252,
    C = 251,
    E = 250,
    Q = 249,
    G = 248,
    H = 247,
    I = 246,
    L = 245,
    K = 244,
    M = 243,
    F = 242,
    P = 241,
    S = 240,
    T = 239,
    W = 238,
    Y = 237,
    V = 236,
    X = 235,
}

impl From<&str> for AminoAcid {
    fn from(name: &str) -> AminoAcid {
        // Translates string representations of the amino acids into char
        match name.to_lowercase().as_str() {
            "alanine" | "ala" => Self::A,
            "arginine" | "arg" => Self::R,
            "asparagine" | "asn" => Self::N,
            "aspartic acid" | "asp" => Self::D,
            "cysteine" | "cys" => Self::C,
            "glutamic acid" | "glu" => Self::E,
            "glutamine" | "gln" => Self::Q,
            "glycine" | "gly" => Self::G,
            "histidine" | "his" => Self::H,
            "isoleucine" | "ile" => Self::I,
            "leucine" | "leu" => Self::L,
            "lysine" | "lys" => Self::K,
            "methionine" | "met" => Self::M,
            "phenylaline" | "phe" => Self::F,
            "proline" | "pro" => Self::P,
            "serine" | "ser" => Self::S,
            "threonine" | "thr" => Self::T,
            "tryptophan" | "trp" => Self::W,
            "tyrosine" | "tyr" => Self::Y,
            "valine" | "val" => Self::V,
            _ => panic!("Unknown amino acid: {}", name),
        }
    }
}

pub fn sequence_array_to_string(input_array: &[Nucleotide]) -> String {
    // Converts a sequence vector into a string representing the DNA sequence
    let mut return_string = String::with_capacity(input_array.len());
    for nuc in input_array {
        // We need the nucleotide the pointer is aimed at
        return_string.push((*nuc).into());
    }
    return_string
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discriminants_follow_the_canonical_acgt_order() {
        for (slot, &base) in ALLOWED_NUCS.iter().enumerate() {
            assert_eq!(base as usize, slot, "{base} discriminant");
            assert_eq!(usize::from(base), slot, "usize::from({base})");
            assert_eq!(Nucleotide::from(slot), base, "Nucleotide::from({slot})");
        }
        for base in [
            Nucleotide::N,
            Nucleotide::Maskeda,
            Nucleotide::Maskedc,
            Nucleotide::Maskedg,
            Nucleotide::Maskedt,
        ] {
            assert_eq!(base as usize, usize::from(base), "{base:?}");
        }
    }

    #[test]
    fn test_get_unmasked_base() {
        assert_eq!(Nucleotide::Maskeda.get_unmasked_base(), Nucleotide::A);
        assert_eq!(Nucleotide::Maskedc.get_unmasked_base(), Nucleotide::C);
        assert_eq!(Nucleotide::Maskedg.get_unmasked_base(), Nucleotide::G);
        assert_eq!(Nucleotide::Maskedt.get_unmasked_base(), Nucleotide::T);
        assert_eq!(Nucleotide::X.get_unmasked_base(), Nucleotide::N);
        assert_eq!(Nucleotide::A.get_unmasked_base(), Nucleotide::A);
    }

    #[test]
    fn test_is_masked() {
        assert!(Nucleotide::Maskeda.is_masked());
        assert!(Nucleotide::Maskedc.is_masked());
        assert!(Nucleotide::Maskedg.is_masked());
        assert!(Nucleotide::Maskedt.is_masked());
        assert!(Nucleotide::X.is_masked());
        assert!(!Nucleotide::A.is_masked());
        assert!(!Nucleotide::C.is_masked());
        assert!(!Nucleotide::G.is_masked());
        assert!(!Nucleotide::T.is_masked());
        assert!(!Nucleotide::N.is_masked());
    }

    #[test]
    fn test_get_masked() {
        assert_eq!(Nucleotide::A.get_masked(), Nucleotide::Maskeda);
        assert_eq!(Nucleotide::C.get_masked(), Nucleotide::Maskedc);
        assert_eq!(Nucleotide::G.get_masked(), Nucleotide::Maskedg);
        assert_eq!(Nucleotide::T.get_masked(), Nucleotide::Maskedt);
        assert_eq!(Nucleotide::N.get_masked(), Nucleotide::X);
    }

    #[test]
    fn test_sequence_array_to_string() {
        let seq = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        assert_eq!(sequence_array_to_string(&seq), "ACGT");
        let masked = vec![Nucleotide::Maskeda, Nucleotide::Maskedc];
        assert_eq!(sequence_array_to_string(&masked), "ac");
        assert_eq!(sequence_array_to_string(&vec![]), "");
    }

    #[test]
    fn test_base_to_nuc() {
        let test_nuc: char = 'A';
        let nuc: Nucleotide = Nucleotide::from(test_nuc);
        let nuc_num: usize = nuc.into();
        assert_eq!(nuc_num, 0);
    }

    #[test]
    fn cast_type() {
        let nuc = Nucleotide::from(0);
        let nuc_char: char = nuc.into();
        assert_eq!(nuc_char, 'A')
    }

    #[test]
    fn test_complement() {
        let nuc1 = Nucleotide::A;
        let nuc2 = Nucleotide::C;
        let nuc3 = Nucleotide::G;
        let nuc4 = Nucleotide::T;
        let nuc5 = Nucleotide::N;

        assert_eq!(nuc1.complement(), nuc4);
        assert_eq!(nuc2.complement(), nuc3);
        assert_eq!(nuc3.complement(), nuc2);
        assert_eq!(nuc4.complement(), nuc1);
        assert_eq!(nuc5.complement(), nuc5);
    }
}
