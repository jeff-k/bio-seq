//! Amino acids
/// TODO Meant to be used with the `alt_chars=[...]` and `alt_bytes=[...]` attributes.
///
/// With the 6-bit representation a `Seq<Dna>` sequence can be read as an amino acid
/// sequence in place.
use crate::codec::{Codec, ParseBioErr};
use core::fmt;
//use core::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Codec)]
#[width = 6]
#[repr(u8)]
pub enum Amino {
    A = 0b000110, // GCA
    C = 0b011011, // TGC
    D = 0b010010, // GAC
    E = 0b000010, // GAA
    F = 0b011111, // TTC
    G = 0b001010, // GGA
    H = 0b010001, // CAC
    I = 0b001100, // ATA
    K = 0b000000, // AAA
    L = 0b001101, // CTA
    M = 0b101100, // ATG
    N = 0b010000, // AAC
    P = 0b000101, // CCA
    Q = 0b000001, // CAA
    R = 0b001000, // AGA
    S = 0b011000, // AGC
    T = 0b000100, // ACA
    V = 0b001110, // GTA
    W = 0b101011, // TGG
    Y = 0b010011, // TAC
    X = 0b000011, // TAA (stop)
}

impl From<Amino> for u8 {
    fn from(amino: Amino) -> Self {
        amino as u8
    }
}

impl fmt::Display for Amino {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[macro_export]
macro_rules! amino {
    ($seq:expr) => {
        match Seq::<Amino>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}
