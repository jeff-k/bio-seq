//! IUPAC nucleotide ambiguity codes

use std::fmt;
use std::str::FromStr;

use crate::codec::dna::Dna;
use crate::codec::{Codec, ParseBioErr};
use bitvec::prelude::*;

#[derive(Debug, PartialEq, Codec)]
pub enum Iupac {
    A = 0b1000,
    C = 0b0100,
    G = 0b0010,
    T = 0b0001,
    R = 0b1010,
    Y = 0b0000,
    S = 0b0000,
    W = 0b0000,
    K = 0b0000,
    M = 0b0000,
    B = 0b0000,
    D = 0b0000,
    H = 0b0000,
    V = 0b0000,
    N = 0b1111,
    X = 0b0000,
}

impl From<Dna> for Iupac {
    fn from(dna: Dna) -> Self {
        match dna {
            Dna::A => Iupac::A,
            Dna::C => Iupac::C,
            Dna::G => Iupac::G,
            Dna::T => Iupac::T,
        }
    }
}

impl FromStr for Iupac {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Iupac::A),
            "C" => Ok(Iupac::C),
            "G" => Ok(Iupac::G),
            "T" => Ok(Iupac::T),
            "R" => Ok(Iupac::R),
            "Y" => Ok(Iupac::Y),
            "S" => Ok(Iupac::S),
            "W" => Ok(Iupac::W),
            "K" => Ok(Iupac::K),
            "M" => Ok(Iupac::M),
            "B" => Ok(Iupac::B),
            "D" => Ok(Iupac::D),
            "H" => Ok(Iupac::H),
            "V" => Ok(Iupac::V),
            "N" => Ok(Iupac::N),
            "." => Ok(Iupac::X),
            "-" => Ok(Iupac::X),
            _ => Err(ParseBioErr),
        }
    }
}

impl fmt::Display for Iupac {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Iupac::X => write!(f, "-"),
            _ => write!(f, "{:?}", self),
        }
    }
}
