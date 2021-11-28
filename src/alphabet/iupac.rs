//! IUPAC nucleotide ambiguity codes

use std::fmt;
use std::str::FromStr;

use crate::alphabet::dna::Dna;
use crate::alphabet::{Alphabet, ParseBioErr};
use bitvec::prelude::*;

#[derive(Debug, PartialEq)]
pub enum Iupac {
    A = 0b0001,
    C = 0b0010,
    G = 0b0100,
    T = 0b1000,
    R = 0b0101,
    Y = 0b1010,
    S = 0b0110,
    W = 0b1001,
    K = 0b1100,
    M = 0b0011,
    B = 0b1110,
    D = 0b1101,
    H = 0b1011,
    V = 0b0111,
    N = 0b1111,
    X = 0b0000,
}

impl Alphabet for Iupac {
    const WIDTH: usize = 4;

    fn to_bits(&self) -> BitVec {
        match &self {
            Iupac::A => bitvec![1, 0, 0, 0],
            Iupac::C => bitvec![0, 1, 0, 0],
            Iupac::G => bitvec![0, 0, 1, 0],
            Iupac::T => bitvec![0, 0, 0, 1],

            Iupac::R => bitvec![1, 0, 1, 0],
            Iupac::Y => bitvec![0, 1, 0, 1],
            Iupac::S => bitvec![0, 1, 1, 0],
            Iupac::W => bitvec![1, 0, 0, 1],

            Iupac::K => bitvec![0, 0, 1, 1],
            Iupac::M => bitvec![1, 1, 0, 0],
            Iupac::B => bitvec![0, 1, 1, 1],
            Iupac::D => bitvec![1, 0, 1, 1],

            Iupac::H => bitvec![1, 1, 0, 1],
            Iupac::V => bitvec![1, 1, 1, 0],
            Iupac::N => bitvec![1, 1, 1, 1],
            Iupac::X => bitvec![0, 0, 0, 0],
        }
    }

    fn from_bits(b: &BitSlice) -> Self {
        let bs: [bool; 4] = [b[0], b[1], b[2], b[3]];
        match bs {
            [true, false, false, false] => Iupac::A,
            [false, true, false, false] => Iupac::C,
            [false, false, true, false] => Iupac::G,
            [false, false, false, true] => Iupac::T,

            [true, false, true, false] => Iupac::R,
            [false, true, false, true] => Iupac::Y,
            [false, true, true, false] => Iupac::S,
            [true, false, false, true] => Iupac::W,

            [false, false, true, true] => Iupac::K,
            [true, true, false, false] => Iupac::M,
            [false, true, true, true] => Iupac::B,
            [true, false, true, true] => Iupac::D,

            [true, true, false, true] => Iupac::H,
            [true, true, true, false] => Iupac::V,
            [true, true, true, true] => Iupac::N,
            [false, false, false, false] => Iupac::X,
        }
    }
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
