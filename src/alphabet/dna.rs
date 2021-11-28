use crate::alphabet::{Alphabet, ParseBioErr};
use bitvec::prelude::*;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl Alphabet for Dna {
    const WIDTH: usize = 2;

    fn to_bits(&self) -> BitVec {
        match &self {
            Dna::A => bitvec![0, 0],
            Dna::C => bitvec![0, 1],
            Dna::G => bitvec![1, 0],
            Dna::T => bitvec![1, 1],
        }
    }

    fn from_bits(b: &BitSlice) -> Self {
        let bs: [bool; 2] = [b[0], b[1]];
        match bs {
            [false, false] => Dna::A,
            [false, true] => Dna::C,
            [true, false] => Dna::G,
            [true, true] => Dna::T,
        }
    }
}

impl FromStr for Dna {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Dna::A),
            "C" => Ok(Dna::C),
            "G" => Ok(Dna::G),
            "T" => Ok(Dna::T),
            _ => Err(ParseBioErr),
        }
    }
}

impl fmt::Display for Dna {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
