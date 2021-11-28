//! IUPAC nucleotide ambiguity codes

use std::fmt;
use std::str::FromStr;

use crate::codec::dna::Dna;
use crate::codec::{Codec, ParseBioErr};
use bitvec::prelude::*;

#[derive(Debug, PartialEq)]
pub enum Iupac {
    A, C, G, T, R, Y, S, W, K, M, B, D, H,V, N, X,
}

impl Codec for Iupac {
    const WIDTH: usize = 4;

    fn to_bits(&self) -> BitArray::<Msb0, u8> {
        match &self {
            Iupac::A => BitArray::new(0b1000),
            Iupac::C => BitArray::new(0b0100),
            Iupac::G => BitArray::new(0b0010),
            Iupac::T => BitArray::new(0b0001),

            Iupac::R => BitArray::new(0b1000),
            Iupac::Y => BitArray::new(0b0100),
            Iupac::S => BitArray::new(0b0010),
            Iupac::W => BitArray::new(0b0001),

            Iupac::K => BitArray::new(0b1000),
            Iupac::M => BitArray::new(0b0100),
            Iupac::B => BitArray::new(0b0010),
            Iupac::D => BitArray::new(0b0001),

            Iupac::H => BitArray::new(0b1000),
            Iupac::V => BitArray::new(0b0100),
            Iupac::N => BitArray::new(0b0010),
            Iupac::X => BitArray::new(0b0001),
        }
    }

    fn from_bits(b: &BitSlice::<Msb0, u8>) -> Self {
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

    fn to_ascii(c: u8) -> Self {
        unimplemented!()
    }

    fn from_ascii(c: u8) -> Self {
        unimplemented!()
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
