//! IUPAC nucleotide ambiguity codes

use crate::codec::dna::Dna;
use crate::codec::{Codec, ParseBioErr};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Codec)]
#[width = 4]
#[repr(u8)]
pub enum Iupac {
    A = 0b1000,
    C = 0b0100,
    G = 0b0010,
    T = 0b0001,
    R = 0b1010,
    Y = 0b0101,
    S = 0b0110,
    W = 0b1001,
    K = 0b0011,
    M = 0b1100,
    B = 0b0111,
    D = 0b1011,
    H = 0b1101,
    V = 0b1110,
    N = 0b1111,
    X = 0b0000,
}

impl TryFrom<char> for Iupac {
    type Error = ParseBioErr;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        Iupac::from_char(c)
    }
}

impl From<Iupac> for char {
    fn from(iupac: Iupac) -> char {
        iupac.to_char()
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

impl From<u8> for Iupac {
    fn from(b: u8) -> Self {
        Iupac::unsafe_from_bits(b)
    }
}

impl From<Iupac> for u8 {
    fn from(iupac: Iupac) -> Self {
        iupac as u8
    }
}

impl FromStr for Iupac {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Iupac::try_from(s.as_bytes()[0] as char)
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
