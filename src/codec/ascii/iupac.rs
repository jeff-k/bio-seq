//! IUPAC nucleotide ambiguity codes

use crate::codec::{Codec, Dna, ParseBioErr};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Codec)]
#[width = 8]
#[repr(u8)]
pub enum Dna {
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
    R = b'R',
    Y = b'Y',
    S = b'S',
    W = b'W',
    K = b'K',
    M = b'M',
    B = b'B',
    D = b'D',
    H = b'H',
    V = b'V',
    N = b'N',
    #[alt = '-']
    X = b'-',
}

impl TryFrom<char> for Iupac {
    type Error = ParseBioErr;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        Dna::from_char(c)
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
