//! IUPAC nucleotide ambiguity codes

use crate::codec;
use crate::codec::dna::Dna;
use crate::codec::{Codec, ParseBioErr};
use std::fmt;
use std::str::FromStr;

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq)]
//#[width(4)]
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

impl Codec for Iupac {
    const WIDTH: u8 = 2;
    fn unsafe_from_bits(b: u8) -> Self {
        match b {
            _ => Self::A,
        }
    }
    fn try_from_bits(b: u8) -> Result<Self, codec::ParseBioErr> {
        match b {
            _ => Err(ParseBioErr),
        }
    }
    fn from_char(c: char) -> Result<Self, codec::ParseBioErr> {
        Ok(Self::A)
    }
    fn to_char(a: Self) -> char {
        match a {
            _ => 'A',
        }
    }
}

impl TryFrom<char> for Iupac {
    type Error = ParseBioErr;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' => Ok(Iupac::A),
            'C' => Ok(Iupac::C),
            'G' => Ok(Iupac::G),
            'T' => Ok(Iupac::T),
            'R' => Ok(Iupac::R),
            'Y' => Ok(Iupac::Y),
            'S' => Ok(Iupac::S),
            'W' => Ok(Iupac::W),
            'K' => Ok(Iupac::K),
            'M' => Ok(Iupac::M),
            'B' => Ok(Iupac::B),
            'D' => Ok(Iupac::D),
            'H' => Ok(Iupac::H),
            'V' => Ok(Iupac::V),
            'N' => Ok(Iupac::N),
            '-' => Ok(Iupac::X),
            '.' => Ok(Iupac::X),
            _ => Err(ParseBioErr),
        }
    }
}

impl From<Iupac> for char {
    fn from(iupac: Iupac) -> char {
        match iupac {
            Iupac::A => 'A',
            Iupac::C => 'C',
            Iupac::G => 'G',
            Iupac::T => 'T',
            Iupac::R => 'R',
            Iupac::Y => 'Y',
            Iupac::S => 'S',
            Iupac::W => 'W',
            Iupac::K => 'K',
            Iupac::M => 'M',
            Iupac::B => 'B',
            Iupac::D => 'D',
            Iupac::H => 'H',
            Iupac::V => 'V',
            Iupac::N => 'N',
            Iupac::X => '-',
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

impl From<u8> for Iupac {
    fn from(b: u8) -> Self {
        match b {
            0b1000 => Iupac::A,
            0b0100 => Iupac::C,
            0b0010 => Iupac::G,
            0b0001 => Iupac::T,
            0b1010 => Iupac::R,
            0b0101 => Iupac::Y,
            0b0110 => Iupac::S,
            0b1001 => Iupac::W,
            0b0011 => Iupac::K,
            0b1100 => Iupac::M,
            0b0111 => Iupac::B,
            0b1011 => Iupac::D,
            0b1101 => Iupac::H,
            0b1110 => Iupac::V,
            0b1111 => Iupac::N,
            0b0000 => Iupac::X,
            _ => Iupac::X,
        }
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
