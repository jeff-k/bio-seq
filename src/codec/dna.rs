use crate::codec::{Codec, ParseBioErr};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl Codec for Dna {
    const WIDTH: usize = 2;

    fn to_bits(&self) -> u8 {
        *self as u8
    }

    fn from_bits(b: &u8) -> Self {
        match b {
            0b00 => Dna::A,
            0b01 => Dna::C,
            0b10 => Dna::G,
            0b11 => Dna::T,
            _ => Dna::A,
        }
    }

    fn to_char(&self) -> char {
        match &self {
            Dna::A => 'A',
            Dna::C => 'C',
            Dna::G => 'G',
            Dna::T => 'T',
        }
    }

    fn from_char(c: &char) -> Result<Self, ParseBioErr> {
        match c {
            'A' => Ok(Dna::A),
            'C' => Ok(Dna::C),
            'G' => Ok(Dna::G),
            'T' => Ok(Dna::T),
            _ => Err(ParseBioErr),
        }
    }
}

impl FromStr for Dna {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Dna::from_char(&(s.as_bytes()[0] as char))
    }
}

impl fmt::Display for Dna {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
