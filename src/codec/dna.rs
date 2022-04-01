use crate::codec::{Codec, ParseBioErr};
use bitvec::prelude::*;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq, Codec)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl FromStr for Dna {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Dna::from_char(&s.as_bytes()[0])
    }
}

impl fmt::Display for Dna {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
