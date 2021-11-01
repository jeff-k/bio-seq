//! Bit packable enums representing biological alphabets

pub mod amino;
pub mod dna;
pub mod iupac;

use bitvec::prelude::*;
use std::fmt;
use std::str::FromStr;

pub trait Alphabet: FromStr + fmt::Display + fmt::Debug {
    fn width() -> usize;
    fn to_bits(&self) -> BitVec;
    //    fn from_u8(&self, b: u8) -> Self;
    fn from_bits(b: &BitSlice) -> Self;
}

pub trait Complement {
    fn complement(base: Self) -> Self;
}

#[derive(Debug, Clone)]
pub struct ParseBioErr;

impl fmt::Display for ParseBioErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid biological sequence")
    }
}
