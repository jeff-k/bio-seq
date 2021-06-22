pub mod dna;
pub mod iupac;

use bitvec::prelude::*;

pub trait Alphabet {
    fn width() -> usize;
    fn to_bits(&self) -> BitVec;
    //    fn from_u8(&self, b: u8) -> Self;
    fn from_bits(b: &BitSlice) -> Self;
}

pub trait Complement {
    fn complement(base: Self) -> Self;
}
