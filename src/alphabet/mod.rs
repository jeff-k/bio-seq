pub mod dna;

use bitvec::prelude::*;

pub trait Alphabet {
    fn width() -> usize;
    fn to_bits(&self) -> BitVec;
    //    fn from_u8(&self, b: u8) -> Self;
    fn from_bits(b: usize) -> Self;
}

pub trait Complement {
    fn complement(base: Self) -> Self;
}
