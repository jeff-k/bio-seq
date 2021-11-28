//! Bit packable enums representing biological alphabets

pub mod amino;
pub mod dna;
pub mod iupac;

use bitvec::prelude::*;
use std::fmt;
use std::str::FromStr;

pub trait Alphabet: FromStr + fmt::Display + fmt::Debug {
    const WIDTH: usize;
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

#[cfg(test)]
mod tests {
    use super::amino::Amino;
    use super::dna::Dna;
    use super::iupac::Iupac;
    use crate::dna;
    use crate::kmer::Kmer;
    use crate::seq::Seq;
    use std::str::FromStr;

    #[test]
    fn dna_to_iupac() {
        assert_eq!(Iupac::from(Dna::A), Iupac::A);
        assert_ne!(Iupac::from(Dna::A), Iupac::T);
        assert_eq!(Iupac::from(Dna::C), Iupac::C);
    }

    #[test]
    fn codon_to_amino() {
        assert_eq!(Amino::from(Kmer::<3>::new(&dna!("AAA").bv)), Amino::K);
    }
}
