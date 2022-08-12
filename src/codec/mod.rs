//! Coding/Decoding trait for bit-packable enums representing biological alphabets

#[macro_use]
pub mod amino;
#[macro_use]
pub mod dna;
#[macro_use]
pub mod iupac;

use std::fmt;

use bio_seq_derive::Codec;

pub trait Codec: Copy + Clone + Into<u8> {
    type Error;
    const WIDTH: u8;

    fn unsafe_from_bits(b: u8) -> Self;
    fn try_from_bits(b: u8) -> Result<Self, Self::Error>;
    fn from_char(c: char) -> Result<Self, Self::Error>;
    fn to_char(self) -> char;
}

pub trait Complement {
    fn comp(self) -> Self;
}

pub trait ReverseComplement {
    /// Reverse complementable
    fn revcomp(self) -> Self;
}

#[derive(Debug, Clone)]
pub struct ParseBioErr;

impl fmt::Display for ParseBioErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Could not encode")
    }
}

/*impl<A: Codec> dyn From<&SeqSlice<A>> {
    fn from(slice: &SeqSlice<A>) {
        A::unsafe_from_bits(slice.load())
    }
}*/

#[cfg(test)]
mod tests {
    use super::dna::Dna;
    use super::iupac::Iupac;

    #[test]
    fn dna_to_iupac() {
        assert_eq!(Iupac::from(Dna::A), Iupac::A);
        assert_ne!(Iupac::from(Dna::A), Iupac::T);
        assert_eq!(Iupac::from(Dna::C), Iupac::C);
    }
}
