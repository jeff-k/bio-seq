//! Bit packable enums representing biological alphabets

pub mod dna;
pub mod iupac;

use std::fmt;
use std::str::FromStr;

pub trait Codec: FromStr + fmt::Display + fmt::Debug + From<u8> + Into<u8> + Copy + Clone {
    const WIDTH: u8;
    //    fn to_bits(&self) -> u8;
    //    fn from_bits(b: &u8) -> Self;
    fn from_char(c: &char) -> Result<Self, ParseBioErr>;
    fn to_char(&self) -> char;
}

#[derive(Debug, Clone)]
pub struct ParseBioErr;

impl fmt::Display for ParseBioErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Could not encode")
    }
}

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
