/// 2-bit DNA representation
///
use std::fmt;
use std::str::FromStr;

use crate::codec::{Codec, ParseBioErr};
use crate::kmer::Kmer;
use crate::Complement;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Codec)]
#[width = 2]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl<const K: usize> Complement for Kmer<Dna, K> {
    fn comp(self: Kmer<Dna, K>) -> Kmer<Dna, K> {
        Kmer::from(self.bs ^ usize::MAX)
    }
}

impl From<Dna> for char {
    fn from(dna: Dna) -> char {
        dna.to_char()
    }
}

impl TryFrom<char> for Dna {
    type Error = ParseBioErr;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        Dna::from_char(c)
    }
}

impl From<u8> for Dna {
    fn from(b: u8) -> Self {
        Dna::unsafe_from_bits(b)
    }
}

impl From<Dna> for u8 {
    fn from(dna: Dna) -> Self {
        dna as u8
    }
}

impl FromStr for Dna {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Dna::try_from(s.as_bytes()[0] as char)
    }
}

impl fmt::Display for Dna {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[macro_export]
macro_rules! dna {
    ($seq:expr) => {
        match Seq::<Dna>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::*;
    use crate::kmer::Kmer;
    use crate::{codec::ParseBioErr, Complement, Seq};
    use std::convert::TryFrom;
    use std::str::FromStr;

    #[test]
    fn dna_kmer_equality() -> Result<(), ParseBioErr> {
        assert_eq!(
            Kmer::<Dna, 8>::try_from(dna!("TGCACATG"))?,
            Kmer::<Dna, 8>::try_from(dna!("TGCACATG"))?
        );
        assert_ne!(
            Kmer::<Dna, 7>::try_from(dna!("GTGACGA"))?,
            Kmer::<Dna, 7>::try_from(dna!("GTGAAGA"))?
        );
        Ok(())
    }

    #[test]
    fn dna_kmer_complement() -> Result<(), ParseBioErr> {
        assert_eq!(
            Kmer::<Dna, 16>::try_from(dna!("AAAATGCACATGTTTT"))?.comp(),
            Kmer::<Dna, 16>::try_from(dna!("TTTTACGTGTACAAAA"))?
        );

        assert_ne!(
            Kmer::<Dna, 8>::try_from(dna!("TGCACATG"))?,
            Kmer::<Dna, 8>::try_from(dna!("ACGTGTAC"))?
        );
        Ok(())
    }
}
