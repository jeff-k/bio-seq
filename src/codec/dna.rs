/// 2-bit DNA representation
///
use core::fmt;
use core::marker::PhantomData;
use core::str::FromStr;

use crate::codec::{Codec, Complement};
use crate::{kmer::Kmer, ParseBioError};

use bitvec::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Codec)]
#[width = 2]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl<const K: usize> Complement for Kmer<Dna, K> {
    fn comp(self: Kmer<Dna, K>) -> Self {
        Kmer {
            _p: PhantomData,
            bs: (self.bs ^ usize::MAX).view_bits::<Lsb0>()[..K * Dna::WIDTH as usize]
                .load_le::<usize>(),
        }
    }
}

impl Complement for Dna {
    fn comp(self: Dna) -> Dna {
        match self {
            Dna::A => Dna::T,
            Dna::C => Dna::G,
            Dna::G => Dna::C,
            Dna::T => Dna::A,
        }
    }
}

impl From<Dna> for char {
    fn from(dna: Dna) -> char {
        dna.to_char()
    }
}

impl TryFrom<char> for Dna {
    type Error = ParseBioError;
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
    type Err = ParseBioError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Dna::try_from(s.as_bytes()[0] as char)
    }
}

impl fmt::Display for Dna {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self:?}")
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

#[macro_export]
macro_rules! kmer {
    ($seq:expr) => {
        match Seq::<Dna>::from_str($seq) {
            Ok(s) => match Kmer::<Dna, { $seq.len() }>::try_from(s) {
                Ok(s) => s,
                Err(_) => panic!(),
            },
            Err(_) => panic!(),
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn dna_kmer_equality() {
        assert_eq!(
            Kmer::<Dna, 8>::try_from(dna!("TGCACATG")).unwrap(),
            Kmer::<Dna, 8>::try_from(dna!("TGCACATG")).unwrap()
        );
        assert_ne!(
            Kmer::<Dna, 7>::try_from(dna!("GTGACGA")).unwrap(),
            Kmer::<Dna, 7>::try_from(dna!("GTGAAGA")).unwrap()
        );
    }

    #[test]
    fn dna_kmer_macro() {
        assert_eq!(
            kmer!("TGCACATG"),
            Kmer::<Dna, 8>::try_from(dna!("TGCACATG")).unwrap()
        );
        assert_ne!(
            kmer!("GTGACGA"),
            Kmer::<Dna, 7>::try_from(dna!("GTGAAGA")).unwrap()
        );
    }

    #[test]
    fn dna_kmer_complement() {
        assert_eq!(
            format!(
                "{:b}",
                Kmer::<Dna, 8>::try_from(dna!("AAAAAAAA"))
                    .unwrap()
                    .comp()
                    .bs
            ),
            format!(
                "{:b}",
                Kmer::<Dna, 8>::try_from(dna!("TTTTTTTT")).unwrap().bs
            )
        );

        assert_eq!(
            Kmer::<Dna, 1>::try_from(dna!("C")).unwrap().comp(),
            Kmer::<Dna, 1>::try_from(dna!("G")).unwrap()
        );

        assert_eq!(
            Kmer::<Dna, 16>::try_from(dna!("AAAATGCACATGTTTT"))
                .unwrap()
                .comp(),
            Kmer::<Dna, 16>::try_from(dna!("TTTTACGTGTACAAAA")).unwrap()
        );
    }
}
