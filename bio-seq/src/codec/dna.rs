//! 2-bit DNA representation
//!
//! Bitwise negation results in the complement.
use core::marker::PhantomData;

use crate::codec::{Codec, Complement};
use crate::kmer::Kmer;

use bitvec::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Codec)]
#[bits(2)]
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
            bs: (self.bs ^ usize::MAX).view_bits::<Lsb0>()[..K * Dna::BITS].load_le::<usize>(),
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
