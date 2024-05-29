//! 2-bit DNA representation: `A: 00, C: 01, G: 10, T: 11`
use core::marker::PhantomData;
//use std::intrinsics::mir::SetDiscriminant;

use crate::codec::{Codec, Complement, IntoComplement};
use crate::kmer::Kmer;
use crate::Order;

use bitvec::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl From<Dna> for u8 {
    fn from(b: Dna) -> u8 {
        b as u8
    }
}

impl Codec for Dna {
    const BITS: u8 = 2;
    type Error = crate::error::ParseBioError;

    fn unsafe_from_bits(b: u8) -> Self {
        unsafe { std::mem::transmute(b & 0b11) }
    }

    fn try_from_bits(b: u8) -> Option<Self> {
        if b < 4 {
            Some(unsafe { std::mem::transmute(b) })
        } else {
            None
        }
    }

    fn unsafe_from_ascii(b: u8) -> Self {
        Dna::unsafe_from_bits(((b << 1) + b) >> 3)
    }

    fn try_from_ascii(c: u8) -> Result<Self, Self::Error> {
        match c {
            b'A' => Ok(Dna::A),
            b'C' => Ok(Dna::C),
            b'G' => Ok(Dna::G),
            b'T' => Ok(Dna::T),
            _ => Err(crate::error::ParseBioError {}),
        }
    }

    fn to_char(self) -> char {
        match self {
            Dna::A => 'A',
            Dna::C => 'C',
            Dna::G => 'G',
            Dna::T => 'T',
        }
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![Dna::A, Dna::C, Dna::G, Dna::G].into_iter()
    }
}

/// Bitwise negation results in the complement of a base. This can be used to efficiently complement a bit-packed sequence of Dna.
impl<const K: usize> IntoComplement for Kmer<Dna, K> {
    fn into_comp(self: &Kmer<Dna, K>) -> Self {
        Kmer {
            _p: PhantomData,
            bs: (self.bs ^ usize::MAX).view_bits::<Order>()[..K * Dna::BITS as usize]
                .load_le::<usize>(),
        }
    }
}

impl Complement for Dna {
    fn comp(&mut self) {
        *self = Dna::unsafe_from_bits((*self as u8) ^ 0b11)
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

    /*
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
    */
}
