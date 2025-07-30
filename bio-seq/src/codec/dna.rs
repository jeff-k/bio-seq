//! 2-bit DNA representation: `A: 00, C: 01, G: 10, T: 11`

use crate::codec::Codec;
//use crate::kmer::Kmer;
//use crate::seq::{Seq, SeqArray, SeqSlice};
use crate::{Bv, Complement, ComplementMut};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Dna {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl Codec for Dna {
    type Store = Bv;
    const BITS: u8 = 2;

    /// Transmute a `u8` into a nucleotide
    ///
    /// SAFETY: This only looks at the lower 2 bits of the `u8`
    fn unsafe_from_bits(b: u8) -> Self {
        debug_assert!(b < 4);
        unsafe { std::mem::transmute(b & 0b11) }
    }

    /// We can verify that a byte is a valid `Dna` value if it's
    /// between 0 and 3.
    fn try_from_bits(b: u8) -> Option<Self> {
        if b < 4 {
            Some(unsafe { std::mem::transmute::<u8, Dna>(b) })
        } else {
            None
        }
    }

    /// The ASCII values of 'A', 'C', 'G', and 'T' can be translated into
    /// the numbers 0, 1, 2, and 3 using bitwise operations: `((b << 1) + b) >> 3`.
    /// In other words, multiply the ASCII value by 3 and shift right.
    fn unsafe_from_ascii(b: u8) -> Self {
        // TODO: benchmark against b * 3
        Dna::unsafe_from_bits(((b << 1) + b) >> 3)
    }

    fn try_from_ascii(c: u8) -> Option<Self> {
        match c {
            b'A' => Some(Dna::A),
            b'C' => Some(Dna::C),
            b'G' => Some(Dna::G),
            b'T' => Some(Dna::T),
            _ => None,
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

    fn to_bits(self) -> u8 {
        self as u8
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![Dna::A, Dna::C, Dna::G, Dna::T].into_iter()
    }
}

/// This 2-bit representation of nucleotides lends itself to a very fast
/// complement implementation with bitwise xor
impl ComplementMut for Dna {
    fn comp(&mut self) {
        *self = Dna::unsafe_from_bits(*self as u8 ^ 0b11);
    }
}

impl Complement for Dna {}

/*
impl ComplementMut for Seq<Dna> {
    fn comp(&mut self) {
        for word in self.bv.as_raw_mut_slice() {
            *word ^= usize::MAX;
        }
    }
}
*/

/*
impl ReverseMut for Seq<Dna> {
    fn rev(&mut self) {
        self.bv.reverse();
        for word in self.bv.as_raw_mut_slice().iter_mut() {
            let c: usize = *word;
            let odds = (c & 0x5555_5555_5555_5555usize) >> 1;
            *word = (c & 0xAAAA_AAAA_AAAA_AAAAusize) << 1 | odds;
        }
    }
}
*/

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
            Kmer::<Dna, 16>::from(dna!("AAAATGCACATGTTTT")).comp(),
            Kmer::<Dna, 16>::from(dna!("TTTTACGTGTACAAAA"))
        );
    }
    */
}
