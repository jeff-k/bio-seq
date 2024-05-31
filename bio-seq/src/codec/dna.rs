//! 2-bit DNA representation: `A: 00, C: 01, G: 10, T: 11`

use crate::codec::{Codec, Complement};

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

    /// Take the two least significant bits of a `u8` and map them to the
    /// corresponding nucleotides.
    fn unsafe_from_bits(b: u8) -> Self {
        unsafe { std::mem::transmute(b & 0b11) }
    }

    /// We can efficient verify that a byte is a valid `Dna` value if it's
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
    fn unsafe_from_ascii(b: u8) -> Self {
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

    fn items() -> impl Iterator<Item = Self> {
        vec![Dna::A, Dna::C, Dna::G, Dna::T].into_iter()
    }
}

impl Complement for Dna {
    /// This 2-bit representation of nucleotides lends itself to a very fast
    /// complement implementation with bitwise xor
    fn comp(&self) -> Self {
        // flip the bits
        let b = *self as u8 ^ 0b11;
        Dna::unsafe_from_bits(b)
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
