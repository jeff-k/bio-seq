//! 4-bit IUPAC nucleotide ambiguity codes
//!
//! IUPAC nucleotide ambiguity codes are represented with 4 bits
//!
//! |   | A | C | G | T |
//! | - | - | - | - | - |
//! | A | 1 | 0 | 0 | 0 |
//! | C | 0 | 1 | 0 | 0 |
//! | G | 0 | 0 | 1 | 0 |
//! | T | 0 | 0 | 0 | 1 |
//! | Y | 0 | 1 | 0 | 1 |
//! | R | 1 | 0 | 1 | 0 |
//! | W | 1 | 0 | 0 | 1 |
//! | S | 0 | 1 | 1 | 0 |
//! | K | 0 | 0 | 1 | 1 |
//! | M | 1 | 1 | 0 | 0 |
//! | D | 1 | 0 | 1 | 1 |
//! | V | 1 | 1 | 1 | 0 |
//! | H | 1 | 1 | 0 | 1 |
//! | B | 0 | 1 | 1 | 1 |
//! | N | 1 | 1 | 1 | 1 |
//! | X/- | 0 | 0 | 0 | 0 |
//!
//! This means that we can treat each symbol as a set and we get meaningful bitwise operations:
//!
//! ```rust
//! use bio_seq::prelude::*;
//!
//! // Set union:
//! let union = iupac!("AS-GYTNAN") | iupac!("ANTGCAT-N");
//! assert_eq!(union, iupac!("ANTGYWNAN"));
//!
//! // Set intersection:
//! let intersection = iupac!("ACGTSWKMN") & iupac!("WKMSTNNAN");
//! assert_eq!(intersection, iupac!("A----WKAN"));
//! ```
//!
//! Which can be used to implement pattern matching:
//!
//! ```rust
//! use bio_seq::prelude::*;
//!
//! let seq = iupac!("AGCTNNCAGTCGACGTATGTA");
//! let pattern = iupac!("AYG");
//!
//! for slice in seq.windows(pattern.len()) {
//!    if pattern.contains(slice) {
//!        println!("{slice} matches pattern");
//!    }
//! }
//!
//! // ACG matches pattern
//! // ATG matches pattern
//! ```
use crate::codec::{dna::Dna, Codec};
use crate::seq::{Seq, SeqArray, SeqSlice};
use crate::{Complement, ComplementMut};

const IUPAC_COMPLEMENT_TABLE: [u8; 16] = {
    let mut table = [0; 16];

    table[Iupac::A as usize] = Iupac::T as u8;
    table[Iupac::C as usize] = Iupac::G as u8;
    table[Iupac::G as usize] = Iupac::C as u8;
    table[Iupac::T as usize] = Iupac::A as u8;
    table[Iupac::Y as usize] = Iupac::R as u8;
    table[Iupac::R as usize] = Iupac::Y as u8;
    table[Iupac::W as usize] = Iupac::W as u8;
    table[Iupac::S as usize] = Iupac::S as u8;
    table[Iupac::K as usize] = Iupac::M as u8;
    table[Iupac::M as usize] = Iupac::K as u8;
    table[Iupac::D as usize] = Iupac::H as u8;
    table[Iupac::V as usize] = Iupac::B as u8;
    table[Iupac::H as usize] = Iupac::D as u8;
    table[Iupac::B as usize] = Iupac::V as u8;
    table[Iupac::N as usize] = Iupac::N as u8;

    table
};

impl From<Iupac> for u8 {
    fn from(b: Iupac) -> u8 {
        b as u8
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Codec)]
#[bits(4)]
#[repr(u8)]
pub enum Iupac {
    A = 0b1000,
    C = 0b0100,
    G = 0b0010,
    T = 0b0001,
    R = 0b1010,
    Y = 0b0101,
    S = 0b0110,
    W = 0b1001,
    K = 0b0011,
    M = 0b1100,
    B = 0b0111,
    D = 0b1011,
    H = 0b1101,
    V = 0b1110,
    N = 0b1111,
    #[display('-')]
    X = 0b0000,
}

impl From<Dna> for Iupac {
    fn from(dna: Dna) -> Self {
        match dna {
            Dna::A => Iupac::A,
            Dna::C => Iupac::C,
            Dna::G => Iupac::G,
            Dna::T => Iupac::T,
        }
    }
}

impl Seq<Iupac> {
    pub fn contains(&self, rhs: &SeqSlice<Iupac>) -> bool {
        if rhs.len() != self.len() {
            return false;
        }

        self.as_ref() & rhs == rhs
    }
}

impl<const N: usize, const W: usize> SeqArray<Iupac, N, W> {
    pub fn contains(&self, rhs: &SeqSlice<Iupac>) -> bool {
        if N != rhs.len() {
            return false;
        }
        self.as_ref() & rhs == rhs
    }
}

impl SeqSlice<Iupac> {
    pub fn contains(&self, rhs: &SeqSlice<Iupac>) -> bool {
        if self.len() != rhs.len() {
            return false;
        }
        self & rhs == rhs
    }
}

/// The complement of an IUPAC base is the reverse of the bit-pattern
impl ComplementMut for Iupac {
    fn comp(&mut self) {
        // Below are two methods:
        // 1. Using a lookup table.
        // 2. The 7 operation from https://graphics.stanford.edu/~seander/bithacks.html
        // See: https://stackoverflow.com/questions/3587826/is-there-a-built-in-function-to-reverse-bit-order

        // Use a lookup table
        *self = Iupac::unsafe_from_bits(IUPAC_COMPLEMENT_TABLE[*self as usize]);

        // Use the The 7 operation from
        // let b = *self as u32;
        // *self = Iupac::unsafe_from_bits(
        //     ((((((b * 0x0802u32) & 0x22110u32) | ((b * 0x8020u32) & 0x88440u32)) * 0x10101u32)
        //         >> 20) // 16 + 4
        //         & 0x0f) as u8,
        // );
    }
}

impl Complement for Iupac {}

/*
impl Complement for Seq<Iupac> {
    type Output = Self;

    fn comp(&mut self) {
        todo!()
    }

    fn to_comp(&self) -> Self::Output {
        todo!()
    }
}

/// Reverse complementing a sequence of IUPAC characters is simply a
/// matter of reversing the entire bit sequence
impl ReverseComplement for Seq<Iupac> {
    type Output = Self;

    fn revcomp(&mut self) {
        // simply do bit-wise reversal
        self.bv.reverse();
    }

    fn to_revcomp(&self) -> Self::Output {
        todo!()
    }
}
*/

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn iupac_ops() {
        let seq = iupac!("AGCTNNCAGTCGACGTATGTA");

        let pattern = iupac!("AYG");

        let matches: Vec<Seq<Iupac>> = seq
            .windows(pattern.len())
            .filter(|w| pattern.contains(w))
            .collect();

        assert_eq!(matches, vec![iupac!("ACG"), iupac!("ATG")]);
    }

    #[test]
    fn iupac_complement() {
        assert_eq!(
            iupac!("AGCTYRWSKMDVHBN").to_comp(),
            iupac!("TCGARYWSMKHBDVN")
        );
    }
}
