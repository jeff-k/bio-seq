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
//! This naturally supports set membership operations:
//!
//! ```rust
//! use bio_seq::prelude::*;
//!
//! // Set union:
//! assert_eq!(iupac!("AS-GYTNA") | iupac!("ANTGCAT-"), iupac!("ANTGYWNA"));
//!
//! // Set intersection:
//! assert_eq!(iupac!("ACGTSWKM") & iupac!("WKMSTNNA"), iupac!("A----WKA"));
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
use crate::seq::{Seq, SeqSlice};

use core::ops::{BitAnd, BitOr};

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

impl BitAnd for Seq<Iupac> {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        self.bit_and(rhs)
    }
}

impl BitOr for Seq<Iupac> {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        self.bit_or(rhs)
    }
}

impl BitAnd for &SeqSlice<Iupac> {
    type Output = Seq<Iupac>;

    fn bitand(self, rhs: Self) -> Self::Output {
        self.bit_and(rhs)
    }
}

impl BitOr for &SeqSlice<Iupac> {
    type Output = Seq<Iupac>;

    fn bitor(self, rhs: Self) -> Self::Output {
        self.bit_or(rhs)
    }
}

impl Seq<Iupac> {
    pub fn contains(&self, rhs: &SeqSlice<Iupac>) -> bool {
        if rhs.len() != self.len() {
            panic!("Cannot compare IUPAC sequences of different length");
        }
        let slice: &SeqSlice<Iupac> = self;
        let intersection: &SeqSlice<Iupac> = &(slice & rhs);
        intersection == rhs
    }
}

#[macro_export]
macro_rules! iupac {
    ($seq:expr) => {
        match Seq::<Iupac>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}
