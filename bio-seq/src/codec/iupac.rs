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

// const IUPAC_COMPLEMENT_TABLE: [u8; 256] = [
//     0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0, 0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0, 0x70, 0xf0,
//     0x08, 0x88, 0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8, 0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8,
//     0x04, 0x84, 0x44, 0xc4, 0x24, 0xa4, 0x64, 0xe4, 0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4,
//     0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac, 0x6c, 0xec, 0x1c, 0x9c, 0x5c, 0xdc, 0x3c, 0xbc, 0x7c, 0xfc,
//     0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2, 0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2,
//     0x0a, 0x8a, 0x4a, 0xca, 0x2a, 0xaa, 0x6a, 0xea, 0x1a, 0x9a, 0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa,
//     0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6, 0x16, 0x96, 0x56, 0xd6, 0x36, 0xb6, 0x76, 0xf6,
//     0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee, 0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe, 0x7e, 0xfe,
//     0x01, 0x81, 0x41, 0xc1, 0x21, 0xa1, 0x61, 0xe1, 0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
//     0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9, 0x19, 0x99, 0x59, 0xd9, 0x39, 0xb9, 0x79, 0xf9,
//     0x05, 0x85, 0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5, 0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5,
//     0x0d, 0x8d, 0x4d, 0xcd, 0x2d, 0xad, 0x6d, 0xed, 0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd,
//     0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3, 0x63, 0xe3, 0x13, 0x93, 0x53, 0xd3, 0x33, 0xb3, 0x73, 0xf3,
//     0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb, 0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb,
//     0x07, 0x87, 0x47, 0xc7, 0x27, 0xa7, 0x67, 0xe7, 0x17, 0x97, 0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7,
//     0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef, 0x1f, 0x9f, 0x5f, 0xdf, 0x3f, 0xbf, 0x7f, 0xff,
// ];

/*
const LTABLE: [u8; 256] = {
    let mut table = [0; 256];

    table[b'A' as usize] = 0b1000;
    table[b'C' as usize] = 0b0100;
    table[b'G' as usize] = 0b0010;
    table[b'T' as usize] = 0b0001;
    table[b'R' as usize] = 0b1010;
    table[b'Y' as usize] = 0b0101;
    table[b'S' as usize] = 0b0110;
    table[b'W' as usize] = 0b1001;
    table[b'K' as usize] = 0b0011;
    table[b'M' as usize] = 0b1100;
    table[b'B' as usize] = 0b0111;
    table[b'D' as usize] = 0b1011;
    table[b'H' as usize] = 0b1101;
    table[b'V' as usize] = 0b1110;
    table[b'N' as usize] = 0b1111;
    table[b'-' as usize] = 0b0000;

    table
};
*/

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
        // *self = Iupac::unsafe_from_bits(IUPAC_COMPLEMENT_TABLE[*self as usize] >> 4 & 0x0f);

        // Use the The 7 operation from
        let b = *self as u32;
        *self = Iupac::unsafe_from_bits(
            ((((((b * 0x0802u32) & 0x22110u32) | ((b * 0x8020u32) & 0x88440u32)) * 0x10101u32)
                >> 20) // 16 + 4
                & 0x0f) as u8,
        );
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
