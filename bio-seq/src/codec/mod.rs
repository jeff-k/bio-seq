//! Coding/Decoding trait for bit-packable enums representing biological alphabets
//!
//! The [dna], [iupac], [text], and [amino] alphabets are built in.
//!
//! This trait implements the translation between the UTF-8 representation of an alphabet and it's efficient bit-packing.
//! The `BITS` attribute stores the number of bits used by the representation.
//! ```
//! use bio_seq::prelude::{Dna, Codec};
//! use bio_seq::codec::text;
//! assert_eq!(Dna::BITS, 2);
//! assert_eq!(text::Dna::BITS, 8);
//! ```
//!
//! ## Deriving custom Codecs
//!
//! Custom encodings can be easily defined on enums using the derivable `Codec` trait.
//!
//! ```ignore
//! use bio_seq::prelude;
//! use bio_seq::prelude::Codec;
//!
//! #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Codec)]
//! pub enum Dna {
//!     A = 0b00,
//!     C = 0b01,
//!     G = 0b10,
//!     T = 0b11,
//! }
//! ```
use core::hash::Hash;

#[macro_use]
pub mod amino;
#[macro_use]
pub mod dna;
#[macro_use]
pub mod iupac;

pub mod text;

pub use bio_seq_derive::Codec;

/// The binary encodings of an alphabet's characters are represented with `u8`s. Encoding from UTF-8 or
/// a raw `u8` will always be fallible but often can be assumed safe.
pub trait Codec: Copy + Clone + Into<u8> + PartialEq + Hash + Eq {
    const BITS: usize;
    type Error: std::error::Error + core::fmt::Display;

    fn unsafe_from_bits(b: u8) -> Self;
    fn try_from_bits(b: u8) -> Result<Self, Self::Error>;
    fn from_char(c: char) -> Result<Self, Self::Error>;
    fn to_char(self) -> char;
}

/// Nucleotide alphabets that can be complemented implement `Complement`
///
/// ```
/// use bio_seq::prelude::{Dna, Complement};
/// assert_eq!(Dna::A.comp(), Dna::T);
/// ````
pub trait Complement {
    fn comp(self) -> Self;
}

#[cfg(test)]
mod tests {
    use super::dna::Dna;
    use super::iupac::Iupac;

    #[test]
    fn dna_to_iupac() {
        assert_eq!(Iupac::from(Dna::A), Iupac::A);
        assert_eq!(Iupac::from(Dna::C), Iupac::C);
        assert_eq!(Iupac::from(Dna::G), Iupac::G);
        assert_eq!(Iupac::from(Dna::T), Iupac::T);

        assert_ne!(Iupac::from(Dna::A), Iupac::T);
        assert_ne!(Iupac::from(Dna::T), Iupac::A);
        assert_ne!(Iupac::from(Dna::C), Iupac::T);
        assert_ne!(Iupac::from(Dna::G), Iupac::T);
    }
}
