//! Coding/Decoding trait for bit-packable enums representing biological alphabets
//!
//! The [dna], [iupac], [text], and [amino] alphabets are built in.
//!
//! ## Deriving custom Codecs
//!
//! ```
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

pub trait Codec: Copy + Clone + Into<u8> + PartialEq + Hash + Eq {
    const WIDTH: u8;
    type Error: std::error::Error + core::fmt::Display;

    fn unsafe_from_bits(b: u8) -> Self;
    fn try_from_bits(b: u8) -> Result<Self, Self::Error>;
    fn from_char(c: char) -> Result<Self, Self::Error>;
    fn to_char(self) -> char;
}

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
        assert_ne!(Iupac::from(Dna::A), Iupac::T);
        assert_eq!(Iupac::from(Dna::C), Iupac::C);
    }
}
