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
//! ## Implementing custom Codecs
//!
//! Custom encodings can be defined on enums by implementing the `Codec` trait.
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

//#[macro_use]
//pub mod amino;
#[macro_use]
pub mod dna;
//#[macro_use]
//pub mod iupac;

//pub mod text;

pub use bio_seq_derive::Codec;

/// The bit encodings of an alphabet's characters can be represented with any type.
/// Encoding from ASCII bytes and decoding the representation is implemented through
/// the `Codec` trait.  
///
/// The intended representation is an `Enum`, transparently represented as a `u8`.
pub trait Codec: Copy + Clone + PartialEq + Hash + Eq + Into<u8> {
    /// Associated parsing error for unknown characters
    type Error;

    /// The number of bits used to encode the characters. e.g. `Dna::BITS` = 2, `Iupac::BITS` = 4.
    const BITS: u8;

    /// Convert raw bits of binary encoding into enum item. Binary values
    /// that don't match an enum member's discriminant will result in panic or random enum
    /// item
    fn unsafe_from_bits(b: u8) -> Self;

    /// Fallibly convert raw bits into enum. If the binary value does not
    /// match a discriminant, return `None`
    fn try_from_bits(b: u8) -> Option<Self>;

    /// Encode an ASCII byte as a codec enum item
    fn unsafe_from_ascii(c: u8) -> Self;

    /// Fallibly encode an ASCII byte as a codec enum item
    fn try_from_ascii(c: u8) -> Result<Self, Self::Error>;

    /// Decode enum item as a UTF-8 character
    fn to_char(self) -> char;

    /// Iterator over the characters of the codec
    fn items() -> impl Iterator<Item = Self>;
}

/// Nucleotide alphabets that can be complemented implement `FromComplement`
///
/// ```
/// use bio_seq::prelude::{Dna, Complement};
/// assert_eq!(Dna::A.from_comp(), Dna::T);
/// ````
pub trait IntoComplement {
    fn into_comp(&self) -> Self;
}

/// `Complement` a value in place
pub trait Complement {
    fn comp(&mut self);
}

/*
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
*/
