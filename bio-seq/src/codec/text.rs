//! 8-bit ASCII representation of nucleotides
//!
//! This encoding is a literal interpretation of bytes of text as DNA
use crate::codec::{dna, Codec, Complement};
use crate::error::ParseBioError;

#[derive(Copy, Clone, PartialEq, Debug, Ord, Eq, PartialOrd, Hash)]
#[repr(transparent)]
pub struct Dna(u8);

impl Codec for Dna {
    const BITS: u8 = 8;

    fn unsafe_from_bits(b: u8) -> Self {
        Self(b)
    }

    fn try_from_bits(b: u8) -> Option<Self> {
        Some(Self(b))
    }

    fn unsafe_from_ascii(c: u8) -> Self {
        Self(c)
    }

    fn try_from_ascii(c: u8) -> Option<Self> {
        //        if c.is_ascii_alphanumeric() {
        match c {
            b'A' | b'C' | b'G' | b'T' | b'N' => Some(Self(c)),
            _ => None,
        }
    }

    fn to_char(self) -> char {
        self.0.into()
    }

    fn to_bits(self) -> u8 {
        self.0
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![Dna(b'A'), Dna(b'C'), Dna(b'G'), Dna(b'T'), Dna(b'N')].into_iter()
    }
}

//impl Eq for Dna {
//    fn eq(self, rhs: Self) {
//        todo!()
//    }
//}

impl Complement for Dna {
    fn comp(&self) -> Self {
        match self {
            Self(b'A') => Self(b'T'),
            Self(b'C') => Self(b'G'),
            Self(b'G') => Self(b'C'),
            Self(b'T') => Self(b'A'),
            _ => Self(b'N'),
        }
    }
}

impl From<Dna> for u8 {
    fn from(val: Dna) -> Self {
        val.0
    }
}

impl From<dna::Dna> for Dna {
    fn from(base: dna::Dna) -> Self {
        match base {
            dna::Dna::A => Dna(b'A'),
            dna::Dna::C => Dna(b'C'),
            dna::Dna::G => Dna(b'G'),
            dna::Dna::T => Dna(b'T'),
        }
    }
}

impl TryFrom<Dna> for dna::Dna {
    type Error = ParseBioError;

    fn try_from(base: Dna) -> Result<Self, Self::Error> {
        match base.0 {
            b'A' => Ok(dna::Dna::A),
            b'C' => Ok(dna::Dna::C),
            b'G' => Ok(dna::Dna::G),
            b'T' => Ok(dna::Dna::T),
            // Todo: decide whether to support lower cases
            /*
            b'a' => Ok(dna::Dna::A),
            b'c' => Ok(dna::Dna::C),
            b'g' => Ok(dna::Dna::G),
            b't' => Ok(dna::Dna::T),
            */
            _ => Err(ParseBioError::UnrecognisedBase(base.0)),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::text;
    use crate::prelude::*;

    #[test]
    fn test_text_dna_encoding() {
        let s: Seq<text::Dna> = "CATCGCGACTGATCACTCGATC".try_into().unwrap();
        println!("{s}");
        assert_eq!(s.nth(0), Dna::C.into());
        assert_eq!(s.nth(4), Dna::G.into());
        assert_eq!(s.nth(6), Dna::G.into());
        assert_eq!(s.nth(8), Dna::C.into());

        assert_eq!(Dna::C, s.nth(0).try_into().unwrap());
        assert_eq!(Dna::G, s.nth(4).try_into().unwrap());
        assert_eq!(Dna::G, s.nth(6).try_into().unwrap());
        assert_eq!(Dna::C, s.nth(8).try_into().unwrap());

        assert_ne!(s.nth(1), Dna::T.into());
        assert_ne!(s.nth(3), Dna::A.into());
    }
}
