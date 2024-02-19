use crate::codec::{dna, Codec, Complement};
use crate::ParseBioError;

#[derive(Copy, Clone, PartialEq, Debug, Ord, Eq, PartialOrd, Hash)]
#[repr(transparent)]
pub struct Dna(u8);

impl Codec for Dna {
    type Error = ParseBioError;

    const WIDTH: u8 = 8;

    fn unsafe_from_bits(b: u8) -> Self {
        Self(b)
    }

    fn try_from_bits(b: u8) -> Result<Self, Self::Error> {
        Ok(Self(b))
    }

    fn from_char(c: char) -> Result<Self, Self::Error> {
        match u8::try_from(c) {
            Ok(b) => Ok(Self(b)),
            _ => Err(Self::Error {}),
        }
    }

    fn to_char(self) -> char {
        self.0.into()
    }
}

//impl Eq for Dna {
//    fn eq(self, rhs: Self) {
//        true
//    }
//}

impl Complement for Dna {
    fn comp(self) -> Self {
        unimplemented!()
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
            b'a' => Ok(dna::Dna::A),
            b'c' => Ok(dna::Dna::C),
            b'g' => Ok(dna::Dna::G),
            b't' => Ok(dna::Dna::T),
            _ => Err(ParseBioError {}),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::text;
    use crate::prelude::*;

    #[test]
    fn test_text_dna_encoding() {
        let v: Vec<u8> = vec![b'A', b'a', b'a', b'c', b'C', b'c', b'G', b'g', b'T', b'a'];
        let s: Seq<text::Dna> = v.into();
        assert_eq!(s.nth(0), Dna::A.into());
        assert_eq!(s.nth(4), Dna::C.into());
        assert_eq!(s.nth(6), Dna::G.into());
        assert_eq!(s.nth(8), Dna::T.into());

        assert_eq!(Dna::A, s.nth(0).try_into().unwrap());
        assert_eq!(Dna::C, s.nth(4).try_into().unwrap());
        assert_eq!(Dna::G, s.nth(6).try_into().unwrap());
        assert_eq!(Dna::T, s.nth(8).try_into().unwrap());

        assert_ne!(s.nth(1), Dna::G.into());
        assert_ne!(s.nth(3), Dna::C.into());
    }
}
