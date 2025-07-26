use crate::ComplementMut;
use crate::codec::Codec;

/// 1-bit encoding for nucleotides with pu**R**ine (`A`/`G`) and p**Y**ramidine (`T`/`C`) structures.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum RY {
    R = 0b0,
    Y = 0b1,
}

impl Codec for RY {
    const BITS: u8 = 1;

    /// Transmute a `u8` into a degenerate 1-bit nucleotide
    ///
    /// SAFETY: This only looks at the lower 2 bits of the `u8`
    fn unsafe_from_bits(b: u8) -> Self {
        debug_assert!(b < 2);
        unsafe { std::mem::transmute(b & 0b1) }
    }

    /// Valid values are `0` and `1`
    fn try_from_bits(b: u8) -> Option<Self> {
        if b < 2 {
            Some(unsafe { std::mem::transmute::<u8, RY>(b) })
        } else {
            None
        }
    }

    /// TODO: fast translation of A, T, W to 0 and C, G, S to 1
    fn unsafe_from_ascii(_b: u8) -> Self {
        todo!()
    }

    fn try_from_ascii(c: u8) -> Option<Self> {
        match c {
            b'R' | b'A' | b'G' => Some(RY::R),
            b'Y' | b'T' | b'C' => Some(RY::Y),
            _ => None,
        }
    }

    fn to_char(self) -> char {
        match self {
            RY::R => 'R',
            RY::Y => 'Y',
        }
    }

    fn to_bits(self) -> u8 {
        self as u8
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![RY::R, RY::Y].into_iter()
    }
}

impl ComplementMut for RY {
    /// This representation preserves complementarity, `R = comp(Y)`
    fn comp(&mut self) {
        *self = unsafe { std::mem::transmute(*self as u8 ^ 1) };
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::degenerate;
    use crate::prelude::*;

    #[test]
    fn test_1bit() {
        let seq = Seq::<degenerate::RY>::from_str("RRYRYRRY").unwrap();
        let seq_rc: Seq<degenerate::RY> = seq.to_revcomp();
        assert_eq!("RYYRYRYY", String::from(seq_rc));
    }
}
