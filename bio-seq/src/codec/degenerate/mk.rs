use crate::ComplementMut;
use crate::codec::Codec;

/// 1-bit encoding for nucleotides with a**M**ino (`A`/`C`) and **K**etone (`T`/`G`) functional groups.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum MK {
    M = 0b0,
    K = 0b1,
}

impl Codec for MK {
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
            Some(unsafe { std::mem::transmute::<u8, MK>(b) })
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
            b'M' | b'A' | b'C' => Some(MK::M),
            b'K' | b'T' | b'G' => Some(MK::K),
            _ => None,
        }
    }

    fn to_char(self) -> char {
        match self {
            MK::M => 'M',
            MK::K => 'K',
        }
    }

    fn to_bits(self) -> u8 {
        self as u8
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![MK::M, MK::K].into_iter()
    }
}

impl ComplementMut for MK {
    /// This representation preserves complementarity, `M = comp(K)`
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
        let seq = Seq::<degenerate::MK>::from_str("MKMKKMKKKKKMMMMKM").unwrap();
        let seq_rc: Seq<degenerate::MK> = seq.to_revcomp();
        assert_eq!("KMKKKKMMMMMKMMKMK", String::from(seq_rc));
    }
}
