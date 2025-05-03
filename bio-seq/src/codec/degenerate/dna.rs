use crate::ComplementMut;
use crate::codec::Codec;

/// 1-bit encoding for `S`trong (`G`/`C`) and `W`eak (`A`/`T`) binding nucleotides
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Dna {
    W = 0b0,
    S = 0b1,
}

impl Codec for Dna {
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
            Some(unsafe { std::mem::transmute::<u8, Dna>(b) })
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
            b'S' | b'C' | b'G' => Some(Dna::S),
            b'W' | b'A' | b'T' => Some(Dna::W),
            _ => None,
        }
    }

    fn to_char(self) -> char {
        match self {
            Dna::S => 'S',
            Dna::W => 'W',
        }
    }

    fn to_bits(self) -> u8 {
        self as u8
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![Dna::S, Dna::W].into_iter()
    }
}

impl ComplementMut for Dna {
    /// This representation erases complements, so this is the identify function
    fn comp(&mut self) {}
}

#[cfg(test)]
mod tests {
    use crate::codec::degenerate;
    use crate::prelude::*;

    #[test]
    fn test_1bit() {
        let seq = Seq::<degenerate::Dna>::from_str("SSSWWWSW").unwrap();
        let seq_rc: Seq<degenerate::Dna> = seq.to_revcomp();
        assert_eq!("WSWWWSSS", String::from(seq_rc));
    }
}
