use crate::codec::Codec;
//use crate::{Complement, Maskable, Reverse, ReverseComplement};
use crate::{ComplementMut, Maskable, MaskableMut}; //, ReverseComplementMut, ReverseMut};

/// 5-bit encoding for maskable IUPAC symbols
/// The middle bit is the mask flag and symbols are complemented by reversing the bit pattern
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Codec)]
#[bits(5)]
#[repr(u8)]
pub enum Iupac {
    A = 0b10000,
    C = 0b01000,
    G = 0b00010,
    T = 0b00001,

    Y = 0b01001,
    R = 0b10010,
    W = 0b10001,
    S = 0b01010,

    K = 0b00011,
    M = 0b11000,
    D = 0b10011,
    V = 0b11010,

    H = 0b11001,
    B = 0b01011,

    N = 0b11011,
    #[display('-')]
    X = 0b00000,

    #[display('a')]
    AMasked = 0b10100,
    #[display('c')]
    CMasked = 0b01100,
    #[display('g')]
    GMasked = 0b00110,
    #[display('t')]
    TMasked = 0b00101,

    #[display('y')]
    YMasked = 0b01101,
    #[display('r')]
    RMasked = 0b10110,
    #[display('w')]
    WMasked = 0b10101,
    #[display('s')]
    SMasked = 0b01110,

    #[display('k')]
    KMasked = 0b00111,
    #[display('m')]
    MMasked = 0b11100,
    #[display('d')]
    DMasked = 0b10111,
    #[display('v')]
    VMasked = 0b11110,

    #[display('h')]
    HMasked = 0b11101,
    #[display('b')]
    BMasked = 0b01111,

    #[display('n')]
    NMasked = 0b11111,
    #[display('.')]
    XMasked = 0b00100,
}

#[allow(clippy::cast_possible_truncation)]
impl ComplementMut for Iupac {
    /// This representation can be complemented by reversing the bit pattern
    fn comp(&mut self) {
        // reverse the bits
        *self = Self::unsafe_from_bits(
            ((self.to_bits() as usize * 0x0202).rotate_right(8) & 0b1111) as u8,
        );
    }
}

impl MaskableMut for Iupac {
    /// Setting the middle bit sets the mask flag
    fn mask(&mut self) {
        let b = *self as u8 | 0b00100;
        *self = Self::unsafe_from_bits(b);
    }

    /// Unsetting the middle bit clears the mask flag
    fn unmask(&mut self) {
        let b = *self as u8 & 0b11011;
        *self = Self::unsafe_from_bits(b);
    }
}

impl Maskable for Iupac {}

#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::prelude::*;

    #[ignore]
    #[test]
    fn mask_iupac_seq() {
        let seq = Seq::<masked::Iupac>::try_from("A.TCGCgtcataN--A").unwrap();

        assert_ne!(seq.to_mask().to_string(), "a.tcgcGTGATAN--a".to_string());
        assert_eq!(seq.to_mask().to_string(), "a.tcgcGTCATAn--a".to_string());
    }
}
