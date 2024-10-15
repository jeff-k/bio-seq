use crate::codec::Codec;
use crate::seq::{ReverseComplement, Seq};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Codec)]
#[bits(5)]
#[repr(u8)]
pub enum Iupac {
    A = 0b10_0_00,
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

    #[display('-')]
    AMasked = 0b10_1_00,
    #[display('-')]
    CMasked = 0b01100,
    #[display('-')]
    GMasked = 0b00110,
    #[display('-')]
    TMasked = 0b00101,

    #[display('-')]
    YMasked = 0b01101,
    #[display('-')]
    RMasked = 0b10110,
    #[display('-')]
    WMasked = 0b10101,
    #[display('-')]
    SMasked = 0b01110,

    #[display('-')]
    KMasked = 0b00111,
    #[display('-')]
    MMasked = 0b11100,
    #[display('-')]
    DMasked = 0b10111,
    #[display('-')]
    VMasked = 0b11110,

    #[display('-')]
    HMasked = 0b11101,
    #[display('-')]
    BMasked = 0b01111,

    #[display('-')]
    NMasked = 0b11111,
    #[display('-')]
    XMasked = 0b00100,

}

/*
impl Complement for Dna {
    /// This representation can be complemented by reversing the bit pattern
    fn comp(&self) -> Self {
        // reverse the bits
        Dna::unsafe_from_bits(b)
    }
}
*/

/*
impl Dna {
    /// Flipping the bit pattern masks/unmasks this representation
    pub fn mask(&self) -> Self {
        let b = *self as u8 ^ 0b1111;
        Dna::unsafe_from_bits(b)
    }
}

impl Seq<Dna> {
    pub fn mask(&self) -> Self {
        Self::from(!self.bv.clone())
    }
}

impl ReverseComplement for Seq<Dna> {
    type Output = Self;

    fn revcomp(&self) -> Self {
        let mut bv = self.bv.clone();
        bv.reverse();
        Self::from(bv)
    }
}
*/
#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::prelude::*;


}
