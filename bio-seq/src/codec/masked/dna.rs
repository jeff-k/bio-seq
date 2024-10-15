use crate::codec::masked::Maskable;
use crate::codec::{Codec, Complement};
use crate::seq::{ReverseComplement, Seq};

/// **Experimental** 4-bit nucleotide encoding with fast reverse complement and toggled mask operation
///
/// Note that masking/unmasking are not idempotent
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Codec)]
#[bits(4)]
#[repr(u8)]
pub enum Dna {
    A = 0b1000,
    C = 0b0100,
    G = 0b0010,
    T = 0b0001,

    #[display('a')]
    AMasked = 0b0111,
    #[display('c')]
    CMasked = 0b1011,
    #[display('g')]
    GMasked = 0b1101,
    #[display('t')]
    TMasked = 0b1110,

    N = 0b0000,
    #[display('n')]
    NMasked = 0b1111,

    #[display('-')]
    #[alt(0b0011)]
    Gap = 0b1100,

    #[display('.')]
    #[alt(0b0101)]
    Pad = 0b1010,

    #[display('?')]
    Unknown1 = 0b0110,

    #[display('!')]
    Unknown2 = 0b1001,
}

impl Complement for Dna {
    /// This representation can be complemented by reversing the bit pattern
    fn comp(&self) -> Self {
        todo!()
    }
}

impl Maskable for Dna {
    /// Inverting the bit pattern masks/unmasks this representation
    fn mask(&self) -> Self {
        let b = *self as u8 ^ 0b1111;
        Dna::unsafe_from_bits(b)
    }

    fn unmask(&self) -> Self {
        let b = *self as u8 ^ 0b1111;
        Dna::unsafe_from_bits(b)
    }
}

impl Maskable for Seq<Dna> {
    fn mask(&self) -> Self {
        Self::from(!self.bv.clone())
    }

    fn unmask(&self) -> Self {
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

#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::codec::masked::Maskable;
    use crate::prelude::*;

    #[test]
    fn mask_sequence() {
        let seq = Seq::<masked::Dna>::try_from("A.TCGCgtcataN--A").unwrap();

        assert_ne!(seq.mask().to_string(), "a.tcgcGTGATAN--a".to_string());
        assert_eq!(seq.mask().to_string(), "a.tcgcGTCATAn--a".to_string());
    }

    #[test]
    fn masked_revcomp() {
        let seq = Seq::<masked::Dna>::try_from("A.TCGCgtcataN--A").unwrap();

        assert_ne!(seq.revcomp().to_string(), "T--NtaagacGCGA.T".to_string());
        assert_eq!(seq.revcomp().to_string(), "T--NtatgacGCGA.T".to_string());
    }
}
