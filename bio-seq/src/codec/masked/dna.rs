use crate::codec::Codec;
//use crate::{Complement, Maskable, Reverse, ReverseComplement};
use crate::{Complement, ComplementMut, MaskableMut}; //, ReverseComplementMut, ReverseMut};

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

impl ComplementMut for Dna {
    /// This representation can be complemented by reversing the bit pattern
    fn comp(&mut self) {
        let bits = self.to_bits();
        *self = Self::unsafe_from_bits(
            ((bits & 0b1000) >> 3)
                | ((bits & 0b0100) >> 1)
                | ((bits & 0b0010) << 1)
                | ((bits & 0b0001) << 3),
        );
    }
}

impl Complement for Dna {}

impl MaskableMut for Dna {
    /// Inverting the bit pattern masks/unmasks this representation
    fn mask(&mut self) {
        let b = *self as u8 ^ 0b1111;
        *self = Dna::unsafe_from_bits(b);
    }

    fn unmask(&mut self) {
        let b = *self as u8 ^ 0b1111;
        *self = Dna::unsafe_from_bits(b);
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::prelude::*;

    #[test]
    fn iupac_comp() {
        let a: masked::Dna = masked::Dna::A;
        let t: masked::Dna = masked::Dna::T;
        let c: masked::Dna = masked::Dna::C;
        let g: masked::Dna = masked::Dna::G;
        assert_eq!(t.to_comp(), a);
        assert_eq!(c.to_comp(), g);

        assert_ne!(t.to_comp(), c);
        assert_ne!(c.to_comp(), c);

        assert_eq!(c.to_comp(), g.to_comp().to_comp());
    }

    #[test]
    fn mask_sequence() {
        let seq = Seq::<masked::Dna>::try_from("A.TCGCgtcataN--A").unwrap();

        assert_ne!(seq.to_mask().to_string(), "a.tcgcGTGATAN--a".to_string());
        assert_eq!(seq.to_mask().to_string(), "a.tcgcGTCATAn--a".to_string());
    }

    #[test]
    fn masked_comp() {
        let seq = Seq::<masked::Dna>::try_from("A.TCGCgtcataN--A").unwrap();

        assert_eq!(seq.to_comp().to_string(), "T.AGCGcagtatN--T".to_string());
    }

    #[test]
    fn masked_revcomp() {
        let seq = Seq::<masked::Dna>::try_from("A.TCGCgtcataN--A").unwrap();

        assert_ne!(seq.to_revcomp().to_string(), "T--NtaagacGCGA.T".to_string());
        assert_eq!(seq.to_revcomp().to_string(), "T--NtatgacGCGA.T".to_string());
    }
}
