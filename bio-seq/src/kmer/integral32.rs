use crate::kmer::{sealed, REV_2BIT};
use bitvec::field::BitField;

use crate::{Ba, Bs};

impl sealed::KmerStorage for usize {
    const BITS: usize = usize::BITS as usize;

    type BaN = Ba<1>;

    fn to_bitarray(self) -> Ba<{ (Self::BITS / usize::BITS) as usize }> {
        Self::BaN::new([self])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        debug_assert!(
            bs.len() <= Self::BITS as usize,
            "bitslice larger than kmer storage type"
        );
        bs.load_le()
    }

    /*
        fn rotate_right(self, n: u32) -> Self {
            self.rotate_right(n)
        }
        fn rotate_left(self, n: u32) -> Self {
            self.rotate_left(n)
        }
    */

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }

    fn complement(&mut self, mask: usize) {
        if mask >= Self::BITS as usize {
            *self ^= Self::MAX;
        } else {
            let mask = (1 << mask) - 1;
            *self ^= mask;
        }
    }

    fn shiftr(&mut self, n: u32) {
        *self >>= n;
    }

    fn shiftl(&mut self, n: u32) {
        *self <<= n;
    }

    fn rev_blocks_2(&mut self) {
        let mut bs = self.swap_bytes().to_le_bytes();

        for b in &mut bs {
            *b = REV_2BIT[*b as usize];
        }

        *self = Self::from_le_bytes(bs);
    }
}

impl sealed::KmerStorage for u64 {
    const BITS: usize = u64::BITS as usize;

    type BaN = Ba<{ (Self::BITS / usize::BITS) as usize }>;

    fn to_bitarray(self) -> Self::BaN {
        Self::BaN::new([
            (self & 0xFFFFFFFF) as usize,
            ((self >> 32) & 0xFFFFFFFF) as usize,
        ])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le::<Self>()
    }

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }

    fn shiftr(&mut self, n: u32) {
        *self >>= n;
    }

    fn shiftl(&mut self, n: u32) {
        *self <<= n;
    }

    fn complement(&mut self, mask: usize) {
        let mask = (1 << mask) - 1;
        *self ^= mask;
    }

    fn rev_blocks_2(&mut self) {
        let mut bs = self.swap_bytes().to_le_bytes();

        for b in &mut bs {
            *b = REV_2BIT[*b as usize];
        }
    }
}

impl sealed::KmerStorage for u128 {
    const BITS: usize = u128::BITS as usize;
    type BaN = Ba<{ (Self::BITS / usize::BITS) as usize }>;

    fn to_bitarray(self) -> Self::BaN {
        Self::BaN::new([
            (self & 0xFFFFFFFF) as usize,
            ((self >> 32) & 0xFFFFFFFF) as usize,
            ((self >> 64) & 0xFFFFFFFF) as usize,
            ((self >> 96) & 0xFFFFFFFF) as usize,
        ])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le::<Self>()
    }

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }

    fn shiftr(&mut self, n: u32) {
        *self >>= n;
    }

    fn shiftl(&mut self, n: u32) {
        *self <<= n;
    }

    fn complement(&mut self, mask: usize) {
        let mask = (1 << mask) - 1;
        *self ^= mask;
    }

    fn rev_blocks_2(&mut self) {
        let mut bs = self.swap_bytes().to_le_bytes();

        for b in &mut bs {
            *b = REV_2BIT[*b as usize];
        }
    }
}
