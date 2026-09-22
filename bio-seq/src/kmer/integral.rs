use crate::kmer::{REV_2BIT, sealed};
use bitvec::field::BitField;

use crate::{Ba, Bs, codec::Codec};

const MASK4: usize = usize::MAX / 0x11;
//const MASK4_32: u32 = 0x0f0f_0f0f;
const MASK4_64: u64 = 0x0f0f_0f0f_0f0f_0f0f;
const MASK4_128: u128 = 0x0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f;

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

    fn rev_blocks<A: Codec, const K: usize>(&mut self) {
        match A::BITS {
            1 => *self = self.reverse_bits(),
            2 => {
                let mut bs = self.swap_bytes().to_le_bytes();

                for b in &mut bs {
                    *b = REV_2BIT[*b as usize];
                }

                *self = Self::from_le_bytes(bs);
            }
            4 => {
                let bs = self.swap_bytes();
                *self = ((bs >> 4) & MASK4) | ((bs & MASK4) << 4);
            }
            _ => todo!(),
        }
    }
}

impl sealed::KmerStorage for u64 {
    const BITS: usize = u64::BITS as usize;

    type BaN = Ba<{ (Self::BITS / usize::BITS) as usize }>;

    #[cfg(target_pointer_width = "64")]
    fn to_bitarray(self) -> Self::BaN {
        Self::BaN::new([self.try_into().unwrap()])
    }

    #[cfg(target_pointer_width = "32")]
    fn to_bitarray(self) -> Self::BaN {
        Self::BaN::new([
            (self & 0xFFFF_FFFF) as usize,
            ((self >> 32) & 0xFFFF_FFFF) as usize,
        ])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le::<Self>()
    }

    fn shiftr(&mut self, n: u32) {
        *self >>= n;
    }

    fn complement(&mut self, mask: usize) {
        if mask >= Self::BITS as usize {
            *self ^= Self::MAX;
        } else {
            *self ^= (1 << mask) - 1;
        }
    }

    fn rev_blocks<A: Codec, const K: usize>(&mut self) {
        match A::BITS {
            1 => *self = self.reverse_bits(),
            2 => {
                let mut bs = self.swap_bytes().to_le_bytes();

                for b in &mut bs {
                    *b = REV_2BIT[*b as usize];
                }

                *self = Self::from_le_bytes(bs);
            }
            4 => {
                let bs = self.swap_bytes();
                *self = ((bs >> 4) & MASK4_64) | ((bs & MASK4_64) << 4);
            }
            _ => todo!(),
        }
    }
}

impl sealed::KmerStorage for u128 {
    const BITS: usize = u128::BITS as usize;
    type BaN = Ba<{ (Self::BITS / usize::BITS) as usize }>;

    #[cfg(target_pointer_width = "64")]
    #[expect(
        clippy::cast_possible_truncation,
        reason = "split the u128 into its low and high 64-bit words"
    )]
    fn to_bitarray(self) -> Self::BaN {
        Self::BaN::new([self as usize, (self >> 64) as usize])
    }

    #[cfg(target_pointer_width = "32")]
    fn to_bitarray(self) -> Self::BaN {
        Self::BaN::new([
            (self & 0xFFFF_FFFF) as usize,
            ((self >> 32) & 0xFFFF_FFFF) as usize,
            ((self >> 64) & 0xFFFF_FFFF) as usize,
            ((self >> 96) & 0xFFFF_FFFF) as usize,
        ])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le::<Self>()
    }

    fn shiftr(&mut self, n: u32) {
        *self >>= n;
    }

    fn complement(&mut self, mask: usize) {
        if mask >= Self::BITS as usize {
            *self ^= Self::MAX;
        } else {
            *self ^= (1 << mask) - 1;
        }
    }

    fn rev_blocks<A: Codec, const K: usize>(&mut self) {
        match A::BITS {
            1 => *self = self.reverse_bits(),
            2 => {
                let mut bs = self.swap_bytes().to_le_bytes();

                for b in &mut bs {
                    *b = REV_2BIT[*b as usize];
                }

                *self = Self::from_le_bytes(bs);
            }
            4 => {
                let bs = self.swap_bytes();
                *self = ((bs >> 4) & MASK4_128) | ((bs & MASK4_128) << 4);
            }
            _ => todo!(),
        }
    }
}
