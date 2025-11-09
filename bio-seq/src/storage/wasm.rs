// Copyright 2025 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! SIMD kmers
//!

use crate::codec::Codec;
use crate::kmer::{sealed, Kmer, KmerStorage};
use crate::{Ba, Bs};
use bitvec::field::BitField;

use core::arch::wasm;
use core::hash::{Hash, Hasher};

impl sealed::KmerStorage for wasm::v128 {
    const BITS: usize = 128;
    type BaN = Ba<2>;

    fn to_bitarray(self) -> Self::BaN {
        let es = self.to_array();
        Ba::<2>::from([es[0] as usize, es[1] as usize])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        debug_assert!(bs.len() == 128);

        Self::from_array([bs[0..64].load_le::<u64>(), bs[64..128].load_le::<u64>()])
    }

    fn mask(&mut self, bits: usize) {
        let mask = if bits >= 128 {
            wasm::v128::splat(u64::MAX)
        } else {
            todo!()
        };

        *self &= mask;
    }
}

impl KmerStorage for wasm::v128 {}

impl<A: Codec, const K: usize> Hash for Kmer<A, K, wasm::v128> {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use core::simd;

    #[test]
    fn kmer_storage_types() {
        todo!()
    }
}
