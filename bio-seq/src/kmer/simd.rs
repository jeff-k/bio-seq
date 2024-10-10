// Copyright 2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! SIMD kmers
//!

use crate::codec::Codec;
use crate::kmer::{sealed, Kmer, KmerStorage};
//use crate::seq::{SeqSlice};
use crate::{Ba, Bs};

use core::hash::{Hash, Hasher};
use std::simd;

impl sealed::KmerStorage for simd::Simd<u64, 4> {
    const BITS: usize = 256;
    type BaN = Ba<4>;

    fn to_bitarray(self) -> Self::BaN {
        todo!()
    }

    fn from_bitslice(_bs: &Bs) -> Self {
        todo!()
    }

    fn mask(&mut self, bits: usize) {
        let mask = if bits >= 256 {
            simd::u64x4::splat(u64::MAX)
        } else {
            todo!()
        };

        *self &= mask;
    }
}

impl KmerStorage for simd::Simd<u64, 4> {}

impl<A: Codec, const K: usize> Hash for Kmer<A, K, simd::Simd<u64, 4>> {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        todo!()
    }
}
