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
use bitvec::field::BitField;

use core::hash::{Hash, Hasher};
use std::simd;

impl sealed::KmerStorage for simd::Simd<u64, 4> {
    const BITS: usize = 256;
    type BaN = Ba<4>;

    fn to_bitarray(self) -> Self::BaN {
        let es = self.to_array();
        Ba::<4>::from([
            es[0] as usize,
            es[1] as usize,
            es[2] as usize,
            es[3] as usize,
        ])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        debug_assert!(bs.len() == 256);

        Self::from_array([
            bs[0..64].load_le::<u64>(),
            bs[64..128].load_le::<u64>(),
            bs[128..192].load_le::<u64>(),
            bs[192..256].load_le::<u64>(),
        ])
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

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::seq::SeqArray;
    use core::simd;

    #[test]
    fn kmer_storage_types() {
        let s1 = "TTCGTAGCCGCGAACTTACGTAGCCGCGAAAAACGTAGCCGCGAACTTACGTAGCCGCGAAAACGTAGCCGCGAACTTACGTAGCCGCGAAAAAACGTAGCACGCGAACTTACGTAGCCGCGCCCCGG";

        let s2 = "TTCGAAGCCGCGAACTTACGTAGCCGCGAAAAACGTAGCCGCGAACTTACGTAGCCGCGAAAACGTAGCCGCGAACTTACGTAGCCGCGAAAAAACGTAGCACGCGAACTTACGTAGCCGCGCCCCGG";

        assert_eq!(s1.len(), 128);

        assert_eq!(s2.len(), 128);

        let kmer = Kmer::<Dna, 128, simd::Simd<u64, 4>>::from_str(&s1).unwrap();
        let kmer2 = Kmer::<Dna, 128, simd::Simd<u64, 4>>::from_str(&s1).unwrap();

        let seq = Seq::<Dna>::from_str(&s1).unwrap();

        let seq2 = Seq::<Dna>::from_str(&s2).unwrap();

        println!("{kmer}");
        println!("{seq}");
        assert_ne!(seq2, seq);

        assert_eq!(kmer, kmer);

        assert_eq!(kmer, kmer2);

        assert_ne!(kmer, &seq2[..]);
        assert_eq!(kmer, &seq[..]);
    }
}
