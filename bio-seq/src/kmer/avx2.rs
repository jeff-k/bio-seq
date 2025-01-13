// Copyright 2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! SIMD kmers

use crate::codec::Codec;
use crate::kmer::{sealed, Kmer, KmerStorage};
use crate::{Ba, Bs};
use bitvec::field::BitField;

use core::arch::x86_64::*;
use core::hash::{Hash, Hasher};
use std::simd;

//impl sealed::KmerStorage for simd::Simd<u64, 4> {
impl sealed::KmerStorage for __m256i {
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

impl<const K: usize> From<Kmer<text::Dna, K, _mm256i>> for Kmer<dna::Dna, K, usize> {
    unsafe fn from(src: Kmer<text::Dna, K, _m256i>) -> Self {
        let src = src.storage();

        // multiply each lane of bytes by 3
        // (shift left 2, add self)
        let x3 = _mm256_add_epi64(_mm2_slli_epi64(src, 1), src);

        // now the 4th and 5th bits are our targets
        const MASK: u64 = 0x1818181818181818;

        // extract lanes and pack masked bits
        let p0 = _pext_u64(_mm256_extract_epi64(x3, 0) as u64, MASK);
        let p1 = _pext_u64(_mm256_extract_epi64(x3, 1) as u64, MASK);
        let p2 = _pext_u64(_mm256_extract_epi64(x3, 2) as u64, MASK);
        let p3 = _pext_u64(_mm256_extract_epi64(x3, 3) as u64, MASK);

        let packed = (p0 << 0) | (p1 << 16) | (p2 << 32) | (p3 << 48);

        Kmer::new(packed)
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
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
