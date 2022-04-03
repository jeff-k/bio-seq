// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub use crate::codec::dna::Dna;
use crate::codec::Codec;
use bitvec::prelude::*;
use std::fmt;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<const K: usize> {
    pub bv: BitVec,
}

impl<const _K: usize> Kmer<_K> {
    pub fn new<const K: usize>(s: &BitSlice) -> Kmer<K> {
        assert_eq!(K, s.len() / Dna::WIDTH);
        Kmer {
            bv: BitVec::from(s),
        }
    }

    pub fn as_usize(&self) -> usize {
        assert_eq!(_K < 33, true);
        self.bv.load()
    }
}

impl<const K: usize> fmt::Display for Kmer<K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let w = Dna::WIDTH;
        for i in 0..(self.bv.len() / Dna::WIDTH) {
            s.push_str(&Dna::from_bits(&self.bv[(i * w)..((i * w) + w)].load()).to_string());
        }
        write!(f, "{}", s,)
    }
}

#[cfg(test)]
mod tests {
    use crate::Dna;
    use crate::Seq;
    use std::str::FromStr;
    #[test]
    fn kmer_to_usize() {
        for (kmer, index) in dna!("AACTT").kmers::<2>().zip([0, 4, 13, 15]) {
            assert_eq!(kmer.as_usize(), index as usize);
        }
    }
}
