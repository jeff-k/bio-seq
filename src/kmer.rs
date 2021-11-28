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
    pub bv: BitVec<Msb0, u8>,
}

impl<const _K: usize> Kmer<_K> {
    pub fn new<const K: usize>(s: &BitSlice::<Msb0, u8>) -> Kmer<K> {
        assert_eq!(K, s.len() / Dna::WIDTH);
        Kmer {
            bv: BitVec::<Msb0, u8>::from(s),
        }
    }
}

impl<const K: usize> fmt::Display for Kmer<K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let w = Dna::WIDTH;
        for i in 0..(self.bv.len() / Dna::WIDTH) {
            s.push_str(&Dna::from_bits(&self.bv[(i * w)..((i * w) + w)]).to_string());
        }
        write!(f, "{}", s,)
    }
}
