// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub use crate::alphabet::dna::Dna;
use crate::alphabet::Alphabet;
use bitvec::prelude::*;
use std::fmt;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<const K: u8> {
    bv: BitVec,
}

impl<const _K: u8> Kmer<_K> {
    pub fn new<const K: u8>(s: &BitSlice) -> Kmer<K> {
        assert_eq!(K as usize, s.len() / Dna::width());
        Kmer {
            bv: BitVec::from(s),
        }
    }
}

impl<const K: u8> fmt::Display for Kmer<K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let w = Dna::width();
        for i in 0..(self.bv.len() / Dna::width()) {
            s.push_str(&Dna::from_bits(&self.bv[(i * w)..((i * w) + w)]).to_string());
        }
        write!(f, "{}", s,)
    }
}
