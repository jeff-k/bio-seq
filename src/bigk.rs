// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;

/// Kmers for large `K`
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<'a, A: Codec, const K: usize> {
    pub bv: &'a BitSlice,
    _p: PhantomData<A>,
}

impl<'a, _A: Codec, const _K: usize> Kmer<'a, _A, _K> {
    pub fn new<A: Codec, const K: usize>(s: &'a BitSlice) -> Kmer<'a, A, K> {
        Kmer {
            bv: s,
            _p: PhantomData,
        }
    }
}

impl<A: Codec, const K: usize> fmt::Display for Kmer<'_, A, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for c in self.bv.chunks_exact(A::WIDTH.into()) {
            s.push_str(&A::unsafe_from_bits(c.load()).to_char().to_string());
        }

        write!(f, "{}", s,)
    }
}
