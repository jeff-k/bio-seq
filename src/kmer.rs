// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
//use crate::Complement;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;

/// ## kmers
///
/// kmers are encoded sequences with a fixed size that can fit into a register. these are implemented with const generics.
///
/// `k * codec::width` must fit in a `usize` (i.e. 64). for larger kmers use `bigk::kmer`: (todo)

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<C: Codec, const K: usize> {
    bs: usize,
    _p: PhantomData<C>,
}

impl<_C: Codec, const _K: usize> Kmer<_C, _K> {
    pub fn new<C: Codec, const K: usize>(s: &BitSlice) -> Kmer<C, K> {
        assert_eq!(K, s.len() / C::WIDTH as usize);
        Kmer {
            bs: s.load::<usize>(),
            _p: PhantomData,
        }
    }
}

impl<A: Codec, const K: usize> From<&Kmer<A, K>> for usize {
    fn from(kmer: &Kmer<A, K>) -> usize {
        kmer.bs
    }
}

impl<A: Codec, const K: usize> fmt::Display for Kmer<A, K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for chunk in BitArray::<usize, Lsb0>::from(self.bs)[..K * A::WIDTH as usize]
            .chunks(A::WIDTH as usize)
        {
            s.push_str(
                &A::unsafe_from_bits(chunk.load::<u8>())
                    .to_char()
                    .to_string(),
            );
        }
        write!(f, "{}", s,)
    }
}

impl<A: Codec, const K: usize> From<usize> for Kmer<A, K> {
    fn from(_uint: usize) -> Self {
        unimplemented!();
    }
}

/*
impl<const K: usize> Complement for Kmer<K> {
    fn complement(_kmer: Kmer<K>) -> Kmer<K> {
        unimplemented!()
    }
}
*/

#[cfg(test)]
mod tests {
    use crate::codec::amino::Amino;
    use crate::codec::dna::Dna;
    use crate::Seq;
    use std::str::FromStr;
    #[test]
    fn kmer_to_usize() {
        for (kmer, index) in dna!("AACTT").kmers::<2>().zip([0, 4, 13, 15]) {
            assert_eq!(index as usize, (&kmer).into());
        }
    }

    #[test]
    fn amino_kmer_to_usize() {
        for (kmer, index) in amino!("SRY")
            .kmers::<2>()
            .zip([0b001000_011000, 0b010011_001000])
        {
            assert_eq!(index as usize, (&kmer).into());
        }
    }

    #[test]
    fn amino_kmer_iter() {
        for (kmer, target) in amino!("SSLMNHKKL")
            .kmers::<3>()
            .zip(["SSL", "SLM", "LMN", "MNH", "NHK", "HKK", "KKL"])
        {
            assert_eq!(format!("{}", kmer), target);
        }
    }
}
