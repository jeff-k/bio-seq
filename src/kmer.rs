// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::{Codec, ParseBioErr};
use crate::seq::Seq;
use bitvec::prelude::*;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

/// ## kmers
///
/// kmers are encoded sequences with a fixed size that can fit into a register. these are implemented with const generics.
///
/// `k * codec::width` must fit in a `usize` (i.e. 64). for larger kmers use `SeqSlice`

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<C: Codec, const K: usize> {
    pub bs: usize,
    pub _p: PhantomData<C>,
}

impl<A: Codec, const K: usize> Kmer<A, K> {
    pub fn new(s: &BitSlice) -> Kmer<A, K> {
        assert_eq!(K, s.len() / A::WIDTH as usize);
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

impl<A: Codec, const K: usize> Hash for Kmer<A, K> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        K.hash(state);
    }
}

impl<A: Codec, const K: usize> From<usize> for Kmer<A, K> {
    fn from(_uint: usize) -> Kmer<A, K> {
        unimplemented!();
    }
}

pub struct KmerIter<A: Codec, const K: usize> {
    bs: BitBox,
    index: usize,
    len: usize,
    _p: PhantomData<A>,
}

impl<A: Codec, const K: usize> Iterator for KmerIter<A, K> {
    type Item = Kmer<A, K>;
    fn next(&mut self) -> Option<Kmer<A, K>> {
        let k = K * A::WIDTH as usize;
        let i = self.index * A::WIDTH as usize;
        if self.index >= self.len - (K - 1) {
            return None;
        }
        self.index += 1;
        Some(Kmer::<A, K>::new(&self.bs[i..k + i]))
        //        Some(Kmer::<A, K>::new(&self.seq.bv[i..k + i]))
    }
}

impl<A: Codec> Seq<A> {
    pub fn kmers<const K: usize>(self) -> KmerIter<A, K> {
        KmerIter::<A, K> {
            bs: BitBox::from_bitslice(&self.bv),
            index: 0,
            len: self.len(),
            _p: PhantomData,
        }
    }
}

impl<A: Codec, const K: usize> TryFrom<Seq<A>> for Kmer<A, K> {
    type Error = ParseBioErr;

    fn try_from(seq: Seq<A>) -> Result<Self, Self::Error> {
        if seq.len() != K {
            Err(ParseBioErr)
        } else {
            Ok(Kmer::<A, K>::new(&seq.bv[0..(K * A::WIDTH as usize)]))
        }
    }
}

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
