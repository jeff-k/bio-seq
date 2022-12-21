// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::{Codec, ParseBioErr};
use crate::{Bound, Seq, SeqSlice, True};
use bitvec::prelude::*;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;

/// ## Kmers
///
/// Encoded sequences of fixed length `k`, known at compile time.
///
/// For this implementation `k * codec::width` must fit in a `usize` (i.e. 64 bits). for larger kmers use `SeqSlice` or
/// `simd::Kmer`
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<C: Codec, const K: usize>
where
    Bound<{ K <= (usize::BITS / C::WIDTH as u32) as usize }>: True,
{
    pub _p: PhantomData<C>,
    pub bs: usize,
}

impl<A: Codec, const K: usize> From<usize> for Kmer<A, K>
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
    fn from(i: usize) -> Kmer<A, K> {
        Kmer {
            _p: PhantomData,
            bs: i,
        }
    }
}

impl<A: Codec, const K: usize> From<&Kmer<A, K>> for usize
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
    fn from(kmer: &Kmer<A, K>) -> usize {
        kmer.bs
    }
}

impl<A: Codec, const K: usize> From<Kmer<A, K>> for usize
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
    fn from(kmer: Kmer<A, K>) -> usize {
        kmer.bs
    }
}

impl<A: Codec, const K: usize> fmt::Display for Kmer<A, K>
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
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
        write!(f, "{s}")
    }
}

/// The value of K is included in the hasher state so that
/// `hash(kmer!("AAA")) != hash(kmer!("AAAA"))
impl<A: Codec, const K: usize> Hash for Kmer<A, K>
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        K.hash(state);
    }
}

impl<A: Codec, const K: usize> TryFrom<Seq<A>> for Kmer<A, K>
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
    type Error = ParseBioErr;

    fn try_from(seq: Seq<A>) -> Result<Self, Self::Error> {
        if seq.len() != K {
            Err(ParseBioErr)
        } else {
            Ok(Kmer::<A, K>::from(&seq[0..K]))
        }
    }
}

impl<A: Codec, const K: usize> From<&SeqSlice<A>> for Kmer<A, K>
where
    Bound<{ K <= (usize::BITS / A::WIDTH as u32) as usize }>: True,
{
    fn from(slice: &SeqSlice<A>) -> Self {
        assert_eq!(K, slice.len());
        Kmer {
            _p: PhantomData,
            bs: slice.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::amino::Amino;
    use crate::codec::dna::Dna;
    use crate::Seq;
    use core::str::FromStr;
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

    #[test]
    fn k_big_enough() {
        // TODO: test instead that this fails to build for K>=11 and K>=33
        for _ in amino!("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS").kmers::<10>() {}
        for _ in
            dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").kmers::<32>()
        {
        }
    }
}
