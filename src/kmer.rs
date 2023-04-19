// Copyright 2021, 2022, 2023 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::seq::{Seq, SeqSlice};
use crate::ParseBioError;
use bitvec::prelude::*;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::str::FromStr;

/// ## Kmers
///
/// Encoded sequences length `k`, fixed at compile time.
///
/// For this implementation `k * codec::WIDTH` must fit in a `usize` (i.e. 64 bits). for larger kmers use `SeqSlice` or
///
/// ```
/// use bio_seq::prelude::*;
///
/// for (amino_kmer, amino_string) in amino!("SSLMNHKKL")
///         .kmers::<3>()
///         .zip(["SSL", "SLM", "LMN", "MNH", "NHK", "HKK", "KKL"])
///     {
///         assert_eq!(amino_kmer, amino_string);
///     }
/// ```
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Copy, Clone)]
pub struct Kmer<C: Codec, const K: usize> {
    pub _p: PhantomData<C>,
    pub bs: usize,
}

impl<A: Codec, const K: usize> Kmer<A, K> {
    /// Kmer uses a run-time test to ensure that K is inside bounds. Compile-time
    /// testing requires nightly.
    #[inline]
    fn check_k() {
        assert!(
            K <= usize::BITS as usize / A::WIDTH as usize,
            "K is too large: it should be <= usize::BITS / A::WIDTH"
        );
    }

    /// Push a base from the right:
    ///
    /// ```
    /// use bio_seq::prelude::*;
    /// use bio_seq::codec::dna::Dna;
    ///
    /// let k = kmer!("ACGAT");
    /// assert_eq!(k.pushr(Dna::T).to_string(), "CGATT");
    /// ```
    pub fn pushr(self, base: A) -> Kmer<A, K> {
        let bs = &BitArray::<usize, Lsb0>::from(base.into() as usize)[..A::WIDTH as usize];
        let ba = &BitArray::<usize, Lsb0>::from(self.bs);

        let mut x: BitVec<usize, Lsb0> = BitVec::new();
        x.extend_from_bitslice(&ba[A::WIDTH as usize..A::WIDTH as usize * K]);
        x.extend_from_bitslice(bs);
        Kmer {
            _p: PhantomData,
            bs: x.load_le::<usize>(),
        }
    }

    /// Push a base from the left
    pub fn pushl(self, base: A) -> Kmer<A, K> {
        let bs = &BitArray::<usize, Lsb0>::from(base.into() as usize)[..A::WIDTH as usize];
        let ba = &BitArray::<usize, Lsb0>::from(self.bs);

        let mut x: BitVec<usize, Lsb0> = BitVec::new();
        x.extend_from_bitslice(bs);
        x.extend_from_bitslice(&ba[..A::WIDTH as usize * K - A::WIDTH as usize]);
        Kmer {
            _p: PhantomData,
            bs: x.load_le::<usize>(),
        }
    }
}

impl<A: Codec, const K: usize> From<usize> for Kmer<A, K> {
    fn from(i: usize) -> Kmer<A, K> {
        Kmer::<A, K>::check_k();
        Kmer {
            _p: PhantomData,
            bs: i,
        }
    }
}

impl<A: Codec, const K: usize> From<&Kmer<A, K>> for usize {
    fn from(kmer: &Kmer<A, K>) -> usize {
        kmer.bs
    }
}

impl<A: Codec, const K: usize> From<Kmer<A, K>> for usize {
    fn from(kmer: Kmer<A, K>) -> usize {
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
                &A::unsafe_from_bits(chunk.load_le::<u8>())
                    .to_char()
                    .to_string(),
            );
        }
        write!(f, "{s}")
    }
}

/// An iterator over all kmers of a sequence with a specified length
pub struct KmerIter<'a, A: Codec, const K: usize> {
    pub slice: &'a SeqSlice<A>,
    pub index: usize,
    pub len: usize,
    pub _p: PhantomData<A>,
}

impl<'a, A: Codec, const K: usize> Iterator for KmerIter<'a, A, K> {
    type Item = Kmer<A, K>;
    fn next(&mut self) -> Option<Kmer<A, K>> {
        let i = self.index;
        if self.index + K > self.len {
            return None;
        }
        self.index += 1;
        Some(Kmer::<A, K>::from(&self.slice[i..i + K]))
    }
}

/// The value of K is included in the hasher state so that
/// `hash(kmer!("AAA")) != hash(kmer!("AAAA"))
impl<A: Codec, const K: usize> Hash for Kmer<A, K> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        K.hash(state);
    }
}

impl<A: Codec, const K: usize> TryFrom<Seq<A>> for Kmer<A, K> {
    type Error = ParseBioError;

    fn try_from(seq: Seq<A>) -> Result<Self, Self::Error> {
        if seq.len() != K {
            Err(ParseBioError {})
        } else {
            Ok(Kmer::<A, K>::from(&seq[0..K]))
        }
    }
}

impl<A: Codec, const K: usize> From<&SeqSlice<A>> for Kmer<A, K> {
    fn from(slice: &SeqSlice<A>) -> Self {
        assert_eq!(K, slice.len());
        Kmer::<A, K>::check_k();
        Kmer {
            _p: PhantomData,
            bs: slice.into(),
        }
    }
}

impl<A: Codec, const K: usize> From<Kmer<A, K>> for String {
    fn from(kmer: Kmer<A, K>) -> Self {
        kmer.to_string()
    }
}

impl<A: Codec, const K: usize> PartialEq<Seq<A>> for Kmer<A, K> {
    fn eq(&self, seq: &Seq<A>) -> bool {
        if seq.len() != K {
            return false;
        }
        &Kmer::<A, K>::from(&seq[..]) == self
    }
}

impl<A: Codec, const K: usize> PartialEq<&str> for Kmer<A, K> {
    fn eq(&self, seq: &&str) -> bool {
        &self.to_string() == seq
    }
}

impl<A: Codec, const K: usize> FromStr for Kmer<A, K> {
    type Err = ParseBioError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != K {
            return Err(ParseBioError {});
        }
        let seq: Seq<A> = Seq::from_str(s)?;
        Kmer::<A, K>::try_from(seq)
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn kmer_to_usize() {
        for (kmer, index) in dna!("AACTT")
            .kmers::<2>()
            .zip([0b00_00, 0b01_00, 0b11_01, 0b11_11])
        {
            assert_eq!(index as usize, (&kmer).into());
        }
    }

    #[test]
    fn pushl_test() {
        let k = kmer!("ACGT");
        let k1 = k.pushl(Dna::G);
        let k2 = k1.pushl(Dna::A);
        let k3 = k2.pushl(Dna::T);
        let k4 = k3.pushl(Dna::C);
        let k5 = k4.pushl(Dna::C);

        assert_eq!(k1, kmer!("GACG"));
        assert_eq!(k2, kmer!("AGAC"));
        assert_eq!(k3, kmer!("TAGA"));
        assert_eq!(k4, kmer!("CTAG"));
        assert_eq!(k5, kmer!("CCTA"));
    }

    #[test]
    fn pushr_test() {
        let k = kmer!("ACGT");
        let k1 = k.pushr(Dna::G);
        let k2 = k1.pushr(Dna::A);
        let k3 = k2.pushr(Dna::T);
        let k4 = k3.pushr(Dna::C);
        let k5 = k4.pushr(Dna::C);

        assert_eq!(k1, kmer!("CGTG"));
        assert_eq!(k2, kmer!("GTGA"));
        assert_eq!(k3, kmer!("TGAT"));
        assert_eq!(k4, kmer!("GATC"));
        assert_eq!(k5, kmer!("ATCC"));
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
    fn big_kmer_shiftr() {
        let mut kmer: Kmer<Dna, 32> = kmer!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTA");
        for base in dna!("TACTATGAGGACGATCAGCACCATAAGAACAAA").into_iter() {
            kmer = kmer.pushr(base);
        }
        assert_eq!(kmer!("ACTATGAGGACGATCAGCACCATAAGAACAAA"), kmer);
    }

    #[test]
    fn big_kmer_shiftl() {
        let mut kmer: Kmer<Dna, 32> = kmer!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTA");
        for base in dna!("GTACTATGAGGACGATCAGCACCATAAGAACAAA").into_iter() {
            kmer = kmer.pushl(base);
        }
        assert_eq!(kmer!("AAACAAGAATACCACGACTAGCAGGAGTATCA"), kmer);
    }

    #[test]
    fn amino_kmer_iter() {
        for (kmer, target) in amino!("SSLMNHKKL")
            .kmers::<3>()
            .zip(["SSL", "SLM", "LMN", "MNH", "NHK", "HKK", "KKL"])
        {
            assert_eq!(kmer, target);
        }
    }

    #[test]
    fn valid_k_check() {
        Kmer::<Dna, 1>::check_k();
        Kmer::<Amino, 1>::check_k();
        Kmer::<Dna, 32>::check_k();
        Kmer::<Amino, 10>::check_k();
    }

    #[test]
    fn eq_functions() {
        assert_eq!(kmer!("ACGT"), dna!("ACGT"));
        assert_ne!(kmer!("ACGT"), dna!("ACGTA"));
        let kmer: Kmer<Iupac, 4> = Kmer::from_str("ACGT").unwrap();
        assert_eq!(kmer, iupac!("ACGT"));
        assert_ne!(kmer, iupac!("NCGT"));
    }

    #[test]
    fn kmer_iter() {
        let seq: Seq<Dna> = dna!("ACTGA");
        let cs: Vec<Kmer<Dna, 3>> = seq.kmers().collect();
        assert_eq!(cs[0], "ACT");
        assert_eq!(cs[1], "CTG");
        assert_eq!(cs[2], "TGA");
        assert_eq!(cs.len(), 3);
    }

    #[test]
    #[should_panic(expected = "K is too large: it should be <= usize::BITS / A::WIDTH")]
    fn invalid_k_check() {
        Kmer::<Dna, 33>::check_k();
    }
    #[test]
    #[should_panic(expected = "K is too large: it should be <= usize::BITS / A::WIDTH")]
    fn invalid_k_check_amino() {
        Kmer::<Amino, 11>::check_k();
    }
    #[test]
    #[should_panic(expected = "K is too large: it should be <= usize::BITS / A::WIDTH")]
    fn invalid_k_check_dna() {
        Kmer::<Dna, 33>::check_k();
    }
}
