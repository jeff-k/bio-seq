// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Encoded sequences of static length
//!
//! Generally, the underlying storage type of `Kmer` should lend itself to optimisation. The default `Kmer` instance is packed into a `usize`, which can be efficiently `Copy`ed on the stack.
//!
//! `k * codec::BITS` must fit in the storage type, e.g. `usize` (64 bits).
//!
//! ```
//! use bio_seq::prelude::*;
//!
//! for (amino_kmer, amino_string) in Seq::<Amino>::try_from("SSLMNHKKL").unwrap()
//!         .kmers::<3>()
//!         .zip(["SSL", "SLM", "LMN", "MNH", "NHK", "HKK", "KKL"])
//!     {
//!         assert_eq!(amino_kmer, amino_string);
//!     }
//! ```
//!
//! Kmers can be copied from other sequence types:
//!
//! ```
//! # use bio_seq::prelude::*;
//! let kmer: Kmer<Dna, 8> = dna!("AGTTGGCA").into();
//! ```

// permit truncations that may happen on 32-bit platforms which are unsupported anyway
#![allow(clippy::cast_possible_truncation)]

use crate::codec::Codec;
use crate::prelude::{Complement, ParseBioError, ReverseComplement};
use crate::seq::{Seq, SeqArray, SeqSlice};
use crate::{Ba, Bs};
use bitvec::field::BitField;
use bitvec::view::BitView;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::ops::Deref;
use core::ptr;
use core::str::FromStr;

#[cfg(feature = "simd")]
pub mod simd;

#[cfg(feature = "serde")]
use serde_derive::{Deserialize, Serialize};

mod sealed {
    use crate::Bs;

    pub trait KmerStorage: Copy + Clone + PartialEq {
        const BITS: usize;
        type BaN: AsRef<Bs> + AsMut<Bs>;

        fn to_bitarray(self) -> Self::BaN;
        fn from_bitslice(bs: &Bs) -> Self;

        //        fn rotate_left(self, n: u32) -> Self;
        //        fn rotate_right(self, n: u32) -> Self;

        fn mask(&mut self, bits: usize);
    }
}

pub trait KmerStorage: sealed::KmerStorage {}

/*
impl sealed::KmerStorage for u32 {
    const BITS: usize = u32::BITS as usize;
    type BaN = Ba<1>;

    fn to_bitarray(self) -> Ba<1> {
        Self::BaN::new([self as usize])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        debug_assert!(bs.len() >= 32, "Target SeqSlice data less than 32 bits");
        bs[..32].load_le()
    }

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }
}

impl KmerStorage for u32 {}
*/

impl sealed::KmerStorage for usize {
    const BITS: usize = usize::BITS as usize;
    type BaN = Ba<1>;

    fn to_bitarray(self) -> Ba<1> {
        Self::BaN::new([self])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le()
    }

    /*
        fn rotate_right(self, n: u32) -> Self {
            self.rotate_right(n)
        }
        fn rotate_left(self, n: u32) -> Self {
            self.rotate_left(n)
        }
    */

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }
}

impl KmerStorage for usize {}

impl sealed::KmerStorage for u64 {
    const BITS: usize = u64::BITS as usize;
    type BaN = Ba<1>;

    fn to_bitarray(self) -> Ba<1> {
        // depends on global assertion that usize::BITS == 64
        Self::BaN::new([self.try_into().unwrap()])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le::<Self>()
    }

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }
}

impl KmerStorage for u64 {}

impl sealed::KmerStorage for u128 {
    const BITS: usize = u128::BITS as usize;
    type BaN = Ba<2>;

    fn to_bitarray(self) -> Ba<2> {
        Self::BaN::new([self as usize, (self >> 64) as usize])
    }

    fn from_bitslice(bs: &Bs) -> Self {
        bs.load_le::<Self>()
    }

    fn mask(&mut self, bits: usize) {
        *self &= (1 << bits) - 1;
    }
}

impl KmerStorage for u128 {}

/// By default k-mers are backed by `usize` and `Codec::BITS` * `K` must be <= 64 on 64-bit platforms
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[repr(transparent)]
pub struct Kmer<C: Codec, const K: usize, S: KmerStorage = usize> {
    pub _p: PhantomData<C>,
    pub bs: S,
}

/*
impl<A: Codec, const K: usize, S1: KmerStorage, S2: KmerStorage> PartialEq<Kmer<A, K, S2>> for Kmer<A, K, S1> {
    fn eq(&self, other: &Kmer<A, K, S2>) -> bool {
        self.bs.to_bitarray() == other.bs.to_bitarray()
    }
}
*/

impl<A: Codec, const K: usize, S: KmerStorage> Kmer<A, K, S> {
    // This error message can be formatted with constants in nightly (const_format)
    const _ASSERT_K: () = assert!(
        K * A::BITS as usize <= S::BITS,
        "`KmerStorage` not large enough for `Kmer`",
    );

    const _ASSERT_K_NONZERO: () = assert!(K > 0, "`K` must be greater than 0");

    const BITS: usize = K * A::BITS as usize;

    pub fn len(&self) -> usize {
        K
    }

    pub fn is_empty(&self) -> bool {
        // This is recommended by clippy since we have `len`
        // Kmers are never empty if K > 0
        false
    }

    pub fn rotated_left(&self, n: u32) -> Self {
        let n: usize = (n as usize % K) * A::BITS as usize;
        let mut ba = self.bs.to_bitarray();
        let bs: &mut Bs = ba.as_mut();
        bs[..Self::BITS].rotate_left(n);

        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(&bs[..Self::BITS]),
        }
    }

    pub fn rotated_right(&self, n: u32) -> Self {
        let n: usize = (n as usize % K) * A::BITS as usize;
        let mut ba = self.bs.to_bitarray();
        let bs: &mut Bs = ba.as_mut();
        bs[..Self::BITS].rotate_right(n);

        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(&bs[..Self::BITS]),
        }
    }

    /// Shift bases to the right and push a base onto the end.
    ///
    /// ```
    /// use bio_seq::prelude::*;
    /// use bio_seq::codec::dna::Dna;
    ///
    /// let k = kmer!("ACGAT");
    /// assert_eq!(k.pushr(Dna::T).to_string(), "CGATT");
    /// ```
    pub fn pushr(self, base: A) -> Self {
        let mut ba = self.rotated_left(1).bs.to_bitarray();
        let bs: &mut Bs = ba.as_mut();

        let start = Self::BITS - A::BITS as usize;
        let end = start + A::BITS as usize;

        bs[start..end].store(base.to_bits());

        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(bs),
        }
    }

    /// Push a base from the left
    pub fn pushl(self, base: A) -> Self {
        let mut ba = self.rotated_right(1).bs.to_bitarray();
        let bs: &mut Bs = ba.as_mut();

        bs[..A::BITS as usize].store(base.to_bits());

        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(bs),
        }
    }

    /// Create Kmer from sequence without checking length
    pub fn unsafe_from_seqslice(seq: &SeqSlice<A>) -> Self {
        debug_assert!(K == seq.len(), "K != seq.len()");
        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(&seq.bs),
        }
    }
    /*
        /// Iterate through all bases of a Kmer
        pub fn iter(self) -> KmerBases<A, K> {
            KmerBases {
                _p: PhantomData,
                bits: Ba::from(self.bs),
                index: 0,
            }
        }
    */
    /*
    /// tail
    pub fn tail(self, base: A) -> Seq<A> {
        let bv: Bv::from(self.bs);

        Seq {
            _p: PhantomData,
            bv: bv[A::BITS ..],
        }
    }
    */
}

impl<A: Codec, const K: usize> From<usize> for Kmer<A, K, usize> {
    fn from(i: usize) -> Kmer<A, K> {
        Kmer {
            _p: PhantomData,
            bs: i,
        }
    }
}

/*
impl<A: Codec, const K: usize, S: KmerStorage> From<&SeqSlice<A>> for Kmer<A, K, S> {
    fn from(seq: &SeqSlice<A>) -> Self {
        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(&seq.bs),
        }
    }
}
*/

impl<A: Codec, const K: usize> From<&Kmer<A, K, usize>> for usize {
    fn from(kmer: &Kmer<A, K, usize>) -> usize {
        kmer.bs
    }
}

/*
impl<A: Codec, const K: usize> Deref for SeqArray<A, K, 1> {
    type Target = Kmer<A, K>;

    fn deref(&self) -> &Self::Target {
        let bs: usize = &self.ba.view_bits()[0..(K * A::BITS as usize)];
        unsafe { &*(bs *const Kmer<A, K>) }
    }
}
*/

impl<A: Codec, const K: usize> Deref for Kmer<A, K> {
    type Target = SeqSlice<A>;

    fn deref(&self) -> &Self::Target {
        let bs: &Bs = &self.bs.view_bits()[0..(K * A::BITS as usize)];
        let bs: *const Bs = ptr::from_ref::<Bs>(bs);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec, const K: usize> AsRef<SeqSlice<A>> for Kmer<A, K> {
    fn as_ref(&self) -> &SeqSlice<A> {
        self
    }
}

impl<A: Codec, const K: usize, S: KmerStorage> fmt::Display for Kmer<A, K, S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let ba = self.bs.to_bitarray();
        let bs: &Bs = &ba.as_ref()[0..(K * A::BITS as usize)];

        bs.chunks(A::BITS as usize).for_each(|chunk| {
            s.push(A::unsafe_from_bits(chunk.load_le::<u8>()).to_char());
        });
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

impl<A: Codec, const K: usize, S: KmerStorage> Kmer<A, K, S> {
    fn unsafe_from(seq: &SeqSlice<A>) -> Self {
        Kmer {
            _p: PhantomData,
            bs: S::from_bitslice(&seq.bs),
        }
    }
}

impl<A: Codec, const K: usize> Iterator for KmerIter<'_, A, K> {
    type Item = Kmer<A, K>;
    fn next(&mut self) -> Option<Kmer<A, K>> {
        let i = self.index;
        if self.index + K > self.len {
            return None;
        }
        self.index += 1;
        Some(Kmer::<A, K>::unsafe_from(&self.slice[i..i + K]))
    }
}

/*
/// An iterator over the bases of a kmer
pub struct KmerBases<A: Codec, const K: usize, S: KmerStorage> {
    pub _p: PhantomData<A>,
    pub bits: S::BaN,
    pub index: usize,
}

/// An iterator over the bases of a kmer
impl<A: Codec, const K: usize> Iterator for KmerBases<A, K> {
    type Item = A;

    fn next(&mut self) -> Option<A> {
        let i = self.index * A::BITS as usize;
        if self.index >= K {
            return None;
        }
        self.index += 1;
        let chunk = &self.bits[i..i + (A::BITS as usize)];
        Some(A::unsafe_from_bits(chunk.load_le::<u8>()))
    }
}
*/

/// ```
/// use bio_seq::prelude::*;
/// use std::hash::{Hash, Hasher, DefaultHasher};
///
/// let mut hasher1 = DefaultHasher::new();
/// kmer!("AAA").hash(&mut hasher1);
/// let hash1 = hasher1.finish();
///
/// let mut hasher2 = DefaultHasher::new();
/// kmer!("AAAA").hash(&mut hasher2);
/// let hash2 = hasher2.finish();
///
/// assert_ne!(hash1, hash2);
/// ```
impl<A: Codec, const K: usize> Hash for Kmer<A, K> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let bs: &SeqSlice<A> = self.as_ref();
        bs.hash(state);
        //K.hash(state);
    }
}

impl<A: Codec, const K: usize, S: KmerStorage> TryFrom<&SeqSlice<A>> for Kmer<A, K, S> {
    type Error = ParseBioError;

    fn try_from(seq: &SeqSlice<A>) -> Result<Self, Self::Error> {
        if seq.len() == K {
            Ok(Kmer::<A, K, S>::unsafe_from(&seq[0..K]))
        } else {
            Err(ParseBioError::MismatchedLength(K, seq.len()))
        }
    }
}

impl<A: Codec, const K: usize> TryFrom<Seq<A>> for Kmer<A, K> {
    type Error = ParseBioError;

    fn try_from(seq: Seq<A>) -> Result<Self, Self::Error> {
        Self::try_from(seq.as_ref())
    }
}

/*
impl<A: Codec, const K: usize> From<Kmer<A, K>> for String {
    fn from(kmer: Kmer<A, K>) -> Self {
        kmer.to_string()
    }
}
*/

impl<A: Codec, const K: usize, S: KmerStorage> PartialEq<SeqArray<A, K, 1>> for Kmer<A, K, S> {
    fn eq(&self, seq: &SeqArray<A, K, 1>) -> bool {
        if seq.len() != K {
            return false;
        }
        &Kmer::<A, K, S>::unsafe_from(seq.as_ref()) == self
    }
}

impl<A: Codec, const K: usize, S: KmerStorage> PartialEq<&SeqArray<A, K, 1>> for Kmer<A, K, S> {
    fn eq(&self, seq: &&SeqArray<A, K, 1>) -> bool {
        if seq.len() != K {
            return false;
        }
        &Kmer::<A, K, S>::unsafe_from(seq.as_ref()) == self
    }
}

impl<A: Codec, const K: usize> PartialEq<Seq<A>> for Kmer<A, K> {
    fn eq(&self, seq: &Seq<A>) -> bool {
        if seq.len() != K {
            return false;
        }
        &Kmer::<A, K>::unsafe_from(seq.as_ref()) == self
    }
}

impl<A: Codec, const K: usize, S: KmerStorage> PartialEq<SeqSlice<A>> for Kmer<A, K, S> {
    fn eq(&self, seq: &SeqSlice<A>) -> bool {
        if seq.len() != K {
            return false;
        }
        &Kmer::<A, K, S>::unsafe_from(seq) == self
    }
}

impl<A: Codec, const K: usize, S: KmerStorage> PartialEq<&SeqSlice<A>> for Kmer<A, K, S> {
    fn eq(&self, seq: &&SeqSlice<A>) -> bool {
        if seq.len() != K {
            return false;
        }
        &Kmer::<A, K, S>::unsafe_from(seq) == self
    }
}

impl<A: Codec, const K: usize> PartialEq<&str> for Kmer<A, K> {
    fn eq(&self, seq: &&str) -> bool {
        &self.to_string() == seq
    }
}

impl<A: Codec, const K: usize, S: KmerStorage> FromStr for Kmer<A, K, S> {
    type Err = ParseBioError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != K {
            return Err(ParseBioError::MismatchedLength(K, s.len()));
        }
        let seq: Seq<A> = Seq::from_str(s)?;
        Kmer::<A, K, S>::try_from(seq.as_ref())
    }
}

impl<A: Codec, const K: usize> From<Kmer<A, K>> for Seq<A> {
    fn from(kmer: Kmer<A, K>) -> Self {
        let mut seq: Seq<A> = Seq::with_capacity(K);
        seq.extend(kmer.iter());
        seq
    }
}

impl<A: Codec + Complement, const K: usize> ReverseComplement for Kmer<A, K> {
    type Output = Self;

    fn revcomp(&self) -> Self {
        let seq: Seq<A> = (*self).into();
        let x: usize = seq.revcomp().into();
        Self::from(x)
    }
}

/// Convenient compile time kmer constructor
///
/// This is a wrapper for the `dna!` macro that returns a `Kmer`:
/// ```
/// # use bio_seq::prelude::*;
/// let kmer: Kmer<Dna, 8> = dna!("ACGTACGT").into();
/// ```
#[macro_export]
macro_rules! kmer {
    ($seq:expr) => {
        Kmer::<Dna, { $seq.len() }>::unsafe_from_seqslice(dna!($seq))
    };
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::seq::SeqArray;

    #[test]
    fn kmer_to_usize() {
        let s: &'static SeqSlice<Dna> = dna!("AACTT");
        println!("{}", s.to_string());

        for (kmer, index) in s.kmers::<2>().zip([0b00_00, 0b01_00, 0b11_01, 0b11_11]) {
            println!("{kmer}");
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

        println!("{}", k1);
        assert_eq!(k1, kmer!("CGTG"));
        assert_eq!(k2, kmer!("GTGA"));
        assert_eq!(k3, kmer!("TGAT"));
        assert_eq!(k4, kmer!("GATC"));
        assert_eq!(k5, kmer!("ATCC"));
    }

    #[test]
    fn amino_kmer_to_usize() {
        for (kmer, index) in Seq::<Amino>::try_from("SRY")
            .unwrap()
            .kmers::<2>()
            .zip([0b001000_011000, 0b010011_001000])
        {
            assert_eq!(index as usize, usize::from(&kmer));
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
        for (kmer, target) in Seq::<Amino>::try_from("SSLMNHKKL")
            .unwrap()
            .kmers::<3>()
            .zip(["SSL", "SLM", "LMN", "MNH", "NHK", "HKK", "KKL"])
        {
            assert_eq!(kmer, target);
        }
    }

    #[test]
    fn test_rotations() {
        let kmer = Kmer::try_from(dna!("ACTGCGATG")).unwrap();

        for (shift, rotation) in vec![
            "ACTGCGATG",
            "CTGCGATGA",
            "TGCGATGAC",
            "GCGATGACT",
            "CGATGACTG",
            "GATGACTGC",
            "ATGACTGCG",
            "TGACTGCGA",
            "GACTGCGAT",
            "ACTGCGATG",
            "CTGCGATGA",
            "TGCGATGAC",
        ]
        .into_iter()
        .enumerate()
        {
            //            println!("{} {} {}", shift, kmer.rotated_left(shift as u32), rotation);
            assert_eq!(kmer.rotated_left(shift as u32), rotation);
        }

        for (shift, rotation) in vec![
            "ACTGCGATG",
            "GACTGCGAT",
            "TGACTGCGA",
            "ATGACTGCG",
            "GATGACTGC",
            "CGATGACTG",
            "GCGATGACT",
            "TGCGATGAC",
            "CTGCGATGA",
            "ACTGCGATG",
            "GACTGCGAT",
        ]
        .into_iter()
        .enumerate()
        {
            //            println!("{} {} {}", shift, kmer.rotated_right(shift as u32), rotation);
            assert_eq!(kmer.rotated_right(shift as u32), rotation);
        }

        let kmer: Kmer<Dna, 8> = Kmer::try_from(dna!("ACTGCGAT")).unwrap().rotated_left(1);

        assert_ne!(kmer.to_string(), "ACTGCGAT");
        assert_eq!(kmer.to_string(), "CTGCGATA");

        let kmer: Kmer<Dna, 8> = kmer!("ACTGCGAT").rotated_right(1);

        assert_ne!(kmer.to_string(), "ACTGCGAT");
        assert_eq!(kmer.to_string(), "TACTGCGA");

        let kmer: Kmer<Dna, 9> = Kmer::from_str("ACTGCGATG").unwrap().rotated_left(0);

        assert_eq!(kmer.to_string(), "ACTGCGATG");
        assert_ne!(kmer.to_string(), "ACTGCGATGA");

        let kmer: Kmer<Dna, 9> = Kmer::try_from(dna!("ACTGCGATG")).unwrap().rotated_right(0);

        assert_eq!(kmer.to_string(), "ACTGCGATG");
        assert_ne!(kmer.to_string(), "ACTGCGATGA");

        let kmer: Kmer<Dna, 9> = Kmer::from_str("ACTGCGATG").unwrap().rotated_left(9 * 3307);

        assert_eq!(kmer.to_string(), "ACTGCGATG");
        assert_ne!(kmer.to_string(), "ACTGCGATGA");

        let kmer: Kmer<Dna, 9> = kmer!("ACTGCGATG").rotated_right(9 * 3307);

        assert_eq!(kmer.to_string(), "ACTGCGATG");
        assert_ne!(kmer.to_string(), "ACTGCGATGA");
    }

    #[test]
    fn eq_functions() {
        assert_eq!(kmer!("ACGT"), dna!("ACGT"));

        // this should be a compiler error:
        // assert_ne!(kmer!("ACGT"), dna!("ACGTA"));

        let kmer: Kmer<Iupac, 4> = Kmer::from_str("ACGT").unwrap();
        assert_eq!(kmer, iupac!("ACGT"));
        assert_ne!(kmer, iupac!("NCGT"));
    }

    #[test]
    fn kmer_iter() {
        //let seq = dna!("ACTGA");
        let cs: Vec<Kmer<Dna, 3>> = dna!("ACTGA").kmers().collect();
        assert_eq!(cs[0], "ACT");
        assert_eq!(cs[1], "CTG");
        assert_eq!(cs[2], "TGA");
        assert_eq!(cs.len(), 3);
    }

    #[test]
    fn k_check() {
        let _kmer = Kmer::<Dna, 32>::from(0);
        let _kmer = Kmer::<Amino, 10>::from(0);
        let _kmer = Kmer::<Iupac, 14>::from(0);
    }

    #[test]
    fn kmer_revcomp() {
        assert_eq!(kmer!("ACGT"), kmer!("ACGT").revcomp());
        assert_ne!(kmer!("GTCGTA"), kmer!("TACGAC"));

        assert_eq!(
            kmer!("ATCGCTATCGATCTGATCGTATATAATATATA"),
            kmer!("TATATATTATATACGATCAGATCGATAGCGAT").revcomp()
        );
    }

    #[test]
    fn kmer_deref() {
        let kmer: Kmer<Dna, 3> = kmer!("ACG");
        let seq: &SeqSlice<Dna> = &kmer;

        assert_eq!(*seq, *kmer);
        assert_eq!(seq.to_string(), "ACG");

        let kmer: Kmer<Dna, 8> = kmer!("AAAAAAAA");
        let seq: &SeqSlice<Dna> = &kmer;

        assert_eq!(seq.to_string(), "AAAAAAAA");

        let kmer: Kmer<Dna, 8> = kmer!("TTTTTTTT");
        let seq: &SeqSlice<Dna> = &kmer;

        assert_eq!(seq.to_string(), "TTTTTTTT");

        let kmer: Kmer<Dna, 32> = kmer!("TTTTTTTTTTTTTTTTAGCTAGCTAGCTAGCT");
        let seq: &SeqSlice<Dna> = &kmer;

        assert_eq!(seq.to_string(), "TTTTTTTTTTTTTTTTAGCTAGCTAGCTAGCT");
    }

    #[test]
    fn kmer_as_ref() {
        let kmer: Kmer<Dna, 4> = kmer!("ACGT");
        let seq: &SeqSlice<Dna> = &kmer.as_ref();

        assert_eq!(seq.to_string(), "ACGT");

        let kmer: Kmer<Dna, 32> = kmer!("TTTTTTTTTTTTTTTTAGCTAGCTAGCTAGCT");
        let seq: &SeqSlice<Dna> = &kmer.as_ref();

        assert_eq!(seq.to_string(), "TTTTTTTTTTTTTTTTAGCTAGCTAGCTAGCT");
    }

    #[test]
    fn kmer_storage_types() {
        let s1 = "AACGTAGCCGCGAACTTACGTAGCCGCGAAAA";
        let s2 = "AACGTAGCCGCGAACTTACGTAGCCGCGAAA";
        let s3 = "ACGTAGCCGCGAACTTACGTAGCCGCGAAAA";

        let s4 = "AACGTAGCCGCGAACTTACGTAGCCGCGAAAAAACGTAGCCGCGAACTTACGTAGCCGCGAAAA";
        let s5 = "AACGTAGCCGCGAACTTACGTAGCCGCGAAAAAACGTAGCCGCGAACTTACGTAGCCGCGAAAAA";

        assert_eq!(s1.len(), 32);
        assert_eq!(s2.len(), 31);
        assert_eq!(s3.len(), 31);
        assert_eq!(s4.len(), 64);
        assert_eq!(s5.len(), 65);

        let kmer1_64 = Kmer::<Dna, 32, u64>::from_str(&s1).unwrap();
        let kmer2_64 = Kmer::<Dna, 31, u64>::from_str(&s2).unwrap();
        let kmer3_64 = Kmer::<Dna, 31, u64>::from_str(&s3).unwrap();

        let kmer1 = Kmer::<Dna, 32>::from_str(&s1).unwrap();
        let kmer2 = Kmer::<Dna, 31>::from_str(&s2).unwrap();
        let kmer3 = Kmer::<Dna, 31>::from_str(&s3).unwrap();

        let kmer4_128 = Kmer::<Dna, 64, u128>::from_str(&s4).unwrap();

        let seq5: Seq<Dna> = s5.try_into().unwrap();

        assert_eq!(kmer4_128, &seq5[..64]);
        assert_ne!(kmer4_128, &seq5[1..]);

        assert_eq!(kmer1, &seq5[..32]);
        assert_eq!(kmer1, &seq5[32..64]);

        assert_eq!(kmer1_64, &seq5[..32]);
        assert_eq!(kmer1_64, &seq5[32..64]);

        assert_ne!(kmer1, &seq5[..31]);
        assert_ne!(kmer1, &seq5[32..]);

        assert_ne!(kmer1_64, &seq5[1..33]);
        assert_ne!(kmer1_64, &seq5[33..]);

        assert_eq!(kmer4_128, kmer4_128);
        assert_eq!(kmer1_64, kmer1_64);
        assert_eq!(kmer2, kmer2);

        assert_ne!(kmer2, kmer3);
        assert_ne!(kmer2_64, kmer3_64);
        // PartialEq is not implemented for different storgage types
        /*
                assert_ne!(kmer2, kmer3_64);
                assert_eq!(kmer2, kmer2_64);
                assert_eq!(kmer1, kmer1_64);
        */
    }

    #[test]
    fn try_from_seq() {
        let seq: Seq<Dna> = Seq::try_from("ACACACACACACGT").unwrap();
        assert_eq!(
            Kmer::<Dna, 8>::try_from(&seq[..8]).unwrap().to_string(),
            "ACACACAC"
        );
        assert_eq!(
            Kmer::<Dna, 8>::try_from(&seq[1..9]).unwrap().to_string(),
            "CACACACA"
        );

        let err: Result<Kmer<Dna, 8>, ParseBioError> = Kmer::try_from(&seq[2..9]);
        assert_eq!(err, Err(ParseBioError::MismatchedLength(8, 7)));

        let err: Result<Kmer<Dna, 8>, ParseBioError> = Kmer::try_from(seq);
        assert_eq!(err, Err(ParseBioError::MismatchedLength(8, 14)));

        let seq: Seq<Dna> = Seq::try_from("ACACACACACACGT").unwrap();

        assert_eq!(
            Kmer::<Dna, 14>::try_from(seq).unwrap().to_string(),
            "ACACACACACACGT"
        );
    }
}
