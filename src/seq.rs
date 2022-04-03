// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

/*!
Sequences of bio alphabet characters. Slicable, Boxable, Iterable.
!*/

pub use crate::codec::dna::Dna;
//pub use crate::codec::iupac::Iupac;
pub use crate::codec::Codec;
pub use crate::kmer::Kmer;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;
//use std::ops::{BitAnd, BitOr};
pub use std::str::FromStr;

/// A sequence of bit packed characters
pub struct Seq<A: Codec> {
    pub bv: BitVec,
    _p: PhantomData<A>,
}

pub struct SeqSlice<'a, A: Codec> {
    pub bv: &'a BitSlice,
    _p: PhantomData<A>,
}

#[derive(Debug, Clone)]
pub struct ParseSeqErr;

impl fmt::Display for ParseSeqErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "could not parse sequence")
    }
}

impl<A: Codec> Seq<A> {
    /// Pack binary representations into a bitvector
    pub fn from_vec(vec: Vec<A>) -> Self {
        let mut bv: BitVec = BitVec::new();
        for b in vec.iter() {
            bv.extend_from_bitslice(&b.to_bits().view_bits::<Lsb0>()[..A::WIDTH]);
        }
        Seq {
            bv,
            _p: PhantomData,
        }
    }

    pub fn raw(&self) -> &[usize] {
        self.bv.as_raw_slice()
    }
}

impl<A: Codec> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for c in self.bv.chunks_exact(A::WIDTH) {
            s.push_str(&A::from_bits(&c.load()).to_string());
        }
        write!(
            f,
            "{}",
            s,
            //&self.bv[(4 * w)..((4 * w) + w)][2],
        )
    }
}

impl<A: Codec> FromStr for Seq<A> {
    type Err = ParseSeqErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut v = Vec::new();
        for i in s.chars() {
            match A::from_char(&i) {
                Ok(b) => v.push(b),
                Err(_) => panic!(),
            }
        }
        Ok(Seq::<A>::from_vec(v))
    }
}

pub struct SeqIter<A: Codec> {
    seq: Seq<A>,
    index: usize,
}

pub struct KmerIter<const K: usize> {
    seq: Seq<Dna>,
    index: usize,
}

impl Seq<Dna> {
    /// Iterate over the k-mers of a sequence.
    ///
    /// K: The number of (characters)[Codec] in the biological sequence
    pub fn kmers<const K: usize>(self) -> KmerIter<K> {
        KmerIter::<K> {
            seq: self,
            index: 0,
        }
    }
}

impl<A: Codec> IntoIterator for Seq<A> {
    type Item = A;
    type IntoIter = SeqIter<A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter::<A> {
            seq: self,
            index: 0,
        }
    }
}

impl<A: Codec> Iterator for SeqIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH;
        let i = self.index;
        if self.index >= (self.seq.bv.len()) {
            return None;
        }
        self.index += w;
        Some(A::from_bits(&self.seq.bv[i..i + w].load()))
    }
}

impl<const K: usize> Iterator for KmerIter<K> {
    type Item = Kmer<K>;
    fn next(&mut self) -> Option<Kmer<K>> {
        let k = K * 2;
        let i = self.index * 2;
        if self.index >= self.seq.len() - (K - 1) {
            return None;
        }
        self.index += 1;
        Some(Kmer::<K>::new(&self.seq.bv[i..k + i]))
    }
}

impl<A: Codec> Seq<A> {
    pub fn len(&self) -> usize {
        self.bv.len() / A::WIDTH
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/*
impl BitAnd for Seq<Iupac> {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        Self {
            bv: self.bv & rhs.bv,
            _p: PhantomData,
        }
    }
}

impl BitOr for Seq<Iupac> {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self {
            bv: self.bv | rhs.bv,
            _p: PhantomData,
        }
    }
}
*/

//impl<Idx> std::ops::Index<Idx> for Seq<Dna>
//where
//    Idx: std::slice::SliceIndex<BitSlice>, //+ std::ops::Mul<Output = usize>,
//{
//    type Output = Idx::Output;
//
//    fn index(&self, i: Idx) -> &Self::Output {
//        let w = Dna::WIDTH as Idx;
//        SeqSlice {
//            bv: &self.bv[i * w..(i * w) + w],
//            _p: PhantomData,
//        }
//    }
//}

#[macro_export]
macro_rules! dna {
    ($seq:expr) => {
        match Seq::<Dna>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}

#[macro_export]
macro_rules! iupac {
    ($seq:expr) => {
        match Seq::<Iupac>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}
