// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

/*!
Sequences of bio alphabet characters. Slicable, Boxable, Iterable.
!*/

pub use crate::alphabet::amino::Amino;
pub use crate::alphabet::dna::Dna;
pub use crate::alphabet::iupac::Iupac;
pub use crate::alphabet::Alphabet;
pub use crate::kmer::Kmer;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;
use std::ops::{BitAnd, BitOr};
pub use std::str::FromStr;

/// A sequence of bit packed characters
pub struct Seq<A: Alphabet> {
    pub bv: BitVec,
    _p: PhantomData<A>,
}

pub struct SeqSlice<'a, A: Alphabet> {
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

impl<A: Alphabet> Seq<A> {
    /// Pack binary representations into a bitvector
    pub fn from_vec(vec: Vec<A>) -> Self {
        let mut bv: BitVec = BitVec::new();
        for b in vec.iter() {
            bv.extend(b.to_bits());
        }
        Seq {
            bv,
            _p: PhantomData,
        }
    }

    pub fn to_usize(&self) -> usize {
        assert!(self.bv.as_raw_slice().len() == 1);
        self.bv.as_raw_slice()[0]
    }
}

impl<A: Alphabet> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let w = A::WIDTH;
        //        for c in self.bv.chunks(2) {
        for i in 0..(self.bv.len() / A::WIDTH) {
            s.push_str(&A::from_bits(&self.bv[(i * w)..((i * w) + w)]).to_string());

            //            v.push(c.as_raw_slice()[0]);
        }
        write!(
            f,
            "{}",
            s,
            //&self.bv[(4 * w)..((4 * w) + w)][2],
        )
    }
}

impl<A: Alphabet> FromStr for Seq<A> {
    type Err = ParseSeqErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut v = Vec::new();
        for i in s.chars() {
            match A::from_str(&i.to_string()) {
                Ok(b) => v.push(b),
                Err(_) => panic!(),
            }
        }
        Ok(Seq::<A>::from_vec(v))
    }
}

pub struct SeqIter<A: Alphabet> {
    seq: Seq<A>,
    index: usize,
}

pub struct KmerIter<const K: u8> {
    seq: Seq<Dna>,
    index: usize,
}

impl Seq<Dna> {
    /// Iterate over the k-mers of a sequence.
    ///
    /// K: The number of (characters)[Alphabet] in the biological sequence
    pub fn kmers<const K: u8>(self) -> KmerIter<K> {
        KmerIter::<K> {
            seq: self,
            index: 0,
        }
    }
}

impl<A: Alphabet> IntoIterator for Seq<A> {
    type Item = A;
    type IntoIter = SeqIter<A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter::<A> {
            seq: self,
            index: 0,
        }
    }
}

impl<A: Alphabet> Iterator for SeqIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH;
        let i = self.index;
        if self.index >= (self.seq.bv.len()) {
            return None;
        }
        self.index += w;
        Some(A::from_bits(&self.seq.bv[i..i + w]))
    }
}

impl<const K: u8> Iterator for KmerIter<K> {
    type Item = Kmer<K>;
    fn next(&mut self) -> Option<Kmer<K>> {
        let k = K as usize * 2;
        let i = self.index * 2;
        if self.index >= self.seq.len() - (K as usize - 1) {
            return None;
        }
        self.index += 1;
        Some(Kmer::<K>::new(&self.seq.bv[i..k + i]))
    }
}

impl<A: Alphabet> Seq<A> {
    pub fn len(&self) -> usize {
        self.bv.len() / A::WIDTH
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

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
macro_rules! amino {
    ($seq:expr) => {
        match Seq::<Amino>::from_str($seq) {
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
