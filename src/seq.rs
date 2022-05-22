// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

///Biological sequence types
use crate::codec::iupac::Iupac;
use crate::codec::Codec;
use crate::kmer::Kmer;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;
use std::ops::{BitAnd, BitOr};
pub use std::str::FromStr;

/// A sequence of bit packed characters
pub struct Seq<A: Codec> {
    pub bv: BitVec,
    _p: PhantomData<A>,
}

pub struct SeqSlice<A: Codec> {
    pub bs: BitBox,
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
            let byte: u8 = (*b).into();
            bv.extend_from_bitslice(&(byte as u8).view_bits::<Lsb0>()[..A::WIDTH as usize]);
        }
        Seq {
            bv,
            _p: PhantomData,
        }
    }

    pub fn nth(&self, i: usize) -> A {
        let w = A::WIDTH as usize;
        A::unsafe_from_bits(self.bv[i * w..(i * w) + w].load())
    }

    pub fn len(&self) -> usize {
        self.bv.len() / A::WIDTH as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over the k-mers of a sequence.
    ///
    /// K: The number of (characters)[Codec] in the biological sequence
    pub fn kmers<const K: usize>(self) -> KmerIter<A, K> {
        KmerIter::<A, K> {
            bs: BitBox::from_bitslice(&self.bv),
            index: 0,
            len: self.len(),
            _p: PhantomData,
        }
    }

    pub fn rev(self) -> RevIter<A> {
        let index = self.bv.len();
        RevIter::<A> { seq: self, index }
    }

    pub fn raw(&self) -> &[usize] {
        self.bv.as_raw_slice()
    }
}

impl<A: Codec> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for c in self.bv.chunks_exact(A::WIDTH.into()) {
            s.push_str(&A::unsafe_from_bits(c.load()).to_char().to_string());
        }
        write!(f, "{}", s,)
    }
}

impl<A: Codec> FromStr for Seq<A> {
    type Err = ParseSeqErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut v = Vec::new();
        for i in s.chars() {
            match A::from_char(i) {
                Ok(b) => v.push(b),
                Err(_) => return Err(ParseSeqErr),
            }
        }
        Ok(Seq::<A>::from_vec(v))
    }
}
pub struct RevIter<A: Codec> {
    seq: Seq<A>,
    index: usize,
}

impl<A: Codec> Iterator for RevIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH as usize;
        if self.index == 0 {
            return None;
        }

        self.index -= w;
        let i = self.index;
        Some(A::unsafe_from_bits(self.seq.bv[i..i + w].load()))
        //Some(self.seq.bv[i..i+w].as_ref())
    }
}

impl<A: Codec> IntoIterator for Seq<A> {
    type Item = A;
    type IntoIter = SeqIter<A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter::<A> {
            seq: SeqSlice {
                bs: BitBox::from_bitslice(&self.bv),
                _p: self._p,
            },
            index: 0,
        }
    }
}

pub struct SeqIter<A: Codec> {
    seq: SeqSlice<A>,
    index: usize,
}

impl<A: Codec> Iterator for SeqIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH as usize;
        let i = self.index;
        if self.index >= (self.seq.bs.len()) {
            return None;
        }
        self.index += w;
        Some(A::unsafe_from_bits(self.seq.bs[i..i + w].load()))
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
