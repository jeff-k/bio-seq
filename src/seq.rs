// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::{Codec, Complement, ReverseComplement};
use crate::Kmer;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;
use std::ops::{Index, Range};
pub use std::str::FromStr;

/// A sequence of bit-packed characters of arbitrary length
///
/// Allocated to the heap
#[derive(Debug, PartialEq, Eq)]
pub struct Seq<A: Codec> {
    pub bv: BitVec,
    _p: PhantomData<A>,
}

pub struct SeqSlice<A: Codec> {
    pub bs: BitBox,
    _p: PhantomData<A>,
}

impl<A: Codec> From<&BitSlice> for SeqSlice<A> {
    fn from(slice: &BitSlice) -> SeqSlice<A> {
        SeqSlice {
            bs: BitBox::from_bitslice(slice),
            _p: PhantomData,
        }
    }
}

impl<A: Codec + Complement> ReverseComplement for Seq<A> {
    fn revcomp(self) -> Self {
        let mut v = vec![];
        for base in self.rev() {
            v.push(base.comp());
        }
        Self::from_vec(v)
    }
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

    pub fn from_bitslice(slice: &BitSlice) -> Self {
        Seq {
            bv: BitVec::from_bitslice(slice),
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

    /*
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
    */

    /// Iterate over the sequence in reverse order
    pub fn rev(self) -> RevIter<A> {
        let index = self.bv.len();
        RevIter::<A> { seq: self, index }
    }

    /// Iterate over sliding windows of size K
    pub fn kmers<const K: usize>(self) -> KmerIter<A, K> {
        KmerIter::<A, K> {
            bs: BitBox::from_bitslice(&self.bv),
            index: 0,
            len: self.len(),
            _p: PhantomData,
        }
    }

    pub fn raw(&self) -> &[usize] {
        self.bv.as_raw_slice()
    }
}

impl<A: Codec> Index<Range<usize>> for Seq<A> {
    type Output = BitSlice;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::WIDTH as usize;
        let e = range.end * A::WIDTH as usize;
        &self.bv[s..e]
    }
}

impl<A: Codec> Index<usize> for Seq<A> {
    type Output = BitSlice;

    fn index(&self, i: usize) -> &Self::Output {
        let s = i * A::WIDTH as usize;
        let e = s + A::WIDTH as usize;
        &self.bv[s..e]
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

impl<A: Codec> fmt::Display for SeqSlice<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for c in self.bs.chunks_exact(A::WIDTH.into()) {
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
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::*;
    use crate::codec::ReverseComplement;
    use crate::seq::{FromStr, Seq, SeqSlice};

    #[test]
    fn test_revcomp() {
        let s1 = dna!("ATGTGTGCGACTGA");
        let s2 = dna!("TCAGTCGCACACAT");
        assert_eq!(s1, s2.revcomp());
    }

    #[test]
    fn test_revcomp_mismatched_sizes() {
        let s1 = dna!("AAAA");
        let s2 = dna!("TTTTT");
        assert_ne!(s1, s2.revcomp());
    }

    #[test]
    fn test_revcomp_idempotence() {
        let s = dna!("AAACGCTACGTACGCGCCTTCGGGGCATCAGCACCAC");
        let sc = dna!("AAACGCTACGTACGCGCCTTCGGGGCATCAGCACCAC");

        assert_eq!(s.revcomp().revcomp(), sc);
    }

    #[test]
    fn slice_indices() {
        let s1 = dna!("ATGTGTGCGACTGATGATCAAACGTAGCTACG");
        let s2 = dna!("ACGTGTGTGCTAGCTAATCGATCAAAAAG");

        assert_eq!(s1[0], s2[0]);
        assert_ne!(s1[1], s2[1]);
        assert_eq!(s1[2..4], s2[2..4]);
        assert_eq!(s1[15..21], s2[19..25]);
        assert_ne!(s1[2..20], s2[2..20]);
    }

    #[test]
    fn from_slice() {
        let s1 = dna!("ATGTGTGCGACTGATGATCAAACGTAGCTACG");
        let s: SeqSlice<Dna> = s1[15..21].into();
        assert_eq!(format!("{}", s), "GATCAA");
    }
}
