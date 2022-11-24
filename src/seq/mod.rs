// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

mod iterators;

use crate::codec::{Codec, Complement, ReverseComplement};

use bitvec::prelude::*;

use core::borrow::Borrow;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::ops::{Deref, Index, Range, RangeFull};
use core::str;
pub use core::str::FromStr;

/// A sequence of bit-packed characters of arbitrary length
///
/// Allocated to the heap
#[derive(Debug, PartialEq, Eq)]
pub struct Seq<A: Codec> {
    _p: PhantomData<A>,
    bv: BitVec,
}

/// A slice of a Seq
#[derive(Debug, Eq)]
#[repr(transparent)]
pub struct SeqSlice<A: Codec> {
    _p: PhantomData<A>,
    bs: BitSlice,
}

impl<A: Codec> From<Seq<A>> for usize {
    fn from(slice: Seq<A>) -> usize {
        assert!(slice.bv.len() <= usize::BITS as usize);
        slice.bv.load::<usize>()
    }
}

impl<A: Codec> From<&SeqSlice<A>> for usize {
    fn from(slice: &SeqSlice<A>) -> usize {
        assert!(slice.bs.len() <= usize::BITS as usize);
        slice.bs.load::<usize>()
    }
}

impl<A: Codec> From<&SeqSlice<A>> for u8 {
    fn from(slice: &SeqSlice<A>) -> u8 {
        slice.bs.load::<u8>()
    }
}

impl<A: Codec + Complement + std::fmt::Debug> ReverseComplement for Seq<A> {
    fn revcomp(self) -> Seq<A> {
        let mut v = vec![];
        for base in self.rev() {
            v.push(base.comp());
        }
        Seq::<A>::from_vec(v)
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
            _p: PhantomData,
            bv,
        }
    }

    pub fn nth(&self, i: usize) -> A {
        A::unsafe_from_bits(self[i].into())
    }

    pub fn len(&self) -> usize {
        self.bv.len() / A::WIDTH as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn bit_and(self, rhs: Seq<A>) -> Seq<A> {
        Seq::<A> {
            _p: PhantomData,
            bv: BitVec::from_bitslice(&(self.bv & rhs.bv)),
        }
    }

    pub fn bit_or(self, rhs: Seq<A>) -> Seq<A> {
        Seq::<A> {
            _p: PhantomData,
            bv: BitVec::from_bitslice(&(self.bv | rhs.bv)),
        }
    }
}

impl<A: Codec> SeqSlice<A> {
    pub fn nth(&self, i: usize) -> A {
        A::unsafe_from_bits(self[i].into())
    }

    pub fn len(&self) -> usize {
        self.bs.len() / A::WIDTH as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /*
    pub fn contains(&self, _other: &SeqSlice<A>) -> bool {
        unimplemented!()
    }
    */

    pub fn bit_and(&self, rhs: &SeqSlice<A>) -> Seq<A> {
        let mut bv: BitVec = BitVec::from_bitslice(&self.bs);
        bv &= BitVec::from_bitslice(&rhs.bs);

        Seq::<A> {
            _p: PhantomData,
            bv,
        }
    }

    pub fn bit_or(&self, rhs: &SeqSlice<A>) -> Seq<A> {
        let mut bv: BitVec = BitVec::from_bitslice(&self.bs);
        bv |= BitVec::from_bitslice(&rhs.bs);

        Seq::<A> {
            _p: PhantomData,
            bv,
        }
    }
}

impl<A: Codec> PartialEq for SeqSlice<A> {
    fn eq(&self, other: &Self) -> bool {
        self.bs == other.bs
    }
}

impl<A: Codec> Hash for SeqSlice<A> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        self.len().hash(state);
    }
}
impl<A: Codec> Borrow<SeqSlice<A>> for Seq<A> {
    fn borrow(&self) -> &SeqSlice<A> {
        &self[..]
    }
}

impl<A: Codec> Index<Range<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::WIDTH as usize;
        let e = range.end * A::WIDTH as usize;
        let bs = &self.bv[s..e] as *const BitSlice as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFull> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs = &self.bv[0..self.bv.len() as usize] as *const BitSlice as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<usize> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        &self[i..i + 1]
    }
}

impl<A: Codec> Index<Range<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::WIDTH as usize;
        let e = range.end * A::WIDTH as usize;
        let bs = &self.bs[s..e] as *const BitSlice as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFull> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs = &self.bs[0..self.bs.len() as usize] as *const BitSlice as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<usize> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        //A::unsafe_from_bits(self.bs[s..e].load());
        //&self.bs[s..e].load::<u8>()
        &self[i..i + 1]
    }
}

impl<A: Codec> Deref for Seq<A> {
    type Target = SeqSlice<A>;
    fn deref(&self) -> &Self::Target {
        &self[..]
    }
}

impl<A: Codec> AsRef<SeqSlice<A>> for Seq<A> {
    fn as_ref(&self) -> &SeqSlice<A> {
        &self[..]
    }
}

impl<A: Codec> AsRef<Seq<A>> for Seq<A> {
    fn as_ref(&self) -> &Self {
        unimplemented!()
    }
}

impl<A: Codec> From<&SeqSlice<A>> for Seq<A> {
    fn from(_slice: &SeqSlice<A>) -> Self {
        unimplemented!()
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

impl<A: Codec> From<Vec<u8>> for Seq<A> {
    fn from(v: Vec<u8>) -> Self {
        Seq::<A>::from_vec(
            v.iter()
                .map(|c| match A::from_char(*c as char) {
                    Ok(base) => base,
                    _ => panic!(),
                })
                .collect(),
        )
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

        assert_eq!(&s1[0], &s2[0]);
        assert_ne!(&s1[1], &s2[1]);
        assert_eq!(&s1[2..4], &s2[2..4]);
        assert_eq!(&s1[15..21], &s2[19..25]);
        assert_ne!(&s1[2..20], &s2[2..20]);
    }

    #[test]
    fn from_slice() {
        let s1 = dna!("ATGTGTGCGACTGATGATCAAACGTAGCTACG");
        let s: &SeqSlice<Dna> = &s1[15..21];
        assert_eq!(format!("{}", s), "GATCAA");
    }
}
