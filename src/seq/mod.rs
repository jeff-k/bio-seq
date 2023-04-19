// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub mod iterators;

use crate::codec::{Codec, Complement};
use crate::error::ParseBioError;

use bitvec::prelude::*;

use core::borrow::Borrow;
use core::convert::TryFrom;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::ops::{
    Deref, Index, IndexMut, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
};
use core::str;
use core::str::FromStr;

/// A sequence of bit-packed characters of arbitrary length
///
/// Allocated to the heap
#[derive(Debug, PartialEq, Eq)]
pub struct Seq<A: Codec> {
    _p: PhantomData<A>,
    bv: BitVec<usize, Lsb0>,
}

/// A slice of a Seq
#[derive(Debug, Eq)]
#[repr(transparent)]
pub struct SeqSlice<A: Codec> {
    _p: PhantomData<A>,
    bs: BitSlice<usize, Lsb0>,
}

impl<A: Codec> From<Seq<A>> for usize {
    fn from(slice: Seq<A>) -> usize {
        assert!(slice.bv.len() <= usize::BITS as usize);
        //let shift = usize::BITS - slice.bv.len() as u32;
        slice.bv.load_le::<usize>() //.wrapping_shr(shift)
    }
}

impl<A: Codec> From<&SeqSlice<A>> for usize {
    fn from(slice: &SeqSlice<A>) -> usize {
        assert!(slice.bs.len() <= usize::BITS as usize);
        //let shift = usize::BITS - slice.bs.len() as u32;
        slice.bs.load_le::<usize>() //.wrapping_shr(shift)
    }
}

impl<A: Codec> From<&SeqSlice<A>> for u8 {
    fn from(slice: &SeqSlice<A>) -> u8 {
        assert!(slice.bs.len() <= u8::BITS as usize);
        //let shift = u8::BITS - slice.bs.len() as u32;
        slice.bs.load_le::<u8>() //.wrapping_shr(shift)
    }
}

pub trait ReverseComplement {
    type A: Codec + Complement;

    /// Reverse complement of a sequence
    fn revcomp(&self) -> Seq<Self::A>;
}

impl<A: Codec + Complement> ReverseComplement for Seq<A> {
    type A = A;

    fn revcomp(&self) -> Seq<A> {
        let mut seq = Seq::<A>::with_capacity(self.len());
        seq.extend(self.rev().map(|base| base.comp()));
        seq
    }
}

impl<A: Codec + Complement> ReverseComplement for SeqSlice<A> {
    type A = A;

    fn revcomp(&self) -> Seq<A> {
        let mut seq = Seq::<A>::with_capacity(self.len());
        seq.extend(self.rev().map(|base| base.comp()));
        seq
    }
}

impl<A: Codec> Default for Seq<A> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: Codec> Seq<A> {
    pub fn new() -> Self {
        Seq {
            _p: PhantomData,
            bv: BitVec::new(),
        }
    }

    pub fn with_capacity(len: usize) -> Self {
        Seq {
            _p: PhantomData,
            bv: BitVec::<usize, Lsb0>::with_capacity(len * A::WIDTH as usize),
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
            bv: BitVec::<usize, Lsb0>::from_bitslice(&(self.bv | rhs.bv)),
        }
    }

    pub fn push(&mut self, item: A) {
        let byte: u8 = item.into();
        self.bv
            .extend_from_bitslice(&byte.view_bits::<Lsb0>()[..A::WIDTH as usize]);
    }

    pub fn clear(&mut self) {
        self.bv.clear();
    }

    pub fn pop(&mut self) -> Option<A> {
        unimplemented!()
    }

    pub fn truncate(&mut self, _len: usize) {
        unimplemented!()
    }

    pub fn remove(&mut self, _index: usize) -> A {
        unimplemented!()
    }

    pub fn insert(&mut self, _index: usize, _element: A) {
        unimplemented!()
    }

    pub fn extend<I: IntoIterator<Item = A>>(&mut self, iter: I) {
        for item in iter {
            self.push(item);
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
        let mut bv: BitVec<usize, Lsb0> = BitVec::from_bitslice(&self.bs);
        bv &= BitVec::from_bitslice(&rhs.bs);

        Seq::<A> {
            _p: PhantomData,
            bv,
        }
    }

    pub fn bit_or(&self, rhs: &SeqSlice<A>) -> Seq<A> {
        let mut bv: BitVec<usize, Lsb0> = BitVec::from_bitslice(&self.bs);
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

impl<A: Codec> From<Seq<A>> for String {
    fn from(seq: Seq<A>) -> Self {
        // potential optimisation: pre-allocate the String upfront:
        // let mut s = String::with_capacity(self.bs.len() / A::WIDTH.into());
        seq.into_iter().map(|base| base.to_char()).collect()
    }
}

impl<A: Codec> PartialEq<&str> for Seq<A> {
    fn eq(&self, other: &&str) -> bool {
        let s: String = self.into_iter().map(|base| base.to_char()).collect();
        s == *other
    }
}

impl<A: Codec> PartialEq<Seq<A>> for SeqSlice<A> {
    fn eq(&self, other: &Seq<A>) -> bool {
        self.bs == other.bv[..]
    }
}

impl<A: Codec> PartialEq<&SeqSlice<A>> for Seq<A> {
    fn eq(&self, other: &&SeqSlice<A>) -> bool {
        self.bv[..] == other.bs
    }
}

impl<A: Codec> Hash for SeqSlice<A> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        // theory: this is prevent Hash(AAAA) from equaling Hash(AAAAA)
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
        let bs = &self.bv[s..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeTo<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        let e = range.end * A::WIDTH as usize;
        let bs = &self.bv[..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeToInclusive<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        let e = (range.end + 1) * A::WIDTH as usize;
        let bs = &self.bv[..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeInclusive<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        let s = range.start() * A::WIDTH as usize;
        let e = (range.end() + 1) * A::WIDTH as usize;
        let bs = &self.bv[s..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFrom<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        let s = range.start * A::WIDTH as usize;
        let bs = &self.bv[s..] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFull> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs = &self.bv[0..self.bv.len()] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
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
        let bs = &self.bs[s..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeTo<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        let e = range.end * A::WIDTH as usize;
        let bs = &self.bs[..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeToInclusive<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        let e = (range.end + 1) * A::WIDTH as usize;
        let bs = &self.bs[..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeInclusive<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        let s = range.start() * A::WIDTH as usize;
        let e = (range.end() + 1) * A::WIDTH as usize;
        let bs = &self.bs[s..e] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFrom<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        let s = range.start * A::WIDTH as usize;
        let bs = &self.bs[s..] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFull> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs = &self.bs[0..self.bs.len()] as *const BitSlice<usize, Lsb0> as *const SeqSlice<A>;
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

impl<A: Codec> ToOwned for SeqSlice<A> {
    type Owned = Seq<A>;

    fn to_owned(&self) -> Self::Owned {
        Seq {
            _p: PhantomData,
            bv: self.bs.into(),
        }
    }
}

impl<A: Codec> Clone for Seq<A> {
    fn clone(&self) -> Self {
        Self {
            _p: PhantomData,
            bv: self.bv.clone(),
        }
    }
}

impl<A: Codec> From<&Vec<A>> for Seq<A> {
    fn from(vec: &Vec<A>) -> Self {
        let mut seq = Seq::<A>::with_capacity(vec.len());
        for c in vec.iter() {
            seq.push(*c);
        }
        seq
    }
}

/*
impl<A: Codec> From<&Vec<u8>> for Seq<A> {
    fn from(vec: &Vec<u8>) -> Self {
        let mut seq = Seq::<A>::with_capacity(vec.len());
        vec.into_iter().map(|c| seq.push(A::unsafe_from_bits(*c)));
        seq
    }
}
*/

impl<A: Codec> From<&SeqSlice<A>> for Seq<A> {
    fn from(slice: &SeqSlice<A>) -> Self {
        Seq {
            _p: PhantomData,
            bv: slice.bs.into(),
        }
    }
}

impl<A: Codec> TryFrom<&str> for Seq<A> {
    type Error = ParseBioError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut seq = Seq::<A>::with_capacity(s.len());
        seq.extend(
            s.chars()
                .map(A::from_char)
                .collect::<Result<Vec<A>, _>>()
                .map_err(|_| ParseBioError {})?,
        );
        Ok(seq)
    }
}

impl<A: Codec> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in self {
            write!(f, "{}", base.to_char())?;
        }
        Ok(())
    }
}

impl<A: Codec> fmt::Display for SeqSlice<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in self {
            write!(f, "{}", base.to_char())?;
        }
        Ok(())
    }
}

impl<A: Codec> FromStr for Seq<A> {
    type Err = ParseBioError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut seq = Seq::with_capacity(s.len());
        for i in s.chars() {
            match A::from_char(i) {
                Ok(b) => seq.push(b),
                Err(_) => return Err(ParseBioError {}),
            }
        }
        Ok(seq)
    }
}

impl<A: Codec> IndexMut<usize> for Seq<A> {
    fn index_mut(&mut self, _index: usize) -> &mut SeqSlice<A> {
        unimplemented!()
    }
}

impl<A: Codec> Extend<A> for Seq<A> {
    fn extend<T: IntoIterator<Item = A>>(&mut self, iter: T) {
        self.extend(iter);
    }
}

impl<A: Codec> FromIterator<A> for Seq<A> {
    fn from_iter<T: IntoIterator<Item = A>>(iter: T) -> Self {
        let mut seq = Seq::new();
        seq.extend(iter);
        seq
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use bitvec::prelude::*;
    use core::marker::PhantomData;

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
    fn slice_index_comparisions() {
        let s1 = dna!("ATGTGTGCGACTGATGATCAAACGTAGCTACG");
        let s2 = dna!("ACGTGTGTGCTAGCTAATCGATCAAAAAG");

        assert_eq!(&s1[0], &s2[0]);
        assert_ne!(&s1[1], &s2[1]);
        assert_eq!(&s1[2..4], &s2[2..4]);
        assert_eq!(&s1[15..21], &s2[19..25]);
        assert_ne!(&s1[2..20], &s2[2..20]);
    }

    #[test]
    fn slice_indexing() {
        let seq = dna!("TGCATCGAT");

        assert_ne!(&seq[..], &dna!("AGCATCGAA")[..]);
        assert_ne!(&seq[3..=6], &dna!("ATC")[..]);
        assert_ne!(&seq[..=6], &dna!("TGCATC")[..]);
        assert_ne!(&seq[4..5], &dna!("TC")[..]);
        assert_ne!(&seq[..6], &dna!("TGCAT")[..]);
        assert_ne!(&seq[5..], &dna!("TCGAT")[..]);

        assert_eq!(&seq[..], &dna!("TGCATCGAT")[..]);
        assert_eq!(&seq[3..=6], &dna!("ATCG")[..]);
        assert_eq!(&seq[..=6], &dna!("TGCATCG")[..]);
        assert_eq!(&seq[4..5], &dna!("T")[..]);
        assert_eq!(&seq[..6], &dna!("TGCATC")[..]);
        assert_eq!(&seq[5..], &dna!("CGAT")[..]);
    }

    #[test]
    fn slice_index_ranges() {
        let s1 = dna!("ACGACTGATCGA");
        let s2 = dna!("TCGAACGACTGA");

        assert_eq!(&s1[..8], &s2[4..]);
        assert_eq!(&s1[8..], &s2[..4]);
        assert_ne!(&s1[8..], &s2[8..]);
        assert_ne!(&s1[..], &s2[..4]);

        assert_eq!(&s1[..=7], &s2[4..]);
        assert_eq!(&s1[8..], &s2[..=3]);
        assert_ne!(&s1[8..11], &s2[8..=11]);
        assert_ne!(&s1[..], &s2[..=4]);
    }

    #[test]
    fn slice_nth() {
        let s = dna!("ATGTGTGCGACTGATGATCAAACGTAGCTACG");

        assert_eq!(s.nth(0), Dna::A);
        assert_ne!(s.nth(0), Dna::G);

        assert_eq!(s.nth(1), Dna::T);
        assert_ne!(s.nth(1), Dna::C);

        assert_eq!(s.nth(s.len() - 1), Dna::G);
        assert_ne!(s.nth(s.len() - 1), Dna::C);
    }

    #[test]
    fn slice_rangeto_and_full() {
        let s1 = dna!("ATCGACTAGCATGCTACG");
        let s2 = dna!("ATCGACTAG");

        assert_eq!(&s1[..s2.len()], &s2[..]);
        assert_ne!(&s2[..s2.len()], &s1[..]);
    }

    #[test]
    fn from_slice() {
        let s1 = dna!("ATGTGTGCGACTGATGATCAAACGTAGCTACG");
        let s: &SeqSlice<Dna> = &s1[15..21];
        assert_eq!(format!("{}", s), "GATCAA");
    }

    #[test]
    fn string_to_seq() {
        let seq_str = "ACTGACTG";
        let seq: Result<Seq<Dna>, _> = seq_str.try_into();
        assert!(seq.is_ok());
        assert_eq!(seq.unwrap().to_string(), seq_str);
    }

    #[test]
    fn invalid_string_to_seq() {
        let invalid_seq_str = "ACUGACTG";
        let seq: Result<Seq<Dna>, _> = invalid_seq_str.try_into();
        assert!(seq.is_err());
    }

    #[test]
    fn seq_to_string() {
        let seq_str = "ACTGACTG";
        let seq: Seq<Dna> = seq_str.try_into().unwrap();
        let result_str: String = seq.into();
        assert_eq!(result_str, seq_str);
    }

    #[test]
    fn seqslice_to_string() {
        let seq_str = "ACTGACTG";
        let seq: Seq<Dna> = seq_str.try_into().unwrap();
        let slice = &seq[1..5];
        let result_str: String = slice.to_owned().into();
        assert_eq!(result_str, "CTGA");
    }

    #[test]
    #[should_panic(expected = "range 2..18 out of bounds: 16")]
    fn invalid_seqslice_to_string() {
        let seq_str = "ACTGACTG";
        let seq: Seq<Dna> = seq_str.try_into().unwrap();
        let _ = &seq[1..9];
    }

    #[test]
    fn test_push() {
        let mut seq = Seq::<Dna>::new();
        seq.push(Dna::A);
        seq.push(Dna::C);
        seq.push(Dna::G);
        seq.push(Dna::T);

        assert_eq!(seq.len(), 4);
        assert_eq!(seq, "ACGT")
    }

    #[test]
    fn test_extend_amino() {
        let mut seq = Seq::<Amino>::new();
        seq.push(Amino::S);
        seq.push(Amino::L);

        seq.extend(vec![Amino::Y, Amino::M].into_iter());

        assert_eq!(seq.len(), 4);
        assert_eq!(seq, "SLYM");
    }

    #[test]
    fn test_extend() {
        let mut seq = Seq::<Dna>::new();
        seq.push(Dna::A);
        seq.push(Dna::C);

        seq.extend(vec![Dna::G, Dna::T].into_iter());

        assert_eq!(seq.len(), 4);
        assert_eq!(seq, "ACGT");
    }

    #[test]
    fn test_eqs() {
        let seq: Seq<Dna> = "ACTAGCATCGA".try_into().unwrap();
        let seq2: Seq<Dna> = "ACTAGCATCGA".try_into().unwrap();
        let slice: &SeqSlice<Dna> = &seq;
        let slice2: &SeqSlice<Dna> = &seq2[..];
        assert_eq!(seq, slice);
        assert_eq!(seq2, slice);
        assert_eq!(slice, slice2);
        assert_eq!(seq, seq2);
    }

    #[test]
    fn test_from_iter() {
        let iter = vec![Dna::A, Dna::C, Dna::G, Dna::T].into_iter();
        let seq: Seq<Dna> = Seq::from_iter(iter);

        assert_eq!(seq.len(), 4);
        assert_eq!(seq, "ACGT");
    }

    #[test]
    fn test_bit_order() {
        let raw: usize = 0b10_11_01_11_10_01_00_01;
        let mut bv: BitVec<usize, Lsb0> = Default::default();
        bv.extend(&raw.view_bits::<Lsb0>()[..(Dna::WIDTH as usize * 8)]);
        let s = Seq::<Dna> {
            bv,
            _p: PhantomData,
        };
        assert_eq!(dna!("CACGTCTG"), "CACGTCTG");
        assert_eq!(s, "CACGTCTG");
        assert_eq!(raw, Kmer::<Dna, 8>::from(&s[..8]).bs);
    }
}
