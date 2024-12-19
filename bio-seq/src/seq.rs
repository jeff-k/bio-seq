// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Arbitrary length sequences of bit-packed genomic data
//!
//! `Seq` and `&SeqSlice` are analogous to `String` and `&str`. A `Seq` owns its data and a `SeqSlice` is a read-only window into a `Seq`.
pub mod index;
pub mod iterators;

mod array;
mod slice;

pub use array::SeqArray;
pub use slice::SeqSlice;

use crate::codec::{text, Codec};
use crate::error::ParseBioError;

use crate::{Bs, Bv, Order};

use bitvec::field::BitField;
use bitvec::view::BitView;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use core::borrow::Borrow;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::ops::{Bound, Deref, RangeBounds};
use core::str::FromStr;
use core::{fmt, ptr, str};

/// A arbitrary length sequence of bit-packed symbols
///
/// Stored on the heap
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[repr(transparent)]
pub struct Seq<A: Codec> {
    pub(crate) _p: PhantomData<A>,
    pub(crate) bv: Bv,
}

impl<A: Codec> From<Seq<A>> for usize {
    fn from(slice: Seq<A>) -> usize {
        debug_assert!(slice.bv.len() <= usize::BITS as usize);
        slice.bv.load_le::<usize>() //.wrapping_shr(shift)
    }
}

impl<A: Codec> Hash for Seq<A> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_ref().hash(state);
    }
}

// Note that we could set a default output type with this feature:
// #![feature(associated_type_defaults)]
/// A reversible sequence of things that can be complemented can be reverse complemented
pub trait ReverseComplement {
    type Output;

    /// Reverse complement of a sequence
    fn revcomp(&self) -> Self::Output;
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
            bv: Bv::new(),
        }
    }

    fn bit_range<R: RangeBounds<usize>>(&self, range: R) -> (usize, usize) {
        let s = match range.start_bound() {
            Bound::Included(&n) => n,
            Bound::Excluded(&n) => n + 1,
            Bound::Unbounded => 0,
        };

        let e = match range.end_bound() {
            Bound::Included(&n) => n + 1,
            Bound::Excluded(&n) => n,
            Bound::Unbounded => self.len(),
        };

        debug_assert!(s <= e, "Start of range must be less than or equal to end");
        debug_assert!(e <= self.len(), "Range out of bounds");

        (s * A::BITS as usize, e * A::BITS as usize)
    }

    /// Trim leading and trailing characters that don't match bases/symbols
    /// ```
    /// # use bio_seq::prelude::*;
    /// let seq = b"NNNNAGAATGATGGGGGGGGGGGCGNNNNNNNNNNN";
    /// let trimmed: Seq<Dna> = Seq::trim_u8(seq).unwrap();
    /// assert_eq!(trimmed, dna!("AGAATGATGGGGGGGGGGGCG"));
    /// ```
    ///
    /// # Errors
    ///
    /// Will return an `UnrecognisedBase` error for non-leading or non-trailing end bases, just as
    /// `TryFrom<&[u8]>` would.
    pub fn trim_u8(v: &[u8]) -> Result<Self, ParseBioError> {
        let start = v
            .iter()
            .position(|&byte| A::try_from_ascii(byte).is_some())
            .unwrap_or(v.len());

        let end = v[start..]
            .iter()
            .rposition(|&byte| A::try_from_ascii(byte).is_some())
            .map_or(start, |pos| start + pos + 1);

        v[start..end]
            .iter()
            .map(|&byte| A::try_from_ascii(byte).ok_or(ParseBioError::UnrecognisedBase(byte)))
            .collect()
    }

    pub fn with_capacity(len: usize) -> Self {
        Seq {
            _p: PhantomData,
            bv: Bv::with_capacity(len * A::BITS as usize),
        }
    }

    pub fn bit_and(self, rhs: Seq<A>) -> Seq<A> {
        Seq::<A> {
            _p: PhantomData,
            bv: Bv::from_bitslice(&(self.bv & rhs.bv)),
        }
    }

    pub fn bit_or(self, rhs: Seq<A>) -> Seq<A> {
        Seq::<A> {
            _p: PhantomData,
            bv: Bv::from_bitslice(&(self.bv | rhs.bv)),
        }
    }

    pub fn push(&mut self, item: A) {
        let byte: u8 = item.to_bits();
        self.bv
            .extend_from_bitslice(&byte.view_bits::<Order>()[..A::BITS as usize]);
    }

    pub fn clear(&mut self) {
        self.bv.clear();
    }

    /// Truncate a sequence
    /// ```
    /// # use bio_seq::prelude::*;
    /// let mut seq: Seq<Dna> = dna!("CCCCC").into();
    /// seq.truncate(2);
    /// assert_eq!(&seq, dna!("CC"));
    /// ```
    pub fn truncate(&mut self, len: usize) {
        self.bv.truncate(len * A::BITS as usize);
    }

    /// Prepend a slice
    /// ```
    /// # use bio_seq::prelude::*;
    /// let mut seq: Seq<Dna> = dna!("CCCCC").into();
    /// seq.prepend(dna!("TTTT"));
    /// assert_eq!(&seq, dna!("TTTTCCCCC"));
    /// ```
    pub fn prepend(&mut self, other: &SeqSlice<A>) {
        let mut bv = Bv::with_capacity(self.bv.len() + other.bs.len());
        bv.extend_from_bitslice(&other.bs);
        bv.extend_from_bitslice(&self.bv);
        self.bv = bv;
    }

    /// Append a slice
    /// ```
    /// # use bio_seq::prelude::*;
    /// let mut seq: Seq<Dna> = dna!("CCCCC").into();
    /// seq.append(dna!("TTTT"));
    /// assert_eq!(&seq, dna!("CCCCCTTTT"));
    /// ```
    pub fn append(&mut self, other: &SeqSlice<A>) {
        self.bv.extend_from_bitslice(&other.bs);
    }

    /// Remove a range and replace it with a slice
    /// ```
    /// # use bio_seq::prelude::*;
    /// let mut seq: Seq<Dna> = dna!("AAAACCAAAA").into();
    /// seq.splice(4..6, dna!("TTTT"));
    /// assert_eq!(&seq, dna!("AAAATTTTAAAA"));
    /// ```
    pub fn splice<R: RangeBounds<usize>>(&mut self, range: R, other: &SeqSlice<A>) {
        let (s, e) = self.bit_range(range);
        self.bv.splice(s..e, other.bs.iter().by_vals());
    }

    /// Insert a slice into a sequence
    /// ```
    /// # use bio_seq::prelude::*;
    /// let mut seq: Seq<Dna> = dna!("AAAAA").into();
    /// seq.insert(3, dna!("TTTT"));
    /// assert_eq!(&seq, dna!("AAATTTTAA"));
    /// ```
    ///
    /// # Panics
    /// Panics if index out of bounds
    pub fn insert(&mut self, index: usize, other: &SeqSlice<A>) {
        assert!(index <= self.len(), "Index out of bounds");

        let i = index * A::BITS as usize;
        let mut bv = Bv::with_capacity(self.bv.len() + other.bs.len());

        bv.extend_from_bitslice(&self.bs[..i]);
        bv.extend_from_bitslice(&other.bs);
        bv.extend_from_bitslice(&self.bs[i..]);

        self.bv = bv;
    }

    /// Remove a region of a sequence
    /// ```
    /// # use bio_seq::prelude::*;
    /// let mut seq: Seq<Dna> = dna!("ACGTACGT").into();
    /// seq.remove(2..5);
    /// assert_eq!(&seq, dna!("ACCGT"));
    /// ```
    pub fn remove<R: RangeBounds<usize>>(&mut self, range: R) {
        let (s, e) = self.bit_range(range);
        self.bv.drain(s..e);
    }

    pub fn extend<I: IntoIterator<Item = A>>(&mut self, iter: I) {
        iter.into_iter().for_each(|base| self.push(base));
    }

    /// **Experimental** Decode sequence from raw `usize` array. This requires the exact length of the target sequence to be known.
    /// ```
    /// # use bio_seq::prelude::*;
    /// let ints: [usize; 2] = [0b0101010101010101010101010101010111111111111111111111111111111111, 0b11100100];
    /// let seq: Option<Seq<Dna>> = Seq::from_raw(36, &ints);
    /// assert_eq!(seq.unwrap(), dna!("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCACGT"));
    /// ```
    pub fn from_raw(len: usize, bits: &[usize]) -> Option<Self> {
        let mut bv: Bv = Bv::from_slice(bits);
        //debug_assert!(len <= bv.len(), "desired length is greater than provided bits string");
        if len > bv.len() {
            None
        } else {
            bv.truncate(len * A::BITS as usize);
            Some(Seq {
                _p: PhantomData,
                bv,
            })
        }
    }

    /// **Experimental** Access raw sequence data as `&[usize]`
    /// ```
    /// # use bio_seq::prelude::*;
    /// let seq: Seq<Dna> = dna!("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCACGT").into();
    /// let ints: Vec<usize> = seq.into_raw().iter().copied().collect();
    /// assert_eq!(ints[0], 0b0101010101010101010101010101010111111111111111111111111111111111);
    /// assert_eq!(ints[1], 0b11100100); // ACGT
    /// ```
    pub fn into_raw(&self) -> &[usize] {
        self.bv.as_raw_slice()
    }
}

impl<A: Codec> PartialEq<SeqSlice<A>> for Seq<A> {
    fn eq(&self, other: &SeqSlice<A>) -> bool {
        self.as_ref() == other
    }
}

impl<A: Codec> PartialEq<&SeqSlice<A>> for Seq<A> {
    fn eq(&self, other: &&SeqSlice<A>) -> bool {
        self.as_ref() == *other
    }
}

impl<A: Codec> PartialEq<Seq<A>> for &Seq<A> {
    fn eq(&self, other: &Seq<A>) -> bool {
        **self == *other
    }
}

impl<A: Codec> PartialEq<&Seq<A>> for Seq<A> {
    fn eq(&self, other: &&Seq<A>) -> bool {
        *self == **other
    }
}

/// Borrow a `Seq<A>` as a `SeqSlice<A>`.
///
/// The `Borrow` trait to is used to obtain a reference to a `SeqSlice` from a `Seq`, allowing it to be used wherever a `SeqSlice` is expected.
/// ```
/// use std::collections::HashMap;
/// use bio_seq::prelude::*;
///
/// let reference = dna!("ACGTTCGCATGCTACGACGATC");
///
/// let mut table: HashMap<Seq<Dna>, usize> = HashMap::new();
///
/// // Associate some kind of count with sequences as keys:
/// table.insert(dna!("ACGTT").into(), 1);
/// table.insert(dna!("ACACCCCC").into(), 0);
///
/// // The query is a short window in the reference `Seq`
/// let query: &SeqSlice<Dna> = &reference[..5];
///
/// if let Some(value) = table.get(query) {
///        // `SeqSlice` implements `Display`
///        println!("{query}: {value}");
/// }
/// ```
impl<A: Codec> Borrow<SeqSlice<A>> for Seq<A> {
    fn borrow(&self) -> &SeqSlice<A> {
        self.as_ref()
    }
}

impl<A: Codec> Borrow<SeqSlice<A>> for &Seq<A> {
    fn borrow(&self) -> &SeqSlice<A> {
        self.as_ref()
    }
}

/// Automatic dereferencing of `Seq<A>` to `SeqSlice<A>`.
///
/// ```
/// # use bio_seq::prelude::*;
/// fn count_bases(s: &SeqSlice<Dna>) -> usize {
///    s.len()
/// }
///
/// let seq: Seq<Dna> = dna!("CATCGATCGATC").into();
/// let count = count_bases(&seq);
/// assert_eq!(count, 12);
/// ```
///
impl<A: Codec> Deref for Seq<A> {
    type Target = SeqSlice<A>;

    fn deref(&self) -> &Self::Target {
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bv);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

/// A Seq can be borrowed as a `SeqSlice` through generic constraints.
///
/// ```
/// # use bio_seq::prelude::*;
/// fn count_bases<S: AsRef<SeqSlice<Dna>>>(s: S) -> usize {
///    s.as_ref().len()
/// }
///
/// let seq: Seq<Dna> = dna!("CATCGATCGATC").into();
/// let count = count_bases(seq); // the AsRef implementation allows us to directly pass a Seq
/// assert_eq!(count, 12);
/// ```
///
impl<A: Codec> AsRef<SeqSlice<A>> for Seq<A> {
    fn as_ref(&self) -> &SeqSlice<A> {
        self
    }
}

/// Creates a deep copy of the sequence.
///
/// ```
/// #[macro_use]
/// # use bio_seq::prelude::*;
/// let mut seq1: Seq<Dna> = dna!("CATCGATCGATC").into();
/// let seq2: Seq<Dna> = seq1.clone();
///
/// seq1.push(Dna::A);
/// assert_ne!(seq1, seq2);
/// ```
///
impl<A: Codec> Clone for Seq<A> {
    fn clone(&self) -> Self {
        Self {
            _p: PhantomData,
            bv: self.bv.clone(),
        }
    }
}

impl<A: Codec> FromIterator<A> for Seq<A> {
    fn from_iter<I: IntoIterator<Item = A>>(iter: I) -> Self {
        let i = iter.into_iter();
        let mut seq = Seq::with_capacity(i.size_hint().0);
        seq.extend(i);
        seq
    }
}

impl<A: Codec> From<&Vec<A>> for Seq<A> {
    fn from(vec: &Vec<A>) -> Self {
        // for a general conversion: vec.iter().copied().map(Into::into).collect()

        vec.iter().copied().collect()
    }
}

impl<A: Codec, B: Codec> From<&SeqSlice<A>> for Seq<B>
where
    A: Into<B>,
{
    fn from(slice: &SeqSlice<A>) -> Self {
        slice.iter().map(Into::into).collect()
    }
}

impl<A: Codec, B: Codec, const N: usize, const W: usize> From<&SeqArray<A, N, W>> for Seq<B>
where
    A: Into<B>,
{
    fn from(slice: &SeqArray<A, N, W>) -> Self {
        slice.iter().map(Into::into).collect()
    }
}

impl<A: Codec, B: Codec, const N: usize, const W: usize> From<SeqArray<A, N, W>> for Seq<B>
where
    A: Into<B>,
{
    fn from(slice: SeqArray<A, N, W>) -> Self {
        slice.iter().map(Into::into).collect()
    }
}

impl<A: Codec> TryFrom<&str> for Seq<A> {
    type Error = ParseBioError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        Seq::<A>::try_from(s.as_bytes())
    }
}

impl<A: Codec> TryFrom<String> for Seq<A> {
    type Error = ParseBioError;

    fn try_from(s: String) -> Result<Self, Self::Error> {
        Seq::<A>::try_from(s.as_str())
    }
}

impl<A: Codec> TryFrom<&String> for Seq<A> {
    type Error = ParseBioError;

    fn try_from(s: &String) -> Result<Self, Self::Error> {
        Seq::<A>::try_from(s.as_str())
    }
}

impl<A: Codec> FromStr for Seq<A> {
    type Err = ParseBioError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Seq::<A>::try_from(s)
    }
}

impl<A: Codec> TryFrom<&[u8]> for Seq<A> {
    type Error = ParseBioError;

    fn try_from(v: &[u8]) -> Result<Self, Self::Error> {
        Self::try_from(v.to_vec())
    }
}

impl<A: Codec> TryFrom<Vec<u8>> for Seq<A> {
    type Error = ParseBioError;

    fn try_from(v: Vec<u8>) -> Result<Self, Self::Error> {
        // potentional optimisation: with an extra allocation we could
        // .collect::<Result<Vec<A>, _>>()
        // .map(|v| { let mut seq = Self::with_capacity etc.
        v.into_iter()
            .map(|byte| A::try_from_ascii(byte).ok_or(ParseBioError::UnrecognisedBase(byte)))
            .collect()
    }
}

impl<A: Codec> From<Seq<A>> for String {
    fn from(seq: Seq<A>) -> Self {
        String::from(seq.as_ref())
    }
}

impl<A: Codec> From<&Seq<A>> for String {
    fn from(seq: &Seq<A>) -> Self {
        String::from(seq.as_ref())
    }
}

impl<A: Codec> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.as_ref(), f)
    }
}

impl<A: Codec> Extend<A> for Seq<A> {
    fn extend<T: IntoIterator<Item = A>>(&mut self, iter: T) {
        self.extend(iter);
    }
}

impl From<Vec<usize>> for Seq<text::Dna> {
    fn from(vec: Vec<usize>) -> Self {
        Seq {
            _p: PhantomData,
            bv: Bv::from_vec(vec),
        }
    }
}

/// **Unstable** construct a `Seq` from a bitslice. This may change in the future.
impl<A: Codec> From<&Bs> for Seq<A> {
    fn from(bs: &Bs) -> Self {
        Seq {
            _p: PhantomData,
            bv: bs.into(),
        }
    }
}

/// **Unstable** construct a `Seq` from a bitvec. This may change in the future.
impl<A: Codec> From<Bv> for Seq<A> {
    fn from(bv: Bv) -> Self {
        Seq {
            _p: PhantomData,
            bv,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::text;
    use crate::prelude::*;
    use crate::{Bv, Order};
    use bitvec::prelude::*;
    use core::borrow::Borrow;
    use core::hash::{Hash, Hasher};
    use core::marker::PhantomData;
    use std::collections::hash_map::DefaultHasher;

    #[test]
    fn test_revcomp() {
        let s1: Seq<Dna> = dna!("ATGTGTGCGACTGA").into();
        let s2: Seq<Dna> = dna!("TCAGTCGCACACAT").into();
        let s3: &SeqSlice<Dna> = &s1;
        assert_eq!(s3.revcomp(), s2.revcomp().revcomp());
        assert_eq!(s3.revcomp(), s2);
        assert_ne!(s3.revcomp(), s2.revcomp());
        assert_eq!(s3, s2.revcomp());
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
    fn slice_index_owned() {
        let seq = dna!("GCTCGATCACT");

        assert_eq!(&seq[..], dna!("GCTCGATCACT"));
        assert_eq!(&seq[..=4], dna!("GCTCG"));
        assert_eq!(&seq[1..=4], dna!("CTCG"));
        assert_eq!(&seq[1..4], dna!("CTC"));
        assert_eq!(&seq[..4], dna!("GCTC"));
        assert_eq!(&seq[1..], dna!("CTCGATCACT"));
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
        let s1: &'static SeqSlice<Dna> = dna!("ACGACTGATCGA");
        let s2: &'static SeqSlice<Dna> = dna!("TCGAACGACTGA");

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
        let result_str: String = slice.into();
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
        assert_eq!(String::from(seq), "ACGT")
    }

    #[test]
    fn test_extend_amino() {
        let mut seq = Seq::<Amino>::new();
        seq.push(Amino::S);
        seq.push(Amino::L);

        seq.extend(vec![Amino::Y, Amino::M].into_iter());

        assert_eq!(seq.len(), 4);
        assert_eq!(String::from(seq), "SLYM");
    }
    #[test]
    fn test_extend() {
        let mut seq = Seq::<Dna>::new();
        seq.push(Dna::A);
        seq.push(Dna::C);

        seq.extend(vec![Dna::G, Dna::T].into_iter());

        assert_eq!(seq.len(), 4);
        assert_eq!(String::from(seq), "ACGT");
    }

    #[test]
    fn test_eqs() {
        let seq: Seq<Dna> = "ACTAGCATCGA".try_into().unwrap();
        let seq2: Seq<Dna> = "ACTAGCATCGA".try_into().unwrap();
        let slice: &SeqSlice<Dna> = &seq;
        let slice2: &SeqSlice<Dna> = &seq2[..];
        assert_eq!(seq, slice);
        assert_eq!(seq2, slice);
        assert_eq!(seq2, slice2);
        assert_eq!(slice, slice2);
        assert_eq!(seq, seq2);
    }

    #[test]
    fn test_str_eqs() {
        let string: String = "ACTAGCATCGA".into();
        let slice: &str = "GCTGCATCGATC";

        let seq1: Seq<Dna> = Seq::<Dna>::try_from(slice).unwrap();
        let seq2: Seq<Dna> = Seq::<Dna>::try_from(string.clone()).unwrap();

        assert_eq!(seq2.to_string(), string);
        assert_eq!(seq1.to_string(), slice);

        assert_ne!(seq2.to_string(), slice);
        assert_ne!(seq1.to_string(), string);
    }

    #[test]
    fn test_from_iter() {
        let iter = vec![Dna::A, Dna::C, Dna::G, Dna::T].into_iter();
        let seq: Seq<Dna> = Seq::from_iter(iter);

        assert_eq!(seq.len(), 4);
        assert_eq!(String::from(seq), "ACGT");
    }

    #[test]
    fn test_bit_order() {
        let raw: usize = 0b10_11_01_11_10_01_00_01;
        let mut bv: Bv = Default::default();
        bv.extend(&raw.view_bits::<Order>()[..(Dna::BITS as usize * 8)]);
        let s = Seq::<Dna> {
            bv,
            _p: PhantomData,
        };
        assert_eq!(dna!("CACGTCTG").to_string(), "CACGTCTG");
        assert_eq!(String::from(s), "CACGTCTG");
        //        assert_eq!(raw, Kmer::<Dna, 8>::from(&s[..8]).bs);
    }

    #[test]
    fn test_borrow() {
        let seq: Seq<Dna> = dna!("ACGACCCCCATAGATGGGCTG").into();
        let slice: &SeqSlice<Dna> = seq.borrow();
        assert_eq!(slice, &seq[..]);
        assert_ne!(slice, &seq[1..]);
    }

    #[test]
    fn test_deref() {
        let seq: Seq<Dna> = dna!("AGAATGATCG").into();
        let slice: &SeqSlice<Dna> = &*seq;

        assert_eq!(slice, &seq[..]);
        assert_ne!(slice, &seq[1..]);
    }

    #[test]
    fn test_asref() {
        let seq: Seq<Dna> = dna!("AGAATGATCAAAATATATATAAAG").into();
        let slice: &SeqSlice<Dna> = seq.as_ref();
        assert_ne!(slice, &seq[2..5]);
        assert_eq!(slice, &seq[..]);
    }

    #[test]
    fn test_to_owned() {
        let seq: Seq<Dna> = dna!("AGAATGAATCG").into();
        let slice: &SeqSlice<Dna> = &seq;
        let owned: Seq<Dna> = slice[2..5].to_owned();
        assert_eq!(&owned, &seq[2..5]);
        assert_eq!(owned, seq[2..5].to_owned());
        assert_ne!(&owned, &seq[..]);
    }

    #[test]
    fn test_clone() {
        let seq: Seq<Dna> = dna!("AGAATGATGGGGGGGGGGGCG").into();
        let cloned = seq.clone();
        assert_eq!(seq, cloned);
    }
    #[test]
    fn test_trim() {
        let seq = b"AGAATGATGGGGGGGGGGGCG";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("AGAATGATGGGGGGGGGGGCG"));

        let seq = b"NNNNAGAATGATGGGGGGGGGGGCGNNNNNNNNNNN";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("AGAATGATGGGGGGGGGGGCG"));

        let seq = b"NNNNAGAATGATGGGGNGGGGGGGCGNNNNNNNNNNN";
        let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
        assert_eq!(s, Err(ParseBioError::UnrecognisedBase(b'N')));

        let seq = b"AGAATGATGGGGGGGGGGGCG";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("AGAATGATGGGGGGGGGGGCG"));

        let seq = b"NNNNAGAATGATGGGGGGGGGGGCGNNNNNNNNNNN";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("AGAATGATGGGGGGGGGGGCG"));

        let seq = b"NNNNAGAATGATGGGGNGGGGGGGCGNNNNNNNNNNN";
        let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
        assert_eq!(s, Err(ParseBioError::UnrecognisedBase(b'N')));

        let seq = b"";
        let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
        assert!(s.is_ok());
        assert_eq!(s.unwrap(), dna!(""));

        let seq = b"XXXX";
        let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
        assert!(s.is_ok());
        assert_eq!(s.unwrap(), dna!(""));

        let seq = b"XXACGT";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("ACGT"));

        let seq = b"ACGTXX";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("ACGT"));

        let seq = b"ACGTACGTACGTACGTACGTACGT";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("ACGTACGTACGTACGTACGTACGT"));

        let seq = b"XXACGTXXACGTXX";
        let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
        assert_eq!(s, Err(ParseBioError::UnrecognisedBase(b'X')));

        let seq = b"A";
        let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
        assert_eq!(s, dna!("A"));

        let seq = b"X";
        let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
        assert!(s.is_ok());
        assert_eq!(s.unwrap(), dna!(""));

        /*
                let seq = b"acgtACGT";
                let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
                assert_eq!(s, dna!("ACGTACGT"));

                let seq = b"AcGtAcGt";
                let s: Seq<Dna> = Seq::trim_u8(seq).unwrap();
                assert_eq!(s, dna!("ACGTACGT"));

                let seq = b"AXCXGXTXaXcXgXt";
                let s: Result<Seq<Dna>, ParseBioError> = Seq::trim_u8(seq);
                assert_eq!(s, Err(ParseBioError::UnrecognisedBase(b'X')));
        */
    }

    #[test]
    fn test_seq_eq_and_hash() {
        let seq1: Seq<Dna> = "ACGT".try_into().unwrap();
        let seq2: Seq<Dna> = "ACGT".try_into().unwrap();

        // Equality checks
        assert_eq!(seq1, seq2);

        // Hash checks
        let mut hasher1 = DefaultHasher::new();
        seq1.hash(&mut hasher1);
        let hash1 = hasher1.finish();

        let mut hasher2 = DefaultHasher::new();
        seq2.hash(&mut hasher2);
        let hash2 = hasher2.finish();

        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_seq_slice_eq() {
        let seq1: Seq<Dna> = "ACGTAAAAAAAAAAAAACGTAAAACCCCGGGGTTTTA".try_into().unwrap();
        let seq2: Seq<Dna> = "ACGTAAAAAAAAAAAAACGTAAAACCCCGGGGTTTTAA".try_into().unwrap();

        let slice1a: &SeqSlice<Dna> = &seq1[..];
        let slice1b: &SeqSlice<Dna> = &seq1[..];
        let slice2a: &SeqSlice<Dna> = &seq2[..seq2.len() - 1];
        let slice2b: &SeqSlice<Dna> = &seq2[..seq2.len() - 1];

        let seq3: Seq<Dna> = slice2b.into();

        assert_eq!(slice1a, slice2b);
        assert_eq!(&slice2a, &slice2b);
        assert_eq!(seq1, slice2a);
        assert_eq!(slice1b, seq3);
        assert_eq!(slice1a, slice1b);
        assert_eq!(seq1, seq3);

        assert_ne!(seq1, seq2);
        assert_ne!(seq2, seq3);
        assert_ne!(seq2, slice2a);
        assert_ne!(seq2, slice1b);

        assert_eq!(&slice1a, &slice1b);

        assert_eq!(seq1, &seq1);

        assert_eq!(seq1, &seq3);
        assert_eq!(&seq1, seq3);
        assert_eq!(&seq1, &seq3);

        assert_eq!(slice1a, seq3);
        //assert_eq!(&slice1a, seq3);

        assert_eq!(seq1, slice2a);
        assert_eq!(&seq1, slice2a);
        //assert_eq!(seq1, &slice2a);

        assert_eq!(&seq1, slice2a);
    }
    #[test]
    fn test_seq_slice_hash() {
        let seq1: Seq<Dna> = "ACGTAAAAAAAAAAAAACGTAAAACCCCGGGGAAAAA".try_into().unwrap();
        let seq2: Seq<Dna> = "ACGTAAAAAAAAAAAAACGTAAAACCCCGGGGAAAAA".try_into().unwrap();

        let seq3: Seq<Dna> = "ACGTAAAAAAAAAAAAACGTAAAACCCCGGGG".try_into().unwrap();

        let slice1 = &seq1[..];
        let slice2 = &seq2[..];

        let slice3 = &seq3[..];

        let slice1_32 = &seq1[..32];

        let mut hasher1 = DefaultHasher::new();
        seq1.hash(&mut hasher1);
        let full1 = hasher1.finish();

        let mut hasher1a = DefaultHasher::new();
        seq1.hash(&mut hasher1a);
        let full1_alt = hasher1a.finish();

        let mut hasher2 = DefaultHasher::new();
        seq2.hash(&mut hasher2);
        let full2 = hasher2.finish();

        let mut hasher3 = DefaultHasher::new();
        slice1.hash(&mut hasher3);
        let full1_slice = hasher3.finish();

        let mut hasher4 = DefaultHasher::new();
        slice2.hash(&mut hasher4);
        let full2_slice = hasher4.finish();

        let mut hasher5 = DefaultHasher::new();
        slice3.hash(&mut hasher5);
        let short1_slice = hasher5.finish();

        let mut hasher6 = DefaultHasher::new();
        slice1_32.hash(&mut hasher6);
        let seq1_short = hasher6.finish();

        assert_eq!(full1, full1_alt);
        assert_eq!(full1, full2);
        assert_ne!(full2_slice, short1_slice);
        assert_eq!(full2, full1_slice);
        assert_eq!(short1_slice, seq1_short);
        assert_ne!(seq1_short, full1_slice);

        assert_ne!(full1, short1_slice);
    }

    #[test]
    fn test_fromstr() {
        let seq: Result<Seq<Dna>, ParseBioError> =
            "ACGATGAGTAGTCGCCATCGTATCTTTGACTGCCGATGCTA".parse();
        assert!(seq.is_ok());

        let seq: Result<Seq<Dna>, ParseBioError> =
            "ACGATGAGTAGBCGCCATCGTATCTTTGACTGCCGATGCTA".parse();
        assert_eq!(seq, Err(ParseBioError::UnrecognisedBase(b'B')));
    }

    #[test]
    fn test_lens() {
        assert_eq!(iupac!("AWANWATNA---SKAGTCAA").len(), 20)
    }

    #[test]
    fn test_unique_bitarray_ident() {
        let s1 = dna!("ATCGACTACGATCGCTACGATCGATCGATCGATCGAATCTCCCGCGCGATCATCGATCATCGCTACGTACGTCGAAAAATATAATGGG");
        let s2 = dna!("ATCGACTACGATCGCTACGATCGATCGATCGATCGAATCTCCCGCGCGATCATCGATCATCGCTACGTACGTCGAAAAATATAATGGG");

        // Changed first base
        let s3 = dna!("CTCGACTACGATCGCTACGATCGATCGATCGATCGAATCTCCCGCGCGATCATCGATCATCGCTACGTACGTCGAAAAATATAATGGG");

        // Added `00` to end
        let s4 = dna!("ATCGACTACGATCGCTACGATCGATCGATCGATCGAATCTCCCGCGCGATCATCGATCATCGCTACGTACGTCGAAAAATATAATGGGA");

        // Changed last base
        let s5 = dna!("ATCGACTACGATCGCTACGATCGATCGATCGATCGAATCTCCCGCGCGATCATCGATCATCGCTACGTACGTCGAAAAATATAATGGC");

        assert_eq!(s1, s2);
        assert_ne!(s1, s3);

        // these two tests won't compile because lengths are different
        assert_ne!(s1.len(), s4.len());
        assert_ne!(s4.len(), s1.len());

        assert_ne!(s5, s1);
    }

    #[test]
    fn test_static() {
        use crate::seq::SeqArray;

        static B: SeqArray<Dna, 5, 1> = SeqArray {
            _p: PhantomData,
            ba: bitarr![const usize, Lsb0; 1,1,0,1,0,0,1,0,1,1],
        };

        let s: &'static SeqSlice<Dna> = &B;

        let x: &'static str = "TGACT";

        assert_eq!(s.to_string(), x);
    }

    #[test]
    fn test_to_from_raw() {
        let s = "TCAGCTAGCTACGACTGATCGATCGACTGATGCCGCGCGCGGCGCCGCGCGCGCGCGCCGCGCGCCCCGCGCGCGGCGCGCGCCGCGCGCGCGCGCGGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";

        let seq: Seq<Dna> = s.try_into().unwrap();
        let raw = seq.into_raw();
        let new = Seq::<Dna>::from_raw(s.len(), raw);
        assert_eq!(new.unwrap(), seq);

        let bad = Seq::<Dna>::from_raw(s.len() + 1, raw);
        assert_ne!(bad.unwrap(), seq);

        let bad = Seq::<Dna>::from_raw(342, raw);
        assert_eq!(bad, None);
    }
    #[test]
    fn test_seq_and_seq_slice_eq_and_hash() {
        let seq: Seq<text::Dna> = Seq::try_from("ACGT".to_string()).unwrap();
        let slice = &seq[..];

        // Equality checks
        assert_eq!(seq, slice);
        assert_eq!(&seq, slice);
        assert_eq!(slice, &seq);

        // Hash checks
        let mut hasher1 = std::collections::hash_map::DefaultHasher::new();
        seq.hash(&mut hasher1);
        let hash1 = hasher1.finish();

        let mut hasher2 = std::collections::hash_map::DefaultHasher::new();
        slice.hash(&mut hasher2);
        let hash2 = hasher2.finish();

        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_prepend() {
        let mut seq1 =
            Seq::<Dna>::from_str("GCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAAT").unwrap();
        seq1.prepend(dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        assert_eq!(seq1.to_string(), "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAAT");

        let mut seq2 =
            Seq::<Dna>::from_str("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        seq2.prepend(dna!("GCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAAT"));
        assert_eq!(seq2.to_string(), "GCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        assert_ne!(seq1, seq2);
    }

    #[test]
    fn test_append() {
        let mut seq1 =
            Seq::<Dna>::from_str("GCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAAT").unwrap();
        seq1.append(dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        assert_eq!(seq1.to_string(), "GCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        let mut seq2 =
            Seq::<Dna>::from_str("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        seq2.append(dna!("GCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAAT"));
        assert_eq!(seq2.to_string(), "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCTCGATCGATCGATCGACTGACTGACGCGCGCATCCGATAAAAAAAAT");

        assert_ne!(seq1, seq2);
    }

    #[test]
    fn test_remove() {
        let mut seq: Seq<Dna> = dna!("TCAGCATCGATCAATCG").into();
        seq.remove(4..6);
        assert_eq!(&seq, dna!("TCAGTCGATCAATCG"));
        seq.remove(..4);
        assert_eq!(&seq, dna!("TCGATCAATCG"));
        seq.remove(6..);
        assert_eq!(&seq, dna!("TCGATC"));
    }

    #[test]
    fn test_insert() {
        let mut seq: Seq<Dna> = dna!("TCAGCATCGATCAATCG").into();
        let insertion = dna!("CCCCC");

        seq.insert(4, insertion);
        assert_eq!(&seq, dna!("TCAGCCCCCCATCGATCAATCG"));

        seq.insert(seq.len(), dna!("AAAAA"));
        assert_eq!(&seq, dna!("TCAGCCCCCCATCGATCAATCGAAAAA"));
    }

    #[test]
    fn test_truncate() {
        let mut seq: Seq<Dna> = dna!("TCAGCATCGATCAATCG").into();

        seq.truncate(seq.len());
        assert_eq!(&seq, dna!("TCAGCATCGATCAATCG"));

        seq.truncate(10);
        assert_eq!(&seq, dna!("TCAGCATCGA"));

        seq.truncate(0);
        assert_eq!(&seq, dna!(""));
    }

    #[test]
    fn test_splice() {
        let mut seq: Seq<Dna> = dna!("TCAGCATCGATCAATCG").into();
        let insertion = dna!("CCCCC");

        seq.splice(4..6, insertion);
        assert_eq!(&seq, dna!("TCAGCCCCCTCGATCAATCG"));

        seq.splice(1..=1, dna!("AAA"));
        assert_eq!(&seq, dna!("TAAAAGCCCCCTCGATCAATCG"));

        seq.splice(10.., dna!("TTTT"));
        assert_eq!(&seq, dna!("TAAAAGCCCCTTTT"));
    }
}
