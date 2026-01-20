// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::error::ParseBioError;
//use crate::seq::Seq;
use crate::storage::{BitSliceStorage, SeqSliceStorage};
//use crate::{
//    Complement, ComplementMut, Reverse, ReverseComplement, ReverseComplementMut, ReverseMut,
//};

use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
//use core::ops::{BitAnd, BitOr};
//use bitvec::prelude::*;

//use core::ops::{BitAnd, BitOr};

/// An unsized, read-only window into part of a sequence
#[derive(Debug)]
#[repr(transparent)]
pub struct SeqSlice<A: Codec, S: SeqSliceStorage + ?Sized = BitSliceStorage> {
    pub(crate) _p: PhantomData<A>,
    pub(crate) bs: S,
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> TryFrom<&SeqSlice<A, S>> for usize {
    type Error = ParseBioError;

    fn try_from(_slice: &SeqSlice<A, S>) -> Result<usize, Self::Error> {
        todo!()
        /*
        if slice.bs.len() <= usize::BITS as usize {
            Ok(slice.bs.load_le::<usize>())
        } else {
            let len: usize = slice.bs.len() / A::BITS as usize;
            let expected: usize = usize::BITS as usize / A::BITS as usize;
            Err(ParseBioError::SequenceTooLong(len, expected))
        }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> PartialEq for SeqSlice<A, S> {
    fn eq(&self, other: &Self) -> bool {
        self.bs == other.bs
    }
}

/*
impl<A: Codec> From<&SeqSlice<A>> for u8 {
    fn from(slice: &SeqSlice<A>) -> u8 {
        todo!()
        /*
        debug_assert!(slice.bs.len() <= u8::BITS as usize);
        slice.bs.load_le::<u8>()
        */
    }
}
*/

impl<A: Codec, S: SeqSliceStorage + ?Sized> SeqSlice<A, S> {
    /// unsafely index into the `i`th position of a sequence
    pub fn nth(&self, i: usize) -> A {
        let s: usize = i * A::BITS as usize;
        let e: usize = s + A::BITS as usize;
        A::unsafe_from_bits(self.bs.get(s, e))
    }

    pub fn len(&self) -> usize {
        self.bs.len() / A::BITS as usize
    }

    /// Get the `i`th element of a `Seq`. Returns `None` if index out of range.
    pub fn get(&self, i: usize) -> Option<A> {
        if i >= self.bs.len() / A::BITS as usize {
            None
        } else {
            Some(self.nth(i))
        }
    }

    pub fn is_empty(&self) -> bool {
        self.bs.len() == 0
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> From<&SeqSlice<A, S>> for String {
    fn from(_seq: &SeqSlice<A, S>) -> Self {
        todo!()
        //        seq.into_iter().map(Codec::to_char).collect()
    }
}

/*
impl<A: Codec, S: SeqStorage> PartialEq<SeqSlice<A, S>> for SeqSlice<A, S> {
    fn eq(&self, other: &SeqSlice<A, S>) -> bool {
        self.bs == other.bs
    }
}
*/

/*
impl<A: Codec> PartialEq<SeqSlice<A>> for &SeqSlice<A> {
    fn eq(&self, other: &SeqSlice<A>) -> bool {
        self.bs == other.bs
    }
}

impl<A: Codec> PartialEq<Seq<A>> for SeqSlice<A> {
    fn eq(&self, other: &Seq<A>) -> bool {
        self == other.as_ref()
    }
}

impl<A: Codec> PartialEq<Seq<A>> for &SeqSlice<A> {
    fn eq(&self, other: &Seq<A>) -> bool {
        *self == other.as_ref()
    }
}

impl<A: Codec> PartialEq<&str> for SeqSlice<A> {
    fn eq(&self, other: &&str) -> bool {
        let bs = other.as_bytes();
        if bs.len() != self.len() {
            return false;
        }
        for (a, c) in self.iter().zip(bs) {
            match A::try_from_ascii(*c) {
                Some(b) => {
                    if a != b {
                        return false;
                    }
                }
                None => return false,
            }
        }
        true
    }
}
*/

/// Warning! hashes are not currently stable between platforms/version
impl<A: Codec, S: SeqSliceStorage + ?Sized> Hash for SeqSlice<A, S> {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        todo!()
        /*
        self.bs.hash(state);
        // prepend length to make robust against matching prefixes
        self.len().hash(state);
        */
    }
}

/*
/// Clone a borrowed slice of a sequence into an owned version.
///
/// ```
/// use bio_seq::prelude::*;
///
/// let seq = dna!("CATCGATCGATCG");
/// let slice = &seq[2..7]; // TCGAT
/// let owned = slice.to_owned();
///
/// assert_eq!(&owned, &seq[2..7]);
/// ```

impl<A: Codec, S: SeqStorage> ToOwned for SeqSlice<A, S> {
    type Owned = Seq<A, S>;

    fn to_owned(&self) -> Self::Owned {
        todo!()
        /*
        Seq {
            _p: PhantomData,
            bv: self.bs.into(),
        }
        */
    }
}
*/

impl<A: Codec, S: SeqSliceStorage + ?Sized> AsRef<SeqSlice<A, S>> for SeqSlice<A, S> {
    fn as_ref(&self) -> &SeqSlice<A, S> {
        self
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> fmt::Display for SeqSlice<A, S> {
    fn fmt(&self, _f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!()
        //write!(f, "{}", String::from(self))
    }
}

/*
impl<A: Codec> BitAnd for &SeqSlice<A> {
    type Output = Seq<A>;

    fn bitand(self, rhs: Self) -> Self::Output {
        let mut bv = self.bs.to_bitvec();
        bv &= &rhs.bs;
        Seq::<A> {
            bv,
            _p: PhantomData,
        }
    }
}

impl<A: Codec> BitOr for &SeqSlice<A> {
    type Output = Seq<A>;

    fn bitor(self, rhs: Self) -> Self::Output {
        let mut bv = self.bs.to_bitvec();
        bv |= &rhs.bs;

        Seq::<A> {
            bv,
            _p: PhantomData,
        }
    }
}
*/

/*
impl<A: Codec> ReverseMut for SeqSlice<A> {
    fn rev(&mut self) {
        self.bs.reverse();
        for chunk in self.bs.rchunks_exact_mut(A::BITS as usize) {
            chunk.reverse();
        }
    }
}

impl<A: Codec + ComplementMut> ComplementMut for SeqSlice<A> {
    fn comp(&mut self) {
        unsafe {
            for base in self.bs.chunks_exact_mut(A::BITS as usize).remove_alias() {
                let mut bc = A::unsafe_from_bits(base.load_le::<u8>());
                bc.comp();
                base.store(bc.to_bits() as usize);
            }
        }
    }
}
*/

/*
impl<A: Codec + ComplementMut> ReverseComplementMut for SeqSlice<A> where
    SeqSlice<A>: ComplementMut + ReverseMut
{
}

impl<A: Codec> Reverse for SeqSlice<A> {}

impl<A: Codec + ComplementMut> Complement for SeqSlice<A> {}

impl<A: Codec + ComplementMut> ReverseComplement for SeqSlice<A> where
    SeqSlice<A>: ComplementMut + ReverseMut
{
}
*/
