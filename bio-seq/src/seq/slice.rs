// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::error::ParseBioError;
use crate::seq::Seq;
use crate::{
    Complement, ComplementMut, Reverse, ReverseComplement, ReverseComplementMut, ReverseMut,
};

use crate::Bs;
use bitvec::field::BitField;

use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::str;

use core::ops::{BitAnd, BitOr};

/// An unsized, read-only window into part of a sequence
#[derive(Debug, Eq)]
#[repr(transparent)]
pub struct SeqSlice<A: Codec> {
    pub(crate) _p: PhantomData<A>,
    pub(crate) bs: Bs,
}

impl<A: Codec> TryFrom<&SeqSlice<A>> for usize {
    type Error = ParseBioError;

    fn try_from(slice: &SeqSlice<A>) -> Result<usize, Self::Error> {
        if slice.bs.len() <= usize::BITS as usize {
            Ok(slice.bs.load_le::<usize>())
        } else {
            let len: usize = slice.bs.len() / A::BITS as usize;
            let expected: usize = usize::BITS as usize / A::BITS as usize;
            Err(ParseBioError::SequenceTooLong(len, expected))
        }
    }
}

impl<A: Codec> From<&SeqSlice<A>> for u8 {
    fn from(slice: &SeqSlice<A>) -> u8 {
        debug_assert!(slice.bs.len() <= u8::BITS as usize);
        slice.bs.load_le::<u8>()
    }
}

impl<A: Codec> SeqSlice<A> {
    /// unsafely index into the `i`th position of a sequence
    pub fn nth(&self, i: usize) -> A {
        A::unsafe_from_bits(self[i].into())
    }

    pub fn len(&self) -> usize {
        self.bs.len() / A::BITS as usize
    }

    /// Get the `i`th element of a `Seq`. Returns `None` if index out of range.
    pub fn get(&self, i: usize) -> Option<A> {
        if i >= self.bs.len() / A::BITS as usize {
            None
        } else {
            Some(A::unsafe_from_bits(self[i].into()))
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<A: Codec> From<&SeqSlice<A>> for String {
    fn from(seq: &SeqSlice<A>) -> Self {
        seq.into_iter().map(Codec::to_char).collect()
    }
}

impl<A: Codec> PartialEq<SeqSlice<A>> for SeqSlice<A> {
    fn eq(&self, other: &SeqSlice<A>) -> bool {
        self.bs == other.bs
    }
}

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

/// Warning! hashes are not currently stable between platforms/version
impl<A: Codec> Hash for SeqSlice<A> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        // prepend length to make robust against matching prefixes
        self.len().hash(state);
    }
}

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
///
impl<A: Codec> ToOwned for SeqSlice<A> {
    type Owned = Seq<A>;

    fn to_owned(&self) -> Self::Owned {
        Seq {
            _p: PhantomData,
            bv: self.bs.into(),
        }
    }
}

impl<A: Codec> AsRef<SeqSlice<A>> for SeqSlice<A> {
    fn as_ref(&self) -> &SeqSlice<A> {
        self
    }
}

impl<A: Codec> fmt::Display for SeqSlice<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", String::from(self))
    }
}

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
