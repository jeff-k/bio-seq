// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//use crate::codec::{Codec, Complement};
use crate::codec::{Codec, Complement};
use crate::error::ParseBioError;
//use crate::seq::array::SeqArray;
use crate::seq::ReverseComplement;
use crate::seq::Seq;

use crate::Order;
use bitvec::field::BitField;
use bitvec::prelude::*;

use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::str;

use core::ops::{BitAnd, BitOr};

/// A lightweight, read-only window into part of a sequence
#[derive(Debug, PartialEq, Eq)]
#[repr(transparent)]
pub struct SeqSlice<A: Codec> {
    pub(crate) _p: PhantomData<A>,
    pub(crate) bs: BitSlice<usize, Order>,
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
        assert!(slice.bs.len() <= u8::BITS as usize);
        slice.bs.load_le::<u8>()
    }
}

impl<A: Codec + Complement> ReverseComplement for SeqSlice<A> {
    type Output = Seq<A>;

    /// The inefficient default complementation of reverse
    fn revcomp(&self) -> Seq<A> {
        let mut seq = Seq::<A>::with_capacity(self.len());
        seq.extend(self.rev().map(|base| base.comp()));
        seq
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

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<A: Codec> From<&SeqSlice<A>> for String {
    fn from(seq: &SeqSlice<A>) -> Self {
        seq.into_iter().map(Codec::to_char).collect()
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

impl<A: Codec> Hash for SeqSlice<A> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bs.hash(state);
        // theory: this is prevent Hash(AAAA) from equaling Hash(AAAAA)
        //self.len().hash(state);
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
