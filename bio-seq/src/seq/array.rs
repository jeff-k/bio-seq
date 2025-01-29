// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;

use crate::seq::slice::SeqSlice;
//use crate::seq::Seq;

use crate::{Ba, Bs};

//use bitvec::field::BitField;

//use core::fmt;
use core::marker::PhantomData;
use core::ops::Deref;
use core::ptr;

//use core::ops::{BitAnd, BitOr};
//use std::hash::{Hash, Hasher};

/// Static bit-packed sequence meant to be accessed as a `&'static SeqSlice`
///
/// ```
/// # use bio_seq::prelude::*;
/// // internally, `dna!` constructs a `SeqArray<Dna>`
/// let seq: &'static SeqSlice<Dna> = dna!("ACGTGT");
/// ```
#[derive(Debug)]
#[repr(transparent)]
pub struct SeqArray<A: Codec, const N: usize, const W: usize> {
    pub _p: PhantomData<A>,
    pub ba: Ba<W>,
}

/*
impl<A: Codec, const K: usize, const W: usize> Hash for SeqArray<A, K, W> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let bs: &SeqSlice<A> = self.as_ref();
        bs.hash(state);
    }
}
*/

impl<A: Codec, const N: usize, const W: usize> Deref for SeqArray<A, N, W> {
    type Target = SeqSlice<A>;

    fn deref(&self) -> &Self::Target {
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.ba[..N * A::BITS as usize]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec, const N: usize, const W: usize> AsRef<SeqSlice<A>> for SeqArray<A, N, W> {
    fn as_ref(&self) -> &SeqSlice<A> {
        self
    }
}
/*
impl<A: Codec, const N: usize> From<&SeqArray<A, N, 1>> for usize {
    fn from(slice: &SeqArray<A, N, 1>) -> usize {
        slice.bs.load_le::<usize>()
    }
}

impl<A: Codec, const N: usize, const W: usize> fmt::Display for SeqArray<A, N, W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.as_ref(), f)
    }
}

impl<A: Codec, const N: usize, const W: usize> BitAnd for &SeqArray<A, N, W> {
    type Output = Seq<A>;

    fn bitand(self, rhs: Self) -> Self::Output {
        let mut bv = self.ba.to_bitvec();
        bv &= &rhs.ba;
        bv.truncate(N * A::BITS as usize);

        Seq::<A> {
            bv,
            _p: PhantomData,
        }
    }
}

impl<A: Codec, const N: usize, const W: usize> BitOr for &SeqArray<A, N, W> {
    type Output = Seq<A>;

    fn bitor(self, rhs: Self) -> Self::Output {
        let mut bv = self.ba.to_bitvec();
        bv |= &rhs.ba;
        bv.truncate(N * A::BITS as usize);

        Seq::<A> {
            bv,
            _p: PhantomData,
        }
    }
}
*/
