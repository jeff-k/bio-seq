// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;

use crate::seq::slice::SeqSlice;

use crate::{Bs, Order};
use bitvec::field::BitField;
use bitvec::prelude::*;

use core::fmt;
use core::marker::PhantomData;
use core::ops::Deref;
use core::ptr;

#[derive(Debug)]
#[repr(transparent)]
pub struct SeqArray<A: Codec, const N: usize, const W: usize> {
    pub _p: PhantomData<A>,
    pub ba: BitArray<[usize; W], Order>,
}

impl<A: Codec, const N: usize, const W: usize> Deref for SeqArray<A, N, W> {
    type Target = SeqSlice<A>;

    fn deref(&self) -> &Self::Target {
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.ba[..N * A::BITS as usize]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec, const N: usize> From<&SeqArray<A, N, 1>> for usize {
    fn from(slice: &SeqArray<A, N, 1>) -> usize {
        slice.bs.load_le::<usize>()
    }
}

impl<A: Codec, const N: usize, const W: usize, const M: usize, const V: usize>
    PartialEq<SeqArray<A, N, W>> for SeqArray<A, M, V>
{
    fn eq(&self, other: &SeqArray<A, N, W>) -> bool {
        if N == M {
            self.ba == other.ba
        } else {
            false
        }
    }
}

impl<A: Codec, const N: usize, const W: usize> PartialEq<SeqSlice<A>> for SeqArray<A, N, W> {
    fn eq(&self, other: &SeqSlice<A>) -> bool {
        self.as_ref() == other
    }
}

impl<A: Codec, const N: usize, const W: usize> PartialEq<SeqSlice<A>> for &SeqArray<A, N, W> {
    fn eq(&self, other: &SeqSlice<A>) -> bool {
        self.as_ref() == other
    }
}

impl<A: Codec, const N: usize, const W: usize> PartialEq<SeqArray<A, N, W>> for &SeqSlice<A> {
    fn eq(&self, other: &SeqArray<A, N, W>) -> bool {
        *self == other.as_ref()
    }
}

impl<A: Codec, const N: usize, const W: usize> PartialEq<SeqArray<A, N, W>> for SeqSlice<A> {
    fn eq(&self, other: &SeqArray<A, N, W>) -> bool {
        self == other.as_ref()
    }
}

impl<A: Codec, const N: usize, const W: usize> AsRef<SeqSlice<A>> for SeqArray<A, N, W> {
    fn as_ref(&self) -> &SeqSlice<A> {
        self
    }
}

impl<A: Codec, const N: usize, const W: usize> fmt::Display for SeqArray<A, N, W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.as_ref(), f)
    }
}
