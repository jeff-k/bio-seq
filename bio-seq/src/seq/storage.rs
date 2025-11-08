// Copyright 2025 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

mod bv;

use core::ops::{Deref, Range, RangeBounds};

pub use bv::{BitSliceStorage, BitVecStorage};

pub trait SeqStorage: Sized + Clone + PartialEq + Deref<Target = Self::Slice> {
    //    type Slice: ?Sized + PartialEq + Eq;
    type Slice: ?Sized + SeqSliceStorage;
    type Array<const N: usize, const W: usize>;

    fn new() -> Self;
    fn len(&self) -> usize;
    fn with_capacity(cap: usize) -> Self;
    fn is_empty(&self) -> bool;

    fn to_usize(&self) -> usize;
    fn push(&mut self, byte: u8, bits: usize);
    fn clear(&mut self);

    fn truncate(&mut self, len: usize);
    fn prepend(&mut self, other: &Self::Slice);
    fn extend(&mut self, other: &Self::Slice);
    fn drain(&mut self, range: Range<usize>);

    fn splice<R: RangeBounds<usize>>(&mut self, range: R, other: &Self::Slice);
    fn insert(&mut self, index: usize, other: &Self::Slice);
}

pub trait SeqSliceStorage: PartialEq + Eq {
    fn get(&self, start: usize, end: usize) -> u8;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
}
