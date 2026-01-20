// Copyright 2025, 2026 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

mod bv;

#[cfg(target_pointer_width = "64")]
pub(crate) mod integral64;

use core::ops::{Deref, Index, Range, RangeBounds};

pub(crate) use bv::{Ba, Bs};
pub use bv::{BitSliceStorage, BitVecStorage};

pub trait PrimitiveStorage: Sized + Clone + PartialEq {
    // Deref<Target = Self::Slice> {
    const BITS: usize;
    type Slice: ?Sized + SeqSliceStorage;

    fn bits(&self) -> usize;
    fn is_empty(&self) -> bool;

    fn to_usize(&self) -> usize;

    fn shiftr(&mut self, n: u32);

    fn shiftl(&mut self, n: u32);

    fn mask(&mut self, bits: Self);

    fn complement(&mut self, mask: Self);

    fn rev_blocks(&mut self, n: usize);

    /// Reverse bits
    fn reverse(&self) -> Self;

    fn unsafe_from_slice(slice: &Self::Slice) -> Self;
}

pub trait SeqStorage: Sized + Clone + PartialEq + Deref<Target = Self::Slice> {
    type Slice: ?Sized + SeqSliceStorage<Owned = Self>;
    type Unit: PrimitiveStorage<Slice = Self::Slice>;

    fn new() -> Self;
    fn with_capacity(cap: usize) -> Self;
    fn bits(&self) -> usize;

    fn push(&mut self, byte: u8, bits: usize);
    fn clear(&mut self);

    fn truncate(&mut self, len: usize);
    fn prepend(&mut self, other: &Self::Slice);
    fn extend(&mut self, other: &Self::Slice);
    fn drain(&mut self, range: Range<usize>);

    fn splice<R: RangeBounds<usize>>(&mut self, range: R, other: &Self::Slice);
    fn insert(&mut self, index: usize, other: &Self::Slice);

    fn pop_unit(&self) -> Self::Unit;
}

pub trait SeqSliceStorage: PartialEq + Eq + Index<Range<usize>, Output = Self> {
    type Owned: SeqStorage<Slice = Self>;
    fn get(&self, start: usize, end: usize) -> u8;
    fn bits(&self) -> usize;
    fn is_empty(&self) -> bool;
}
