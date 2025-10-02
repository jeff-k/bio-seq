// Copyright 2025 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

mod bv;

use core::ops::Deref;

pub use bv::{BitSliceStorage, BitVecStorage};

pub trait SeqStorage: Sized + Clone + PartialEq + Deref {
    type Array<const N: usize, const W: usize>;

    fn new() -> Self;
    fn len(&self) -> usize;
    fn with_capacity(cap: usize) -> Self;
    fn is_empty(&self) -> bool;

    fn to_usize(&self) -> usize;
    fn push(&mut self, bits: u8);
    fn clear(&mut self);
}
