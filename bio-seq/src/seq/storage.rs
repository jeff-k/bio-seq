// Copyright 2025 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::{Bs, Bv};
use bitvec::field::BitField;

pub trait SeqStorage: Sized + Clone + PartialEq {
    type Slice<'a>: 'a + ?Sized
    where
        Self: 'a;

    fn new() -> Self;
    fn len(&self) -> usize;
    fn with_capacity(cap: usize) -> Self;
    fn is_empty(&self) -> bool;

    fn as_slice(&self) -> &Self::Slice<'_>;
    fn to_usize(&self) -> usize;
    fn push(&mut self, bits: u8);
    fn clear(&mut self);
}

impl SeqStorage for Bv {
    type Slice<'a> = Bs;

    fn new() -> Self {
        todo!()
    }

    fn len(&self) -> usize {
        todo!()
    }

    fn with_capacity(_len: usize) -> Self {
        todo!()
    }

    fn is_empty(&self) -> bool {
        todo!()
    }

    fn as_slice(&self) -> &Self::Slice<'_> {
        todo!()
    }

    fn push(&mut self, _bits: u8) {
        todo!()
    }

    fn to_usize(&self) -> usize {
        self.load_le::<usize>()
    }

    fn clear(&mut self) {
        todo!()
    }

}

impl SeqStorage for Vec<u8> {
    type Slice<'a> = &'a [u8];

    fn new() -> Self {
        todo!()
    }

    fn len(&self) -> usize {
        todo!()
    }

    fn with_capacity(_len: usize) -> Self {
        todo!()
    }

    fn is_empty(&self) -> bool {
        todo!()
    }

    fn as_slice(&self) -> &Self::Slice<'_> {
        todo!()
    }

    fn push(&mut self, _bits: u8) {
        todo!()
    }

   fn to_usize(&self) -> usize {
       todo!()
    }

    fn clear(&mut self) {
        todo!()
    }
}
