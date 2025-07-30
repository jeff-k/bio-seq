// Copyright 2025 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::{Bs, Bv};

pub trait SeqStorage: Sized {
    type Slice<'a>: 'a + ?Sized
    where
        Self: 'a;

    fn new() -> Self;
    fn len(&self) -> usize;
    fn with_capacity(cap: usize) -> Self;
    fn is_empty(&self) -> bool;
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
}
