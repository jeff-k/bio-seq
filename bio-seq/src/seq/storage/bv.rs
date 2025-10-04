use bitvec::prelude::*;

use core::ops::{Deref, Range};
use core::ptr;

type Order = Lsb0;
type Bs = BitSlice<usize, Order>;
type Bv = BitVec<usize, Order>;
type Ba<const W: usize> = BitArray<[usize; W], Order>;

use bitvec::field::BitField;
use bitvec::view::BitView;

use crate::seq::SeqStorage;

#[derive(Clone, PartialEq, Eq)]
#[repr(transparent)]
pub struct BitVecStorage {
    pub(crate) bv: Bv,
}

#[repr(transparent)]
pub struct BitSliceStorage {
    pub(crate) bs: Bs,
}

impl PartialEq for BitSliceStorage {
    fn eq(&self, other: &Self) -> bool {
        self.bs == other.bs
    }
}

impl Eq for BitSliceStorage {}

impl SeqStorage for BitVecStorage {
    //    type Slice<'a> = BitSliceStorage where Self: 'a;
    type Slice = BitSliceStorage;
    type Array<const A: usize, const B: usize> = ();

    fn new() -> Self {
        BitVecStorage { bv: Bv::new() }
    }

    fn len(&self) -> usize {
        self.bv.len()
    }

    fn with_capacity(len: usize) -> Self {
        BitVecStorage {
            bv: Bv::with_capacity(len),
        }
    }

    fn is_empty(&self) -> bool {
        self.bv.is_empty()
    }

    fn push(&mut self, bit: u8) {
        Bv::push(&mut self.bv, bit & 1 != 0);
    }

    fn to_usize(&self) -> usize {
        self.bv.load_le::<usize>()
    }

    fn clear(&mut self) {
        self.bv.clear()
    }

    fn truncate(&mut self, len: usize) {
        todo!()
    }

    fn extend(&mut self, other: &Self::Slice) {
        todo!()
    }

    fn prepend(&mut self, other: &Self::Slice) {
        todo!()
    }

    fn drain(&mut self, range: Range<usize>) {
        todo!()
    }
}

impl Deref for BitVecStorage {
    type Target = BitSliceStorage;

    fn deref(&self) -> &Self::Target {
        let bs: *const Bs = ptr::from_ref(self.bv.as_bitslice());
        unsafe { &*(bs as *const BitSliceStorage) }
    }
}
