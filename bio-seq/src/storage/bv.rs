use bitvec::prelude::*;

use core::ops::{Deref, Index, Range, RangeBounds};
use core::ptr;

type Order = Lsb0;
pub(crate) type Bs = BitSlice<usize, Order>;
type Bv = BitVec<usize, Order>;
pub(crate) type Ba<const W: usize> = BitArray<[usize; W], Order>;

use bitvec::field::BitField;
use bitvec::view::BitView;

use crate::storage::{SeqSliceStorage, SeqStorage};

const _: () = {
    assert!(core::mem::size_of::<&BitSliceStorage>() == core::mem::size_of::<&Bs>());
};

#[derive(Clone, Debug, PartialEq, Eq)]
#[repr(transparent)]
pub struct BitVecStorage {
    pub(crate) bv: Bv,
}

#[derive(Debug)]
#[repr(transparent)]
pub struct BitSliceStorage {
    pub(crate) bs: Bs,
}

impl Index<Range<usize>> for BitSliceStorage {
    type Output = Self;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        let bs = &self.bs[index];
        // SAFETY: BitSliceStorage is a repr(transparent) of BitSlice
        unsafe { &*(ptr::from_ref(bs) as *const BitSliceStorage) }
    }
}

impl SeqSliceStorage for BitSliceStorage {
    type Owned = BitVecStorage;

    fn get(&self, start: usize, end: usize) -> u8 {
        self.bs[start..end].load_le::<u8>()
    }

    fn bits(&self) -> usize {
        self.bs.len()
    }

    fn is_empty(&self) -> bool {
        self.bs.is_empty()
    }
}

impl PartialEq for BitSliceStorage {
    fn eq(&self, other: &Self) -> bool {
        self.bs == other.bs
    }
}

impl Eq for BitSliceStorage {}

impl SeqStorage for BitVecStorage {
    //    type Array<const A: usize, const B: usize> = ();
    type Unit = usize;
    type Slice = BitSliceStorage;

    fn new() -> Self {
        BitVecStorage { bv: Bv::new() }
    }

    fn bits(&self) -> usize {
        todo!()
    }

    fn with_capacity(len: usize) -> Self {
        BitVecStorage {
            bv: Bv::with_capacity(len),
        }
    }

    fn push(&mut self, byte: u8, bits: usize) {
        let bs = &byte.view_bits::<Order>()[..bits];
        self.bv.extend_from_bitslice(&bs[..bits]);
    }

    fn clear(&mut self) {
        self.bv.clear();
    }

    fn truncate(&mut self, len: usize) {
        self.bv.truncate(len);
    }

    fn extend(&mut self, other: &Self::Slice) {
        self.bv.extend_from_bitslice(&other.bs);
    }

    fn prepend(&mut self, other: &Self::Slice) {
        let mut bv = Bv::with_capacity(self.bits() + other.bs.len());
        bv.extend(&other.bs);
        bv.extend_from_bitslice(&self.bv);
        self.bv = bv;
    }

    fn drain(&mut self, range: Range<usize>) {
        self.bv.drain(range);
    }

    fn splice<R: RangeBounds<usize>>(&mut self, range: R, other: &Self::Slice) {
        self.bv.splice(range, other.bs.iter().by_vals());
    }

    fn insert(&mut self, index: usize, other: &Self::Slice) {
        assert!(index <= self.bits(), "Index out of bounds");

        let mut bv = Bv::with_capacity(self.bv.len() + other.bs.len());

        let bs = self.bv.as_bitslice();
        bv.extend_from_bitslice(&bs[..index]);
        bv.extend_from_bitslice(&other.bs);
        bv.extend_from_bitslice(&bs[index..]);

        self.bv = bv;
    }

    fn pop_unit(&self) -> Self::Unit {
        self.bs.load_le::<usize>()
    }
}

impl Deref for BitVecStorage {
    type Target = BitSliceStorage;

    fn deref(&self) -> &Self::Target {
        let bs: *const Bs = ptr::from_ref(self.bv.as_bitslice());
        unsafe { &*(bs as *const BitSliceStorage) }
    }
}
