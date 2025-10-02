use bitvec::prelude::*;

type Order = Lsb0;
type Bs = BitSlice<usize, Order>;
type Bv = BitVec<usize, Order>;
type Ba<const W: usize> = BitArray<[usize; W], Order>;

use bitvec::field::BitField;
use bitvec::view::BitView;

use crate::seq::SeqStorage;

#[derive(Clone, PartialEq, Eq)]
pub struct BitVecStorage {
    pub(crate) bv: Bv,
}

pub struct BitSliceStorage<'a> {
    pub(crate) bs: &'a Bs,
}


impl SeqStorage for BitVecStorage {
    type Slice<'a> = BitSliceStorage<'a>;
    type Array<const A: usize, const B: usize> = ();

    fn new() -> Self {
        BitVecStorage { bv: Bv::new() }
    }

    fn len(&self) -> usize {
        self.bv.len()
    }

    fn with_capacity(len: usize) -> Self {
        BitVecStorage { bv: Bv::with_capacity(len) }
    }

    fn is_empty(&self) -> bool {
        self.bv.is_empty()
    }

    fn as_slice(&self) -> &Self::Slice<'_> {
        todo!()
        //&Self::Slice { bs: self.bv.as_slice() }
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

}

