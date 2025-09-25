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
        self.bv.load_le::<usize>()
    }

    fn clear(&mut self) {
        todo!()
    }

}

