//use crate::bitops::REV_2BIT;
use crate::storage::PrimitiveStorage;
use crate::storage::bv;
//use bitvec::field::BitField;

impl PrimitiveStorage for usize {
    //const BITS: usize = usize::BITS as usize;

    const BITS: usize = usize::BITS as usize;
    type Slice = bv::BitSliceStorage;

    fn len(&self) -> usize {
        Self::BITS as usize
    }

    fn is_empty(&self) -> bool {
        false
    }

    fn to_usize(&self) -> usize {
        *self
    }

    fn reverse(&self) -> Self {
        todo!()
    }

    /*
        fn unsafe_from_slice<Q: SeqSliceStorage + ?Sized>(bs: &Q) -> Self {
            debug_assert!(
                bs.len() <=  usize::BITS as usize,
                "bitslice larger than kmer storage type"
            );

            todo!()
        }
    */

    fn shiftr(&mut self, _n: u32) {
        todo!()
    }

    fn shiftl(&mut self, _n: u32) {
        todo!()
    }

    fn mask(&mut self, _bits: usize) {
        todo!()
    }

    fn complement(&mut self, _mask: usize) {
        todo!()
    }

    fn rev_blocks_2(&mut self) {
        //REV_2BIT(&self);
        todo!()
    }

    fn unsafe_from_slice(_slice: &Self::Slice) -> Self {
        todo!()
    }
}
