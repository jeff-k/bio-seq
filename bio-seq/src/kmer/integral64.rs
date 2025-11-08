
use crate::kmer::storage::{REV_2BIT, sealed};
use crate::seq::storage::SeqSliceStorage;

impl sealed::KmerStorage for usize {
    const BITS: usize = usize::BITS as usize;

    fn from_slice<Q: SeqSliceStorage + ?Sized>(bs: &Q) -> Self {
        debug_assert!(
            bs.len() <= Self::BITS as usize,
            "bitslice larger than kmer storage type"
        );

        todo!()
    }

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
}
