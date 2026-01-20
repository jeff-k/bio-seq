use core::ops::{Index, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};
use core::ptr;

use crate::codec::Codec;
use crate::seq::SeqSlice;
use crate::storage::SeqSliceStorage;

fn slice_from_storage<A: Codec, S: SeqSliceStorage + ?Sized>(slice: &S) -> &SeqSlice<A, S> {
    // SAFETY: Seqslice is repr(transparent) for the backing storage type
    unsafe { &*(ptr::from_ref(slice) as *const SeqSlice<A, S>) }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<Range<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::BITS;
        let e = range.end * A::BITS;
        slice_from_storage::<A, S>(&self.bs[s..e])
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeTo<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        let e = range.end * A::BITS;
        slice_from_storage::<A, S>(&self.bs[0..e])
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeToInclusive<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        let e = (range.end + 1) * A::BITS;
        slice_from_storage::<A, S>(&self.bs[0..e])
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeInclusive<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        let s = range.start() * A::BITS;
        let e = (range.end() + 1) * A::BITS;
        slice_from_storage::<A, S>(&self.bs[s..e])
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeFrom<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        let s = range.start * A::BITS;
        let e = self.bs.bits();
        slice_from_storage::<A, S>(&self.bs[s..e])
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeFull> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let e = self.bs.bits();
        slice_from_storage::<A, S>(&self.bs[0..e])
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<usize> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, i: usize) -> &Self::Output {
        let s = i * A::BITS;
        let e = s + A::BITS;
        slice_from_storage::<A, S>(&self.bs[s..e])
    }
}
