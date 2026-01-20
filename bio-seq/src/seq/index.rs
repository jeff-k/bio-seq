use core::ops::{Index, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};
//use core::ptr;

use crate::storage::SeqSliceStorage;
//use crate::storage::BitSliceStorage;

//use crate::Bs;
use crate::codec::Codec;
use crate::seq::SeqSlice;

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<Range<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: Range<usize>) -> &Self::Output {
        todo!()
        /*
        let s = range.start * A::BITS as usize;
        let e = range.end * A::BITS as usize;
        let bs = &self.bs.bs[s..e];
        unsafe { &*(ptr::from_ref(bs) as *const SeqSlice<A, S>) }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeTo<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: RangeTo<usize>) -> &Self::Output {
        todo!()
        /*
        let e = range.end * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeToInclusive<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: RangeToInclusive<usize>) -> &Self::Output {
        todo!()
        /*
        let e = (range.end + 1) * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeInclusive<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: RangeInclusive<usize>) -> &Self::Output {
        todo!()
        /*
        let s = range.start() * A::BITS as usize;
        let e = (range.end() + 1) * A::BITS as usize;

        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeFrom<usize>> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: RangeFrom<usize>) -> &Self::Output {
        todo!()
        /*
        let s = range.start * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..]);
        unsafe { &*(bs as *const SeqSlice<A>) }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<RangeFull> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        todo!()
        /*
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[..]);
        unsafe { &*(bs as *const SeqSlice<A>) }
        */
    }
}

impl<A: Codec, S: SeqSliceStorage + ?Sized> Index<usize> for SeqSlice<A, S> {
    type Output = SeqSlice<A, S>;

    fn index(&self, _i: usize) -> &Self::Output {
        todo!()
        /*
        let s = i * A::BITS as usize;
        let e = s + A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
        */
    }
}
