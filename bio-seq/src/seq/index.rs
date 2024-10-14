use core::ops::{Index, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};
use core::ptr;

use crate::codec::Codec;
use crate::seq::SeqSlice;
use crate::Bs;

impl<A: Codec> Index<Range<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::BITS as usize;
        let e = range.end * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<RangeTo<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        let e = range.end * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<RangeToInclusive<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        let e = (range.end + 1) * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<RangeInclusive<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        let s = range.start() * A::BITS as usize;
        let e = (range.end() + 1) * A::BITS as usize;

        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<RangeFrom<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        let s = range.start * A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<RangeFull> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[..]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<usize> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        let s = i * A::BITS as usize;
        let e = s + A::BITS as usize;
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs[s..e]);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}
