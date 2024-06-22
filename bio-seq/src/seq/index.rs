use core::ops::{
    Index, IndexMut, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
};
use core::ptr;

use crate::codec::Codec;
use crate::seq::{Seq, SeqSlice};
use crate::Bs;

impl<A: Codec> Index<Range<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        self.as_ref().index(range)
    }
}

impl<A: Codec> Index<RangeTo<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        self.as_ref().index(range)
    }
}

impl<A: Codec> Index<RangeToInclusive<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        self.as_ref().index(range)
    }
}

impl<A: Codec> Index<RangeInclusive<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        self.as_ref().index(range)
    }
}

impl<A: Codec> Index<RangeFrom<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        self.as_ref().index(range)
    }
}

impl<A: Codec> Index<RangeFull> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFull) -> &Self::Output {
        self.as_ref().index(range)
    }
}

impl<A: Codec> Index<usize> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        self.as_ref().index(i)
    }
}

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
        let bs: *const Bs = ptr::from_ref::<Bs>(&self.bs);
        unsafe { &*(bs as *const SeqSlice<A>) }
    }
}

impl<A: Codec> Index<usize> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        &self[i..=i]
    }
}

impl<A: Codec> IndexMut<usize> for Seq<A> {
    fn index_mut(&mut self, _index: usize) -> &mut SeqSlice<A> {
        unimplemented!()
    }
}
