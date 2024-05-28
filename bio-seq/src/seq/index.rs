use core::ops::{
    Index, IndexMut, Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
};

use crate::codec::Codec;
use crate::seq::{Seq, SeqSlice};
use crate::Bs;

impl<A: Codec> Index<Range<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::BITS as usize;
        let e = range.end * A::BITS as usize;
        let bs = &self.bv[s..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeTo<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        let e = range.end * A::BITS as usize;
        let bs = &self.bv[..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeToInclusive<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        let e = (range.end + 1) * A::BITS as usize;
        let bs = &self.bv[..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeInclusive<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        let s = range.start() * A::BITS as usize;
        let e = (range.end() + 1) * A::BITS as usize;
        let bs = &self.bv[s..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFrom<usize>> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        let s = range.start * A::BITS as usize;
        let bs = &self.bv[s..] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFull> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs = &self.bv[0..self.bv.len()] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<usize> for Seq<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        &self[i..i + 1]
    }
}

impl<A: Codec> Index<Range<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: Range<usize>) -> &Self::Output {
        let s = range.start * A::BITS as usize;
        let e = range.end * A::BITS as usize;
        let bs = &self.bs[s..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeTo<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeTo<usize>) -> &Self::Output {
        let e = range.end * A::BITS as usize;
        let bs = &self.bs[..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeToInclusive<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeToInclusive<usize>) -> &Self::Output {
        let e = (range.end + 1) * A::BITS as usize;
        let bs = &self.bs[..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeInclusive<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeInclusive<usize>) -> &Self::Output {
        let s = range.start() * A::BITS as usize;
        let e = (range.end() + 1) * A::BITS as usize;
        let bs = &self.bs[s..e] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFrom<usize>> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, range: RangeFrom<usize>) -> &Self::Output {
        let s = range.start * A::BITS as usize;
        let bs = &self.bs[s..] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<RangeFull> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, _range: RangeFull) -> &Self::Output {
        let bs = &self.bs[0..self.bs.len()] as *const Bs as *const SeqSlice<A>;
        unsafe { &*bs }
    }
}

impl<A: Codec> Index<usize> for SeqSlice<A> {
    type Output = SeqSlice<A>;

    fn index(&self, i: usize) -> &Self::Output {
        //A::unsafe_from_bits(self.bs[s..e].load());
        //&self.bs[s..e].load::<u8>()
        &self[i..i + 1]
    }
}

impl<A: Codec> IndexMut<usize> for Seq<A> {
    fn index_mut(&mut self, _index: usize) -> &mut SeqSlice<A> {
        unimplemented!()
    }
}
