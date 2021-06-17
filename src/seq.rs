use crate::alphabet::Alphabet;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;

pub struct Seq<A: Alphabet> {
    //<A: Alphabet> {
    bv: BitVec,
    _p: PhantomData<A>,
    //width: usize,
}

impl<A: Alphabet> Seq<A> {
    pub fn from_vec(vec: Vec<A>) -> Self {
        let mut bv: BitVec = BitVec::new();
        for b in vec.iter() {
            bv.extend(b.to_bits());
        }
        Seq {
            bv: bv,
            _p: PhantomData,
        }
    }

    //    pub fn from_string(s: String) -> Self { }

    pub fn to_usize(&self) -> usize {
        self.bv.as_raw_slice()[0]
    }
}

impl<A: Alphabet + fmt::Debug> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut v = Vec::new();
        let w = A::width();
        for i in 0..(self.bv.len() / A::width()) {
            v.push(A::from_bits(&self.bv[(i * w)..((i * w) + w)]));
        }
        write!(
            f,
            "{:?} {:?}",
            v,
            A::from_bits(&self.bv[(4 * w)..((4 * w) + w)])
        )
    }
}

macro_rules! dna {
    [$seq:literal] => {
       Seq::from_vec($seq)
    };
}

macro_rules! iupac {
    [$seq:literal] => {
       Seq::from_vec($seq)
    };
}

macro_rules! amino {
    [$seq:literal] => {
       Seq::from_vec($seq)
    };
}

//impl Iterator for Seq<A> {
//    fn next(&mut self) -> Option<A> {
//        A::from_bits(self.bv.pop(A::width()))
//    }
//}
