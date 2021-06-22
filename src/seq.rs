/*!
Sequences of bio alphabet characters. Slicable, Boxable, Iterable.
!*/

use crate::alphabet::Alphabet;
use bitvec::prelude::*;
use std::fmt;
use std::marker::PhantomData;
use std::str::FromStr;

pub struct Seq<A: Alphabet> {
    //<A: Alphabet> {
    bv: BitVec,
    _p: PhantomData<A>,
    //width: usize,
}

pub struct SeqSlice<'a, A: Alphabet> {
    //<A: Alphabet> {
    bv: &'a BitSlice,
    _p: PhantomData<A>,
    //width: usize,
}

#[derive(Debug, Clone)]
pub struct ParseSeqErr;

impl fmt::Display for ParseSeqErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "could not parse sequence")
    }
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

impl<A: Alphabet> fmt::Display for Seq<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let w = A::width();
        //        for c in self.bv.chunks(2) {
        for i in 0..(self.bv.len() / A::width()) {
            s.push_str(&A::from_bits(&self.bv[(i * w)..((i * w) + w)]).to_string());

            //            v.push(c.as_raw_slice()[0]);
        }
        write!(
            f,
            "{}",
            s,
            //&self.bv[(4 * w)..((4 * w) + w)][2],
        )
    }
}

impl<A: Alphabet> FromStr for Seq<A> {
    type Err = ParseSeqErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut v = Vec::new();
        let w = A::width();
        for i in s.chars() {
            match A::from_str(&i.to_string()) {
                Ok(b) => v.push(b),
                Err(_) => panic!(),
            }
            //v.push(A::from_str(&i.to_string())?);
        }
        //Ok(Self::from_vec(v))
        Ok(Seq::<A>::from_vec(v))
    }
}
//impl<A> std::ops::Index<usize> for Seq<A>
//where
//    //Idx: std::ops::Index<usize>,
//    A: Alphabet,
//{
//    type Output = SeqSlice<A>;
//
//    fn index(&self, i: usize) -> Self::Output {
//        let w = A::width();
//        SeqSlice { bv: &self.bv[i*w..(i*w)+w], _p: PhantomData }
//    }
//}

//impl Iterator for Seq<A> {
//    fn next(&mut self) -> Option<A> {
//        A::from_bits(self.bv.pop(A::width()))
//    }
//}

macro_rules! dna {
    ($seq:expr) => {
        Seq::<Dna>::from_str($seq)
    };
}

macro_rules! iupac {
    ($seq:expr) => {
        Seq::from_str($seq)
    };
}

macro_rules! amino {
    ($seq:expr) => {
        Seq::from_str($seq)
    };
}

//impl Iterator for Seq<A> {
//    fn next(&mut self) -> Option<A> {
//        A::from_bits(self.bv.pop(A::width()))
//    }
//}
