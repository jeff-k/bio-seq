use crate::alphabet::Alphabet;
use bitvec::prelude::*;

#[derive(Debug)]
pub enum Dna {
    A,
    C,
    G,
    T,
}

impl Alphabet for Dna {
    fn width() -> usize {
        2
    }

    fn to_bits(&self) -> BitVec {
        match &self {
            Dna::A => bitvec![0, 0],
            Dna::C => bitvec![0, 1],
            Dna::G => bitvec![1, 0],
            Dna::T => bitvec![1, 1],
        }
    }

    fn from_bits(b: &BitSlice) -> Self {
        let bs: [bool; 2] = [b[0], b[1]];
        match bs {
            [false, false] => Dna::A,
            [false, true] => Dna::C,
            [true, false] => Dna::G,
            [true, true] => Dna::T,
            _ => panic!("got bitvector: {:?}", b),
        }
    }
}
