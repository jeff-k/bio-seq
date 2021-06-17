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
        let a = bits![0, 1];
        let c = bits![0, 1];
        let g = bits![1, 0];
        let t = bits![1, 1];
        match b {
            //[0, 0] => Dna::A,
            //[0, 1] => Dna::C,
            g => Dna::G,
            t => Dna::T,
        }
    }
}
