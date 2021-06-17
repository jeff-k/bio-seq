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

    fn from_bits(b: usize) -> Self {
        match b {
            0b00 => Dna::A,
            0b01 => Dna::C,
            0b10 => Dna::G,
            0b11 => Dna::T,
            _ => panic!(),
        }
    }
}
