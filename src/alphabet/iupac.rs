use crate::alphabet::Alphabet;
use bitvec::prelude::*;

#[derive(Debug)]
pub enum Iupac {
    A,
    C,
    G,
    T,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    N,
    X,
}

impl Alphabet for Iupac {
    fn width() -> usize {
        4
    }

    fn to_bits(&self) -> BitVec {
        match &self {
            Iupac::A => bitvec![1, 0, 0, 0],
            Iupac::C => bitvec![0, 1, 0, 0],
            Iupac::G => bitvec![0, 0, 1, 0],
            Iupac::T => bitvec![0, 0, 0, 1],

            Iupac::R => bitvec![1, 0, 1, 0],
            Iupac::Y => bitvec![0, 1, 0, 1],
            Iupac::S => bitvec![0, 1, 1, 0],
            Iupac::W => bitvec![1, 0, 0, 1],

            Iupac::K => bitvec![0, 0, 1, 1],
            Iupac::M => bitvec![1, 1, 0, 0],
            Iupac::B => bitvec![0, 1, 1, 1],
            Iupac::D => bitvec![1, 0, 1, 1],

            Iupac::H => bitvec![1, 1, 0, 1],
            Iupac::V => bitvec![1, 1, 1, 0],
            Iupac::N => bitvec![1, 1, 1, 1],
            Iupac::X => bitvec![0, 0, 0, 0],
        }
    }

    fn from_bits(b: &BitSlice) -> Self {
        let bs: [bool; 4] = [b[0], b[1], b[2], b[3]];
        match bs {
            [true, false, false, false] => Iupac::A,
            [false, true, false, false] => Iupac::C,
            [false, false, true, false] => Iupac::G,
            [false, false, false, true] => Iupac::T,

            [true, false, true, false] => Iupac::R,
            [false, true, false, true] => Iupac::Y,
            [false, true, true, false] => Iupac::S,
            [true, false, false, true] => Iupac::W,

            [false, false, true, true] => Iupac::K,
            [true, true, false, false] => Iupac::M,
            [false, true, true, true] => Iupac::B,
            [true, false, true, true] => Iupac::D,

            [true, true, false, true] => Iupac::H,
            [true, true, true, false] => Iupac::V,
            [true, true, true, true] => Iupac::N,
            [false, false, false, false] => Iupac::X,
        }
    }
}
