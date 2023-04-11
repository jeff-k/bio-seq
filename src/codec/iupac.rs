//! IUPAC nucleotide ambiguity codes

/// IUPAC  nucleotide ambiguity codes are represented with 4 bits
///
/// |   | A | C | G | T |
/// | - | - | - | - | - |
/// | S | 0 | 1 | 1 | 0 |
/// | - | 0 | 0 | 0 | 0 |
/// | C | 0 | 1 | 0 | 0 |
/// | N | 1 | 1 | 1 | 1 |
/// | B | 0 | 1 | 1 | 1 |
///
/// This supports membership resolution with bitwise operations:
///
/// ```rust
/// use bio_seq::prelude::*;
///
/// // Set union:
/// assert_eq!(iupac!("AS-GYTNA") | iupac!("ANTGCAT-"), iupac!("ANTGYWNA"));
///
/// // Set intersection:
/// assert_eq!(iupac!("ACGTSWKM") & iupac!("WKMSTNNA"), iupac!("A----WKA"));
/// ```
use crate::codec::{dna::Dna, Codec};
use crate::seq::{Seq, SeqSlice};
use crate::ParseBioError;

use core::fmt;
use core::ops::{BitAnd, BitOr};
use core::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Codec)]
#[width = 4]
#[repr(u8)]
pub enum Iupac {
    A = 0b1000,
    C = 0b0100,
    G = 0b0010,
    T = 0b0001,
    R = 0b1010,
    Y = 0b0101,
    S = 0b0110,
    W = 0b1001,
    K = 0b0011,
    M = 0b1100,
    B = 0b0111,
    D = 0b1011,
    H = 0b1101,
    V = 0b1110,
    N = 0b1111,
    #[display = '-']
    X = 0b0000,
}

impl TryFrom<char> for Iupac {
    type Error = ParseBioError;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        Iupac::from_char(c)
    }
}

impl From<Iupac> for char {
    fn from(iupac: Iupac) -> char {
        iupac.to_char()
    }
}

impl From<Dna> for Iupac {
    fn from(dna: Dna) -> Self {
        match dna {
            Dna::A => Iupac::A,
            Dna::C => Iupac::C,
            Dna::G => Iupac::G,
            Dna::T => Iupac::T,
        }
    }
}

impl From<u8> for Iupac {
    fn from(b: u8) -> Self {
        Iupac::unsafe_from_bits(b)
    }
}

impl From<Iupac> for u8 {
    fn from(iupac: Iupac) -> Self {
        iupac as u8
    }
}

impl BitAnd for Seq<Iupac> {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        self.bit_and(rhs)
    }
}

impl BitOr for Seq<Iupac> {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        self.bit_or(rhs)
    }
}

impl BitAnd for &SeqSlice<Iupac> {
    type Output = Seq<Iupac>;

    fn bitand(self, rhs: Self) -> Self::Output {
        self.bit_and(rhs)
    }
}

impl BitOr for &SeqSlice<Iupac> {
    type Output = Seq<Iupac>;

    fn bitor(self, rhs: Self) -> Self::Output {
        self.bit_or(rhs)
    }
}

impl FromStr for Iupac {
    type Err = ParseBioError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Iupac::try_from(s.as_bytes()[0] as char)
    }
}

impl fmt::Display for Iupac {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Iupac::X => write!(f, "-"),
            _ => write!(f, "{self:?}"),
        }
    }
}

impl Seq<Iupac> {
    pub fn contains(&self, rhs: &SeqSlice<Iupac>) -> bool {
        let slice: &SeqSlice<Iupac> = self;
        let intersection: &SeqSlice<Iupac> = &(slice & rhs);
        intersection == rhs
    }
}

#[macro_export]
macro_rules! iupac {
    ($seq:expr) => {
        match Seq::<Iupac>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}
