//! Amino acids

use std::fmt;
use std::str::FromStr;

use crate::codec::dna::Dna;
use crate::codec::{Codec, ParseBioErr};
use crate::kmer::Kmer;
use crate::Seq;
use bitvec::prelude::*;

#[derive(Debug, PartialEq)]
pub enum Amino {
    A = 0b000110,    // GCA
    C = 0b011011,    // TGC
    D = 0b010010,    // GAC
    E = 0b000010,    // GAA
    F = 0b011111,    // TTC
    G = 0b001010,    // GGA
    H = 0b010001,    // CAC
    I = 0b001100,    // ATA
    K = 0b000000,    // AAA
    L = 0b001101,    // CTA
    M = 0b101100,    // ATG
    N = 0b010000,    // AAC
    P = 0b000101,    // CCA
    Q = 0b000001,    // CAA
    R = 0b001000,    // AGA
    S = 0b011000,    // AGC
    T = 0b000100,    // ACA
    V = 0b001110,    // GTA
    W = 0b101011,    // TGG
    Y = 0b010011,    // TAC
    Stop = 0b000011, // TAA
}

impl Codec for Amino {
    const WIDTH: usize = 6;

    fn to_bits(&self) -> BitArray::<Msb0, u8> {
        match &self {
            Amino::K => BitArray::new(0b000000),
            _ => unimplemented!(),
        }
    }

    fn from_bits(b: &BitSlice::<Msb0, u8>) -> Self {
        let bs: [bool; 6] = [b[0], b[1], b[2], b[3], b[4], b[5]];
        match bs {
            [false, false, false, false, false, false] => Amino::K,
            _ => unimplemented!(),
        }
    }

    fn from_ascii(c: u8) -> Self {
        unimplemented!()
    }

    fn to_ascii(c: u8) -> Self {
        unimplemented!()
    }
}

impl From<Seq<Dna>> for Amino {
    fn from(_item: Seq<Dna>) -> Self {
        unimplemented!()
    }
}

impl From<[Dna; 3]> for Amino {
    fn from(codon: [Dna; 3]) -> Self {
        match codon {
            [Dna::A, Dna::A, Dna::A] => Amino::K,
            _ => unimplemented!(),
        }
    }
}

impl From<Kmer<3>> for Amino {
    fn from(codon: Kmer<3>) -> Self {
        Amino::from_bits(&codon.bv[0..5])
    }
}

impl FromStr for Amino {
    type Err = ParseBioErr;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Amino::A),
            "K" => Ok(Amino::K),
            _ => Err(ParseBioErr),
        }
    }
}

impl fmt::Display for Amino {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Amino::Stop => write!(f, "*"),
            _ => write!(f, "{:?}", self),
        }
    }
}
