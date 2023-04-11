//! Amino acids
///
/// With the 6-bit representation a `Seq<Dna>` sequence can be read as an amino acid
/// sequence in place.
use crate::codec::Codec;
use core::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Codec)]
#[width = 6]
#[repr(u8)]
pub enum Amino {
    #[alt(0b110110, 0b010110, 0b100110)]
    A = 0b000110, // GCA
    #[alt(0b111011)]
    C = 0b011011, // TGC
    #[alt(0b110010)]
    D = 0b010010, // GAC
    #[alt(0b100010)]
    E = 0b000010, // GAA
    #[alt(0b111111)]
    F = 0b011111, // TTC
    #[alt(0b101010, 0b011010, 0b111010)]
    G = 0b001010, // GGA
    #[alt(0b110001)]
    H = 0b010001, // CAC
    #[alt(0b011100, 0b111100)]
    I = 0b001100, // ATA
    #[alt(0b100000)]
    K = 0b000000, // AAA
    #[alt(0b001111, 0b101111, 0b111101, 0b011101, 0b101101)]
    L = 0b001101, // CTA
    M = 0b101100, // ATG
    #[alt(0b110000)]
    N = 0b010000, // AAC
    #[alt(0b010101, 0b100101, 0b110101)]
    P = 0b000101, // CCA
    #[alt(0b100001)]
    Q = 0b000001, // CAA
    #[alt(0b101000, 0b111001, 0b011001, 0b001001, 0b101001)]
    R = 0b001000, // AGA
    #[alt(0b110111, 0b010111, 0b000111, 0b100111, 0b111000)]
    S = 0b011000, // AGC
    #[alt(0b110100, 0b010100, 0b100100)]
    T = 0b000100, // ACA
    #[alt(0b011110, 0b111110, 0b101110)]
    V = 0b001110, // GTA
    W = 0b101011, // TGG
    #[alt(0b110011)]
    Y = 0b010011, // TAC
    #[display = '*']
    #[alt(0b001011, 0b100011)]
    X = 0b000011, // TAA (stop)
}

impl From<Amino> for u8 {
    fn from(amino: Amino) -> Self {
        amino as u8
    }
}

impl fmt::Display for Amino {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

#[macro_export]
macro_rules! amino {
    ($seq:expr) => {
        match Seq::<Amino>::from_str($seq) {
            Ok(s) => s,
            Err(_) => panic!(),
        }
    };
}
