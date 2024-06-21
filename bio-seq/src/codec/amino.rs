//! 6-bit representation of amino acids
//!
//! A `Seq<Dna>` sequence of length 3 can be read as an amino acid
//! sequence in place.
use crate::codec::Codec;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Codec)]
#[bits(6)]
#[repr(u8)]
pub enum Amino {
    #[alt(0b11_01_10, 0b01_01_10, 0b10_01_10)]
    A = 0b00_01_10, // GCA
    #[alt(0b11_10_11)]
    C = 0b01_10_11, // TGC
    #[alt(0b11_00_10)]
    D = 0b01_00_10, // GAC
    #[alt(0b10_00_10)]
    E = 0b00_00_10, // GAA
    #[alt(0b11_11_11)]
    F = 0b01_11_11, // TTC
    #[alt(0b10_10_10, 0b01_10_10, 0b11_10_10)]
    G = 0b00_10_10, // GGA
    #[alt(0b11_00_01)]
    H = 0b01_00_01, // CAC
    #[alt(0b01_11_00, 0b11_11_00)]
    I = 0b00_11_00, // ATA
    #[alt(0b10_00_00)]
    K = 0b00_00_00, // AAA
    #[alt(0b00_11_11, 0b10_11_11, 0b11_11_01, 0b01_11_01, 0b10_11_01)]
    L = 0b00_11_01, // CTA
    M = 0b10_11_00, // ATG
    #[alt(0b11_00_00)]
    N = 0b01_00_00, // AAC
    #[alt(0b01_01_01, 0b10_01_01, 0b11_01_01)]
    P = 0b00_01_01, // CCA
    #[alt(0b10_00_01)]
    Q = 0b00_00_01, // CAA
    #[alt(0b10_10_00, 0b11_10_01, 0b01_10_01, 0b00_10_01, 0b10_10_01)]
    R = 0b00_10_00, // AGA
    #[alt(0b11_01_11, 0b01_01_11, 0b00_01_11, 0b10_01_11, 0b11_10_00)]
    S = 0b01_10_00, // AGC
    #[alt(0b11_01_00, 0b01_01_00, 0b10_01_00)]
    T = 0b00_01_00, // ACA
    #[alt(0b01_11_10, 0b11_11_10, 0b10_11_10)]
    V = 0b00_11_10, // GTA
    W = 0b10_10_11, // TGG
    #[alt(0b11_00_11)]
    Y = 0b01_00_11, // TAC
    #[display('*')]
    #[alt(0b00_10_11, 0b10_00_11)]
    X = 0b00_00_11, // TAA (stop)
}

impl From<Amino> for u8 {
    fn from(b: Amino) -> u8 {
        b as u8
    }
}

impl core::fmt::Display for Amino {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
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
