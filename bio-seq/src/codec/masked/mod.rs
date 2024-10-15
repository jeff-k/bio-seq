//! Experimental encodings with maskable bases

pub mod dna;
pub mod iupac;

pub use dna::Dna;
pub use iupac::Iupac;

pub trait Maskable {
    fn mask(&self) -> Self;
    fn unmask(&self) -> Self;
}

#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::prelude::*;
}
