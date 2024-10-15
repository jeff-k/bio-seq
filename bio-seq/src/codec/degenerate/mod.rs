//! Experimental encodings with degenerate representations
//! Includes `N`, `n`, `.`, `-`

pub mod dna;

// TODO
//pub mod iupac;

pub use dna::Dna;
//pub use iupac::Iupac;

#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::prelude::*;
}
