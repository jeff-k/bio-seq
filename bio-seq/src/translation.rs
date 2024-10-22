//! # Amino acid translation tables
//!
//! This module provides traits for implementing amino acid translation tables.
//!
//! Enable the translation feature in `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! bio-seq = { version="0.13", features=["translation"] }
//! ```
//!
//! ## Examples
//!
//! The standard genetic code is provided as a `translation::STANDARD` constant:
//!
//! ```rust
//! use bio_seq::prelude::*;
//! use bio_seq::translation::STANDARD;
//! use bio_seq::translation::TranslationTable;
//!
//! let seq = dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
//!
//! let aminos: Seq<Amino> = seq
//!     .windows(3)
//!     .map(|codon| STANDARD.to_amino(&codon))
//!     .collect::<Seq<Amino>>();
//!
//! assert_eq!(
//!     aminos,
//!     Seq::<Amino>::try_from("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK").unwrap()
//! );
//!
//! ```
//!
//! Custom translation tables can be implemented from associative datastructures:
//!
//! ```
//! use bio_seq::prelude::*;
//! use bio_seq::translation::{TranslationTable, TranslationError};
//!
//! struct Mitochondria;
//! impl TranslationTable<Dna, Amino> for Mitochondria {
//!     fn to_amino(&self, codon: &SeqSlice<Dna>) -> Amino {
//!         if codon == dna!("AGA") {
//!             Amino::X
//!         } else if codon == dna!("AGG") {
//!             Amino::X
//!         } else if codon == dna!("ATA") {
//!             Amino::M
//!        } else if codon == dna!("TGA") {
//!             Amino::W
//!         } else {
//!                 Amino::unsafe_from_bits(Into::<u8>::into(codon))
//!               }
//!           }
//!
//!          fn to_codon(&self, _amino: Amino) -> Result<Seq<Dna>, TranslationError> {
//!               unimplemented!()
//!           }
//!       }
//!
//!        let seq: Seq<Dna> =
//!            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA").into();
//!        let aminos: Seq<Amino> = seq
//!            .windows(3)
//!            .map(|codon| Mitochondria.to_amino(&codon))
//!            .collect::<Seq<Amino>>();
//!        assert_eq!(seq.len() - 2, aminos.len());
//!
//!        for (x, y) in aminos.into_iter().zip(
//!            Seq::<Amino>::try_from(
//!                "NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYMWE*GDTRDISQSAHTPHM*K*ENTQK",
//!            )
//!            .unwrap()
//!            .into_iter(),
//!        ) {
//!            assert_eq!(x, y)
//!        }
//! ```
//!
//! ## Errors
//!
//! Translation tables may not be complete or they may be ambiguous
//!
use core::cmp::Eq;
use core::fmt;
use std::collections::HashMap;

use crate::codec::Codec;
use crate::prelude::{Amino, Dna, Seq, SeqSlice};

mod standard;

pub use crate::translation::standard::STANDARD;

/// Error conditions for codon/amino acid translation
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum TranslationError<A: Codec = Dna, B: Codec = Amino> {
    /// Amino acid can be translation from multiple codons
    AmbiguousCodon(B),
    /// Codon sequence maps to multiple amino acids
    AmbiguousTranslation(Seq<A>),
    /// Codon sequence does not map to an amino acid
    InvalidCodon(Seq<A>),
    /// Amino acid symbol is not valid (i.e. `X`)
    InvalidAmino(B),
}

impl<A: Codec, B: Codec> fmt::Display for TranslationError<A, B> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TranslationError::AmbiguousCodon(amino) => {
                let amino = amino.to_char();
                write!(f, "Multiple codon sequences: {amino}")
            }
            TranslationError::AmbiguousTranslation(codon) => {
                write!(f, "Ambiguous translations for codon: {codon}")
            }
            TranslationError::InvalidCodon(codon) => write!(f, "Invalid codon sequence: {codon}"),
            TranslationError::InvalidAmino(amino) => {
                let amino = amino.to_char();
                write!(f, "Invalid amino acid character: {amino}")
            }
        }
    }
}

// #![feature(error_in_core)
impl<A: Codec, B: Codec> std::error::Error for TranslationError<A, B> {}

/// A codon translation table where all codons map to amino acids
pub trait TranslationTable<A: Codec, B: Codec> {
    fn to_amino(&self, codon: &SeqSlice<A>) -> B;

    /// # Errors
    ///
    /// Will return `Err` when an amino acid has multiple codons (most cases)
    fn to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError<A, B>>;
}

/// A partial translation table where not all triples of characters map to amino acids
pub trait PartialTranslationTable<A: Codec, B: Codec> {
    /// # Errors
    ///
    /// Will return an `Err` if a codon does not map to an amino acid. This would be
    /// the case for a translation table from codons with ambiguous nucleotide codes such as `ANC`, `SWS`, `NNN`, etc.
    fn try_to_amino(&self, codon: &SeqSlice<A>) -> Result<B, TranslationError<A, B>>;
    /// # Errors
    ///
    /// Will return an `Err` if the amino acid can be translated from different codons
    fn try_to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError<A, B>>;
}

/// A customisable translation table
pub struct CodonTable<A: Codec, B: Codec> {
    // I'm open to using a better bidirectional mapping datastructure
    table: HashMap<Seq<A>, B>,
    inverse_table: HashMap<B, Option<Seq<A>>>,
}

impl<A: Codec, B: Codec> CodonTable<A, B> {
    pub fn from_map<T>(table: T) -> Self
    where
        T: Into<HashMap<Seq<A>, B>>,
    {
        let table: HashMap<Seq<A>, B> = table.into();
        let mut inverse_table = HashMap::new();
        for (codon, amino) in &table {
            if inverse_table.contains_key(amino) {
                inverse_table.insert(*amino, None);
            } else {
                inverse_table.insert(*amino, Some(codon.clone()));
            }
        }
        CodonTable {
            table,
            inverse_table,
        }
    }
}

impl<A: Codec, B: Codec> PartialTranslationTable<A, B> for CodonTable<A, B> {
    fn try_to_amino(&self, codon: &SeqSlice<A>) -> Result<B, TranslationError<A, B>> {
        self.table
            .get(codon)
            .ok_or_else(|| TranslationError::InvalidCodon(codon.into()))
            .copied()
    }

    fn try_to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError<A, B>> {
        if let Some(codon) = self.inverse_table.get(&amino) {
            match codon {
                Some(codon) => Ok(codon.clone()),
                None => Err(TranslationError::AmbiguousCodon(amino)),
            }
        } else {
            Err(TranslationError::InvalidAmino(amino))
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::{
        CodonTable, PartialTranslationTable, TranslationError, TranslationTable,
    };

    #[test]
    fn custom_codon_table() {
        let mito: [(Seq<Dna>, Amino); 6] = [
            (dna!("AAA").into(), Amino::A),
            (dna!("ATG").into(), Amino::A),
            (dna!("CCC").into(), Amino::C),
            (dna!("GGG").into(), Amino::E),
            (dna!("TTT").into(), Amino::D),
            (dna!("TTA").into(), Amino::F),
        ];

        let table = CodonTable::from_map(mito);

        let seq: Seq<Dna> = dna!("AAACCCGGGTTTTTATTAATG").into();
        let mut amino_seq: Seq<Amino> = Seq::new();
        for codon in seq.chunks(3) {
            amino_seq.push(table.try_to_amino(codon).unwrap());
        }
        assert_eq!(amino_seq, Seq::<Amino>::try_from("ACEDFFA").unwrap());

        assert_ne!(table.try_to_codon(Amino::E), Ok(dna!("CCC").into()));
        assert_eq!(table.try_to_codon(Amino::C), Ok(dna!("CCC").into()));
        assert_eq!(
            table.try_to_codon(Amino::A),
            Err(TranslationError::AmbiguousCodon(Amino::A))
        );
        assert_eq!(
            table.try_to_codon(Amino::X),
            Err(TranslationError::InvalidAmino(Amino::X))
        );
    }

    #[test]
    fn mitochondrial_coding_table() {
        struct Mitochondria;

        impl TranslationTable<Dna, Amino> for Mitochondria {
            fn to_amino(&self, codon: &SeqSlice<Dna>) -> Amino {
                if codon == dna!("AGA") {
                    Amino::X
                } else if codon == dna!("AGG") {
                    Amino::X
                } else if codon == dna!("ATA") {
                    Amino::M
                } else if codon == dna!("TGA") {
                    Amino::W
                } else {
                    Amino::unsafe_from_bits(Into::<u8>::into(codon))
                }
            }

            fn to_codon(&self, _amino: Amino) -> Result<Seq<Dna>, TranslationError> {
                unimplemented!()
            }
        }

        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA").into();
        let aminos: Seq<Amino> = seq
            .windows(3)
            .map(|codon| Mitochondria.to_amino(&codon))
            .collect::<Seq<Amino>>();
        assert_eq!(seq.len() - 2, aminos.len());

        for (x, y) in aminos.into_iter().zip(
            Seq::<Amino>::try_from(
                "NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYMWE*GDTRDISQSAHTPHM*K*ENTQK",
            )
            .unwrap()
            .into_iter(),
        ) {
            assert_eq!(x, y)
        }
    }
}
