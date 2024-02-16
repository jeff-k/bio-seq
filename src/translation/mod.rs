//! # Genetic Code Translation
//!
//! This module provides traits
//!
//! ## Examples
//!
//! ## Errors
//!
use core::fmt::Display;
use core::iter::FromIterator;
use std::collections::HashMap;

use crate::codec::Codec;
use crate::error::TranslationError;
use crate::prelude::{Seq, SeqSlice};

pub mod standard;

/// A codon translation table where all codons map to amino acids
trait TranslationTable<A: Codec, B: Codec> {
    fn to_amino(&self, codon: &SeqSlice<A>) -> B;
    fn to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError>;
}

/// A partial translation table where not all triples of characters map to amino acids
trait PartialTranslationTable<A: Codec, B: Codec> {
    fn try_to_amino(&self, codon: &SeqSlice<A>) -> Result<B, TranslationError>;
    fn try_to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError>;
}

/// A customisable translation table
pub struct CodonTable<A: Codec, B: Codec> {
    table: HashMap<Seq<A>, B>,
    inverse_table: HashMap<B, Option<Seq<A>>>,
}

impl<A: Codec, B: Codec + Display> CodonTable<A, B> {
    pub fn from_map(table: HashMap<Seq<A>, B>) -> Self {
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

impl<A: Codec, B: Codec> FromIterator<(Seq<A>, B)> for CodonTable<A, B> {
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = (Seq<A>, B)>,
    {
        let mut table: HashMap<Seq<A>, B> = HashMap::new();
        let mut inverse_table: HashMap<B, Option<Seq<A>>> = HashMap::new();

        for (codon, amino) in iter {
            table.insert(codon.clone(), amino);
            inverse_table.insert(amino, Some(codon.clone()));
        }

        CodonTable {
            table,
            inverse_table,
        }
    }
}

impl<A: Codec, B: Codec + Display> PartialTranslationTable<A, B> for CodonTable<A, B> {
    fn try_to_amino(&self, codon: &SeqSlice<A>) -> Result<B, TranslationError> {
        match self.table.get(&Seq::from(codon)) {
            Some(amino) => Ok(*amino),
            None => Err(TranslationError::InvalidCodon),
        }
    }

    fn try_to_codon(&self, amino: B) -> Result<Seq<A>, TranslationError> {
        match self.inverse_table.get(&amino) {
            Some(Some(codon)) => Ok(codon.clone()),
            Some(None) => Err(TranslationError::AmbiguousCodon),
            None => Err(TranslationError::InvalidCodon),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::{CodonTable, PartialTranslationTable, TranslationError};
    use std::collections::HashMap;

    #[test]
    fn custom_codon_table() {
        let mito: [(Seq<Dna>, Amino); 6] = [
            (dna!("AAA"), Amino::A),
            (dna!("ATG"), Amino::A),
            (dna!("CCC"), Amino::C),
            (dna!("GGG"), Amino::E),
            (dna!("TTT"), Amino::D),
            (dna!("TTA"), Amino::F),
        ];

        let codons = HashMap::from(mito);
        let table = CodonTable::from_map(codons);

        let seq: Seq<Dna> = dna!("AAACCCGGGTTTTTATTAATG");
        let mut amino_seq: Seq<Amino> = Seq::new();
        for codon in seq.chunks(3) {
            amino_seq.push(table.try_to_amino(codon).unwrap());
        }
        assert_eq!(amino_seq, amino!("ACEDFFA"));

        assert_ne!(table.try_to_codon(Amino::E), Ok(dna!("CCC")));
        assert_eq!(table.try_to_codon(Amino::C), Ok(dna!("CCC")));
        assert_eq!(
            table.try_to_codon(Amino::A),
            Err(TranslationError::AmbiguousCodon)
        );
        assert_eq!(
            table.try_to_codon(Amino::X),
            Err(TranslationError::InvalidCodon)
        );
    }
}
