//! Translation tables for coding sequences
use std::collections::HashMap;
use std::iter::FromIterator;

use crate::codec::Codec;
use crate::prelude::{Seq, SeqSlice};

pub mod standard;

#[derive(Debug, PartialEq, Eq, Clone)]
enum TranslationError {
    AmbiguousCodon,
    InvalidCodon,
}

trait TranslationTable<A: Codec, B: Codec> {
    fn to_amino(self, codon: &SeqSlice<A>) -> B;
    fn to_codon(self, amino: B) -> Result<Seq<A>, TranslationError>;
}

trait PartialTranslationTable<A: Codec, B: Codec> {
    fn try_to_amino(self, codon: &SeqSlice<A>) -> Result<B, TranslationError>;
    fn try_to_codon(self, amino: B) -> Result<Seq<A>, TranslationError>;
}

pub struct CodonTable<A: Codec, B: Codec> {
    table: HashMap<Seq<A>, B>,
    inverse_table: HashMap<B, Seq<A>>,
}

impl<A: Codec, B: Codec> CodonTable<A, B> {
    pub fn from_map(table: HashMap<Seq<A>, B>) -> Self {
        let mut inverse_table = HashMap::new();
        for (codon, amino) in &table {
            inverse_table.insert(*amino, codon.clone());
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
        let mut inverse_table: HashMap<B, Seq<A>> = HashMap::new();

        for (codon, amino) in iter {
            table.insert(codon.clone(), amino);
            inverse_table.insert(amino, codon.clone());
        }

        CodonTable {
            table,
            inverse_table,
        }
    }
}

impl<A: Codec, B: Codec> PartialTranslationTable<A, B> for CodonTable<A, B> {
    fn try_to_amino(self, codon: &SeqSlice<A>) -> Result<B, TranslationError> {
        match self.table.get(codon) {
            Some(amino) => Ok(*amino),
            None => Err(TranslationError::InvalidCodon),
        }
    }

    fn try_to_codon(self, amino: B) -> Result<Seq<A>, TranslationError> {
        match self.inverse_table.get(&amino) {
            Some(codon) => Ok(codon.clone()),
            None => Err(TranslationError::AmbiguousCodon),
        }
    }
}

/*
#[cfg(test)]
mod tests {
}
*/
