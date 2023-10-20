//use std::lazy::OnceCell;
use once_cell::sync::OnceCell;
use std::collections::HashMap;
use std::str::FromStr;

use crate::codec::Codec;
use crate::codec::{amino::Amino, dna::Dna, iupac::Iupac};
use crate::prelude::{Seq, SeqSlice};
use crate::translation::{PartialTranslationTable, TranslationError, TranslationTable};

struct Standard;

static AMINO_TO_IUPAC: OnceCell<HashMap<Amino, Option<Seq<Iupac>>>> = OnceCell::new();

static IUPAC_TO_AMINO: OnceCell<[(Seq<Iupac>, Amino); 29]> = OnceCell::new();

fn initialise_iupac_to_amino() -> [(Seq<Iupac>, Amino); 29] {
    [
        (iupac!("GCN"), Amino::A),
        (iupac!("TGY"), Amino::C),
        (iupac!("GAY"), Amino::D),
        (iupac!("GAR"), Amino::E),
        (iupac!("TTY"), Amino::F),
        (iupac!("GGN"), Amino::G),
        (iupac!("CAY"), Amino::H),
        (iupac!("ATH"), Amino::I),
        (iupac!("AAR"), Amino::K),
        (iupac!("CTN"), Amino::L),
        (iupac!("TTR"), Amino::L),
        (iupac!("CTY"), Amino::L),
        (iupac!("YTR"), Amino::L),
        (iupac!("ATG"), Amino::M),
        (iupac!("AAY"), Amino::N),
        (iupac!("CCN"), Amino::P),
        (iupac!("CAR"), Amino::Q),
        (iupac!("CGN"), Amino::R),
        (iupac!("AGR"), Amino::R),
        (iupac!("CGY"), Amino::R),
        (iupac!("MGR"), Amino::R),
        (iupac!("TCN"), Amino::S),
        (iupac!("AGY"), Amino::S),
        (iupac!("ACN"), Amino::T),
        (iupac!("GTN"), Amino::V),
        (iupac!("TGG"), Amino::W),
        (iupac!("TAY"), Amino::Y),
        (iupac!("TAR"), Amino::X),
        (iupac!("TRA"), Amino::X),
    ]
}

fn initialise_amino_to_iupac() -> HashMap<Amino, Option<Seq<Iupac>>> {
    let mut amino_to_iupac = HashMap::new();
    for (iupac_set, amino) in IUPAC_TO_AMINO.get_or_init(initialise_iupac_to_amino) {
        if amino_to_iupac.contains_key(amino) {
            amino_to_iupac.insert(*amino, None);
        } else {
            amino_to_iupac.insert(*amino, Some(iupac_set.clone()));
        }
    }

    amino_to_iupac
}

impl TranslationTable<Dna, Amino> for Standard {
    fn to_amino(self, codon: &SeqSlice<Dna>) -> Amino {
        Amino::unsafe_from_bits(Into::<u8>::into(codon))
    }

    /// There are no unambiguous translations from amino acid to DNA codon
    fn to_codon(self, _amino: Amino) -> Result<Seq<Dna>, TranslationError> {
        Err(TranslationError::AmbiguousCodon)
    }
}

impl PartialTranslationTable<Iupac, Amino> for Standard {
    fn try_to_amino(self, codon: &SeqSlice<Iupac>) -> Result<Amino, TranslationError> {
        for (iupac_set, amino) in IUPAC_TO_AMINO.get_or_init(initialise_iupac_to_amino) {
            if iupac_set.contains(codon) {
                return Ok(*amino);
            }
        }

        Err(TranslationError::AmbiguousCodon)
    }

    fn try_to_codon(self, amino: Amino) -> Result<Seq<Iupac>, TranslationError> {
        let amino_to_iupac = AMINO_TO_IUPAC.get_or_init(initialise_amino_to_iupac);

        match amino_to_iupac.get(&amino) {
            Some(None) => Err(TranslationError::AmbiguousCodon),
            Some(Some(codon)) => Ok(codon.clone()),
            None => Err(TranslationError::AmbiguousCodon),
        }
    }
}

#[allow(dead_code)]
const STANDARD: Standard = Standard;

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::standard::STANDARD;
    use crate::translation::{PartialTranslationTable, TranslationTable};

    #[test]
    fn dna_to_amino() {
        let seq: Seq<Dna> = dna!("GCATGCGACGAATTCGGACACATAAAACTAATGAACCCACAAAGAAGCACAGTATGGTACTAA");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.to_amino(codon.into()))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("ACDEFGHIKLMNPQRSTVWY*"));
        assert_ne!(aminos, amino!("ACDEFGHIKLMNPQRSTVWY*A"));
        assert_ne!(aminos, amino!("CDEFGHIKLMNPQRSTVWY*A"));
    }

    #[test]
    fn alternate_codons() {
        let seq: Seq<Dna> = dna!("AGCTCGTCATCCTCTAGTTGATAATAG");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.to_amino(codon.into()))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("SSSSSS***"));
    }

    #[test]
    fn test_debruin_sequence() {
        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq
            .windows(3)
            .map(|codon| STANDARD.to_amino(&codon))
            .collect::<Seq<Amino>>();
        assert_eq!(seq.len() - 2, aminos.len());
        assert_eq!(
            aminos,
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
        );
    }

    #[test]
    fn test_iupac_codons() {
        let seq: Seq<Iupac> = iupac!("GCTTGYGCAGCN");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.try_to_amino(&codon).unwrap())
            .collect::<Seq<Amino>>();
        println!("{}", aminos);
        assert_eq!(aminos, amino!("ACAA"));
    }
}
