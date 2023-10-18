use std::str::FromStr;

use crate::codec::Codec;
use crate::codec::{amino::Amino, dna::Dna, iupac::Iupac};
use crate::prelude::{Seq, SeqSlice};
use crate::translation::{PartialTranslationTable, TranslationError, TranslationTable};

struct Standard;

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
    fn try_to_amino(self, _codon: &SeqSlice<Iupac>) -> Result<Amino, TranslationError> {
        unimplemented!();
        //        Err(TranslationError::InvalidCodon)
    }

    fn try_to_codon(self, amino: Amino) -> Result<Seq<Iupac>, TranslationError> {
        match amino {
            Amino::A => Ok(iupac!("GCN")),
            Amino::C => Ok(iupac!("UGY")),
            Amino::D => Ok(iupac!("GAY")),
            Amino::E => Ok(iupac!("GAR")),
            Amino::F => Ok(iupac!("UUY")),
            Amino::G => Ok(iupac!("GGN")),
            Amino::H => Ok(iupac!("CAY")),
            Amino::I => Ok(iupac!("AUH")),
            Amino::K => Ok(iupac!("AAR")),
            Amino::L => Err(TranslationError::AmbiguousCodon),
            Amino::M => Ok(iupac!("AUG")),
            Amino::N => Ok(iupac!("AAY")),
            Amino::P => Ok(iupac!("CCN")),
            Amino::Q => Ok(iupac!("CAR")),
            Amino::R => Err(TranslationError::AmbiguousCodon),
            Amino::S => Err(TranslationError::AmbiguousCodon),
            Amino::T => Ok(iupac!("ACN")),
            Amino::V => Ok(iupac!("GUN")),
            Amino::W => Ok(iupac!("UGG")),
            Amino::Y => Ok(iupac!("UAY")),
            Amino::X => Err(TranslationError::AmbiguousCodon),
        }
    }
}

#[allow(dead_code)]
const STANDARD: Standard = Standard;

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::standard::STANDARD;
    use crate::translation::TranslationTable;

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
}
