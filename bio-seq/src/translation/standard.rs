//! Standard amino acid translation table

use core::str::FromStr;
use std::collections::HashMap;
use std::sync::OnceLock;

use crate::codec::Codec;
use crate::codec::{amino::Amino, dna::Dna, iupac::Iupac};
use crate::prelude::{Seq, SeqSlice};
use crate::translation::{PartialTranslationTable, TranslationError, TranslationTable};

pub struct Standard;

static AMINO_TO_IUPAC: OnceLock<HashMap<Amino, Option<Seq<Iupac>>>> = OnceLock::new();

static IUPAC_TO_AMINO: OnceLock<[(Seq<Iupac>, Amino); 29]> = OnceLock::new();

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
    fn to_amino(&self, codon: &SeqSlice<Dna>) -> Amino {
        assert!(codon.len() == 3, "Invalid codon length {}", codon.len());
        // It should be possible to assert that this is safe at compile time
        Amino::unsafe_from_bits(Into::<u8>::into(codon))
    }

    /// There are no unambiguous translations from amino acid to DNA codon except M and W
    fn to_codon(&self, amino: Amino) -> Result<Seq<Dna>, TranslationError<Dna, Amino>> {
        Err(TranslationError::AmbiguousCodon(amino))
    }
}

impl PartialTranslationTable<Iupac, Amino> for Standard {
    fn try_to_amino(
        &self,
        codon: &SeqSlice<Iupac>,
    ) -> Result<Amino, TranslationError<Iupac, Amino>> {
        if codon.len() != 3 {
            return Err(TranslationError::InvalidCodon(codon.into()));
        }
        for (iupac_set, amino) in IUPAC_TO_AMINO.get_or_init(initialise_iupac_to_amino) {
            if iupac_set.contains(codon) {
                return Ok(*amino);
            }
        }

        Err(TranslationError::AmbiguousTranslation(codon.into()))
    }

    fn try_to_codon(&self, amino: Amino) -> Result<Seq<Iupac>, TranslationError<Iupac, Amino>> {
        let amino_to_iupac = AMINO_TO_IUPAC.get_or_init(initialise_amino_to_iupac);

        match amino_to_iupac.get(&amino) {
            Some(Some(codon)) => Ok(codon.clone()),
            None | Some(None) => Err(TranslationError::AmbiguousCodon(amino)),
        }
    }
}

pub const STANDARD: Standard = Standard;

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::STANDARD;
    use crate::translation::{PartialTranslationTable, TranslationError, TranslationTable};

    #[test]
    fn dna_to_amino() {
        let seq: Seq<Dna> =
            dna!("GCATGCGACGAATTCGGACACATAAAACTAATGAACCCACAAAGAAGCACAGTATGGTACTAA").into();
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
        let seq: Seq<Dna> = dna!("AGCTCGTCATCCTCTAGTTGATAATAG").into();
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.to_amino(codon.into()))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("SSSSSS***"));
    }

    #[test]
    fn test_debruin_sequence() {
        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA").into();
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
        let seq: Seq<Iupac> =
            iupac!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq
            .windows(3)
            .map(|codon| STANDARD.try_to_amino(&codon).unwrap())
            .collect::<Seq<Amino>>();
        assert_eq!(
            aminos,
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
        );
    }

    #[test]
    fn test_iupac_aminos() {
        let seq: Seq<Amino> =
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK");
        let mut iupacs: Seq<Iupac> = Seq::new();
        let mut ambs: Seq<Amino> = Seq::new();

        for amino in &seq {
            match STANDARD.try_to_codon(amino) {
                Ok(codon) => iupacs.extend(&codon),
                Err(TranslationError::AmbiguousCodon(amino)) => ambs.push(amino),
                _ => panic!(),
            }
        }
        assert_eq!(iupacs, iupac!("AAYATHTTYTGYGTNTGGGGNGGNGTNTTYGTNTGYGCNGGNGCNCCNGCNCCNCCNGTNTAYACNTAYATGGARGGNGAYACNGAYATHCARGCNCAYACNCCNCAYATHAARGARAAYACNCARAAR"));
        assert_eq!(ambs, amino!("LSRSLRLSRLL*SL*RRSS*R"));
    }

    #[test]
    fn test_iupac_ambiguous_codons() {
        let seq: Seq<Iupac> = iupac!("AAYATHTTYTGYGTNTGGGGNGGNGTNTTYGTNTGYGCNGGNGCNCCNGCNCCNCCNGTNTAYACNTAYATGGARGGNGAYACNGAYATHCARGCNCAYACNCCNCAYATHAARGARAAYACNCARAAR");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.try_to_amino(&codon).unwrap())
            .collect::<Seq<Amino>>();
        assert_eq!(
            aminos,
            amino!("NIFCVWGGVFVCAGAPAPPVYTYMEGDTDIQAHTPHIKENTQK")
        );
    }

    #[test]
    fn iupac_to_amino_errors() {
        assert_eq!(
            STANDARD.try_to_amino(&iupac!("TN")),
            Err(TranslationError::InvalidCodon(
                Seq::<Iupac>::try_from("TN").unwrap()
            ))
        );

        assert_eq!(
            STANDARD.try_to_amino(&iupac!("NYTN")),
            Err(TranslationError::InvalidCodon("NYTN".try_into().unwrap()))
        );
        assert_eq!(
            STANDARD.try_to_amino(&iupac!("YTN")),
            Err(TranslationError::AmbiguousTranslation(
                "YTN".try_into().unwrap()
            ))
        );
    }

    #[test]
    fn ambiguous_amino_to_iupac() {
        assert_eq!(
            STANDARD.try_to_codon(Amino::L),
            Err(TranslationError::AmbiguousCodon(Amino::L))
        );
        assert_eq!(
            STANDARD.try_to_codon(Amino::S),
            Err(TranslationError::AmbiguousCodon(Amino::S))
        );
        assert_eq!(
            STANDARD.try_to_codon(Amino::R),
            Err(TranslationError::AmbiguousCodon(Amino::R))
        );
        assert_eq!(
            STANDARD.try_to_codon(Amino::X),
            Err(TranslationError::AmbiguousCodon(Amino::X))
        );

        assert_ne!(
            STANDARD.try_to_codon(Amino::A),
            Err(TranslationError::AmbiguousCodon(Amino::A))
        );
        assert_ne!(
            STANDARD.try_to_codon(Amino::N),
            Err(TranslationError::AmbiguousCodon(Amino::N))
        );
        assert_ne!(
            STANDARD.try_to_codon(Amino::W),
            Err(TranslationError::AmbiguousCodon(Amino::W))
        );
        assert_ne!(
            STANDARD.try_to_codon(Amino::M),
            Err(TranslationError::AmbiguousCodon(Amino::M))
        );
    }
}
