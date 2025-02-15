//! Standard amino acid translation table
use std::collections::HashMap;
use std::sync::OnceLock;

use crate::codec::Codec;
use crate::codec::{amino::Amino, dna::Dna, iupac::Iupac};
use crate::iupac;
use crate::prelude::{Seq, SeqArray, SeqSlice};
use crate::prelude::{__bio_seq_Lsb0, __bio_seq_bitarr, __bio_seq_count_words};
use crate::translation::{PartialTranslationTable, TranslationError, TranslationTable};

pub struct Standard;

static AMINO_TO_IUPAC: OnceLock<HashMap<Amino, Option<Seq<Iupac>>>> = OnceLock::new();

/// We can also translate from degenerate nucleotide symbols to amino acids where the coding isn't ambiguous
///
/// ```
/// use bio_seq::prelude::*;
/// use bio_seq::translation::PartialTranslationTable;
/// use bio_seq::translation::STANDARD;
///
/// let seq: Seq<Iupac> = iupac!("AAYATHTTYTGYGTNTGGGGNGGNGTNTTYGTNTGYGCNGGNGCNCCNGCNCCNCCNGTNTAYACNTAYATGGARGGNGAYACNGAYATHCARGCNCAYACNCCNCAYATHAARGARAAYACNCARAAR").into();
///     let aminos: Seq<Amino> = seq
///     .chunks(3)
///     .map(|codon| STANDARD.try_to_amino(&codon).unwrap())
///     .collect::<Seq<Amino>>();
/// assert_eq!(
///     aminos,
///     Seq::<Amino>::try_from("NIFCVWGGVFVCAGAPAPPVYTYMEGDTDIQAHTPHIKENTQK").unwrap()
/// );
/// ```
static IUPAC_TO_AMINO: OnceLock<[(Seq<Iupac>, Amino); 29]> = OnceLock::new();

fn initialise_iupac_to_amino() -> [(Seq<Iupac>, Amino); 29] {
    [
        (iupac!("GCN").into(), Amino::A),
        (iupac!("TGY").into(), Amino::C),
        (iupac!("GAY").into(), Amino::D),
        (iupac!("GAR").into(), Amino::E),
        (iupac!("TTY").into(), Amino::F),
        (iupac!("GGN").into(), Amino::G),
        (iupac!("CAY").into(), Amino::H),
        (iupac!("ATH").into(), Amino::I),
        (iupac!("AAR").into(), Amino::K),
        (iupac!("CTN").into(), Amino::L),
        (iupac!("TTR").into(), Amino::L),
        (iupac!("CTY").into(), Amino::L),
        (iupac!("YTR").into(), Amino::L),
        (iupac!("ATG").into(), Amino::M),
        (iupac!("AAY").into(), Amino::N),
        (iupac!("CCN").into(), Amino::P),
        (iupac!("CAR").into(), Amino::Q),
        (iupac!("CGN").into(), Amino::R),
        (iupac!("AGR").into(), Amino::R),
        (iupac!("CGY").into(), Amino::R),
        (iupac!("MGR").into(), Amino::R),
        (iupac!("TCN").into(), Amino::S),
        (iupac!("AGY").into(), Amino::S),
        (iupac!("ACN").into(), Amino::T),
        (iupac!("GTN").into(), Amino::V),
        (iupac!("TGG").into(), Amino::W),
        (iupac!("TAY").into(), Amino::Y),
        (iupac!("TAR").into(), Amino::X),
        (iupac!("TRA").into(), Amino::X),
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
        assert_eq!(
            aminos,
            Seq::<Amino>::try_from("ACDEFGHIKLMNPQRSTVWY*").unwrap()
        );
        assert_ne!(
            aminos,
            Seq::<Amino>::try_from("ACDEFGHIKLMNPQRSTVWY*A").unwrap()
        );
        assert_ne!(
            aminos,
            Seq::<Amino>::try_from("CDEFGHIKLMNPQRSTVWY*A").unwrap()
        );
    }

    #[test]
    fn alternate_codons() {
        let seq: Seq<Dna> = dna!("AGCTCGTCATCCTCTAGTTGATAATAG").into();
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.to_amino(codon.into()))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, Seq::<Amino>::try_from("SSSSSS***").unwrap());
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
            Seq::<Amino>::try_from(
                "NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_iupac_codons() {
        let seq = iupac!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq
            .windows(3)
            .map(|codon| STANDARD.try_to_amino(&codon).unwrap())
            .collect::<Seq<Amino>>();
        assert_eq!(
            aminos,
            Seq::<Amino>::try_from(
                "NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_iupac_aminos() {
        let seq: Seq<Amino> = Seq::<Amino>::try_from(
            "NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK",
        )
        .unwrap();
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
        assert_eq!(
            ambs,
            Seq::<Amino>::try_from("LSRSLRLSRLL*SL*RRSS*R").unwrap()
        );
    }

    #[test]
    fn test_iupac_ambiguous_codons() {
        let seq: Seq<Iupac> = iupac!("AAYATHTTYTGYGTNTGGGGNGGNGTNTTYGTNTGYGCNGGNGCNCCNGCNCCNCCNGTNTAYACNTAYATGGARGGNGAYACNGAYATHCARGCNCAYACNCCNCAYATHAARGARAAYACNCARAAR").into();
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| STANDARD.try_to_amino(&codon).unwrap())
            .collect::<Seq<Amino>>();
        assert_eq!(
            aminos,
            Seq::<Amino>::try_from("NIFCVWGGVFVCAGAPAPPVYTYMEGDTDIQAHTPHIKENTQK").unwrap()
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
