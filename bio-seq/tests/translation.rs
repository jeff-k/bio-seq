use bio_seq::prelude::*;
use bio_seq::translation::{
    CodonTable, PartialTranslationTable, STANDARD, TranslationError, TranslationTable,
};

fn message(error: impl std::error::Error) -> String {
    error.to_string()
}

#[test]
fn errors_display_the_offending_codon_or_amino_acid() {
    let table = CodonTable::from_map([(Seq::<Dna>::from(dna!("ATG")), Amino::M)]);

    assert_eq!(
        message(STANDARD.to_codon(Amino::L).unwrap_err()),
        "Multiple codon sequences: L"
    );
    assert_eq!(
        message(STANDARD.try_to_amino(iupac!("NNN")).unwrap_err()),
        "Ambiguous translations for codon: NNN"
    );
    assert_eq!(
        message(STANDARD.try_to_amino(iupac!("A-G")).unwrap_err()),
        "Invalid codon sequence: A-G"
    );
    assert_eq!(
        message(table.try_to_codon(Amino::K).unwrap_err()),
        "Invalid amino acid character: K"
    );
    assert_eq!(
        message(table.try_to_codon(Amino::X).unwrap_err()),
        "Invalid amino acid character: *"
    );
}

#[test]
fn codon_tables_reject_what_they_do_not_contain() {
    let table = CodonTable::from_map([
        (Seq::<Dna>::from(dna!("ATG")), Amino::M),
        (dna!("TAA").into(), Amino::X),
        (dna!("TAG").into(), Amino::X),
    ]);

    assert_eq!(table.try_to_amino(dna!("ATG")), Ok(Amino::M));
    for codon in [dna!("CCC"), dna!("AT"), dna!("ATGA")] {
        assert_eq!(
            table.try_to_amino(codon),
            Err(TranslationError::InvalidCodon(codon.into())),
            "{codon}"
        );
    }

    assert_eq!(table.try_to_codon(Amino::M), Ok(dna!("ATG").into()));
    assert_eq!(
        table.try_to_codon(Amino::X),
        Err(TranslationError::AmbiguousCodon(Amino::X))
    );
    assert_eq!(
        table.try_to_codon(Amino::K),
        Err(TranslationError::InvalidAmino(Amino::K))
    );
}

#[test]
fn only_methionine_and_tryptophan_have_a_single_standard_codon() {
    for amino in Amino::items() {
        let expected = match amino {
            Amino::M => Ok(dna!("ATG").to_owned()),
            Amino::W => Ok(dna!("TGG").to_owned()),
            _ => Err(TranslationError::AmbiguousCodon(amino)),
        };
        assert_eq!(STANDARD.to_codon(amino), expected, "{amino}");
    }
}
