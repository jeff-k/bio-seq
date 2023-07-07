//! Translation tables for coding sequences

use crate::codec::{amino::Amino, dna::Dna, Codec};
use crate::kmer::Kmer;

struct CustomCode {
    table: Vec<Amino>,
}

impl CustomCode {
    fn new(_table: Vec<(Kmer<Dna, 3>, Amino)>) -> Self {
        unimplemented!()
    }
}

impl TranslationTable for CustomCode {
    fn translate(_codon: Kmer<Dna, 3>) -> Amino {
        unimplemented!()
    }

    fn back_translate(_amino: Amino) -> Option<Kmer<Dna, 3>> {
        unimplemented!()
    }
}

trait TranslationTable {
    fn translate(codon: Kmer<Dna, 3>) -> Amino;
    fn back_translate(amino: Amino) -> Option<Kmer<Dna, 3>>;
}

struct StandardCode;

impl TranslationTable for StandardCode {
    fn translate(codon: Kmer<Dna, 3>) -> Amino {
        Amino::unsafe_from_bits(Into::<usize>::into(codon) as u8)
    }

    fn back_translate(_amino: Amino) -> Option<Kmer<Dna, 3>> {
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::StandardCode;
    use crate::translation::TranslationTable;

    #[test]
    fn dna_to_amino() {
        let seq: Seq<Dna> = dna!("GCATGCGACGAATTCGGACACATAAAACTAATGAACCCACAAAGAAGCACAGTATGGTACTAA");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| StandardCode::translate(codon.into()))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("ACDEFGHIKLMNPQRSTVWY*"));
        assert_ne!(aminos, amino!("ACDEFGHIKLMNPQRSTVWY*A"));
        assert_ne!(aminos, amino!("CDEFGHIKLMNPQRSTVWY*A"));
    }

    #[test]
    fn alternate_codons() {
        let seq: Seq<Dna> = dna!("AGCTCGTCATCCTCTAGTTGATAATAG");
        let codons: Vec<Kmer<Dna, 3>> = seq.chunks(3).map(|codon| codon.into()).collect();
        let aminos: Seq<Amino> = codons
            .iter()
            .map(|kmer| StandardCode::translate(*kmer))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("SSSSSS***"));
    }

    #[test]
    fn test_debruin_sequence() {
        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq
            .kmers()
            .map(|kmer| StandardCode::translate(kmer))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos.len(), 64);
        assert_eq!(
            aminos,
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
        );
    }
}
