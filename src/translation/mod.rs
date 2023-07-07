//! Translation tables for coding sequences

use crate::codec::{amino::Amino, dna::Dna, Codec};
use crate::kmer::Kmer;

struct Custom {
    table: Vec<Amino>,
}

impl Custom {
    fn from_vec(table: Vec<Amino>) -> Self {
        Custom { table }
    }
}

impl TranslationTable for Custom {
    fn to_amino(self, _codon: Kmer<Dna, 3>) -> Amino {
        unimplemented!()
    }

    fn from_amino(self, _amino: Amino) -> Option<Kmer<Dna, 3>> {
        unimplemented!()
    }
}

trait TranslationTable {
    fn to_amino(self, codon: Kmer<Dna, 3>) -> Amino;
    fn from_amino(self, amino: Amino) -> Option<Kmer<Dna, 3>>;
}

struct Standard;

impl TranslationTable for Standard {
    fn to_amino(self, codon: Kmer<Dna, 3>) -> Amino {
        Amino::unsafe_from_bits(Into::<usize>::into(codon) as u8)
    }

    fn from_amino(self, _amino: Amino) -> Option<Kmer<Dna, 3>> {
        None
    }
}

const standard: Standard = Standard;

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::standard;
    use crate::translation::TranslationTable;

    #[test]
    fn dna_to_amino() {
        let seq: Seq<Dna> = dna!("GCATGCGACGAATTCGGACACATAAAACTAATGAACCCACAAAGAAGCACAGTATGGTACTAA");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| standard.to_amino(codon.into()))
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
            .map(|kmer| standard.to_amino(*kmer))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("SSSSSS***"));
    }

    #[test]
    fn test_debruin_sequence() {
        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq
            .kmers()
            .map(|kmer| standard.to_amino(kmer))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos.len(), 64);
        assert_eq!(
            aminos,
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
        );
    }
}
