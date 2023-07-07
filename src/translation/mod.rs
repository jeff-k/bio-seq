//! Translation tables for coding sequences

use crate::codec::{amino::Amino, dna::Dna, Codec};
use crate::kmer::Kmer;

struct Custom {
    table: Vec<Amino>,
}

impl Custom {
    pub fn from_vec(table: Vec<Amino>) -> Self {
        Custom { table }
    }
}

impl<A: Codec> TranslationTable<A> for Custom {
    fn to_amino(self, codon: Kmer<A, 3>) -> Amino {
        self.table[Into::<usize>::into(codon)]
    }

    fn to_codon(self, _amino: Amino) -> Option<Kmer<A, 3>> {
        unimplemented!()
    }
}

trait TranslationTable<A: Codec> {
    fn to_amino(self, codon: Kmer<A, 3>) -> Amino;
    fn to_codon(self, amino: Amino) -> Option<Kmer<A, 3>>;
}

struct Standard;

impl TranslationTable<Dna> for Standard {
    fn to_amino(self, codon: Kmer<Dna, 3>) -> Amino {
        Amino::unsafe_from_bits(Into::<usize>::into(codon) as u8)
    }

    fn to_codon(self, _amino: Amino) -> Option<Kmer<Dna, 3>> {
        None
    }
}

const STANDARD: Standard = Standard;

#[cfg(test)]
mod tests {
    use crate::prelude::*;
    use crate::translation::TranslationTable;
    use crate::translation::STANDARD;

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
        let codons: Vec<Kmer<Dna, 3>> = seq.chunks(3).map(|codon| codon.into()).collect();
        let aminos: Seq<Amino> = codons
            .iter()
            .map(|kmer| STANDARD.to_amino(*kmer))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("SSSSSS***"));
    }

    #[test]
    fn test_debruin_sequence() {
        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq
            .kmers()
            .map(|kmer| STANDARD.to_amino(kmer))
            .collect::<Seq<Amino>>();
        assert_eq!(aminos.len(), 64);
        assert_eq!(
            aminos,
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
        );
    }
}
