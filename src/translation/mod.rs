//! Translation tables for coding sequences

use crate::codec::{amino::Amino, dna::Dna, Codec};
use crate::kmer::Kmer;
use crate::seq::SeqSlice;

impl From<Kmer<Dna, 3>> for Amino {
    fn from(k: Kmer<Dna, 3>) -> Amino {
        let x: usize = k.into();
        Amino::unsafe_from_bits(x as u8)
    }
}

impl From<&Kmer<Dna, 3>> for Amino {
    fn from(k: &Kmer<Dna, 3>) -> Amino {
        let x: usize = k.into();
        Amino::unsafe_from_bits(x as u8)
    }
}

impl From<&SeqSlice<Dna>> for Amino {
    fn from(seq: &SeqSlice<Dna>) -> Amino {
        assert_eq!(seq.len(), 3);
        Amino::unsafe_from_bits(seq.into())
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn dna_to_amino() {
        let seq: Seq<Dna> = dna!("GCATGCGACGAATTCGGACACATAAAACTAATGAACCCACAAAGAAGCACAGTATGGTACTAA");
        let aminos: Seq<Amino> = seq
            .chunks(3)
            .map(|codon| codon.into())
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("ACDEFGHIKLMNPQRSTVWY*"));
    }

    #[test]
    fn alternate_codons() {
        let seq: Seq<Dna> = dna!("AGCTCGTCATCCTCTAGTTGATAATAG");
        let codons: Vec<Kmer<Dna, 3>> = seq.chunks(3).map(|codon| codon.into()).collect();
        let aminos: Seq<Amino> = codons
            .iter()
            .map(|kmer| kmer.into())
            .collect::<Seq<Amino>>();
        assert_eq!(aminos, amino!("SSSSSS***"));
    }

    #[test]
    fn test_debruin_sequence() {
        let seq: Seq<Dna> =
            dna!("AATTTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATGAGGACGATCAGCACCATAAGAACAAA");
        let aminos: Seq<Amino> = seq.kmers().map(|kmer| kmer.into()).collect::<Seq<Amino>>();
        assert_eq!(aminos.len(), 64);
        assert_eq!(
            aminos,
            amino!("NIFLCVWGGVFSRVSLCARGALSPRAPPLL*SVYTLYM*ERGDTRDISQSAHTPHI*KRENTQK")
        );
    }
}
