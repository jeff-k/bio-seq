// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

/*! # `bioseq`

Bit packed biological sequences

!*/

pub mod alphabet;

#[macro_use]
pub mod seq;

pub mod kmer;

pub use alphabet::dna::Dna;
pub use seq::Seq;
pub use std::str::FromStr;

#[cfg(test)]
mod tests {
    use crate::alphabet::dna::Dna;
    use crate::alphabet::dna::Dna::{A, C, G, T};
    use crate::seq::Seq;
    use std::str::FromStr;

    #[test]
    fn make_from_vector() {
        assert_eq!(Seq::from_vec(vec![A, C, G, T]).to_usize(), 0b11011000);
        assert_eq!(Seq::from_vec(vec![C, G, C, G]).to_usize(), 0b01100110);
        assert_eq!(Seq::from_vec(vec![T, T]).to_usize(), 0b1111);
        assert_eq!(
            Seq::from_vec(vec![C, G, T, A, C, G, A, T]).to_usize(),
            0b1100011000110110
        );
        assert_eq!(Seq::from_vec(vec![A,]).to_usize(), 0b00);
    }

    #[test]
    fn test_display_dna() {
        let seq = Seq::from_vec(vec![A, C, G, T, T, A, T, C]);
        assert_eq!(format!("{}", seq), "ACGTTATC");
        assert_eq!(format!("{}", dna!("ACGT")), "ACGT");
    }

    #[test]
    fn iterate_bases() {
        let seq = dna!("ACGTACGT");
        assert_eq!(
            seq.into_iter().collect::<Vec<Dna>>(),
            vec![A, C, G, T, A, C, G, T]
        );
    }

    #[test]
    fn iterate_kmers() {
        let seq = dna!("ACGTAAGGGG");
        for (kmer, answer) in seq
            .kmers::<4>()
            .zip(["ACGT", "CGTA", "GTAA", "TAAG", "AAGG", "AGGG", "GGGG"])
        {
            assert_eq!(format!("{}", kmer), answer);
        }
    }

    #[test]
    fn iterate_kmer8() {
        let seq = dna!("AAAACCCCGGGG");
        for (kmer, answer) in seq
            .kmers::<8>()
            .zip(["AAAACCCC", "AAACCCCG", "AACCCCGG", "ACCCCGGG", "CCCCGGGG"])
        {
            assert_eq!(format!("{}", kmer), answer);
        }
    }

    #[test]
    fn iterate_kmer4() {
        let seq = dna!("AAAACCCCGGGGTTTT");
        for (kmer, answer) in seq.kmers::<4>().zip([
            "AAAA", "AAAC", "AACC", "ACCC", "CCCC", "CCCG", "CCGG", "CGGG", "GGGG", "GGGT", "GGTT",
            "GTTT", "TTTT",
        ]) {
            assert_eq!(format!("{}", kmer), answer);
        }
    }
}
