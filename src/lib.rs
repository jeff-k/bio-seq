// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bit-packed and well-typed biological sequences
//!
//! - [Seq] heap allocated sequences of variable length
//! - [Kmer] short fixed length sequences
//! - [Codec] coder/decoder implementations
//!
//! This crate is designed to facilitate common bioinformatics tasks,
//! incuding amino acid translation, k-mer minimisation and hashing, and
//! nucleotide sequence manipulation.
//!
//! ```rust
//! use bio_seq::prelude::*;
//!
//! let seq = dna!("ATACGATCGATCGATCGATCCGT");
//!
//! // iterate over the 8-mers of the reverse complement
//! for kmer in seq.revcomp().kmers::<8>() {
//!     println!("{}", kmer);
//! }
//! ```
//!
//! Custom encodings are supported with the help of the `bio-seq-derive`
//! crate.
//!
//! ```rust
//! use bio_seq::prelude::*;
//!
//! let seq = iupac!("AGCTNNCAGTCGACGTATGTA");
//! let pattern = Seq::<Iupac>::from_str("AYG").unwrap();
//!
//! for slice in seq.windows(pattern.len()) {
//!     if pattern.contains(slice) {
//!         println!("{} matches pattern", slice);
//!     }
//! }
//! ```

#[macro_use]
pub mod codec;
pub mod kmer;
pub mod seq;
pub mod translation;

pub use crate::kmer::Kmer;
pub use crate::seq::{Seq, SeqSlice};

pub use core::str::FromStr;

#[cfg(test)]
mod tests {
    use crate::codec::amino::Amino;
    use crate::codec::dna::Dna;
    use crate::codec::dna::Dna::{A, C, G, T};
    use crate::codec::iupac::Iupac;
    use crate::seq::Seq;
    use core::str::FromStr;

    #[test]
    fn alt_repr() {
        assert_eq!(iupac!("-").nth(0), Iupac::X);
    }

    #[test]
    fn into_usize() {
        let a: usize = dna!("ACGT").into();
        assert_eq!(a, 0b11_10_01_00);

        let b: usize = dna!("CGCG").into();
        assert_eq!(b, 0b10_01_10_01);

        let c: usize = Seq::from_vec(vec![T, T]).into();
        assert_eq!(c, 0b11_11);

        let d: usize = Seq::<Dna>::from_str("TCA").unwrap().into();
        assert_eq!(d, 0b00_01_11);

        let e: usize = Seq::<Dna>::from_str("TGA").unwrap().into();
        assert_eq!(e, 0b00_10_11);

        let f: usize = Seq::from_vec(vec![C, G, T, A, C, G, A, T]).into();
        assert_eq!(f, 0b11_00_10_01_00_11_10_01);

        let g: usize = Seq::from_vec(vec![A]).into();
        assert_eq!(g, 0b00);
    }

    #[test]
    fn test_display_aminos() {
        let a: Seq<Amino> = Seq::from_str("DCMNLKG*HI").unwrap();
        assert_eq!(format!("{}", a), "DCMNLKG*HI");
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
    fn from_string() {
        let seq = Seq::<Dna>::from_str("ACGTACGT").unwrap();
        assert_eq!(
            seq.into_iter().collect::<Vec<Dna>>(),
            vec![A, C, G, T, A, C, G, T]
        );
    }
    #[test]
    fn rev_seq() {
        let seq = dna!("ACGTACGT");
        assert_eq!(
            seq.rev().collect::<Vec<Dna>>(),
            vec![T, G, C, A, T, G, C, A]
        );

        assert_eq!(
            iupac!("GN-").rev().collect::<Vec<Iupac>>(),
            vec![Iupac::X, Iupac::N, Iupac::G]
        );

        assert_eq!(
            amino!("DCMNLKGHI").rev().collect::<Vec<Amino>>(),
            vec![
                Amino::I,
                Amino::H,
                Amino::G,
                Amino::K,
                Amino::L,
                Amino::N,
                Amino::M,
                Amino::C,
                Amino::D
            ]
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

    #[test]
    fn iupac_bitwise_ops() {
        assert_eq!(iupac!("AS-GYTNA") | iupac!("ANTGCAT-"), iupac!("ANTGYWNA"));
        assert_eq!(iupac!("ACGTSWKM") & iupac!("WKMSTNNA"), iupac!("A----WKA"));
    }

    #[test]
    fn nth_chars() {
        assert_eq!(iupac!("ACGTRYSWKMBDHVN-").nth(0), Iupac::A);
        assert_ne!(iupac!("ACGTRYSWKMBDHVN-").nth(0), Iupac::C);
        assert_eq!(iupac!("ACGTRYSWKMBDHVN-").nth(15), Iupac::X);
        assert_eq!(iupac!("ACGTRYSWKMBDHVN-").nth(3), Iupac::from(Dna::T));
        assert_ne!(iupac!("ACGTRYSWKMBDHVN-").nth(3), Iupac::from(Dna::G));

        assert_eq!(amino!("DCMNLKGHI").nth(1), Amino::C);
        assert_ne!(amino!("DCMNLKGHI").nth(7), Amino::I);
    }
}
