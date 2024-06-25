// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bit-packed and well-typed biological sequences
//!
//! A [`Seq`](seq::Seq) is a heap allocated [sequence](seq) of variable length that owns its data. A [`SeqSlice`](seq::SeqSlice) is a read-only window into a `Seq`.
//!
//! [`Kmer`](mod@kmer)s are short, fixed-length sequences. They generally implement `Copy` and are used for optimised algorithms on sequences. The default implementation uses a `usize` for storage.
//!
//! Binary encodings of genomic data types are implemented as "[`codec`]s." Custom codecs can be defined, and this crate has four built in:
//!   - [`codec::dna`]: 2-bit encoding of the four nucleotides
//!   - [`codec::text`]: 8-bit ASCII encoding of nucleotides, meant to be compatible with plaintext sequencing data formats
//!   - [`codec::iupac`]: 4-bit encoding of ambiguous nucleotide identities (the IUPAC ambiguity codes)
//!   - [`codec::amino`]: 6-bit encoding of amino acids
//!
//! Each of these encodings is designed to facilitate common bioinformatics tasks, such as minimising k-mers and implementing succinct datastructures. The [translation] module provides traits and methods for translating between nucleotide and amino acid sequences.
//!
//! Add `bio-seq` to `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! bio-seq = "0.12"
//! ```
//!
//! ```rust
//! use bio_seq::prelude::*;
//!
//! let seq = dna!("ATACGATCGATCGATCGATCCGT");
//!
//! // iterate over the 8-mers of the reverse complement
//! for kmer in seq.revcomp().kmers::<8>() {
//!     println!("{kmer}");
//! }
//!
//! // ACGGATCG
//! // CGGATCGA
//! // GGATCGAT
//! // GATCGATC
//! // ATCGATCG
//! // ...
//! ```
//!
#![warn(clippy::pedantic)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::module_name_repetitions)]

use bitvec::prelude::*;

type Order = Lsb0;
type Bs = BitSlice<usize, Order>;
type Bv = BitVec<usize, Order>;
type Ba = BitArray<usize, Order>;

#[macro_use]
pub mod codec;
pub mod error;
pub mod kmer;
pub mod seq;

//#[macro_use]
pub use bio_seq_derive::dna;

pub use bitvec::bitarr;
pub use bitvec::prelude::Lsb0;

#[cfg(feature = "translation")]
pub mod translation;

pub mod prelude {
    pub use crate::codec::amino::Amino;
    pub use crate::codec::dna::Dna;
    pub use crate::codec::iupac::Iupac;
    pub use crate::codec::{Codec, Complement};

    pub use crate::kmer::Kmer;
    pub use crate::seq::{ReverseComplement, Seq, SeqArray, SeqSlice};
    pub use crate::{amino, iupac, kmer};

    #[cfg(feature = "translation")]
    pub use crate::translation;

    pub use core::str::FromStr;

    pub use crate::error::ParseBioError;

    pub use crate::dna;

    pub use crate::bitarr as __bio_seq_bitarr;
    pub use crate::Lsb0 as __bio_seq_Lsb0;
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::Dna::{A, C, G, T};
    use crate::prelude::*;

    #[test]
    fn alt_repr() {
        assert_eq!(iupac!("-").nth(0), Iupac::X);
    }

    /*
    #[test]
    fn into_usize() {
        let a: usize = dna!("ACGT").into();
        assert_eq!(a, 0b11_10_01_00);

        let b: usize = dna!("CGCG").into();
        assert_eq!(b, 0b10_01_10_01);

        let c: usize = Seq::from(&vec![T, T]).into();
        assert_eq!(c, 0b11_11);

        let d: usize = Seq::<Dna>::from_str("TCA").unwrap().into();
        assert_eq!(d, 0b00_01_11);

        let e: usize = Seq::<Dna>::from_str("TGA").unwrap().into();
        assert_eq!(e, 0b00_10_11);

        let f: usize = Seq::from(&vec![C, G, T, A, C, G, A, T]).into();
        assert_eq!(f, 0b11_00_10_01_00_11_10_01);

        let g: usize = Seq::from(&vec![A]).into();
        assert_eq!(g, 0b00);
    }
    */
    #[test]
    fn test_display_aminos() {
        let a: Seq<Amino> = Seq::from_str("DCMNLKG*HI").unwrap();
        assert_eq!(format!("{a}"), "DCMNLKG*HI");
    }
    #[test]
    fn test_display_dna() {
        let seq = Seq::from(&vec![A, C, G, T, T, A, T, C]);
        assert_eq!(format!("{}", &seq), "ACGTTATC");
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

    #[test]
    fn colexicographic_order() {
        for (i, e) in ["AA", "CA", "GA", "TA", "AC", "CC", "GC", "TC"]
            .iter()
            .enumerate()
        {
            assert_eq!(format!("{}", Kmer::<Dna, 2>::from(i)), format!("{}", e));
            assert_eq!(Kmer::<Dna, 2>::from(i), *e);
        }
    }
}
