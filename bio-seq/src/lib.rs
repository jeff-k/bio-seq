// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Bit-packed and well-typed biological sequences
//!
//! The strength of rust is that we can safely separate the science (well-typed) and the engineering (bit-packed) of bioinformatics. An incremental benchmark improvement in the reverse complement algorithm should benefit the user of a succinct datastructure without anyone unwillingly learning about endianess.
//!
//! Contributions are very welcome. There's lots of low hanging fruit for optimisation and ideally we should only have to write them once!
//!
//! ## Sequences
//!
//! A [`Seq`](seq::Seq) is a heap allocated [sequence](seq) of symbols that owns its data. A [`SeqSlice`](seq::SeqSlice) is a read-only window into a `Seq`. Static [`fixed length sequences`](seq::SeqArray) can be declared with the [`dna!`](macro@dna) and [`iupac!`](macro@iupac) macros. Generally these should be dereferenced as `&'static SeqSlice`s or kmers.
//!
//! [`Kmer`](mod@kmer)s are shorter, fixed-length sequences. They generally fit in a single register and implement `Copy`. They are used for optimised algorithms on sequences and succinct datastructures. The default implementation uses a `usize` for storage. Using the 2-bit `Dna` encoding a `Kmer<Dna, 32>` occupies 64 bits.
//!
//! These sequence types are parameterised with [encodings](`codec`) (e.g. `Seq<Dna>`, `Seq<Amino>`, etc.) that define how symbols are encoded into strings of bits and decoded as readable strings.
//!
//! ## Quick start
//!
//! Add `bio-seq` to `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! bio-seq = "0.13"
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
//! Sequences are analogous to rust's string types and follow similar dereferencing conventions.
//!
//! ```rust
//! # use bio_seq::prelude::*;
//! // The `dna!` macro packs a static sequence with 2-bits per symbol at compile time.
//! // This is analogous to rust's string literals:
//! let s: &'static str = "hello!";
//! let seq: &'static SeqSlice<Dna> = dna!("CGCTAGCTACGATCGCAT");
//!
//! // Static sequences can also be copied as kmers
//! let kmer: Kmer<Dna, 18> = dna!("CGCTAGCTACGATCGCAT").into();
//! // or with the kmer! macro:
//! let kmer = kmer!("CGCTAGCTACGATCGCAT")
//!
//! // `Seq`s can be allocated on the heap like `String`s are:
//! let s: String = "hello!".into();
//! let seq: Seq<Dna> = dna!("CGCTAGCTACGATCGCAT").into();
//!
//! // Alternatively, a `Seq` can be fallibly encoded at runtime:
//! let seq: Seq<Dna> = "CGCTAGCTACGATCGCAT".try_into().unwrap();
//!
//! // &SeqSlices are analogous to &str, String slices:
//! let slice: &str = &s[1..3];
//! let seqslice: &SeqSlice<Dna> = &seq[2..4];
//! ```
//!
//! ## Bit-packed encodings
//!
//! Encodings of genomic symbols are implemented as "[`Codecs`](codec)." This crate provides four common ones:
//!   - [`codec::dna`]: 2-bit encoding of the four nucleotides
//!   - [`codec::text`]: 8-bit ASCII encoding of nucleotides, meant to be compatible with plaintext sequencing data formats
//!   - [`codec::iupac`]: 4-bit encoding of ambiguous nucleotide identities (the IUPAC ambiguity codes)
//!   - [`codec::amino`]: 6-bit encoding of amino acids
//!
//! Each of these encodings is designed to facilitate common bioinformatics tasks, such as minimising k-mers and implementing succinct datastructures. The [translation] module provides traits and methods for translating between nucleotide and amino acid sequences.
//!
//! Custom codecs can also be implemented with the `Codec` trait and derived on well crafted enums.
//!

#![warn(clippy::pedantic)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::module_name_repetitions)]
// the lint doesn't seem to recognise our implementations
#![allow(clippy::into_iter_without_iter)]

use bitvec::prelude::*;

type Order = Lsb0;
type Bs = BitSlice<usize, Order>;
type Bv = BitVec<usize, Order>;
type Ba = BitArray<usize, Order>;

pub mod codec;
pub mod error;
#[macro_use]
pub mod kmer;
pub mod seq;

//#[macro_use]
pub use bio_seq_derive::{dna, iupac};

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

    #[cfg(feature = "translation")]
    pub use crate::translation;

    pub use core::str::FromStr;

    pub use crate::error::ParseBioError;

    pub use crate::{dna, iupac, kmer};

    pub use crate::bitarr as __bio_seq_bitarr;
    pub use crate::Lsb0 as __bio_seq_Lsb0;
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::Dna::{A, C, G, T};
    use crate::prelude::*;
    use std::hash::{DefaultHasher, Hash, Hasher};

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
            Seq::<Amino>::try_from("DCMNLKGHI")
                .unwrap()
                .rev()
                .collect::<Vec<Amino>>(),
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
        let s1: &SeqSlice<Iupac> = iupac!("AS-GYTNA");
        let s2: &SeqSlice<Iupac> = iupac!("ANTGCAT-");

        let s3: &SeqSlice<Iupac> = iupac!("ACGTSWKM");
        let s4: &SeqSlice<Iupac> = iupac!("WKMSTNNA");

        assert_eq!(s1 | s2, iupac!("ANTGYWNA"));
        assert_eq!(s3 & s4, iupac!("A----WKA"));
    }
    #[test]
    fn min_sequence() {
        let seq = dna!("GCTCGATCGTAAAAAATCGTATT");

        let minimised = seq.kmers::<8>().min().unwrap();
        assert_eq!(minimised, Kmer::from(dna!("GTAAAAAA")));
    }

    #[test]
    fn hash_minimiser() {
        use core::cmp::min;

        fn hash<T: Hash>(seq: T) -> u64 {
            let mut hasher = DefaultHasher::new();
            seq.hash(&mut hasher);
            hasher.finish()
        }

        let seq =
            dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCTAAAAAAAAAAAAAAAAGGGGTGTGTGGGTTGTGGAGGAGAGAGAGCC");

        //        let minimised = seq.kmers::<16>().map(hash).min().unwrap();

        let (minimiser_rc, min_hash_rc) = seq
            .revcomp()
            .kmers::<16>()
            .map(|kmer| (kmer, hash(&kmer)))
            .min_by_key(|&(_, hash)| hash)
            .unwrap();

        let (minimiser, min_hash) = seq
            .kmers::<16>()
            .map(|kmer| (kmer, hash(&kmer)))
            .min_by_key(|&(_, hash)| hash)
            .unwrap();

        //        let x = min(min_hash, min_hash_rc);

        let (canonical_minimiser, canonical_hash) = seq
            .kmers::<16>()
            .map(|kmer| {
                let canonical_hash = min(hash(&kmer), hash(&kmer.revcomp()));
                (kmer, canonical_hash)
            })
            .min_by_key(|&(_, hash)| hash)
            .unwrap();

        println!("{minimiser_rc} {min_hash_rc}\n{minimiser} {min_hash}\n{canonical_minimiser} {canonical_hash}");
        assert_eq!(min_hash_rc, canonical_hash);
        assert_eq!(minimiser_rc, canonical_minimiser.revcomp());
    }

    #[test]
    fn hash_characteristics() {
        fn hash<T: Hash>(chunk: &T) -> u64 {
            let mut hasher = DefaultHasher::new();
            chunk.hash(&mut hasher);
            hasher.finish()
        }

        let s1 = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCT");
        let s2 = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCTA");

        let q1: Seq<Dna> = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCT").into();
        let q2: Seq<Dna> = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCTA").into();

        let s3 = dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let s4 = dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        let q3 = dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let q4 = dna!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        let l3: &SeqSlice<Dna> = &q3;
        let l3_a: &SeqSlice<Dna> = &q4[1..];
        let l3_b: &SeqSlice<Dna> = &q4[..32];
        let l4: &SeqSlice<Dna> = &q4;

        let k1: Kmer<Dna, 32> = s1.into();
        let k1_a: Kmer<Dna, 32> = s1.into();

        let k3: Kmer<Dna, 32> = s3.into();

        assert_eq!(hash(&l3), hash(q3));
        assert_eq!(hash(&l3), hash(&l3_a));
        assert_eq!(hash(&l3_a), hash(&l3_b));

        assert_eq!(hash(&s2), hash(&q2));

        assert_eq!(hash(&s1), hash(s1));
        assert_eq!(hash(s2), hash(&s2));
        assert_ne!(hash(&s4), hash(&s3));

        assert_ne!(hash(&l3), hash(&l4));
        assert_ne!(hash(&l3_a), hash(&l4));

        assert_ne!(hash(&q2), hash(&q1));

        assert_eq!(hash(q3), hash(s3));
        assert_eq!(hash(s1), hash(&q1));
        assert_ne!(hash(s3), hash(s4));

        assert_ne!(hash(&k3), hash(&k1));
        assert_eq!(hash(&k1_a), hash(&k1));
        assert_eq!(hash(s1), hash(&k1));
    }

    #[test]
    fn sequence_type_hashes() {
        fn hash<T: Hash>(chunk: &T) -> u64 {
            let mut hasher = DefaultHasher::new();
            chunk.hash(&mut hasher);
            hasher.finish()
        }

        let seq_arr: &SeqArray<Dna, 32, 1> = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCT");
        let seq: Seq<Dna> = seq_arr.into();
        let seq_slice: &SeqSlice<Dna> = &seq;
        let kmer: Kmer<Dna, 32> = seq_arr.into();

        assert_eq!(hash(seq_arr), hash(&seq));
        assert_eq!(hash(&seq), hash(&seq_slice));
        assert_eq!(hash(&seq_slice), hash(&kmer));
    }

    #[test]
    fn nth_chars() {
        assert_eq!(iupac!("ACGTRYSWKMBDHVN-").nth(0), Iupac::A);
        assert_ne!(iupac!("ACGTRYSWKMBDHVN-").nth(0), Iupac::C);
        assert_eq!(iupac!("ACGTRYSWKMBDHVN-").nth(15), Iupac::X);
        assert_eq!(iupac!("ACGTRYSWKMBDHVN-").nth(3), Iupac::from(Dna::T));
        assert_ne!(iupac!("ACGTRYSWKMBDHVN-").nth(3), Iupac::from(Dna::G));

        assert_eq!(
            Seq::<Amino>::try_from("DCMNLKGHI").unwrap().nth(1),
            Amino::C
        );
        assert_ne!(
            Seq::<Amino>::try_from("DCMNLKGHI").unwrap().nth(7),
            Amino::I
        );
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
