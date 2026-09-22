// Copyright 2021-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![doc = include_str!("../README.md")]
#![warn(clippy::pedantic)]
//#[cfg(not(target_pointer_width = "64"))]
//compile_error!("bio-seq currently only supports 64-bit platforms");
//#![feature(simd_wasm64)]
//#![feature(portable_simd)]

use bitvec::prelude::*;

type Order = Lsb0;
type Bs = BitSlice<usize, Order>;
type Bv = BitVec<usize, Order>;
type Ba<const W: usize> = BitArray<[usize; W], Order>;

mod hash;

pub mod codec;
pub mod error;
#[macro_use]
pub mod kmer;
pub mod seq;

pub use bio_seq_derive::{dna, iupac};

#[doc(hidden)]
pub use bitvec::bitarr as __bio_seq_bitarr;

#[doc(hidden)]
pub use bitvec::prelude::Lsb0 as __bio_seq_Lsb0;

#[cfg(feature = "translation")]
pub mod translation;

pub mod prelude {
    pub use crate::codec::Codec;
    pub use crate::codec::amino::Amino;
    pub use crate::codec::dna::Dna;
    pub use crate::codec::iupac::Iupac;
    pub use crate::{
        Complement, ComplementMut, Maskable, MaskableMut, Reverse, ReverseComplement,
        ReverseComplementMut, ReverseMut,
    };

    pub use crate::kmer::Kmer;
    pub use crate::seq::{Seq, SeqArray, SeqSlice};

    #[cfg(feature = "translation")]
    pub use crate::translation;

    pub use core::str::FromStr;

    pub use crate::error::ParseBioError;

    pub use crate::{dna, iupac, kmer};

    #[doc(hidden)]
    pub use crate::__bio_seq_Lsb0;
    #[doc(hidden)]
    pub use crate::__bio_seq_bitarr;
    #[doc(hidden)]
    pub use crate::__bio_seq_count_words;
}

/// Nucleotide bases and sequences can be complemented
pub trait ComplementMut {
    fn comp(&mut self);
}

pub trait Complement: ComplementMut + ToOwned
where
    <Self as ToOwned>::Owned: ComplementMut,
{
    /// ```
    /// use bio_seq::prelude::{Dna, Complement};
    /// assert_eq!(Dna::A.to_comp(), Dna::T);
    /// ````
    #[must_use]
    fn to_comp(&self) -> <Self as ToOwned>::Owned {
        let mut owned = self.to_owned();
        owned.comp();
        owned
    }
}

//impl<T: ?Sized + ComplementMut + ToOwned> Complement for T where <T as ToOwned>::Owned: ComplementMut {}

/// A reversible sequence
pub trait ReverseMut {
    /// Reverse sequence in place
    fn rev(&mut self);
}

pub trait Reverse: ReverseMut + ToOwned
where
    <Self as ToOwned>::Owned: ReverseMut,
{
    #[must_use]
    fn to_rev(&self) -> <Self as ToOwned>::Owned {
        let mut owned = self.to_owned();
        owned.rev();
        owned
    }
}

//impl<T: ?Sized + ReverseMut + ToOwned> Reverse for T where <T as ToOwned>::Owned: ReverseMut {}

/// A reversible sequence that can be complemented can be reverse complemented
pub trait ReverseComplementMut: ComplementMut + ReverseMut {
    /// Reverse complement a sequence in place
    fn revcomp(&mut self) {
        self.comp();
        self.rev();
    }
}

//impl<T: ?Sized + ReverseMut + ComplementMut> ReverseComplementMut for T {}

pub trait ReverseComplement: ReverseComplementMut + ToOwned
where
    <Self as ToOwned>::Owned: ReverseComplementMut,
{
    #[must_use]
    fn to_revcomp(&self) -> <Self as ToOwned>::Owned {
        let mut owned = self.to_owned();
        owned.revcomp();
        owned
    }
}

// TODO: Marker trait to allow overriding this blanket impl
//impl<T: ReverseComplementMut + ToOwned> ReverseComplement for T where
//    <T as ToOwned>::Owned: ReverseComplementMut
//{
//}

/// Some sequence types may be maskable
pub trait MaskableMut {
    fn mask(&mut self);
    fn unmask(&mut self);
}

pub trait Maskable: MaskableMut + ToOwned
where
    <Self as ToOwned>::Owned: MaskableMut,
{
    #[must_use]
    fn to_mask(&self) -> <Self as ToOwned>::Owned {
        let mut owned = self.to_owned();
        owned.mask();
        owned
    }
    #[must_use]
    fn to_unmask(&self) -> <Self as ToOwned>::Owned {
        let mut owned = self.to_owned();
        owned.unmask();
        owned
    }
}

//impl<T: MaskableMut + ToOwned> Maskable for T where <T as ToOwned>::Owned: MaskableMut {}

#[macro_export]
macro_rules! __bio_seq_count_words {
    ($len:expr) => {{ $len.div_ceil(usize::BITS) as usize }};
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

    #[test]
    fn into_usize() {
        let bits: usize = dna!("ACGT").to_owned().into_raw()[0];
        assert_eq!(bits, 0b11_10_01_00);

        let bits: usize = dna!("CGCG").to_owned().into_raw()[0];
        assert_eq!(bits, 0b10_01_10_01);

        let bits: usize = Seq::from(&vec![T, T]).into();
        assert_eq!(bits, 0b11_11);

        let bits: usize = Seq::<Dna>::from_str("TCA").unwrap().into();
        assert_eq!(bits, 0b00_01_11);

        let bits: usize = Seq::<Dna>::from_str("TGA").unwrap().into();
        assert_eq!(bits, 0b00_10_11);

        let bits: usize = Seq::from(&vec![C, G, T, A, C, G, A, T]).into();
        assert_eq!(bits, 0b11_00_10_01_00_11_10_01);

        let bits: usize = Seq::from(&vec![A]).into();
        assert_eq!(bits, 0b00);
    }

    #[test]
    fn test_display_aminos() {
        let a: Seq<Amino> = Seq::from_str("DCMNLKG*HI").unwrap();
        assert_eq!(format!("{a}"), "DCMNLKG*HI");
    }
    #[test]
    fn test_display_dna() {
        let seq = Seq::from(&vec![A, C, G, T, T, A, T, C]);
        assert_eq!(format!("{seq}"), "ACGTTATC");
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
            seq.rev_iter().collect::<Vec<Dna>>(),
            vec![T, G, C, A, T, G, C, A]
        );
        assert_eq!(
            seq.to_rev().into_iter().collect::<Vec<Dna>>(),
            vec![T, G, C, A, T, G, C, A]
        );
        assert_eq!(
            iupac!("GN-").rev_iter().collect::<Vec<Iupac>>(),
            vec![Iupac::X, Iupac::N, Iupac::G]
        );

        assert_eq!(
            Seq::<Amino>::try_from("DCMNLKGHI")
                .unwrap()
                .to_rev()
                .into_iter()
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
            assert_eq!(format!("{kmer}"), answer);
        }
    }

    #[test]
    fn iterate_kmer8() {
        let seq = dna!("AAAACCCCGGGG");
        for (kmer, answer) in seq
            .kmers::<8>()
            .zip(["AAAACCCC", "AAACCCCG", "AACCCCGG", "ACCCCGGG", "CCCCGGGG"])
        {
            assert_eq!(format!("{kmer}"), answer);
        }
    }

    #[test]
    fn iterate_kmer4() {
        let seq = dna!("AAAACCCCGGGGTTTT");
        for (kmer, answer) in seq.kmers::<4>().zip([
            "AAAA", "AAAC", "AACC", "ACCC", "CCCC", "CCCG", "CCGG", "CGGG", "GGGG", "GGGT", "GGTT",
            "GTTT", "TTTT",
        ]) {
            assert_eq!(format!("{kmer}"), answer);
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
        assert_eq!(minimised, Kmer::try_from(dna!("GTAAAAAA")).unwrap());
    }

    #[test]
    fn hash_minimiser() {
        use core::cmp::min;

        fn hash(seq: &SeqSlice<Dna>) -> u64 {
            u64::from(seq != dna!("GGCTCTCTCTCCTCCA"))
        }

        let seq =
            dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCTAAAAAAAAAAAAAAAAGGGGTGTGTGGGTTGTGGAGGAGAGAGAGCC");

        //        let minimised = seq.kmers::<16>().map(hash).min().unwrap();

        let (minimiser_rc, min_hash_rc) = seq
            .to_revcomp()
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
                let canonical_hash = min(hash(&kmer), hash(&kmer.to_revcomp()));
                (kmer, canonical_hash)
            })
            .min_by_key(|&(_, hash)| hash)
            .unwrap();

        println!(
            "{minimiser_rc} {min_hash_rc}\n{minimiser} {min_hash}\n{canonical_minimiser} {canonical_hash}"
        );
        assert_eq!(min_hash_rc, canonical_hash);
        assert_ne!(min_hash, canonical_hash);
        assert_eq!(minimiser_rc, canonical_minimiser.to_revcomp());
    }

    #[test]
    #[expect(
        clippy::needless_borrows_for_generic_args,
        reason = "exercise hashing both values and references"
    )]
    fn hash_characteristics() {
        fn hash<T: Hash>(chunk: T) -> u64 {
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

        let l3: &SeqSlice<Dna> = q3;
        let l3_a: &SeqSlice<Dna> = &q4[1..];
        let l3_b: &SeqSlice<Dna> = &q4[..32];
        let l4: &SeqSlice<Dna> = q4;

        let k1: Kmer<Dna, 32, u64> = s1.try_into().unwrap();
        let k1_a: Kmer<Dna, 32, u64> = s1.try_into().unwrap();

        let k3: Kmer<Dna, 32, u64> = s3.try_into().unwrap();

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

        let seq_arr: &SeqSlice<Dna> = dna!("AGCGCTAGTCGTACTGCCGCATCGCTAGCGCT");
        let seq: Seq<Dna> = seq_arr.into();
        let seq_slice: &SeqSlice<Dna> = &seq;
        let kmer: Kmer<Dna, 32, u64> = seq_arr.try_into().unwrap();

        assert_eq!(hash(&seq_arr), hash(&seq));
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
            assert_eq!(format!("{}", Kmer::<Dna, 2>::from(i)), e.to_string());
            assert_eq!(Kmer::<Dna, 2>::from(i), *e);
        }
    }

    #[test]
    #[expect(
        clippy::redundant_slicing,
        reason = "exercise equality with full-range indexing"
    )]
    #[expect(
        clippy::similar_names,
        reason = "names identify sequence variants and k-mer storage widths"
    )]
    fn sequence_type_equality() {
        let raw_a = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAA";
        let raw_b = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA";
        let raw_c = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA";
        let raw_d = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAAA";

        assert_eq!(raw_a.len(), 63);
        assert_eq!(raw_b.len(), 64);
        assert_eq!(raw_d.len(), 65);

        assert_eq!(raw_c, raw_b);
        assert_eq!(raw_c, &raw_b[..]);

        assert_ne!(raw_b, raw_d);
        assert_ne!(raw_a, raw_b);

        // Seq

        let seq_a: Seq<Dna> = raw_a.try_into().unwrap();
        let seq_b: Seq<Dna> = raw_b.try_into().unwrap();
        let seq_c: Seq<Dna> = raw_c.try_into().unwrap();
        let seq_d: Seq<Dna> = raw_d.try_into().unwrap();

        assert_eq!(seq_a.len(), raw_a.len());
        assert_eq!(seq_d.len(), raw_d.len());

        assert_eq!(seq_c, seq_b);
        assert_eq!(seq_c, &seq_b);

        assert_ne!(seq_a, &seq_b);
        assert_ne!(seq_a, seq_b);
        assert_ne!(seq_c, seq_d);

        // SeqSlice

        let slice_a: &SeqSlice<Dna> = &seq_a;
        let slice_b: &SeqSlice<Dna> = &seq_b;
        let slice_c: &SeqSlice<Dna> = &seq_c;
        let slice_d: &SeqSlice<Dna> = &seq_d;

        assert_eq!(slice_a.len(), raw_a.len());
        assert_eq!(slice_d.len(), raw_d.len());

        assert_eq!(slice_c, slice_b);
        assert_eq!(slice_c, &slice_b[..]);

        assert_ne!(slice_a, slice_b);
        assert_ne!(slice_c, slice_d);
        assert_ne!(slice_c, &slice_d[..]);

        // SeqArray references

        let array_a = dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAA");
        let array_b = dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA");
        let array_c: &'static SeqSlice<Dna> =
            dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA");
        let array_d: &'static SeqSlice<Dna> =
            dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAAA");

        assert_eq!(array_a.len(), raw_a.len());
        assert_eq!(array_d.len(), raw_d.len());

        assert_eq!(array_c, array_b);

        assert_ne!(array_a, array_b);
        assert_ne!(array_c, array_d);

        // Kmers

        let kmer_ax_32: Kmer<Dna, 32, u64> = kmer!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAG", u64);
        let kmer_bx_32 = Kmer::<Dna, 32, u64>::from_str(&raw_b[..32]).unwrap();

        let kmer_x_32: Kmer<Dna, 32, u64> = kmer!("AATTGTGGGTTCGTCTGCGCCTCCGCCCTTAG", u64);

        assert_eq!(kmer_ax_32.len(), 32);

        assert_eq!(kmer_ax_32, kmer_bx_32);
        assert_ne!(kmer_ax_32, kmer_x_32);

        let kmer_b_64 = Kmer::<Dna, 64, u128>::from_str(raw_b).unwrap();
        let kmer_cx_64 = Kmer::<Dna, 64, u128>::from_str(&raw_d[..64]).unwrap();
        let kmer_dx_64 = Kmer::<Dna, 64, u128>::from_str(&raw_d[1..]).unwrap();

        assert_eq!(kmer_cx_64.len(), 64);

        assert_eq!(kmer_b_64, kmer_cx_64);
        assert_ne!(kmer_b_64, kmer_dx_64);

        // Cross-type equality:

        assert_eq!(seq_c, slice_b);
        assert_eq!(seq_c, *array_b);
        //        assert_eq!(seq_c, kmer_b_64);

        assert_eq!(&seq_c, slice_b);
        assert_eq!(&seq_c, array_b);
        //        assert_eq!(&seq_c, kmer_b_64);

        assert_eq!(slice_c, seq_b);
        assert_eq!(slice_c, &seq_b);
        assert_eq!(&slice_c, array_b);
        //        assert_eq!(slice_c, kmer_b_64);

        assert_eq!(array_c, &seq_b);
        assert_eq!(array_c, seq_b);
        assert_eq!(array_c, slice_b);
        //        assert_eq!(array_c, kmer_b_64);

        //        assert_eq!(kmer_b_64, &seq_c);
        //        assert_eq!(kmer_b_64, seq_c);
        //        assert_eq!(kmer_b_64, slice_c);
        //        assert_eq!(kmer_b_64, array_c);

        // Cross-type inequality (shorter):

        assert_ne!(&seq_a, slice_b);
        assert_ne!(&seq_a, array_b);
        assert_ne!(seq_a, slice_b);
        assert_ne!(seq_a, array_b);
        //        assert_ne!(seq_a, kmer_b_64);
        //        assert_ne!(&seq_a, kmer_b_64);

        assert_ne!(slice_a, &seq_b);
        assert_ne!(slice_a, seq_b);
        assert_ne!(&slice_a, array_b);
        //        assert_ne!(slice_a, kmer_b_64);

        assert_ne!(array_a, &seq_b);
        assert_ne!(array_a, seq_b);
        assert_ne!(array_a, slice_b);
        //        assert_ne!(array_a, kmer_b_64);

        //        assert_ne!(kmer_b_64, &seq_a);
        //        assert_ne!(kmer_b_64, seq_a);
        //        assert_ne!(kmer_b_64, slice_a);
        //        assert_ne!(kmer_b_64, array_a);

        // Cross-type inequality (longer):

        assert_ne!(seq_d, slice_b);
        assert_ne!(seq_d, array_b);
        //        assert_ne!(seq_d, kmer_b_64);
        assert_ne!(&seq_d, slice_b);
        assert_ne!(&seq_d, array_b);
        //        assert_ne!(&seq_d, kmer_b_64);

        assert_ne!(slice_d, &seq_b);
        assert_ne!(slice_d, seq_b);
        assert_ne!(&slice_d, array_b);
        //        assert_ne!(slice_d, kmer_b_64);

        assert_ne!(array_d, &seq_b);
        assert_ne!(array_d, seq_b);

        assert_ne!(slice_b, array_d);
        assert_ne!(array_d, slice_b);
        //        assert_ne!(array_d, kmer_b_64);

        //        assert_ne!(kmer_b_64, &seq_d);
        //        assert_ne!(kmer_b_64, seq_d);
        //        assert_ne!(kmer_b_64, slice_d);
        //        assert_ne!(kmer_b_64, array_d);
    }
}

#[cfg(test)]
#[cfg(target_arch = "wasm32")]
mod wasm_tests {
    use crate::prelude::*;
    use wasm_bindgen_test::*;

    #[wasm_bindgen_test]
    fn sequence_type_equality() {
        let raw_a = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAA";
        let raw_b = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA";
        let raw_c = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA";
        let raw_d = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAAA";

        assert_eq!(raw_a.len(), 63);
        assert_eq!(raw_b.len(), 64);
        assert_eq!(raw_d.len(), 65);

        assert_eq!(raw_c, raw_b);
        assert_eq!(raw_c, &raw_b[..]);

        assert_ne!(raw_b, raw_d);
        assert_ne!(raw_a, raw_b);

        // Seq

        let seq_a: Seq<Dna> = raw_a.try_into().unwrap();
        let seq_b: Seq<Dna> = raw_b.try_into().unwrap();
        let seq_c: Seq<Dna> = raw_c.try_into().unwrap();
        let seq_d: Seq<Dna> = raw_d.try_into().unwrap();

        assert_eq!(seq_a.len(), raw_a.len());
        assert_eq!(seq_d.len(), raw_d.len());

        assert_eq!(seq_c, seq_b);
        assert_eq!(seq_c, &seq_b);

        assert_ne!(seq_a, &seq_b);
        assert_ne!(seq_a, seq_b);
        assert_ne!(seq_c, seq_d);

        // SeqSlice

        let slice_a: &SeqSlice<Dna> = &seq_a;
        let slice_b: &SeqSlice<Dna> = &seq_b;
        let slice_c: &SeqSlice<Dna> = &seq_c;
        let slice_d: &SeqSlice<Dna> = &seq_d;

        assert_eq!(slice_a.len(), raw_a.len());
        assert_eq!(slice_d.len(), raw_d.len());

        assert_eq!(slice_c, slice_b);
        assert_eq!(slice_c, &slice_b[..]);

        assert_ne!(slice_a, slice_b);
        assert_ne!(slice_c, slice_d);
        assert_ne!(slice_c, &slice_d[..]);

        // SeqArray references

        let array_a = dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAA");
        let array_b = dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA");
        let array_c: &'static SeqSlice<Dna> =
            dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA");
        let array_d: &'static SeqSlice<Dna> =
            dna!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAAA");

        assert_eq!(array_a.len(), raw_a.len());
        assert_eq!(array_d.len(), raw_d.len());

        assert_eq!(array_c, array_b);

        assert_ne!(array_a, array_b);
        assert_ne!(array_c, array_d);

        assert_eq!(seq_c, slice_b);
        assert_eq!(seq_c, *array_b);

        assert_eq!(&seq_c, slice_b);
        assert_eq!(&seq_c, array_b);

        assert_eq!(slice_c, seq_b);
        assert_eq!(slice_c, &seq_b);
        assert_eq!(&slice_c, array_b);

        assert_eq!(array_c, &seq_b);
        assert_eq!(array_c, seq_b);
        assert_eq!(array_c, slice_b);
        // Cross-type inequality (shorter):

        assert_ne!(&seq_a, slice_b);
        assert_ne!(&seq_a, array_b);
        assert_ne!(seq_a, slice_b);
        assert_ne!(seq_a, array_b);
        assert_ne!(slice_a, &seq_b);
        assert_ne!(slice_a, seq_b);
        assert_ne!(&slice_a, array_b);

        assert_ne!(array_a, &seq_b);
        assert_ne!(array_a, seq_b);
        assert_ne!(array_a, slice_b);
        // Cross-type inequality (longer):

        assert_ne!(seq_d, slice_b);
        assert_ne!(seq_d, array_b);
        assert_ne!(&seq_d, slice_b);
        assert_ne!(&seq_d, array_b);

        assert_ne!(slice_d, &seq_b);
        assert_ne!(slice_d, seq_b);
        assert_ne!(&slice_d, array_b);

        assert_ne!(array_d, &seq_b);
        assert_ne!(array_d, seq_b);

        assert_ne!(slice_b, array_d);
        assert_ne!(array_d, slice_b);
    }

    #[wasm_bindgen_test]
    fn wasm_kmers() {
        //let raw_a = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAA";
        let raw_b = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA";
        //let raw_c = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAA";
        let raw_d = "AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAGTACTATAGGACGATCAGCACCATAAGAACAAAA";

        let kmer_ax_32: Kmer<Dna, 32, u64> = kmer!("AATTGTGGGTTCGTCTGCGGCTCCGCCCTTAG", u64);
        let kmer_bx_32 = Kmer::<Dna, 32, u64>::from_str(&raw_b[..32]).unwrap();

        let kmer_x_32: Kmer<Dna, 32, u64> = kmer!("AATTGTGGGTTCGTCTGCGCCTCCGCCCTTAG", u64);

        assert_eq!(kmer_ax_32.len(), 32);

        assert_eq!(kmer_ax_32, kmer_bx_32);
        assert_ne!(kmer_ax_32, kmer_x_32);

        let kmer_b_64 = Kmer::<Dna, 64, u128>::from_str(&raw_b).unwrap();
        let kmer_cx_64 = Kmer::<Dna, 64, u128>::from_str(&raw_d[..64]).unwrap();
        let kmer_dx_64 = Kmer::<Dna, 64, u128>::from_str(&raw_d[1..]).unwrap();

        assert_eq!(kmer_cx_64.len(), 64);

        assert_eq!(kmer_b_64, kmer_cx_64);
        assert_ne!(kmer_b_64, kmer_dx_64);
    }

    /*
        #[wasm_bindgen_test]
        fn test_splice() {
            let mut seq: Seq<Dna> = dna!("TCAGCATCGATCAATCG").into();
            let insertion = dna!("CCCCC");

            seq.splice(4..6, insertion);
            assert_eq!(&seq, dna!("TCAGCCCCCTCGATCAATCG"));

            seq.splice(1..=1, dna!("AAA"));
            assert_eq!(&seq, dna!("TAAAAGCCCCCTCGATCAATCG"));

            seq.splice(10.., dna!("TTTT"));
            assert_eq!(&seq, dna!("TAAAAGCCCCTTTT"));
        }
    */
}
