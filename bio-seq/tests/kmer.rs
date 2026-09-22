use bio_seq::prelude::*;
use bitvec::array::BitArray;
use std::marker::PhantomData;

/// `ACGT` packed with the layout `dna!` generates.
fn acgt_array() -> SeqArray<Dna, 4, 1> {
    SeqArray {
        _p: PhantomData,
        ba: BitArray::new([0b11_10_01_00]),
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Codec)]
#[repr(u8)]
#[bits(3)]
enum Sparse {
    A = 1,
    C = 2,
    X = 5,
}

macro_rules! storage_tests {
    ($module:ident, $storage:ty) => {
        mod $module {
            use super::*;

            #[cfg(feature = "extra_codecs")]
            #[test]
            fn reverse_ws_kmer() {
                use bio_seq::codec::degenerate::WS;
                let k: Kmer<WS, 2, $storage> = "WS".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "SW");
            }

            #[test]
            fn reverse_dna_kmer() {
                let k: Kmer<Dna, 3, $storage> = "ACG".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "GCA");
                assert_eq!(k.to_rev().to_comp().to_string(), "CGT");

                const FULL_K: usize = <$storage>::BITS as usize / 2;
                let fk: Kmer<Dna, FULL_K, $storage> = "ACGT".repeat(FULL_K / 4).parse().unwrap();
                assert_eq!(fk.to_comp().to_string(), "TGCA".repeat(FULL_K / 4));
            }

            #[test]
            #[should_panic]
            fn reverse_sparse_kmer() {
                let k: Kmer<Sparse, 3, $storage> = "ACX".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "XCA");
            }

            #[test]
            fn reverse_iupac_kmer() {
                let k: Kmer<Iupac, 2, $storage> = "AC".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "CA");
            }

            #[test]
            fn parse_checks_length_before_symbols() {
                type K4 = Kmer<Dna, 4, $storage>;

                for (text, len) in [("", 0), ("ACG", 3), ("ACGTA", 5), ("NNN", 3)] {
                    let expected = Err(ParseBioError::MismatchedLength(4, len));
                    assert_eq!(text.parse::<K4>(), expected, "{text:?}");
                }
                assert_eq!(
                    "ACGN".parse::<K4>(),
                    Err(ParseBioError::UnrecognisedBase(b'N'))
                );
                assert_eq!("ACGT".parse::<K4>().unwrap().to_string(), "ACGT");
            }

            #[test]
            fn equals_slices_of_the_same_sequence() {
                let seq: Seq<Dna> = "TACGTA".parse().unwrap();
                let k: Kmer<Dna, 4, $storage> = "ACGT".parse().unwrap();

                assert_eq!(k, seq[1..5]);
                assert_eq!(k, &seq[1..5]);
                assert_ne!(k, seq[0..4]);
                assert_ne!(k, seq[1..4]);
                assert_ne!(k, seq[1..]);
                assert_ne!(k, &seq[1..]);
            }

            #[test]
            fn equals_seq_arrays() {
                let array = acgt_array();
                let k: Kmer<Dna, 4, $storage> = "ACGT".parse().unwrap();
                let other: Kmer<Dna, 4, $storage> = "ACGA".parse().unwrap();

                assert_eq!(k, array);
                assert_eq!(k, &array);
                assert_ne!(other, array);
                assert_ne!(other, &array);
            }
        }
    };
}

storage_tests!(storage_usize, usize);
storage_tests!(storage_u64, u64);
storage_tests!(storage_u128, u128);

#[test]
fn hand_packed_array_matches_macro() {
    assert_eq!(*acgt_array(), *dna!("ACGT"));
}

#[test]
fn u64_kmers_from_integers() {
    // Two bits per base, first base in the lowest bits.
    let k = Kmer::<Dna, 4, u64>::from(0b11_10_01_00_u64);
    assert_eq!(k.to_string(), "ACGT");
    assert_eq!(k, "ACGT".parse::<Kmer<Dna, 4, u64>>().unwrap());
    assert_eq!(Kmer::<Dna, 4, u64>::from(0b11_10_01_00_usize), k);
    assert_eq!(
        Kmer::<Dna, 4>::from(0b11_10_01_00_usize).to_string(),
        "ACGT"
    );

    // 0x9c is ATCG.
    let full = Kmer::<Dna, 32, u64>::from(0x9c9c_9c9c_9c9c_9c9c_u64);
    assert_eq!(full.to_string(), "ATCG".repeat(8));
    let half = Kmer::<Dna, 16, u64>::from(0x9c9c_9c9c_usize);
    assert_eq!(half.to_string(), "ATCG".repeat(4));
}

#[test]
fn kmers_convert_to_sequences() {
    let k = kmer!("GATTACA");
    let seq = Seq::<Dna>::from(k);
    assert_eq!(seq, dna!("GATTACA"));
    assert_eq!(Kmer::<Dna, 7>::try_from(seq), Ok(k));

    const FULL_K: usize = usize::BITS as usize / 2;
    let full: Kmer<Dna, FULL_K> = "ACGT".repeat(FULL_K / 4).parse().unwrap();
    assert_eq!(Seq::from(full).to_string(), "ACGT".repeat(FULL_K / 4));

    let amino: Kmer<Amino, 5> = "MWKL*".parse().unwrap();
    assert_eq!(Seq::from(amino).to_string(), "MWKL*");
}

#[test]
fn kmers_equal_owned_sequences_of_the_same_length() {
    let k = kmer!("ACGT");

    for (text, equal) in [
        ("ACGT", true),
        ("ACGA", false),
        ("", false),
        ("ACG", false),
        ("ACGTA", false),
    ] {
        let seq: Seq<Dna> = text.parse().unwrap();
        assert_eq!(k == seq, equal, "{text:?}");
    }
}

#[test]
fn kmer_iterator_agrees_with_windows() {
    let seq: Seq<Dna> = "GATTACACGTTGCAAGGCTTCATGACCGTAAGTCCATGCG".parse().unwrap();

    let kmers: Vec<Seq<Dna>> = seq.kmers::<5>().map(Seq::from).collect();
    let windows: Vec<Seq<Dna>> = seq.windows(5).collect();
    assert_eq!(kmers.len(), seq.len() - 4);
    assert_eq!(kmers, windows);
}

#[test]
fn kmer_iterator_reports_its_remaining_length() {
    let seq = dna!("ACGTACGTAC");
    let mut kmers = seq.kmers::<4>();

    assert_eq!(kmers.len(), 7);
    kmers.next();
    assert_eq!(kmers.len(), 6);
    assert_eq!(kmers.by_ref().count(), 6);
    assert_eq!(kmers.len(), 0);
    assert_eq!(kmers.next(), None);

    for (text, count) in [("", 0), ("ACG", 0), ("ACGT", 1)] {
        let seq: Seq<Dna> = text.parse().unwrap();
        assert_eq!(seq.kmers::<4>().len(), count, "{text:?}");
        assert_eq!(seq.kmers::<4>().count(), count, "{text:?}");
    }
}
