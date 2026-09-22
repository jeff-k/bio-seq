//! Public API tests for owned sequences.

use bio_seq::codec::text;
use bio_seq::prelude::*;
use bitvec::prelude::{BitArray, BitVec};
use std::collections::HashMap;
use std::marker::PhantomData;
use std::ops::Bound::{self, Excluded, Included, Unbounded};

const LONG: &str = "GATTACACGTTGCAAGGCTTCATGACCGTAAGTCCATGCG";

#[test]
fn default_is_an_empty_sequence() {
    let seq = Seq::<Dna>::default();

    assert!(seq.is_empty());
    assert_eq!(seq, Seq::new());
    assert_eq!(seq, Seq::with_capacity(LONG.len()));
    assert_eq!(seq.to_string(), "");
}

#[test]
fn clear_empties_the_sequence_for_reuse() {
    let mut seq: Seq<Dna> = LONG.parse().unwrap();

    seq.clear();
    assert!(seq.is_empty());
    assert_eq!(seq, Seq::new());

    seq.clear();
    assert!(seq.is_empty());

    seq.extend([Dna::G, Dna::A]);
    assert_eq!(seq, dna!("GA"));
}

#[test]
fn integer_conversion_agrees_with_slice_conversion() {
    let capacity = usize::BITS as usize / usize::from(Dna::BITS);
    let text = "TGCA".repeat(capacity / 4 + 1);

    assert_eq!(usize::from(Seq::<Dna>::new()), 0);
    for len in [0, 1, 5, capacity - 1, capacity] {
        let seq: Seq<Dna> = text[..len].parse().unwrap();
        assert_eq!(
            Ok(usize::from(seq.clone())),
            usize::try_from(&seq[..]),
            "{len}"
        );
    }
}

#[test]
#[should_panic]
fn integer_conversion_rejects_more_than_a_word() {
    let capacity = usize::BITS as usize / usize::from(Dna::BITS);
    let seq: Seq<Dna> = "A".repeat(capacity + 1).parse().unwrap();
    let _ = usize::from(seq);
}

fn bounds(len: usize) -> Vec<Bound<usize>> {
    (0..=len)
        .flat_map(|i| [Included(i), Excluded(i)])
        .chain([Unbounded])
        .collect()
}

fn is_valid_range((start, end): (Bound<usize>, Bound<usize>), len: usize) -> bool {
    let start = match start {
        Included(i) => i,
        Excluded(i) => i + 1,
        Unbounded => 0,
    };
    let end = match end {
        Included(i) => i + 1,
        Excluded(i) => i,
        Unbounded => len,
    };
    start <= end && end <= len
}

#[test]
fn remove_and_splice_agree_with_vec_for_every_kind_of_bound() {
    for text in ["", "G", LONG] {
        let seq: Seq<Dna> = text.parse().unwrap();
        let symbols: Vec<Dna> = seq.iter().collect();
        let bounds = bounds(seq.len());
        let ranges = bounds
            .iter()
            .flat_map(|&start| bounds.iter().map(move |&end| (start, end)))
            .filter(|&range| is_valid_range(range, seq.len()));

        for range in ranges {
            let mut removed = seq.clone();
            removed.remove(range);
            let mut expected = symbols.clone();
            expected.drain(range);
            assert_eq!(
                removed,
                Seq::from(&expected),
                "remove {range:?} of {text:?}"
            );

            for patch in [dna!(""), dna!("CAT")] {
                let mut spliced = seq.clone();
                spliced.splice(range, patch);
                let mut expected = symbols.clone();
                expected.splice(range, patch);
                assert_eq!(
                    spliced,
                    Seq::from(&expected),
                    "splice {range:?} of {text:?}"
                );
            }
        }
    }
}

#[test]
#[should_panic(expected = "Start of range must be less than or equal to end")]
fn remove_rejects_a_start_after_the_end() {
    let mut seq: Seq<Dna> = "ACGT".parse().unwrap();
    seq.remove((Excluded(2), Excluded(2)));
}

#[test]
#[should_panic(expected = "Range out of bounds")]
fn splice_rejects_a_range_past_the_end() {
    let mut seq: Seq<Dna> = "ACGT".parse().unwrap();
    seq.splice(2..=4, dna!("A"));
}

#[test]
#[should_panic(expected = "bound overflow")]
fn remove_rejects_an_overflowing_start_bound() {
    let mut seq: Seq<Dna> = "ACGT".parse().unwrap();
    seq.remove((Excluded(usize::MAX), Unbounded));
}

#[test]
#[should_panic(expected = "bound overflow")]
fn splice_rejects_an_overflowing_end_bound() {
    let mut seq: Seq<Dna> = "ACGT".parse().unwrap();
    seq.splice(..=usize::MAX, dna!("A"));
}

#[test]
#[should_panic(expected = "Index out of bounds")]
fn insert_rejects_a_position_past_the_end() {
    let mut seq: Seq<Dna> = "ACGT".parse().unwrap();
    seq.insert(5, dna!("A"));
}

#[test]
fn bitwise_operations_combine_iupac_symbols_pairwise() {
    // Every ordered pair of symbols: 256 of them, spanning several storage words.
    let (lhs, rhs): (Vec<Iupac>, Vec<Iupac>) = Iupac::items()
        .flat_map(|l| Iupac::items().map(move |r| (l, r)))
        .unzip();
    let pairwise = |op: fn(u8, u8) -> u8| -> Seq<Iupac> {
        lhs.iter()
            .zip(&rhs)
            .map(|(l, r)| Iupac::unsafe_from_bits(op(l.to_bits(), r.to_bits())))
            .collect()
    };
    let (left, right) = (Seq::from(&lhs), Seq::from(&rhs));

    let and = left.clone().bit_and(right.clone());
    assert_eq!(and, pairwise(|l, r| l & r));
    assert_eq!(and, &left[..] & &right[..]);

    let or = left.clone().bit_or(right.clone());
    assert_eq!(or, pairwise(|l, r| l | r));
    assert_eq!(or, &left[..] | &right[..]);

    // IUPAC codes are sets of bases, so these are intersection and union.
    let ambiguous: Seq<Iupac> = "RYSN".parse().unwrap();
    let other: Seq<Iupac> = "YRNS".parse().unwrap();
    assert_eq!(ambiguous.clone().bit_and(other.clone()), iupac!("--SS"));
    assert_eq!(ambiguous.bit_or(other), iupac!("NNNN"));
}

#[test]
fn borrowed_sequence_keys_can_be_found_by_slice() {
    let keys: Vec<Seq<Dna>> = ["ACGT", "GGCC", "CCAA"]
        .iter()
        .map(|key| key.parse::<Seq<Dna>>().unwrap())
        .collect();
    let counts: HashMap<&Seq<Dna>, usize> = keys.iter().zip(1..).collect();
    let reference: Seq<Dna> = "TTACGTGGCCAA".parse().unwrap();

    assert_eq!(counts.get(&reference[2..6]), Some(&1));
    assert_eq!(counts.get(&reference[6..10]), Some(&2));
    assert_eq!(counts.get(&reference[8..]), Some(&3));
    assert_eq!(counts.get(&reference[1..5]), None);
    assert_eq!(counts.get(&reference[2..5]), None);
}

/// Pack `text` into a `SeqArray` by hand, with the layout `dna!` generates.
fn seq_array<const N: usize, const W: usize>(text: &str) -> SeqArray<Dna, N, W> {
    let seq: Seq<Dna> = text.parse().unwrap();
    assert_eq!(seq.len(), N);
    let raw = seq.into_raw();
    SeqArray {
        _p: PhantomData,
        ba: BitArray::new(std::array::from_fn(|i| raw.get(i).copied().unwrap_or(0))),
    }
}

#[test]
fn seq_arrays_convert_to_sequences_of_any_compatible_codec() {
    const LEN: usize = LONG.len();
    const WORDS: usize = (LEN * 2).div_ceil(usize::BITS as usize);
    let array: SeqArray<Dna, LEN, WORDS> = seq_array(LONG);

    let dna: Seq<Dna> = Seq::from(&array);
    assert_eq!(dna.to_string(), LONG);

    let iupac: Seq<Iupac> = Seq::from(&array);
    assert_eq!(iupac.to_string(), LONG);

    let owned: Seq<Dna> = array.into();
    assert_eq!(owned, dna);
}

#[test]
fn every_parsing_path_agrees() {
    let cases = [
        ("", Ok(Seq::new())),
        ("GATTACA", Ok(dna!("GATTACA").to_owned())),
        ("ACGTN", Err(ParseBioError::UnrecognisedBase(b'N'))),
        ("acgt", Err(ParseBioError::UnrecognisedBase(b'a'))),
        ("AC GT", Err(ParseBioError::UnrecognisedBase(b' '))),
        ("ACGTé", Err(ParseBioError::UnrecognisedBase(0xC3))),
    ];

    for (text, expected) in cases {
        let owned = text.to_owned();
        let bytes = text.as_bytes().to_vec();
        assert_eq!(text.parse::<Seq<Dna>>(), expected, "{text:?}");
        assert_eq!(Seq::<Dna>::try_from(text), expected, "{text:?}");
        assert_eq!(Seq::<Dna>::try_from(&owned), expected, "{text:?}");
        assert_eq!(Seq::<Dna>::try_from(owned), expected, "{text:?}");
        assert_eq!(Seq::<Dna>::try_from(text.as_bytes()), expected, "{text:?}");
        assert_eq!(Seq::<Dna>::try_from(bytes), expected, "{text:?}");
    }

    let not_utf8 = vec![b'A', 0xFF];
    assert_eq!(
        Seq::<Dna>::try_from(not_utf8),
        Err(ParseBioError::UnrecognisedBase(0xFF))
    );
}

#[test]
fn string_conversions_match_display() {
    for text in ["", "A", "GATTACA-NRYKMSWBDHV"] {
        let seq: Seq<Iupac> = text.parse().unwrap();
        assert_eq!(seq.to_string(), text);
        assert_eq!(String::from(&seq), text);
        assert_eq!(String::from(seq), text);
    }
}

#[test]
fn text_sequences_reinterpret_words_as_bytes() {
    let bytes: Vec<u8> = b"GATTACANNACGT"
        .iter()
        .copied()
        .cycle()
        .take(3 * size_of::<usize>())
        .collect();
    // Lsb0 order: the first symbol of each word is its least significant byte.
    let words: Vec<usize> = bytes
        .chunks(size_of::<usize>())
        .map(|chunk| usize::from_le_bytes(chunk.try_into().unwrap()))
        .collect();

    let seq = Seq::<text::Dna>::from(words);
    assert_eq!(seq.len(), bytes.len());
    assert_eq!(seq, Seq::<text::Dna>::try_from(&bytes[..]).unwrap());
    assert!(Seq::<text::Dna>::from(Vec::<usize>::new()).is_empty());
}

#[test]
fn bit_vectors_convert_symbol_by_symbol() {
    let seq: Seq<Dna> = LONG.parse().unwrap();
    // Each base contributes its two bits, least significant first.
    let bits: BitVec = seq
        .iter()
        .flat_map(|base| {
            let b = base.to_bits();
            [b & 1 == 1, b & 2 == 2]
        })
        .collect();

    assert_eq!(Seq::<Dna>::from(bits.as_bitslice()), seq);
    // Starts part-way through the first storage word.
    assert_eq!(Seq::<Dna>::from(&bits[6..]), seq[3..]);
    assert_eq!(Seq::<Dna>::from(bits), seq);
}

#[test]
fn from_raw_checks_capacity_without_overflowing() {
    let per_word = usize::BITS as usize / usize::from(Dna::BITS);
    let words = [usize::MAX, 0b11_10_01_00];

    let seq = Seq::<Dna>::from_raw(per_word + 4, &words).unwrap();
    assert_eq!(seq.to_string(), "T".repeat(per_word) + "ACGT");
    assert_eq!(
        Seq::<Dna>::from_raw(2 * per_word, &words).map(|seq| seq.len()),
        Some(2 * per_word)
    );
    assert_eq!(Seq::<Dna>::from_raw(2 * per_word + 1, &words), None);
    // Lengths whose size in bits overflows a usize, including one that would wrap to zero.
    assert_eq!(Seq::<Dna>::from_raw(usize::MAX, &words), None);
    assert_eq!(Seq::<Dna>::from_raw(usize::MAX / 2 + 1, &words), None);
    assert_eq!(Seq::<Dna>::from_raw(0, &[]), Some(Seq::new()));
}
