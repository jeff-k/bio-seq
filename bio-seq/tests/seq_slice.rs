//! Public API tests for borrowed sequence slices.

use bio_seq::prelude::*;
use std::fmt::{self, Write};
use std::hash::{Hash, Hasher};

#[test]
fn empty_slice_access_and_conversions() {
    let empty = &dna!("TACG")[2..2];

    assert_eq!(empty.len(), 0);
    assert!(empty.is_empty());
    assert_eq!(empty.get(0), None);
    assert_eq!(empty.get(usize::MAX), None);
    assert_eq!(usize::try_from(empty), Ok(0));
    assert_eq!(u8::from(empty), 0);
    assert_eq!(String::from(empty), "");
    assert_eq!(empty.to_string(), "");
    assert_eq!(*empty, "");
    assert_eq!(empty.to_owned(), Seq::<Dna>::new());
    assert!(std::ptr::eq(empty, SeqSlice::as_ref(empty)));
}

#[test]
fn access_uses_symbol_indices_relative_to_the_slice() {
    // Six-bit symbols include a symbol spanning a storage-word boundary.
    let seq: Seq<Amino> = "AMWLPQHIKRST*V".parse().unwrap();
    let slice = &seq[1..13];
    let expected = [
        Amino::M,
        Amino::W,
        Amino::L,
        Amino::P,
        Amino::Q,
        Amino::H,
        Amino::I,
        Amino::K,
        Amino::R,
        Amino::S,
        Amino::T,
        Amino::X,
    ];

    assert_eq!(slice.len(), expected.len());
    assert!(!slice.is_empty());
    for (i, symbol) in expected.into_iter().enumerate() {
        assert_eq!(slice.get(i), Some(symbol), "symbol {i}");
        assert_eq!(slice.nth(i), symbol, "symbol {i}");
    }
    assert_eq!(slice.get(slice.len()), None);
    assert_eq!(slice.get(usize::MAX), None);
    assert!(std::ptr::eq(slice, SeqSlice::as_ref(slice)));
}

#[test]
#[should_panic]
fn nth_past_the_slice_end_panics() {
    let slice = &dna!("TACGTG")[1..5];
    let _ = slice.nth(slice.len());
}

#[test]
fn integer_conversions_respect_slice_offsets() {
    let seq = dna!("TACGTG");
    let byte = &seq[1..5];
    let partial_byte = &seq[2..5];

    assert_eq!(u8::from(byte), 0b11_10_01_00);
    assert_eq!(usize::try_from(byte), Ok(0b11_10_01_00));
    assert_eq!(u8::from(partial_byte), 0b11_10_01);
    assert_eq!(usize::try_from(partial_byte), Ok(0b11_10_01));
}

#[test]
fn usize_conversion_checks_capacity_in_symbols() {
    let capacity = usize::BITS as usize / usize::from(Dna::BITS);
    let seq: Seq<Dna> = format!("A{}C", "T".repeat(capacity)).parse().unwrap();

    // A full integer's worth of symbols, starting two bits into the backing storage.
    let full = &seq[1..capacity + 1];
    assert_eq!(usize::try_from(full), Ok(usize::MAX));
    assert_eq!(
        usize::try_from(&seq[1..]),
        Err(ParseBioError::SequenceTooLong(capacity, capacity + 1))
    );
}

#[test]
#[should_panic]
fn u8_conversion_rejects_more_than_eight_bits() {
    let _ = u8::from(dna!("ACGTA"));
}

#[test]
fn string_equality_checks_length_and_every_symbol() {
    let slice = &dna!("TACGTG")[1..5];

    assert_eq!(*slice, "ACGT");
    for different in ["", "ACG", "ACGTA", "TCGT", "ACGA", "acgt", "ACNT", "ACé"] {
        assert_ne!(*slice, different, "{different:?}");
    }
    assert_eq!(*iupac!("R-N"), "R-N");
}

#[test]
fn sequence_equality_ignores_storage_offset() {
    let slice = &dna!("TACGTG")[1..5];
    let same = &dna!("CCACGTA")[2..6];

    for (other, equal) in [
        (same, true),
        (dna!("ACGA"), false),
        (dna!("ACGTA"), false),
        (&same[..0], false),
    ] {
        let owned = other.to_owned();
        // Exercise each comparison implemented for SeqSlice and &SeqSlice.
        assert_eq!(*slice == *other, equal);
        assert_eq!(slice == *other, equal);
        assert_eq!(*slice == owned, equal);
        assert_eq!(slice == owned, equal);
    }
}

#[test]
fn owning_a_slice_copies_its_symbols_independently() {
    let capacity = usize::BITS as usize / usize::from(Dna::BITS);
    let text = "TACG".repeat(capacity / 4 + 2);
    let seq: Seq<Dna> = text.parse().unwrap();
    let slice = &seq[1..capacity + 2];
    let expected = &text[1..capacity + 2];
    let mut owned = slice.to_owned();

    assert_eq!(owned.to_string(), expected);
    owned.comp();
    assert_ne!(owned.as_ref(), slice);
    assert_eq!(String::from(slice), expected);
    assert_eq!(seq.to_string(), text);
}

#[test]
fn iupac_set_operations_on_offset_slices() {
    let left = &iupac!("NAS-GYTNAN-")[1..10];
    let right = &iupac!("NNANTGCAT-NC")[2..11];

    assert_eq!(left & right, iupac!("AS-GC-T-N"));
    assert_eq!(left | right, iupac!("ANTGYWNAN"));
    assert_eq!(left.to_string(), "AS-GYTNAN");
    assert_eq!(right.to_string(), "ANTGCAT-N");

    let empty = &left[..0];
    assert!((empty & empty).is_empty());
    assert!((empty | empty).is_empty());
}

#[test]
fn display_propagates_writer_errors() {
    struct RejectWrites;

    impl fmt::Write for RejectWrites {
        fn write_str(&mut self, _: &str) -> fmt::Result {
            Err(fmt::Error)
        }
    }

    let slice = &dna!("TACGTG")[1..5];
    let mut writer = RejectWrites;
    assert_eq!(write!(writer, "{slice}"), Err(fmt::Error));
}

#[derive(Default)]
struct RecordingHasher(Vec<u8>);

impl Hasher for RecordingHasher {
    fn finish(&self) -> u64 {
        0
    }

    fn write(&mut self, bytes: &[u8]) {
        self.0.extend_from_slice(bytes);
    }
}

#[test]
fn hash_bytes_ignore_surrounding_symbols_and_alignment() {
    let first = &dna!("TACGTATT")[1..6];
    let second = &dna!("CCACGTAGG")[2..7];

    // Eight-byte length prefix, followed by ACGT and a partial byte containing A.
    let expected = [5, 0, 0, 0, 0, 0, 0, 0, 0xe4, 0];
    for slice in [first, second] {
        let mut hasher = RecordingHasher::default();
        slice.hash(&mut hasher);
        assert_eq!(hasher.0, expected);
    }
}
