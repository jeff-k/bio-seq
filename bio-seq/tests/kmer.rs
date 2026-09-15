use bio_seq::codec::degenerate::WS;
use bio_seq::prelude::*;
/*
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Codec)]
#[repr(u8)]
enum Sparse {
    A = 1,
    C = 2,
}
*/

#[test]
fn reverse_iupac_kmer() {
    let k: Kmer<Iupac, 2> = "AC".parse().unwrap();
    assert_eq!(k.to_rev().to_string(), "CA");
}

#[cfg(feature = "extra_codecs")]
#[test]
fn reverse_ws_kmer() {
    let k: Kmer<WS, 2> = "WS".parse().unwrap();
    assert_eq!(k.to_rev().to_string(), "SW");
}
