use bitvec::prelude::*;

pub struct Kmer<const K: u8> {
    bv: BitArray,
}
