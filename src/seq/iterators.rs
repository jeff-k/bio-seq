// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::seq::{Seq, SeqSlice};
use crate::Kmer;
use bitvec::prelude::*;
use std::marker::PhantomData;

pub struct RevIter<A: Codec> {
    pub seq: Seq<A>,
    pub index: usize,
}

impl<A: Codec> Iterator for RevIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH as usize;
        if self.index == 0 {
            return None;
        }

        self.index -= w;
        let i = self.index;
        Some(A::unsafe_from_bits(self.seq.bv[i..i + w].load()))
        //Some(self.seq.bv[i..i+w].as_ref())
    }
}

pub struct SeqChunks<A: Codec> {
    pub seq: SeqSlice<A>,
    pub width: usize,
    pub skip: usize,
    pub index: usize,
}

impl<A: Codec> Iterator for SeqChunks<A> {
    type Item = SeqSlice<A>;
    fn next(&mut self) -> Option<Self::Item> {
        let w = A::WIDTH as usize * self.width;

        if self.index + w >= self.seq.bs.len() {
            return None;
        }
        let i = self.index;
        self.index += w + self.skip;
        Some(SeqSlice {
            bs: BitBox::from_bitslice(&self.seq.bs[i..i + w]),
            _p: self.seq._p,
        })
    }
}

impl<A: Codec> IntoIterator for Seq<A> {
    type Item = A;
    type IntoIter = SeqIter<A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter::<A> {
            seq: SeqSlice {
                bs: BitBox::from_bitslice(&self.bv),
                _p: self._p,
            },
            index: 0,
        }
    }
}

pub struct SeqIter<A: Codec> {
    seq: SeqSlice<A>,
    index: usize,
}

impl<A: Codec> Iterator for SeqIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH as usize;
        let i = self.index;
        if self.index >= (self.seq.bs.len()) {
            return None;
        }
        self.index += w;
        Some(A::unsafe_from_bits(self.seq.bs[i..i + w].load()))
    }
}

pub struct KmerIter<A: Codec, const K: usize> {
    pub bs: BitBox,
    pub index: usize,
    pub len: usize,
    pub _p: PhantomData<A>,
}

impl<A: Codec, const K: usize> Iterator for KmerIter<A, K> {
    type Item = Kmer<A, K>;
    fn next(&mut self) -> Option<Kmer<A, K>> {
        let k = K * A::WIDTH as usize;
        let i = self.index * A::WIDTH as usize;
        if self.index >= self.len - (K - 1) {
            return None;
        }
        self.index += 1;
        Some(Kmer::<A, K>::new(&self.bs[i..k + i]))
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::*;
    use crate::seq::{FromStr, Seq, SeqSlice};

    #[test]
    fn chunks() {
        let cs: Vec<SeqSlice<Dna>> = dna!("ACTGATACGTA").chunks(5).collect();
        assert_eq!(format!("{}", cs[0]), "ACTGA");
    }

    #[test]
    fn windows() {
        let cs: Vec<SeqSlice<Dna>> = dna!("ACTGATACGTA").windows(5).collect();
        assert_eq!(format!("{}", cs[0]), "ACTGA");
    }
}
