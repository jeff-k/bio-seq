// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::seq::SeqSlice;
use crate::Kmer;
use core::marker::PhantomData;

pub struct SeqChunks<'a, A: Codec> {
    slice: &'a SeqSlice<A>,
    width: usize,
    skip: usize,
    index: usize,
}

pub struct SeqIter<'a, A: Codec> {
    slice: &'a SeqSlice<A>,
    index: usize,
}

impl<A: Codec> SeqSlice<A> {
    /// Iterate over sliding windows of size K
    pub fn kmers<const K: usize>(&self) -> KmerIter<A, K> {
        KmerIter::<A, K> {
            slice: self,
            index: 0,
            len: self.len(),
            _p: PhantomData,
        }
    }

    /// Iterate over the sequence in reverse order
    pub fn rev(&self) -> RevIter<A> {
        RevIter {
            slice: self,
            index: self.len(),
        }
    }

    pub fn windows(&self, width: usize) -> SeqChunks<A> {
        SeqChunks {
            slice: self,
            width,
            skip: 1,
            index: 0,
        }
    }

    pub fn chunks(&self, width: usize) -> SeqChunks<A> {
        SeqChunks {
            slice: self,
            width,
            skip: width,
            index: 0,
        }
    }
}

pub struct RevIter<'a, A: Codec> {
    pub slice: &'a SeqSlice<A>,
    pub index: usize,
}

impl<'a, A: Codec> Iterator for RevIter<'a, A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let i = self.index;
        if self.index == 0 {
            return None;
        }
        self.index -= 1;
        Some(A::unsafe_from_bits(self.slice[i].into()))
    }
}

impl<'a, A: Codec> Iterator for SeqChunks<'a, A> {
    type Item = &'a SeqSlice<A>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.width >= self.slice.len() {
            return None;
        }
        let i = self.index;
        self.index += self.width + self.skip;
        Some(&self.slice[i..i + self.width])
    }
}

impl<'a, A: Codec> IntoIterator for &'a SeqSlice<A> {
    type Item = A;
    type IntoIter = SeqIter<'a, A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter {
            slice: self,
            index: 0,
        }
    }
}

// TODO
// IntoIter for Seq should box the contained slice
// IntoIter for &'a Seq

impl<'a, A: Codec> Iterator for SeqIter<'a, A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let i = self.index;
        if self.index >= (self.slice.len()) {
            return None;
        }
        self.index += 1;
        Some(A::unsafe_from_bits(self.slice[i].into()))
    }
}

pub struct KmerIter<'a, A: Codec, const K: usize> {
    pub slice: &'a SeqSlice<A>,
    pub index: usize,
    pub len: usize,
    pub _p: PhantomData<A>,
}

impl<'a, A: Codec, const K: usize> Iterator for KmerIter<'a, A, K> {
    type Item = Kmer<A, K>;
    fn next(&mut self) -> Option<Kmer<A, K>> {
        let i = self.index;
        if self.index >= self.len - (K - 1) {
            return None;
        }
        self.index += 1;
        Some(Kmer::<A, K>::from(&self.slice[i..K + i]))
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::*;
    use crate::seq::{FromStr, Seq, SeqSlice};

    /*
    #[test]
    fn chunks() {
        let cs: Vec<&SeqSlice<Dna>> = dna!("ACTGATACGTA").chunks(5).collect();
        assert_eq!(format!("{}", cs[0]), "ACTGA");
    }

    #[test]
    fn windows() {
        let cs: Vec<Seq<Dna>> = dna!("ACTGATACGTA").windows(5).collect();
        assert_eq!(format!("{}", cs[0]), "ACTGA");
    }
    */
}
