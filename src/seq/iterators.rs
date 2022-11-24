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
        Some(A::unsafe_from_bits(self.slice[i - 1].into()))
    }
}

impl<'a, A: Codec + std::fmt::Debug> Iterator for SeqChunks<'a, A> {
    type Item = &'a SeqSlice<A>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.width > self.slice.len() {
            return None;
        }
        let i = self.index;
        self.index += self.skip;
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
        if self.index >= self.slice.len() {
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
        if self.index + K > self.len {
            return None;
        }
        self.index += 1;
        Some(Kmer::<A, K>::from(&self.slice[i..i + K]))
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::*;
    use crate::kmer::Kmer;
    use crate::seq::{FromStr, Seq, SeqSlice};

    #[test]
    fn chunks() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC");
        let cs: Vec<&SeqSlice<Dna>> = seq.chunks(5).collect();
        assert_eq!(format!("{}", cs[0]), "ACTGA");
        assert_eq!(format!("{}", cs[1]), "TCGAT");
        assert_eq!(cs.len(), 2);
    }

    #[test]
    fn kmer_iter() {
        let seq: Seq<Dna> = dna!("ACTGA");
        let cs: Vec<Kmer<Dna, 3>> = seq.kmers().collect();
        assert_eq!(format!("{}", cs[0]), "ACT");
        assert_eq!(format!("{}", cs[1]), "CTG");
        assert_eq!(format!("{}", cs[2]), "TGA");
        assert_eq!(cs.len(), 3);
    }

    #[test]
    fn windows() {
        let seq: Seq<Dna> = dna!("ACTGATACG");
        let cs: Vec<&SeqSlice<Dna>> = seq.windows(5).collect();
        assert_eq!(format!("{}", cs[0]), "ACTGA");
        assert_eq!(format!("{}", cs[1]), "CTGAT");
        assert_eq!(format!("{}", cs[2]), "TGATA");
        assert_eq!(format!("{}", cs[3]), "GATAC");
        assert_eq!(format!("{}", cs[4]), "ATACG");
        assert_eq!(cs.len(), 5);
    }
}
