// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::seq::SeqSlice;
use crate::Kmer;
use std::marker::PhantomData;

pub struct SeqChunks<'a, A: Codec> {
    slice: &'a SeqSlice<A>,
    width: usize,
    skip: usize,
    index: usize,
}

impl<A: Codec> SeqSlice<A> {
    /// Iterate over sliding windows of size K
    pub fn kmers<'a, const K: usize>(&'a self) -> KmerIter<'a, A, K> {
        KmerIter::<A, K> {
            slice: self,
            index: 0,
            len: self.len(),
            _p: PhantomData,
        }
    }

    /*
        /// Iterate over the sequence in reverse order
        pub fn rev(self) -> RevIter<A> {
            let index = self.bs.len();
            RevIter::<A> { slice: self, index }
        }
    */
    pub fn windows<'a>(&'a self, width: usize) -> SeqChunks<'a, A> {
        SeqChunks {
            slice: self,
            width,
            skip: 1,
            index: 0,
        }
    }

    pub fn chunks<'a>(&'a self, width: usize) -> SeqChunks<'a, A> {
        SeqChunks {
            slice: self,
            width,
            skip: width,
            index: 0,
        }
    }
}

/*
pub struct RevIter<'a, A: Codec> {
    pub slice: &'a SeqSlice<A>,
    pub index: usize,
}

impl<'a, A: Codec> Iterator for RevIter<'a, A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH as usize;
        if self.index == 0 {
            return None;
        }

        self.index -= w;
        let i = self.index;
        Some(A::unsafe_from_bits(self.slice[i]))
    }
}
*/

impl<'a, A: Codec> Iterator for SeqChunks<'a, A> {
    type Item = &'a SeqSlice<A>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.width + self.skip >= self.slice.len() {
            return None;
        }
        let i = self.index;
        self.index += self.width + self.skip;
        Some(&self.slice[i..self.width])
    }
}

/*
impl<A: Codec> IntoIterator for SeqSlice<A> {
    type Item = A;
    type IntoIter = SeqIter<A>;

    fn into_iter(&self) -> Self::IntoIter {
        SeqIter::<A> {
            slice: Box::new(self),
            index: 0,
        }
    }
}

pub struct SeqIter<A: Codec> {
    slice: Box<SeqSlice<A>>,
    index: usize,
}

impl<A: Codec> Iterator for SeqIter<A> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let w = A::WIDTH as usize;
        let i = self.index;
        if self.index >= (self.slice.len()) {
            return None;
        }
        self.index += w;
        Some(A::unsafe_from_bits(self.slice[i]))
    }
}
*/

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
