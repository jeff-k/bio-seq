// Copyright 2021, 2022, 2025 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
//use crate::kmer::KmerIter;
use crate::seq::storage::{SeqSliceStorage, SeqStorage};
use crate::seq::{Seq, SeqSlice};
use core::iter::Chain;
use core::marker::PhantomData;

/// An iterator over fixed-size non-overlapping chunks of a sequence
pub struct SeqChunks<'a, A: Codec, S: SeqStorage> {
    slice: &'a SeqSlice<A, S>,
    width: usize,
    skip: usize,
    index: usize,
}

/// An iterator over the elements of a sequence
pub struct SeqIter<'a, A: Codec, S: SeqStorage> {
    slice: &'a SeqSlice<A, S>,
    index: usize,
}

impl<'a, A: Codec, S: SeqStorage> SeqSlice<A, S> {
    pub fn chain(
        self: &'a SeqSlice<A, S>,
        second: &'a SeqSlice<A, S>,
    ) -> Chain<SeqIter<'a, A, S>, SeqIter<'a, A, S>> {
        self.into_iter().chain(second)
    }

    pub fn iter(&'a self) -> SeqIter<'a, A, S> {
        <&Self as IntoIterator>::into_iter(self)
    }

    /*
        /// Iterate over sliding windows of length K
        pub fn kmers<const K: usize>(&self) -> KmerIter<'_, A, K> {
            KmerIter::<A, K> {
                slice: self,
                index: 0,
                len: self.len(),
                _p: PhantomData,
            }
        }
    */

    /// Iterate over the sequence in reverse order
    pub fn rev_iter(&self) -> RevIter<'_, A, S> {
        RevIter {
            slice: self,
            index: self.len(),
        }
    }

    /// Iterate over the sequence in overlapping windows of a specified width
    ///
    /// ```
    /// use bio_seq::prelude::*;
    ///
    /// let seq: Seq<Dna> = "ACTGATCG".try_into().unwrap();
    /// let windows: Vec<String> = seq.windows(3).map(String::from).collect();
    /// assert_eq!(windows, vec!["ACT", "CTG", "TGA", "GAT", "ATC", "TCG"]);
    /// ```
    pub fn windows(&self, width: usize) -> SeqChunks<'_, A, S> {
        SeqChunks {
            slice: self,
            width,
            skip: 1,
            index: 0,
        }
    }

    /// Iterate over the sequence in non-overlapping chunks of a specified width
    ///
    /// The last incomplete chunk will be excluded if the sequence length is not divisible by the specified
    /// width.
    ///
    /// ```
    /// use bio_seq::prelude::*;
    ///
    /// let seq: Seq<Dna> = "ACTGATCG".try_into().unwrap();
    /// let chunks: Vec<Seq<Dna>> = seq.chunks(3).collect();
    /// assert_eq!(chunks, vec![dna!("ACT"), dna!("GAT")]);
    /// ```
    pub fn chunks(&self, width: usize) -> SeqChunks<'_, A, S> {
        SeqChunks {
            slice: self,
            width,
            skip: width,
            index: 0,
        }
    }
}

/// An iterator over the elements of a sequence in reverse order
pub struct RevIter<'a, A: Codec, S: SeqStorage> {
    pub slice: &'a SeqSlice<A, S>,
    pub index: usize,
}

impl<A: Codec, S: SeqStorage> Iterator for RevIter<'_, A, S> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let i = self.index;

        if self.index == 0 {
            return None;
        }
        self.index -= 1;
        Some(self.slice.nth(i - 1))
    }
}

impl<'a, A: Codec + core::fmt::Debug, S: SeqStorage> Iterator for SeqChunks<'a, A, S> {
    type Item = &'a SeqSlice<A, S>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index + self.width > self.slice.len() {
            return None;
        }
        let i = self.index;
        self.index += self.skip;
        if i + self.width > self.slice.len() {
            return None;
        }
        Some(&self.slice[i..i + self.width])
    }
}

impl<'a, A: Codec, S: SeqStorage> IntoIterator for &'a Seq<A, S> {
    type Item = A;
    type IntoIter = SeqIter<'a, A, S>;

    fn into_iter(self) -> Self::IntoIter {
        self.as_ref().into_iter()
    }
}

impl<'a, A: Codec, S: SeqStorage> IntoIterator for &'a SeqSlice<A, S> {
    type Item = A;
    type IntoIter = SeqIter<'a, A, S>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter {
            slice: self,
            index: 0,
        }
    }
}

impl<'a, A: Codec, S: SeqStorage> FromIterator<&'a SeqSlice<A, S>> for Vec<Seq<A, S>> {
    fn from_iter<T: IntoIterator<Item = &'a SeqSlice<A, S>>>(iter: T) -> Self {
        iter.into_iter().collect()
    }
}

impl<A: Codec, S: SeqStorage> Iterator for SeqIter<'_, A, S> {
    type Item = A;
    fn next(&mut self) -> Option<A> {
        let i = self.index;
        if self.index >= self.slice.len() {
            return None;
        }
        self.index += 1;
        Some(self.slice.nth(i))
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::dna::Dna::*;
    use crate::prelude::*;

    #[test]
    fn seq_iter() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC").into();
        let elements: Vec<Dna> = seq.into_iter().collect();
        assert_eq!(elements, vec![A, C, T, G, A, T, C, G, A, T, A, C]);
        assert_ne!(elements, vec![A, C, T, G, A, T, C, G, A, T, A, C, A]);
        assert_ne!(elements, vec![C, A, T, A, G, C, T, A, G, T, C, A]);
    }

    #[test]
    fn rev_iter() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC").into();
        let rev_elements: Vec<Dna> = seq.rev_iter().collect();
        assert_ne!(rev_elements, vec![A, C, T, G, A, T, C, G, A, T, A, C]);
        assert_eq!(rev_elements, vec![C, A, T, A, G, C, T, A, G, T, C, A]);
    }

    #[test]
    fn iterators() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC").into();
        let slice = &seq[2..9];
        let elements: Vec<Dna> = slice.into_iter().collect();
        assert_eq!(elements, vec![T, G, A, T, C, G, A]);
    }

    #[test]
    fn chunks() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC").into();
        let cs: Vec<Seq<Dna>> = seq.chunks(5).collect();
        assert_eq!(cs[0], dna!("ACTGA"));
        assert_eq!(cs[1], dna!("TCGAT"));
        assert_eq!(cs.len(), 2);
    }

    #[test]
    fn test_chain() {
        let seq1 = Seq::<Dna>::try_from("ATG").unwrap();
        let seq2 = Seq::<Dna>::try_from("TAC").unwrap();

        let chained = seq1.chain(&seq2);

        let expected_seq = Seq::<Dna>::try_from("ATGTAC").unwrap();
        for (a, b) in chained.zip(expected_seq.into_iter()) {
            assert_eq!(a, b);
        }

        let chained = seq1.chain(&seq2);
        for (a, b) in chained.map(|b| b.to_comp()).zip(expected_seq.into_iter()) {
            assert_ne!(a, b);
        }
    }

    #[test]
    fn windows() {
        let seq: Seq<Dna> = dna!("ACTGATACG").into();
        let windows: Vec<String> = seq.windows(5).map(String::from).collect();
        assert_eq!(windows, vec!["ACTGA", "CTGAT", "TGATA", "GATAC", "ATACG"]);
    }
}
