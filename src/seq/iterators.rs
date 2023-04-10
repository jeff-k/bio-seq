// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::kmer::Kmer;
use crate::seq::{Seq, SeqSlice};
use core::marker::PhantomData;

/// An iterator over fixed-size non-overlapping chunks of a sequence
pub struct SeqChunks<'a, A: Codec> {
    slice: &'a SeqSlice<A>,
    width: usize,
    skip: usize,
    index: usize,
}

/// An iterator over the elements of a sequence
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

    /// Iterate over the sequence in overlapping windows of a specified width
    ///
    /// This will panic if the length of the window is greater than the length of the sequence.
    ///
    /// Example:
    ///
    /// ```
    /// use bio_seq::prelude::*;
    ///
    /// let seq: Seq<Dna> = "ACTGATCG".try_into().unwrap();
    /// let windows: Vec<Seq<Dna>> = seq.windows(3).collect();
    /// assert_eq!(windows, vec!["ACT", "CTG", "TGA", "GAT", "ATC", "TCG"]);
    /// ```
    pub fn windows(&self, width: usize) -> SeqChunks<A> {
        SeqChunks {
            slice: self,
            width,
            skip: 1,
            index: 0,
        }
    }

    /// Iterate over the sequence in non-overlapping chunks of a specified width
    ///
    /// The last chunk may be smaller if the sequence length is not divisible by the specified
    /// width.
    ///
    /// Example:
    ///
    /// ```
    /// use bio_seq::prelude::*;
    ///
    /// let seq: Seq<Dna> = "ACTGATCG".try_into().unwrap();
    /// let chunks: Vec<Seq<Dna>> = seq.chunks(3).collect();
    /// assert_eq!(chunks, vec!["ACT", "GAT", "CG"]);
    /// ```
    pub fn chunks(&self, width: usize) -> SeqChunks<A> {
        SeqChunks {
            slice: self,
            width,
            skip: width,
            index: 0,
        }
    }
}

/// An iterator over the elements of a sequence in reverse order
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

impl<'a, A: Codec> IntoIterator for &'a Seq<A> {
    type Item = A;
    type IntoIter = SeqIter<'a, A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter {
            slice: self,
            index: 0,
        }
    }
}

impl<'a, A: Codec> IntoIterator for &'a &Seq<A> {
    type Item = A;
    type IntoIter = SeqIter<'a, A>;

    fn into_iter(self) -> Self::IntoIter {
        SeqIter {
            slice: self,
            index: 0,
        }
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

/// An iterator over all kmers of a sequence with a specified length
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
    use crate::codec::dna::{Dna, Dna::*};
    use crate::kmer::Kmer;
    use crate::seq::{FromStr, Seq, SeqSlice};

    #[test]
    fn seq_iter() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC");
        let elements: Vec<Dna> = seq.into_iter().collect();
        assert_eq!(elements, vec![A, C, T, G, A, T, C, G, A, T, A, C]);
        assert_ne!(elements, vec![A, C, T, G, A, T, C, G, A, T, A, C, A]);
        assert_ne!(elements, vec![C, A, T, A, G, C, T, A, G, T, C, A]);
    }

    #[test]
    fn rev_iter() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC");
        let rev_elements: Vec<Dna> = seq.rev().collect();
        assert_ne!(rev_elements, vec![A, C, T, G, A, T, C, G, A, T, A, C]);
        assert_eq!(rev_elements, vec![C, A, T, A, G, C, T, A, G, T, C, A]);
    }

    #[test]
    fn iterators() {
        let seq: Seq<Dna> = dna!("ACTGATCGATAC");
        let slice = &seq[2..9];
        let elements: Vec<Dna> = slice.into_iter().collect();
        assert_eq!(elements, vec![T, G, A, T, C, G, A]);
    }

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
