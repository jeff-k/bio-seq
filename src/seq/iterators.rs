// Copyright 2021, 2022 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::codec::Codec;
use crate::kmer::KmerIter;
use crate::seq::{Seq, SeqSlice};
use core::iter::Chain;
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

impl<'a, A: Codec> SeqSlice<A> {
    pub fn chain(
        self: &'a SeqSlice<A>,
        second: &'a SeqSlice<A>,
    ) -> Chain<SeqIter<'a, A>, SeqIter<'a, A>> {
        self.into_iter().chain(second.into_iter())
    }
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
    /// The last incomplete chunk will be excluded if the sequence length is not divisible by the specified
    /// width.
    ///
    /// Example:
    ///
    /// ```
    /// use bio_seq::prelude::*;
    ///
    /// let seq: Seq<Dna> = "ACTGATCG".try_into().unwrap();
    /// let chunks: Vec<Seq<Dna>> = seq.chunks(3).collect();
    /// assert_eq!(chunks, vec!["ACT", "GAT"]);
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
        if i + self.width > self.slice.len() {
            return None;
        }
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

impl<'a, A: Codec> FromIterator<&'a SeqSlice<A>> for Vec<Seq<A>> {
    fn from_iter<T: IntoIterator<Item = &'a SeqSlice<A>>>(iter: T) -> Self {
        iter.into_iter().map(|slice| slice.to_owned()).collect()
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

#[cfg(test)]
mod tests {
    use crate::codec::dna::{Dna, Dna::*};
    use crate::codec::Complement;
    use crate::seq::{FromStr, Seq};

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
        for (a, b) in chained.map(|b| b.comp()).zip(expected_seq.into_iter()) {
            assert_ne!(a, b);
        }
    }

    #[test]
    fn windows() {
        let seq: Seq<Dna> = dna!("ACTGATACG");
        let windows: Vec<Seq<Dna>> = seq.windows(5).collect();
        assert_eq!(windows, vec!["ACTGA", "CTGAT", "TGATA", "GATAC", "ATACG"]);
    }
}
