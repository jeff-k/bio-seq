// Copyright 2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! SIMD kmers
//!

#![feature(portable_simd)]

use crate::codec::Codec;
use crate::kmer::{Kmer, KmerStorage};
use crate::prelude::{Complement, ParseBioError, ReverseComplement};
use crate::seq::{Seq, SeqArray, SeqSlice};
use crate::{Ba, Bs, Bv};

use bitvec::field::BitField;
use bitvec::view::BitView;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::ops::Deref;
use core::ptr;
use core::str::FromStr;
use std::simd;

impl KmerStorage for simd::Simd<u64, 4> {
    const BITS: usize = todo!();

    fn new() -> Self {
        todo!()
    }
}

impl<A: Codec, const K: usize> Kmer<A, K, simd::Simd<u64, 4>> {
    fn from_seq(seq: &SeqSlice<A>) -> Self {
        todo!()
    }
}

impl<A: Codec, const K: usize> Hash for Kmer<A, K, simd::Simd<u64, 4>> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        todo!()
    }
}

impl<A: Codec, const K: usize> TryFrom<&SeqSlice<A>> for Kmer<A, K, simd::Simd<u64, 4>> {
    type Error = ParseBioError;

    fn try_from(seq: &SeqSlice<A>) -> Result<Self, Self::Error> {
        todo!()
    }
}
