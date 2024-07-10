// Copyright 2023-2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! `bio-seq-derive` is a procedural macro crate that provides the `Codec` derive macro for the `bio-seq` library.
//! It allows users to define custom bit-packed encodings from an enum. The representation of the enum is derived from the discriminants.
//! Please refer to the `bio-seq` [documentation](https://github.com/jeff-k/bio-seq) for a complete guide on defining custom alphabets.

use quote::quote;
use syn::LitStr;

pub(crate) fn gen_seqarray(
    codec: &syn::Ident,
    seq_name: &str,
    len: usize,
    bits: &[u8],
) -> proc_macro2::TokenStream {
    let num_words: usize = (bits.len() + (usize::BITS as usize - 1)) / usize::BITS as usize;

    let seq_name_ident = syn::Ident::new(seq_name, proc_macro2::Span::call_site());

    quote! {
    {
    type Lsb0 = __bio_seq_Lsb0;
    static #seq_name_ident: SeqArray<#codec, #len, #num_words> = SeqArray {
        _p: core::marker::PhantomData,
        ba: __bio_seq_bitarr![const usize, Lsb0; #(#bits),*]
    };
    & #seq_name_ident
    }
    }
}

pub(crate) fn dna_seq(seq: &LitStr) -> Result<(usize, Vec<u8>), syn::Error> {
    let mut bits: Vec<u8> = Vec::new();
    let mut bases: usize = 0;

    for (i, c) in seq.value().char_indices() {
        match c {
            'A' => bits.extend([0, 0]),
            'C' => bits.extend([1, 0]),
            'G' => bits.extend([0, 1]),
            'T' => bits.extend([1, 1]),
            _ => {
                return Err(syn::Error::new_spanned(
                    seq.value(),
                    format!("Invalid DNA base at position {i}: {c}"),
                ))
            }
        }
        bases += 1;
    }

    Ok((bases, bits))
}

pub(crate) fn iupac_seq(seq: &LitStr) -> Result<(usize, Vec<u8>), syn::Error> {
    let mut bits: Vec<u8> = Vec::new();
    let mut bases: usize = 0;

    for (i, c) in seq.value().char_indices() {
        match c {
            'A' => bits.extend([0, 0, 0, 1]),
            'C' => bits.extend([0, 0, 1, 0]),
            'G' => bits.extend([0, 1, 0, 0]),
            'T' => bits.extend([1, 0, 0, 0]),

            'R' => bits.extend([0, 1, 0, 1]),
            'Y' => bits.extend([1, 0, 1, 0]),
            'S' => bits.extend([0, 1, 1, 0]),
            'W' => bits.extend([1, 0, 0, 1]),
            'K' => bits.extend([1, 1, 0, 0]),
            'M' => bits.extend([0, 0, 1, 1]),

            'B' => bits.extend([1, 1, 1, 0]),
            'D' => bits.extend([1, 1, 0, 1]),
            'H' => bits.extend([1, 0, 1, 1]),
            'V' => bits.extend([0, 1, 1, 1]),

            'N' => bits.extend([1, 1, 1, 1]),
            'X' | '-' => bits.extend([0, 0, 0, 0]),
            _ => {
                return Err(syn::Error::new_spanned(
                    seq.value(),
                    format!("Invalid IUPAC code at position {i}: {c}"),
                ))
            }
        }
        bases += 1;
    }

    Ok((bases, bits))
}
