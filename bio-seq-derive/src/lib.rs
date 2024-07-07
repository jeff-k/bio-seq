// Copyright 2023 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! `bio-seq-derive` is a procedural macro crate that provides the `Codec` derive macro for the `bio-seq` library.
//! It allows users to define custom bit-packed encodings from an enum. The representation of the enum is derived from the discriminants.
//! Please refer to the `bio-seq` [documentation](https://github.com/jeff-k/bio-seq) for a complete guide on defining custom alphabets.

#![warn(clippy::pedantic)]

mod codec;
mod seqarray;

use crate::codec::{parse_variants, parse_width, test_repr, CodecVariants};
use crate::seqarray::{dna_seq, gen_seqarray, iupac_seq};
use quote::quote;
use std::hash::{DefaultHasher, Hash, Hasher};
use syn::parse_macro_input;
use syn::LitStr;

#[proc_macro_derive(Codec, attributes(bits, display, alt))]
pub fn codec_derive(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input = parse_macro_input!(input as syn::Item);

    // Test for correct usage
    let syn::Item::Enum(enum_ast) = input else {
        return syn::Error::new_spanned(input, "Codec can only be derived for enums")
            .to_compile_error()
            .into();
    };

    // Test whether enum is #[repr(u8)]
    let _is_repr8 = test_repr(&enum_ast);

    let variants = match parse_variants(&enum_ast.variants) {
        Ok(variants) => variants,
        Err(err) => return err.to_compile_error().into(),
    };

    let enum_ident = enum_ast.ident;

    let CodecVariants {
        idents,
        to_chars,
        from_chars,
        unsafe_alts,
        alts,
        max_discriminant,
    } = variants;

    let width = match parse_width(&enum_ast.attrs, max_discriminant) {
        Ok(width) => width,
        Err(err) => return err.into_compile_error().into(),
    };

    // Generate the implementation
    let output = quote! {
        impl Codec for #enum_ident {
            const BITS: u8 = #width;

            fn unsafe_from_bits(b: u8) -> Self {
                match b {
                    #(#unsafe_alts),*,
                    _ => panic!("Unrecognised bit pattern: {:08b}", b),
                }
            }

            fn try_from_bits(b: u8) -> Option<Self> {
                match b {
                    #(#alts),*,
                    _ => None,
                }
            }

            fn unsafe_from_ascii(c: u8) -> Self {
                match c {
                    #(#from_chars),*,
                    _ => {
                        if c.is_ascii_alphanumeric() {
                            panic!("Unrecognised character: {} ({:#04X?})", c as char, c);
                        } else {
                            panic!("Unrecognised character: {:#04X?}", c);
                        }
                    },
                }.unwrap()
            }

            fn try_from_ascii(c: u8) -> Option<Self> {
                match c {
                    #(#from_chars),*,
                    _ => None,
                }
            }

            fn to_char(self) -> char {
                match self {
                    #(#to_chars),*,
                }.into()
            }

            fn to_bits(self) -> u8 {
                self as u8
            }

            fn items() -> impl Iterator<Item = Self> {
                vec![ #(Self::#idents,)* ].into_iter()
            }
        }

    };
    output.into()
}

#[proc_macro]
pub fn dna(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let seq: LitStr = parse_macro_input!(input as LitStr);

    if !seq.value().is_ascii() {
        return syn::Error::new_spanned(seq, "Non-ASCII characters in DNA string")
            .to_compile_error()
            .into();
    }

    let seq_name = {
        let mut hasher = DefaultHasher::new();
        seq.value().hash(&mut hasher);
        format!("DNA_SEQ_{:0X}", hasher.finish())
    };

    match dna_seq(&seq) {
        Ok((len, bits)) => {
            let encoding: syn::Ident = syn::Ident::new("Dna", proc_macro2::Span::call_site());
            gen_seqarray(&encoding, &seq_name, len, &bits).into()
        }
        Err(e) => e.to_compile_error().into(),
    }
}

#[proc_macro]
pub fn iupac(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let seq: LitStr = parse_macro_input!(input as LitStr);

    if !seq.value().is_ascii() {
        return syn::Error::new_spanned(seq, "Non-ASCII characters in IUPAC string")
            .to_compile_error()
            .into();
    }

    let seq_name = {
        let mut hasher = DefaultHasher::new();
        seq.value().hash(&mut hasher);
        format!("IUPAC_SEQ_{:0X}", hasher.finish())
    };

    match iupac_seq(&seq) {
        Ok((len, bits)) => {
            let encoding: syn::Ident = syn::Ident::new("Iupac", proc_macro2::Span::call_site());
            gen_seqarray(&encoding, &seq_name, len, &bits).into()
        }
        Err(e) => e.to_compile_error().into(),
    }
}
