// Copyright 2023 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! `bio-seq-derive` is a procedural macro crate that provides the `Codec` derive macro for the `bio-seq` library.
//! It allows users to define custom bit-packed encodings from an enum. The representation of the enum is derived from the discriminants.
//! Please refer to the `bio-seq` [documentation](https://github.com/jeff-k/bio-seq) for a complete guide on defining custom alphabets.

extern crate proc_macro;

use crate::proc_macro::TokenStream;

use quote::quote;

use syn::punctuated::Punctuated;
use syn::{parse_macro_input, Token};

#[proc_macro_derive(Codec, attributes(width, display, alt))]
pub fn codec_derive(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as syn::Item);

    // Test for correct usage
    let enum_ast = match input {
        syn::Item::Enum(e) => e,
        _ => {
            return syn::Error::new_spanned(input, "Codec can only be derived for enums")
                .to_compile_error()
                .into()
        }
    };

    // Test that enum is repr(u8)
    let test_repr_u8 = enum_ast.attrs.iter().any(|attr| {
        attr.path().is_ident("repr")
            && match attr.parse_args::<syn::Ident>() {
                Ok(ident) => ident == "u8",
                Err(_) => false,
            }
    });

    if !test_repr_u8 {
        return syn::Error::new_spanned(
            &enum_ast.ident,
            "Enums deriving Codec must be annotated with #[repr(u8)]",
        )
        .to_compile_error()
        .into();
    }

    let variants = enum_ast.variants;
    let enum_ident = enum_ast.ident;
    let mut max_variant = 0u8;
    let mut variant_idents = Vec::new();

    let mut variants_to_char = Vec::new();
    let mut chars_to_variant = Vec::new();
    let mut alt_discriminants = Vec::new();
    let mut unsafe_alts = Vec::new();

    for variant in variants.iter() {
        let ident = &variant.ident;
        variant_idents.push(ident.clone());
        let discriminant = &variant.discriminant;

        if let Some((_, syn::Expr::Lit(expr_lit))) = discriminant {
            let value = match &expr_lit.lit {
                // discriminants must be either integers or byte literals
                syn::Lit::Byte(lit_byte) => lit_byte.value(),
                syn::Lit::Int(lit_int) => lit_int.base10_parse::<u8>().unwrap(),
                _ => {
                    return syn::Error::new_spanned(
                        ident,
                        "Codec derivations require byte or integer discriminants",
                    )
                    .to_compile_error()
                    .into();
                }
            };

            alt_discriminants.push(quote! { #value => Ok(Self::#ident) });
            unsafe_alts.push(quote! { #value => Self::#ident });

            max_variant = max_variant.max(value);
        } else {
            return syn::Error::new_spanned(ident, "Codec derivations require discriminants")
                .to_compile_error()
                .into();
        }

        let mut char_repr = ident.to_string().chars().next().unwrap();

        for attr in &variant.attrs {
            if attr.path().is_ident("display") {
                let alt_attr: syn::LitChar = match attr.parse_args() {
                    Ok(attr) => attr,
                    Err(err) => return err.to_compile_error().into(),
                };
                char_repr = alt_attr.value();
            } else if attr.path().is_ident("alt") {
                let discs: Punctuated<syn::ExprLit, Token![,]> =
                    match attr.parse_args_with(Punctuated::parse_terminated) {
                        Ok(discs) => discs,
                        Err(err) => return err.to_compile_error().into(),
                    };

                for d in discs.into_iter() {
                    alt_discriminants.push(quote! { #d => Ok(Self::#ident) });
                    unsafe_alts.push(quote! { #d => Self::#ident });
                }
            };
        }

        variants_to_char.push(quote! { Self::#ident => #char_repr });
        chars_to_variant.push(quote! { #char_repr => Ok(Self::#ident) });
    }

    // default width is the log2 of the max_variant
    let mut width = f32::ceil(f32::log2(max_variant as f32)) as u8;

    for attr in &enum_ast.attrs {
        if attr.path().is_ident("width") {
            width = match attr.parse_args::<syn::LitInt>() {
                Ok(w) => {
                    let chosen_width = w.base10_parse::<u8>().unwrap();
                    // test whether the specified width is too small
                    if chosen_width < width {
                        return syn::Error::new_spanned(
                            attr,
                            format!("Width is not large enough encode all variants (max: {width})"),
                        )
                        .to_compile_error()
                        .into();
                    }
                    chosen_width
                }
                Err(err) => return err.to_compile_error().into(),
            }
        };
    }

    let parse_error = quote! { crate::ParseBioError };

    // Generate the implementation
    let output = quote! {
        impl Codec for #enum_ident {
            type Error = #parse_error;
            const WIDTH: u8 = #width;
            fn unsafe_from_bits(b: u8) -> Self {
                //debug_assert!(false, "Invalid encoding: {b:?}");
                match b {
                    #(#unsafe_alts),*,
                    _ => unreachable!(),
                }
            }

            fn try_from_bits(b: u8) -> Result<Self, Self::Error> {
                match b {
                    #(#alt_discriminants),*,
                    _ => Err(#parse_error {}),
                }
            }

            fn from_char(c: char) -> Result<Self, Self::Error> {
                match c {
                    #(#chars_to_variant),*,
                    _ => Err(#parse_error {}),
                }
            }

            fn to_char(self) -> char {
                match self {
                    #(#variants_to_char),*,
                }
            }
        }

        impl #enum_ident {
            pub fn iter() -> impl Iterator<Item = Self> {
                vec![ #(Self::#variant_idents,)* ].into_iter()
            }
        }

        impl From<#enum_ident> for char {
            fn from(item: #enum_ident) -> char {
                item.to_char()
            }
        }

        impl core::convert::TryFrom<char> for #enum_ident {
            type Error = #parse_error;
            fn try_from(c: char) -> Result<Self, Self::Error> {
                Self::from_char(c)
            }
        }

        impl From<u8> for #enum_ident {
            fn from(b: u8) -> Self {
                Self::unsafe_from_bits(b)
            }
        }

        impl From<#enum_ident> for u8 {
            fn from(item: #enum_ident) -> Self {
                item as u8
            }
        }

        impl core::str::FromStr for #enum_ident {
            type Err = #parse_error;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                Self::try_from(s.chars().next().unwrap())
            }
        }

        impl core::fmt::Display for #enum_ident {
            fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
                write!(f, "{:?}", self)
            }
        }
    };
    output.into()
}
