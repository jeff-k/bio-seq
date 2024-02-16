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

/// The width attribute should take the form #[width = 2]
struct WidthAttr {
    width: syn::LitInt,
}

impl syn::parse::Parse for WidthAttr {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        let _: syn::Token![=] = input.parse()?;
        let width: syn::LitInt = input.parse()?;
        Ok(Self { width })
    }
}

/// Alternate representation attributes look like `#[alt(0b00, 0b11)]`
struct AltAttr {
    chr: syn::LitChar,
}

impl syn::parse::Parse for AltAttr {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        let _: syn::Token![=] = input.parse()?;
        let chr: syn::LitChar = input.parse()?;
        Ok(Self { chr })
    }
}

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
            if attr.path.is_ident("display") {
                let alt_attr: AltAttr = match syn::parse2(attr.tokens.clone()) {
                    Ok(attr) => attr,
                    Err(err) => return err.to_compile_error().into(),
                };
                char_repr = alt_attr.chr.value();
            } else if attr.path.is_ident("alt") {
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
        if attr.path.is_ident("width") {
            width = match syn::parse2::<WidthAttr>(attr.tokens.clone()) {
                Ok(w) => {
                    let chosen_width = w.width.base10_parse::<u8>().unwrap();
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
                match b {
                    #(#unsafe_alts),*,
                    _ => panic!(),
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
    };
    output.into()
}
