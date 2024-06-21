// Copyright 2024 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use proc_macro2::TokenStream;

use quote::quote;

use syn::punctuated::Punctuated;
use syn::Token;

/// Allow the user to request more bits than used by their encodings
pub(crate) fn parse_width(attrs: &Vec<syn::Attribute>, max_variant: u8) -> Result<u8, syn::Error> {
    // minimum width is the log2 of the max_variant
    let min_width = f32::ceil(f32::log2(max_variant as f32)) as u8;

    for attr in attrs {
        if attr.path().is_ident("bits") {
            return match attr.parse_args::<syn::LitInt>() {
                Ok(w) => {
                    let chosen_width = w.base10_parse::<u8>().unwrap();
                    // test whether the specified width is too small
                    if chosen_width < min_width {
                        Err(syn::Error::new_spanned(
                            attr,
                            format!(
                                "Bit width is not large enough encode all variants (min: {min_width})"
                            ),
                        ))
                    } else {
                        Ok(chosen_width)
                    }
                }
                Err(err) => Err(err),
            };
        };
    }
    Ok(min_width)
}

pub(crate) fn test_repr(enum_ast: &syn::ItemEnum) -> Result<(), syn::Error> {
    // Test that enum is repr(u8)
    for attr in enum_ast.attrs.iter() {
        if attr.path().is_ident("repr") {
            match attr.parse_args::<syn::Ident>() {
                Ok(ident) => {
                    if ident == "u8" {
                        return Ok(());
                    } else {
                        return Err(syn::Error::new_spanned(
                            &enum_ast.ident,
                            "Enums deriving Codec must be annotated with #[repr(u8)]",
                        ));
                    }
                }
                Err(err) => return Err(err),
            }
        }
    }

    Err(syn::Error::new_spanned(
        &enum_ast.ident,
        "Enums deriving Codec must be annotated with #[repr(u8)]",
    ))
}

pub(crate) struct CodecVariants {
    pub(crate) idents: Vec<syn::Ident>,
    pub(crate) to_chars: Vec<TokenStream>,
    pub(crate) from_chars: Vec<TokenStream>,
    pub(crate) alts: Vec<TokenStream>,
    pub(crate) unsafe_alts: Vec<TokenStream>,
    pub(crate) max_discriminant: u8,
}

pub(crate) fn parse_variants(
    variants: &syn::punctuated::Punctuated<syn::Variant, syn::token::Comma>,
) -> Result<CodecVariants, syn::Error> {
    let mut max_discriminant = 0u8;
    let mut idents = Vec::new();

    let mut to_chars = Vec::new();
    let mut from_chars = Vec::new();
    let mut alts = Vec::new();
    let mut unsafe_alts = Vec::new();

    for variant in variants.iter() {
        let ident = &variant.ident;
        idents.push(ident.clone());
        let discriminant = &variant.discriminant;

        if let Some((_, syn::Expr::Lit(expr_lit))) = discriminant {
            let value = match &expr_lit.lit {
                // discriminants must be either integers or byte literals
                syn::Lit::Byte(lit_byte) => lit_byte.value(),
                syn::Lit::Int(lit_int) => lit_int.base10_parse::<u8>().unwrap(),
                _ => {
                    return Err(syn::Error::new_spanned(
                        ident,
                        "Codec derivations require byte or integer discriminants",
                    ))
                }
            };

            alts.push(quote! { #value => Some(Self::#ident) });
            unsafe_alts.push(quote! { #value => Self::#ident });

            max_discriminant = max_discriminant.max(value);
        } else {
            return Err(syn::Error::new_spanned(
                ident,
                "Codec derivations require discriminants",
            ));
        }

        //let mut char_repr = ident.to_string().chars().next().unwrap();

        let mut char_repr = ident.to_string().bytes().next().unwrap();

        for attr in &variant.attrs {
            if attr.path().is_ident("display") {
                let alt_attr: syn::LitChar = match attr.parse_args() {
                    Ok(attr) => attr,
                    Err(err) => return Err(err),
                };
                char_repr = alt_attr.value() as u8;
            } else if attr.path().is_ident("alt") {
                let discs: Punctuated<syn::ExprLit, Token![,]> =
                    match attr.parse_args_with(Punctuated::parse_terminated) {
                        Ok(discs) => discs,
                        Err(err) => return Err(err),
                    };

                for d in discs.into_iter() {
                    alts.push(quote! { #d => Some(Self::#ident) });
                    unsafe_alts.push(quote! { #d => Self::#ident });
                }
            };
        }

        to_chars.push(quote! { Self::#ident => #char_repr });
        from_chars.push(quote! { #char_repr => Some(Self::#ident) });
    }

    Ok(CodecVariants {
        idents,
        to_chars,
        from_chars,
        alts,
        unsafe_alts,
        max_discriminant,
    })
}
