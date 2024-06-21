// Copyright 2023 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! `bio-seq-derive` is a procedural macro crate that provides the `Codec` derive macro for the `bio-seq` library.
//! It allows users to define custom bit-packed encodings from an enum. The representation of the enum is derived from the discriminants.
//! Please refer to the `bio-seq` [documentation](https://github.com/jeff-k/bio-seq) for a complete guide on defining custom alphabets.

use proc_macro2::TokenStream;

use quote::quote;

use syn::punctuated::Punctuated;
use syn::{parse_macro_input, Token};

#[proc_macro_derive(Codec, attributes(bits, display, alt))]
pub fn codec_derive(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
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

            fn items() -> impl Iterator<Item = Self> {
                vec![ #(Self::#idents,)* ].into_iter()
            }
        }

    };
    output.into()
}

/// Allow the user to request more bits than used by their encodings
fn parse_width(attrs: &Vec<syn::Attribute>, max_variant: u8) -> Result<u8, syn::Error> {
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

fn test_repr(enum_ast: &syn::ItemEnum) -> Result<(), syn::Error> {
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

struct CodecVariants {
    idents: Vec<syn::Ident>,
    to_chars: Vec<TokenStream>,
    from_chars: Vec<TokenStream>,
    alts: Vec<TokenStream>,
    unsafe_alts: Vec<TokenStream>,
    max_discriminant: u8,
}

fn parse_variants(
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
