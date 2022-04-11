extern crate proc_macro;

use crate::proc_macro::TokenStream;

use quote::quote;

use syn::{parse_macro_input, DeriveInput, Item};

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

fn codec_width(attrs: &Vec<syn::Attribute>) -> u8 {
    let mut width: u8 = 4;
    for attr in attrs {
        if attr.path.is_ident("width") {
            width = match syn::parse2::<WidthAttr>(attr.tokens.clone()) {
                Ok(w) => w.width.base10_parse::<u8>().unwrap(),
                _ => 4,
            };

            //panic!("width {}", width);
        };
    }
    width
}

#[proc_macro_derive(Codec, attributes(width))]
pub fn codec_derive(input: TokenStream) -> TokenStream {
    //    let ast = parse_macro_input!(input as DeriveInput);

    let ast = parse_macro_input!(input as syn::Item);

    //let width = codec_width(&ast.into());

    //let ident = ast.ident;

    match ast {
        //let variants = match ast.data {
        syn::Item::Enum(e) => {
            let variants = e.variants;
            let ident = e.ident;
            let width: u8 = codec_width(&e.attrs);
            // TODO
            // We can check whether the max value of the variants fits into the
            // specified bit width
            // We want to run this check on the alternative variants as well
            let max_variant = 255;
            //let width: u8 = 4;

            for v in variants.iter() {
                let ident = &v.ident;
                let discriminant = &v.discriminant;

                match discriminant {
                    // get max value out of the discriminant Expr (should be LitByte() with value member)
                    Some((_, d)) => (),
                    _ => panic!(),
                    //None => return syn::Error::new_spanned(ident, "Codec derivations require discriminants").to_compile_error().into(),
                }
            }

            if max_variant < u64::pow(2, width as u32) {
                return syn::Error::new_spanned(
                    ident,
                    "Codec specifier does not fit inside bit width",
                )
                .to_compile_error()
                .into();
            }

            let variants_to_char = variants.iter().map(|v| {
                let ident = &v.ident;
                let discriminant = &v.discriminant;

                match discriminant {
                    Some(d) => {
                        let char_repr: char = ident.to_string().chars().next().unwrap();
                        quote! { Self::#ident => #char_repr }
                    }
                    None => panic!(),
                }
            });
            let chars_to_variant = variants.iter().map(|v| {
                let ident = &v.ident;
                let discriminant = &v.discriminant;

                match discriminant {
                    Some(d) => {
                        let char_repr: char = ident.to_string().chars().next().unwrap();
                        quote! { #char_repr => Ok(Self::#ident) }
                    }
                    None => panic!(),
                }
            });
            let bits_to_variant = variants.iter().map(|v| {
                let ident = &v.ident;
                let discriminant = &v.discriminant;

                match discriminant {
                    Some((_, d)) => {
                        quote! { #d => Ok(Self::#ident) }
                    }
                    None => panic!(),
                }
            });
            let unsafe_bit_variants = variants.iter().map(|v| {
                let ident = &v.ident;
                let discriminant = &v.discriminant;

                match discriminant {
                    Some((_, d)) => {
                        quote! { #d => Self::#ident }
                    }
                    None => panic!(),
                }
            });

            let output = quote! {
                impl Codec for #ident {
                    const WIDTH: u8 = #width;
                    fn unsafe_from_bits(b: u8) -> Self {
                        match b {
                            #(#unsafe_bit_variants),*,
                            _ => panic!(),
                        }
                    }

                    fn try_from_bits(b: u8) -> Result<Self, ParseBioErr> {
                        match b {
                            #(#bits_to_variant),*,
                            _ => Err(ParseBioErr),
                        }
                    }

                    fn from_char(c: char) -> Result<Self, ParseBioErr> {
                        match c {
                            #(#chars_to_variant),*,
                            _ => Err(ParseBioErr),
                        }
                    }

                    fn to_char(a: Self) -> char {
                        match a {
                            #(#variants_to_char),*,
                        }
                    }
                }
            };
            output.into()
        }
        _ => panic!("not enum"),
    }

    //    panic!("{}", output);
}
