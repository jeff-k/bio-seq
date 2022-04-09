extern crate proc_macro;

use crate::proc_macro::TokenStream;

use quote::quote;

use syn::{parse_macro_input, DeriveInput};

#[proc_macro_derive(Codec)]
pub fn codec_derive(input: TokenStream) -> TokenStream {
    let DeriveInput { ident, data, .. } = parse_macro_input!(input);

    /*
    for variant in &data_enum.variants {
        let ref variant_name = variant.ident;
    }
    */

    let output = quote! {
        impl Codec for #ident {
            const WIDTH: u8 = 2;
            fn unsafe_from_bits(b: u8) -> Self {
                match b {
                    _ => Self::A,
                    //#(#repr => #ident::#name)*,
                }
            }

            fn try_from_bits(b: u8) -> Result<Self, ParseBioErr> {
                match b {
                    //#repr => Ok(#ident::#variant_name),
                    _ => Err(ParseBioErr),
                }
            }

            fn from_char(c: char) -> Result<Self, ParseBioErr> {
                Ok(Self::A)
                //Ok(Self::#variant_name)
            }

            fn to_char(a: Self) -> char {
                match a {
                    _ => 'A',
                    //#variant_name => #repr,
                }
            }
        }
    };
    //    panic!("{}", output);
    output.into()
}
