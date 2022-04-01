extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn;

#[proc_macro_derive(Codec)]
pub fn codec_derive(input: TokenStream) -> TokenStream {
    let ast = syn::parse(input).unwrap();

    impl_codec(&ast)
}

fn impl_codec(ast: &syn::DeriveInput) -> TokenStream {
/*
    for variant in &ast.data.variants {

        quote! {
                #name::#key => BitArray::new(#value),
        }
    }
*/

    let name = &ast.ident;
    let gen = quote! {
        impl Codec for #name {
            const WIDTH: usize = 2;

            fn to_bits(&self) -> BitArray<Msb0, u8> {
                match &self { 
                    _ => unimplemented!()
                }
            }
            
            fn from_bits(b: &BitSlice<Msb0, u8>) -> Self {
                if b == BitArray::<Msb0, u8>::new(0b00)[8 - Self::WIDTH..] {
                    #name::A
                } else if b == BitArray::<Msb0, u8>::new(0b00)[8 - Self::WIDTH..] {
                    #name::C
                } else {
                    panic!()
                }
            }

            fn from_char(c: &u8) -> Result<Self, ParseBioErr> {
                match c {
                    b'A' => Ok(#name::A),
                    _ => Err(ParseBioErr),
                }
            }

            fn to_char(c: Self) -> u8 {
                match c {
                    #name::A => b'A',
                }
            }
        }
    };
    gen.into()
}
