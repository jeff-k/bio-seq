extern crate proc_macro;

use crate::proc_macro::TokenStream;

use quote::quote;

use syn::{parse_macro_input, DeriveInput};

#[proc_macro_derive(Codec)]
pub fn codec_derive(input: TokenStream) -> TokenStream {
    let ast = parse_macro_input!(input as DeriveInput);
    let ident = ast.ident;

    /*
             let ident = &s.ident;
             let fields_ty = s.fields.iter().map(|field| &field.ty);

             let attr_check = s.fields.iter().filter_map(|field| {
                 let ty = &field.ty;
                 let attrs = &field.attrs;
                 for attr in attrs {
                     if attr.path.is_ident("width") {
                         let width = syn::parse2::<WidthAttribute>(attr.tokens.clone()).ok()?.width;
                         quote_spanned!(width.span() => const _: [(); #width] = [(); #ty::WIDTH];),
                         );
                     }
                 }
                 None
             });
    */

    let variants = match ast.data {
        syn::Data::Enum(e) => e.variants,
        _ => {
            return syn::Error::new_spanned(ident, "Can only derive bio-seq codecs on enum")
                .to_compile_error()
                .into()
        }
    };

    // TODO
    // We can check whether the max value of the variants fits into the
    // specified bit width
    // We want to run this check on the alternative variants as well
    let max_variant = 15;
    let width: u8 = 4;
    if max_variant >= u64::pow(2, width as u32) {
        return syn::Error::new_spanned(ident, "Codec specifier does not fit inside bit width")
            .to_compile_error()
            .into();
    }

    let output = quote! {
        impl Codec for #ident {
            const WIDTH: u8 = #width;
            fn unsafe_from_bits(b: u8) -> Self {
                match b {
                    //#(#bit_variants),*,
                    _ => Self::A,
                    //#(#repr => #ident::#name)*,
                }
            }

            fn try_from_bits(b: u8) -> Result<Self, ParseBioErr> {
                match b {
                    //#(#bit_variants),*,
                    //#repr => Ok(#ident::#variant_name),
                    _ => Err(ParseBioErr),
                }
            }

            fn from_char(c: char) -> Result<Self, ParseBioErr> {
                match c {
                    //#(#match_variants),*,
                    _ => Err(ParseBioErr),
                }
                //Ok(Self::#variant_name)
            }

            fn to_char(a: Self) -> char {
                match a {
                    _ => 'B',
                    //#(#match_variants),*,
                    //#variant_name => #repr,
                }
            }
        }
    };
    //    panic!("{}", output);
    output.into()
}
