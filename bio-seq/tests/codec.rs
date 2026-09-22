use bio_seq::prelude::*;

macro_rules! codec_tests {
    ($module:ident, $codec:ty) => {
        mod $module {
            use super::*;

            #[test]
            #[should_panic]
            fn from_invalid_bits() {
                let _: $codec = <$codec>::unsafe_from_bits(255u8);
            }

            #[test]
            fn from_invalid_bits_checked() {
                assert_eq!(None, <$codec>::try_from_bits(255u8));
            }

            #[test]
            fn from_invalid_ascii_checked() {
                assert_eq!(None, <$codec>::try_from_ascii(255u8));
            }
        }
    };
}

codec_tests!(dna_from, Dna);
codec_tests!(iupac_from, Iupac);
codec_tests!(amino_from, Amino);

#[cfg(feature = "extra_codecs")]
mod extras {
    use bio_seq::codec::Codec;
    use bio_seq::codec::degenerate::{MK, RY, WS};

    codec_tests!(ws_from, WS);
    codec_tests!(mk_from, MK);
    codec_tests!(ry_from, RY);
}
