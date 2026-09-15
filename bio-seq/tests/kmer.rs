use bio_seq::prelude::*;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Codec)]
#[repr(u8)]
enum Sparse {
    A = 1,
    C = 2,
}

macro_rules! storage_tests {
    ($module:ident, $storage:ty) => {
        mod $module {
            use super::*;

            #[cfg(feature = "extra_codecs")]
            #[test]
            fn reverse_ws_kmer() {
                use bio_seq::codec::degenerate::WS;
                let k: Kmer<WS, 2> = "WS".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "SW");
            }

            #[test]
            fn reverse_dna_kmer() {
                let k: Kmer<Dna, 2, $storage> = "AC".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "CA")
            }

            #[test]
            fn reverse_sparse_kmer() {
                let k: Kmer<Sparse, 2, $storage> = "AC".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "CA")
            }

            #[test]
            fn reverse_iupac_kmer() {
                let k: Kmer<Iupac, 2, $storage> = "AC".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "CA")
            }
        }
    };
}

storage_tests!(storage_usize, usize);
storage_tests!(storage_u64, u64);
storage_tests!(storage_u128, u128);
