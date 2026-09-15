use bio_seq::prelude::*;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Codec)]
#[repr(u8)]
#[bits(3)]
enum Sparse {
    A = 1,
    C = 2,
    X = 5,
}

macro_rules! storage_tests {
    ($module:ident, $storage:ty) => {
        mod $module {
            use super::*;

            #[cfg(feature = "extra_codecs")]
            #[test]
            fn reverse_ws_kmer() {
                use bio_seq::codec::degenerate::WS;
                let k: Kmer<WS, 2, $storage> = "WS".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "SW");
            }

            #[test]
            fn reverse_dna_kmer() {
                let k: Kmer<Dna, 3, $storage> = "ACG".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "GCA");
                assert_eq!(k.to_rev().to_comp().to_string(), "CGT");

                const FULL_K: usize = <$storage>::BITS as usize / 2;
                let fk: Kmer<Dna, FULL_K, $storage> = "ACGT".repeat(FULL_K / 4).parse().unwrap();
                assert_eq!(fk.to_comp().to_string(), "TGCA".repeat(FULL_K / 4));
            }

            #[test]
            #[should_panic]
            fn reverse_sparse_kmer() {
                let k: Kmer<Sparse, 3, $storage> = "ACX".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "XCA");
            }

            #[test]
            fn reverse_iupac_kmer() {
                let k: Kmer<Iupac, 2, $storage> = "AC".parse().unwrap();
                assert_eq!(k.to_rev().to_string(), "CA");
            }
        }
    };
}

storage_tests!(storage_usize, usize);
storage_tests!(storage_u64, u64);
storage_tests!(storage_u128, u128);
