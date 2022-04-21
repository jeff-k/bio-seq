// Copyright 2021 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

/*! # bio-seq

Bit packed types for biological [sequences](seq) and [k-mers](kmer).

Add `bio-seq` to your project's `Cargo.toml`:

```toml
[dependencies]
bio-seq = "*"
```

## Example: Kmers

```rust
use bio_seq::codec::Dna;
use bio_seq::*;

let seq = dna!("ACTGCTAGCA");

for kmer in seq.kmers::<8>() {
        println!("{}", kmer);
}
```

!*/

pub mod bigk;
pub mod codec;

#[macro_use]
pub mod seq;

pub mod kmer;

//pub use codec::iupac::Iupac;
pub use seq::Seq;
pub use std::str::FromStr;

pub trait ReverseComplement {
    fn rc(s: Self) -> Self;
}

#[cfg(test)]
mod tests {
    use crate::codec::amino::Amino;
    use crate::codec::dna::Dna;
    use crate::codec::dna::Dna::{A, C, G, T};
    use crate::codec::iupac::Iupac;
    use crate::seq::Seq;
    use std::str::FromStr;

    #[test]
    fn make_from_vector() {
        assert_eq!(Seq::from_vec(vec![A, C, G, T]).raw(), &[0b11_10_01_00]);
        assert_eq!(Seq::from_vec(vec![C, G, C, G]).raw(), &[0b10_01_10_01]);
        assert_eq!(Seq::from_vec(vec![T, T]).raw(), &[0b11_11]);
        assert_eq!(Seq::from_vec(vec![T, C, A]).raw(), &[0b00_01_11]);
        assert_eq!(Seq::from_vec(vec![T, G, A]).raw(), &[0b00_10_11]);

        assert_eq!(
            Seq::from_vec(vec![C, G, T, A, C, G, A, T]).raw(),
            &[0b11_00_10_01_00_11_10_01]
        );
        assert_eq!(Seq::from_vec(vec![A,]).raw(), &[0b00]);
    }

    #[test]
    fn test_display_dna() {
        let seq = Seq::from_vec(vec![A, C, G, T, T, A, T, C]);
        assert_eq!(format!("{}", seq), "ACGTTATC");
        assert_eq!(format!("{}", dna!("ACGT")), "ACGT");
    }

    #[test]
    fn iterate_bases() {
        let seq = dna!("ACGTACGT");
        assert_eq!(
            seq.into_iter().collect::<Vec<Dna>>(),
            vec![A, C, G, T, A, C, G, T]
        );
    }

    #[test]
    fn rev_seq() {
        let seq = dna!("ACGTACGT");
        assert_eq!(
            seq.rev().collect::<Vec<Dna>>(),
            vec![T, G, C, A, T, G, C, A]
        );

        assert_eq!(
            iupac!("GNX").rev().collect::<Vec<Iupac>>(),
            vec![Iupac::X, Iupac::N, Iupac::G]
        );

        assert_eq!(
            amino!("DCMNLKGHI").rev().collect::<Vec<Amino>>(),
            vec![
                Amino::I,
                Amino::H,
                Amino::G,
                Amino::K,
                Amino::L,
                Amino::N,
                Amino::M,
                Amino::C,
                Amino::D
            ]
        );
    }

    #[test]
    fn iterate_kmers() {
        let seq = dna!("ACGTAAGGGG");
        for (kmer, answer) in seq
            .kmers::<4>()
            .zip(["ACGT", "CGTA", "GTAA", "TAAG", "AAGG", "AGGG", "GGGG"])
        {
            assert_eq!(format!("{}", kmer), answer);
        }
    }

    #[test]
    fn iterate_kmer8() {
        let seq = dna!("AAAACCCCGGGG");
        for (kmer, answer) in seq
            .kmers::<8>()
            .zip(["AAAACCCC", "AAACCCCG", "AACCCCGG", "ACCCCGGG", "CCCCGGGG"])
        {
            assert_eq!(format!("{}", kmer), answer);
        }
    }

    #[test]
    fn iterate_kmer4() {
        let seq = dna!("AAAACCCCGGGGTTTT");
        for (kmer, answer) in seq.kmers::<4>().zip([
            "AAAA", "AAAC", "AACC", "ACCC", "CCCC", "CCCG", "CCGG", "CGGG", "GGGG", "GGGT", "GGTT",
            "GTTT", "TTTT",
        ]) {
            assert_eq!(format!("{}", kmer), answer);
        }
    }

    #[test]
    fn iupac_bitwise_ops() {
        assert_eq!(
            format!("{}", iupac!("ASXGYTNA") | iupac!("ANTGCATX")),
            "ANTGYWNA"
        );
        assert_eq!(
            format!("{}", iupac!("ACGTSWKM") & iupac!("WKMSTNNA")),
            "AXXXXWKA"
        );
    }

    #[test]
    fn nth_chars() {
        assert_eq!(iupac!("ACGTRYSWKMBDHVNX").nth(0), Iupac::A);
        assert_ne!(iupac!("ACGTRYSWKMBDHVNX").nth(0), Iupac::C);
        assert_eq!(iupac!("ACGTRYSWKMBDHVNX").nth(15), Iupac::X);
        assert_eq!(iupac!("ACGTRYSWKMBDHVNX").nth(3), Iupac::from(Dna::T));
        assert_ne!(iupac!("ACGTRYSWKMBDHVNX").nth(3), Iupac::from(Dna::G));

        assert_eq!(amino!("DCMNLKGHI").nth(1), Amino::C);
        assert_ne!(amino!("DCMNLKGHI").nth(7), Amino::I);
    }
}
