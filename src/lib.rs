/*! # `bioseq`

Bit packed biological sequences

!*/

mod alphabet;

#[macro_use]
mod seq;

mod kmer;

#[cfg(test)]
mod tests {
    use crate::alphabet::dna::Dna;
    use crate::alphabet::dna::Dna::{A, C, G, T};
    use crate::seq::Seq;

    #[test]
    fn make_from_vector() {
        assert_eq!(Seq::from_vec(vec![A, C, G, T]).to_usize(), 0b11011000);
        assert_eq!(Seq::from_vec(vec![C, G, C, G]).to_usize(), 0b01100110);
        assert_eq!(Seq::from_vec(vec![T, T]).to_usize(), 0b1111);
        assert_eq!(
            Seq::from_vec(vec![C, G, T, A, C, G, A, T]).to_usize(),
            0b1100011000110110
        );
        assert_eq!(Seq::from_vec(vec![A,]).to_usize(), 0b00);
    }

    #[test]
    fn test_display_dna() {
        let seq = Seq::from_vec(vec![A, C, G, T, T, A, T, C]);
        assert_eq!(format!("{}", seq), "[A, C, G, T, T, A, T, C]");
        assert_eq!(format!("{}", dna!("ACGT")), "[A, C, G, T]");
    }

    //    #[test]
    //    fn complement_iupac() {
    //        assert_eq!((iupac!([])).complement(), iupac!("TAGC"));
    //    }
}
