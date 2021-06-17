mod alphabet;

#[macro_use]
mod seq;

#[cfg(test)]
mod tests {
    use crate::alphabet::dna::Dna::*;
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
    fn test_display() {
        let seq = Seq::from_vec(vec![T, A, C, G, T, T]);
        assert_eq!(format!("{}", seq), "[A, C, G, T, T, T, T]");
    }

    //    #[test]
    //    fn make_dnaseq_from_macro() {
    //        assert_eq!(dna!([A, T, C, G]), Seq::<Dna>::from_vec([Dna::A, Dna::T, Dna::C, Dna::G]));
    //    }

    //    #[test]
    //    fn complement_iupac() {
    //        assert_eq!((iupac!([])).complement(), iupac!("TAGC"));
    //    }
}
