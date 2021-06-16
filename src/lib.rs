mod alphabets;

use alphabets::{Amino, Dna, Iupac};
use bitvec::prelude::*;

#[cfg(test)]
mod tests {
    #[test]
    fn make_dnaseq_from_macro() {
        assert_eq!(dna!"ATCG", [Dna::A, Dna::T, Dna::C, Dna::G]);
    }
}
