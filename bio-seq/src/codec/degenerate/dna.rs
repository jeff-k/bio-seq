use crate::codec::Codec;
use crate::codec::Complement;
use crate::seq::{ReverseComplement, Seq};

/// 1-bit encoding for `S`trong (`G`/`C`) and `W`eak (`A`/`T`) nucleotides
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Codec)]
#[bits(1)]
#[repr(u8)]
pub enum Dna {
    W = 0b0,
    S = 0b1,
}

impl Complement for Dna {
    /// This representation collapses complements; `comp(x) == x`
    fn comp(&self) -> Self {
        *self
    }
}

impl ReverseComplement for Seq<Dna> {
    type Output = Self;

    fn revcomp(&self) -> Self {
        let mut bv = self.bv.clone();
        bv.reverse();
        Self::from(bv)
    }
}

#[cfg(test)]
mod tests {
    use crate::codec::masked;
    use crate::prelude::*;
}
