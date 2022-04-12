//! Translation tables for coding sequences

//pub mod default;

pub trait Translate {
    fn kmer_to_amino();
    fn amino_to_iupac();
}

#[cfg(test)]
mod tests {
    #[test]
    fn kmer_to_amino() {
    }
}
