use core::fmt;
use core::result;

#[derive(Debug, PartialEq)]
pub enum ParseBioError {
    UnrecognisedBase(u8),
    MismatchedLength(usize, usize),
    SequenceTooLong(usize, usize),
}

pub type Result<T> = result::Result<T, ParseBioError>;

impl fmt::Display for ParseBioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseBioError::UnrecognisedBase(byte) => {
                if byte.is_ascii_alphanumeric() {
                    write!(
                        f,
                        "Unrecognised character: '{}' ({byte:#04X?})",
                        *byte as char,
                    )
                } else {
                    write!(f, "Unrecognised character: {byte:#04X?}")
                }
            }

            ParseBioError::MismatchedLength(got, expected) => {
                write!(f, "Expected length {expected}, got {got}")
            }

            ParseBioError::SequenceTooLong(got, expected) => {
                write!(f, "Expected length <= {expected}, got {got}")
            }
        }
    }
}

impl From<ParseBioError> for std::io::Error {
    fn from(err: ParseBioError) -> Self {
        std::io::Error::other(err)
    }
}

// #![feature(error_in_core)
impl std::error::Error for ParseBioError {}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn test_parse_error() {
        let seq = Seq::<Dna>::from_str("ACGTx").unwrap_err();
        assert_eq!(
            format!("{seq}"),
            "Unrecognised character: 'x' (0x78)".to_string()
        );

        let seq = Seq::<Dna>::from_str("ACGT\n").unwrap_err();
        assert_eq!(format!("{seq}"), "Unrecognised character: 0x0A".to_string());

        let err: std::io::Error = seq.into();
        assert_eq!(err.kind(), std::io::ErrorKind::Other);

        let seq = Seq::<Dna>::from_str("ACGTT");
        assert!(seq.is_ok());
    }

    #[test]
    fn test_mismatched_length() {
        let seq: &SeqSlice<Dna> = dna!("ACGTGAT");
        let kmer = Kmer::<Dna, 14>::try_from(seq).unwrap_err();
        assert_eq!(format!("{kmer}"), "Expected length 7, got 14".to_string());

        let kmer = Kmer::<Dna, 6>::try_from(seq).unwrap_err();
        assert_eq!(format!("{kmer}"), "Expected length 7, got 6".to_string());

        let kmer = Kmer::<Dna, 7>::try_from(seq);
        assert!(kmer.is_ok());
    }

    #[test]
    fn test_size_error() {
        let seq = dna!("ACGTACGTACGTACGTACGTACGTACGTACGTA");
        let int = TryInto::<usize>::try_into(seq).unwrap_err();

        let expected_size = usize::BITS / Dna::BITS as u32;

        assert_eq!(
            format!("{int}"),
            format!("Expected length <= {expected_size}, got 33")
        );

        let good = TryInto::<usize>::try_into(&seq[..expected_size as usize]);
        assert!(good.is_ok());
    }
}
