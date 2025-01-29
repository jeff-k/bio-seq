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
        std::io::Error::new(std::io::ErrorKind::Other, err)
    }
}

// #![feature(error_in_core)
impl std::error::Error for ParseBioError {}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn test_parse_error() {
        let seq = Seq::<Dna>::from_str("ACGTx");
        assert_eq!(seq, Err(ParseBioError::UnrecognisedBase(b'x')));

        let seq = Seq::<Dna>::from_str("ACGT\n");
        assert_eq!(seq, Err(ParseBioError::UnrecognisedBase(10)));

        let err: std::io::Error = seq.unwrap_err().into();
        assert_eq!(err.kind(), std::io::ErrorKind::Other);

        let seq = Seq::<Dna>::from_str("ACGTT");
        assert!(seq.is_ok());
    }

    #[test]
    fn test_mismatched_length() {
        let seq: &SeqSlice<Dna> = dna!("ACGTGAT");
        let kmer = Kmer::<Dna, 14>::try_from(seq);
        assert_eq!(kmer, Err(ParseBioError::MismatchedLength(14, 7)));

        let kmer = Kmer::<Dna, 6>::try_from(seq);
        assert_eq!(kmer, Err(ParseBioError::MismatchedLength(6, 7)));

        let kmer = Kmer::<Dna, 7>::try_from(seq);
        assert!(kmer.is_ok());
    }

    #[test]
    fn test_size_error() {
        let seq = dna!("ACGTACGTACGTACGTACGTACGTACGTACGTA");
        let int = TryInto::<usize>::try_into(seq);

        assert_eq!(int, Err(ParseBioError::SequenceTooLong(33, 32)));

        let good = TryInto::<usize>::try_into(&seq[..32]);
        assert!(good.is_ok());
    }
}
