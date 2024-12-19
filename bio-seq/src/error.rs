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
