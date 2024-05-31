use core::fmt;

#[derive(Debug, PartialEq)]
pub enum ParseBioError {
    UnrecognisedBase(u8),
    MismatchedLength(usize, usize),
}

impl fmt::Display for ParseBioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseBioError::UnrecognisedBase(byte) => {
                write!(f, "Unrecognised character: {byte}")
            }
            ParseBioError::MismatchedLength(got, expected) => {
                write!(f, "Expected length {expected}, got {got}")
            }
        }
    }
}

// #![feature(error_in_core)
impl std::error::Error for ParseBioError {}

/// Error conditions for codon/amino acid translation
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum TranslationError {
    AmbiguousCodon,
    InvalidCodon,
}

impl fmt::Display for TranslationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            TranslationError::AmbiguousCodon => write!(f, ""),
            TranslationError::InvalidCodon => write!(f, ""),
        }
    }
}

// #![feature(error_in_core)
impl std::error::Error for TranslationError {}
