use core::fmt;

#[derive(Debug)]
pub struct ParseBioError {}

impl fmt::Display for ParseBioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Error parsing sequence")
    }
}

// #![feature(error_in_core)
impl std::error::Error for ParseBioError {}
