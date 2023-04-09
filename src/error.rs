use std::fmt;

#[derive(Debug)]
pub struct ParseBioError {}

impl fmt::Display for ParseBioError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Error parsing sequence")
    }
}

impl std::error::Error for ParseBioError {}
