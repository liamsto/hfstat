use core::fmt::{Display, Formatter};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Input(&'static str),
    Data(&'static str),
    Calculation(&'static str),
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        let (kind, message) = match self {
            Self::Input(message) => ("invalid input", message),
            Self::Data(message) => ("invalid model data", message),
            Self::Calculation(message) => ("calculation failed", message),
        };
        write!(f, "{kind}: {message}")
    }
}

impl std::error::Error for Error {}

pub type Result<T> = core::result::Result<T, Error>;
