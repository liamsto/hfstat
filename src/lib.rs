#![doc = include_str!("../README.md")]

pub mod antenna;
mod data;
pub mod error;
pub mod input;
pub mod output;

mod p533;

pub use antenna::{Antenna, Orientation, Pattern};
pub use error::{Error, Result};
pub use input::{DigitalConfig, Location, ManMadeNoise, Modulation, Month, Path, Request};
pub use output::Prediction;
pub use p533::predict;
