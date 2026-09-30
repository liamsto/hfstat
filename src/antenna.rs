use core::fmt::{Debug, Formatter};

use crate::error::{Error, Result};

pub const AZIMUTHS: usize = 360;
pub const ELEVATIONS: usize = 91;
pub const GAINS_PER_FREQUENCY: usize = AZIMUTHS * ELEVATIONS;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Orientation {
    #[default]
    Path,
    BearingDeg(f64),
}

#[derive(Clone, Copy, PartialEq)]
pub struct Pattern<'a> {
    frequencies_mhz: &'a [f64],
    gains_dbi: &'a [f64],
}

impl Debug for Pattern<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Pattern")
            .field("frequency_count", &self.frequencies_mhz.len())
            .field("gain_count", &self.gains_dbi.len())
            .finish()
    }
}

impl<'a> Pattern<'a> {
    pub fn new(frequencies_mhz: &'a [f64], gains_dbi: &'a [f64]) -> Result<Self> {
        if frequencies_mhz.is_empty() {
            return Err(Error::Data("antenna pattern has no frequencies"));
        }
        let Some(expected) = frequencies_mhz.len().checked_mul(GAINS_PER_FREQUENCY) else {
            return Err(Error::Data("antenna pattern is too large"));
        };
        if gains_dbi.len() != expected {
            return Err(Error::Data(
                "antenna gains must contain 360 by 91 values per frequency",
            ));
        }
        if frequencies_mhz
            .iter()
            .any(|frequency| !frequency.is_finite() || *frequency <= 0.0)
        {
            return Err(Error::Data(
                "antenna frequencies must be finite and positive",
            ));
        }
        if frequencies_mhz.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(Error::Data(
                "antenna frequencies must be strictly increasing",
            ));
        }
        if gains_dbi.iter().any(|gain| !gain.is_finite()) {
            return Err(Error::Data("antenna gains must be finite"));
        }
        Ok(Self {
            frequencies_mhz,
            gains_dbi,
        })
    }

    #[must_use]
    pub const fn frequencies_mhz(self) -> &'a [f64] {
        self.frequencies_mhz
    }

    #[must_use]
    pub const fn gains_dbi(self) -> &'a [f64] {
        self.gains_dbi
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Model<'a> {
    Isotropic(f64),
    Pattern(Pattern<'a>),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Antenna<'a> {
    model: Model<'a>,
    orientation: Orientation,
}

impl<'a> Antenna<'a> {
    #[must_use]
    pub const fn isotropic(gain_dbi: f64) -> Self {
        Self {
            model: Model::Isotropic(gain_dbi),
            orientation: Orientation::Path,
        }
    }

    #[must_use]
    pub const fn pattern(pattern: Pattern<'a>, orientation: Orientation) -> Self {
        Self {
            model: Model::Pattern(pattern),
            orientation,
        }
    }

    pub(crate) const fn model(self) -> Model<'a> {
        self.model
    }

    pub(crate) const fn orientation(self) -> Orientation {
        self.orientation
    }

    pub(crate) fn validate(self, frequency_mhz: f64) -> Result<()> {
        if let Orientation::BearingDeg(bearing) = self.orientation
            && !bearing.is_finite()
        {
            return Err(Error::Input("antenna bearing must be finite"));
        }

        match self.model {
            Model::Isotropic(gain) if gain.is_finite() => Ok(()),
            Model::Isotropic(_) => Err(Error::Input("isotropic antenna gain must be finite")),
            Model::Pattern(pattern) => {
                let frequencies = pattern.frequencies_mhz;
                if frequency_mhz < frequencies[0]
                    || frequency_mhz > frequencies[frequencies.len() - 1]
                {
                    return Err(Error::Input(
                        "operating frequency is outside the antenna pattern range",
                    ));
                }
                Ok(())
            }
        }
    }
}

impl Default for Antenna<'_> {
    fn default() -> Self {
        Self::isotropic(0.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::antenna::{GAINS_PER_FREQUENCY, Pattern};

    #[test]
    fn pattern_validates_shape_and_frequency_order() {
        let gains = vec![0.0; GAINS_PER_FREQUENCY * 2];
        assert!(Pattern::new(&[10.0, 20.0], &gains).is_ok());
        assert!(Pattern::new(&[20.0, 10.0], &gains).is_err());
        assert!(Pattern::new(&[10.0], &gains).is_err());
    }

    #[test]
    fn pattern_rejects_nonfinite_values() {
        let mut gains = vec![0.0; GAINS_PER_FREQUENCY];
        gains[0] = f64::NAN;
        assert!(Pattern::new(&[10.0], &gains).is_err());
        assert!(Pattern::new(&[f64::INFINITY], &gains).is_err());
    }
}
