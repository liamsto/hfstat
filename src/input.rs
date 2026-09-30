use crate::{
    antenna::Antenna,
    error::{Error, Result},
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum Month {
    #[default]
    Jan,
    Feb,
    Mar,
    Apr,
    May,
    Jun,
    Jul,
    Aug,
    Sep,
    Oct,
    Nov,
    Dec,
}

impl Month {
    #[must_use]
    pub const fn index(self) -> usize {
        self as usize
    }
}

impl TryFrom<u8> for Month {
    type Error = Error;

    fn try_from(month: u8) -> Result<Self> {
        match month {
            1 => Ok(Self::Jan),
            2 => Ok(Self::Feb),
            3 => Ok(Self::Mar),
            4 => Ok(Self::Apr),
            5 => Ok(Self::May),
            6 => Ok(Self::Jun),
            7 => Ok(Self::Jul),
            8 => Ok(Self::Aug),
            9 => Ok(Self::Sep),
            10 => Ok(Self::Oct),
            11 => Ok(Self::Nov),
            12 => Ok(Self::Dec),
            _ => Err(Error::Input("month must be in 1..=12")),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Location {
    pub(crate) lat: f64,
    pub(crate) lng: f64,
}

impl Location {
    pub fn new(latitude_deg: f64, longitude_deg: f64) -> Result<Self> {
        if !latitude_deg.is_finite() || !(-90.0..=90.0).contains(&latitude_deg) {
            return Err(Error::Input(
                "latitude must be finite and in -90..=90 degrees",
            ));
        }
        if !longitude_deg.is_finite() || !(-180.0..=180.0).contains(&longitude_deg) {
            return Err(Error::Input(
                "longitude must be finite and in -180..=180 degrees",
            ));
        }
        Ok(Self {
            lat: latitude_deg,
            lng: longitude_deg,
        })
    }

    #[must_use]
    pub const fn latitude_deg(self) -> f64 {
        self.lat
    }

    #[must_use]
    pub const fn longitude_deg(self) -> f64 {
        self.lng
    }

    pub(crate) fn radians(self) -> Self {
        Self {
            lat: self.lat.to_radians(),
            lng: self.lng.to_radians(),
        }
    }

    pub(crate) const fn from_radians(latitude: f64, longitude: f64) -> Self {
        Self {
            lat: latitude,
            lng: longitude,
        }
    }
}

impl Default for Location {
    fn default() -> Self {
        Self::from_radians(0.0, 0.0)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Path {
    #[default]
    Short,
    Long,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ManMadeNoise {
    #[default]
    City,
    Residential,
    Rural,
    QuietRural,
    Noisy,
    Quiet,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DigitalConfig {
    pub sir_db: f64,
    pub f0_hz: f64,
    pub t0_ms: f64,
    pub a_db: f64,
    pub tw_ms: f64,
    pub fw_hz: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Modulation {
    #[default]
    Analog,
    Digital(DigitalConfig),
}

#[derive(Clone, Copy, Debug)]
pub struct Request<'a> {
    pub month: Month,
    /// UTC hour, with midnight represented by 0.
    pub utc_hour: u8,
    /// Twelve-month smoothed sunspot number R12, using the version 1 convention
    /// (observatory factor 0.6) specified by P.1239-4, in 0..=311.
    pub ssn: u16,
    pub modulation: Modulation,
    pub path: Path,
    pub freq_mhz: f64,
    pub bandwidth_hz: f64,
    pub power_dbkw: f64,
    pub noise: ManMadeNoise,
    pub snr_percent: u8,
    pub snr_db: f64,
    pub tx: Location,
    pub rx: Location,
    pub tx_antenna: Antenna<'a>,
    pub rx_antenna: Antenna<'a>,
}

impl Request<'_> {
    pub(crate) fn validate(&self) -> Result<()> {
        if self.utc_hour > 23 {
            return Err(Error::Input("UTC hour must be in 0..=23"));
        }
        if self.ssn > 311 {
            return Err(Error::Input("sunspot number must be in 0..=311"));
        }
        check_range(self.freq_mhz, 2.0, 30.0, "frequency must be in 2..=30 MHz")?;
        check_range(
            self.bandwidth_hz,
            0.005,
            3_000_000.0,
            "bandwidth must be in 0.005..=3000000 Hz",
        )?;
        check_range(
            self.power_dbkw,
            -30.0,
            60.0,
            "transmit power must be in -30..=60 dB(1 kW)",
        )?;
        check_range(
            self.snr_db,
            -30.0,
            200.0,
            "required SNR must be in -30..=200 dB",
        )?;
        if !(1..=99).contains(&self.snr_percent) {
            return Err(Error::Input("SNR percentile must be in 1..=99"));
        }

        if let Modulation::Digital(config) = self.modulation {
            check_range(
                config.sir_db,
                -30.0,
                200.0,
                "required SIR must be in -30..=200 dB",
            )?;
            check_range(
                config.f0_hz,
                0.0,
                1000.0,
                "frequency dispersion must be in 0..=1000 Hz",
            )?;
            check_range(
                config.t0_ms,
                0.0,
                1000.0,
                "time spread must be in 0..=1000 ms",
            )?;
            check_range(
                config.a_db,
                0.0,
                1000.0,
                "amplitude ratio must be in 0..=1000 dB",
            )?;
            check_range(config.tw_ms, 0.0, 50.0, "time window must be in 0..=50 ms")?;
            check_range(
                config.fw_hz,
                0.0,
                1000.0,
                "frequency window must be in 0..=1000 Hz",
            )?;
        }

        self.tx_antenna.validate(self.freq_mhz)?;
        self.rx_antenna.validate(self.freq_mhz)
    }
}

fn check_range(value: f64, min: f64, max: f64, message: &'static str) -> Result<()> {
    if value.is_finite() && (min..=max).contains(&value) {
        Ok(())
    } else {
        Err(Error::Input(message))
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        antenna::Antenna,
        input::{Location, ManMadeNoise, Modulation, Month, Path, Request},
    };

    fn request() -> Request<'static> {
        Request {
            month: Month::Jan,
            utc_hour: 0,
            ssn: 100,
            modulation: Modulation::Analog,
            path: Path::Short,
            freq_mhz: 10.0,
            bandwidth_hz: 3000.0,
            power_dbkw: 0.0,
            noise: ManMadeNoise::Rural,
            snr_percent: 90,
            snr_db: 10.0,
            tx: Location::new(0.0, 0.0).unwrap(),
            rx: Location::new(1.0, 1.0).unwrap(),
            tx_antenna: Antenna::isotropic(0.0),
            rx_antenna: Antenna::isotropic(0.0),
        }
    }

    #[test]
    fn month_conversion_is_calendar_based_and_fallible() {
        assert_eq!(Month::try_from(1).unwrap(), Month::Jan);
        assert_eq!(Month::try_from(12).unwrap(), Month::Dec);
        assert!(Month::try_from(0).is_err());
        assert!(Month::try_from(13).is_err());
    }

    #[test]
    fn location_rejects_nonfinite_and_out_of_range_values() {
        assert!(Location::new(f64::NAN, 0.0).is_err());
        assert!(Location::new(91.0, 0.0).is_err());
        assert!(Location::new(0.0, 181.0).is_err());
    }

    #[test]
    fn request_rejects_invalid_hours_and_nan() {
        let mut request = request();
        request.utc_hour = 24;
        assert!(request.validate().is_err());

        request = self::request();
        request.freq_mhz = f64::NAN;
        assert!(request.validate().is_err());
    }
}
