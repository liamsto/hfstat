use hfstat::{
    Antenna, DigitalConfig, Location, ManMadeNoise, Modulation, Month, Orientation, Path, Pattern,
    Request,
    antenna::{ELEVATIONS, GAINS_PER_FREQUENCY},
};

pub fn gains() -> Vec<f64> {
    (0..3 * GAINS_PER_FREQUENCY)
        .map(|i| {
            let az = ((i / ELEVATIONS) % 360) as f64;
            let el = (i % ELEVATIONS) as f64;
            4.0 * az.to_radians().cos()
                + 2.0 * el.to_radians().sin()
                + (i / GAINS_PER_FREQUENCY) as f64
        })
        .collect()
}

pub fn cases(gains: &[f64]) -> Vec<(&'static str, Request<'_>)> {
    let base = Request {
        month: Month::Jan,
        utc_hour: 12,
        ssn: 100,
        modulation: Modulation::Analog,
        path: Path::Short,
        freq_mhz: 10.0,
        bandwidth_hz: 3000.0,
        power_dbkw: 0.0,
        noise: ManMadeNoise::Rural,
        snr_percent: 90,
        snr_db: 10.0,
        tx: Location::new(51.5, -0.1).unwrap(),
        rx: Location::new(40.7, -74.0).unwrap(),
        tx_antenna: Antenna::isotropic(0.0),
        rx_antenna: Antenna::isotropic(0.0),
    };
    let mut cases = Vec::new();
    for (name, lat, lng) in [
        ("local", 52.0, 0.5),
        ("regional", 41.9, 12.5),
        ("3000km", 30.0, 15.0),
        ("atlantic", 40.7, -74.0),
        ("8000km", 19.1, 72.9),
        ("12000km", -1.3, 103.8),
        ("polar", 71.3, -156.8),
    ] {
        cases.push((
            name,
            Request {
                rx: Location::new(lat, lng).unwrap(),
                ..base
            },
        ));
    }
    cases.push((
        "long",
        Request {
            path: Path::Long,
            ..base
        },
    ));
    let digital = Modulation::Digital(DigitalConfig {
        sir_db: 10.0,
        f0_hz: 10.0,
        t0_ms: 3.0,
        a_db: 6.0,
        tw_ms: 1.0,
        fw_hz: 3.0,
    });
    for (name, index) in [
        ("digital_local", 0),
        ("digital_atlantic", 3),
        ("digital_8000km", 4),
        ("digital_long", 7),
    ] {
        cases.push((
            name,
            Request {
                modulation: digital,
                ..cases[index].1
            },
        ));
    }
    let pattern = Pattern::new(&[2.0, 15.0, 30.0], gains).unwrap();
    for (name, path) in [
        ("pattern_atlantic", Path::Short),
        ("pattern_long", Path::Long),
    ] {
        cases.push((
            name,
            Request {
                path,
                tx_antenna: Antenna::pattern(pattern, Orientation::BearingDeg(123.25)),
                rx_antenna: Antenna::pattern(pattern, Orientation::Path),
                ..base
            },
        ));
    }
    cases
}
