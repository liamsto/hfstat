use hfstat::{Antenna, Location, ManMadeNoise, Modulation, Month, Path, Request, predict};

#[test]
fn predicts_with_embedded_data() {
    let mut request = Request {
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
    for month in 1..=12 {
        request.month = Month::try_from(month).unwrap();
        for hour in [0, 6, 12, 18, 23] {
            request.utc_hour = hour;
            for ssn in [0, 50, 100, 160, 200] {
                request.ssn = ssn;
                for path in [Path::Short, Path::Long] {
                    request.path = path;
                    let result = predict(request).unwrap_or_else(|error| {
                        panic!("month={month}, hour={hour}, ssn={ssn}, path={path:?}: {error}")
                    });
                    assert!(result.rx_power_dbw.is_finite());
                    assert!(result.noise.total.median_db.is_finite());
                    assert!((0.0..=100.0).contains(&result.reliability.basic_percent));
                    assert_eq!(result, predict(request).unwrap());
                }
            }
        }
    }
    request.ssn = 312;
    assert!(predict(request).is_err());
}
