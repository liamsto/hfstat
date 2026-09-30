mod common;

use hfstat::{Location, ManMadeNoise, Month, Request, predict};

// Don't touch
#[test]
fn matches_baseline() {
    let gains = common::gains();
    let mut requests = Vec::new();
    for (_, base) in common::cases(&gains) {
        for month in 1..=12 {
            for phase in 0..2 {
                let i = (month as usize - 1) * 2 + phase;
                requests.push(Request {
                    month: Month::try_from(month).unwrap(),
                    utc_hour: i as u8,
                    ssn: [0, 50, 100, 160, 200, 311][i % 6],
                    freq_mhz: [2.0, 7.0, 10.0, 14.0, 21.0, 30.0][(i / 2) % 6],
                    noise: [
                        ManMadeNoise::City,
                        ManMadeNoise::Residential,
                        ManMadeNoise::Rural,
                        ManMadeNoise::QuietRural,
                        ManMadeNoise::Noisy,
                        ManMadeNoise::Quiet,
                    ][i % 6],
                    snr_percent: [1, 10, 50, 90, 99][i % 5],
                    ..base
                });
            }
        }
    }
    let base = requests[0];
    for distance in [
        0.001, 1999.999, 2000.001, 3999.999, 4000.001, 6999.999, 7000.001, 8999.999, 9000.001,
        19999.0,
    ] {
        requests.push(Request {
            tx: Location::new(0.0, 0.0).unwrap(),
            rx: Location::new(0.0, (distance / 6371.009_f64).to_degrees()).unwrap(),
            ..base
        });
    }
    for (lat1, lng1, lat2, lng2) in [
        (90.0, 0.0, 70.0, 90.0),
        (-90.0, 0.0, -70.0, -90.0),
        (0.0, 179.0, 10.0, -179.0),
        (-33.9, 151.2, -33.9, 18.4),
    ] {
        requests.push(Request {
            tx: Location::new(lat1, lng1).unwrap(),
            rx: Location::new(lat2, lng2).unwrap(),
            ..base
        });
    }
    let actual: Vec<_> = requests
        .iter()
        .map(|&r| format!("{:?}", predict(r)))
        .collect();
    if let Ok(path) = std::env::var("HFSTAT_CAPTURE") {
        std::fs::write(path, actual.join("\n") + "\n").unwrap();
        return;
    }
    let expected = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/predictions.txt"
    ))
    .unwrap();
    assert_eq!(actual.len(), expected.lines().count());
    let tokens = |line: &str| {
        line.split(|c: char| c.is_ascii_whitespace() || "{}[](),:".contains(c))
            .filter(|s| !s.is_empty())
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let mut max_error = 0.0_f64;
    let mut max_scaled = 0.0_f64;
    for (i, (actual, expected)) in actual.iter().zip(expected.lines()).enumerate() {
        let actual = tokens(actual);
        let expected = tokens(expected);
        assert_eq!(actual.len(), expected.len(), "case {i}: {:?}", requests[i]);
        for (a, e) in actual.iter().zip(&expected) {
            if a == e {
                continue;
            }
            match (a.parse::<f64>(), e.parse::<f64>()) {
                (Ok(a), Ok(e)) => {
                    let error = (a - e).abs();
                    max_error = max_error.max(error);
                    max_scaled = max_scaled.max(error / e.abs().max(1.0));
                    assert!(
                        a.is_finite() && e.is_finite() && error <= 1e-11 + 1e-12 * e.abs(),
                        "case {i}: {a} != {e}, error={error:e}; {:?}",
                        requests[i]
                    );
                }
                _ => assert_eq!(a, e, "case {i}: {:?}", requests[i]),
            }
        }
    }
    eprintln!(
        "{} predictions; maximum absolute difference: {max_error:e}; normalized: {max_scaled:e}",
        requests.len()
    );
}
