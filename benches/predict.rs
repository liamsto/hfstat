#[path = "../tests/common/mod.rs"]
mod common;

use hfstat::{Month, predict};
use std::{
    hint::black_box,
    time::{Duration, Instant},
};

fn main() {
    let gains = common::gains();
    let filter = std::env::var("HFSTAT_BENCH_FILTER").unwrap_or_default();
    println!("case,median_ns,p10_ns,p90_ns,iterations");
    for (name, base) in common::cases(&gains) {
        if !name.contains(&filter) {
            continue;
        }
        let requests: Vec<_> = (0..24)
            .map(|hour| hfstat::Request {
                month: Month::try_from(hour / 2 + 1).unwrap(),
                utc_hour: hour,
                ssn: [0, 50, 100, 160, 200, 311][hour as usize % 6],
                freq_mhz: [3.5, 7.0, 10.0, 14.0, 21.0, 28.0][hour as usize % 6],
                ..base
            })
            .collect();
        for &request in &requests {
            predict(request).expect(name);
        }
        let run = |iterations: usize| {
            let start = Instant::now();
            for i in 0..iterations {
                black_box(predict(black_box(requests[i % requests.len()]))).unwrap();
            }
            start.elapsed()
        };
        let mut iterations = requests.len();
        while run(iterations) < Duration::from_millis(30) {
            iterations *= 2;
        }
        let mut samples: Vec<_> = (0..21)
            .map(|_| run(iterations).as_nanos() as f64 / iterations as f64)
            .collect();
        samples.sort_by(f64::total_cmp);
        println!(
            "{name},{:.1},{:.1},{:.1},{iterations}",
            samples[10], samples[2], samples[18]
        );
    }
}
