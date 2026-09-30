# hfstat

Rust library for predicting HF circuit performance. Designed to be small, self-contained, and high-performance. The Rust standard library is the only dependency. 

Based on ITU-R P.533-14, P.1239-4 and P.372-17. Call `predict(request)` to get MUF, field strength, received power, noise, SNR and reliability. Prediction does no I/O and uses immutable embedded reference data and stack storage. 

Portions of this implementation were derived from the ITU-R Study Group 3 ITU-R-HF reference implementation. The ITU-R-HF project states that the ITURHFProp, P533 and P372 software may be used by implementers in implementations of the relevant Recommendation free from copyright assertions. See [NOTICE](NOTICE) for the upstream copyright notice, permission and disclaimer.

## Basic Usage

```rust
use hfstat::{Antenna, Location, ManMadeNoise, Modulation, Month, Path, Request, predict};

fn main() -> hfstat::Result<()> {
    let result = predict(Request {
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
        tx: Location::new(51.5, -0.1)?,
        rx: Location::new(40.7, -74.0)?,
        tx_antenna: Antenna::isotropic(0.0),
        rx_antenna: Antenna::isotropic(0.0),
    })?;
    assert!(result.snr.median_db.is_finite());
    Ok(())
}
```

Replace parameters as needed. It is a good idea to read the standards to understand what the parameters and outputs mean.

## Notes

- `utc_hour` uses 0 for midnight. 
- `ssn` is the twelve-month smoothed R12 sunspot number, values from 0 to 311 are accepted. 
- Antenna patterns are optional, not needed for an isotropic antenna

foF2 and M(3000)F2 use the calculations found in P.1239-4 Annex 1, equations (1)-(4) and Table 1, using the specified 1960 magnetic field at 300 km. All twelve months of coefficients occupy 268 KiB as `f64`, replacing 128 MiB of sampled `f32` grids. Noise coefficients and MUF variability factors are also embedded, they're about 638 KiB total. The coefficients come from the [ITU-R-HF data files](https://github.com/ITU-R-Study-Group-3/ITU-R-HF/tree/82017594a1c6cacfaa7e86954c4ae7b3a5825a3d/P533/Data), at the pinned commit shown in that link. `src/data.rs` records the input SHA-256 hashes.

To regenerate the arrays, download `COEFF01W.txt` through `COEFF12W.txt` and `P1239-3 Decile Factors.txt` from that commit into a directory, then run `python3 tools/import_data.py DIRECTORY`. 

`cargo test` runs a couple reference comparisons, solar checks, variability-table checks, and predictions across all months on short and long paths. The grid comparisons allow 0.0001 MHz for foF2 and 0.000002 for M(3000)F2, which are tolerances that I just made up. Run `cargo bench --bench predict` for 14 end-to-end prediction workloads. Set `HFSTAT_BENCH_FILTER=long` to run a subset.
