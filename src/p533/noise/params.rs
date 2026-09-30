use crate::input::{ManMadeNoise, Month};

/// Represents noise parameters including atmospheric, man-made, and galactic components.
#[derive(Default)]
pub(in crate::p533) struct NoiseParams {
    /// Atmospheric noise.
    pub fa_a: f64,
    /// Atmospheric noise upper decile.
    pub du_a: f64,
    /// Atmospheric noise lower decile.
    pub dl_a: f64,
    /// Man-made noise.
    pub fa_m: f64,
    /// Man-made noise upper decile.
    pub du_m: f64,
    /// Man-made noise lower decile.
    pub dl_m: f64,
    /// Galactic noise.
    pub fa_g: f64,
    /// Galactic noise upper decile.
    pub du_g: f64,
    /// Galactic noise lower decile.
    pub dl_g: f64,
    /// Total noise upper decile.
    pub du_t: f64,
    /// Total noise lower decile.
    pub dl_t: f64,
    /// Total noise.
    pub fam_t: f64,
    pub man_made_noise: ManMadeNoise,
    pub month: Month,
}
