#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Muf {
    pub basic_mhz: f64,
    pub median_mhz: f64,
    pub p10_mhz: f64,
    pub p90_mhz: f64,
    pub operational_mhz: f64,
    pub operational_p10_mhz: f64,
    pub operational_p90_mhz: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NoiseComponent {
    pub median_db: f64,
    pub upper_decile_db: f64,
    pub lower_decile_db: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Noise {
    pub atmospheric: NoiseComponent,
    pub man_made: NoiseComponent,
    pub galactic: NoiseComponent,
    pub total: NoiseComponent,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Snr {
    pub median_db: f64,
    pub upper_decile_db: f64,
    pub lower_decile_db: f64,
    pub percentile_db: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Reliability {
    pub basic_percent: f64,
    pub overall_percent: Option<f64>,
    pub scatter_percent: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DigitalResult {
    pub sir_db: Option<f64>,
    pub sir_upper_decile_db: Option<f64>,
    pub sir_lower_decile_db: Option<f64>,
    pub snr_percent: Option<f64>,
    pub time_percent: Option<f64>,
    pub frequency_percent: Option<f64>,
    pub multimode_percent: Option<f64>,
    pub scatter_probability: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Mode {
    pub hops: u8,
    pub basic_muf_mhz: f64,
    pub median_muf_mhz: f64,
    pub muf_p10_mhz: f64,
    pub muf_p90_mhz: f64,
    pub operational_muf_mhz: f64,
    pub operational_muf_p10_mhz: f64,
    pub operational_muf_p90_mhz: f64,
    pub support_probability: f64,
    pub reflection_height_km: f64,
    pub screening_frequency_mhz: f64,
    pub field_strength_dbuv_m: f64,
    pub elevation_deg: f64,
    pub rx_power_dbw: f64,
    pub tx_gain_dbi: f64,
    pub rx_gain_dbi: f64,
    pub delay_ms: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Prediction {
    pub distance_km: f64,
    pub max_hop_km: f64,
    pub elevation_deg: f64,
    pub muf: Muf,
    pub field_strength_dbuv_m: f64,
    pub short_field_strength_dbuv_m: f64,
    pub intermediate_field_strength_dbuv_m: f64,
    pub long_field_strength_dbuv_m: f64,
    pub rx_power_dbw: f64,
    pub tx_eirp_dbkw: f64,
    pub tx_gain_dbi: f64,
    pub rx_gain_dbi: f64,
    pub noise: Noise,
    pub snr: Snr,
    pub reliability: Reliability,
    pub digital: Option<DigitalResult>,
    pub e_modes: [Option<Mode>; 3],
    pub f2_modes: [Option<Mode>; 6],
}
