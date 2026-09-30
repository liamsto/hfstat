/// Represents the properties of a radio mode, including MUFs and other propagation parameters.
///
#[derive(Clone, Copy, Default)]
pub(super) struct Mode {
    /// Basic MUF (MHz). Typically there is no difference between the basic and the 50% MUF.
    pub(crate) bmuf: Option<f64>,
    /// MUF exceeded for 90% of the days of the month (MHz).
    pub(crate) muf90: f64,
    /// MUF exceeded for 50% of the days of the month (MHz).
    pub(crate) muf50: f64,
    /// MUF exceeded for 10% of the days of the month (MHz).
    pub(crate) muf10: f64,
    /// Operational MUF (MHz).
    pub(crate) opmuf: f64,
    /// Operational MUF exceeded 10% of the days of the month (MHz).
    pub(crate) opmuf10: f64,
    /// Operational MUF exceeded 90% of the days of the month (MHz).
    pub(crate) opmuf90: f64,
    /// Probability that the mode is supported at the frequency of interest.
    pub(crate) fprob: f64,
    /// Lower decile for the MUF calculations.
    pub(crate) deltal: f64,
    /// Upper decile for the MUF calculations.
    pub(crate) deltau: f64,
    /// Reflection height for the mode.
    pub(crate) hr: f64,
    /// E-Layer screening frequency for F2 modes only (MHz).
    pub(crate) fs: f64,
    /// < 9000 km path basic loss.
    pub(crate) lb: f64,
    /// < 9000 km field strength (dB(1 μV/m)).
    pub(crate) ew: f64,
    /// Elevation angle.
    pub(crate) ele: f64,
    /// Receiver power (dBW).
    pub(crate) prw: f64,
    /// Transmit antenna gain (dBi).
    pub(crate) gtw: f64,
    /// Receive antenna gain (dBi).
    pub(crate) grw: f64,
    /// Time delay.
    pub(crate) tau: f64,
}
