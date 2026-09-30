#[derive(Default)]
pub(in crate::p533) struct FamStats {
    /// Timeblock
    pub(in crate::p533) tmblk: i32,
    /// Atmospheric noise in dB above KT0b at 1 MHz
    pub(in crate::p533) fa: f64,
    /// Standard deviation of Fam values
    pub(in crate::p533) sigma_fam: f64,
    /// Radio of upper decile to median Fam value
    pub(in crate::p533) du: f64,
    /// Standard deviations of values of Du
    pub(in crate::p533) sigma_du: f64,
    /// Ratio of median value, Fam, to lower decile
    pub(in crate::p533) dl: f64,
    /// Standard deviation of values of Dl
    pub(in crate::p533) sigma_dl: f64,
}
