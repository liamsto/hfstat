#[derive(Default, Clone, Copy)]
pub(super) struct SolarParameters {
    pub(crate) sza: f64,
    pub(crate) decl: f64,
    pub(crate) lsr: f64,
    pub(crate) lsn: f64,
    pub(crate) lss: f64,
}
