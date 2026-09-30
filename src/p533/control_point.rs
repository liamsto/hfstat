use crate::{input::Location, p533::solar::SolarParameters};
use std::ops::{Index, IndexMut};

#[derive(Default, Clone, Copy)]
pub(super) struct ControlPt {
    pub(crate) l: Location,
    pub(crate) foe: f64,
    pub(crate) fof2: f64,
    pub(crate) m3kf2: f64,
    pub(crate) dip: (f64, f64),
    /// 100km and 300km gyrofrequency/magdip parameters.
    pub(crate) fh: (f64, f64),
    pub(crate) ltime: f64,
    pub(crate) hr: f64,
    pub(crate) x: f64,
    pub(crate) sun: SolarParameters,
}

#[derive(Debug, Clone, Copy)]
#[repr(usize)]
/// Control point indices for readability.
pub(super) enum ControlPointId {
    T1k,
    R1k,
    MP,
    Td02,
    Rd02,
}

impl IndexMut<ControlPointId> for [ControlPt; 5] {
    fn index_mut(&mut self, index: ControlPointId) -> &mut Self::Output {
        &mut self[index as usize]
    }
}

impl Index<ControlPointId> for [ControlPt; 5] {
    type Output = ControlPt;

    fn index(&self, index: ControlPointId) -> &Self::Output {
        &self[index as usize]
    }
}
