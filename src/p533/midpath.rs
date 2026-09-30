use crate::p533::control_point::ControlPointId::{Rd02, Td02};
use crate::p533::muf_basic::{calc_dmax, calc_f2dmuf};
use crate::{p533::muf_basic::calc_b, p533::state::State};

/// Interpolates between the path-distances of 7000 and 9000 km.
/// The interpolation method is given in ITU-R P.533-14 Section 5.4 "Paths between 7000 and 9000 km"
pub(super) fn midpath(path: &mut State) {
    if !((7000.0 < path.distance) && (path.distance < 9000.0)) {
        return;
    }

    let Some(n0_f2) = path.n0_f2 else {
        return;
    };

    let xl = 10.0_f64.powf(path.el / 100.0);
    let xs = 10.0_f64.powf(path.es / 100.0);
    let xi = ((path.distance - 7000.0) / 2000.0).mul_add(xl - xs, xs);
    path.ei = 100.0 * xi.log10();

    let b_td02 = calc_b(&mut path.controlpt[Td02 as usize]);
    let dmax_td02 = calc_dmax(&mut path.controlpt[Td02 as usize]).min(4000.0);

    let b_rd02 = calc_b(&mut path.controlpt[Rd02 as usize]);
    let dmax_rd02 = calc_dmax(&mut path.controlpt[Rd02 as usize]).min(4000.0);

    let bmuf1 = calc_f2dmuf(
        &mut path.controlpt[Td02 as usize],
        path.distance / (n0_f2 as f64 + 1.0),
        dmax_td02,
        b_td02,
    );
    let bmuf2 = calc_f2dmuf(
        &mut path.controlpt[Rd02 as usize],
        path.distance / (n0_f2 as f64 + 1.0),
        dmax_rd02,
        b_rd02,
    );

    path.bmuf = bmuf1.min(bmuf2);
}
