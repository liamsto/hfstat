use std::f64::consts::PI;

use crate::p533::{
    consts::{D2R, MAXEMDS, MAXF2MDS, MINELEANGLE, R0},
    control_point::{
        ControlPointId::{MP, R1k, Rd02, T1k, Td02},
        ControlPt,
    },
    controlpt_params::calculate_ctrlpt_param,
    geometry::great_circle_point,
    ionosphere::Maps,
    screening_freq::{elevation_angle, incidence_angle},
    state::State,
};

const C_SUB_3000: f64 = 3000.0;

/// Compute and store basic MUF values for the provided propagation path.
///
/// # Errors
/// Returns an error if a computed MUF value exceeds supported bounds.
pub(super) fn muf_basic(path: &mut State, maps: &Maps) -> Result<(), &'static str> {
    // Only run the subroutine if the path is less than or equal to 9000 km.
    if path.distance > 9000.0 {
        return Ok(());
    }

    // Determine the F2 layer basic MUF
    // First, determine mirror reflection height (hr) at midpoint
    let hr = (1490.0 / path.controlpt[MP as usize].m3kf2 - 176.0).min(500.0);
    path.controlpt[MP as usize].hr = hr;

    // The mode is determined by geometry. For simplicity, assume the earth is a sphere. Find the lowest mode by geometry:
    // Mirror reflection height tells you if multiple hops are needed.
    // Determine where the horizon is - if less than half the distance, you must have multiple hops.
    // Try the other nodes with a max of 6.
    // Minimum elevation angle is MINELEANGLE degrees.

    let minele = MINELEANGLE * D2R;
    // Get hop distance to the horizon
    let aoi = incidence_angle(minele, hr);
    // Now get the arc length of the half hop length times the earth's radius. Can be no more than 4000 km.
    let dh = ((PI - aoi - ((PI / 2.0) + minele)) * R0 * 2.0).min(4000.0);

    let mut mode_idx = 0;

    while mode_idx < MAXF2MDS {
        // Is mirror reflection height horizon less than n0 hop distance?
        if dh > path.distance / (mode_idx + 1) as f64 {
            path.n0_f2 = Some(mode_idx);
            break;
        }

        mode_idx += 1;
    }

    // Check to see if a low order mode exists. If it does not, skip sections 3.5.1 and 3.5.2.
    if let Some(n0) = path.n0_f2 {
        path.dmax = calc_dmax(&mut path.controlpt[MP as usize]).min(4000.0);
        let dmax = path.dmax;

        // Lowest order mode
        if path.distance <= path.dmax {
            let intermediate = calc_b(&mut path.controlpt[MP as usize]);
            let n0f2dmuf = calc_f2dmuf(
                &mut path.controlpt[MP as usize],
                path.distance / (n0 as f64 + 1.0),
                dmax,
                intermediate,
            );
            // The lowest order F2 mode has been determined and can be loaded into the path structure.
            path.md_f2[n0].bmuf = Some(n0f2dmuf);
            // n0F2DMUF is basic F2 MUF of the path.
            path.bmuf = n0f2dmuf;
        } else {
            // Paths longer than dmax km.
            let mut frac = 1.0 / (2.0 * (n0 as f64 + 1.0));
            great_circle_point(
                path.location.0,
                path.location.1,
                &mut path.controlpt[Td02 as usize],
                path.distance,
                frac,
            );
            frac = 1.0 - (1.0 / (2.0 * (n0 as f64 + 1.0)));
            great_circle_point(
                path.location.0,
                path.location.1,
                &mut path.controlpt[Rd02 as usize],
                path.distance,
                frac,
            );

            // Find foF2, M3kF2, and foE for these control points.
            calculate_ctrlpt_param(
                &mut path.controlpt[Td02],
                path.month,
                path.hour,
                path.ssn,
                maps,
            );
            calculate_ctrlpt_param(
                &mut path.controlpt[Rd02],
                path.month,
                path.hour,
                path.ssn,
                maps,
            );

            // Get the F2 basic MUF at each control point
            // For these control points get the basic MUF at F2(dmax)MUF

            let b_td02 = calc_b(&mut path.controlpt[Td02 as usize]);
            let b_rd02 = calc_b(&mut path.controlpt[Rd02 as usize]);

            let f2dmuf_td02 = calc_f2dmuf(
                &mut path.controlpt[Td02 as usize],
                path.distance / (n0 as f64 + 1.0),
                dmax,
                b_td02,
            );
            let f2dmuf_rd02 = calc_f2dmuf(
                &mut path.controlpt[Rd02 as usize],
                path.distance / (n0 as f64 + 1.0),
                dmax,
                b_rd02,
            );

            path.md_f2[n0].bmuf = Some(f2dmuf_td02.min(f2dmuf_rd02));
            path.bmuf = f2dmuf_td02.min(f2dmuf_rd02);
        }

        if path.distance <= 9000.0 {
            for mode in n0 + 1..MAXF2MDS {
                if path.distance <= path.dmax {
                    // Paths up to dmax
                    let b = calc_b(&mut path.controlpt[MP as usize]);
                    path.md_f2[mode].bmuf = Some(calc_f2dmuf(
                        &mut path.controlpt[MP as usize],
                        path.distance / (mode as f64 + 1.0),
                        dmax,
                        b,
                    ));
                } else {
                    // Higher order modes calculated here imply the lowest order mode has already been set up by the two
                    // new control points and the basic MUF has been determined. The path.bmuf is then n0F2muf. We now
                    // find the mode basic MUF as a function of F2(dmax[controlpt[MP]])MUF and a scaling factor.
                    let td02_dmax = calc_dmax(&mut path.controlpt[Td02 as usize]);
                    let td02_b = calc_b(&mut path.controlpt[Td02 as usize]);

                    let rd02_dmax = calc_dmax(&mut path.controlpt[Rd02 as usize]);
                    let rd02_b = calc_b(&mut path.controlpt[Rd02 as usize]);

                    let mn_0_td02 = calc_f2dmuf(
                        &mut path.controlpt[Td02 as usize],
                        path.distance / (n0 as f64 + 1.0),
                        td02_dmax,
                        td02_b,
                    );
                    let mn_0_rd02 = calc_f2dmuf(
                        &mut path.controlpt[Rd02 as usize],
                        path.distance / (n0 as f64 + 1.0),
                        rd02_dmax,
                        rd02_b,
                    );

                    let mn_td02 = calc_f2dmuf(
                        &mut path.controlpt[Td02 as usize],
                        path.distance / (mode as f64 + 1.0),
                        td02_dmax,
                        td02_b,
                    );
                    let mn_rd02 = calc_f2dmuf(
                        &mut path.controlpt[Rd02 as usize],
                        path.distance / (mode as f64 + 1.0),
                        rd02_dmax,
                        rd02_b,
                    );

                    path.md_f2[mode].bmuf =
                        Some(path.bmuf * (mn_td02 / mn_0_td02).min(mn_rd02 / mn_0_rd02));
                }
            }
        }
    }

    // E layer basic MUF calculation
    if path.distance < 4000.0 {
        let hr = 110.0;
        let delta = minele;
        let aoi = incidence_angle(delta, hr);

        let dh = ((PI - aoi - ((PI / 2.0) + minele)) * R0 * 2.0).min(4000.0);

        for mode_idx in 0..MAXEMDS {
            if dh > path.distance / (mode_idx as f64 + 1.0) {
                path.n0_e = Some(mode_idx);
                break;
            }
        }

        if let Some(n0_e) = path.n0_e {
            for mode in n0_e..MAXEMDS {
                path.md_e[mode].hr = hr;
                let dh = (path.distance / (mode as f64 + 1.0)).min(4000.0);
                let delta = elevation_angle(dh, hr);
                let i110 = incidence_angle(delta, hr);

                if path.distance < 2000.0 {
                    path.md_e[mode].bmuf = Some(path.controlpt[MP as usize].foe / i110.cos());
                } else if path.distance >= 2000.0 && path.distance <= 4000.0 {
                    path.md_e[mode].bmuf = Some(
                        (path.controlpt[R1k as usize].foe).min(path.controlpt[T1k as usize].foe)
                            / i110.cos(),
                    );
                }
            }
        } else {
            // There are no E layer modes for path.distance >= 4000.0
        }
    }

    if let Some(n0_e) = path.n0_e {
        if let Some(n0_f2) = path.n0_f2 {
            // F2 and E modes exist.
            let e_bmuf = path.md_e[n0_e].bmuf.ok_or("n0_e BMUF missing")?;
            let f2_bmuf = path.md_f2[n0_f2].bmuf.ok_or("n0_f2 BMUF missing")?;
            path.bmuf = e_bmuf.max(f2_bmuf);
        } else {
            // E modes exist. F2 modes do not.
            path.bmuf = path.md_e[n0_e].bmuf.ok_or("n0_e BMUF missing")?;
        }
    } else if let Some(n0_f2) = path.n0_f2 {
        path.bmuf = path.md_f2[n0_f2].bmuf.ok_or("n0_f2 BMUF missing")?;
    } else {
        // This is an error.
        return Err("BMUF too big!");
    }

    Ok(())
}

/// Determines the intermediate values B from Eqn (6) P.533-14
pub(super) fn calc_b(ctrlpt: &mut ControlPt) -> f64 {
    if ctrlpt.foe != 0.0 {
        ctrlpt.x = f64::max(ctrlpt.fof2 / ctrlpt.foe, 2.0);
    } else {
        ctrlpt.x = 2.0;
    }

    ctrlpt.m3kf2.mul_add(ctrlpt.m3kf2, -4.0).mul_add(
        0.005f64.mul_add(((7.854 / ctrlpt.x) - 1.9635).sin(), 0.0215),
        ctrlpt.m3kf2 - 0.124,
    )
}

///Finds the `dmax` for the path. The `dmax` can be calculated to be greater than 4000 km.
///
/// `dmax` can be greater than 4000 km for the calculation of higher order MUFs.
/// all other locations it will typically be restricted to 4000 km.
pub(super) fn calc_dmax(ctrlpt: &mut ControlPt) -> f64 {
    let b = calc_b(ctrlpt);
    (12610.0 + (2140.0 / ctrlpt.x.powi(2)) - (49720.0 / ctrlpt.x.powi(4))
        + (688900.0 / ctrlpt.x.powi(6)))
    .mul_add((1.0 / b) - 0.303, 4780.0)
}

/// Determines the F2 Layer MUF from Eqn (3) P.533-14
pub(super) fn calc_f2dmuf(ctrlpt: &mut ControlPt, distance: f64, dmax: f64, b: f64) -> f64 {
    let d = if distance <= dmax { distance } else { dmax };

    let cd = calc_cd(d, dmax);

    let c3k = calc_cd(C_SUB_3000, dmax);
    (ctrlpt.fh.1 / 2.0).mul_add(
        1.0 - (distance / dmax),
        (cd / c3k).mul_add(b - 1.0, 1.0) * ctrlpt.fof2,
    )
}

#[must_use]
pub(super) fn calc_cd(d: f64, dmax: f64) -> f64 {
    let z = 1.0 - 2.0 * d / dmax;
    0.096f64.mul_add(
        z.powi(6),
        0.181f64.mul_add(
            z.powi(5),
            0.088f64.mul_add(
                z.powi(4),
                0.090f64.mul_add(
                    -z.powi(3),
                    0.424f64.mul_add(-z.powi(2), 0.591f64.mul_add(-z, 0.74)),
                ),
            ),
        ),
    )
}
