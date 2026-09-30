use crate::p533::{
    consts::{MAXF2MDS, R0},
    control_point::{ControlPointId::*, ControlPt},
    state::State,
};

/// Calculates the E-layer screening frequency for the F2 modes.
pub(super) fn e_layer_screening_frequency(path: &mut State) {
    // E layer scrrening is restricted to paths shorter than 4000 km.
    if path.distance > 4000.0 {
        return;
    }

    let foe = if path.distance <= 2000.0 {
        path.controlpt[MP as usize].foe
    } else {
        (path.controlpt[T1k as usize].foe).max(path.controlpt[R1k as usize].foe)
    };

    let Some(n0_f2) = path.n0_f2 else {
        return;
    };

    for i in n0_f2..MAXF2MDS {
        let dh = path.distance / ((i + 1) as f64);
        if path.distance <= path.dmax {
            path.md_f2[i].hr = mirror_reflection_height(path, &path.controlpt[MP as usize], dh);
        } else if path.distance > path.dmax {
            // In this case, the mean of the mirror reflection heights of all control points must be taken.
            // Assume the hop distance is path.dmax
            path.md_f2[i].hr = (mirror_reflection_height(path, &path.controlpt[Td02 as usize], dh)
                + mirror_reflection_height(path, &path.controlpt[MP as usize], dh)
                + mirror_reflection_height(path, &path.controlpt[Rd02 as usize], dh))
                / 3.0;
        }

        let deltaf = elevation_angle(dh, path.md_f2[i].hr);
        let incidence = incidence_angle(deltaf, 110.0);

        path.md_f2[i].fs = (1.05 * foe) / incidence.cos();
    }
}

#[must_use]
pub(super) fn elevation_angle(dh: f64, hr: f64) -> f64 {
    ((dh / (2.0 * R0)).tan().recip() - (R0 / (R0 + hr)) / (dh / (2.0 * R0)).sin()).atan()
}

#[must_use]
pub(super) fn incidence_angle(deltaf: f64, hr: f64) -> f64 {
    (R0 * deltaf.cos() / (R0 + hr)).asin()
}

fn mirror_reflection_height(path: &State, controlpt: &ControlPt, d: f64) -> f64 {
    let hr;

    // Calculate critical frequency ratio
    let x = controlpt.fof2 / controlpt.foe;
    let y = x.max(1.8);
    let delta_m = 0.096f64.mul_add(
        (f64::from(path.ssn).min(160.0) - 25.0) / 150.0,
        0.18 / (y - 1.4),
    );
    let xr = path.frequency / controlpt.fof2;

    let h = 1490.0 / (controlpt.m3kf2 + delta_m) - 316.0;

    if x > 3.33 {
        if xr >= 1.0 {
            // Case a)
            let e1 =
                0.7506f64.mul_add(-xr, 0.6870f64.mul_add(xr.powi(2), -0.09707 * xr.powi(3))) + 0.6;

            let f1 = if xr <= 1.71 {
                33.50f64.mul_add(
                    xr,
                    32.03f64.mul_add(
                        -xr.powi(2),
                        12.95f64.mul_add(xr.powi(3), -1.862 * xr.powi(4)),
                    ),
                ) - 10.91
            } else {
                0.2f64.mul_add(xr, 1.21)
            };

            let g = if xr <= 3.7 {
                63.15f64.mul_add(
                    -xr.powi(2),
                    19.50f64.mul_add(xr.powi(3), -2.102 * xr.powi(4)),
                ) - 44.73
            } else {
                19.25
            };

            let ds = (h + 43.0).mul_add(g, 160.0);
            let a = (d - ds) / (h + 140.0);

            let a1 = (h - 47.0).mul_add(e1, 140.0);
            let b1 = (h - 17.0).mul_add(f1, 150.0) - a1;

            let h_val = if b1 >= 0.0 && a >= 0.0 {
                a1 + b1 * 2.4_f64.powf(-a)
            } else {
                a1 + b1
            };

            hr = h_val.min(800.0);
        } else {
            // Case b)
            let z = xr.max(0.1);
            let e2 = 0.00583f64.mul_add(z, 0.1906 * z.powi(2)) + 0.1936;
            let a2 = (h - 47.0).mul_add(e2, 151.0);

            let f2 = 0.883f64.mul_add(z, 0.645 * z.powi(2)) + 0.162;
            let b2 = (h - 24.0).mul_add(f2, 141.0) - a2;

            let df = (0.115 * d) / (z * (h + 140.0)).min(0.65);
            let b = 0.378f64.mul_add(
                -df,
                8.834f64.mul_add(
                    -df.powi(2),
                    15.75f64.mul_add(df.powi(3), -7.535 * df.powi(4)),
                ),
            ) + 1.0;

            let h_val = if b2 >= 0.0 { a2 + b2 * b } else { a2 + b2 };

            hr = h_val.min(800.0);
        }
    } else {
        // Case c)
        let j = 16.13f64.mul_add(-y, 5.863f64.mul_add(y.powi(2), -0.7126 * y.powi(3)));
        let u = (1.2e-3 * h).mul_add(
            y.powf(-3.6),
            8.0e-5 * (h - 80.0) * 11.0f64.mul_add(y.powf(-2.2), 1.0),
        );
        hr = (115.0 + h * j + u * d).min(800.0);
    }

    hr
}
