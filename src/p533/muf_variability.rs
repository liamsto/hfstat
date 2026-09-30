use crate::data::DECILES;
use crate::p533::{
    consts::{D2R, MAXEMDS, MAXF2MDS},
    control_point::ControlPointId::*,
    controlpt_params::bilinear_interpolation,
    state::State,
};

/// Calculates the F2 and E Fprob, 50% MUF (`muf50`), 90% MUF, (`muf90`), and 10% MUF (`muf10`).
pub(super) const fn muf_variability(path: &mut State) {
    if path.distance > 9000.0 {
        return;
    }

    path.muf50 = path.bmuf;

    // There are 3 possible E layer modes and 6 F2 layer modes
    let mut i = 0;
    while i < MAXF2MDS {
        if let Some(mdf2_bmuf) = path.md_f2[i].bmuf {
            path.md_f2[i].muf50 = mdf2_bmuf;

            // Determine the decile factors
            // Lower decile
            path.md_f2[i].deltal = find_fof2_var(
                path,
                path.controlpt[MP as usize].ltime,
                path.controlpt[MP as usize].l.lat,
                Decile::DL,
            );
            // Upper decile
            path.md_f2[i].deltau = find_fof2_var(
                path,
                path.controlpt[MP as usize].ltime,
                path.controlpt[MP as usize].l.lat,
                Decile::DU,
            );

            // Get other MUFs
            path.md_f2[i].muf10 = path.md_f2[i].deltau * path.md_f2[i].muf50;
            path.md_f2[i].muf90 = path.md_f2[i].deltal * path.md_f2[i].muf50;

            // Now determine the probability the mode can be supported
            if path.frequency < path.md_f2[i].muf50 {
                path.md_f2[i].fprob = 1.3
                    - (0.8
                        / (1.0
                            + ((1.0 - (path.frequency / path.md_f2[i].muf50))
                                / (1.0 - path.md_f2[i].deltal))))
                        .min(1.0);
            } else {
                path.md_f2[i].fprob = ((0.8
                    / (1.0
                        + (((path.frequency / path.md_f2[i].muf50) - 1.0)
                            / (path.md_f2[i].deltau - 1.0))))
                    - 0.3)
                    .max(0.0);
            }
        }
        i += 1;
    }

    // Same loop as F2, but decile values are fixed to between 0.95 and 1.05.
    // Find the E layer MUFs
    let mut i = 0;
    while i < MAXEMDS {
        if let Some(mde_bmuf) = path.md_e[i].bmuf {
            path.md_e[i].muf50 = mde_bmuf;
            //deltal for E layer
            path.md_e[i].deltal = 0.95;
            //deltau for E layer
            path.md_e[i].deltau = 1.05;

            // Find other MUFs
            path.md_e[i].muf10 = path.md_e[i].deltau * path.md_e[i].muf50;
            path.md_e[i].muf90 = path.md_e[i].deltal * path.md_e[i].muf50;

            if path.frequency < path.md_e[i].muf50 {
                path.md_e[i].fprob = (1.3
                    - (0.8
                        / (1.0
                            + ((1.0 - (path.frequency / path.md_e[i].muf50))
                                / (1.0 - path.md_e[i].deltal)))))
                    .min(1.0);
            } else {
                path.md_e[i].fprob = ((0.8
                    / (1.0
                        + (((path.frequency / path.md_e[i].muf50) - 1.0)
                            / (path.md_e[i].deltau - 1.0))))
                    - 0.3)
                    .max(0.0);
            }
        }
        i += 1;
    }

    // As done in MUFBasic, determine the largest of the MUF90 and MUF10
    // Pick the MUF extrema as the path for the E and F2 layer.

    let mut emuf_10 = 0.0;
    let mut emuf_90 = 0.0;

    let mut i = 0;
    while i < MAXEMDS {
        if path.md_e[i].bmuf.is_some() {
            if path.md_e[i].muf90 > emuf_90 {
                emuf_90 = path.md_e[i].muf90;
            }
            if path.md_e[i].muf10 > emuf_10 {
                emuf_10 = path.md_e[i].muf10;
            }
        }
        i += 1;
    }

    let mut f2muf_10 = 0.0;
    let mut f2muf_90 = 0.0;

    let mut i = 0;
    while i < MAXF2MDS {
        if path.md_f2[i].bmuf.is_some() {
            if path.md_f2[i].muf90 > f2muf_90 {
                f2muf_90 = path.md_f2[i].muf90;
            }
            if path.md_f2[i].muf10 > f2muf_10 {
                f2muf_10 = path.md_f2[i].muf10;
            }
        }
        i += 1;
    }

    path.muf90 = emuf_90.max(f2muf_90);
    path.muf10 = emuf_10.max(f2muf_10);
}

#[must_use]
pub(super) const fn find_fof2_var(path: &State, hour: f64, lat: f64, decile: Decile) -> f64 {
    let lat = (lat.abs() / (5.0 * D2R)).clamp(0.0, 18.0);
    let hour = hour - (hour / 24.0).floor() * 24.0;
    let lat_l = lat as usize;
    let lat_u = lat.ceil() as usize;
    let hour_l = hour as usize;
    let hour_u = (hour_l + 1) % 24;

    // Sunspot index
    let ssn = match path.ssn {
        s if s < 50 => 0,
        s if s <= 100 => 1,
        _ => 2,
    };

    let table = &DECILES[decile as usize][path.season as usize][ssn];
    let ll = table[lat_l][hour_l];
    let lr = table[lat_l][hour_u];
    let ul = table[lat_u][hour_l];
    let ur = table[lat_u][hour_u];

    bilinear_interpolation(ll, lr, ul, ur, lat.fract(), hour.fract())
}

#[derive(Default, Clone, Copy)]
pub(super) enum Decile {
    #[default]
    // Lower decile
    DL,
    // Upper decile
    DU,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        input::{Location, Month},
        p533::{control_point::ControlPt, controlpt_params::solar_parameters, state::Season},
    };

    #[test]
    fn deciles_use_local_time_and_interpolate_the_reference_table() {
        let path = State {
            ssn: 20,
            season: Season::Winter,
            ..State::default()
        };
        let mut pt = ControlPt {
            l: Location::new(57.5, 22.5).unwrap().radians(),
            ..ControlPt::default()
        };
        solar_parameters(&mut pt, Month::Jan, 12.0);
        assert_eq!(pt.ltime, 13.5);
        // P.1239-4 Table 2a, latitudes 55/60 and local hours 13/14:
        // (0.82 + 0.82 + 0.80 + 0.81) / 4.
        let factor = find_fof2_var(&path, pt.ltime, pt.l.lat, Decile::DL);
        assert!((factor - 0.8125).abs() < 1e-12);
        assert_eq!(
            find_fof2_var(&path, -0.5, pt.l.lat, Decile::DL),
            find_fof2_var(&path, 23.5, pt.l.lat, Decile::DL),
        );
        assert_eq!(find_fof2_var(&path, 24.0, 90.0 * D2R, Decile::DL), 0.67);
    }
}
