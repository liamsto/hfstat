use std::f64::consts::PI;

use crate::{
    input::{Location, Month},
    p533::{
        antenna::{Direction, antenna_gain},
        consts::{ATNO, D2R, LH, MAXEMDS, MAXF2MDS, MAXSSN, PPT, PVAL1, PVAL2, R0, R2D, TINY_DB},
        control_point::{
            ControlPointId::{MP, R1k, Rd02, T1k, Td02},
            ControlPt,
        },
        controlpt_params::{find_foe, solar_parameters, solar_time, solar_zenith},
        geometry::{geomagnetic_coordinates, great_circle_point},
        magfit::magfit,
        mode::Mode,
        screening_freq::{elevation_angle, incidence_angle},
        state::{Season, State},
    },
};

const fn idx3(i: usize, j: usize, k: usize) -> usize {
    i * 14 + j * 7 + k
}

const fn idx_lh(txrange: usize, season: usize, gmlat: usize, mplt: usize) -> usize {
    (((txrange * 3 + season) * 8 + gmlat) * 8) + mplt
}

const fn deg_limit() -> f64 {
    70.0 * D2R
}

const PEN: bool = true;
const NOIL: f64 = 8.72;

/// Calculate the median skywave field strength as described in ITU-R P.533-14 Section 5.2 paths up to 7000 km.
///
/// Compared to the original C code, this is intentionally less readable but more performant. We use `collect_modes` to merge both
/// the E and F2 layer mode arrays into a single array of size `MAXF2MDS + MAXEMDS`, which is better for performance but slightly
/// confusing. This allows us to avoid duplication of effort.
pub(super) fn median_skywave_field_strength_short(path: &mut State) {
    if path.distance > 9000.0 {
        return;
    }

    let ssn = path.ssn.min(MAXSSN);
    let hr_e = 110.0;
    let hr_f2 = if path.distance > path.dmax {
        // Calculate hr_F2 like in C
        let smallest = smallest_cp_fo_f2(path);
        (1490.0 / path.controlpt[smallest].m3kf2).min(500.0) - 176.0
    } else {
        path.controlpt[MP].hr
    };

    // Local time at midpath
    let tz = (path.controlpt[MP].l.lng / (15.0 * D2R)) as i32;
    let mpltime = ((path.controlpt[MP].ltime + f64::from(tz)) % 24.0) as i32;

    // Combined F2 + E modes
    let mut combined = collect_modes(&path.md_f2, &path.md_e);

    let mut etw = 0.0;

    // ---- F2 Layer ----
    if let Some(n0_f2) = path.n0_f2 {
        for n in n0_f2..MAXF2MDS {
            let idx = n;
            let mode = &mut combined[idx];

            let fs = mode.fs;
            let hop_length = path.distance / (n as f64 + 1.0);
            let is_first_mode = n == n0_f2;
            let valid_hop = hop_length <= path.dmax;
            let has_valid_bmuf = mode.bmuf.is_some();
            let below_fs = fs < path.frequency;

            let passes_filter =
                below_fs && ((is_first_mode && valid_hop) || (!is_first_mode && has_valid_bmuf));

            if !passes_filter {
                continue;
            }

            compute_median_mode_fields(path, mode, n, hr_f2, mpltime, f64::from(ssn), true);

            etw += 10f64.powf(mode.ew / 10.0);
        }
    }

    // ---- E Layer ----
    if let Some(n0_e) = path.n0_e {
        for n in n0_e..MAXEMDS {
            let idx = MAXF2MDS + n;
            let mode = &mut combined[idx];

            let hop_length = path.distance / (n as f64 + 1.0);

            let passes_filter =
                (n == n0_e && hop_length <= 2000.0) || (n != n0_e && mode.bmuf.is_some());

            if !passes_filter {
                continue;
            }

            compute_median_mode_fields(path, mode, n, hr_e, mpltime, f64::from(ssn), false);

            etw += 10f64.powf(mode.ew / 10.0);
        }
    }

    // ---- Final Result ----
    path.es = if etw != 0.0 {
        10.0 * etw.log10()
    } else {
        TINY_DB
    };

    path.md_f2.copy_from_slice(&combined[..MAXF2MDS]);
    path.md_e.copy_from_slice(&combined[MAXF2MDS..]);
}

fn compute_median_mode_fields(
    path: &mut State,
    mode: &mut Mode,
    n: usize,
    hr: f64,
    mpltime: i32,
    ssn: f64,
    is_f2: bool,
) {
    let delta = elevation_angle(path.distance / (n as f64 + 1.0), hr);
    mode.ele = delta;

    let aoi110 = incidence_angle(delta, if is_f2 { 110.0 } else { hr });

    let fv = path.frequency * aoi110.cos();
    let dh = path.distance / (n as f64 + 1.0);
    let psi = dh / (2.0 * R0);
    path.ptick = (2.0 * R0 * (psi.sin() / (delta + psi).cos())).abs() * (n as f64 + 1.0);

    // --- Absorption term ---
    let at = if PEN {
        penetration_points(path, n as f64, hr, fv)
    } else if path.distance <= 2000.0 {
        absorption_term(&path.controlpt[MP], path.month, fv)
    } else if path.distance <= path.dmax {
        (absorption_term(&path.controlpt[MP], path.month, fv)
            + absorption_term(&path.controlpt[T1k], path.month, fv)
            + absorption_term(&path.controlpt[R1k], path.month, fv))
            / 3.0
    } else {
        (absorption_term(&path.controlpt[MP], path.month, fv)
            + absorption_term(&path.controlpt[T1k], path.month, fv)
            + absorption_term(&path.controlpt[R1k], path.month, fv)
            + absorption_term(&path.controlpt[Td02], path.month, fv)
            + absorption_term(&path.controlpt[Rd02], path.month, fv))
            / 5.0
    };

    // --- Gyrofrequency ---
    let fl = if path.distance <= 2000.0 {
        (path.controlpt[MP].fh.0 * path.controlpt[MP].dip.0.sin()).abs()
    } else if path.distance <= path.dmax {
        ((path.controlpt[MP].fh.0 * path.controlpt[MP].dip.0.sin()).abs()
            + (path.controlpt[T1k].fh.0 * path.controlpt[T1k].dip.0.sin()).abs()
            + (path.controlpt[R1k].fh.0 * path.controlpt[R1k].dip.0.sin()).abs())
            / 3.0
    } else {
        ((path.controlpt[MP].fh.0 * path.controlpt[MP].dip.0.sin()).abs()
            + (path.controlpt[T1k].fh.0 * path.controlpt[T1k].dip.0.sin()).abs()
            + (path.controlpt[R1k].fh.0 * path.controlpt[R1k].dip.0.sin()).abs()
            + (path.controlpt[Td02].fh.0 * path.controlpt[Td02].dip.0.sin()).abs()
            + (path.controlpt[Rd02].fh.0 * path.controlpt[Rd02].dip.0.sin()).abs())
            / 5.0
    };

    // --- Other signal losses ---
    let lh = if path.distance <= 2000.0 {
        find_lh(&path.controlpt[MP], dh, mpltime, path.month)
    } else if path.distance <= path.dmax {
        (find_lh(&path.controlpt[MP], dh, mpltime, path.month)
            + find_lh(&path.controlpt[T1k], dh, mpltime, path.month)
            + find_lh(&path.controlpt[R1k], dh, mpltime, path.month))
            / 3.0
    } else {
        (find_lh(&path.controlpt[MP], dh, mpltime, path.month)
            + find_lh(&path.controlpt[T1k], dh, mpltime, path.month)
            + find_lh(&path.controlpt[R1k], dh, mpltime, path.month)
            + find_lh(&path.controlpt[Td02], dh, mpltime, path.month)
            + find_lh(&path.controlpt[Rd02], dh, mpltime, path.month))
            / 5.0
    };

    // --- Li ---
    let li = ((n as f64 + 1.0) * 0.0067f64.mul_add(ssn, 1.0) * at)
        / ((path.frequency + fl).powi(2) * aoi110.cos());

    // --- Lm ---
    if let Some(bmuf) = mode.bmuf {
        let lm = above_muf_loss(path.frequency, bmuf, is_f2);

        let lg = 2.0 * (n as f64); // C-based n+1 - 1
        let lz = NOIL;

        let lb = 20.0f64.mul_add(
            path.ptick.log10(),
            20.0f64.mul_add(path.frequency.log10(), 32.45),
        ) + li
            + lm
            + lg
            + lh
            + lz;

        let gt = antenna_gain(path, path.antenna.0, delta, Direction::TxToRx);
        let gr = antenna_gain(path, path.antenna.1, delta, Direction::RxToTx);
        let pt = path.txpower;
        let ew = 20.0f64.mul_add(path.frequency.log10(), 136.6 + pt + gt) - lb;

        // Store results back
        mode.lb = lb;
        mode.ew = ew;
        mode.gtw = gt;
        mode.grw = gr;
        mode.prw = 20.0f64.mul_add(-path.frequency.log10(), ew + gr) - 107.2;
    }
}

fn above_muf_loss(frequency: f64, basic_muf: f64, is_f2: bool) -> f64 {
    if frequency <= basic_muf {
        0.0
    } else if is_f2 {
        (36.0 * (frequency / basic_muf - 1.0).sqrt()).min(62.0)
    } else {
        (130.0 * (frequency / basic_muf - 1.0).powi(2)).min(81.0)
    }
}

fn absorption_term(pt: &ControlPt, month: Month, fv: f64) -> f64 {
    let cp = *pt;
    let p = diurnal_absorption_exponent(&cp, month);
    let chij = cp.sun.sza.min(102.0 * D2R);
    let fchij = (((0.881 * chij).cos()).powf(p)).max(0.02);

    let chijnoon = solar_zenith(cp.l, cp.sun.lsn, solar_time(month, cp.sun.lsn));
    let fchijnoon = (((0.881 * chijnoon).cos()).powf(p)).max(0.02);

    let at_noon = absorption_factor(&cp, month as u8);
    let phin = absorption_layer_penetration_factor(fv / cp.foe);

    at_noon * phin * fchij / fchijnoon
}

const fn absorption_factor(pt: &ControlPt, month: u8) -> f64 {
    let i = month_to_index(month);

    let mut x = (pt.l.lat.abs() * R2D).min(69.99) / 2.5;
    let j = x.floor() as usize;
    x -= j as f64;

    let at = ATNO[i];

    at[j].mul_add(1.0 - x, at[j + 1] * x)
}

#[must_use]
pub(super) const fn absorption_layer_penetration_factor(t: f64) -> f64 {
    let phi = if t < 0.0 {
        0.0
    } else if t <= 1.0 {
        let x = (t - 0.475) / 0.475;
        let phi = (-0.093f64)
            .mul_add(x, 0.04)
            .mul_add(x, 0.127)
            .mul_add(x, -0.027)
            .mul_add(x, 0.044)
            .mul_add(x, 0.159)
            .mul_add(x, 0.225);
        phi.min(0.53)
    } else if t <= 2.2 {
        let x = (t - 1.65) / 0.55;
        let phi = 0.043f64
            .mul_add(x, -0.07)
            .mul_add(x, -0.027)
            .mul_add(x, 0.034)
            .mul_add(x, 0.054)
            .mul_add(x, -0.049)
            .mul_add(x, 0.375);
        phi.min(0.53)
    } else if t <= 10.0 {
        0.34 + (((10.0 - t) * 0.02) / 7.8)
    } else {
        0.34
    };

    // Normalize by the base scaling factor
    phi / 0.34
}

/// Finds the value of `Lh` from Table 2 ITU-R P.533-14 - "values of Lh giving auroral and other signal losses".
#[must_use]
pub(super) fn find_lh(pt: &ControlPt, dh: f64, hour: i32, month: Month) -> f64 {
    let mut gn = Location { lat: 0.0, lng: 0.0 };

    // Geomagnetic conversion
    geomagnetic_coordinates(&pt.l, &mut gn);
    let season = what_season_for_lh(&pt.l, month);

    let txrange = if dh <= 2500.0 { 0 } else { 1 };

    let mplt = match hour {
        1..=3 => 0,
        4..=6 => 1,
        7..=9 => 2,
        10..=12 => 3,
        13..=15 => 4,
        16..=18 => 5,
        19..=21 => 6,
        _ => 7, // handles 22..=24 and 0
    };

    let lat = gn.lat.abs();
    const D2R: f64 = std::f64::consts::PI / 180.0;

    let gmlat = if lat >= 77.5 * D2R {
        0
    } else if lat >= 72.5 * D2R {
        1
    } else if lat >= 67.5 * D2R {
        2
    } else if lat >= 62.5 * D2R {
        3
    } else if lat >= 57.5 * D2R {
        4
    } else if lat >= 52.5 * D2R {
        5
    } else if lat >= 47.5 * D2R {
        6
    } else if lat >= 42.5 * D2R {
        7
    } else {
        return 0.0;
    };

    LH[idx_lh(txrange, season as usize, gmlat, mplt)]
}

const fn what_season_for_lh(location: &Location, month: Month) -> Season {
    if location.lat >= 0.0 {
        match month {
            Month::Dec | Month::Jan | Month::Feb => Season::Winter,
            Month::Mar | Month::Apr | Month::May | Month::Sep | Month::Oct | Month::Nov => {
                Season::Equinox
            }
            Month::Jun | Month::Jul | Month::Aug => Season::Summer,
        }
    } else {
        match month {
            Month::Dec | Month::Jan | Month::Feb => Season::Summer,
            Month::Mar | Month::Apr | Month::May | Month::Sep | Month::Oct | Month::Nov => {
                Season::Equinox
            }
            Month::Jun | Month::Jul | Month::Aug => Season::Winter,
        }
    }
}

/// Determines the smallest Control point foF2. The original C function skips the smallest if it is at index 0, which sounds like a bug but may be intentional. Check later.
fn smallest_cp_fo_f2(path: &State) -> usize {
    path.controlpt
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a.fof2
                .partial_cmp(&b.fof2)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(index, _)| index)
        .unwrap_or(0)
}

fn penetration_points(path: &mut State, noh: f64, hr: f64, fv: f64) -> f64 {
    let mut cpts = (ControlPt::default(), ControlPt::default());
    // Hop distance
    let dh = path.distance / (noh + 1.0);
    // Elevation angle
    let delta = elevation_angle(dh, hr);
    // Angle of incidence 90 km penetration points
    let aoi90 = incidence_angle(delta, 90.0);

    // Determine where the rays penetrate the 90 km height to calculate the dips.
    let phi = PI / 2.0 - delta - aoi90;

    // Hop distance
    let dh90 = R0 * phi;

    // Vertical incidence frequency
    let mut at_sum = 0.0;

    for i in 0..=(noh as i64) {
        let mut fracd = (i as f64).mul_add(dh, dh90) / path.distance;
        great_circle_point(
            path.location.0,
            path.location.1,
            &mut cpts.0,
            path.distance,
            fracd,
        );
        // Absorption uses the 100 km dip, Sun and foE, not the F2 maps.
        magfit::<100>(&mut cpts.0);
        solar_parameters(&mut cpts.0, path.month, path.hour as f64);
        find_foe(&mut cpts.0, path.month, path.hour, path.ssn);

        at_sum += absorption_term(&cpts.0, path.month, fv);

        fracd = (i as f64 + 1.0).mul_add(dh, -dh90) / path.distance;

        great_circle_point(
            path.location.0,
            path.location.1,
            &mut cpts.1,
            path.distance,
            fracd,
        );
        magfit::<100>(&mut cpts.1);
        solar_parameters(&mut cpts.1, path.month, path.hour as f64);
        find_foe(&mut cpts.1, path.month, path.hour, path.ssn);

        at_sum += absorption_term(&cpts.1, path.month, fv);
    }

    at_sum / (2.0 * (noh + 1.0))
}

/// Calculates the diural absorption exponent, a function of the month and magnetic dip angle.
fn diurnal_absorption_exponent(cp: &ControlPt, month: Month) -> f64 {
    let mut p = 0.0;

    let mut month = month.index() as i32;

    let mut moddip = cp.dip.0.atan2(cp.l.lat.cos().sqrt()).abs();
    if moddip > deg_limit() {
        moddip = deg_limit();
    }

    if cp.l.lat < 0.0 {
        month += 6;
        if month > 11 {
            month -= 12;
        }
    }

    let pp = PPT[month as usize] * D2R;
    let (i, moddip) = if moddip > pp {
        let norm = 2.0 * (moddip - pp) / (deg_limit() - pp) - 1.0;
        (1, norm)
    } else {
        let norm = 2.0 * moddip / pp - 1.0;
        (0, norm)
    };

    let mut sx = 1.0;

    let mut j = 0;
    while j < 7 {
        let a = if month <= 5 {
            PVAL1[idx3(month as usize, i, j)]
        } else {
            PVAL2[idx3(month as usize - 6, i, j)]
        };

        p += a * sx;
        sx *= moddip;

        j += 1;
    }

    p
}

/// Merges the two mode layer fields, F2 and E, into one statically sized array. By using a merged array, we can avoid some of the redundancy that occurs in the original C code.
fn collect_modes(
    f2_modes: &[Mode; MAXF2MDS],
    e_modes: &[Mode; MAXEMDS],
) -> [Mode; MAXF2MDS + MAXEMDS] {
    let mut combined = [Mode::default(); MAXF2MDS + MAXEMDS];
    combined[..MAXF2MDS].copy_from_slice(f2_modes);
    combined[MAXF2MDS..].copy_from_slice(e_modes);
    combined
}

/// Converts calendar month (1–12) to the internal index used in ATNO
const fn month_to_index(month: u8) -> usize {
    match month {
        0 | 11 => 0, // Jan or Dec
        1 => 1,
        2 => 2,
        3 => 3,
        4 | 7 => 4, // May or Aug
        5 | 6 => 5, // Jun or Jul
        8 => 6,
        9 => 7,
        10 => 8,
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use crate::p533::skywave_short::above_muf_loss;

    #[test]
    fn above_muf_loss_uses_p533_14_equations() {
        assert_eq!(above_muf_loss(10.0, 10.0, true), 0.0);
        assert!((above_muf_loss(12.5, 10.0, true) - 18.0).abs() < 1e-12);
        assert!((above_muf_loss(12.0, 10.0, false) - 5.2).abs() < 1e-12);
        assert_eq!(above_muf_loss(30.0, 2.0, true), 62.0);
        assert_eq!(above_muf_loss(30.0, 2.0, false), 81.0);
    }
}
