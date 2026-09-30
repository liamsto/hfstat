use std::f64::consts::PI;

use crate::{
    input::{Month, Path},
    p533::{
        antenna::{Direction, antenna_gain_08},
        consts::{D2R, MINELEANGLE, R0, R2D, WINTER_ANOMALY},
        control_point::{ControlPointId::*, ControlPt},
        controlpt_params::{calculate_ctrlpt_param, solar_parameters, solar_time, solar_zenith},
        geometry::{bearing, great_circle_point},
        ionosphere::{Maps, basis},
        magfit::magfit,
        muf_variability::find_fof2_var,
        screening_freq::{elevation_angle, incidence_angle},
        state::State,
    },
};

use crate::p533::muf_variability::Decile::{DL, DU};

const HOURS: usize = 24;

const TX: usize = 0;
const RX: usize = 1;
const NOIL_LONG: f64 = -0.14;

// Interestingly, inlining this 24-hour model slowed short paths.
#[inline(never)]
pub(super) fn median_skywave_field_strength_long(path: &mut State, maps: &Maps) {
    if path.distance >= 7000.0 {
        let mut ends = [[ControlPt::default(); HOURS]; 2];
        let mut sum_cos_chi = [0.0; HOURS];
        let mut tx_1k = ControlPt::default();
        let mut rx_1k = ControlPt::default();
        let hr = 300.0;

        // fL Calculation
        let mut n_l = 0;
        while path.distance / (n_l as f64 + 1.0) > 3000.0 {
            n_l += 1;
        }

        let d_l = path.distance / (n_l as f64 + 1.0);
        let delta_l = elevation_angle(d_l, hr);

        // fM Calculation
        let mut n_m = 0;
        while path.distance / (f64::from(n_m) + 1.0) > 4000.0 {
            n_m += 1;
        }

        let mut d_m = path.distance / (f64::from(n_m) + 1.0);
        let mut delta_m = elevation_angle(d_m, hr);

        if delta_m < MINELEANGLE.to_radians() {
            n_m += 1;
            d_m = path.distance / (f64::from(n_m) + 1.0);
            delta_m = elevation_angle(d_m, hr);
        }

        let i90 = incidence_angle(delta_l, 90.0);
        let phi = PI / 2.0 - delta_l - i90;
        let dh90 = R0 * phi;

        let sun: [_; HOURS] = core::array::from_fn(|hour| solar_time(path.month, hour as f64));

        // Locations and magnetic fields do not depend on the hour. Keep each
        // hour's start/end summation order unchanged when interchanging loops.
        for i in 0..=n_l {
            let mut start = ControlPt::default();
            let mut end = ControlPt::default();
            great_circle_point(
                path.location.0,
                path.location.1,
                &mut start,
                path.distance,
                (i as f64).mul_add(d_l, dh90) / path.distance,
            );
            great_circle_point(
                path.location.0,
                path.location.1,
                &mut end,
                path.distance,
                ((i + 1) as f64).mul_add(d_l, -dh90) / path.distance,
            );
            for (hour, sum) in sum_cos_chi.iter_mut().enumerate() {
                add_sun(sum, solar_zenith(start.l, hour as f64, sun[hour]));
                add_sun(sum, solar_zenith(end.l, hour as f64, sun[hour]));
            }
            if i == 0 || i == n_l {
                calculate_ctrlpt_param(&mut start, path.month, path.hour, path.ssn, maps);
                if i == 0 {
                    tx_1k = start;
                }
                if i == n_l {
                    rx_1k = start;
                }
            }
        }

        let fracd_tx = 1.0 / (2.0 * (f64::from(n_m) + 1.0));
        for (points, fraction) in ends.iter_mut().zip([fracd_tx, 1.0 - fracd_tx]) {
            let mut point = ControlPt::default();
            great_circle_point(
                path.location.0,
                path.location.1,
                &mut point,
                path.distance,
                fraction,
            );
            magfit::<100>(&mut point);
            magfit::<300>(&mut point);
            point.hr = 300.0;
            points.fill(point);
        }

        let bases = [basis(&ends[TX][0]), basis(&ends[RX][0])];
        for hour in 0..HOURS {
            let maps = Maps::new(path.month, hour, path.ssn);
            for (points, basis) in ends.iter_mut().zip(&bases) {
                let point = &mut points[hour];
                (point.fof2, point.m3kf2) = maps.at(basis);
                solar_parameters(point, path.month, hour as f64);
            }
        }

        let psi = d_m / (2.0 * R0);
        let tick = 2.0 * R0 * (psi.sin() / (delta_m + psi).cos()) * (f64::from(n_m) + 1.0);
        path.ptick = tick.abs();

        path.e0 = 20.0f64.mul_add(-path.ptick.log10(), 139.6);
        let mut tx_elevation = 0.0;
        path.gtl = antenna_gain_08(path, path.antenna.0, Direction::TxToRx, &mut tx_elevation);

        let d = path.distance;
        path.gap = (10.0 * (d / (R0 * (d / R0).sin().abs())).log10()).min(15.0);
        path.ly = NOIL_LONG;
        let pt = path.txpower;

        path.fh = (ends[TX][path.hour].fh.1 + ends[RX][path.hour].fh.1) / 2.0;

        find_mufs_and_fm(path, &ends, d_m);
        find_fl(path, &sum_cos_chi, path.ptick, path.fh, i90);

        let f = path.frequency;
        let flh = path.fl + path.fh;
        let fmh = path.fm + path.fh;
        let fh = path.fh;

        let etl = ((flh * flh) / ((f + fh).powi(2)) + ((f + fh).powi(2)) / (fmh * fmh))
            * (fmh * fmh / flh.mul_add(flh, fmh * fmh));
        let mut etl = path.e0 * (1.0 - etl);
        etl = etl - 30.0 + pt + path.gtl + path.gap - path.ly;
        path.el = etl;

        if path.distance > 9000.0 {
            path.ele = delta_m;

            path.controlpt[Rd02] = ends[RX][path.hour];
            path.controlpt[Td02] = ends[TX][path.hour];
            path.controlpt[T1k] = tx_1k;
            path.controlpt[R1k] = rx_1k;

            path.dmax = 4000.0;
        }
    }
}

fn add_sun(sum: &mut f64, zenith: f64) {
    if zenith > 0.0 && zenith < PI / 2.0 {
        *sum += zenith.cos().sqrt();
    }
}

fn find_fl(path: &mut State, sum_cos_chi: &[f64; HOURS], ptick: f64, fh: f64, i90: f64) {
    let mut fl = [0.0f64; 24];

    let month = path.month;

    // Precompute cosine of i90
    let cos_i90 = i90.cos();

    // === Compute winter anomaly factor at midpoint CP ===
    let aw = winter_anomaly(path.controlpt[MP].l.lat, month);

    let fln = (path.distance / 3000.0).sqrt();

    // === Calculate 24 values of fL ===
    for i in 0..24 {
        let sum_cos = sum_cos_chi[i];
        let numerator = 0.009f64.mul_add(f64::from(path.ssn), 1.0) * sum_cos;
        let denom = cos_i90 * (9.5e6 / ptick).ln();

        let base = 5.3f64.mul_add((numerator / denom).sqrt(), -fh) * (aw + 1.0);
        fl[i] = base.max(fln);
    }

    // === Find transition time from day to night ===
    let mut tr = None;

    for now in 0..24 {
        let prev = (now + 23) % 24; // equivalent to (now - 1 + 24) % 24

        if tr.is_none() && fl[prev] >= 2.0 * fln && fl[now] <= 2.0 * fln {
            tr = Some(now);

            if let Some(tr) = tr {
                let dt = 2.0f64.mul_add(fln, -fl[tr]) / (fl[prev] - fl[tr]);
                fl[tr] = 0.7945 * fl[prev] * (dt * (1.0 - 0.7945) + 0.7945);

                if fl[now] < fl[tr] {
                    fl[now] = fl[tr];
                }
            }
        }
    }

    // === Apply smoothing over next 3 hours if transition time was found ===
    if let Some(tr) = tr {
        for i in 1..4 {
            let now = (tr + i) % 24;
            let prev = (now + 23) % 24;

            fl[now] = (fl[prev] * 0.7945).max(fl[now]);
        }
    }

    path.fl = fl[(path.hour + 1) % HOURS];
}

fn find_mufs_and_fm(path: &mut State, points: &[[ControlPt; HOURS]; 2], d_m: f64) {
    let mut f_bm = [[0.0; 24]; 2];
    let mut f_bm_min: [f64; 2] = [100.0, 100.0];

    let f_d = (-2.40074637494790e-24f64)
        .mul_add(d_m, 25.8520201885984e-21)
        .mul_add(d_m, -92.4986988833091e-18)
        .mul_add(d_m, 102.342990689362e-15)
        .mul_add(d_m, 22.0776941764705e-12)
        .mul_add(d_m, 87.4376851991085e-9)
        .mul_add(d_m, 29.1996868566837e-6)
        * d_m;

    // Calculate noon index for each control point (0 = TdM2, 1 = RdM2)
    let noon = [
        ((12.0 - points[TX][1].l.lng / (15.0 * D2R)).floor() as i32 - 1).rem_euclid(24) as usize,
        ((12.0 - points[RX][1].l.lng / (15.0 * D2R)).floor() as i32 - 1).rem_euclid(24) as usize,
    ];

    for t in 0..24 {
        for i in 0..2 {
            let point = &points[i][t];
            let f4 = 1.1 * point.fof2 * point.m3kf2;
            let fz = 0.5f64.mul_add(point.fh.1, point.fof2);
            let f_bm_val = (f4 - fz).mul_add(f_d, fz);
            // println!("Storing f_bm value: {} = {} + ({} - {}) * {}", f_bm_val, fz, f4, fz, f_d);

            f_bm[i][t] = f_bm_val;
            f_bm_min[i] = f_bm_min[i].min(f_bm_val);
        }
    }

    // Forward azimuth at midpoint
    let mut a = bearing(&path.controlpt[MP].l, &path.location.1, Path::Short);
    if a > PI {
        a -= PI;
    }
    if a >= PI / 2.0 {
        a -= PI / 2.0;
    } else {
        a = PI / 2.0 - a;
    }

    let ew = a / (PI / 2.0);
    let iw = (1.0 - ew).mul_add(0.1, ew * 0.2);
    let ix = (1.0 - ew).mul_add(1.2, ew * 0.2);
    let iy = (1.0 - ew).mul_add(0.6, ew * 0.4);

    let (mut k0, mut k1) = path.k;

    for i in 0..2 {
        let hour_val = f_bm[i][path.hour];
        let noon_val = f_bm[i][noon[i]];
        let min_val = f_bm_min[i];

        let k_val = 1.2
            + iw * (hour_val / noon_val)
            + ix * ((noon_val / hour_val).cbrt() - 1.0)
            + iy * (min_val / noon_val).powi(2);

        match i {
            0 => k0 = k_val,
            1 => k1 = k_val,
            _ => unreachable!(),
        }
    }

    path.k = (k0, k1);

    let f_m0 = path.k.0 * f_bm[0][path.hour];
    let f_m1 = path.k.1 * f_bm[1][path.hour];
    path.fm = f_m0.min(f_m1);

    if path.distance > 9000.0 {
        // Determine which control point has the smaller basic MUF at current hour
        let (smaller_cp_idx, b_muf) = if f_bm[0][path.hour] < f_bm[1][path.hour] {
            (TX, f_bm[0][path.hour])
        } else {
            (RX, f_bm[1][path.hour])
        };

        path.bmuf = b_muf;

        let cp_point = &points[smaller_cp_idx][path.hour];
        let decile_l = find_fof2_var(path, cp_point.ltime, cp_point.l.lat, DL);
        let decile_u = find_fof2_var(path, cp_point.ltime, cp_point.l.lat, DU);

        path.muf50 = path.bmuf;
        path.muf10 = path.muf50 * decile_u;
        path.muf90 = path.muf50 * decile_l;

        path.opmuf = path.fm;
        path.opmuf10 = path.opmuf * decile_u;
        path.opmuf90 = path.opmuf * decile_l;
    }
}

/// Calculates the winter anomaly factor for long paths ( > 9000 km). Values other than 60 degrees latitude are determined by interpolation.
fn winter_anomaly(lat: f64, month: Month) -> f64 {
    let cardinality = usize::from(lat < 0.0);
    let f_lat = lat.abs();
    if (0.0..=30.0 * D2R).contains(&f_lat) || (f_lat >= 90.0 * D2R) {
        0.0
    } else if f_lat < 60.0 * D2R {
        WINTER_ANOMALY[month as usize][cardinality] * f_lat.mul_add(R2D, -30.0) / 30.0
    } else {
        WINTER_ANOMALY[month as usize][cardinality] * f_lat.mul_add(-R2D, 90.0) / 30.0
    }
}
