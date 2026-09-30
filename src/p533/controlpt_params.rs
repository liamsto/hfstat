use crate::{
    input::{Location, Month},
    p533::{
        consts::{D2R, MAXSSN, R2D},
        control_point::ControlPt,
        ionosphere::{Maps, basis},
        magfit::magfit,
    },
};
use std::f64::consts::PI;

/// Average angle per day
const A: f64 = 0.98565327;
/// Minutes per degree of Earth's rotation
const B: f64 = 3.98891967;
/// Earth's tile sine
const S: f64 = 0.397948631;
/// Earth's tile cosine
const C: f64 = 0.917407699;
/// Value of nu on March 21st
const V: f64 = 78.746118 * D2R;
/// To determine the day count of the day of interest.
const DOTY: [i32; 12] = [0, 31, 59, 90, 120, 152, 181, 212, 243, 273, 304, 334];
const DAY_CONST: i32 = 15;
const SHAPE_FACTOR: f64 = 1.91516873;

/// Evaluate the ionosphere, magnetic field and Sun at a path control point.
pub(super) fn calculate_ctrlpt_param(
    pt: &mut ControlPt,
    month: Month,
    hour: usize,
    ssn: u16,
    maps: &Maps,
) {
    magfit::<100>(pt);
    magfit::<300>(pt);
    (pt.fof2, pt.m3kf2) = maps.at(&basis(pt));
    solar_parameters(pt, month, hour as f64);
    find_foe(pt, month, hour, ssn);
}

/// Equation of time and declination depend on the date/hour, not the location.
pub(super) fn solar_time(month: Month, hour: f64) -> (f64, f64) {
    let d = f64::from(DOTY[month as usize]) + f64::from(DAY_CONST) + hour / 24.0;
    let lambda = A * D2R * (d - 2.0);
    let nu = (SHAPE_FACTOR * D2R).mul_add(lambda.sin(), lambda);
    let mut epsilon = A * D2R * (d - 80.0);

    if epsilon >= 270.0 * D2R {
        epsilon = 2.0f64.mul_add(-PI, epsilon);
    } else if epsilon >= 90.0 * D2R {
        epsilon -= PI;
    }

    let beta = (C * epsilon.tan()).atan();

    // Equation of time: tilt effect + elliptic effect
    let eot = B * ((epsilon - beta) + (lambda - nu)) * R2D;

    // Solar declination, radians
    let decl = (S
        * ((A * (d - 2.0)).mul_add(D2R, (A * (d - 2.0) * D2R).sin() * 0.016713) - V).sin())
    .asin();

    (eot, decl)
}

pub(super) fn solar_zenith(location: Location, hour: f64, (eot, decl): (f64, f64)) -> f64 {
    // Integer time zone from longitude (matches C's (int) cast: trunc toward zero)
    let tzone = (location.lng / (15.0 * D2R)).trunc() as i32;

    // Local time using integer time zone only
    let ltime = hour + f64::from(tzone);

    // Hour angle from solar time corrected for local longitude and EOT
    let toffset = ((location.lng / (15.0 * D2R)) - f64::from(tzone)).mul_add(60.0, eot);

    // True solar time in minutes
    let tst = (ltime * 60.0) + toffset;
    let ha = ((tst / 4.0) - 180.0) * D2R;

    // Cosine of the solar zenith angle (correct product form)
    let mut cosphi =
        (location.lat.cos() * decl.cos()).mul_add(ha.cos(), location.lat.sin() * decl.sin());

    if cosphi.abs() > 1.0 {
        if cosphi >= 0.0 {
            cosphi = 1.0;
        } else {
            cosphi = -1.0;
        }
    }

    // Solar zenith angle which is positive for all latitudes.
    cosphi.acos()
}

pub(super) fn solar_parameters(pt: &mut ControlPt, month: Month, hour: f64) {
    let (eot, decl) = solar_time(month, hour);
    pt.sun.decl = decl;
    pt.sun.sza = solar_zenith(pt.l, hour, (eot, decl));
    // Hour angle at sunrise and sunset in radians
    let sha =
        pt.l.lat
            .tan()
            .mul_add(
                -decl.tan(),
                (90.833 * D2R).cos() / (pt.l.lat.cos() * decl.cos()),
            )
            .acos();

    // Local sunrise/sunset/noon relative to UTC in fractional hours
    pt.sun.lsr = (((-pt.l.lng - sha) * R2D).mul_add(4.0, 720.0) - eot) / 60.0;
    pt.sun.lss = (((-pt.l.lng + sha) * R2D).mul_add(4.0, 720.0) - eot) / 60.0;
    pt.sun.lsn = (((-pt.l.lng) * R2D).mul_add(4.0, 720.0) - eot) / 60.0;

    // Wrap to [0, 24)
    pt.sun.lsr = (pt.sun.lsr + 24.0).rem_euclid(24.0);
    pt.sun.lss = (pt.sun.lss + 24.0).rem_euclid(24.0);
    pt.sun.lsn = (pt.sun.lsn + 24.0).rem_euclid(24.0);

    pt.ltime = (hour + pt.l.lng / (15.0 * D2R)).rem_euclid(24.0);
}

#[must_use]
pub(super) const fn bilinear_interpolation(
    ll: f64,
    lr: f64,
    ul: f64,
    ur: f64,
    row: f64,
    column: f64,
) -> f64 {
    ur.mul_add(
        row * column,
        ll * ((1.0 - row) * (1.0 - column))
            + ul * (row * (1.0 - column))
            + lr * ((1.0 - row) * column),
    )
}

pub(super) fn find_foe(pt: &mut ControlPt, month: Month, mut hour: usize, ssn: u16) {
    // Limit SSN
    let ssn = ssn.min(MAXSSN);

    // A: Solar activity factor
    let phi = 0.00089f64.mul_add(
        f64::from(ssn).powi(2),
        0.728f64.mul_add(f64::from(ssn), 63.7),
    );
    let a = 0.0094f64.mul_add(phi - 66.0, 1.0);

    // B: Seasonal factor
    let lat = pt.l.lat;
    let decl = pt.sun.decl;

    let m = if lat.abs() < 32.0 * D2R {
        1.92f64.mul_add(lat.cos(), -1.93)
    } else {
        0.49f64.mul_add(-lat.cos(), 0.11)
    };

    let n = if (lat - decl).abs() < 80.0 * D2R {
        lat - decl
    } else {
        80.0 * D2R
    };

    let b = n.cos().powf(m);

    // C: Main latitude factor
    let (x, y): (f64, f64) = if lat.abs() < 32.0 * D2R {
        (23.0, 116.0)
    } else {
        (92.0, 35.0)
    };

    let c = y.mul_add(lat.cos(), x);

    // D: Time-of-day factor
    let sza = pt.sun.sza;
    let p = if lat.abs() <= 12.0 * D2R { 1.31 } else { 1.2 };

    let d = if sza <= 73.0 * D2R {
        sza.cos().powf(p)
    } else if sza < std::f64::consts::FRAC_PI_2 {
        let dsza = 6.27e-13 * (sza.mul_add(R2D, -50.0).powi(8)) * D2R;
        (sza - dsza).cos().powf(p)
    } else {
        // Adjust time
        hour = (hour + 1) % 24;

        let lss = pt.sun.lss;
        let lsr = pt.sun.lsr;

        let hour = hour as f64;
        let h = if hour >= lss && (lss >= lsr || hour < lsr) {
            hour - lss
        } else if lss >= lsr && hour < lsr {
            24.0 - lss + hour
        } else {
            0.0
        };

        // Polar winter test
        let polar_lat = 72.5622 * D2R;
        let is_polar_winter = (lat > polar_lat
            && matches!(month, Month::Nov | Month::Dec | Month::Jan))
            || (lat < -polar_lat && matches!(month, Month::May | Month::Jun | Month::Jul));

        if is_polar_winter {
            (0.072f64).powf(p) * (0.28 * sza).mul_add(-R2D, 25.2).exp()
        } else {
            let d1 = (0.072f64).powf(p) * (-1.4 * h).exp();
            let d2 = (0.072f64).powf(p) * (0.28 * sza).mul_add(-R2D, 25.2).exp();
            d1.max(d2)
        }
    };

    // Final foE
    let foe1 = (a * b * c * d).powf(0.25);
    let foe2 = (0.004 * 0.021f64.mul_add(phi, 1.0).powi(2)).powf(0.25);
    pt.foe = foe1.max(foe2);
}
