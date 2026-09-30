use crate::{
    input::{Location, Path},
    p533::consts::{D2R, R0},
    p533::control_point::ControlPt,
};
use std::f64::consts::PI;

const GEOMAGNETIC_LATITUDE: f64 = 78.5 * D2R;
const GEOMAGNETIC_LONGITUDE: f64 = -68.2 * D2R;

/// Determines the distance in km between `source` and `destination` on the great circle.
#[must_use]
pub(super) fn great_circle_distance(source: Location, destination: Location) -> f64 {
    2.0 * R0
        * (((source.lat - destination.lat) / 2.0).sin())
            .mul_add(
                ((source.lat - destination.lat) / 2.0).sin(),
                source.lat.cos()
                    * destination.lat.cos()
                    * (((source.lng - destination.lng) / 2.0).sin()).powi(2),
            )
            .sqrt()
            .asin()
}

/// Determines the lat and long of a point on a great circle path of length distance from `source` to `destination`. The point is determined at `fraction` * `distance` from `source`.
pub(super) fn great_circle_point(
    source: Location,
    destination: Location,
    midpoint: &mut ControlPt,
    distance: f64,
    fraction: f64,
) {
    if distance == 0.0 {
        midpoint.l.lat = source.lat;
        midpoint.l.lng = source.lng;
    } else {
        let d = distance / R0;
        let a = ((1.0 - fraction) * d).sin() / d.sin();
        let b = (fraction * d).sin() / d.sin();

        let x = (b * destination.lat.cos()).mul_add(
            destination.lng.cos(),
            a * source.lat.cos() * source.lng.cos(),
        );
        let y = (b * destination.lat.cos()).mul_add(
            destination.lng.sin(),
            a * source.lat.cos() * source.lng.sin(),
        );
        let z = b.mul_add(destination.lat.sin(), a * source.lat.sin());
        midpoint.l.lat = z.atan2(x.hypot(y));
        midpoint.l.lng = y.atan2(x);
        debug_assert!(!midpoint.l.lat.is_nan());
        debug_assert!(!midpoint.l.lat.is_nan());
    }
}

/// Converts from geographic coordinates to geomagnetic coordinates
pub(super) fn geomagnetic_coordinates(source: &Location, destination: &mut Location) {
    destination.lat = (source.lat.cos() * GEOMAGNETIC_LATITUDE.cos())
        .mul_add(
            (source.lng - GEOMAGNETIC_LONGITUDE).cos(),
            source.lat.sin() * GEOMAGNETIC_LATITUDE.sin(),
        )
        .asin();
    destination.lng = (source.lat.cos() * (source.lng - GEOMAGNETIC_LONGITUDE).sin()
        / destination.lat.cos())
    .asin();
}

#[must_use]
pub(super) fn bearing(source: &Location, destination: &Location, direction: Path) -> f64 {
    let numerator = (destination.lng - source.lng).sin() * destination.lat.cos();
    let denominator = (source.lat.sin() * destination.lat.cos()).mul_add(
        -(destination.lng - source.lng).cos(),
        source.lat.cos() * destination.lat.sin(),
    );
    let mut bearing = numerator.atan2(denominator);
    bearing = 2.0f64.mul_add(PI, bearing).rem_euclid(2.0 * PI);
    if direction == Path::Long {
        bearing = 2.0f64.mul_add(PI, bearing + PI).rem_euclid(2.0 * PI);
    }

    bearing
}
