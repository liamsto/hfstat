use crate::{
    antenna::{AZIMUTHS, Antenna, ELEVATIONS, GAINS_PER_FREQUENCY, Model, Orientation},
    p533::{
        consts::{D2R, R2D, TINY_DB},
        controlpt_params::bilinear_interpolation,
        geometry::bearing,
        state::State,
    },
};

#[derive(Clone, Copy)]
pub(super) enum Direction {
    TxToRx,
    RxToTx,
}

const fn pattern_idx(frequency: usize, azimuth: usize, elevation: usize) -> usize {
    (frequency * AZIMUTHS + azimuth) * ELEVATIONS + elevation
}

fn closest_frequency(frequencies: &[f64], target: f64) -> usize {
    let mut best = 0;
    let mut delta = f64::INFINITY;
    for (index, &frequency) in frequencies.iter().enumerate() {
        let next = (frequency - target).abs();
        if next < delta {
            best = index;
            delta = next;
        }
    }
    best
}

fn path_bearing_deg(path: &State<'_>, direction: Direction) -> f64 {
    let bearing = match direction {
        Direction::TxToRx => bearing(&path.location.0, &path.location.1, path.sorl),
        Direction::RxToTx => bearing(&path.location.1, &path.location.0, path.sorl),
    };
    (bearing * R2D).rem_euclid(360.0)
}

fn gain_from_pattern(
    path: &State<'_>,
    antenna: Antenna<'_>,
    frequencies: &[f64],
    gains: &[f64],
    elevation_rad: f64,
    direction: Direction,
) -> f64 {
    debug_assert_eq!(gains.len(), frequencies.len() * GAINS_PER_FREQUENCY);
    let frequency = closest_frequency(frequencies, path.frequency);
    let azimuth = match antenna.orientation() {
        Orientation::Path => 0.0,
        Orientation::BearingDeg(value) => {
            (path_bearing_deg(path, direction) - value).rem_euclid(360.0)
        }
    };

    let elevation = (elevation_rad * R2D).clamp(0.0, 90.0);
    let elevation_low = elevation.floor() as usize;
    let elevation_high = elevation.ceil() as usize;
    let azimuth_low = azimuth.floor() as usize;
    let azimuth_high = azimuth.ceil() as usize % AZIMUTHS;

    let low_low = gains[pattern_idx(frequency, azimuth_low, elevation_low)];
    let low_high = gains[pattern_idx(frequency, azimuth_high, elevation_low)];
    let high_low = gains[pattern_idx(frequency, azimuth_low, elevation_high)];
    let high_high = gains[pattern_idx(frequency, azimuth_high, elevation_high)];

    bilinear_interpolation(
        low_low,
        low_high,
        high_low,
        high_high,
        elevation.fract(),
        azimuth.fract(),
    )
}

pub(super) fn antenna_gain(
    path: &State<'_>,
    antenna: Antenna<'_>,
    elevation_rad: f64,
    direction: Direction,
) -> f64 {
    match antenna.model() {
        Model::Isotropic(gain_dbi) => gain_dbi,
        Model::Pattern(pattern) => gain_from_pattern(
            path,
            antenna,
            pattern.frequencies_mhz(),
            pattern.gains_dbi(),
            elevation_rad,
            direction,
        ),
    }
}

pub(super) fn antenna_gain_08(
    path: &State<'_>,
    antenna: Antenna<'_>,
    direction: Direction,
    elevation: &mut f64,
) -> f64 {
    if let Model::Isotropic(gain) = antenna.model() {
        if gain > TINY_DB {
            *elevation = 0.0;
            return gain;
        }
        return TINY_DB;
    }
    let mut max_gain = TINY_DB;
    for degree in 0..=8 {
        let angle = f64::from(degree) * D2R;
        let gain = antenna_gain(path, antenna, angle, direction);
        if gain > max_gain {
            max_gain = gain;
            *elevation = angle;
        }
    }
    max_gain
}

#[cfg(test)]
mod tests {
    use crate::{
        antenna::{AZIMUTHS, Antenna, ELEVATIONS, GAINS_PER_FREQUENCY, Orientation, Pattern},
        input::{Location, Path},
        p533::{
            antenna::{Direction, antenna_gain, antenna_gain_08},
            consts::{D2R, TINY_DB},
            state::State,
        },
    };

    fn index(frequency: usize, azimuth: usize, elevation: usize) -> usize {
        (frequency * AZIMUTHS + azimuth) * ELEVATIONS + elevation
    }

    #[test]
    fn isotropic_gain_is_constant() {
        let path = State::default();
        let gain = antenna_gain(
            &path,
            Antenna::isotropic(3.5),
            12.0 * D2R,
            Direction::TxToRx,
        );
        assert_eq!(gain, 3.5);
    }

    #[test]
    fn isotropic_max_preserves_the_floor_and_first_elevation() {
        for (gain, expected, angle) in [
            (-400.0, TINY_DB, 1.0),
            (TINY_DB, TINY_DB, 1.0),
            (3.5, 3.5, 0.0),
        ] {
            let mut elevation = 1.0;
            assert_eq!(
                antenna_gain_08(
                    &State::default(),
                    Antenna::isotropic(gain),
                    Direction::TxToRx,
                    &mut elevation,
                ),
                expected,
            );
            assert_eq!(elevation, angle);
        }
    }

    #[test]
    fn fixed_orientation_rotates_pattern() {
        let mut gains = vec![0.0; GAINS_PER_FREQUENCY];
        gains[index(0, 0, 0)] = 42.0;
        let pattern = Pattern::new(&[10.0], &gains).unwrap();
        let antenna = Antenna::pattern(pattern, Orientation::BearingDeg(90.0));
        let path = State {
            sorl: Path::Short,
            location: (
                Location::from_radians(0.0, 0.0),
                Location::from_radians(0.0, 90.0 * D2R),
            ),
            frequency: 10.0,
            ..State::default()
        };

        assert_eq!(antenna_gain(&path, antenna, 0.0, Direction::TxToRx), 42.0);
    }

    #[test]
    fn path_orientation_stays_on_pattern_boresight() {
        let mut gains = vec![0.0; GAINS_PER_FREQUENCY];
        gains[index(0, 0, 0)] = 42.0;
        let pattern = Pattern::new(&[10.0], &gains).unwrap();
        let path = State {
            sorl: Path::Short,
            location: (
                Location::from_radians(0.0, 0.0),
                Location::from_radians(1.0 * D2R, 90.5 * D2R),
            ),
            frequency: 10.0,
            ..State::default()
        };

        assert_eq!(
            antenna_gain(
                &path,
                Antenna::pattern(pattern, Orientation::Path),
                0.0,
                Direction::TxToRx,
            ),
            42.0
        );
    }

    #[test]
    fn fractional_bearing_is_interpolated() {
        let mut gains = vec![0.0; GAINS_PER_FREQUENCY];
        gains[index(0, 1, 0)] = 4.0;
        let pattern = Pattern::new(&[10.0], &gains).unwrap();
        let antenna = Antenna::pattern(pattern, Orientation::BearingDeg(89.75));
        let path = State {
            sorl: Path::Short,
            location: (
                Location::from_radians(0.0, 0.0),
                Location::from_radians(0.0, 90.5 * D2R),
            ),
            frequency: 10.0,
            ..State::default()
        };

        assert!((antenna_gain(&path, antenna, 0.0, Direction::TxToRx) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn nearest_frequency_is_selected() {
        let mut gains = vec![0.0; GAINS_PER_FREQUENCY * 2];
        gains[..GAINS_PER_FREQUENCY].fill(1.0);
        gains[GAINS_PER_FREQUENCY..].fill(5.0);
        let pattern = Pattern::new(&[10.0, 20.0], &gains).unwrap();
        let path = State {
            frequency: 18.0,
            ..State::default()
        };

        assert_eq!(
            antenna_gain(
                &path,
                Antenna::pattern(pattern, Orientation::Path),
                0.0,
                Direction::TxToRx,
            ),
            5.0
        );
    }
}
