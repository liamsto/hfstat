use crate::{
    antenna::Antenna,
    input::{DigitalConfig, Location, Modulation, Month, Path, Request},
    p533::{
        consts::{MAXEMDS, MAXF2MDS, R0, TINY_DB},
        control_point::{ControlPointId::*, ControlPt},
        controlpt_params::calculate_ctrlpt_param,
        geometry::{great_circle_distance, great_circle_point},
        ionosphere::Maps,
        mode::Mode,
        noise::params::NoiseParams,
    },
};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum ModulationType {
    #[default]
    Analog,
    Digital,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) enum Season {
    #[default]
    Winter,
    Equinox,
    Summer,
}

#[derive(Default)]
pub(super) struct State<'a> {
    pub(crate) month: Month,
    pub(crate) hour: usize,
    pub(crate) ssn: u16,
    pub(crate) modulation: ModulationType,
    pub(crate) sorl: Path,
    pub(crate) frequency: f64,
    pub(crate) bw: f64,
    pub(crate) txpower: f64,
    pub(crate) snrxxp: i32,
    pub(crate) snrr: f64,
    pub(crate) sirr: f64,
    pub(crate) f0: f64,
    pub(crate) t0: f64,
    pub(crate) a: f64,
    pub(crate) tw: f64,
    pub(crate) fw: f64,
    pub(crate) location: (Location, Location),
    pub(crate) antenna: (Antenna<'a>, Antenna<'a>),
    pub(crate) season: Season,
    pub(crate) distance: f64,
    pub(crate) ptick: f64,
    pub(crate) dmax: f64,
    pub(crate) ele: f64,
    pub(crate) bmuf: f64,
    pub(crate) muf50: f64,
    pub(crate) muf90: f64,
    pub(crate) muf10: f64,
    pub(crate) opmuf: f64,
    pub(crate) opmuf90: f64,
    pub(crate) opmuf10: f64,
    pub(crate) n0_f2: Option<usize>,
    pub(crate) n0_e: Option<usize>,
    pub(crate) es: f64,
    pub(crate) el: f64,
    pub(crate) ei: f64,
    pub(crate) ep: f64,
    pub(crate) pr: f64,
    pub(crate) e0: f64,
    pub(crate) gap: f64,
    pub(crate) ly: f64,
    pub(crate) fm: f64,
    pub(crate) fl: f64,
    pub(crate) fh: f64,
    pub(crate) gtl: f64,
    pub(crate) k: (f64, f64),
    pub(crate) snr: f64,
    pub(crate) dusn: f64,
    pub(crate) dlsn: f64,
    pub(crate) snrxx: f64,
    pub(crate) sir: f64,
    pub(crate) dusi: f64,
    pub(crate) dlsi: f64,
    pub(crate) has_interference: bool,
    pub(crate) rsn: f64,
    pub(crate) rt: f64,
    pub(crate) rf: f64,
    pub(crate) bcr: f64,
    pub(crate) ocr: f64,
    pub(crate) ocrs: f64,
    pub(crate) mir: f64,
    pub(crate) probocc: f64,
    pub(crate) grw: f64,
    pub(crate) controlpt: [ControlPt; 5],
    pub(crate) md_f2: [Mode; MAXF2MDS],
    pub(crate) md_e: [Mode; MAXEMDS],
    pub(crate) noisep: NoiseParams,
}

impl<'a> State<'a> {
    pub(crate) fn new(request: Request<'a>, maps: &Maps) -> Self {
        let (modulation, digital) = match request.modulation {
            Modulation::Analog => (
                ModulationType::Analog,
                DigitalConfig {
                    sir_db: 0.0,
                    f0_hz: 0.0,
                    t0_ms: 0.0,
                    a_db: 0.0,
                    tw_ms: 0.0,
                    fw_hz: 0.0,
                },
            ),
            Modulation::Digital(config) => (ModulationType::Digital, config),
        };

        let mut path = Self {
            month: request.month,
            hour: usize::from(request.utc_hour),
            ssn: request.ssn,
            modulation,
            sorl: request.path,
            frequency: request.freq_mhz,
            bw: request.bandwidth_hz,
            txpower: request.power_dbkw,
            snrxxp: i32::from(request.snr_percent),
            snrr: request.snr_db,
            sirr: digital.sir_db,
            f0: digital.f0_hz,
            t0: digital.t0_ms,
            a: digital.a_db,
            tw: digital.tw_ms,
            fw: digital.fw_hz,
            location: (request.tx.radians(), request.rx.radians()),
            antenna: (request.tx_antenna, request.rx_antenna),
            ei: TINY_DB,
            el: TINY_DB,
            noisep: NoiseParams {
                man_made_noise: request.noise,
                month: request.month,
                ..NoiseParams::default()
            },
            ..Self::default()
        };

        path.distance = great_circle_distance(path.location.0, path.location.1);
        if path.distance == 0.0 {
            path.distance = f64::EPSILON;
        }
        if path.sorl == Path::Long {
            path.distance = 2.0 * R0 * PI - path.distance;
        }

        path.initialize_ctrlpt(maps);
        path.season = season(path.controlpt[MP].l, path.month);
        path
    }

    fn initialize_ctrlpt(&mut self, maps: &Maps) {
        great_circle_point(
            self.location.0,
            self.location.1,
            &mut self.controlpt[MP],
            self.distance,
            0.5,
        );
        calculate_ctrlpt_param(
            &mut self.controlpt[MP],
            self.month,
            self.hour,
            self.ssn,
            maps,
        );

        if self.distance >= 2000.0 {
            great_circle_point(
                self.location.0,
                self.location.1,
                &mut self.controlpt[R1k],
                self.distance,
                (self.distance - 1000.0) / self.distance,
            );
            great_circle_point(
                self.location.0,
                self.location.1,
                &mut self.controlpt[T1k],
                self.distance,
                1000.0 / self.distance,
            );
            calculate_ctrlpt_param(
                &mut self.controlpt[T1k],
                self.month,
                self.hour,
                self.ssn,
                maps,
            );
            calculate_ctrlpt_param(
                &mut self.controlpt[R1k],
                self.month,
                self.hour,
                self.ssn,
                maps,
            );
        }
    }
}

const fn season(location: Location, month: Month) -> Season {
    if location.lat >= 0.0 {
        match month {
            Month::Nov | Month::Dec | Month::Jan | Month::Feb => Season::Winter,
            Month::Mar | Month::Apr | Month::Sep | Month::Oct => Season::Equinox,
            Month::May | Month::Jun | Month::Jul | Month::Aug => Season::Summer,
        }
    } else {
        match month {
            Month::Nov | Month::Dec | Month::Jan | Month::Feb => Season::Summer,
            Month::Mar | Month::Apr | Month::Sep | Month::Oct => Season::Equinox,
            Month::May | Month::Jun | Month::Jul | Month::Aug => Season::Winter,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seasons_follow_hemisphere() {
        let north = Location::from_radians(0.5, 0.0);
        let south = Location::from_radians(-0.5, 0.0);

        assert!(matches!(season(north, Month::Jan), Season::Winter));
        assert!(matches!(season(south, Month::Jan), Season::Summer));
        assert!(matches!(season(north, Month::Apr), Season::Equinox));
    }
}
