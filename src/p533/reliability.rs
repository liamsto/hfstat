use crate::{
    input::Location,
    p533::{
        consts::{
            D2R, FS2_INV, MAXEMDS, MAXF2MDS, MAXMDS, NORM, R0, R2D, TABLE2LD, TABLE2UD, TINY_DB,
            TS2_INV, VOFL,
        },
        control_point::{ControlPointId, ControlPt},
        geometry::geomagnetic_coordinates,
        mode::Mode,
        screening_freq::elevation_angle,
        state::{ModulationType, State},
    },
};

use crate::p533::control_point::ControlPointId::{MP, R1k, T1k};

#[inline]
fn mode_exists(mode: &Mode) -> bool {
    mode.bmuf.is_some()
}

pub(super) fn circuit_reliability(path: &mut State) {
    let mut i_s = [None; MAXMDS];
    let mut i_i = [None; MAXMDS];
    let mut geomag = Location::default();

    let s = if path.modulation == ModulationType::Analog {
        path.pr
    } else {
        digital_modulation_signal_and_interferers(path, &mut i_s, &mut i_i)
    };

    // SNR calculation
    path.snr = 10.0f64.mul_add(
        -path.bw.log10(),
        10.0f64.mul_add(
            -((10.0_f64.powf(path.noisep.fa_a / 10.0)
                + 10.0_f64.powf(path.noisep.fa_m / 10.0)
                + 10.0_f64.powf(path.noisep.fa_g / 10.0))
            .log10()),
            s,
        ),
    ) + 204.0;

    // Determine if path crosses 60 degrees
    let mut gt60deg = false;
    if path.distance > 2000.0 {
        geomagnetic_coordinates(&path.controlpt[MP].l, &mut geomag);
        if geomag.lat >= 60.0 * D2R {
            gt60deg = true;
        }
        geomagnetic_coordinates(&path.controlpt[T1k].l, &mut geomag);
        if geomag.lat >= 60.0 * D2R {
            gt60deg = true;
        }
        geomagnetic_coordinates(&path.controlpt[R1k].l, &mut geomag);
        if geomag.lat >= 60.0 * D2R {
            gt60deg = true;
        }
    }

    let fbmuf_r = path.frequency / path.bmuf;

    let bmuf = match fbmuf_r {
        x if x <= 0.8 => 0,
        x if x <= 1.0 => 1,
        x if x <= 1.2 => 2,
        x if x <= 1.4 => 3,
        x if x <= 1.6 => 4,
        x if x <= 1.8 => 5,
        x if x <= 2.0 => 6,
        x if x <= 3.0 => 7,
        x if x <= 4.0 => 8,
        _ => 9,
    };

    let dlsd = TABLE2LD[usize::from(gt60deg)][bmuf];
    let dusd = TABLE2UD[usize::from(gt60deg)][bmuf];

    let dush: f64 = 5.0;
    let dlsh: f64 = 8.0;

    let x = 10.0_f64.powf(path.noisep.fa_a / 10.0)
        + 10.0_f64.powf(path.noisep.fa_m / 10.0)
        + 10.0_f64.powf(path.noisep.fa_g / 10.0);

    let y = 10.0_f64.powf((path.noisep.fa_a - path.noisep.dl_a) / 10.0)
        + 10.0_f64.powf((path.noisep.fa_m - path.noisep.dl_m) / 10.0)
        + 10.0_f64.powf((path.noisep.fa_g - path.noisep.dl_g) / 10.0);

    path.dusn = dush
        .mul_add(dush, dusd.mul_add(dusd, (10.0 * (x / y).log10()).powi(2)))
        .sqrt();

    let y = 10.0_f64.powf((path.noisep.fa_a + path.noisep.du_a) / 10.0)
        + 10.0_f64.powf((path.noisep.fa_m + path.noisep.du_m) / 10.0)
        + 10.0_f64.powf((path.noisep.fa_g + path.noisep.du_g) / 10.0);

    path.dlsn = dlsh
        .mul_add(dlsh, dlsd.mul_add(dlsd, (10.0 * (y / x).log10()).powi(2)))
        .sqrt();

    if path.snr >= path.snrr {
        path.bcr = (130.0 - 80.0 / (1.0 + ((path.snr - path.snrr) / path.dlsn))).min(100.0);
    } else {
        path.bcr = (80.0 / (1.0 + ((path.snrr - path.snr) / path.dusn)) - 30.0).max(0.0);
    }

    if path.modulation == ModulationType::Digital && path.t0 != 0.0 && path.f0 != 0.0 {
        let d = path.distance;
        let tm = if d <= 2000.0 {
            (2.5e7
                * (path.frequency / path.bmuf).mul_add(-(path.frequency / path.bmuf), 1.0)
                * d.powi(-2))
            .min(0.00175f64.mul_add(-d, 7.0))
        } else {
            (4.27e-2
                * (path.frequency / path.bmuf).mul_add(-(path.frequency / path.bmuf), 1.0)
                * d.powf(0.65))
            .min(3.5)
        };

        let fm = 0.02 * path.frequency * tm;
        let dtu = 0.15 * tm;
        let dtl = 0.15 * tm;
        let dfu = 0.1 * fm;
        let dfl = 0.1 * fm;

        if path.snr >= path.snrr {
            path.rsn = (130.0 - 80.0 / (1.0 + ((path.snr - path.snrr) / path.dlsn))).min(100.0);
        } else {
            path.rsn = (80.0 / (1.0 + ((path.snrr - path.snr) / path.dusn)) - 30.0).max(0.0);
        }

        if tm >= path.t0 {
            path.rt = (130.0 - 80.0 / (1.0 + (path.t0 - tm) / dtl)).min(100.0);
        } else {
            path.rt = (80.0 / (1.0 + (tm - path.t0) / dtu) - 30.0).max(0.0);
        }

        if fm >= path.f0 {
            path.rf = (130.0 - 80.0 / (1.0 + (path.f0 - fm) / dfl)).min(100.0);
        } else {
            path.rf = (80.0 / (1.0 + (fm - path.f0) / dfu) - 30.0).max(0.0);
        }
    }

    if path.modulation == ModulationType::Digital && path.distance <= 9000.0 {
        let dl_ih = 8.0;
        let du_ih = 5.0;
        let mut i_sum = 0.0;
        let mut i_sumu = 0.0;
        let mut i_suml = 0.0;

        for idx in i_i.iter().take(MAXMDS).flatten().copied() {
            if idx < MAXEMDS {
                // E-mode
                if let Some(e) = path.md_e.get(idx) {
                    i_sum += 10.0_f64.powf((e.prw - path.a) / 10.0);
                    i_sumu += 10.0_f64.powf((e.prw - path.a + du_ih) / 10.0);
                    i_suml += 10.0_f64.powf((e.prw - path.a - dl_ih) / 10.0);
                }
            } else {
                // F2-mode
                let f2_idx = idx - MAXEMDS;
                if let Some(f2) = path.md_f2.get(f2_idx) {
                    i_sum += 10.0_f64.powf((f2.prw - path.a) / 10.0);
                    i_sumu += 10.0_f64.powf((f2.prw - path.a + du_ih) / 10.0);
                    i_suml += 10.0_f64.powf((f2.prw - path.a - dl_ih) / 10.0);
                }
            }
        }
        // With no interfering modes the multimode reliability is 100%.
        if i_sum == 0.0 {
            path.dusi = dush;
            path.dlsi = dlsh;
            path.mir = 100.0;
        } else {
            path.has_interference = true;
            path.sir = 10.0f64.mul_add(-i_sum.log10(), s);
            path.dusi = dush.hypot(10.0 * (i_sum / i_suml).log10());
            path.dlsi = dlsh.hypot(10.0 * (i_sumu / i_sum).log10());

            if path.sir >= path.sirr {
                path.mir = (130.0 - 80.0 / (1.0 + ((path.sir - path.sirr) / path.dlsi))).min(100.0);
            } else {
                path.mir = ((80.0 / (1.0 + ((path.sirr - path.sir) / path.dusi))) - 30.0).max(0.0);
            }
        }

        path.ocr = path.bcr * path.mir / 100.0;
        equatorial_scattering(path, i_s);
    }

    if path.snrxxp < 50 {
        path.snrxx = path.snr + path.dusn * NORM[50 - path.snrxxp as usize] / NORM[40];
    } else {
        path.snrxx = path.snr - path.dlsn * NORM[path.snrxxp as usize - 50] / NORM[40];
    }
}

#[derive(Clone, Copy)]
enum Criteria {
    Field,
    Delay,
}

fn digital_modulation_signal_and_interferers(
    path: &mut State,
    i_s: &mut [Option<usize>; MAXMDS],
    i_i: &mut [Option<usize>; MAXMDS],
) -> f64 {
    if path.distance <= 9000.0 {
        // --- Calculate delays for E layer
        let hr_default = 110.0;
        for (n, mode) in path.md_e.iter_mut().enumerate() {
            if mode_exists(mode) {
                let dh = path.distance / (n as f64 + 1.0);
                let delta = elevation_angle(dh, hr_default);
                let psi = dh / (2.0 * R0);
                let ptick = 2.0 * R0 * (psi.sin() / (delta - psi).cos());
                mode.tau = (n as f64 + 1.0) * (ptick / VOFL) * 1000.0;
            }
        }

        // --- Calculate delays for F2 layer
        for (n, mode) in path.md_f2.iter_mut().enumerate() {
            if mode_exists(mode) {
                let hr = mode.hr;
                let dh = path.distance / (n as f64 + 1.0);
                let delta = elevation_angle(dh, hr);
                let psi = dh / (2.0 * R0);
                let ptick = 2.0 * R0 * (psi.sin() / (delta - psi).cos());
                mode.tau = (n as f64 + 1.0) * (ptick / VOFL) * 1000.0;
            }
        }

        if number_of_modes(path) >= 2 {
            // --- Flatten all modes into a single array of references
            let mut m: [Option<&Mode>; MAXMDS] = [None; MAXMDS];
            for (n, m_n) in m.iter_mut().enumerate().take(MAXEMDS) {
                *m_n = Some(&path.md_e[n]);
            }
            for n in 0..MAXF2MDS {
                m[n + 3] = Some(&path.md_f2[n]);
            }

            // --- Mode sorting (indices into m[])
            let i_ew = mode_sort(&m, Criteria::Field);
            let i_tau = mode_sort(&m, Criteria::Delay);

            let Some(tau_idx) = i_tau.first() else {
                return path.pr;
            };
            let Some(ew_idx) = i_ew.first() else {
                return path.pr;
            };

            let Some(tau_mode) = m.get(tau_idx).and_then(|mode| *mode) else {
                return path.pr;
            };
            let Some(ew_mode) = m.get(ew_idx).and_then(|mode| *mode) else {
                return path.pr;
            };

            let deltat = tau_mode.tau + path.tw;
            let delta_a = ew_mode.ew - path.a;

            let mut power_sum = 0.0;

            for (n, mode_opt) in m.iter().enumerate() {
                if let Some(mode) = mode_opt
                    && mode_exists(mode)
                {
                    if mode.ew >= delta_a && mode.tau <= deltat {
                        power_sum += 10.0f64.powf(mode.prw / 10.0);
                        i_s[n] = Some(n);
                    } else {
                        i_i[n] = Some(n);
                    }
                }
            }

            let signal = if power_sum > 0.0 {
                10.0 * power_sum.log10()
            } else {
                TINY_DB
            };

            sort_mode_indices(i_s);
            sort_mode_indices(i_i);
            signal
        } else {
            for (index, mode) in path.md_e.iter().enumerate() {
                if mode_exists(mode) {
                    i_s[0] = Some(index);
                    break;
                }
            }
            if i_s[0].is_none() {
                for (index, mode) in path.md_f2.iter().enumerate() {
                    if mode_exists(mode) {
                        i_s[0] = Some(index + MAXEMDS);
                        break;
                    }
                }
            }
            path.pr
        }
    } else {
        // Distance > 9000
        path.pr
    }
}

/// Sorts the indices of modes based on the provided criteria.
/// Filters out modes with BMUF == 0.0 (nonexistent modes).
#[derive(Default)]
struct ModeOrder {
    indices: [usize; MAXMDS],
    len: usize,
}

impl ModeOrder {
    fn first(&self) -> Option<usize> {
        (self.len != 0).then_some(self.indices[0])
    }
}

fn mode_sort(modes: &[Option<&Mode>; MAXMDS], criteria: Criteria) -> ModeOrder {
    let mut order = ModeOrder::default();
    for (index, mode) in modes.iter().enumerate() {
        if mode.is_some_and(mode_exists) {
            order.indices[order.len] = index;
            order.len += 1;
        }
    }

    for index in 1..order.len {
        let value = order.indices[index];
        let mut insert = index;
        while insert > 0 && mode_before(modes, value, order.indices[insert - 1], criteria) {
            order.indices[insert] = order.indices[insert - 1];
            insert -= 1;
        }
        order.indices[insert] = value;
    }
    order
}

fn mode_before(
    modes: &[Option<&Mode>; MAXMDS],
    left: usize,
    right: usize,
    criteria: Criteria,
) -> bool {
    let left = modes[left].expect("mode order contains only existing modes");
    let right = modes[right].expect("mode order contains only existing modes");
    match criteria {
        Criteria::Field => left.ew > right.ew,
        Criteria::Delay => left.tau < right.tau,
    }
}

fn number_of_modes(path: &State) -> usize {
    path.md_e.iter().filter(|m| mode_exists(m)).count()
        + path.md_f2.iter().filter(|m| mode_exists(m)).count()
}

fn sort_mode_indices(indices: &mut [Option<usize>; MAXMDS]) {
    let mut write = 0;
    for read in 0..MAXMDS {
        if let Some(index) = indices[read] {
            indices[write] = Some(index);
            write += 1;
        }
    }
    indices[write..].fill(None);
}

/// Applies the equatorial-scattering model in P.533-14 § 10.3.
fn equatorial_scattering(path: &mut State, is_indices: [Option<usize>; MAXMDS]) {
    if path.distance > 9000.0 || is_indices[0].is_none() {
        return;
    }

    let mut earliest = f64::INFINITY;
    let mut dominant_f2: Option<usize> = None;
    for index in is_indices.iter().flatten().copied() {
        let mode = if index < MAXEMDS {
            &path.md_e[index]
        } else {
            &path.md_f2[index - MAXEMDS]
        };
        earliest = earliest.min(mode.tau);
        if index < MAXEMDS {
            continue;
        }
        let f2 = index - MAXEMDS;
        if dominant_f2.is_none_or(|best| path.md_f2[best].ew < mode.ew) {
            dominant_f2 = Some(f2);
        }
    }

    let Some(dominant_f2) = dominant_f2 else {
        path.ocrs = path.ocr;
        return;
    };
    let dominant_level = path.md_f2[dominant_f2].ew;
    let threshold = dominant_level - path.a;
    let edge = earliest + path.tw;
    let mut affected = [false; MAXF2MDS];

    for index in is_indices.iter().flatten().copied() {
        if index >= MAXEMDS {
            let f2 = index - MAXEMDS;
            let mode = path.md_f2[f2];
            let spread = scattered_level(mode.ew, edge - mode.tau, TS2_INV);
            affected[f2] = spread >= threshold;
        }
    }
    if scattered_level(dominant_level, path.fw, FS2_INV) >= threshold {
        affected[dominant_f2] = true;
    }
    if !affected.iter().any(|affected| *affected) {
        path.ocrs = path.ocr;
        return;
    }

    let Some(n0_f2) = path.n0_f2 else {
        path.ocrs = path.ocr;
        return;
    };
    let sunspot = sunspot_factor(path.ssn);
    let seasonal = 0.45f64.mul_add(
        (60.0 * D2R * (path.month.index() as f64 + 1.0 - 1.5)).sin(),
        0.55,
    );

    let mut probability: f64 = 0.0;
    for (mode, affected) in affected.into_iter().enumerate() {
        if affected {
            probability = probability.max(mode_probability(path, mode, n0_f2, sunspot, seasonal));
        }
    }

    path.probocc = probability;
    path.ocrs = path.ocr * (1.0 - probability);
}

fn scattered_level(specular_db: f64, delta: f64, inverse_variance: f64) -> f64 {
    specular_db + 10.0 * (0.056 * (-delta * delta * inverse_variance).exp()).log10()
}

fn mode_probability(
    path: &State<'_>,
    mode: usize,
    lowest: usize,
    sunspot: f64,
    seasonal: f64,
) -> f64 {
    let mut probability: f64 = 0.0;
    if mode == lowest {
        if path.distance <= path.dmax {
            probability = occurrence(&path.controlpt[ControlPointId::MP], sunspot, seasonal);
        } else {
            for point in [ControlPointId::Td02, ControlPointId::Rd02] {
                probability =
                    probability.max(occurrence(&path.controlpt[point], sunspot, seasonal));
            }
        }
    } else if path.distance <= path.dmax {
        for point in [ControlPointId::T1k, ControlPointId::R1k, ControlPointId::MP] {
            probability = probability.max(occurrence(&path.controlpt[point], sunspot, seasonal));
        }
    } else {
        for point in [
            ControlPointId::T1k,
            ControlPointId::R1k,
            ControlPointId::MP,
            ControlPointId::Td02,
            ControlPointId::Rd02,
        ] {
            probability = probability.max(occurrence(&path.controlpt[point], sunspot, seasonal));
        }
    }
    probability
}

fn occurrence(point: &ControlPt, sunspot: f64, seasonal: f64) -> f64 {
    find_f_lambda(point) * find_ftl(point) * sunspot * seasonal
}

///Determines the F sub T sub l parameter in in P.533-14 Appendix 1 to Annex 1 "A model for scattering of HF signals".
const fn find_ftl(cpt: &ControlPt) -> f64 {
    let tl = cpt.ltime;

    if tl >= 0.0 && tl <= 3.0 {
        1.0
    } else if (tl > 3.0) && (tl <= 7.0) {
        let t1 = (7.0 - tl) / 4.0;
        t1 * t1 * ((tl - 1.0) / 2.0)
    } else if (tl > 7.0) && (tl <= 19.0) {
        0.0
    } else if (tl > 19.0) && (tl <= 20.0) {
        let t2 = tl - 19.0;
        t2 * t2 * 2.0f64.mul_add(-tl, 41.0)
    } else if (tl > 20.0) && (tl <= 24.0) {
        1.0
    } else {
        0.0
    }
}

/// Determines F sub lambda sub d in P.533-14 Appendix 1 to Annex 1 "A model for scattering of HF signals"
fn find_f_lambda(cpt: &ControlPt) -> f64 {
    let dip_deg = cpt.dip.0.abs() * R2D;
    if dip_deg <= 15.0 {
        1.0
    } else if dip_deg < 25.0 {
        ((25.0 - dip_deg) / 10.0).powi(2) * ((dip_deg - 10.0) / 5.0)
    } else {
        0.0
    }
}

fn sunspot_factor(ssn: u16) -> f64 {
    0.008f64.mul_add(f64::from(ssn), 0.1).min(1.0)
}

#[cfg(test)]
mod tests {
    use crate::p533::{
        consts::D2R,
        control_point::ControlPt,
        reliability::{find_f_lambda, find_ftl, scattered_level, sunspot_factor},
    };

    #[test]
    fn scattering_latitude_uses_degrees() {
        let mut point = ControlPt::default();

        point.dip.0 = 15.0 * D2R;
        assert_eq!(find_f_lambda(&point), 1.0);
        point.dip.0 = 20.0 * D2R;
        assert!((find_f_lambda(&point) - 0.5).abs() < 1e-12);
        point.dip.0 = 25.0 * D2R;
        assert_eq!(find_f_lambda(&point), 0.0);
    }

    #[test]
    fn scattering_sunspot_factor_is_capped() {
        assert!((sunspot_factor(100) - 0.9).abs() < 1e-12);
        assert_eq!(sunspot_factor(160), 1.0);
        assert_eq!(sunspot_factor(311), 1.0);
    }

    #[test]
    fn scattering_windows_use_relative_offsets() {
        let at_center = scattered_level(50.0, 0.0, 1.0 / 18.0);
        let at_three_hz = scattered_level(50.0, 3.0, 1.0 / 18.0);

        assert!((at_center - (50.0 + 10.0 * 0.056f64.log10())).abs() < 1e-12);
        assert!(at_three_hz < at_center);
    }

    #[test]
    fn midnight_has_full_local_time_factor() {
        let point = ControlPt::default();
        assert_eq!(find_ftl(&point), 1.0);
    }
}
