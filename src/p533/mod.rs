mod antenna;
mod consts;
mod control_point;
mod controlpt_params;
mod geometry;
mod ionosphere;
mod magfit;
mod midpath;
mod mode;
mod muf_basic;
mod muf_operational;
mod muf_variability;
mod noise;
mod receiver_power;
mod reliability;
mod screening_freq;
mod skywave_long;
mod skywave_short;
mod solar;
mod state;

use crate::{
    error::{Error, Result},
    input::Request,
    output::{
        DigitalResult, Mode as ModeResult, Muf, Noise, NoiseComponent, Prediction, Reliability, Snr,
    },
    p533::{
        ionosphere::Maps,
        midpath::midpath,
        mode::Mode,
        muf_basic::muf_basic,
        muf_operational::muf_operational,
        muf_variability::muf_variability,
        noise::model::noise,
        receiver_power::median_available_receiver_power,
        reliability::circuit_reliability,
        screening_freq::e_layer_screening_frequency,
        skywave_long::median_skywave_field_strength_long,
        skywave_short::median_skywave_field_strength_short,
        state::{ModulationType, State},
    },
};

/// Predict an HF circuit using the embedded ITU reference data, without I/O.
pub fn predict(request: Request<'_>) -> Result<Prediction> {
    request.validate()?;
    let maps = Maps::new(request.month, usize::from(request.utc_hour), request.ssn);
    let mut path = State::new(request, &maps);
    muf_basic(&mut path, &maps).map_err(Error::Calculation)?;
    muf_variability(&mut path);
    muf_operational(&mut path);
    e_layer_screening_frequency(&mut path);
    median_skywave_field_strength_short(&mut path);
    median_skywave_field_strength_long(&mut path, &maps);
    midpath(&mut path);
    median_available_receiver_power(&mut path);
    noise(
        &mut path.noisep,
        path.hour,
        path.location.1.lng,
        path.location.1.lat,
        path.frequency,
    );
    circuit_reliability(&mut path);
    if !finite(&path) {
        return Err(Error::Calculation(
            "model data produced a non-finite prediction",
        ));
    }
    Ok(output(&path))
}

fn finite(path: &State<'_>) -> bool {
    let main = [
        path.distance,
        path.dmax,
        path.ele,
        path.bmuf,
        path.muf50,
        path.muf90,
        path.muf10,
        path.opmuf,
        path.opmuf90,
        path.opmuf10,
        path.es,
        path.ei,
        path.el,
        path.ep,
        path.pr,
        path.gtl,
        path.grw,
        path.noisep.fa_a,
        path.noisep.du_a,
        path.noisep.dl_a,
        path.noisep.fa_m,
        path.noisep.du_m,
        path.noisep.dl_m,
        path.noisep.fa_g,
        path.noisep.du_g,
        path.noisep.dl_g,
        path.noisep.fam_t,
        path.noisep.du_t,
        path.noisep.dl_t,
        path.snr,
        path.dusn,
        path.dlsn,
        path.snrxx,
        path.sir,
        path.dusi,
        path.dlsi,
        path.rsn,
        path.rt,
        path.rf,
        path.bcr,
        path.ocr,
        path.ocrs,
        path.mir,
        path.probocc,
    ];
    main.into_iter().all(f64::is_finite)
        && path.md_e.iter().chain(path.md_f2.iter()).all(mode_finite)
}

fn mode_finite(mode: &Mode) -> bool {
    mode.bmuf.is_none_or(f64::is_finite)
        && [
            mode.muf90,
            mode.muf50,
            mode.muf10,
            mode.opmuf,
            mode.opmuf10,
            mode.opmuf90,
            mode.fprob,
            mode.hr,
            mode.fs,
            mode.ew,
            mode.ele,
            mode.prw,
            mode.gtw,
            mode.grw,
            mode.tau,
        ]
        .into_iter()
        .all(f64::is_finite)
}

fn output(path: &State<'_>) -> Prediction {
    let digital = path.modulation == ModulationType::Digital;
    let spread = digital && path.t0 != 0.0 && path.f0 != 0.0;
    let short_digital = digital && path.distance <= 9000.0;

    Prediction {
        distance_km: path.distance,
        max_hop_km: path.dmax,
        elevation_deg: path.ele.to_degrees(),
        muf: Muf {
            basic_mhz: path.bmuf,
            median_mhz: path.muf50,
            p10_mhz: path.muf10,
            p90_mhz: path.muf90,
            operational_mhz: path.opmuf,
            operational_p10_mhz: path.opmuf10,
            operational_p90_mhz: path.opmuf90,
        },
        field_strength_dbuv_m: path.ep,
        short_field_strength_dbuv_m: path.es,
        intermediate_field_strength_dbuv_m: path.ei,
        long_field_strength_dbuv_m: path.el,
        rx_power_dbw: path.pr,
        tx_eirp_dbkw: path.txpower + path.gtl,
        tx_gain_dbi: path.gtl,
        rx_gain_dbi: path.grw,
        noise: Noise {
            atmospheric: NoiseComponent {
                median_db: path.noisep.fa_a,
                upper_decile_db: path.noisep.du_a,
                lower_decile_db: path.noisep.dl_a,
            },
            man_made: NoiseComponent {
                median_db: path.noisep.fa_m,
                upper_decile_db: path.noisep.du_m,
                lower_decile_db: path.noisep.dl_m,
            },
            galactic: NoiseComponent {
                median_db: path.noisep.fa_g,
                upper_decile_db: path.noisep.du_g,
                lower_decile_db: path.noisep.dl_g,
            },
            total: NoiseComponent {
                median_db: path.noisep.fam_t,
                upper_decile_db: path.noisep.du_t,
                lower_decile_db: path.noisep.dl_t,
            },
        },
        snr: Snr {
            median_db: path.snr,
            upper_decile_db: path.dusn,
            lower_decile_db: path.dlsn,
            percentile_db: path.snrxx,
        },
        reliability: Reliability {
            basic_percent: path.bcr,
            overall_percent: short_digital.then_some(path.ocr),
            scatter_percent: short_digital.then_some(path.ocrs),
        },
        digital: digital.then(|| DigitalResult {
            sir_db: path.has_interference.then_some(path.sir),
            sir_upper_decile_db: path.has_interference.then_some(path.dusi),
            sir_lower_decile_db: path.has_interference.then_some(path.dlsi),
            snr_percent: spread.then_some(path.rsn),
            time_percent: spread.then_some(path.rt),
            frequency_percent: spread.then_some(path.rf),
            multimode_percent: short_digital.then_some(path.mir),
            scatter_probability: short_digital.then_some(path.probocc),
        }),
        e_modes: core::array::from_fn(|index| mode(path.md_e[index], index)),
        f2_modes: core::array::from_fn(|index| mode(path.md_f2[index], index)),
    }
}

fn mode(mode: Mode, index: usize) -> Option<ModeResult> {
    let basic_muf_mhz = mode.bmuf?;
    Some(ModeResult {
        hops: index as u8 + 1,
        basic_muf_mhz,
        median_muf_mhz: mode.muf50,
        muf_p10_mhz: mode.muf10,
        muf_p90_mhz: mode.muf90,
        operational_muf_mhz: mode.opmuf,
        operational_muf_p10_mhz: mode.opmuf10,
        operational_muf_p90_mhz: mode.opmuf90,
        support_probability: mode.fprob,
        reflection_height_km: mode.hr,
        screening_frequency_mhz: mode.fs,
        field_strength_dbuv_m: mode.ew,
        elevation_deg: mode.ele.to_degrees(),
        rx_power_dbw: mode.prw,
        tx_gain_dbi: mode.gtw,
        rx_gain_dbi: mode.grw,
        delay_ms: mode.tau,
    })
}
