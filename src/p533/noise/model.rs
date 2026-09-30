use crate::{
    data::{DUD, FAKABP, FAKP, FAM},
    input::ManMadeNoise,
    p533::noise::params::NoiseParams,
    p533::{consts::D2R, noise::types::FamStats},
};

use std::f64::consts::{LN_10, PI};

/// Determines atmospheric, man-made, and galactic noise from the associated decile values.
///
pub(in crate::p533) fn noise(
    noisep: &mut NoiseParams,
    hour: usize,
    rlng: f64,
    rlat: f64,
    frequency: f64,
) {
    // Compute each noise component
    atmospheric_noise(noisep, hour, rlng, rlat, frequency);
    galactic_noise(noisep, frequency);
    man_made_noise(noisep, frequency);

    // Constant
    let c = 10.0 / LN_10;
    let sigma_g = 1.56;

    // ----- Upper Decile -----
    let sigma_a = noisep.du_a / 1.282;
    let sigma_m = noisep.du_m / 1.282;

    let alpha_t_u = exp_component(noisep.fa_a, sigma_a, c)
        + exp_component(noisep.fa_g, sigma_g, c)
        + exp_component(noisep.fa_m, sigma_m, c);

    let beta_t_u = beta_component(noisep.fa_a, sigma_a, c)
        + beta_component(noisep.fa_g, sigma_g, c)
        + beta_component(noisep.fa_m, sigma_m, c);

    let gamma_t = (noisep.fa_a / c).exp() + (noisep.fa_g / c).exp() + (noisep.fa_m / c).exp();

    let sigma_t_u = if noisep.du_a > 12.0 || noisep.du_g > 12.0 || noisep.du_m > 12.0 {
        c * (2.0 * (alpha_t_u / gamma_t).ln()).sqrt()
    } else {
        c * (beta_t_u / alpha_t_u.powi(2)).ln_1p().sqrt()
    };

    let fam_tu = c * (alpha_t_u.ln() - (sigma_t_u.powi(2) / (2.0 * c.powi(2))));

    noisep.du_t = 1.282 * sigma_t_u;

    // ----- Lower Decile -----
    let sigma_a = noisep.dl_a / 1.282;
    let sigma_m = noisep.dl_m / 1.282;

    let alpha_t_l = exp_component(noisep.fa_a, sigma_a, c)
        + exp_component(noisep.fa_g, sigma_g, c)
        + exp_component(noisep.fa_m, sigma_m, c);

    let beta_t_l = beta_component(noisep.fa_a, sigma_a, c)
        + beta_component(noisep.fa_g, sigma_g, c)
        + beta_component(noisep.fa_m, sigma_m, c);

    let sigma_t_l = if noisep.dl_a > 12.0 || noisep.dl_g > 12.0 || noisep.dl_m > 12.0 {
        c * (2.0 * (alpha_t_l / gamma_t).ln()).sqrt()
    } else {
        c * (beta_t_l / alpha_t_l.powi(2)).ln_1p().sqrt()
    };

    let fam_tl = c * (alpha_t_l.ln() - (sigma_t_l.powi(2) / (2.0 * c.powi(2))));

    noisep.dl_t = 1.282 * sigma_t_l;

    // Worst-case (minimum) total noise
    noisep.fam_t = fam_tu.min(fam_tl);
}

fn atmospheric_noise(noisep: &mut NoiseParams, hour: usize, rlng: f64, rlat: f64, frequency: f64) {
    let mut fs_now = FamStats::default();
    let mut fs_adj = FamStats::default();

    // Find the local time from the clock UTC, hour, and longitude
    let lrxmt = (hour as i32 + (rlng / (15.0 * D2R)) as i32).rem_euclid(24);

    // Atmospheric noise is determined by finding the atmospheric noise at the current time block at the local mean time,
    // then finding the noise for the adjacent time block, and then iterating between the two noise values
    fs_now.tmblk = (lrxmt / 4) % 6;
    fs_adj.tmblk = (fs_now.tmblk + 1) % 6;

    get_fam_parameters(noisep, &mut fs_now, rlng, rlat, frequency);
    get_fam_parameters(noisep, &mut fs_adj, rlng, rlat, frequency);

    let slp = (f64::from(lrxmt) % 4.0) / 4.0;

    let mut fa = (10.0_f64.powf(fs_adj.fa / 10.0) - 10.0_f64.powf(fs_now.fa / 10.0))
        .mul_add(slp, 10.0_f64.powf(fs_now.fa / 10.0));

    noisep.fa_a = 10.0 * fa.log10();

    fa = (10.0_f64.powf(fs_adj.du / 10.0) - 10.0_f64.powf(fs_now.du / 10.0))
        .mul_add(slp, 10.0_f64.powf(fs_now.du / 10.0));
    noisep.du_a = 10.0 * fa.log10();

    fa = (10.0_f64.powf(fs_adj.dl / 10.0) - 10.0_f64.powf(fs_now.dl / 10.0))
        .mul_add(slp, 10.0_f64.powf(fs_now.dl / 10.0));
    noisep.dl_a = 10.0 * fa.log10();
}

/// Determines the galactic noise in accordance with Section 5 "Man-made noise"
fn galactic_noise(noisep: &mut NoiseParams, frequency: f64) {
    let c = 52.0;
    let d = 23.0;
    noisep.fa_g = c - d * frequency.log10();

    noisep.du_g = 2.0;
    noisep.dl_g = 2.0;
}

fn man_made_noise(noisep: &mut NoiseParams, frequency: f64) {
    let c;
    let d: f64;
    match noisep.man_made_noise {
        ManMadeNoise::City => {
            c = 76.8;
            d = 27.7;
            noisep.du_m = 11.0;
            noisep.dl_m = 6.7;
        }
        ManMadeNoise::Residential => {
            c = 72.5;
            d = 27.7;
            noisep.du_m = 10.6;
            noisep.dl_m = 5.3;
        }
        ManMadeNoise::Rural => {
            c = 67.2;
            d = 27.7;
            noisep.du_m = 9.2;
            noisep.dl_m = 4.6;
        }
        ManMadeNoise::QuietRural => {
            c = 53.6;
            d = 28.6;
            noisep.du_m = 9.2;
            noisep.dl_m = 4.6;
        }
        ManMadeNoise::Noisy => {
            c = 83.2;
            d = 37.5;
            noisep.du_m = 11.0;
            noisep.dl_m = 6.7;
        }
        ManMadeNoise::Quiet => {
            c = 65.2;
            d = 29.1;
            noisep.du_m = 9.2;
            noisep.dl_m = 4.6;
        }
    }

    noisep.fa_m = d.mul_add(-frequency.log10(), c);
}

/// Computes atmospheric noise parameters for the given coordinates and frequency.
fn get_fam_parameters(noisep: &NoiseParams, fs: &mut FamStats, lng: f64, lat: f64, frequency: f64) {
    let month = noisep.month.index();
    let block = fs.tmblk as usize;
    let fakp = &FAKP[month][block];
    let mut v = [0.0_f64; 5];
    let mut u = [0.0_f64; 2];
    let mut cz = 0.0;
    let mut pz = 0.0;
    let mut px = 0.0;
    let mut x;
    let mut y;
    let mut zz = [0.0_f64; 29];

    let lm = 29;
    let ln = 15;

    // Convert longitude to east longitude in [0, 2PI), and halve it
    let mut q = if lng < 0.0 {
        f64::midpoint(lng, 2.0 * PI)
    } else {
        lng / 2.0
    };

    // Longitude Fourier series (longitude component of noise model)
    for (j, zz_j) in zz.iter_mut().enumerate().take(lm) {
        let mut r = 0.0;
        for (k, row) in fakp.iter().enumerate().take(ln) {
            let angle = (k as f64 + 1.0) * q;
            r = angle.sin().mul_add(row[j], r);
        }
        *zz_j = r + fakp[15][j];
    }

    // Latitude Fourier series
    q = lat + PI / 2.0;
    let mut r = 0.0;
    for (j, zz_j) in zz.iter().enumerate().take(lm) {
        r = ((j + 1) as f64 * q).sin().mul_add(*zz_j, r);
    }

    let fakabp = &FAKABP[month][block];
    let fam_1mhz = fakabp[1].mul_add(q, r + fakabp[0]);

    // Determine timeblock index offset based on latitude
    let i = if lat < 0.0 {
        fs.tmblk as usize + 6
    } else {
        fs.tmblk as usize
    };

    // Compute u[1] using: u = (8 * 2^log10(f) - 11) / 4
    u[0] = -0.75;
    u[1] = 8.0f64.mul_add(frequency.log10().exp2(), -11.0) / 4.0;

    for (k, u_k) in u.iter().enumerate() {
        let fam = &FAM[month][i];
        pz = (*u_k).mul_add(fam[0], fam[1]);
        px = (*u_k).mul_add(fam[7], fam[8]);

        for j in 2..7 {
            pz = (*u_k).mul_add(pz, fam[j]);
            px = (*u_k).mul_add(px, fam[j + 7]);
        }

        if k == 0 {
            cz = fam_1mhz * (2.0 - pz) - px;
        }
    }

    fs.fa = cz * pz + px;

    // Limit x to log10(20.0) for Du, Dl, SigmaDu, SigmaDl
    x = frequency.log10();
    if frequency > 20.0 {
        x = 20.0_f64.log10();
    }

    for (j, v_j) in v.iter_mut().enumerate() {
        // For SigmaFam, clamp to 10 MHz
        if j == 4 && frequency > 10.0 {
            x = 10.0_f64.log10();
        }

        let dud = &DUD[month][j][i];
        y = dud[0];
        for coefficient in &dud[1..] {
            y = y.mul_add(x, *coefficient);
        }
        *v_j = y;
    }

    fs.du = v[0];
    fs.dl = v[1];
    fs.sigma_du = v[2];
    fs.sigma_dl = v[3];
    fs.sigma_fam = v[4];
}

// --- Helper functions ---
fn exp_component(fa: f64, sigma: f64, c: f64) -> f64 {
    ((fa / c) + (sigma.powi(2) / (2.0 * c.powi(2)))).exp()
}

fn beta_component(fa: f64, sigma: f64, c: f64) -> f64 {
    let term = ((fa / c) + (sigma.powi(2) / (2.0 * c.powi(2)))).exp();
    term.powi(2) * ((sigma / c).powi(2)).exp_m1()
}
