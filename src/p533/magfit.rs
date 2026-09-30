use crate::p533::control_point::ControlPt;

const G: [f64; 49] = [
    0.000000, 0.304112, 0.024035, -0.031518, -0.041794, 0.016256, -0.019523, 0.000000, 0.021474,
    -0.051253, 0.062130, -0.045298, -0.034407, -0.004853, 0.000000, 0.000000, -0.013381, -0.024898,
    -0.021795, -0.019447, 0.003212, 0.000000, 0.000000, 0.000000, -0.006496, 0.007008, -0.000608,
    0.021413, 0.000000, 0.000000, 0.000000, 0.000000, -0.002044, 0.002775, 0.001051, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000697, 0.000227, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.001115,
];

const H: [f64; 49] = [
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.057989,
    0.033124, 0.014870, -0.011825, -0.000796, -0.005758, 0.000000, 0.000000, -0.001579, -0.004075,
    0.010006, -0.002000, -0.008735, 0.000000, 0.000000, 0.000000, 0.000210, 0.000430, 0.004597,
    -0.003406, 0.000000, 0.000000, 0.000000, 0.000000, 0.001385, 0.002421, -0.000118, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, -0.001218, -0.001116, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, -0.000325,
];

const CT: [f64; 49] = [
    0.0000000,
    0.0000000,
    0.33333333,
    0.266666666,
    0.25714286,
    0.25396825,
    0.25252525,
    0.0000000,
    0.0000000,
    0.00000000,
    0.200000000,
    0.22857142,
    0.23809523,
    0.24242424,
    0.0000000,
    0.0000000,
    0.00000000,
    0.000000000,
    0.14285714,
    0.19047619,
    0.21212121,
    0.0000000,
    0.0000000,
    0.00000000,
    0.000000000,
    0.00000000,
    0.11111111,
    0.16161616,
    0.0000000,
    0.0000000,
    0.00000000,
    0.000000000,
    0.00000000,
    0.00000000,
    0.09090909,
    0.0000000,
    0.0000000,
    0.00000000,
    0.000000000,
    0.00000000,
    0.00000000,
    0.00000000,
    0.0000000,
    0.0000000,
    0.00000000,
    0.000000000,
    0.00000000,
    0.00000000,
    0.00000000,
];

/// For easy access to flattened 2D arrays from the original C code.
const fn idx(row: usize, col: usize) -> usize {
    row * 7 + col
}

/// Calculates the magnetic dip and gyrofrequency, as described in Section 2 of Rec P.1239. Some readability is sacrificed here for the sake of better performance.
pub(super) fn magfit<const HEIGHT: usize>(here: &mut ControlPt) {
    let mut p = [[0.0; 7]; 7];
    let mut dp = [[0.0; 7]; 7];

    p[0][0] = 1.0;

    // P.1239-4 equation (8): radius used by the 1960 magnetic field model.
    let ar = 6371.2 / (6371.2 + HEIGHT as f64);
    let mut ar_pow = ar * ar * ar;

    let sin_lat = here.l.lat.sin();
    let cos_lat = here.l.lat.cos();
    let cos_lat_inv = 1.0 / cos_lat;

    // Each longitude harmonic is shared by every degree with n >= m.
    let harmonics: [(f64, f64); 7] = core::array::from_fn(|m| (m as f64 * here.l.lng).sin_cos());

    let mut fx = 0.0;
    let mut fy = 0.0;
    let mut fz = 0.0;
    let mut n = 1;
    while n <= 6 {
        let mut sumz = 0.0;
        let mut sumx = 0.0;
        let mut sumy = 0.0;
        let mut m = 0;
        while m <= n {
            let m_f64 = m as f64;

            if n == m {
                p[m][n] = cos_lat * p[m - 1][n - 1];
                dp[m][n] = sin_lat.mul_add(p[m - 1][n - 1], cos_lat * dp[m - 1][n - 1]);
            } else if n != 1 {
                let i = idx(m, n);
                p[m][n] = CT[i].mul_add(-p[m][n - 2], sin_lat * p[m][n - 1]);
                dp[m][n] = CT[i].mul_add(
                    -dp[m][n - 2],
                    cos_lat.mul_add(-p[m][n - 1], sin_lat * dp[m][n - 1]),
                );
            } else {
                p[m][n] = sin_lat * p[m][n - 1];
                dp[m][n] = cos_lat.mul_add(-p[m][n - 1], sin_lat * dp[m][n - 1]);
            }

            let (sin_m_lng, cos_m_lng) = harmonics[m];

            let i = idx(m, n);
            let g = G[i];
            let h = H[i];

            let common = h.mul_add(sin_m_lng, g * cos_m_lng);
            let common_y = h.mul_add(-cos_m_lng, g * sin_m_lng);

            sumz = p[m][n].mul_add(common, sumz);
            sumx = dp[m][n].mul_add(common, sumx);
            sumy = (m_f64 * p[m][n]).mul_add(common_y, sumy);
            m += 1;
        }

        let n_f64 = n as f64;
        fz = (ar_pow * (n_f64 + 1.0)).mul_add(sumz, fz);
        fx = ar_pow.mul_add(-sumx, fx);
        fy = ar_pow.mul_add(sumy, fy);

        ar_pow *= ar;
        n += 1;
    }

    let fx_sq = fx * fx;
    let fy_adj = fy * cos_lat_inv;
    let horiz_mag = (fx_sq + fy_adj * fy_adj).sqrt();

    let dip_val = (fz / horiz_mag).atan();
    let fh_val = 2.8 * (fx_sq + fy_adj * fy_adj + fz * fz).sqrt();

    match HEIGHT {
        100 => {
            here.dip.0 = dip_val;
            here.fh.0 = fh_val;
        }
        300 => {
            here.dip.1 = dip_val;
            here.fh.1 = fh_val;
        }
        _ => unreachable!("magnetic field height must be 100 or 300 km"),
    }
}
