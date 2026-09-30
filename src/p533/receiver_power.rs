use crate::p533::{
    antenna::{Direction, antenna_gain_08},
    consts::{MAXEMDS, MAXF2MDS, TINY_DB},
    state::State,
};

/// Calculates the median available receiver power by the method in ITU-R P.533-14 Section 6.
pub(super) fn median_available_receiver_power(path: &mut State) {
    let mut sum_pr = 0.0;
    let mut max_power = TINY_DB;
    let elevation = 2.0 * std::f64::consts::PI;

    if path.distance <= 7000.0 {
        if let Some(n0_e) = path.n0_e {
            for i in n0_e..MAXEMDS {
                let mode = path.md_e[i];

                if (i == n0_e && path.distance / (n0_e + 1) as f64 <= 2000.0)
                    || (i != n0_e && mode.bmuf.is_some())
                {
                    if max_power < mode.prw {
                        max_power = mode.prw;
                        path.ele = mode.ele;
                        path.gtl = mode.gtw;
                        path.grw = mode.grw;
                    }
                    sum_pr += 10f64.powf(mode.prw / 10.0);
                }
            }
        }

        if let Some(n0_f2) = path.n0_f2 {
            for i in n0_f2..MAXF2MDS {
                let mode = path.md_f2[i];
                let is_first_mode = i == n0_f2;
                let valid_hop = path.distance / (n0_f2 + 1) as f64 <= path.dmax;
                let has_valid_bmuf = mode.bmuf.is_some();
                let below_fs = mode.fs < path.frequency;

                if below_fs && ((is_first_mode && valid_hop) || (!is_first_mode && has_valid_bmuf))
                {
                    if max_power < mode.prw {
                        max_power = mode.prw;
                        path.ele = mode.ele;
                        path.gtl = mode.gtw;
                        path.grw = mode.grw;
                    }
                    sum_pr += 10f64.powf(mode.prw / 10.0);
                }
            }
        }

        // Final total power
        path.pr = if sum_pr != 0.0 {
            10.0 * sum_pr.log10()
        } else {
            TINY_DB
        };

        // Use Es for Ep at short distances
        path.ep = path.es;
    } else if path.distance < 9000.0 {
        let mut long_path_elevation = elevation;
        let grw = antenna_gain_08(
            path,
            path.antenna.1,
            Direction::RxToTx,
            &mut long_path_elevation,
        );
        path.pr = 20.0f64.mul_add(-path.frequency.log10(), path.ei + grw) - 107.2;
        path.grw = grw;
        path.ep = path.ei;
        path.ele = long_path_elevation;
    } else {
        let mut long_path_elevation = elevation;
        let grw = antenna_gain_08(
            path,
            path.antenna.1,
            Direction::RxToTx,
            &mut long_path_elevation,
        );
        path.pr = 20.0f64.mul_add(-path.frequency.log10(), path.el + grw) - 107.2;
        path.grw = grw;
        path.ep = path.el;
        path.ele = long_path_elevation;
    }
}
