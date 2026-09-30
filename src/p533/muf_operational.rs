use crate::p533::{
    consts::{MAXEMDS, MAXF2MDS},
    control_point::ControlPointId::MP,
    state::State,
};

const ROP: [f64; 6] = [1.20, 1.30, 1.15, 1.25, 1.10, 1.20];

const fn idx(season: usize, time: usize) -> usize {
    season * 2 + time
}

/// Calculates the Operation MUF according to P.533-14 Section 3.7.
pub(super) const fn muf_operational(path: &mut State) {
    // Only do this subroutine for paths less than or equal to 9000 km.
    if path.distance > 9000.0 {
        return;
    }

    // To determine if the time of interest is night or day. It is either day when sunset is greater than sunrise with local time in between and it is not night, otherwise it must be night.
    let time = if (path.controlpt[MP as usize].ltime < path.controlpt[MP as usize].sun.lss)
        && (path.controlpt[MP as usize].ltime > path.controlpt[MP as usize].sun.lsr)
        && !((path.controlpt[MP as usize].ltime > path.controlpt[MP as usize].sun.lss)
            && (path.controlpt[MP as usize].ltime < path.controlpt[MP as usize].sun.lsr))
    {
        Time::Day
    } else {
        Time::Night
    };

    let mut opf2_muf = 0.0;
    let mut opf2_muf10 = 0.0;
    let mut opf2_muf90 = 0.0;

    // Set the OPMUF for the ith F2 mode and determine if it is the largest for all existant F2 modes
    let mut i = 0;
    while i < MAXF2MDS {
        if path.md_f2[i].bmuf.is_some() {
            path.md_f2[i].opmuf =
                path.md_f2[i].muf50 * ROP[idx(path.season as usize, time as usize)];
            path.md_f2[i].opmuf10 = path.md_f2[i].opmuf * path.md_f2[i].deltau;
            path.md_f2[i].opmuf90 = path.md_f2[i].opmuf * path.md_f2[i].deltal;

            // Determine the largest of the OPMUFs.
            if path.md_f2[i].opmuf > opf2_muf {
                opf2_muf = path.md_f2[i].opmuf;
            }
            if path.md_f2[i].opmuf90 > opf2_muf90 {
                opf2_muf90 = path.md_f2[i].opmuf90;
            }
            if path.md_f2[i].opmuf10 > opf2_muf10 {
                opf2_muf10 = path.md_f2[i].opmuf10;
            }
        }
        i += 1;
    }

    let mut ope_muf = 0.0;
    let mut ope_muf10 = 0.0;
    let mut ope_muf90 = 0.0;

    // Set the OPMUF for the ith F2 mode and determine if it is the largest for all existant F2 modes
    let mut i = 0;
    while i < MAXEMDS {
        if let Some(mde_opmuf) = path.md_e[i].bmuf {
            path.md_e[i].opmuf = mde_opmuf;
            path.md_e[i].opmuf10 = path.md_e[i].opmuf * path.md_e[i].deltau;
            path.md_e[i].opmuf90 = path.md_e[i].opmuf * path.md_e[i].deltal;

            // Determine the largest of the OPMUFs.
            if path.md_e[i].opmuf > ope_muf {
                ope_muf = path.md_e[i].opmuf;
            }
            if path.md_e[i].opmuf90 > ope_muf90 {
                ope_muf90 = path.md_e[i].opmuf90;
            }
            if path.md_e[i].opmuf10 > ope_muf10 {
                ope_muf10 = path.md_e[i].opmuf10;
            }
        }
        i += 1;
    }

    path.opmuf = ope_muf.max(opf2_muf);
    path.opmuf90 = ope_muf90.max(opf2_muf90);
    path.opmuf10 = ope_muf10.max(opf2_muf10);
}

#[derive(Clone, Copy)]
enum Time {
    Day,
    Night,
}
