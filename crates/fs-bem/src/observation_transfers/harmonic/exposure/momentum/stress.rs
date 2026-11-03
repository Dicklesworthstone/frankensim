//! Compensated spatial bilinear wave-stress contraction, before tone storage.

use super::*;

pub(super) fn spatial_pair(
    a: &Line,
    b: &Line,
    reference: &PreparedHarmonicObservation,
    quadrature: &[HarmonicSurfaceQuadrature],
    origin: [f64; 3],
    real: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WaveMomentumPair, TransferError> {
    let density = Dd::from_f64(reference.medium.density);
    let speed = Dd::from_f64(reference.medium.sound_speed);
    let half = Dd::from_f64(0.5);
    let mut difference = [[Dd::ZERO; 2]; 6];
    let mut sum = [[Dd::ZERO; 2]; 6];
    for sample in quadrature {
        poll(cancelled)?;
        // Exact zero areas still passed admission and count against budgets.
        // Their field products and moment arms cannot contribute to this flux.
        if sample.oriented_area_m2.iter().all(|&area| area == 0.0) {
            continue;
        }
        let af = &a.fields[sample.receiver];
        let bf = &b.fields[sample.receiver];
        // Divide pressure before squaring, avoiding a needless c² intermediate.
        let ap = component(af, 0).map(|x| x / speed);
        let bp = component(bf, 0).map(|x| x / speed);
        let (mut isotropic_d, mut isotropic_s) = products(ap, bp, real);
        for value in isotropic_d.iter_mut().chain(&mut isotropic_s) {
            *value = *value / density;
        }
        let area = sample.oriented_area_m2.map(Dd::from_f64);
        let mut projected_b = [Dd::ZERO; 2];
        for i in 0..3 {
            let au = component(af, i + 1);
            let bu = component(bf, i + 1);
            let (d, s) = products(au, bu, real);
            for part in 0..2 {
                isotropic_d[part] = isotropic_d[part] - density * d[part];
                if real {
                    isotropic_s[part] = isotropic_s[part] - density * s[part];
                }
                projected_b[part] = projected_b[part] + area[i] * bu[part];
            }
        }
        let mut traction_d = [[Dd::ZERO; 2]; 3];
        let mut traction_s = [[Dd::ZERO; 2]; 3];
        for i in 0..3 {
            let (d, s) = products(component(af, i + 1), projected_b, real);
            for part in 0..2 {
                traction_d[i][part] = density * d[part] + half * isotropic_d[part] * area[i];
                if real {
                    traction_s[i][part] = density * s[part] + half * isotropic_s[part] * area[i];
                }
            }
        }
        let position = reference.receivers[sample.receiver].position_m;
        let arm: [Dd; 3] =
            core::array::from_fn(|i| Dd::from_f64(position[i]) - Dd::from_f64(origin[i]));
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            for part in 0..2 {
                difference[i][part] = difference[i][part] + traction_d[i][part];
                difference[i + 3][part] = difference[i + 3][part] + arm[j] * traction_d[k][part]
                    - arm[k] * traction_d[j][part];
                if real {
                    sum[i][part] = sum[i][part] + traction_s[i][part];
                    sum[i + 3][part] = sum[i + 3][part] + arm[j] * traction_s[k][part]
                        - arm[k] * traction_s[j][part];
                }
            }
        }
    }
    // Re(a e^-iwt) Re(b e^-ivt) gives half a*conj(b) at w-v,
    // and half a*b at w+v. A proper factor has no unconjugated term.
    for value in difference.iter_mut().chain(&mut sum).flatten() {
        *value = *value * half;
        if !finite(*value) {
            return Err(bad("acoustic wave-stress contraction is unrepresentable"));
        }
    }
    Ok((difference, sum))
}

fn products(a: [Dd; 2], b: [Dd; 2], real: bool) -> ([Dd; 2], [Dd; 2]) {
    let rr = a[0] * b[0];
    let ii = a[1] * b[1];
    let ri = a[0] * b[1];
    let ir = a[1] * b[0];
    let difference = [rr + ii, ir - ri];
    let sum = if real {
        [rr - ii, ri + ir]
    } else {
        [Dd::ZERO; 2]
    };
    (difference, sum)
}
