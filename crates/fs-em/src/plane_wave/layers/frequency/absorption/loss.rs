//! Analytic tangents of the original slab's positive even/odd loss channels.
//!
//! D denotes omega*d/domega. Numerical common scaling is held fixed; physical
//! attenuation, phase, port admittance and material loss all retain their jets.
//! The native loss eigenvalues remain authoritative. No small loss is recovered
//! by subtracting reflected/transmitted power or dividing by a loss eigenvalue.

use fs_math::{c64::C64, dd::Dd, det};

#[derive(Clone, Copy)]
pub(super) struct SlabLossDirection {
    pub(super) normal: C64,
    pub(super) normal_scaled_derivative: C64,
    pub(super) permittivity_ratio: C64,
    pub(super) permittivity_ratio_scaled_derivative: C64,
    pub(super) permeability_ratio: f64,
    pub(super) cosine: f64,
    pub(super) cosine_scaled_derivative: f64,
    pub(super) depth: f64,
    pub(super) depth_scaled_derivative: f64,
    pub(super) phase_parts: Dd,
    pub(super) loss: [f64; 2],
    pub(super) transverse_electric: bool,
}

pub(super) fn scaled_loss_derivatives(input: SlabLossDirection) -> Option<[f64; 2]> {
    let SlabLossDirection {
        normal: q,
        normal_scaled_derivative: dq,
        permittivity_ratio: epsilon,
        permittivity_ratio_scaled_derivative: depsilon,
        permeability_ratio: mu,
        cosine,
        cosine_scaled_derivative: dcosine,
        depth,
        depth_scaled_derivative: ddepth,
        phase_parts,
        loss,
        transverse_electric: te,
    } = input;
    if [q, dq, epsilon, depsilon].iter().any(|&z| !finite(z))
        || [
            mu,
            cosine,
            dcosine,
            depth,
            ddepth,
            phase_parts.hi,
            phase_parts.lo,
        ]
        .iter()
        .any(|x| !x.is_finite())
        || mu <= 0.0
        || cosine <= 0.0
        || cosine > 1.0
        || depth < 0.0
        || q.re < 0.0
        || q.im > 0.0
        || epsilon.im > 0.0
        || loss.iter().any(|x| !x.is_finite() || *x < 0.0)
    {
        return None;
    }
    if depth == 0.0 {
        return (ddepth == 0.0 && loss == [0.0; 2]).then_some([0.0; 2]);
    }
    if q == C64::ZERO {
        return None;
    }
    let (s, ds) = if te {
        (C64::from_re(mu * cosine), C64::from_re(mu * dcosine))
    } else {
        (
            epsilon.scale(cosine),
            depsilon.scale(cosine) + epsilon.scale(dcosine),
        )
    };
    let scale = maximum(s).max(maximum(q));
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    let (s, ds, normal, dnormal) = (
        divide(s, scale),
        divide(ds, scale),
        divide(q, scale),
        divide(dq, scale),
    );
    let (a, b, da, db) = if te {
        (s, normal, ds, dnormal)
    } else {
        (normal, s, dnormal, ds)
    };
    if [a, b, da, db].iter().any(|&z| !finite(z)) {
        return None;
    }
    let phase = q.re * depth;
    let decay = q.im * depth;
    let dphase = sum([(dq.re, depth), (q.re, ddepth)]);
    let ddecay = sum([(dq.im, depth), (q.im, ddepth)]);
    if [phase, decay, dphase, ddecay]
        .iter()
        .any(|x| !x.is_finite())
    {
        return None;
    }
    let attenuation = det::exp(decay);
    // The caller can still publish the opaque native carrier. Its missing
    // transmitted frequency jet is not silently replaced by zero here.
    if attenuation == 0.0 || !attenuation.is_finite() {
        return None;
    }
    let (sine, cosine_phase) = det::sin_cos_dd(phase_parts);
    let (sine, cosine_phase) = (sine.to_f64(), cosine_phase.to_f64());
    let p = C64::new(attenuation * cosine_phase, -attenuation * sine);
    let dattenuation = attenuation * ddecay;
    let phase_rotation = attenuation * dphase;
    let dp = C64::new(
        sum([(dattenuation, cosine_phase), (-phase_rotation, sine)]),
        -sum([(dattenuation, sine), (phase_rotation, cosine_phase)]),
    );
    let dwork = if phase.abs().max(decay.abs()) < 0.1 {
        // Differentiate the same small-depth polynomial before multiplying its
        // positive loss factor. In particular the odd O(depth^3) work survives.
        let x = phase * phase;
        let y = decay * decay;
        let dx = (phase * dphase) * 2.0;
        let dy = (decay * ddecay) * 2.0;
        let f = x * (1.0 / 6.0 - x * (1.0 / 120.0 - x * (1.0 / 5040.0 - x / 362880.0)));
        let g = y * (1.0 / 6.0 + y * (1.0 / 120.0 + y * (1.0 / 5040.0 + y / 362880.0)));
        let df = (1.0 / 6.0 - x * (1.0 / 60.0 - x * (1.0 / 1680.0 - x / 90720.0))) * dx;
        let dg = (1.0 / 6.0 + y * (1.0 / 60.0 + y * (1.0 / 1680.0 + y / 90720.0))) * dy;
        let plus = (2.0 + g - f, dg - df);
        let minus = (g + f, dg + df);
        if te {
            let material = (-epsilon.im * mu / scale, -depsilon.im * mu / scale);
            let factor = product(
                product(
                    product((attenuation, dattenuation), (a.re, da.re)),
                    material,
                ),
                (depth, ddepth),
            );
            [product(factor, plus).1, product(factor, minus).1]
        } else {
            let material = (
                -epsilon.im * cosine / scale,
                -sum([(depsilon.im, cosine), (epsilon.im, dcosine)]) / scale,
            );
            let factor = product(
                product((attenuation, dattenuation), material),
                (depth, ddepth),
            );
            let weight = (
                sum([(q.re, normal.re), (q.im, normal.im)]),
                2.0 * sum([(dq.re, normal.re), (dq.im, normal.im)]),
            );
            let lateral = (
                (1.0 - cosine) * (1.0 + cosine) / scale,
                -2.0 * cosine * dcosine / scale,
            );
            let even = add(product(weight, plus), product(lateral, minus));
            let odd = add(product(weight, minus), product(lateral, plus));
            [product(factor, even).1, product(factor, odd).1]
        }
    } else {
        let cross = a * b.conj();
        let dcross = complex_sum(&[(da, b.conj()), (a, db.conj())]);
        let power_loss = -det::expm1(2.0 * decay);
        // Form rho*(rho*Ddecay), not exp(2*decay)*Ddecay: the latter
        // can underflow before a large finite source derivative restores it.
        let dpower_loss = -2.0 * attenuation * dattenuation;
        let common = sum([(dcross.re, power_loss), (cross.re, dpower_loss)]);
        let interference = 2.0 * sum([(dcross.im, p.im), (cross.im, dp.im)]);
        [common - interference, common + interference]
    };
    let plus = C64::ONE + p;
    let half_sine = det::sin_cos_dd(phase_parts * Dd::from_f64(0.5)).0.to_f64();
    let minus = C64::new(
        -det::expm1(decay) + 2.0 * attenuation * half_sine * half_sine,
        attenuation * sine,
    );
    let z = [a * plus + b * minus, a * minus + b * plus];
    let dz = [
        complex_sum(&[(da, plus), (a, dp), (db, minus), (-b, dp)]),
        complex_sum(&[(da, minus), (-a, dp), (db, plus), (b, dp)]),
    ];
    let mut output = [0.0; 2];
    for i in 0..2 {
        let denominator = z[i].norm_sq();
        let ddenominator = 2.0 * sum([(z[i].re, dz[i].re), (z[i].im, dz[i].im)]);
        if !denominator.is_finite()
            || denominator <= 0.0
            || !ddenominator.is_finite()
            || !dwork[i].is_finite()
        {
            return None;
        }
        let quarter = 0.25 * denominator;
        let work_term = if quarter > 0.0 {
            dwork[i] / quarter
        } else {
            4.0 * (dwork[i] / denominator)
        };
        output[i] = sum([(1.0, work_term), (-loss[i], ddenominator / denominator)]);
        if !output[i].is_finite() {
            return None;
        }
    }
    Some(output)
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn maximum(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}
fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}
fn sum<const N: usize>(terms: [(f64, f64); N]) -> f64 {
    terms
        .into_iter()
        .fold(Dd::ZERO, |value, (a, b)| {
            value + Dd::from_f64(a) * Dd::from_f64(b)
        })
        .to_f64()
}
fn product(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0, sum([(a.1, b.0), (a.0, b.1)]))
}
fn add(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (sum([(a.0, 1.0), (b.0, 1.0)]), sum([(a.1, 1.0), (b.1, 1.0)]))
}
fn complex_sum(terms: &[(C64, C64)]) -> C64 {
    let (mut real, mut imaginary) = (Dd::ZERO, Dd::ZERO);
    for &(a, b) in terms {
        real = real + Dd::from_f64(a.re) * Dd::from_f64(b.re)
            - Dd::from_f64(a.im) * Dd::from_f64(b.im);
        imaginary = imaginary
            + Dd::from_f64(a.re) * Dd::from_f64(b.im)
            + Dd::from_f64(a.im) * Dd::from_f64(b.re);
    }
    C64::new(real.to_f64(), imaginary.to_f64())
}

#[cfg(test)]
mod tests;
