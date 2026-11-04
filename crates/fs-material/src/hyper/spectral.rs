//! Second-order spectral lift of tr g(C). The Daleckii–Krein trace
//! derivative is D tr g(C)[H] = sum g'(d_i) H_ii; its Hessian is
//! sum g'[d_i,d_j] H_ij H_ji in an orthonormal eigenbasis. At equal roots
//! the divided difference is g''(d), with no eigenvector derivative.
//! See Simo & Taylor, CMAME 85 (1991), 273–310, doi:10.1016/0045-7825(91)90100-F.
//!
//! C-C.value() has exactly zero primal. Its degree-two lift therefore
//! carries the exact value/gradient/Hessian, including the chain rule for
//! C(F). It is not an approximation in stretch and has no finite-difference
//! step. Third and higher AD derivatives are outside this interface's claim.

use super::OgdenTerm;
use crate::MaterialError;
use fs_ad::Real;
use fs_math::det;

pub(super) fn principal_stretches(f: &[f64; 9]) -> Result<[f64; 3], MaterialError> {
    let c: [f64; 9] =
        core::array::from_fn(|ij| (0..3).map(|r| f[3 * r + ij / 3] * f[3 * r + ij % 3]).sum());
    let (d, _) = fs_la::eigen::jacobi_eigh(&c, 3);
    if d.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(MaterialError::State {
            what: "non-positive/nonfinite Cauchy–Green eigenvalue".to_string(),
        });
    }
    Ok(core::array::from_fn(|i| det::sqrt(d[i])))
}

// ln(1+t) without losing the small difference when a and b coalesce.
// For |t|<1/8, |s|=|t/(2+t)|<1/15 and the omitted atanh tail is <1e-25.
fn log1p(t: f64) -> f64 {
    if t.abs() >= 0.125 {
        return det::ln(1.0 + t);
    }
    let s = t / (2.0 + t);
    let mut power = s;
    let mut sum = s;
    for n in 1..10 {
        power *= s * s;
        sum += power / f64::from(2 * n + 1);
    }
    2.0 * sum
}

fn derivative_difference(a: f64, b: f64, r: f64) -> f64 {
    let (a, b) = (a.max(b), a.min(b));
    let t = (a - b) / b;
    let base = r * det::exp((r - 2.0) * det::ln(b));
    if t == 0.0 {
        base * (r - 1.0)
    } else {
        base * det::expm1((r - 1.0) * log1p(t)) / t
    }
}

// Original isochoric Ogden potential near a rigid motion, without forming
// unit-size eigenvalues. Write C=m(I+B), tr(B)=0, then
// L=dev(log(I+B))=log(Cbar). The energy is
// sum mu/alpha * tr(exp(alpha L/2)-I-alpha L/2).
// Every retained term is at least quadratic in strain, including for alpha
// near zero. Both bounded series are evaluated on Real, so the same potential
// supplies the stress and consistent tangent; this is not linearization.
pub(super) fn small_ogden_energy<T: Real>(a: &[T; 9], terms: &[OgdenTerm]) -> Option<T> {
    let three = T::from_f64(3.0);
    let mean = T::one() + (a[0] + a[4] + a[8]) / three;
    let mut b = a.map(|entry| entry / mean);
    let denominator = three * mean;
    b[0] = ((a[0] - a[4]) + (a[0] - a[8])) / denominator;
    b[4] = ((a[4] - a[0]) + (a[4] - a[8])) / denominator;
    b[8] = -b[0] - b[4];
    let norm = (0..3)
        .map(|i| (0..3).map(|j| b[3 * i + j].value().abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    let exponent = terms
        .iter()
        .map(|p| 0.5 * p.alpha.abs())
        .fold(0.0_f64, f64::max);
    // ||dev log(I+B)||_inf <= 2 ||B||_inf/(1-||B||_inf).
    // The bounds apply to the entire original material, not a sampled
    // eigenvalue; large exponents retain the existing spectral fallback.
    if norm > 1.0 / 16.0 || exponent * norm / (1.0 - norm) > 1.0 / 16.0 {
        return None;
    }
    let square: [T; 9] = core::array::from_fn(|ij| {
        (0..3).fold(T::zero(), |sum, k| {
            b[3 * (ij / 3) + k].mul_add(b[3 * k + ij % 3], sum)
        })
    });
    let q = -T::from_f64(0.5) * (square[0] + square[4] + square[8]);
    let cubic = super::det3(&b);
    // Cayley-Hamilton: B^3=-q B+det(B) I. Retain only the two
    // coefficients needed by dev(log(I+B)); no eigensolver or allocation.
    let mut power = [T::zero(), T::one(), T::zero()];
    let (mut linear, mut quadratic) = (T::zero(), T::zero());
    for degree in 1..=24 {
        let coefficient = T::from_f64(if degree % 2 == 0 { -1.0 } else { 1.0 } / f64::from(degree));
        linear = coefficient.mul_add(power[1], linear);
        quadratic = coefficient.mul_add(power[2], quadratic);
        power = [cubic * power[2], power[0] - q * power[2], power[1]];
    }
    // The log tail and its first two derivatives are below 1e-26 for the
    // admitted norm. Subtract the trace by differences, not unit quantities.
    let mut log: [T; 9] = core::array::from_fn(|ij| linear * b[ij] + quadratic * square[ij]);
    for (i, j, k) in [(0, 4, 8), (4, 0, 8)] {
        log[i] =
            linear * b[i] + quadratic * ((square[i] - square[j]) + (square[i] - square[k])) / three;
    }
    log[8] = -log[0] - log[4];
    let trace2 = log[0] * log[0]
        + log[4] * log[4]
        + log[8] * log[8]
        + T::from_f64(2.0) * (log[1] * log[1] + log[2] * log[2] + log[5] * log[5]);
    let trace3 = three * super::det3(&log);
    let mut energy = T::zero();
    for term in terms {
        let r = T::from_f64(0.5 * term.alpha);
        let h = log.map(|entry| r * entry);
        let q = -T::from_f64(0.5)
            * (h[0] * h[0]
                + h[4] * h[4]
                + h[8] * h[8]
                + T::from_f64(2.0) * (h[1] * h[1] + h[2] * h[2] + h[5] * h[5]));
        let cubic = super::det3(&h);
        // s_n=r^(n-2) tr(L^n), s_1=0. Factoring out mu*alpha/4
        // preserves the alpha->0 limit without subtracting O(alpha) terms
        // or first squaring a possibly very small alpha.
        let mut history = [T::zero(), trace2, r * trace3];
        let mut sum = trace2 / T::from_f64(2.0) + history[2] / T::from_f64(6.0);
        let mut coefficient = 1.0 / 6.0;
        for degree in 4..=16 {
            let next = -q * history[1] + cubic * history[0];
            coefficient /= f64::from(degree);
            sum = T::from_f64(coefficient).mul_add(next, sum);
            history = [history[1], history[2], next];
        }
        // ||alpha L/2||_inf <= 1/8; the exponential tail and its first
        // two derivatives in that argument are below 1e-25 before the modulus.
        energy = T::from_f64(0.25 * (term.mu_pa * term.alpha)).mul_add(sum, energy);
    }
    Some(energy)
}

pub(super) fn ogden_energy<T: Real>(c: &[T; 9], terms: &[OgdenTerm]) -> T {
    let primal = c.map(Real::value);
    let (d, v) = fs_la::eigen::jacobi_eigh(&primal, 3);
    let delta: [T; 9] = core::array::from_fn(|ij| c[ij] - T::from_f64(primal[ij]));
    let h: [T; 9] = core::array::from_fn(|ij| {
        let mut value = T::zero();
        for p in 0..3 {
            for q in 0..3 {
                value =
                    value + T::from_f64(v[3 * p + ij / 3] * v[3 * q + ij % 3]) * delta[3 * p + q];
            }
        }
        value
    });
    let mut w = T::zero();
    for term in terms {
        let r = 0.5 * term.alpha;
        // This algebraic identity also preserves the NH reduction without
        // eigensolver roundoff entering its energy or derivatives.
        let trace = if r == 1.0 {
            c[0] + c[4] + c[8]
        } else {
            let mut trace = T::zero();
            for i in 0..3 {
                trace = trace
                    + T::from_f64(det::exp(r * det::ln(d[i])))
                    + T::from_f64(r * det::exp((r - 1.0) * det::ln(d[i]))) * h[3 * i + i];
                for j in 0..3 {
                    trace = trace
                        + T::from_f64(0.5 * derivative_difference(d[i], d[j], r))
                            * h[3 * i + j]
                            * h[3 * j + i];
                }
            }
            trace
        };
        w = w + T::from_f64(term.mu_pa / term.alpha) * (trace - T::from_f64(3.0));
    }
    w
}
