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
