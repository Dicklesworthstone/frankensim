//! Dimensionless Maugis equations, parameterized by t=a/sqrt(c*c-a*a).
//! At fixed t, cohesive-zone closure is a positive quadratic in c. There is
//! consequently one scalar branch inversion, not a nested two-variable solve.
use super::MaugisError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct State {
    pub t: f64,
    pub a: f64,
    pub c: f64,
    pub l: f64,
    pub d: f64,
    pub p: f64,
    // d'(t)/cos(theta), p'(t)/cos(theta), with finite limits at t=0.
    pub ds: f64,
    pub ps: f64,
    pub dd: f64,
}

/// No subtraction of O(theta) or O(1) terms in a thin cohesive annulus.
fn closure_terms(theta: f64) -> (f64, f64) {
    if theta < 0.125 {
        let z = theta * theta;
        let a = theta
            * z
            * (4.0 / 3.0
                + z * (-8.0 / 15.0
                    + z * (8.0 / 105.0 + z * (-16.0 / 2835.0 + z * (8.0 / 31185.0)))));
        let b = z
            * (0.5
                + z * (-0.125
                    + z * (1.0 / 144.0
                        + z * (-1.0 / 5760.0 + z * (1.0 / 403200.0 - z / 43545600.0)))));
        (a, b)
    } else {
        let (s, q) = theta.sin_cos();
        ((s * s - q * q) * theta + s * q, s * theta - 1.0 + q)
    }
}

pub(super) fn state(lambda: f64, t: f64) -> Result<State, MaugisError> {
    if !t.is_finite() || t < 0.0 {
        return Err(MaugisError::Unrepresentable {
            field: "contact branch coordinate",
        });
    }
    let h = t.hypot(1.0);
    let s = 1.0 / h;
    let q = t / h;
    let theta = 1.0_f64.atan2(t);
    let (aa, bb) = closure_terms(theta);
    let alpha = 0.5 * lambda * aa;
    let beta = (4.0 / 3.0) * lambda * lambda * bb;
    let c = 2.0 / (beta + beta.hypot(2.0 * alpha.sqrt()));
    let a = c * q;
    let l = c * s;
    let d = a * a - (4.0 / 3.0) * lambda * l;
    let p = a * a * a - lambda * c * c * (q * s + theta);
    let denom = 2.0 * alpha * c + beta;
    let z = (2.0 * lambda * theta * s * s * s * c * c
        + (4.0 / 3.0) * lambda * lambda * theta * s * s * c)
        / denom;
    let da = q * q * z + c * s * s * s;
    let ds = 2.0 * c * da - (4.0 / 3.0) * lambda * s * (z - c * s);
    let ps = 3.0 * a * c * da - 2.0 * lambda * c * z * (q * s + theta)
        + 2.0 * lambda * c * c * q * s * s;
    if [c, a, l, d, p, ds, ps].iter().any(|v| !v.is_finite()) || c <= 0.0 || l <= 0.0 {
        return Err(MaugisError::Unrepresentable {
            field: "dimensionless Maugis response",
        });
    }
    Ok(State {
        t,
        a,
        c,
        l,
        d,
        p,
        ds,
        ps,
        dd: q * ds,
    })
}

pub(super) fn fold(lambda: f64, force: bool) -> Result<f64, MaugisError> {
    let slope = |s: State| if force { s.ps } else { s.ds };
    if slope(state(lambda, 0.0)?) >= 0.0 {
        return Ok(0.0);
    }
    let mut hi = 1.0;
    for _ in 0..128 {
        if slope(state(lambda, hi)?) > 0.0 {
            break;
        }
        hi *= 2.0;
    }
    if slope(state(lambda, hi)?) <= 0.0 {
        return Err(MaugisError::Unrepresentable {
            field: "stable branch fold",
        });
    }
    let mut lo = 0.0;
    for _ in 0..128 {
        let mid = 0.5 * lo + 0.5 * hi;
        if mid == lo || mid == hi {
            break;
        }
        if slope(state(lambda, mid)?) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // The stable side, including at a representationally unresolved fold.
    Ok(hi)
}

pub(super) fn invert(
    lambda: f64,
    mut lo: f64,
    mut hi: f64,
    target: f64,
    force: bool,
) -> Result<State, MaugisError> {
    let value = |s: State| if force { s.p } else { s.d };
    let left = state(lambda, lo)?;
    let right = state(lambda, hi)?;
    if target == value(left) {
        return Ok(left);
    }
    if target == value(right) {
        return Ok(right);
    }
    for _ in 0..128 {
        let mid = 0.5 * lo + 0.5 * hi;
        if mid == lo || mid == hi {
            break;
        }
        let s = state(lambda, mid)?;
        if value(s) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let a = state(lambda, lo)?;
    let b = state(lambda, hi)?;
    Ok(if (value(a) - target).abs() < (value(b) - target).abs() {
        a
    } else {
        b
    })
}
