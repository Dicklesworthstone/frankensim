//! One compression realization for positive caloric capacity/conductance laws.
//!
//! For capacity-normalized conductance A and pressure forcing b, elimination
//! of the ideal-gas equation of state gives
//! K(s)=p+p*alpha*s*b' [s*(I-alpha*b*b')+A]^-1 b.
//! The rank-one capacity correction is positive when alpha*b'b < 1.
//! A symmetric congruence then gives positive Maxwell arms, without a fit or
//! an independent thermal integrator. Both prescribed and material walls use it.

use fs_la::eigen::jacobi_eigh;
use fs_math::{c64::C64, det};

pub(super) const MAX_COORDINATES: usize = 128;

#[derive(Clone, Copy, Debug)]
pub(super) struct BulkArm {
    pub increment_pa: f64,
    pub relaxation_time_s: f64,
}

#[derive(Debug)]
pub(super) struct Spectrum {
    pub arms: Vec<BulkArm>,
    pub frozen_bulk_modulus_pa: f64,
}

impl Spectrum {
    pub(super) fn bulk_at(&self, pressure: f64, omega: f64) -> C64 {
        let s = C64::new(0.0, omega);
        self.arms.iter().fold(C64::from_re(pressure), |sum, arm| {
            sum + (s / (s + C64::from_re(arm.relaxation_time_s.recip()))).scale(arm.increment_pa)
        })
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), &'static str> {
    if cancelled() {
        Err("cavity caloric preparation cancelled")
    } else {
        Ok(())
    }
}

/// All coordinates are retained. This cold solve has a hard dimension bound;
/// positivity and eigen residuals are checked before any runtime is published.
pub(super) fn realize(
    pressure: f64,
    alpha: f64,
    mut a: Vec<f64>,
    b: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Spectrum, &'static str> {
    poll(cancelled)?;
    let n = b.len();
    if n == 0
        || n > MAX_COORDINATES
        || a.len() != n * n
        || !positive(pressure)
        || !positive(alpha)
        || b.iter().any(|v| !v.is_finite())
    {
        return Err("cavity caloric operator requires finite bounded coordinates");
    }
    let delta = alpha * b.iter().map(|v| v * v).sum::<f64>();
    if !delta.is_finite() || delta <= 0.0 || delta >= 1.0 {
        return Err("cavity effective heat capacity is not positive");
    }
    let t = det::sqrt(1.0 - delta);
    let eta = alpha / (t * (1.0 + t));
    let scale = a.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if !positive(scale) || a.iter().any(|v| !v.is_finite()) {
        return Err("cavity diffusion rates are unrepresentable");
    }
    for value in &mut a {
        *value /= scale;
    }
    let ab: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * b[j]).sum())
        .collect();
    let bab: f64 = b.iter().zip(&ab).map(|(b, ab)| b * ab).sum();
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            h[i * n + j] =
                a[i * n + j] + eta * (ab[i] * b[j] + b[i] * ab[j]) + eta * eta * bab * b[i] * b[j];
        }
    }
    let norm = h.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if !positive(norm) || h.iter().any(|v| !v.is_finite()) {
        return Err("cavity capacity normalization is unrepresentable");
    }
    poll(cancelled)?;
    let (rates, vectors) = jacobi_eigh(&h, n);
    poll(cancelled)?;
    let mut arms = Vec::with_capacity(n);
    for k in 0..n {
        poll(cancelled)?;
        let mut projection = 0.0;
        let mut unit = 0.0;
        let mut residual = 0.0_f64;
        for i in 0..n {
            let v = vectors[i * n + k];
            projection += b[i] * v;
            unit += v * v;
            let hv: f64 = (0..n).map(|j| h[i * n + j] * vectors[j * n + k]).sum();
            residual = residual.max((hv - rates[k] * v).abs());
        }
        let rate = rates[k] * scale;
        let tau = rate.recip();
        let increment = pressure * alpha * (projection / t).powi(2);
        if !positive(rate)
            || !positive(tau)
            || !increment.is_finite()
            || increment < 0.0
            || (unit - 1.0).abs() > 2e-12
            || residual > 2e-12 * norm
        {
            return Err("cavity caloric spectrum is unresolved");
        }
        // Exactly unobservable coordinates have no compression residue.
        // Every positive residue survives, without an arbitrary cutoff.
        if increment > 0.0 {
            arms.push(BulkArm {
                increment_pa: increment,
                relaxation_time_s: tau,
            });
        }
    }
    let total_increment: f64 = arms.iter().map(|a| a.increment_pa).sum();
    let expected_increment = pressure * delta / (1.0 - delta);
    if !positive(total_increment)
        || !positive(expected_increment)
        || (total_increment - expected_increment).abs() > 2e-11 * expected_increment
    {
        return Err("cavity spectrum lost frozen compression storage");
    }
    Ok(Spectrum {
        arms,
        frozen_bulk_modulus_pa: pressure + total_increment,
    })
}
