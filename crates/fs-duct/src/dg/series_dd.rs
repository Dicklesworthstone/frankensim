//! Double-double fallback for the dense series propagator.
//!
//! `DgDuct::build_series_propagator` certifies its f64 matrices by `max(carried, halving)`,
//! where `halving` is the Frobenius distance between the series taken with `s` and with `s + 1`
//! squarings. On a stiff chimney (wall memory out to 20 kHz at a 96 kHz clock, 13 squarings)
//! the Taylor remainder is far below the tolerance (`carried` 7.5e-15) but `halving` is a
//! rounding estimate that doubles with every squaring (6.07e-12 at 12 squarings, 1.136e-11 at 13
//! against a tolerance of 1e-11), so the build refuses a duct whose f64 entries are good to a
//! few 1e-12.
//!
//! This module runs ONLY where that certificate has refused with a finite value, for ducts of at
//! most [`SERIES_DD_MAX_DOFS`] DOFs. It evaluates the same series in double-double (about 106
//! bits) and re-certifies it with the same test, never loosened:
//!
//! - `A = L dt / 2^s` is formed EXACTLY in double-double (the product of two f64 is exact as a
//!   pair), so the result is the exponential of the operator as given, not of its f64-rounded
//!   scaling;
//! - the Taylor series runs to a scaled 1-norm remainder of 1e-32 (the f64 path stops at 1e-18),
//!   and `carried = 2^s sqrt(n) remainder` is formed as before;
//! - `halving` is the Frobenius distance between the double-double `phi_0` at `s` and at
//!   `s + 1` squarings, differenced and summed in double-double;
//! - the result is rounded to f64, and that rounding is counted: `rounding` is one ulp per
//!   entry (`2^-52` times the Frobenius norm, of the largest of `phi_0`, `phi_1`, `phi_2`),
//!   added to the certificate. The certificate is `max(carried, halving) + rounding`.
//!
//! A duct that passes the f64 certificate never reaches this module, so its bits are the f64
//! path's bits. A duct this module cannot certify returns the f64 path's own refusal.

use fs_math::det;
use fs_math::dd::Dd;

/// Largest dense duct (total DOFs) the fallback will take. A double-double product costs about
/// 35 f64 operations, and the series needs about 50 dense products of an `n x n` matrix; 512
/// DOFs is minutes, 2048 (the dense cap) would be hours. Above this the f64 refusal stands.
pub(super) const SERIES_DD_MAX_DOFS: usize = 512;

/// Scaled 1-norm Taylor remainder at which the double-double series stops.
const DD_REMAINDER: f64 = 1e-32;

/// The double-double series rounded to f64, with the three terms of its certificate.
#[derive(Debug, Clone)]
pub(super) struct DdSeries {
    pub(super) phi0: Vec<f64>,
    pub(super) phi1: Vec<f64>,
    pub(super) phi2: Vec<f64>,
    /// `2^s sqrt(n) *` the double-double Taylor remainder bound.
    pub(super) carried: f64,
    /// Frobenius distance of the double-double `phi_0` at `s` and `s + 1` squarings.
    pub(super) halving: f64,
    /// One ulp per entry of the rounding to f64.
    pub(super) rounding: f64,
}

impl DdSeries {
    /// `max(carried, halving) + rounding`: the f64 certificate's own form, plus the rounding.
    pub(super) fn certificate(&self) -> f64 {
        self.carried.max(self.halving) + self.rounding
    }
}

/// Whether the fallback runs: the duct is small enough, and the f64 certificate refused with a
/// finite value (a NaN or infinite certificate is a bad operator, not a roundoff floor).
#[allow(clippy::float_cmp)]
pub(super) fn engages(n: usize, f64_certificate: f64, tol: f64) -> bool {
    n <= SERIES_DD_MAX_DOFS && f64_certificate.is_finite() && f64_certificate > tol
}

#[cfg(test)]
thread_local! {
    /// Refinements run on this thread (the identity test reads it to show the fallback is not
    /// reached on a duct that passes today).
    pub(super) static REFINEMENTS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

fn identity(n: usize) -> Vec<Dd> {
    let mut out = vec![Dd::ZERO; n * n];
    for i in 0..n {
        out[i * n + i] = Dd::ONE;
    }
    out
}

/// `a b`, skipping zero entries of `a` (put the sparse factor on the left).
#[allow(clippy::float_cmp)]
fn mat_mul(a: &[Dd], b: &[Dd], n: usize) -> Vec<Dd> {
    let mut c = vec![Dd::ZERO; n * n];
    for i in 0..n {
        let row = &mut c[i * n..(i + 1) * n];
        for k in 0..n {
            let aik = a[i * n + k];
            if aik.hi == 0.0 {
                continue;
            }
            let b_row = &b[k * n..(k + 1) * n];
            for (cj, bj) in row.iter_mut().zip(b_row) {
                *cj = *cj + aik * *bj;
            }
        }
    }
    c
}

fn halve(x: Dd) -> Dd {
    Dd { hi: 0.5 * x.hi, lo: 0.5 * x.lo }
}

fn double(x: Dd) -> Dd {
    Dd { hi: 2.0 * x.hi, lo: 2.0 * x.lo }
}

struct Phis {
    phi0: Vec<Dd>,
    phi1: Vec<Dd>,
    phi2: Vec<Dd>,
    remainder: f64,
}

/// `phi_0(L dt)` (and with `all`, `phi_1`, `phi_2`) by the Taylor series of `L dt / 2^shift`
/// and the exact doubling recurrences, in double-double. `phi_1` and `phi_2` are left empty
/// without `all`: the `phi_0` recurrence does not use them.
fn phis(l: &[f64], n: usize, dt: f64, shift: u32, all: bool) -> Phis {
    let dt_sub = dt / f64::from(1u32 << shift);
    let a: Vec<Dd> = l
        .iter()
        .map(|&v| Dd::from_f64(v) * Dd::from_f64(dt_sub))
        .collect();
    let mut a_norm = 0.0f64;
    for j in 0..n {
        a_norm = a_norm.max((0..n).map(|i| a[i * n + j].hi.abs()).sum());
    }
    // The hi parts carry the norm to a few 1e-14 (n terms of 2^-53); a margin covers it.
    let a_norm = a_norm * (1.0 + 1e-9);
    let mut phi0 = identity(n);
    let (mut phi1, mut phi2): (Vec<Dd>, Vec<Dd>) = if all {
        (identity(n), identity(n).into_iter().map(halve).collect())
    } else {
        (Vec::new(), Vec::new())
    };
    let mut term = identity(n);
    let mut remainder = f64::INFINITY;
    let mut power = 1.0;
    let mut factorial = 1.0;
    for k in 1..=64usize {
        // A is the sparse factor: A * (A^{k-1}/(k-1)!) = A^k/(k-1)!, divided by k below.
        term = mat_mul(&a, &term, n);
        let inv_k = Dd::ONE / Dd::from_f64(k as f64);
        let inv_k1 = Dd::ONE / Dd::from_f64((k + 1) as f64);
        let inv_k2 = Dd::ONE / Dd::from_f64(((k + 1) * (k + 2)) as f64);
        for i in 0..(n * n) {
            term[i] = term[i] * inv_k;
            phi0[i] = phi0[i] + term[i];
            if all {
                phi1[i] = phi1[i] + term[i] * inv_k1;
                phi2[i] = phi2[i] + term[i] * inv_k2;
            }
        }
        power *= a_norm;
        factorial *= k as f64;
        let next = power * a_norm / (factorial * (k + 1) as f64);
        remainder = next / (1.0 - a_norm / (k + 2) as f64);
        if k >= 12 && remainder <= DD_REMAINDER {
            break;
        }
    }
    for _ in 0..shift {
        let next_phi0 = mat_mul(&phi0, &phi0, n);
        if all {
            let mut i_plus_phi0 = phi0.clone();
            for d in 0..n {
                i_plus_phi0[d * n + d] = i_plus_phi0[d * n + d] + Dd::ONE;
            }
            let next_phi1: Vec<Dd> = mat_mul(&i_plus_phi0, &phi1, n)
                .into_iter()
                .map(halve)
                .collect();
            let phi1_sq = mat_mul(&phi1, &phi1, n);
            let next_phi2: Vec<Dd> = phi1_sq
                .iter()
                .zip(&phi2)
                .map(|(a, b)| halve(halve(*a + double(*b))))
                .collect();
            phi1 = next_phi1;
            phi2 = next_phi2;
        }
        phi0 = next_phi0;
    }
    Phis { phi0, phi1, phi2, remainder }
}

fn frobenius(a: &[f64]) -> f64 {
    det::sqrt(a.iter().map(|x| x * x).sum::<f64>())
}

/// The double-double series of the dense operator `l` (`n x n`, row-major) at step `dt` with
/// `shift` squarings, rounded to f64, with its certificate terms. `None` when `shift + 1` would
/// pass the 30 halvings the f64 path allows, or any term is not finite.
pub(super) fn refine(l: &[f64], n: usize, dt: f64, shift: u32) -> Option<DdSeries> {
    if shift + 1 > 30 {
        return None;
    }
    #[cfg(test)]
    REFINEMENTS.with(|count| count.set(count.get() + 1));
    let fine = phis(l, n, dt, shift, true);
    let finer = phis(l, n, dt, shift + 1, false);
    let mut sum = Dd::ZERO;
    for (a, b) in fine.phi0.iter().zip(&finer.phi0) {
        let d = *a - *b;
        sum = sum + d * d;
    }
    let halving = sum.sqrt().to_f64();
    let carried = f64::from(1u32 << shift) * det::sqrt(n as f64) * fine.remainder;
    let round = |m: &[Dd]| -> Vec<f64> { m.iter().map(|v| v.to_f64()).collect() };
    let (phi0, phi1, phi2) = (round(&fine.phi0), round(&fine.phi1), round(&fine.phi2));
    let rounding =
        f64::EPSILON * frobenius(&phi0).max(frobenius(&phi1)).max(frobenius(&phi2));
    if !(carried.is_finite() && halving.is_finite() && rounding.is_finite()) {
        return None;
    }
    Some(DdSeries { phi0, phi1, phi2, carried, halving, rounding })
}

#[cfg(test)]
mod tests;
