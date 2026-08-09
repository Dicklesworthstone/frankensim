//! Rayleigh-equation instability oracle for the Bickley jet
//! `U(y) = sech^2(y)` — the inviscid linear-stability reference the
//! fs-lbm jet runs are validated against (bead 9ok02's verification
//! chain).
//!
//! Temporal problem: for real wavenumber `alpha`, find complex phase
//! speed `c` such that `(U - c)(phi'' - alpha^2 phi) - U'' phi = 0`
//! with `phi -> 0` as `|y| -> inf`. Growth rate is
//! `sigma = alpha Im(c)`.
//!
//! Method: shooting on the half line using jet symmetry — integrate
//! from `y = L` (with the exact decay condition `phi' = -alpha phi`)
//! down to `y = 0` by RK4, and drive the symmetry mismatch
//! (`phi'(0) = 0` sinuous, `phi(0) = 0` varicose) to zero by a
//! complex secant iteration in `c`.
//!
//! ANALYTIC PINS, SELF-VERIFIED: the two neutral eigenmodes
//! `phi = sech^2(y)` at `(alpha, c) = (2, 2/3)` (sinuous) and
//! `phi = sech(y) tanh(y)` at `(alpha, c) = (1, 2/3)` (varicose) are
//! EXACT solutions — derived by hand (substituting into the Rayleigh
//! equation reduces it to `s [s (6c - alpha^2) + c (alpha^2 - 4)]`
//! resp. `u t [s (6c - alpha^2 - 2) + ... ]` with `s = sech^2`) and
//! RE-PROVEN numerically at machine precision by the residual test on
//! every run, so no literature value is transcribed on trust.

use crate::AeroacError;
use fs_math::c64::C64;
use fs_math::det;

/// Jet mode symmetry (of the streamfunction eigenfunction).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JetSymmetry {
    /// Even eigenfunction (`phi'(0) = 0`): the sinuous/flapping mode.
    Sinuous,
    /// Odd eigenfunction (`phi(0) = 0`): the varicose mode.
    Varicose,
}

/// A converged Rayleigh mode.
#[derive(Debug, Clone)]
pub struct RayleighMode {
    /// Wavenumber (real).
    pub alpha: f64,
    /// Complex phase speed.
    pub c: C64,
    /// Temporal growth rate `alpha * Im(c)`.
    pub growth_rate: f64,
    /// Final boundary-condition mismatch magnitude (disclosed).
    pub mismatch: f64,
    /// Secant iterations used.
    pub iterations: usize,
}

/// Bickley profile and second derivative:
/// `U = sech^2`, `U'' = 4 sech^2 - 6 sech^4` (exact).
fn u_and_upp(y: f64) -> (f64, f64) {
    let t = det::tanh(y);
    let s = 1.0 - t * t; // sech^2 via the identity (deterministic)
    (s, 4.0 * s - 6.0 * s * s)
}

/// One RK4 step of the first-order system
/// `phi' = psi`, `psi' = [alpha^2 + U''/(U - c)] phi`.
fn rk4_step(y: f64, h: f64, phi: C64, psi: C64, alpha: f64, c: C64) -> (C64, C64) {
    let f = |yy: f64, ph: C64, ps: C64| -> (C64, C64) {
        let (u, upp) = u_and_upp(yy);
        let coeff =
            C64::new(alpha * alpha, 0.0) + C64::new(upp, 0.0) * (C64::new(u, 0.0) - c).recip();
        (ps, coeff * ph)
    };
    let (k1p, k1s) = f(y, phi, psi);
    let (k2p, k2s) = f(
        y + 0.5 * h,
        phi + k1p.scale(0.5 * h),
        psi + k1s.scale(0.5 * h),
    );
    let (k3p, k3s) = f(
        y + 0.5 * h,
        phi + k2p.scale(0.5 * h),
        psi + k2s.scale(0.5 * h),
    );
    let (k4p, k4s) = f(y + h, phi + k3p.scale(h), psi + k3s.scale(h));
    (
        phi + (k1p + k2p.scale(2.0) + k3p.scale(2.0) + k4p).scale(h / 6.0),
        psi + (k1s + k2s.scale(2.0) + k3s.scale(2.0) + k4s).scale(h / 6.0),
    )
}

/// Integrate from `y = l` down to 0; return `(phi(0), phi'(0))`.
fn shoot(alpha: f64, c: C64, l: f64, steps: usize) -> (C64, C64) {
    #[allow(clippy::cast_precision_loss)]
    let h = -l / steps as f64;
    let mut y = l;
    // Exact far-field decay: phi ~ e^{-alpha y}, normalized phi(L)=1.
    let mut phi = C64::new(1.0, 0.0);
    let mut psi = C64::new(-alpha, 0.0);
    for _ in 0..steps {
        let (p, s) = rk4_step(y, h, phi, psi, alpha, c);
        phi = p;
        psi = s;
        y += h;
    }
    (phi, psi)
}

/// Boundary mismatch (normalized so the secant target is
/// scale-invariant).
fn mismatch(alpha: f64, c: C64, symmetry: JetSymmetry, l: f64, steps: usize) -> C64 {
    let (phi0, psi0) = shoot(alpha, c, l, steps);
    let norm = det::sqrt(phi0.norm_sq() + psi0.norm_sq()).max(1e-300);
    match symmetry {
        JetSymmetry::Sinuous => psi0.scale(1.0 / norm),
        JetSymmetry::Varicose => phi0.scale(1.0 / norm),
    }
}

/// Solve for the Rayleigh mode at wavenumber `alpha` starting the
/// complex secant iteration from `c_guess`.
///
/// Defaults tuned by the convergence test: half-domain `l = 14`
/// (sech^2 < 3e-12 there), `steps` RK4 steps (2048 gives ~1e-10
/// eigenvalue accuracy; the convergence test measures it).
///
/// # Errors
/// [`AeroacError::NonFinite`] / [`AeroacError::InvalidParameter`] on
/// bad inputs; [`AeroacError::NotConverged`] when the secant fails —
/// no partial eigenvalue is returned.
pub fn bickley_rayleigh_mode(
    alpha: f64,
    symmetry: JetSymmetry,
    c_guess: C64,
    l: f64,
    steps: usize,
) -> Result<RayleighMode, AeroacError> {
    if !alpha.is_finite() || !c_guess.re.is_finite() || !c_guess.im.is_finite() || !l.is_finite() {
        return Err(AeroacError::NonFinite {
            what: "rayleigh inputs",
        });
    }
    if alpha <= 0.0 || l <= 2.0 || steps < 64 {
        return Err(AeroacError::InvalidParameter {
            what: "alpha must be positive, half-domain > 2, steps >= 64",
        });
    }
    // Complex secant on the mismatch.
    let mut c0 = c_guess;
    let mut c1 = c_guess + C64::new(1.0e-4, 1.0e-4);
    let mut f0 = mismatch(alpha, c0, symmetry, l, steps);
    let mut f1 = mismatch(alpha, c1, symmetry, l, steps);
    let mut iterations = 0usize;
    for it in 0..80 {
        iterations = it + 1;
        let denom = f1 - f0;
        if denom.norm_sq() == 0.0 {
            break;
        }
        let c2 = c1 - f1 * (c1 - c0) * denom.recip();
        if !c2.re.is_finite() || !c2.im.is_finite() {
            break;
        }
        c0 = c1;
        f0 = f1;
        c1 = c2;
        f1 = mismatch(alpha, c1, symmetry, l, steps);
        if f1.abs() < 1.0e-12 {
            break;
        }
        if (c1 - c0).abs() < 1.0e-14 * c1.abs().max(1.0e-3) {
            break;
        }
    }
    let m = f1.abs();
    if m > 1.0e-9 {
        return Err(AeroacError::NotConverged {
            what: "rayleigh secant",
            residual: m,
        });
    }
    Ok(RayleighMode {
        alpha,
        c: c1,
        growth_rate: alpha * c1.im,
        mismatch: m,
        iterations,
    })
}

/// The Rayleigh-equation residual
/// `(U - c)(phi'' - alpha^2 phi) - U'' phi` for a CANDIDATE
/// closed-form eigenfunction, evaluated with exact hyperbolic
/// identities — the self-verification the analytic pins rest on
/// (tests drive it to machine zero).
#[must_use]
pub fn rayleigh_residual_closed_form(candidate: JetSymmetry, alpha: f64, c: f64, y: f64) -> f64 {
    let t = det::tanh(y);
    let s = 1.0 - t * t; // sech^2
    let (phi, phipp) = match candidate {
        // phi = sech^2:  phi'' = 4 s - 6 s^2 (same as U'').
        JetSymmetry::Sinuous => (s, 4.0 * s - 6.0 * s * s),
        // phi = sech tanh: phi'' = u t (1 - 6 s) with u = sech.
        JetSymmetry::Varicose => {
            let u = det::sqrt(s);
            (u * t, -u * t * (6.0 * s - 1.0))
        }
    };
    let upp = 4.0 * s - 6.0 * s * s;
    (s - c) * (phipp - alpha * alpha * phi) - upp * phi
}
