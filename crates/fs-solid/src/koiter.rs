//! Koiter-style initial post-buckling asymptotics [F] (feature
//! `koiter-asymptotics`, off by default per the Ambition-Tag gating
//! rule): the a/b coefficients that classify a bifurcation as
//! imperfection-tolerant or imperfection-sensitive.
//!
//! v1 computes the coefficients by finite-difference expansion of the
//! exact potential energy along the (normalized) buckling mode at the
//! critical state — `Π(t) = Π₀ + ½Π″t² + (a)t³ + (b)t⁴ + …` with the
//! quadratic term stationary at λ_cr — and classifies:
//! - `|a| > tol` → ASYMMETRIC (imperfection-sensitive, λ_s knockdown
//!   O(ξ^1/2));
//! - `a ≈ 0, b > 0` → SYMMETRIC STABLE (imperfection-tolerant);
//! - `a ≈ 0, b < 0` → SYMMETRIC UNSTABLE (knockdown O(ξ^2/3)).
//!
//! The SAMPLED-CONTINUATION FALLBACK ORACLE (the bead's blessed
//! fallback) cross-checks the classification: imperfect paths of a
//! symmetric-stable structure show no limit point below the critical
//! load; sensitive ones snap early. The battery runs both.

use crate::hyper2d::HyperProblem;

/// Bifurcation classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bifurcation {
    /// a ≈ 0, b > 0: rising symmetric post-buckling path.
    SymmetricStable,
    /// a ≈ 0, b < 0: falling symmetric path (imperfection-sensitive).
    SymmetricUnstable,
    /// a ≠ 0: one-sided sensitivity.
    Asymmetric,
}

/// Initial post-buckling coefficients from the energy expansion.
#[derive(Debug, Clone, Copy)]
pub struct KoiterCoefficients {
    /// Cubic coefficient (energy expansion / 3!).
    pub a: f64,
    /// Quartic coefficient (energy expansion / 4!).
    pub b: f64,
    /// The classification at tolerance `a_tol`.
    pub class: Bifurcation,
}

/// FD energy expansion along the buckling mode at the critical state.
/// `scale` sets the FD amplitude (mode is max-normalized internally);
/// `a_tol` is the asymmetry gate relative to |b|·scale.
#[must_use]
pub fn koiter_coefficients(
    problem: &HyperProblem<'_>,
    u_cr: &[f64],
    lambda_cr: f64,
    mode: &[[f64; 2]],
    scale: f64,
    a_tol: f64,
) -> Option<KoiterCoefficients> {
    let n = problem.mesh.node_count();
    let mmax = mode
        .iter()
        .fold(0.0f64, |m, v| m.max(v[0].abs()).max(v[1].abs()));
    if mmax <= 0.0 {
        return None;
    }
    let pi_at = |t: f64| -> Option<f64> {
        let mut u = u_cr.to_vec();
        for node in 0..n {
            u[2 * node] += t * mode[node][0] / mmax;
            u[2 * node + 1] += t * mode[node][1] / mmax;
        }
        problem.potential_energy(&u, lambda_cr)
    };
    let h = scale;
    let (p2m, p1m, p0, p1p, p2p) = (
        pi_at(-2.0 * h)?,
        pi_at(-h)?,
        pi_at(0.0)?,
        pi_at(h)?,
        pi_at(2.0 * h)?,
    );
    let d3 = (p2p - 2.0 * p1p + 2.0 * p1m - p2m) / (2.0 * h * h * h);
    let d4 = (p2p - 4.0 * p1p + 6.0 * p0 - 4.0 * p1m + p2m) / (h * h * h * h);
    let a = d3 / 6.0;
    let b = d4 / 24.0;
    let class = if a.abs() > a_tol * b.abs().max(1e-30) {
        Bifurcation::Asymmetric
    } else if b > 0.0 {
        Bifurcation::SymmetricStable
    } else {
        Bifurcation::SymmetricUnstable
    };
    Some(KoiterCoefficients { a, b, class })
}
