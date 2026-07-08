//! Stage 5: CVaR-constrained MASS MINIMIZATION in the
//! Rockafellar–Uryasev form: CVaR_β(L) = min_t t + E[(L−t)₊]/(1−β),
//! evaluated over the motion ensemble with the section scale s as the
//! design variable. Peak drift decreases monotonically in s at smoke
//! scale (bigger sections, stiffer/stronger hinges), so the minimal
//! feasible scale is found by bisection — deterministic and honest;
//! the multi-variable trust-region tier is the recorded successor.
//! The chosen scale then snaps UP to the section catalog and the
//! snapped design is INDEPENDENTLY re-checked.

use crate::history::{StoryFrame, StoryParams, peak_drift};
use fs_scenario::ensemble::StochasticEnsemble;

/// The CVaR design record.
pub struct CvarDesign {
    /// Minimal feasible section scale (continuous).
    pub scale_star: f64,
    /// Catalog-snapped scale (≥ scale_star).
    pub scale_snapped: f64,
    /// CVaR at the snapped design (re-checked, must pass).
    pub cvar_snapped: f64,
    /// CVaR at the continuous optimum.
    pub cvar_star: f64,
    /// Mass proxy at the snapped design (scale × member count — the
    /// smoke-tier stand-in for Σ ρAL).
    pub mass: f64,
    /// Bisection iterations.
    pub iters: u32,
}

/// Empirical CVaR_β by the Rockafellar–Uryasev minimization (exact
/// for empirical measures: t* is the β-quantile).
#[must_use]
pub fn empirical_cvar(losses: &[f64], beta: f64) -> f64 {
    let mut sorted = losses.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    let q_idx = ((beta * n as f64).floor() as usize).min(n - 1);
    let t = sorted[q_idx];
    let tail: f64 = sorted.iter().map(|&l| (l - t).max(0.0)).sum();
    t + tail / ((1.0 - beta) * n as f64)
}

/// CVaR of the peak-drift loss over the ensemble at section scale
/// `s` — the battery's monotonicity probe and limit-bracketing tool.
#[must_use]
pub fn ensemble_cvar(ensemble: &StochasticEnsemble, base: StoryParams, s: f64, beta: f64) -> f64 {
    empirical_cvar(&losses(ensemble, base, s), beta)
}

/// Peak-drift losses over the whole ensemble at section scale `s`.
fn losses(ensemble: &StochasticEnsemble, base: StoryParams, s: f64) -> Vec<f64> {
    let dt = ensemble.dt.value;
    let mut out = Vec::with_capacity(ensemble.members as usize);
    for member in 0..ensemble.members {
        let real = ensemble.realize(member).expect("ensemble realizes");
        let params = StoryParams { scale: s, ..base };
        let mut frame = StoryFrame::new(params);
        let drifts = frame.run(&real.values, dt);
        out.push(peak_drift(&drifts, base.h));
    }
    out
}

/// Minimize mass (∝ scale) subject to CVaR_β(peak drift) ≤ `limit` by
/// bisection on the scale, then snap UP to `catalog` and re-check.
///
/// # Panics
/// If even the largest catalog scale is infeasible (the drill fixture
/// checks the diagnostics path instead).
#[must_use]
pub fn cvar_mass_min(
    ensemble: &StochasticEnsemble,
    base: StoryParams,
    beta: f64,
    limit: f64,
    catalog: &[f64],
) -> CvarDesign {
    let cvar_at = |s: f64| empirical_cvar(&losses(ensemble, base, s), beta);
    let (mut lo, mut hi) = (0.25f64, 4.0f64);
    assert!(
        cvar_at(hi) <= limit,
        "even scale {hi} violates the CVaR limit — infeasible study"
    );
    let mut iters = 0u32;
    // If the smallest scale is already feasible, take it.
    if cvar_at(lo) <= limit {
        hi = lo;
    }
    while hi - lo > 0.02 {
        let mid = f64::midpoint(lo, hi);
        if cvar_at(mid) <= limit {
            hi = mid;
        } else {
            lo = mid;
        }
        iters += 1;
    }
    let scale_star = hi;
    let cvar_star = cvar_at(scale_star);
    // Snap UP to the catalog (feasibility preserved by monotonicity —
    // and re-checked anyway).
    let scale_snapped = catalog
        .iter()
        .copied()
        .filter(|&c| c >= scale_star)
        .fold(f64::INFINITY, f64::min);
    assert!(
        scale_snapped.is_finite(),
        "catalog has no section above the optimum — infeasible snap"
    );
    let cvar_snapped = cvar_at(scale_snapped);
    CvarDesign {
        scale_star,
        scale_snapped,
        cvar_snapped,
        cvar_star,
        mass: scale_snapped * 2.0,
        iters,
    }
}
