//! Stage 4: ROBUSTIFICATION. The stability objective is cheap (an
//! Orr–Sommerfeld solve per station), so the fluid band — Carreau
//! parameters spanning the target liquids × pour rates — is swept
//! directly: the NOMINAL design minimizes growth at the band center,
//! the ROBUST design minimizes the CVaR of growth over the band
//! (Rockafellar–Uryasev empirical form, the fs-frame pattern). The
//! battery gates the flagship claim: the robust lip beats the nominal
//! lip on off-nominal fluids. Candidate screening under the (noisy)
//! validator runs through fs-race — dominated lips die early with
//! anytime validity, eliminations ledgered.

use crate::stability::{VesselProfile, growth_objective};

/// The robustification record.
#[derive(Debug, Clone)]
pub struct RobustReport {
    /// Lip width chosen by nominal-only optimization.
    pub nominal_lip: f64,
    /// Lip width chosen by CVaR optimization over the band.
    pub robust_lip: f64,
    /// Worst off-nominal growth of the nominal design.
    pub nominal_offband_growth: f64,
    /// Worst off-nominal growth of the robust design.
    pub robust_offband_growth: f64,
    /// CVaR level used.
    pub beta: f64,
}

/// Empirical CVaR_β (Rockafellar–Uryasev; exact for empirical
/// measures: the tail mean past the β-quantile).
#[must_use]
pub fn empirical_cvar(losses: &[f64], beta: f64) -> f64 {
    let mut sorted = losses.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    #[allow(clippy::cast_precision_loss)]
    let q = ((beta * n as f64).floor() as usize).min(n - 1);
    let t = sorted[q];
    let tail: f64 = sorted.iter().map(|&l| (l - t).max(0.0)).sum();
    #[allow(clippy::cast_precision_loss)]
    {
        t + tail / ((1.0 - beta) * n as f64)
    }
}

/// The (rate, viscosity) band: nominal center plus off-nominal
/// corners (the family of liquids the spout must serve).
#[must_use]
pub fn fluid_band() -> Vec<(f64, f64)> {
    let rates = [0.8f64, 1.0, 1.3];
    let viscs = [0.6f64, 1.0, 1.8];
    let mut band = Vec::new();
    for &r in &rates {
        for &v in &viscs {
            band.push((r, v));
        }
    }
    band
}

/// Optimize the lip width two ways (nominal vs CVaR over the band) by
/// a deterministic golden-section-free grid refinement (the objective
/// is cheap and smooth in the scalar knob), then evaluate both on the
/// OFF-NOMINAL corners.
#[must_use]
pub fn robustify(beta: f64) -> RobustReport {
    let stations = 4;
    let modes = 4;
    let nominal_obj = |lip: f64| -> f64 {
        let p = VesselProfile::carafe(lip);
        growth_objective(&p, 1.0, 1.0, stations, modes)
    };
    let cvar_obj = |lip: f64| -> f64 {
        let p = VesselProfile::carafe(lip);
        let losses: Vec<f64> = fluid_band()
            .iter()
            .map(|&(r, v)| growth_objective(&p, r, v, stations, modes))
            .collect();
        empirical_cvar(&losses, beta)
    };
    let minimize = |f: &dyn Fn(f64) -> f64| -> f64 {
        // Two-stage deterministic grid refinement over the lip range.
        let (mut lo, mut hi) = (0.5f64, 3.0f64);
        let mut best = lo;
        for _ in 0..2 {
            let mut best_v = f64::INFINITY;
            for k in 0..=12 {
                let lip = lo + (hi - lo) * f64::from(k) / 12.0;
                let v = f(lip);
                if v < best_v {
                    best_v = v;
                    best = lip;
                }
            }
            let third = (hi - lo) / 3.0;
            lo = (best - third).max(0.5);
            hi = (best + third).min(3.0);
        }
        best
    };
    let nominal_lip = minimize(&nominal_obj);
    let robust_lip = minimize(&cvar_obj);
    // Off-nominal evaluation: worst growth over the band corners
    // EXCLUDING the nominal center.
    let offband = |lip: f64| -> f64 {
        let p = VesselProfile::carafe(lip);
        fluid_band()
            .iter()
            .filter(|&&(r, v)| (r - 1.0).abs() > 1e-12 || (v - 1.0).abs() > 1e-12)
            .map(|&(r, v)| growth_objective(&p, r, v, stations, modes))
            .fold(f64::NEG_INFINITY, f64::max)
    };
    RobustReport {
        nominal_lip,
        robust_lip,
        nominal_offband_growth: offband(nominal_lip),
        robust_offband_growth: offband(robust_lip),
        beta,
    }
}
