//! fs-flowcert-e2e — FlowCert: a certified credibility map for a lattice-
//! Boltzmann channel flow. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! A CFD run gives you a number; it does not tell you whether to believe it.
//! This illuminates the operating space of a lattice-Boltzmann channel and
//! certifies credibility at every point, composing crates never designed to meet:
//!
//! - **The simulation + its exact solution** ([`fs_lbm`]): each channel is
//!   marched to STEADY STATE (not a fixed step budget — otherwise a slow
//!   high-Reynolds transient would masquerade as inaccuracy), then compared to
//!   the ANALYTIC Poiseuille solution: a manufactured-solution accuracy check
//!   that reflects the inherent `O(1/ny²)` discretization error.
//! - **The scaling certificate** ([`fs_lbm::plan_scaling`]): the lattice scaling
//!   planner derives `ν`, `τ`, and the Mach number for the target Reynolds and
//!   flags the regime `Verified` only when comfortably stable (positive viscosity,
//!   low Mach, safe `τ` margin) — the operating-envelope certificate.
//! - **Illumination** ([`fs_archive`]): MAP-Elites over (Reynolds × resolution)
//!   keeps the most-accurate operating point in every niche — the credibility
//!   atlas, not a single run.
//! - **Honest colors** ([`fs_evidence`]): once converged, EVERY point matches the
//!   analytic solution, so the credibility differentiation is the REGIME — a
//!   point accurate AND comfortably stable is `Verified`; a near-`τ=½` point is
//!   flagged `Estimated` as risky even where it is (currently) accurate.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_archive::MapElites;
use fs_evidence::{Color, ColorRank};
use fs_lbm::{Lbm, plan_scaling, poiseuille_analytic};

/// One certified operating point.
#[derive(Debug, Clone, Copy)]
pub struct OperatingPoint {
    /// Target Reynolds number.
    pub reynolds: f64,
    /// Channel resolution (rows).
    pub ny: usize,
    /// Relaxation time.
    pub tau: f64,
    /// Lattice viscosity.
    pub viscosity: f64,
    /// Relative max-profile error vs the analytic Poiseuille solution, at
    /// STEADY STATE (so this reflects inherent discretization error, not the
    /// step budget).
    pub profile_error: f64,
    /// Did the flow reach steady state within the step cap?
    pub converged: bool,
    /// Steps actually run to reach steady state.
    pub steps_run: usize,
    /// Is the (converged) profile accurate within tolerance?
    pub accurate: bool,
    /// Is the lattice scaling in a `Verified` (stable) regime?
    pub regime_stable: bool,
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct FlowReport {
    /// Every operating point probed.
    pub points: Vec<OperatingPoint>,
    /// Fraction of (Reynolds × resolution) niches filled.
    pub coverage: f64,
    /// QD score (Σ accuracy over the atlas).
    pub qd_score: f64,
    /// Number of filled niches.
    pub num_niches: usize,
    /// Best (smallest) profile error.
    pub best_error: f64,
    /// Are all points accurate within tolerance?
    pub all_accurate: bool,
    /// Fraction of points in a certified-stable regime.
    pub stable_fraction: f64,
    /// The credibility color: `Verified` iff every point is accurate & stable.
    pub credibility_color: Color,
}

/// Simulate one channel operating point (to steady state) and certify it.
/// `max_steps` caps the transient; the run stops early once the profile stops
/// changing, so `profile_error` measures the INHERENT accuracy rather than how
/// far the transient happened to decay in a fixed budget.
#[must_use]
pub fn certify_point(
    reynolds: f64,
    ny: usize,
    u_lattice: f64,
    max_steps: usize,
    tol: f64,
) -> OperatingPoint {
    let plan = plan_scaling(reynolds, ny as f64, u_lattice);
    let nu = plan.viscosity;
    let tau = plan.tau;
    // Body force sized so the centerline velocity ≈ u_lattice.
    let gx = 8.0 * nu * u_lattice / (ny as f64).powi(2);

    // March to steady state in fixed chunks; stop when the profile stabilizes.
    let mut lbm = Lbm::channel(4, ny, tau, gx);
    let chunk = 2000usize;
    let mut profile = lbm.x_velocity_profile();
    let mut steps_run = 0usize;
    let mut converged = false;
    while steps_run < max_steps {
        lbm.run(chunk);
        steps_run += chunk;
        let next = lbm.x_velocity_profile();
        let (mut delta, mut scale) = (0.0_f64, 1e-12_f64);
        for (a, b) in next.iter().zip(&profile) {
            delta = delta.max((a - b).abs());
            scale = scale.max(a.abs());
        }
        profile = next;
        // Steady once the per-chunk change is far below the O(1/ny²)
        // discretization floor — tighter would only burn steps.
        if delta / scale < 1e-4 {
            converged = true;
            break;
        }
    }

    let mut peak = 0.0_f64;
    let mut max_err = 0.0_f64;
    for (y, &u) in profile.iter().enumerate() {
        let exact = poiseuille_analytic(gx, nu, ny, y);
        peak = peak.max(exact.abs());
        max_err = max_err.max((u - exact).abs());
    }
    let profile_error = if peak > 1e-12 {
        max_err / peak
    } else {
        max_err
    };

    OperatingPoint {
        reynolds,
        ny,
        tau,
        viscosity: nu,
        profile_error,
        converged,
        steps_run,
        // Accuracy is only claimed for a converged, tolerance-matching profile.
        accurate: converged && profile_error <= tol,
        regime_stable: plan.color().rank() == ColorRank::Verified,
    }
}

/// Run the FlowCert campaign over the Reynolds × resolution grid.
///
/// # Panics
/// If `reynolds` or `resolutions` is empty.
#[must_use]
pub fn run_campaign(
    reynolds: &[f64],
    resolutions: &[usize],
    max_steps: usize,
    tol: f64,
) -> FlowReport {
    assert!(
        !reynolds.is_empty() && !resolutions.is_empty(),
        "empty sweep"
    );
    let re_lo = reynolds.iter().copied().fold(f64::INFINITY, f64::min);
    let re_hi = reynolds.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let ny_lo = *resolutions.iter().min().unwrap() as f64;
    let ny_hi = *resolutions.iter().max().unwrap() as f64;
    let mut archive = MapElites::new(
        vec![re_lo - 1.0, ny_lo - 1.0],
        vec![re_hi + 1.0, ny_hi + 1.0],
        vec![reynolds.len().max(2), resolutions.len().max(2)],
    );

    let mut points = Vec::new();
    for &re in reynolds {
        for &ny in resolutions {
            // A low lattice velocity keeps the Mach number inside the strict
            // low-Mach band the scaling planner certifies.
            let p = certify_point(re, ny, 0.05, max_steps, tol);
            // Fitness = accuracy in (0, 1]; the archive requires ≥ 0.
            let fitness = 1.0 / (1.0 + p.profile_error);
            archive.add(vec![re, ny as f64], vec![re, ny as f64], fitness);
            points.push(p);
        }
    }

    let all_accurate = points.iter().all(|p| p.accurate);
    let stable_count = points.iter().filter(|p| p.regime_stable).count();
    let stable_fraction = stable_count as f64 / points.len() as f64;
    let best_error = points
        .iter()
        .map(|p| p.profile_error)
        .fold(f64::INFINITY, f64::min);
    let credibility_color = if all_accurate && stable_count == points.len() {
        // declared-color-ok: demo credibility candidate derived from local profile-error bounds; admitted only at a consumer's authority boundary (6pf9)
        Color::Verified {
            lo: 0.0,
            hi: points.iter().map(|p| p.profile_error).fold(0.0, f64::max),
        }
    } else {
        Color::Estimated {
            estimator: "lbm-credibility".to_string(),
            dispersion: best_error,
        }
    };

    FlowReport {
        coverage: archive.coverage(),
        qd_score: archive.qd_score(),
        num_niches: archive.num_elites(),
        best_error,
        all_accurate,
        stable_fraction,
        credibility_color,
        points,
    }
}

/// The default sweep.
#[must_use]
pub fn default_sweep() -> (Vec<f64>, Vec<usize>) {
    // Reynolds spans a comfortable regime (τ well above ½) to a near-boundary
    // one (τ→½); resolutions kept modest so every point reaches steady state fast.
    (vec![20.0, 50.0, 90.0], vec![12, 16, 24])
}
