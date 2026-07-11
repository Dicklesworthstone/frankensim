//! fs-adaptbo-e2e — AnytimeBO: Bayesian optimization with anytime-valid
//! evidence for a declared stall-rate null. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! Every practical Bayesian optimizer faces "have we searched enough?" — and
//! peeking at the best-so-far after each iteration and stopping on a threshold
//! silently inflates the chance of stopping too early. This answers the question
//! with a guarantee, composing crates never designed to meet:
//!
//! - **The optimizer** ([`fs_bo`]): a Matérn-5⁄2 Gaussian process with closed-form
//!   Expected Improvement drives a deterministic minimization loop over a
//!   candidate grid.
//! - **The stopping certificate** ([`fs_eproc`]): a betting e-process watches a
//!   per-iteration STALL indicator (did the best improve by at least `δ`?). Its
//!   null is "we are still improving often enough"; when the e-value crosses
//!   `1/α` the search stops — an ANYTIME-VALID decision (Ville's inequality), so
//!   testing after every iteration never inflates the false-stop rate beyond `α`.
//! - **A trace diagnostic** ([`fs_eproc::GaussianMixtureCs`]): the mixture
//!   boundary is evaluated on the adaptive best-value trace to visualize its
//!   contraction. No fixed-mean coverage claim is made for that trace.
//! - **Honest evidence** ([`fs_evidence`]): the stop carries an e-value candidate
//!   for its declared statistical null; the observed incumbent and GP remain
//!   `Estimated` and do not masquerade as global-optimum enclosures.
//!
//! Deterministic (a fixed grid, a polynomial objective — no RNG, no libm). No
//! dependencies beyond the composed crates.

use fs_bo::{Gp, Kernel, Matern, expected_improvement};
use fs_eproc::{BettingEProcess, GaussianMixtureCs};
use fs_evidence::{Color, StatisticalCertificate};

const MAX_CAMPAIGN_ITERS: usize = 64;

/// The black-box objective: a tilted double well on `[0, 4]`, minimized near
/// `x ≈ 3` (polynomial ⇒ bit-identical on every ISA).
#[must_use]
pub fn objective(x: f64) -> f64 {
    (x - 1.0).powi(2) * (x - 3.0).powi(2) - 0.15 * x
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct AdaptBoReport {
    /// The best design found.
    pub best_x: f64,
    /// Its objective value.
    pub best_value: f64,
    /// BO iterations run.
    pub iterations: usize,
    /// Total objective evaluations (initial design + iterations).
    pub evaluations: usize,
    /// Did the e-process stop the search before the iteration cap?
    pub stopped_early: bool,
    /// The final e-process log e-value (evidence the search has stalled).
    pub log_e_value: f64,
    /// Center of the anytime-valid confidence sequence on the best-value TRACE
    /// (a running-mean diagnostic, not a bound on the true optimum).
    pub ci_center: f64,
    /// Its (shrinking) radius.
    pub ci_radius: f64,
    /// Statistical evidence against the declared stall-rate null. This is a
    /// certificate candidate, not an admitted numerical optimum enclosure.
    pub stop_evidence: StatisticalCertificate,
    /// Epistemic color of the observed incumbent value.
    pub incumbent_color: Color,
    /// The surrogate color (`Estimated` — a GP).
    pub surrogate_color: Color,
}

fn argmin(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let mut bi = 0usize;
    for i in 1..ys.len() {
        if ys[i] < ys[bi] {
            bi = i;
        }
    }
    (xs[bi], ys[bi])
}

/// Run the AnytimeBO campaign. Stops when the e-process rejects the declared
/// stall-rate null at `alpha`, or after `max_iters` iterations.
///
/// # Panics
/// If `max_iters` exceeds the bounded native work envelope, `delta` is not
/// finite and non-negative, or `alpha` is outside `(0, 1)`.
#[must_use]
pub fn run_campaign(max_iters: usize, delta: f64, alpha: f64) -> AdaptBoReport {
    assert!(
        max_iters <= MAX_CAMPAIGN_ITERS,
        "max_iters exceeds the bounded campaign envelope ({MAX_CAMPAIGN_ITERS})"
    );
    assert!(
        delta.is_finite() && delta >= 0.0,
        "improvement threshold must be finite and non-negative"
    );
    assert!(
        alpha.is_finite() && alpha > 0.0 && alpha < 1.0,
        "alpha must be finite and lie in (0,1)"
    );
    // Candidate grid over [0, 4].
    let grid: Vec<f64> = (0..=80).map(|i| 4.0 * f64::from(i) / 80.0).collect();
    // Deterministic spread-out initial design.
    let mut xs: Vec<Vec<f64>> = vec![vec![0.4], vec![2.6], vec![3.6]];
    let mut ys: Vec<f64> = xs.iter().map(|x| objective(x[0])).collect();

    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.5],
    };
    // H₀: "still improving ≥ 70% of the time" (stall rate ≤ 0.3).
    let mut eproc = BettingEProcess::new(0.3);
    let mut cs = GaussianMixtureCs::new(1.0, 4.0, alpha);
    let (_, mut best_value) = argmin_flat(&xs, &ys);
    cs.observe(best_value);

    let mut iterations = 0usize;
    let mut stopped_early = false;
    for _ in 0..max_iters {
        let Some(gp) = Gp::try_fit(&xs, &ys, kernel.clone(), 1e-6) else {
            break;
        };
        let f_best = ys.iter().copied().fold(f64::INFINITY, f64::min);
        // Pick the grid point of maximum Expected Improvement.
        let mut best_ei = f64::NEG_INFINITY;
        let mut x_next = grid[0];
        for &g in &grid {
            let ei = expected_improvement(&gp, &[g], f_best, 0.01);
            if ei > best_ei {
                best_ei = ei;
                x_next = g;
            }
        }
        let y_next = objective(x_next);
        xs.push(vec![x_next]);
        ys.push(y_next);
        iterations += 1;

        let new_best = best_value.min(y_next);
        let improvement = best_value - new_best;
        best_value = new_best;
        cs.observe(new_best);
        // Stall indicator ∈ {0,1}: 1 = this step did NOT improve by ≥ δ.
        let stall = if improvement < delta { 1.0 } else { 0.0 };
        eproc.observe(stall);
        if eproc.rejects_at(alpha) {
            stopped_early = true;
            break;
        }
    }

    let (best_x, best_value) = argmin_flat(&xs, &ys);
    let (ci_center, ci_radius) = cs.interval().unwrap_or((best_value, f64::INFINITY));
    let stop_evidence = StatisticalCertificate::EValue {
        e: eproc.e_value(),
        alpha,
    };
    let incumbent_color = Color::Estimated {
        estimator: if stopped_early {
            "observed-incumbent-at-e-stopping-time"
        } else {
            "observed-incumbent-at-budget-end"
        }
        .to_string(),
        dispersion: ci_radius,
    };

    AdaptBoReport {
        best_x,
        best_value,
        iterations,
        evaluations: ys.len(),
        stopped_early,
        log_e_value: eproc.log_e_value(),
        ci_center,
        ci_radius,
        stop_evidence,
        incumbent_color,
        surrogate_color: Color::Estimated {
            estimator: "gp-matern52".to_string(),
            dispersion: ci_radius,
        },
    }
}

fn argmin_flat(xs: &[Vec<f64>], ys: &[f64]) -> (f64, f64) {
    let flat: Vec<f64> = xs.iter().map(|x| x[0]).collect();
    argmin(&flat, ys)
}
