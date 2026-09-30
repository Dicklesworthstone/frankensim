//! Noisy outcome-constrained design search (plan sections 9.4 and 9.9).
//!
//! Each output has an independent GP, but baseline/candidate correlations
//! WITHIN each output are retained. Feasibility is evaluated on joint latent
//! draws, not on noisy observations or posterior means. This is a model-based
//! acquisition, not a safety guarantee or an anytime-valid probability bound.

use crate::gp::Gp;

mod bo;
pub use bo::{ConstrainedBoConfig, ConstrainedBoError, ConstrainedBoReport,
    ConstrainedObservation, ConstrainedRecommendation, OutcomeConstraint,
    OutputFitRecord, minimize_constrained};

/// One modeled upper constraint, `g(x) <= upper_bound`.
///
/// `gp` predicts centered outputs. Supply the bound in those same centered
/// units; the closed-loop driver handles prior-mean subtraction for callers.
pub struct ConstraintModel<'a> {
    /// Independent latent GP for this outcome.
    pub gp: &'a Gp,
    /// Finite upper bound in the GP's output units.
    pub upper_bound: f64,
}

/// Expected improvement in reference-capped feasible utility, for minimization.
///
/// For each joint draw define
/// `U(S) = max({0} union {reference - f(x): x in S, all g_j(x) <= b_j})`.
/// Return the sample mean of `U(baseline union candidates) - U(baseline)`.
/// The explicit finite reference defines zero utility when no point is feasible
/// and gives no reward to objectives above it. It is a preference, NOT a bound
/// inferred from the data. This differs from uncapped q-NEI on such draws.
///
/// `bank` is row-major. Each row has `(baseline.len() + candidates.len()) *
/// (1 + constraints.len())` columns, POINT first, OUTPUT second: the objective
/// then the declared constraints. Reuse row-wise prefixes of the same bank
/// while growing a batch. Under independent standard normal columns this
/// represents independent outputs with each output's joint spatial posterior.
///
/// Equal coordinates, including signed zero, share one latent draw per output.
/// Only their first occurrence consumes bank columns; a baseline duplicate
/// cannot manufacture improvement through Cholesky jitter. Constraints use an
/// inclusive upper bound. An empty constraint list is permitted and yields
/// reference-capped unconstrained improvement.
///
/// Inherits `Gp::predict_joint` adaptive jitter and degenerate covariance
/// fallback. No output cross-covariance, global optimum, safe-evaluation policy,
/// exact integration, or cross-ISA replay guarantee is asserted.
///
/// # Panics
/// Empty baseline/candidates, incompatible dimensions, non-finite inputs or
/// bounds, malformed banks, shape overflow and non-finite posterior arithmetic
/// refuse instead of returning an apparent acquisition value.
#[must_use]
pub fn q_feasible_noisy_improvement(
    objective: &Gp,
    constraints: &[ConstraintModel<'_>],
    baseline: &[Vec<f64>],
    candidates: &[Vec<f64>],
    reference: f64,
    bank: &[f64],
) -> f64 {
    assert!(!baseline.is_empty(), "constrained q-NEI needs a baseline");
    assert!(!candidates.is_empty(), "constrained q-NEI needs candidates");
    assert!(reference.is_finite(), "constrained q-NEI reference must be finite");
    let dim = objective.kernel.lengthscales.len();
    assert!(dim > 0, "constrained q-NEI needs input dimensions");
    for constraint in constraints {
        assert_eq!(constraint.gp.kernel.lengthscales.len(), dim, "constraint GP dimension mismatch");
        assert!(constraint.upper_bound.is_finite(), "constraint bound must be finite");
    }
    for point in baseline.iter().chain(candidates) {
        assert_eq!(point.len(), dim, "constrained q-NEI input dimension mismatch");
        assert!(point.iter().all(|v| v.is_finite()), "constrained q-NEI coordinates must be finite");
    }
    let points = baseline.len().checked_add(candidates.len()).expect("point count overflow");
    let outputs = constraints.len().checked_add(1).expect("output count overflow");
    let width = points.checked_mul(outputs).expect("constrained bank width overflow");
    points.checked_mul(points).expect("constrained posterior shape overflow");
    assert!(!bank.is_empty() && bank.len() % width == 0, "constrained bank must be nonempty and rectangular");
    assert!(bank.iter().all(|v| v.is_finite()), "constrained bank must be finite");

    let mut unique: Vec<Vec<f64>> = Vec::new();
    let mut columns = Vec::new();
    let mut baseline_len = 0;
    for (index, point) in baseline.iter().chain(candidates).enumerate() {
        if !unique.iter().any(|previous| previous == point) {
            unique.push(point.clone());
            columns.push(index);
        }
        if index + 1 == baseline.len() { baseline_len = unique.len(); }
    }
    if baseline_len == unique.len() { return 0.0; }
    let posteriors: Vec<_> = std::iter::once(objective).chain(constraints.iter().map(|c| c.gp))
        .map(|gp| gp.predict_joint(&unique)).collect();
    assert!(posteriors.iter().all(|(mean, lower)| mean.iter().chain(lower).all(|v| v.is_finite())),
        "constrained posterior must be finite");
    let n = unique.len();
    let mut average = 0.0;
    for (sample, row) in bank.chunks_exact(width).enumerate() {
        // Minima capped by the same reference are numerically equivalent to
        // utility maxima, without subtracting two large reference-offset gains.
        let mut old_best = reference;
        let mut new_best = reference;
        for i in 0..n {
            let mut feasible = true;
            let mut f = 0.0;
            for (output, (mean, lower)) in posteriors.iter().enumerate() {
                let mut value = mean[i];
                for j in 0..=i {
                    value = lower[i * n + j].mul_add(row[columns[j] * outputs + output], value);
                }
                assert!(value.is_finite(), "constrained posterior draw overflow");
                if output == 0 { f = value; }
                else { feasible &= value <= constraints[output - 1].upper_bound; }
            }
            if feasible {
                if i < baseline_len { old_best = old_best.min(f); }
                else { new_best = new_best.min(f); }
            }
        }
        let gain = (old_best - new_best).max(0.0);
        assert!(gain.is_finite(), "constrained improvement overflow");
        average += (gain - average) / (sample + 1) as f64;
    }
    average
}

#[cfg(test)]
mod tests;
