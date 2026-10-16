//! CG with a recomputed-residual acceptance gate and bounded defect correction.
//!
//! This is an opt-in driver over the original [`CgState`], not a new Krylov
//! recurrence. A small recursive residual is only a request to check `b - Ax`.
//! When that check fails, remaining work may solve for the defect instead of
//! discarding an otherwise useful iterate. Every attempt shares one iteration
//! budget. The caller still owns SPD admission and the cancellation scope.
//!
//! "Checked" means an explicit floating-point operator application. It is not
//! an interval enclosure, a condition-number bound, or a forward-error claim.

use crate::{CgState, LinearOp, ResidualClaim, SolveReport, StallDiagnosis, norm2};
use fs_sparse::precond::Precond;

/// Explicit accuracy and work limits for [`checked_cg`].
#[derive(Debug, Clone, Copy)]
pub struct CheckedCgConfig {
    /// Strict relative Euclidean residual target, finite and in `(0, 1)`.
    pub tolerance: f64,
    /// Total CG iterations, including every defect-correction solve.
    pub max_iterations: usize,
    /// Additional solves after the initial CG attempt. Zero disables repair,
    /// but never disables the recomputed-residual acceptance gate.
    pub max_corrections: usize,
}

/// An input refusal or the caller's unchanged interruption reason.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckedCgError<E> {
    /// Invalid dimensions, non-finite RHS, or invalid tolerance.
    InvalidInput(&'static str),
    /// The caller declined further work. No partial solution is published.
    Interrupted(E),
}

impl<E: std::fmt::Display> std::fmt::Display for CheckedCgError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(what) => write!(f, "checked CG: {what}"),
            Self::Interrupted(error) => write!(f, "checked CG interrupted: {error}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for CheckedCgError<E> {}

/// Best recomputed iterate, including a truthful report on nonconvergence.
#[derive(Debug, Clone)]
pub struct CheckedCgSolution {
    /// The iterate whose residual is carried by `report`.
    pub x: Vec<f64>,
    /// Always `TrueEuclidean`; `converged_euclidean()` is the acceptance test.
    /// History contains the initial residual and accepted outer corrections,
    /// not the inner CG recurrence's estimates.
    pub report: SolveReport,
    /// Number of additional defect solves actually attempted.
    pub corrections: usize,
}

/// Solve an admitted SPD system with bounded, cancellable defect correction.
///
/// The preconditioner is reused. Each residual RHS is scaled before entering
/// CG, avoiding squared-norm underflow/overflow for very small/large loads.
/// Acceptance is checked on the *returned*, rescaled iterate in the original
/// system. A correction that makes no measured progress is rejected, and the
/// previous best iterate is retained. No tolerance or iteration budget expands.
///
/// `checkpoint` receives cumulative CG iterations before/after operator work
/// and between batches of at most 16 CG iterations. Operator/preconditioner
/// applications themselves must supply any finer-grained cancellation needed.
/// An interruption returns the caller's error, never a partially accepted solve.
/// Ordinary nonconvergence returns a solution with a nonconverged report.
///
/// # Errors
/// Returns [`CheckedCgError::InvalidInput`] before numerical work for invalid
/// inputs, or [`CheckedCgError::Interrupted`] when the callback refuses work.
pub fn checked_cg<A: LinearOp, P: Precond, E>(
    a: &A,
    preconditioner: &P,
    b: &[f64],
    config: CheckedCgConfig,
    mut checkpoint: impl FnMut(usize) -> Result<(), E>,
) -> Result<CheckedCgSolution, CheckedCgError<E>> {
    if b.len() != a.n() {
        return Err(CheckedCgError::InvalidInput("RHS length does not match operator"));
    }
    if !b.iter().all(|v| v.is_finite()) {
        return Err(CheckedCgError::InvalidInput("RHS must be finite"));
    }
    if !(config.tolerance.is_finite() && config.tolerance > 0.0 && config.tolerance < 1.0) {
        return Err(CheckedCgError::InvalidInput("tolerance must lie in (0, 1)"));
    }
    checkpoint(0).map_err(CheckedCgError::Interrupted)?;
    let mut x = vec![0.0; b.len()];
    let mut residual = b.to_vec();
    let bscale = max_abs(b);
    let bnorm = if bscale == 0.0 { 1.0 } else { scaled_norm(b, bscale) };
    let mut relative = if bscale == 0.0 { 0.0 } else { 1.0 };
    let mut history = vec![relative];
    let mut iterations = 0;
    let mut corrections = 0;
    let mut diagnosis = StallDiagnosis::BudgetExhausted;

    while relative >= config.tolerance && iterations < config.max_iterations {
        checkpoint(iterations).map_err(CheckedCgError::Interrupted)?;
        let scale = max_abs(&residual);
        let rhs: Vec<f64> = residual.iter().map(|v| v / scale).collect();
        // The first attempt uses the requested tolerance. A defect solve must
        // reduce its own RHS by at least a factor of two; asking it to meet a
        // relative tolerance near one can merely reproduce the failed check.
        let inner_tolerance = if iterations == 0 {
            config.tolerance
        } else {
            (config.tolerance / relative).min(0.5)
        };
        let mut state = CgState::new(a, preconditioner, &rhs);
        let mut broken = false;
        while state.rel_residual() >= inner_tolerance && iterations < config.max_iterations {
            checkpoint(iterations).map_err(CheckedCgError::Interrupted)?;
            let before = state.iters;
            let report = state.run(a, preconditioner, inner_tolerance,
                (config.max_iterations - iterations).min(16));
            iterations += state.iters - before;
            checkpoint(iterations).map_err(CheckedCgError::Interrupted)?;
            if !report.rel_residual.is_finite() || !state.x.iter().all(|v| v.is_finite()) {
                broken = true;
                break;
            }
            if state.iters == before { break; }
        }
        if broken {
            diagnosis = StallDiagnosis::Breakdown;
            break;
        }
        let candidate: Vec<f64> = x.iter().zip(&state.x)
            .map(|(&old, &delta)| scale.mul_add(delta, old)).collect();
        if !candidate.iter().all(|v| v.is_finite()) {
            diagnosis = StallDiagnosis::Breakdown;
            break;
        }
        let mut applied = vec![0.0; b.len()];
        checkpoint(iterations).map_err(CheckedCgError::Interrupted)?;
        a.apply(&candidate, &mut applied);
        checkpoint(iterations).map_err(CheckedCgError::Interrupted)?;
        let defect: Vec<f64> = b.iter().zip(applied).map(|(&bi, ai)| bi - ai).collect();
        if !defect.iter().all(|v| v.is_finite()) {
            diagnosis = StallDiagnosis::Breakdown;
            break;
        }
        let rscale = max_abs(&defect);
        let measured = if rscale == 0.0 { 0.0 } else {
            (rscale / bscale) * (scaled_norm(&defect, rscale) / bnorm)
        };
        if !measured.is_finite() {
            diagnosis = StallDiagnosis::Breakdown;
            break;
        }
        if measured >= relative {
            diagnosis = StallDiagnosis::Plateau;
            break;
        }
        x = candidate;
        residual = defect;
        relative = measured;
        history.push(relative);
        if relative < config.tolerance || iterations == config.max_iterations
            || corrections == config.max_corrections
        { break; }
        corrections += 1;
    }
    checkpoint(iterations).map_err(CheckedCgError::Interrupted)?;
    Ok(CheckedCgSolution {
        x,
        report: SolveReport::from_claim_with_diagnosis(iterations,
            ResidualClaim::TrueEuclidean(relative), config.tolerance, history, diagnosis),
        corrections,
    })
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|v| v.abs()).fold(0.0, f64::max)
}

fn scaled_norm(values: &[f64], scale: f64) -> f64 {
    norm2(&values.iter().map(|v| v / scale).collect::<Vec<_>>())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_sparse::precond::IdentityPrecond;
    use std::convert::Infallible;

    struct Dense(Vec<Vec<f64>>);
    impl LinearOp for Dense {
        fn n(&self) -> usize { self.0.len() }
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            for (yi, row) in y.iter_mut().zip(&self.0) {
                *yi = crate::dot(row, x);
            }
        }
    }
    fn config() -> CheckedCgConfig {
        CheckedCgConfig { tolerance: 1e-12, max_iterations: 64, max_corrections: 2 }
    }
    fn run(a: &Dense, b: &[f64], cfg: CheckedCgConfig) -> CheckedCgSolution {
        checked_cg(a, &IdentityPrecond, b, cfg, |_| Ok::<_, Infallible>(())).unwrap()
    }

    #[test]
    fn checked_cg_measures_returned_solution_at_extreme_load_scales() {
        let a = Dense(vec![vec![4.0, 1.0], vec![1.0, 3.0]]);
        for scale in [1.0, 1e-250, 1e250] {
            let b = [6.0 * scale, 7.0 * scale];
            let solved = run(&a, &b, config());
            assert!(solved.report.converged_euclidean(), "{:?}", solved.report);
            assert!((solved.x[0] / scale - 1.0).abs() < 1e-12);
            assert!((solved.x[1] / scale - 2.0).abs() < 1e-12);
            let mut applied = vec![0.0; 2];
            a.apply(&solved.x, &mut applied);
            let independent = ((b[0]-applied[0])/scale).hypot((b[1]-applied[1])/scale)
                / 6.0_f64.hypot(7.0);
            assert!(independent < config().tolerance);
            assert!((solved.report.rel_residual-independent).abs() < 1e-28);
        }
    }

    #[test]
    fn checked_cg_repairs_a_real_recursive_residual_gap_without_relaxing_tolerance() {
        // A = G^T G + I, G = [[-39,21,46],[-28,48,36],[-5,-36,4]].
        // Hence SPD independently of the solver. Its cancellation-prone rows
        // expose an actual residual gap, not a mocked report or perturbed state.
        let a = Dense(vec![vec![2331.0, -1983.0, -2822.0],
            vec![-1983.0, 4042.0, 2550.0], vec![-2822.0, 2550.0, 3429.0]]);
        let b = [1.0, 0.0, 0.0];
        let mut cfg = config(); cfg.tolerance = 1e-13; cfg.max_corrections = 0;
        let mut original = CgState::new(&a, &IdentityPrecond, &b);
        let recurrence = original.run(&a, &IdentityPrecond, cfg.tolerance, cfg.max_iterations);
        assert!(recurrence.converged);
        assert!(!recurrence.converged_euclidean());
        let unchecked = run(&a, &b, cfg);
        assert!(!unchecked.report.converged_euclidean(), "{:?}", unchecked.report);
        assert_eq!(unchecked.corrections, 0);
        assert_eq!(unchecked.report.iters, recurrence.iters);

        cfg.max_corrections = 2;
        let repaired = run(&a, &b, cfg);
        assert!(repaired.report.converged_euclidean(), "{:?}", repaired.report);
        assert!(repaired.corrections > 0 && repaired.corrections <= 2);
        assert!(repaired.report.iters > recurrence.iters);
        assert!(repaired.report.iters <= cfg.max_iterations);
        let mut applied = vec![0.0; 3]; a.apply(&repaired.x, &mut applied);
        let truth = (1.0-applied[0]).hypot(applied[1]).hypot(applied[2]);
        assert!(truth < cfg.tolerance, "actual residual {truth}");
        assert!(repaired.report.history.windows(2).all(|w| w[1] < w[0]));

        cfg.max_iterations = recurrence.iters;
        let exhausted = run(&a, &b, cfg);
        assert_eq!(exhausted.report.iters, cfg.max_iterations);
        assert!(!exhausted.report.converged_euclidean());
        assert_eq!(exhausted.x, unchecked.x);
        cfg.max_iterations = 64;
        let interrupted = checked_cg(&a, &IdentityPrecond, &b, cfg, |iters|
            if iters > recurrence.iters { Err("stop during repair") } else { Ok(()) });
        assert!(matches!(interrupted, Err(CheckedCgError::Interrupted("stop during repair"))));
        assert_eq!(run(&a, &b, cfg).x, repaired.x);
    }

    #[test]
    fn checked_cg_zero_rhs_and_zero_budget_do_no_operator_work() {
        struct Never;
        impl LinearOp for Never {
            fn n(&self) -> usize { 2 }
            fn apply(&self, _: &[f64], _: &mut [f64]) { panic!("unfunded work"); }
        }
        let zero = checked_cg(&Never, &IdentityPrecond, &[0.0, 0.0], config(),
            |_| Ok::<_, Infallible>(())).unwrap();
        assert!(zero.report.converged_euclidean());
        assert_eq!(zero.report.iters, 0);
        let mut cfg = config(); cfg.max_iterations = 0;
        let unfunded = checked_cg(&Never, &IdentityPrecond, &[1.0, 2.0], cfg,
            |_| Ok::<_, Infallible>(())).unwrap();
        assert!(!unfunded.report.converged);
        assert_eq!(unfunded.report.rel_residual, 1.0);
        assert_eq!(unfunded.report.iters, 0);
    }

    #[test]
    fn checked_cg_shares_budget_and_preserves_interruption_reason() {
        let a = Dense(vec![vec![4.0, 1.0], vec![1.0, 3.0]]);
        let mut cfg = config(); cfg.max_iterations = 1;
        let short = run(&a, &[1.0, 0.0], cfg);
        assert_eq!(short.report.iters, 1);
        assert!(!short.report.converged);
        let cancelled = checked_cg(&a, &IdentityPrecond, &[1.0, 0.0], config(),
            |iters| if iters > 0 { Err("stop inside CG") } else { Ok(()) });
        assert!(matches!(cancelled, Err(CheckedCgError::Interrupted("stop inside CG"))));
        let retry = run(&a, &[1.0, 0.0], config());
        let fresh = run(&a, &[1.0, 0.0], config());
        assert_eq!(retry.x, fresh.x);
        assert_eq!(retry.report.history, fresh.report.history);
    }

    #[test]
    fn checked_cg_refuses_invalid_inputs_and_singular_breakdown() {
        let a = Dense(vec![vec![1.0]]);
        for tolerance in [0.0, -1.0, 1.0, f64::NAN, f64::INFINITY] {
            let mut cfg = config(); cfg.tolerance = tolerance;
            assert!(matches!(checked_cg(&a, &IdentityPrecond, &[1.0], cfg,
                |_| Ok::<_, Infallible>(())), Err(CheckedCgError::InvalidInput(_))));
        }
        assert!(matches!(checked_cg(&a, &IdentityPrecond, &[f64::NAN], config(),
            |_| Ok::<_, Infallible>(())), Err(CheckedCgError::InvalidInput(_))));
        assert!(matches!(checked_cg(&a, &IdentityPrecond, &[], config(),
            |_| Ok::<_, Infallible>(())), Err(CheckedCgError::InvalidInput(_))));
        let singular = run(&Dense(vec![vec![0.0]]), &[1.0], config());
        assert!(!singular.report.converged_euclidean());
        assert_eq!(singular.report.diagnosis, Some(StallDiagnosis::Breakdown));
        assert_eq!(singular.x, [0.0]);
    }
}
