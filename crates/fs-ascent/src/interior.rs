//! Interior-point option (bead ijil, §9.2): the LOG-BARRIER method —
//! the bead's sanctioned simple variant. Inequalities become barrier
//! terms −μ·Σ ln(−cᵢ), equalities ride the augmented-Lagrangian term
//! (λᵀc + ρ/2‖c‖²), and each barrier subproblem is minimized by the
//! existing resumable L-BFGS. μ decreases geometrically; the barrier
//! multiplier estimates νⱼ = μ/(−cᵢⱼ) feed the SAME
//! [`crate::auglag::kkt_residual`] certificate the AL path uses —
//! converged and stalled stay distinguishable outcomes, and parity
//! with AL on the landed KKT fixtures is gated in the battery.
//!
//! Strictly feasible interior start required for the inequalities (the
//! method's nature); infeasible starts are nudged by a feasibility
//! phase-1 (minimize max(cᵢ) by L-BFGS on a smooth softmax proxy)
//! before the barrier loop.

use crate::auglag::{
    ConstrainedProblem, KktResidual, assert_finite, checked_constraints, checked_fg, checked_jt,
    kkt_residual, validate_problem_at_start, validate_tolerance,
};
use crate::lbfgs::LbfgsState;
use crate::stop::StopRule;

/// Outcome of an interior-point solve (mirrors [`crate::AugLagReport`]).
#[derive(Debug, Clone)]
pub struct InteriorReport {
    /// Final iterate.
    pub x: Vec<f64>,
    /// Final objective.
    pub f: f64,
    /// The certificate.
    pub kkt: KktResidual,
    /// Equality multipliers λ.
    pub lambda: Vec<f64>,
    /// Barrier inequality multipliers ν = μ/(−cᵢ) ≥ 0.
    pub nu: Vec<f64>,
    /// Barrier outer iterations (μ steps).
    pub outer_iters: usize,
    /// Total inner evaluations.
    pub evals: usize,
    /// All four KKT residuals below tolerance.
    pub converged: bool,
}

/// Phase-1: from `x0`, reduce max(cᵢ) below −margin with L-BFGS on the
/// smooth log-sum-exp proxy. Returns a strictly feasible point or the
/// best found (the barrier loop then reports honest non-convergence).
fn phase1(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    margin: f64,
    ni: usize,
) -> (Vec<f64>, usize) {
    let ci0 = checked_constraints("inequality", problem.ci, x0, Some(ni));
    if ci0.iter().all(|&c| c < -margin) {
        return (x0.to_vec(), 0);
    }
    let beta = 30.0; // softmax sharpness
    let mut evals = 0usize;
    let ce = problem.ce;
    let ci = problem.ci;
    let ci_jt = problem.ci_jt;
    let _ = ce;
    let mut fg = |x: &[f64]| -> (f64, Vec<f64>) {
        evals += 1;
        let c = checked_constraints("inequality", ci, x, Some(ni));
        let m = c.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let ws: Vec<f64> = c
            .iter()
            .map(|&v| fs_math::det::exp(beta * (v - m)))
            .collect();
        let sum: f64 = ws.iter().sum();
        let val = m + fs_math::det::ln(sum) / beta;
        let w: Vec<f64> = ws.iter().map(|v| v / sum).collect();
        let g = checked_jt("inequality", ci_jt, x, &w);
        assert!(val.is_finite(), "phase-1 objective must remain finite");
        assert_finite("phase-1 gradient", &g);
        (val, g)
    };
    let mut st = LbfgsState::new(x0, 8, &mut fg);
    let rule = StopRule::Any(vec![
        StopRule::ObjectiveBelow(-2.0 * margin),
        StopRule::GradNorm(1e-12),
    ]);
    st.run(&mut fg, &rule, 300);
    (st.x.clone(), evals)
}

/// Run the log-barrier interior-point loop from `x0`.
///
/// `tol` gates the KKT certificate; `max_outer` caps μ reductions.
#[allow(clippy::too_many_lines)] // One ordered barrier solve and certificate protocol.
pub fn interior_point(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    tol: f64,
    max_outer: usize,
) -> InteriorReport {
    validate_tolerance(tol);
    let (ne, ni) = validate_problem_at_start(problem, x0);
    let (mut x, mut evals) = if ni > 0 {
        phase1(problem, x0, 1e-9, ni)
    } else {
        (x0.to_vec(), 0)
    };
    let mut mu = 1.0f64;
    let mut rho = 10.0f64;
    let mut lambda = vec![0.0f64; ne];
    let mut nu = vec![0.0f64; ni];
    let mut outer = 0usize;
    for _ in 0..max_outer {
        outer += 1;
        // Barrier + AL subproblem.
        let (lam, m, r) = (lambda.clone(), mu, rho);
        let fg_cb = &mut *problem.fg;
        let ce = problem.ce;
        let ce_jt = problem.ce_jt;
        let ci = problem.ci;
        let ci_jt = problem.ci_jt;
        let mut inner_evals = 0usize;
        let mut inner = |xv: &[f64]| -> (f64, Vec<f64>) {
            inner_evals += 1;
            let (f, mut g) = checked_fg(fg_cb, xv);
            let mut val = f;
            let cev = checked_constraints("equality", ce, xv, Some(ne));
            if !cev.is_empty() {
                let w: Vec<f64> = cev
                    .iter()
                    .zip(&lam)
                    .map(|(c, l)| r.mul_add(*c, *l))
                    .collect();
                for (c, l) in cev.iter().zip(&lam) {
                    val += l * c + 0.5 * r * c * c;
                }
                let pull = checked_jt("equality", ce_jt, xv, &w);
                for i in 0..g.len() {
                    g[i] += pull[i];
                }
            }
            let civ = checked_constraints("inequality", ci, xv, Some(ni));
            if !civ.is_empty() {
                // Log barrier; infeasible samples get +inf (the line
                // search backtracks into the interior).
                let mut w = vec![0.0f64; civ.len()];
                for (j, &c) in civ.iter().enumerate() {
                    if c >= 0.0 {
                        return (f64::INFINITY, vec![0.0; xv.len()]);
                    }
                    val -= m * fs_math::det::ln(-c);
                    w[j] = m / (-c);
                }
                let pull = checked_jt("inequality", ci_jt, xv, &w);
                for i in 0..g.len() {
                    g[i] += pull[i];
                }
            }
            assert!(val.is_finite(), "barrier objective must remain finite");
            assert_finite("barrier gradient", &g);
            (val, g)
        };
        let mut st = LbfgsState::new(&x, 10, &mut inner);
        let rule = StopRule::GradNorm((mu * 0.1).max(tol * 0.1));
        st.run(&mut inner, &rule, 400);
        x.clone_from(&st.x);
        evals += inner_evals;
        // Multiplier updates.
        let cev = checked_constraints("equality", problem.ce, &x, Some(ne));
        for i in 0..lambda.len() {
            lambda[i] += rho * cev[i];
        }
        assert_finite("equality multiplier", &lambda);
        let civ = checked_constraints("inequality", problem.ci, &x, Some(ni));
        nu = civ.iter().map(|&c| mu / (-c).max(1e-300)).collect();
        assert_finite("inequality multiplier", &nu);
        // Certificate check at the CURRENT multipliers.
        let kkt = kkt_residual(problem, &x, &lambda, &nu);
        evals += 1;
        if kkt.within_tolerance(tol) {
            let (f, _) = checked_fg(&mut *problem.fg, &x);
            return InteriorReport {
                x,
                f,
                kkt,
                lambda,
                nu,
                outer_iters: outer,
                evals,
                converged: true,
            };
        }
        mu *= 0.2;
        rho *= 2.0;
    }
    // `nu` belongs to the last barrier subproblem actually solved. The
    // loop reduces `mu` only in preparation for a possible next solve;
    // recomputing here with that unsolved parameter can false-certify.
    assert_finite("inequality multiplier", &nu);
    let kkt = kkt_residual(problem, &x, &lambda, &nu);
    let (f, _) = checked_fg(&mut *problem.fg, &x);
    let converged = kkt.within_tolerance(tol);
    InteriorReport {
        x,
        f,
        kkt,
        lambda,
        nu,
        outer_iters: outer,
        evals,
        converged,
    }
}
