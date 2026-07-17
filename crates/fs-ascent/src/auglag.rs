//! Augmented Lagrangian — the robust constrained default: minimize
//! f(x) s.t. c_e(x) = 0, c_i(x) ≤ 0 by an outer multiplier loop over
//! L-BFGS inner solves of
//! L_μ(x) = f + λᵀc_e + (μ/2)‖c_e‖² + (μ/2)‖max(0, c_i + s/μ)‖²-style
//! terms (the standard PHR form for inequalities). Every returned
//! optimum carries a KKT-RESIDUAL CERTIFICATE — stationarity, primal
//! feasibility, dual feasibility, and complementarity — so
//! "converged" and "stalled" are distinguishable outcomes, not vibes.

use crate::lbfgs::LbfgsState;
use crate::stop::StopRule;

/// The KKT residuals of a returned point (the certificate).
#[derive(Debug, Clone, PartialEq)]
pub struct KktResidual {
    /// ‖∇f + Σλ∇c_e + Σν∇c_i‖∞ (stationarity of the Lagrangian).
    pub stationarity: f64,
    /// max(‖c_e‖∞, ‖max(0, c_i)‖∞) (feasibility).
    pub feasibility: f64,
    /// max max(0, -ν_j) (inequality dual feasibility).
    pub dual_feasibility: f64,
    /// max |ν_j · c_i_j| (complementary slackness).
    pub complementarity: f64,
}

impl KktResidual {
    /// Whether every KKT residual is strictly below `tol`.
    ///
    /// # Panics
    ///
    /// Panics when `tol` is not finite and positive. Non-finite residuals
    /// never satisfy the predicate.
    #[must_use]
    pub fn within_tolerance(&self, tol: f64) -> bool {
        validate_tolerance(tol);
        [
            self.stationarity,
            self.feasibility,
            self.dual_feasibility,
            self.complementarity,
        ]
        .into_iter()
        .all(|residual| residual.is_finite() && residual < tol)
    }
}

/// Outcome of an augmented-Lagrangian solve.
#[derive(Debug, Clone)]
pub struct AugLagReport {
    /// Final iterate.
    pub x: Vec<f64>,
    /// Final objective (of f, not the Lagrangian).
    pub f: f64,
    /// The certificate.
    pub kkt: KktResidual,
    /// Equality multipliers λ.
    pub lambda: Vec<f64>,
    /// Inequality multipliers ν ≥ 0.
    pub nu: Vec<f64>,
    /// Outer iterations.
    pub outer_iters: usize,
    /// Total inner evaluations.
    pub evals: usize,
    /// All four KKT residuals below the tolerance.
    pub converged: bool,
}

/// Problem callbacks: objective+gradient, equality constraints and
/// their Jacobian-transpose action, inequalities likewise.
pub struct ConstrainedProblem<'a> {
    /// (f, ∇f).
    pub fg: crate::FnGrad<'a>,
    /// c_e(x) (empty vec for none).
    pub ce: crate::FnConstraints<'a>,
    /// (∂c_e/∂x)ᵀ·w.
    pub ce_jt: crate::FnJacobianTranspose<'a>,
    /// c_i(x) (≤ 0 feasible; empty vec for none).
    pub ci: crate::FnConstraints<'a>,
    /// (∂c_i/∂x)ᵀ·w.
    pub ci_jt: crate::FnJacobianTranspose<'a>,
}

fn inf_norm(v: &[f64]) -> f64 {
    assert_finite("norm input", v);
    v.iter().map(|x| x.abs()).fold(0.0f64, f64::max)
}

pub(crate) fn validate_tolerance(tol: f64) {
    assert!(
        tol.is_finite() && tol > 0.0,
        "KKT tolerance must be finite and positive"
    );
}

pub(crate) fn validate_point(x: &[f64]) {
    assert!(
        !x.is_empty(),
        "constrained decision vector must be non-empty"
    );
    assert_finite("constrained decision vector", x);
}

pub(crate) fn assert_finite(label: &str, values: &[f64]) {
    assert!(
        values.iter().all(|value| value.is_finite()),
        "{label} entries must be finite"
    );
}

pub(crate) fn checked_fg(fg: crate::FnGrad<'_>, x: &[f64]) -> (f64, Vec<f64>) {
    validate_point(x);
    let (f, gradient) = fg(x);
    assert!(f.is_finite(), "objective value must be finite");
    assert_eq!(
        gradient.len(),
        x.len(),
        "objective gradient length must match the decision dimension"
    );
    assert_finite("objective gradient", &gradient);
    (f, gradient)
}

pub(crate) fn checked_constraints(
    label: &str,
    callback: crate::FnConstraints<'_>,
    x: &[f64],
    expected: Option<usize>,
) -> Vec<f64> {
    validate_point(x);
    let values = callback(x);
    if let Some(expected) = expected {
        assert_eq!(
            values.len(),
            expected,
            "{label} constraint count changed across evaluations"
        );
    }
    assert_finite(&format!("{label} constraint"), &values);
    values
}

pub(crate) fn checked_jt(
    label: &str,
    callback: crate::FnJacobianTranspose<'_>,
    x: &[f64],
    weights: &[f64],
) -> Vec<f64> {
    validate_point(x);
    assert_finite(&format!("{label} multiplier"), weights);
    let pullback = callback(x, weights);
    assert_eq!(
        pullback.len(),
        x.len(),
        "{label} Jacobian-transpose output length must match the decision dimension"
    );
    assert_finite(&format!("{label} Jacobian-transpose output"), &pullback);
    pullback
}

fn validate_jt_zero(
    label: &str,
    callback: crate::FnJacobianTranspose<'_>,
    x: &[f64],
    constraint_count: usize,
) {
    if constraint_count == 0 {
        return;
    }
    let zero_pullback = checked_jt(label, callback, x, &vec![0.0; constraint_count]);
    assert!(
        zero_pullback.iter().all(|value| *value == 0.0),
        "{label} Jacobian-transpose action must map zero weights to zero"
    );
}

pub(crate) fn validate_problem_at_start(
    problem: &mut ConstrainedProblem<'_>,
    x: &[f64],
) -> (usize, usize) {
    validate_point(x);
    let _ = checked_fg(&mut *problem.fg, x);
    let ce = checked_constraints("equality", problem.ce, x, None);
    let ci = checked_constraints("inequality", problem.ci, x, None);
    validate_jt_zero("equality", problem.ce_jt, x, ce.len());
    validate_jt_zero("inequality", problem.ci_jt, x, ci.len());
    (ce.len(), ci.len())
}

/// Run the PHR augmented-Lagrangian loop from `x0` with cold (zero)
/// multipliers.
pub fn augmented_lagrangian(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    tol: f64,
    max_outer: usize,
) -> AugLagReport {
    let (ne, ni) = validate_problem_at_start(problem, x0);
    augmented_lagrangian_warm(problem, x0, tol, max_outer, &vec![0.0; ne], &vec![0.0; ni])
}

/// Run the PHR augmented-Lagrangian loop from `x0` with WARM-STARTED
/// multipliers (continuation solves, e.g. ε-constraint sweeps). Cold
/// starts on an objective whose unconstrained descent runs to an
/// asymptote can escape before the penalty activates; carrying the
/// previous solve's multipliers keeps constraint pressure alive from
/// the first inner iteration.
pub fn augmented_lagrangian_warm(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    tol: f64,
    max_outer: usize,
    lambda0: &[f64],
    nu0: &[f64],
) -> AugLagReport {
    validate_tolerance(tol);
    let (ne, ni) = validate_problem_at_start(problem, x0);
    assert_eq!(
        lambda0.len(),
        ne,
        "warm equality multipliers must match the constraint dimension"
    );
    assert_eq!(
        nu0.len(),
        ni,
        "warm inequality multipliers must match the constraint dimension"
    );
    assert_finite("warm equality multiplier", lambda0);
    assert_finite("warm inequality multiplier", nu0);
    assert!(
        nu0.iter().all(|v| *v >= 0.0),
        "warm inequality multipliers must be nonnegative"
    );
    let mut lambda = lambda0.to_vec();
    let mut nu = nu0.to_vec();
    let mut mu = 10.0f64;
    let mut x = x0.to_vec();
    let mut evals = 0usize;
    let mut outer = 0usize;
    let mut prev_feas = f64::INFINITY;
    for _ in 0..max_outer {
        outer += 1;
        // Inner minimization of the augmented Lagrangian.
        let (lam, nuv, m) = (lambda.clone(), nu.clone(), mu);
        let mut inner = |xv: &[f64]| -> (f64, Vec<f64>) {
            let (f, mut g) = checked_fg(&mut *problem.fg, xv);
            let cev = checked_constraints("equality", problem.ce, xv, Some(ne));
            let civ = checked_constraints("inequality", problem.ci, xv, Some(ni));
            let mut val = f;
            // Equalities: λᵀc + (μ/2)‖c‖²; gradient (λ + μc)ᵀ∇c.
            if !cev.is_empty() {
                let w: Vec<f64> = cev
                    .iter()
                    .zip(&lam)
                    .map(|(c, l)| m.mul_add(*c, *l))
                    .collect();
                for (c, l) in cev.iter().zip(&lam) {
                    val += l * c + 0.5 * m * c * c;
                }
                let pull = checked_jt("equality", problem.ce_jt, xv, &w);
                for i in 0..g.len() {
                    g[i] += pull[i];
                }
            }
            // Inequalities (PHR): (1/2μ)·Σ [max(0, ν + μc)² − ν²].
            if !civ.is_empty() {
                let w: Vec<f64> = civ
                    .iter()
                    .zip(&nuv)
                    .map(|(c, v)| m.mul_add(*c, *v).max(0.0))
                    .collect();
                for (wi, v) in w.iter().zip(&nuv) {
                    val += (wi * wi - v * v) / (2.0 * m);
                }
                let pull = checked_jt("inequality", problem.ci_jt, xv, &w);
                for i in 0..g.len() {
                    g[i] += pull[i];
                }
            }
            assert!(
                val.is_finite(),
                "augmented-Lagrangian objective must remain finite"
            );
            assert_finite("augmented-Lagrangian gradient", &g);
            (val, g)
        };
        let mut st = LbfgsState::new(&x, 10, &mut inner);
        let rep = st.run(&mut inner, &StopRule::GradNorm(0.1 * tol), 300);
        evals += rep.evals;
        x.clone_from(&st.x);
        // Multiplier updates.
        let cev = checked_constraints("equality", problem.ce, &x, Some(ne));
        let civ = checked_constraints("inequality", problem.ci, &x, Some(ni));
        for i in 0..lambda.len() {
            lambda[i] = mu.mul_add(cev[i], lambda[i]);
        }
        for i in 0..nu.len() {
            nu[i] = mu.mul_add(civ[i], nu[i]).max(0.0);
        }
        assert_finite("equality multiplier", &lambda);
        assert_finite("inequality multiplier", &nu);
        let feas = inf_norm(&cev).max(civ.iter().map(|c| c.max(0.0)).fold(0.0f64, f64::max));
        // Penalty growth when feasibility stalls (classical schedule).
        if feas > 0.25 * prev_feas {
            mu = (mu * 10.0).min(1e10);
        }
        prev_feas = feas;
        let kkt = kkt_residual(problem, &x, &lambda, &nu);
        evals += 1;
        if kkt.within_tolerance(tol) {
            let (f, _) = checked_fg(&mut *problem.fg, &x);
            return AugLagReport {
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
    }
    let kkt = kkt_residual(problem, &x, &lambda, &nu);
    let (f, _) = checked_fg(&mut *problem.fg, &x);
    AugLagReport {
        x,
        f,
        kkt,
        lambda,
        nu,
        outer_iters: outer,
        evals,
        converged: false,
    }
}

/// Compute the KKT residuals at (x, λ, ν) — the certificate builder.
pub fn kkt_residual(
    problem: &mut ConstrainedProblem<'_>,
    x: &[f64],
    lambda: &[f64],
    nu: &[f64],
) -> KktResidual {
    validate_point(x);
    let (_, mut g) = checked_fg(&mut *problem.fg, x);
    let cev = checked_constraints("equality", problem.ce, x, None);
    let civ = checked_constraints("inequality", problem.ci, x, None);
    assert_eq!(
        lambda.len(),
        cev.len(),
        "equality multiplier length must match the equality constraint count"
    );
    assert_eq!(
        nu.len(),
        civ.len(),
        "inequality multiplier length must match the inequality constraint count"
    );
    assert_finite("equality multiplier", lambda);
    assert_finite("inequality multiplier", nu);
    validate_jt_zero("equality", problem.ce_jt, x, cev.len());
    validate_jt_zero("inequality", problem.ci_jt, x, civ.len());
    if !cev.is_empty() {
        let equality_pull = checked_jt("equality", problem.ce_jt, x, lambda);
        for i in 0..g.len() {
            g[i] += equality_pull[i];
        }
    }
    if !civ.is_empty() {
        let inequality_pull = checked_jt("inequality", problem.ci_jt, x, nu);
        for i in 0..g.len() {
            g[i] += inequality_pull[i];
        }
    }
    let feasibility = inf_norm(&cev).max(civ.iter().map(|c| c.max(0.0)).fold(0.0f64, f64::max));
    let dual_feasibility = nu
        .iter()
        .map(|multiplier| (-multiplier).max(0.0))
        .fold(0.0f64, f64::max);
    let complementarity = civ
        .iter()
        .enumerate()
        .map(|(i, c)| (c * nu[i]).abs())
        .fold(0.0f64, f64::max);
    KktResidual {
        stationarity: inf_norm(&g),
        feasibility,
        dual_feasibility,
        complementarity,
    }
}
