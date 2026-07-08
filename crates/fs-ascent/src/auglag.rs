//! Augmented Lagrangian — the robust constrained default: minimize
//! f(x) s.t. c_e(x) = 0, c_i(x) ≤ 0 by an outer multiplier loop over
//! L-BFGS inner solves of
//! L_μ(x) = f + λᵀc_e + (μ/2)‖c_e‖² + (μ/2)‖max(0, c_i + s/μ)‖²-style
//! terms (the standard PHR form for inequalities). Every returned
//! optimum carries a KKT-RESIDUAL CERTIFICATE — stationarity,
//! feasibility, complementarity — so "converged" and "stalled" are
//! distinguishable outcomes, not vibes.

use crate::lbfgs::LbfgsState;
use crate::stop::StopRule;

/// The KKT residuals of a returned point (the certificate).
#[derive(Debug, Clone)]
pub struct KktResidual {
    /// ‖∇f + Σλ∇c_e + Σν∇c_i‖∞ (stationarity of the Lagrangian).
    pub stationarity: f64,
    /// max(‖c_e‖∞, ‖max(0, c_i)‖∞) (feasibility).
    pub feasibility: f64,
    /// max |ν_j · c_i_j| (complementary slackness).
    pub complementarity: f64,
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
    /// All three KKT residuals below the tolerance.
    pub converged: bool,
}

/// Problem callbacks: objective+gradient, equality constraints and
/// their Jacobian-transpose action, inequalities likewise.
#[allow(clippy::type_complexity)]
pub struct ConstrainedProblem<'a> {
    /// (f, ∇f).
    pub fg: crate::FnGrad<'a>,
    /// c_e(x) (empty vec for none).
    pub ce: &'a dyn Fn(&[f64]) -> Vec<f64>,
    /// (∂c_e/∂x)ᵀ·w.
    pub ce_jt: &'a dyn Fn(&[f64], &[f64]) -> Vec<f64>,
    /// c_i(x) (≤ 0 feasible; empty vec for none).
    pub ci: &'a dyn Fn(&[f64]) -> Vec<f64>,
    /// (∂c_i/∂x)ᵀ·w.
    pub ci_jt: &'a dyn Fn(&[f64], &[f64]) -> Vec<f64>,
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0.0f64, f64::max)
}

/// Run the PHR augmented-Lagrangian loop from `x0`.
pub fn augmented_lagrangian(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    tol: f64,
    max_outer: usize,
) -> AugLagReport {
    let ne = (problem.ce)(x0).len();
    let ni = (problem.ci)(x0).len();
    let mut lambda = vec![0.0f64; ne];
    let mut nu = vec![0.0f64; ni];
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
            let (f, mut g) = (problem.fg)(xv);
            let cev = (problem.ce)(xv);
            let civ = (problem.ci)(xv);
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
                let pull = (problem.ce_jt)(xv, &w);
                for (gi, pi) in g.iter_mut().zip(&pull) {
                    *gi += pi;
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
                let pull = (problem.ci_jt)(xv, &w);
                for (gi, pi) in g.iter_mut().zip(&pull) {
                    *gi += pi;
                }
            }
            (val, g)
        };
        let mut st = LbfgsState::new(&x, 10, &mut inner);
        let rep = st.run(&mut inner, &StopRule::GradNorm(0.1 * tol), 300);
        evals += rep.evals;
        x.clone_from(&st.x);
        // Multiplier updates.
        let cev = (problem.ce)(&x);
        let civ = (problem.ci)(&x);
        for (l, c) in lambda.iter_mut().zip(&cev) {
            *l = mu.mul_add(*c, *l);
        }
        for (v, c) in nu.iter_mut().zip(&civ) {
            *v = mu.mul_add(*c, *v).max(0.0);
        }
        let feas = inf_norm(&cev).max(civ.iter().map(|c| c.max(0.0)).fold(0.0f64, f64::max));
        // Penalty growth when feasibility stalls (classical schedule).
        if feas > 0.25 * prev_feas {
            mu = (mu * 10.0).min(1e10);
        }
        prev_feas = feas;
        let kkt = kkt_residual(problem, &x, &lambda, &nu);
        evals += 1;
        if kkt.stationarity < tol && kkt.feasibility < tol && kkt.complementarity < tol {
            let (f, _) = (problem.fg)(&x);
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
    let (f, _) = (problem.fg)(&x);
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
    let (_, mut g) = (problem.fg)(x);
    let cev = (problem.ce)(x);
    let civ = (problem.ci)(x);
    if !cev.is_empty() {
        let pull = (problem.ce_jt)(x, lambda);
        for (gi, pi) in g.iter_mut().zip(&pull) {
            *gi += pi;
        }
    }
    if !civ.is_empty() {
        let pull = (problem.ci_jt)(x, nu);
        for (gi, pi) in g.iter_mut().zip(&pull) {
            *gi += pi;
        }
    }
    let feasibility = inf_norm(&cev).max(civ.iter().map(|c| c.max(0.0)).fold(0.0f64, f64::max));
    let complementarity = civ
        .iter()
        .zip(nu)
        .map(|(c, v)| (c * v).abs())
        .fold(0.0f64, f64::max);
    KktResidual {
        stationarity: inf_norm(&g),
        feasibility,
        complementarity,
    }
}
