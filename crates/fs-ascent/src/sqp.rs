//! SQP for tightly-constrained SMALL-DIMENSION polish (bead ijil,
//! §9.2): active-set sequential quadratic programming with a damped
//! BFGS Lagrangian-Hessian approximation. Each iteration solves the
//! equality-constrained QP on the current working set through a dense
//! KKT factorization (fs-la LU — the "QP subproblems via fs-la dense
//! factors" the bead names), takes a backtracking step on the ℓ1
//! merit, and updates the working set from QP multiplier signs and
//! constraint violations.
//!
//! Scope: n and the constraint counts are SMALL (the polish regime —
//! warm starts near an optimum converge in a handful of iterations,
//! gated). Large-scale SQP (sparse KKT, trust-region globalization)
//! is recorded follow-up, not claimed.

use crate::auglag::{
    ConstrainedProblem, KktResidual, assert_finite, checked_constraints, checked_fg, checked_jt,
    kkt_residual, validate_problem_at_start, validate_tolerance,
};
use fs_la::factor::lu;

type JtAction<'a> = dyn Fn(&[f64], &[f64]) -> Vec<f64> + 'a;

/// Outcome of an SQP solve.
#[derive(Debug, Clone)]
pub struct SqpReport {
    /// Final iterate.
    pub x: Vec<f64>,
    /// Final objective.
    pub f: f64,
    /// The certificate.
    pub kkt: KktResidual,
    /// Equality multipliers.
    pub lambda: Vec<f64>,
    /// Inequality multipliers (≥ 0; zero off the working set).
    pub nu: Vec<f64>,
    /// SQP iterations.
    pub iters: usize,
    /// Objective/constraint evaluations.
    pub evals: usize,
    /// Certificate below tolerance.
    pub converged: bool,
}

/// Reconstruct a constraint Jacobian (rows = constraints) by probing
/// the Jᵀ·w action with unit vectors — fixture-scale by design.
fn jacobian(label: &str, jt: &JtAction<'_>, x: &[f64], m: usize, n: usize) -> Vec<f64> {
    let mut j = vec![0.0f64; m * n];
    for k in 0..m {
        let mut w = vec![0.0f64; m];
        w[k] = 1.0;
        let row = checked_jt(label, jt, x, &w);
        j[k * n..(k + 1) * n].copy_from_slice(&row);
    }
    j
}

/// Damped BFGS update (Powell damping keeps B positive definite even
/// when the Lagrangian curvature is indefinite near the boundary).
fn bfgs_update(b: &mut [f64], n: usize, s: &[f64], y: &[f64]) {
    let mut bs = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            bs[i] += b[i * n + j] * s[j];
        }
    }
    let sbs: f64 = s.iter().zip(&bs).map(|(a, c)| a * c).sum();
    let sy: f64 = s.iter().zip(y).map(|(a, c)| a * c).sum();
    if sbs <= 0.0 {
        return;
    }
    // Powell damping.
    let theta = if sy >= 0.2 * sbs {
        1.0
    } else {
        0.8 * sbs / (sbs - sy)
    };
    let r: Vec<f64> = y
        .iter()
        .zip(&bs)
        .map(|(yi, bsi)| theta * yi + (1.0 - theta) * bsi)
        .collect();
    let sr: f64 = s.iter().zip(&r).map(|(a, c)| a * c).sum();
    if sr.abs() < 1e-300 {
        return;
    }
    for i in 0..n {
        for j in 0..n {
            b[i * n + j] += r[i] * r[j] / sr - bs[i] * bs[j] / sbs;
        }
    }
}

/// Solve the working-set QP: min ½dᵀBd + gᵀd s.t. A d = −c (rows =
/// equalities + active inequalities) through the dense KKT system.
/// Returns (d, multipliers) or None on a singular KKT (degenerate set).
fn solve_qp(
    b: &[f64],
    g: &[f64],
    a: &[f64],
    c: &[f64],
    n: usize,
    m: usize,
) -> Option<(Vec<f64>, Vec<f64>)> {
    let dim = n + m;
    let mut kkt = vec![0.0f64; dim * dim];
    for i in 0..n {
        for j in 0..n {
            kkt[i * dim + j] = b[i * n + j];
        }
    }
    for r in 0..m {
        for j in 0..n {
            kkt[(n + r) * dim + j] = a[r * n + j];
            kkt[j * dim + n + r] = a[r * n + j];
        }
    }
    let fact = lu(&kkt, dim).ok()?;
    let mut rhs = vec![0.0f64; dim];
    for i in 0..n {
        rhs[i] = -g[i];
    }
    for r in 0..m {
        rhs[n + r] = -c[r];
    }
    fact.solve(&mut rhs);
    if rhs.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let d = rhs[..n].to_vec();
    let mult = rhs[n..].to_vec();
    Some((d, mult))
}

/// ℓ1 merit: f + w·(‖c_e‖₁ + ‖max(0, c_i)‖₁).
fn merit_of(
    problem: &mut ConstrainedProblem<'_>,
    xx: &[f64],
    w: f64,
    ev: &mut usize,
    ne: usize,
    ni: usize,
) -> f64 {
    let (fv, _) = checked_fg(&mut *problem.fg, xx);
    *ev += 1;
    let ce_v = checked_constraints("equality", problem.ce, xx, Some(ne));
    let ci_v = checked_constraints("inequality", problem.ci, xx, Some(ni));
    let viol: f64 =
        ce_v.iter().map(|c| c.abs()).sum::<f64>() + ci_v.iter().map(|c| c.max(0.0)).sum::<f64>();
    let merit = w.mul_add(viol, fv);
    assert!(merit.is_finite(), "SQP merit value must remain finite");
    merit
}

fn active_constraints(problem: &ConstrainedProblem<'_>, x: &[f64], ni: usize) -> Vec<usize> {
    let civ = checked_constraints("inequality", problem.ci, x, Some(ni));
    (0..ni).filter(|&j| civ[j] > -1e-8).collect()
}

fn working_set_block(
    problem: &ConstrainedProblem<'_>,
    x: &[f64],
    active: &[usize],
    cev: &[f64],
    civ: &[f64],
    dims: (usize, usize, usize),
) -> (Vec<f64>, Vec<f64>) {
    let (ne, ni, n) = dims;
    let je = jacobian("equality", problem.ce_jt, x, ne, n);
    let ji = jacobian("inequality", problem.ci_jt, x, ni, n);
    let m = ne + active.len();
    let mut a = vec![0.0f64; m * n];
    let mut c = vec![0.0f64; m];
    a[..ne * n].copy_from_slice(&je);
    c[..ne].copy_from_slice(cev);
    for (r, &j) in active.iter().enumerate() {
        a[(ne + r) * n..(ne + r + 1) * n].copy_from_slice(&ji[j * n..(j + 1) * n]);
        c[ne + r] = civ[j];
    }
    (a, c)
}

fn update_multipliers(
    active: &mut Vec<usize>,
    lambda: &mut [f64],
    nu: &mut [f64],
    mult: &[f64],
    ne: usize,
) -> bool {
    assert_eq!(
        mult.len(),
        ne + active.len(),
        "SQP multiplier dimension must match the working set"
    );
    assert_finite("SQP multiplier", mult);
    lambda.copy_from_slice(&mult[..ne]);
    nu.fill(0.0);
    for (r, &j) in active.iter().enumerate() {
        nu[j] = mult[ne + r].max(0.0);
    }
    let mut drop_idx = None;
    let mut worst = -1e-10f64;
    for r in 0..active.len() {
        if mult[ne + r] < worst {
            worst = mult[ne + r];
            drop_idx = Some(r);
        }
    }
    if let Some(r) = drop_idx {
        let j = active.remove(r);
        nu[j] = 0.0;
        true
    } else {
        false
    }
}

struct MeritStep<'a> {
    x: &'a [f64],
    d: &'a [f64],
    g: &'a [f64],
    lambda: &'a [f64],
    nu: &'a [f64],
}

fn accept_merit_step(
    problem: &mut ConstrainedProblem<'_>,
    step: &MeritStep<'_>,
    b: &mut [f64],
    evals: &mut usize,
    ne: usize,
    ni: usize,
) -> Option<Vec<f64>> {
    let n = step.x.len();
    let m0 = merit_of(problem, step.x, 10.0, evals, ne, ni);
    let mut alpha = 1.0f64;
    for _ in 0..40 {
        let xt: Vec<f64> = step
            .x
            .iter()
            .zip(step.d)
            .map(|(xi, di)| alpha.mul_add(*di, *xi))
            .collect();
        if merit_of(problem, &xt, 10.0, evals, ne, ni) < m0 - 1e-12 {
            let (_, gt) = checked_fg(&mut *problem.fg, &xt);
            *evals += 1;
            let pull_e = checked_jt("equality", problem.ce_jt, &xt, step.lambda);
            let pull_i = checked_jt("inequality", problem.ci_jt, &xt, step.nu);
            let pull_e0 = checked_jt("equality", problem.ce_jt, step.x, step.lambda);
            let pull_i0 = checked_jt("inequality", problem.ci_jt, step.x, step.nu);
            let s: Vec<f64> = xt.iter().zip(step.x).map(|(a2, b2)| a2 - b2).collect();
            let y: Vec<f64> = (0..n)
                .map(|i| (gt[i] + pull_e[i] + pull_i[i]) - (step.g[i] + pull_e0[i] + pull_i0[i]))
                .collect();
            assert_finite("SQP curvature vector", &y);
            bfgs_update(b, n, &s, &y);
            assert_finite("SQP Hessian approximation", b);
            return Some(xt);
        }
        alpha *= 0.5;
    }
    None
}

fn activate_violated(active: &mut Vec<usize>, civ: &[f64]) {
    for (j, &cj) in civ.iter().enumerate() {
        if cj > -1e-10 && !active.contains(&j) {
            active.push(j);
        }
    }
    active.sort_unstable();
}

/// Run active-set SQP from `x0`.
pub fn sqp(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    tol: f64,
    max_iters: usize,
) -> SqpReport {
    validate_tolerance(tol);
    let n = x0.len();
    let (ne, ni) = validate_problem_at_start(problem, x0);
    let mut x = x0.to_vec();
    let mut b = vec![0.0f64; n * n];
    for i in 0..n {
        b[i * n + i] = 1.0;
    }
    let mut evals = 0usize;
    let mut iters = 0usize;
    // Working set: active inequality indices (violated-or-near ones).
    let mut active = active_constraints(problem, &x, ni);
    let mut lambda = vec![0.0f64; ne];
    let mut nu = vec![0.0f64; ni];
    for _ in 0..max_iters {
        iters += 1;
        let (f, g) = checked_fg(&mut *problem.fg, &x);
        let cev = checked_constraints("equality", problem.ce, &x, Some(ne));
        let civ = checked_constraints("inequality", problem.ci, &x, Some(ni));
        evals += 1;
        let m = ne + active.len();
        let (a, c) = working_set_block(problem, &x, &active, &cev, &civ, (ne, ni, n));
        let Some((d, mult)) = solve_qp(&b, &g, &a, &c, n, m) else {
            break; // degenerate working set — report honestly below
        };
        let dropped = update_multipliers(&mut active, &mut lambda, &mut nu, &mult, ne);
        // Convergence: small step + certificate.
        let dnorm = d.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        let kkt = kkt_residual(problem, &x, &lambda, &nu);
        evals += 1;
        if !dropped && dnorm < tol && kkt.within_tolerance(tol) {
            return SqpReport {
                x,
                f,
                kkt,
                lambda,
                nu,
                iters,
                evals,
                converged: true,
            };
        }
        if dropped {
            continue; // re-solve with the reduced set before stepping
        }
        let step = MeritStep {
            x: &x,
            d: &d,
            g: &g,
            lambda: &lambda,
            nu: &nu,
        };
        let Some(accepted_x) = accept_merit_step(problem, &step, &mut b, &mut evals, ne, ni) else {
            break; // merit stall — certificate below tells the truth
        };
        x = accepted_x;
        // Activate violated inequalities at the new point.
        let civ_new = checked_constraints("inequality", problem.ci, &x, Some(ni));
        activate_violated(&mut active, &civ_new);
    }
    let (f, _) = checked_fg(&mut *problem.fg, &x);
    let kkt = kkt_residual(problem, &x, &lambda, &nu);
    let converged = kkt.within_tolerance(tol);
    SqpReport {
        x,
        f,
        kkt,
        lambda,
        nu,
        iters,
        evals,
        converged,
    }
}
