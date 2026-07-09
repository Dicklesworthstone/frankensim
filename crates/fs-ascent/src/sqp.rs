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

use crate::auglag::{ConstrainedProblem, KktResidual, kkt_residual};
use fs_la::factor::lu;

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
fn jacobian(
    jt: &dyn Fn(&[f64], &[f64]) -> Vec<f64>,
    x: &[f64],
    m: usize,
    n: usize,
) -> Vec<f64> {
    let mut j = vec![0.0f64; m * n];
    for k in 0..m {
        let mut w = vec![0.0f64; m];
        w[k] = 1.0;
        let row = jt(x, &w);
        j[k * n..(k + 1) * n].copy_from_slice(&row[..n]);
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
fn solve_qp(b: &[f64], g: &[f64], a: &[f64], c: &[f64], n: usize, m: usize) -> Option<(Vec<f64>, Vec<f64>)> {
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
    let d = rhs[..n].to_vec();
    let mult = rhs[n..].to_vec();
    Some((d, mult))
}

/// ℓ1 merit: f + w·(‖c_e‖₁ + ‖max(0, c_i)‖₁).
fn merit_of(problem: &mut ConstrainedProblem<'_>, xx: &[f64], w: f64, ev: &mut usize) -> f64 {
    let (fv, _) = (problem.fg)(xx);
    *ev += 1;
    let ce_v = (problem.ce)(xx);
    let ci_v = (problem.ci)(xx);
    let viol: f64 = ce_v.iter().map(|c| c.abs()).sum::<f64>()
        + ci_v.iter().map(|c| c.max(0.0)).sum::<f64>();
    w.mul_add(viol, fv)
}

/// Run active-set SQP from `x0`.
pub fn sqp(
    problem: &mut ConstrainedProblem<'_>,
    x0: &[f64],
    tol: f64,
    max_iters: usize,
) -> SqpReport {
    let n = x0.len();
    let ne = (problem.ce)(x0).len();
    let ni = (problem.ci)(x0).len();
    let mut x = x0.to_vec();
    let mut b = vec![0.0f64; n * n];
    for i in 0..n {
        b[i * n + i] = 1.0;
    }
    let mut evals = 0usize;
    let mut iters = 0usize;
    // Working set: active inequality indices (violated-or-near ones).
    let mut active: Vec<usize> = {
        let civ = (problem.ci)(&x);
        (0..ni).filter(|&j| civ[j] > -1e-8).collect()
    };
    let mut lambda = vec![0.0f64; ne];
    let mut nu = vec![0.0f64; ni];
    let merit_w = 10.0f64;
    for _ in 0..max_iters {
        iters += 1;
        let (f, g) = (problem.fg)(&x);
        let cev = (problem.ce)(&x);
        let civ = (problem.ci)(&x);
        evals += 1;
        // Assemble the working-set constraint block.
        let je = jacobian(problem.ce_jt, &x, ne, n);
        let ji = jacobian(problem.ci_jt, &x, ni, n);
        let m = ne + active.len();
        let mut a = vec![0.0f64; m * n];
        let mut c = vec![0.0f64; m];
        a[..ne * n].copy_from_slice(&je);
        c[..ne].copy_from_slice(&cev);
        for (r, &j) in active.iter().enumerate() {
            a[(ne + r) * n..(ne + r + 1) * n].copy_from_slice(&ji[j * n..(j + 1) * n]);
            c[ne + r] = civ[j];
        }
        let Some((d, mult)) = solve_qp(&b, &g, &a, &c, n, m) else {
            break; // degenerate working set — report honestly below
        };
        lambda.copy_from_slice(&mult[..ne]);
        nu.fill(0.0);
        for (r, &j) in active.iter().enumerate() {
            nu[j] = mult[ne + r];
        }
        // Working-set update: drop the most negative multiplier, add
        // violated inactive constraints.
        let mut drop_idx: Option<usize> = None;
        let mut worst = -1e-10f64;
        for r in 0..active.len() {
            if mult[ne + r] < worst {
                worst = mult[ne + r];
                drop_idx = Some(r);
            }
        }
        let mut dropped = false;
        if let Some(r) = drop_idx {
            let j = active.remove(r);
            nu[j] = 0.0;
            dropped = true;
        }
        // Convergence: small step + certificate.
        let dnorm = d.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        let kkt = kkt_residual(problem, &x, &lambda, &nu);
        evals += 1;
        if !dropped
            && dnorm < tol
            && kkt.stationarity < tol
            && kkt.feasibility < tol
            && kkt.complementarity < tol
        {
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
        // ℓ1-merit backtracking line search.
        let m0 = merit_of(problem, &x, merit_w, &mut evals);
        let mut alpha = 1.0f64;
        let mut accepted = false;
        for _ in 0..40 {
            let xt: Vec<f64> = x.iter().zip(&d).map(|(xi, di)| alpha.mul_add(*di, *xi)).collect();
            if merit_of(problem, &xt, merit_w, &mut evals) < m0 - 1e-12 {
                // BFGS curvature pair on the Lagrangian gradient.
                let (_, gt) = (problem.fg)(&xt);
                evals += 1;
                let pull_e = (problem.ce_jt)(&xt, &lambda);
                let pull_i = (problem.ci_jt)(&xt, &nu);
                let pull_e0 = (problem.ce_jt)(&x, &lambda);
                let pull_i0 = (problem.ci_jt)(&x, &nu);
                let s: Vec<f64> = xt.iter().zip(&x).map(|(a2, b2)| a2 - b2).collect();
                let y: Vec<f64> = (0..n)
                    .map(|i| (gt[i] + pull_e[i] + pull_i[i]) - (g[i] + pull_e0[i] + pull_i0[i]))
                    .collect();
                bfgs_update(&mut b, n, &s, &y);
                x = xt;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break; // merit stall — certificate below tells the truth
        }
        // Activate violated inequalities at the new point.
        let civ_new = (problem.ci)(&x);
        for (j, &cj) in civ_new.iter().enumerate() {
            if cj > -1e-10 && !active.contains(&j) {
                active.push(j);
            }
        }
        active.sort_unstable();
    }
    let (f, _) = (problem.fg)(&x);
    let kkt = kkt_residual(problem, &x, &lambda, &nu);
    let converged = kkt.stationarity < tol && kkt.feasibility < tol && kkt.complementarity < tol;
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
