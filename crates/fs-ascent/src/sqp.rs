//! SQP for tightly-constrained SMALL-DIMENSION polish (bead ijil,
//! §9.2): active-set sequential quadratic programming with a damped
//! BFGS Lagrangian-Hessian approximation. Each iteration solves the
//! complete inequality-constrained QP through a dual working-set method
//! and dense KKT factorizations (fs-la LU), then globalizes the step
//! with multiplier-adaptive exact ℓ1 merit and Armijo backtracking.
//! Constraint addition and release happen inside the QP, before any
//! nonlinear objective trial.
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
    /// Objective+gradient callback calls, including validation and KKT checks.
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
    let original_rhs = rhs.clone();
    fact.solve(&mut rhs);
    // Reuse the dense factor for two residual corrections. Large physical
    // multipliers otherwise lose the much smaller feasibility step through
    // cancellation, especially when constraint units rescale the KKT rows.
    for _ in 0..2 {
        let mut correction: Vec<f64> = kkt.chunks_exact(dim)
            .zip(&original_rhs)
            .map(|(row, initial)| row.iter().zip(&rhs)
                .fold(*initial, |r, (a, x)| (-*a).mul_add(*x, r)))
            .collect();
        fact.solve(&mut correction);
        for (value, change) in rhs.iter_mut().zip(correction) {
            *value += change;
        }
    }
    if rhs.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let d = rhs[..n].to_vec();
    let mult = rhs[n..].to_vec();
    Some((d, mult))
}

/// Unweighted violation used by the exact ℓ1 merit.
fn violation(ce: &[f64], ci: &[f64]) -> f64 {
    ce.iter().map(|c| c.abs()).sum::<f64>()
        + ci.iter().map(|c| c.max(0.0)).sum::<f64>()
}

struct QpSolution {
    d: Vec<f64>,
    lambda: Vec<f64>,
    nu: Vec<f64>,
}

fn constraint_block(
    je: &[f64],
    ce: &[f64],
    ji: &[f64],
    ci: &[f64],
    active: &[usize],
    n: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut a = je.to_vec();
    let mut c = ce.to_vec();
    for &j in active {
        a.extend_from_slice(&ji[j * n..(j + 1) * n]);
        c.push(ci[j]);
    }
    (a, c)
}

/// Residual and a scale-relative rounding allowance for one linear row.
fn linear_residual(row: &[f64], d: &[f64], c: f64) -> Option<(f64, f64)> {
    let mut value = c;
    let mut scale = c.abs();
    for (a, x) in row.iter().zip(d) {
        value = a.mul_add(*x, value);
        scale += (a * x).abs();
    }
    let allowance = 64.0 * f64::EPSILON * scale.max(f64::MIN_POSITIVE);
    (value.is_finite() && scale.is_finite()).then_some((value, allowance))
}

/// Solve the complete convex QP, including currently inactive inequalities.
///
/// Start at the equality-constrained minimum (dual feasible with zero
/// inequality multipliers). A violated row raises its multiplier along a
/// projected search direction until it binds, or an old multiplier reaches
/// zero and releases that row. Thus violated/redundant inequalities are not
/// all imposed as simultaneous equalities. Dense KKT solves deliberately keep
/// this a small-problem path. Singular equality blocks, incompatible linear
/// constraints, arithmetic failure, or a bounded pivot limit return None;
/// none of these is a nonlinear infeasibility certificate.
fn solve_inequality_qp(
    b: &[f64],
    g: &[f64],
    je: &[f64],
    ce: &[f64],
    ji: &[f64],
    ci: &[f64],
) -> Option<QpSolution> {
    let (n, ne, ni) = (g.len(), ce.len(), ci.len());
    let (mut d, mut mult) = solve_qp(b, g, je, ce, n, ne)?;
    let mut active = Vec::<usize>::new();
    let dimension = n.checked_add(ne)?.checked_add(ni)?.checked_add(1)?;
    let mut pivots_left = dimension.checked_mul(dimension)?.checked_mul(8)?;
    loop {
        let mut selected = None;
        let mut worst = 0.0f64;
        for j in 0..ni {
            if active.contains(&j) {
                continue;
            }
            let (value, allowance) = linear_residual(&ji[j * n..(j + 1) * n], &d, ci[j])?;
            if value > allowance && value > worst {
                selected = Some(j);
                worst = value;
            }
        }
        let Some(selected) = selected else {
            let mut nu = vec![0.0; ni];
            for (r, &j) in active.iter().enumerate() {
                nu[j] = mult[ne + r];
            }
            return Some(QpSolution { d, lambda: mult[..ne].to_vec(), nu });
        };
        let row = &ji[selected * n..(selected + 1) * n];
        let mut pending_multiplier = 0.0f64;
        loop {
            pivots_left = pivots_left.checked_sub(1)?;
            let (a, _) = constraint_block(je, ce, ji, ci, &active, n);
            let m = ne + active.len();
            // B z + Aᵀ r = -a_new; A z = 0. Increasing the pending
            // multiplier by t moves d by t*z and the old duals by t*r.
            let (z, r) = solve_qp(b, row, &a, &vec![0.0; m], n, m)?;
            let (value, _) = linear_residual(row, &d, ci[selected])?;
            let rate = -row.iter().zip(&z).map(|(a, z)| a * z).sum::<f64>();
            if !rate.is_finite() {
                return None;
            }
            let full = if rate > 0.0 { value / rate } else { f64::INFINITY };
            let mut partial = f64::INFINITY;
            let mut blocking = None;
            for j in 0..active.len() {
                if r[ne + j] < 0.0 {
                    let bound = mult[ne + j] / -r[ne + j];
                    if bound < partial {
                        partial = bound;
                        blocking = Some(j);
                    }
                }
            }
            let t = full.min(partial);
            if !t.is_finite() || t < 0.0 {
                return None;
            }
            for (value, z) in d.iter_mut().zip(z) {
                *value = t.mul_add(z, *value);
            }
            for (j, (value, r)) in mult.iter_mut().zip(r).enumerate() {
                let updated = t.mul_add(r, *value);
                if j >= ne {
                    let roundoff = 64.0 * f64::EPSILON * (value.abs() + (t * r).abs());
                    if !updated.is_finite() || !roundoff.is_finite() || updated < -roundoff {
                        return None;
                    }
                    *value = updated.max(0.0);
                } else {
                    *value = updated;
                }
            }
            pending_multiplier += t;
            if !pending_multiplier.is_finite() || d.iter().chain(&mult).any(|x| !x.is_finite()) {
                return None;
            }
            if full <= partial {
                active.push(selected);
                // Re-solve the final block instead of retaining cancellation
                // in d + t*z when the equality-only minimizer was far away.
                let (a, c) = constraint_block(je, ce, ji, ci, &active, n);
                (d, mult) = solve_qp(b, g, &a, &c, n, ne + active.len())?;
                let scale = mult.iter().map(|v| v.abs()).fold(1.0f64, f64::max);
                for value in &mut mult[ne..] {
                    if *value < -64.0 * f64::EPSILON * scale {
                        return None;
                    }
                    *value = (*value).max(0.0);
                }
                break;
            }
            let j = blocking?;
            active.remove(j);
            mult.remove(ne + j);
        }
    }
}

struct MeritStep<'a> {
    f: f64,
    violation: f64,
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
    penalty: &mut f64,
) -> Option<Vec<f64>> {
    let n = step.x.len();
    // Exact-penalty descent requires a weight above the QP multiplier norm.
    // A fixed weight can forbid every feasibility-restoring step when the
    // physical objective or the constraint units change. Never lower a weight
    // already established by earlier iterations.
    let multiplier_norm = step.lambda.iter().chain(step.nu)
        .map(|v| v.abs()).fold(0.0f64, f64::max);
    *penalty = (*penalty).max(1.1 * multiplier_norm);
    let gd: f64 = step.g.iter().zip(step.d).map(|(g, d)| g * d).sum();
    // The QP step satisfies the linearized constraints, so this is an upper
    // bound on the one-sided directional derivative of the ℓ1 merit.
    let slope = gd - *penalty * step.violation;
    if !penalty.is_finite() || !slope.is_finite() || slope >= 0.0 {
        return None;
    }
    let mut alpha = 1.0f64;
    for _ in 0..40 {
        let xt: Vec<f64> = step
            .x
            .iter()
            .zip(step.d)
            .map(|(xi, di)| alpha.mul_add(*di, *xi))
            .collect();
        if xt.as_slice() == step.x || xt.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let (ft, gt) = checked_fg(&mut *problem.fg, &xt);
        *evals += 1;
        let cet = checked_constraints("equality", problem.ce, &xt, Some(ne));
        let cit = checked_constraints("inequality", problem.ci, &xt, Some(ni));
        // Compare changes rather than subtracting two large penalized totals.
        // Armijo scales with the proposed improvement: there is no absolute
        // 1e-12 floor that strands small-amplitude physical objectives.
        let change = (ft - step.f) + *penalty * (violation(&cet, &cit) - step.violation);
        if change.is_finite() && change < 0.0 && change <= 1e-4 * alpha * slope {
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
    let mut evals = 1usize; // validate_problem_at_start evaluates fg once.
    let mut penalty = 10.0f64;
    let mut iters = 0usize;
    let mut lambda = vec![0.0f64; ne];
    let mut nu = vec![0.0f64; ni];
    for _ in 0..max_iters {
        iters += 1;
        let (f, g) = checked_fg(&mut *problem.fg, &x);
        let cev = checked_constraints("equality", problem.ce, &x, Some(ne));
        let civ = checked_constraints("inequality", problem.ci, &x, Some(ni));
        evals += 1;
        let je = jacobian("equality", problem.ce_jt, &x, ne, n);
        let ji = jacobian("inequality", problem.ci_jt, &x, ni, n);
        let Some(qp) = solve_inequality_qp(&b, &g, &je, &cev, &ji, &civ) else {
            break; // no admissible linearized step; report the actual KKT below
        };
        let d = qp.d;
        lambda = qp.lambda;
        nu = qp.nu;
        // Convergence: small step + certificate.
        let dnorm = d.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        let kkt = kkt_residual(problem, &x, &lambda, &nu);
        evals += 1;
        if dnorm < tol && kkt.within_tolerance(tol) {
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
        let step = MeritStep {
            f,
            violation: violation(&cev, &civ),
            x: &x,
            d: &d,
            g: &g,
            lambda: &lambda,
            nu: &nu,
        };
        let Some(accepted_x) =
            accept_merit_step(problem, &step, &mut b, &mut evals, ne, ni, &mut penalty)
        else {
            break; // merit stall — certificate below tells the truth
        };
        x = accepted_x;
    }
    let (f, _) = checked_fg(&mut *problem.fg, &x);
    let kkt = kkt_residual(problem, &x, &lambda, &nu);
    evals += 2;
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

#[cfg(test)]
mod qp_tests {
    use super::solve_inequality_qp;

    #[test]
    fn inactive_inequalities_enter_the_qp_before_a_trial_is_taken() {
        let q = solve_inequality_qp(
            &[1.0, 0.0, 0.0, 1.0], &[-3.0, -4.0], &[], &[],
            &[1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, 1.0],
            &[-1.0, -2.0, 0.0, 0.0, -4.0],
        ).expect("feasible box QP");
        assert!((q.d[0] - 1.0).abs() < 1e-12);
        assert!((q.d[1] - 2.0).abs() < 1e-12);
        for (actual, expected) in q.nu.iter().zip([2.0, 2.0, 0.0, 0.0, 0.0]) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn tighter_redundant_bound_releases_the_looser_dual_constraint() {
        for (rows, rhs, expected) in [
            ([-100.0, -1.0], [100.0, 2.0], [0.0, 2.0]),
            ([-1.0, -100.0], [2.0, 100.0], [2.0, 0.0]),
        ] {
            let q = solve_inequality_qp(&[1.0], &[0.0], &[], &[], &rows, &rhs)
                .expect("compatible redundant bounds");
            assert!((q.d[0] - 2.0).abs() < 1e-12);
            for (actual, expected) in q.nu.iter().zip(expected) {
                assert!((actual - expected).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn duplicate_active_bounds_do_not_make_the_equality_block_singular() {
        let q = solve_inequality_qp(
            &[1.0, 0.0, 0.0, 1.0], &[-1.6, -0.4], &[1.0, 1.0], &[0.0],
            &[1.0, 0.0, 2.0, 0.0], &[0.0, 0.0],
        ).expect("duplicate inequality rows need only one active representative");
        assert!(q.d.iter().all(|v| v.abs() < 1e-12));
        assert!((q.lambda[0] - 0.4).abs() < 1e-12);
        assert!((q.nu[0] + 2.0 * q.nu[1] - 1.2).abs() < 1e-12);
        assert!(q.nu.iter().all(|v| *v >= 0.0));
    }

    #[test]
    fn equality_rows_are_preserved_while_an_inequality_is_added() {
        let q = solve_inequality_qp(
            &[1.0, 0.0, 0.0, 1.0], &[-3.0, -4.0], &[1.0, 1.0], &[-1.0],
            &[1.0, 0.0], &[0.5],
        ).expect("compatible equality and inequality");
        assert!((q.d[0] + 0.5).abs() < 1e-12);
        assert!((q.d[1] - 1.5).abs() < 1e-12);
        assert!((q.lambda[0] - 2.5).abs() < 1e-12);
        assert!((q.nu[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn inconsistent_linear_constraints_are_refused_instead_of_clipped() {
        assert!(solve_inequality_qp(
            &[1.0], &[0.0], &[], &[], &[1.0, -1.0], &[0.0, 1.0],
        ).is_none());
        assert!(solve_inequality_qp(
            &[1.0], &[0.0], &[], &[], &[0.0], &[1.0],
        ).is_none());
    }
}
