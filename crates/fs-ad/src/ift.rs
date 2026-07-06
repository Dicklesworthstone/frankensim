//! Implicit-function-theorem adjoints (plan §6.6 regime 2): differentiate
//! THROUGH A SOLUTION instead of through the iteration that found it.
//!
//! Given a residual F(u, p) = 0 defining u(p) implicitly, and a scalar
//! objective J(u, p):
//!
//!   dJ/dp = ∂J/∂p − λᵀ·∂F/∂p,   where  (∂F/∂u)ᵀ·λ = ∂J/∂u.
//!
//! One adjoint solve replaces differentiating every solver iteration —
//! and, crucially, the gradient does not degrade when the solver stops
//! early (the report carries the primal residual so callers KNOW how well
//! u* actually solves F).
//!
//! v1 builds the Jacobians DENSELY, column by column, with single-lane
//! forward duals: deterministic (fixed seeding order, fused arithmetic),
//! exact to floating point, O((N+M)·cost(F)). Matrix-free adjoints for
//! large N join the solver-stack work; the FrankenTorch reverse bridge is
//! the recorded follow-up.

use crate::Real;
use crate::dual::Dual64;
use fs_la::factor::{FactorError, lu};

/// Evidence attached to an IFT gradient: how good the primal was, and how
/// well the adjoint system was solved.
#[derive(Debug, Clone)]
pub struct IftReport {
    /// ‖F(u*, p)‖∞ — how much of a "solution" u* really is. The gradient
    /// formula is exact only at F = 0; this is the caller's honesty check.
    pub primal_residual: f64,
    /// ‖(∂F/∂u)ᵀλ − ∂J/∂u‖∞ achieved by the adjoint solve.
    pub adjoint_residual: f64,
}

/// dJ/dp at an (approximate) solution u* of F(u, p) = 0.
///
/// `f` maps (u, p) → residual vector (length = u.len()), generic over
/// [`Real`]; `j` maps (u, p) → scalar objective. Both must be pure.
///
/// # Errors
/// [`FactorError`] if ∂F/∂u is singular (the IFT hypothesis fails — no
/// implicit function exists at this point).
pub fn ift_gradient<F, J>(
    f: &F,
    j: &J,
    u_star: &[f64],
    p: &[f64],
) -> Result<(Vec<f64>, IftReport), FactorError>
where
    F: Fn(&[Dual64<1>], &[Dual64<1>]) -> Vec<Dual64<1>>,
    J: Fn(&[Dual64<1>], &[Dual64<1>]) -> Dual64<1>,
{
    let n = u_star.len();
    let m = p.len();
    let lift = |xs: &[f64]| -> Vec<Dual64<1>> {
        xs.iter().map(|&v| Dual64::<1>::constant(v)).collect()
    };
    // Primal residual (dual constants; primal channel is bit-exact).
    let u0 = lift(u_star);
    let p0 = lift(p);
    let r0 = f(&u0, &p0);
    assert_eq!(r0.len(), n, "F must return u.len() = {n} residuals");
    let primal_residual = r0.iter().fold(0.0f64, |a, d| a.max(d.re.abs()));

    // ∂F/∂u column by column (seed u_k), row-major N×N.
    let mut dfdu = vec![0.0f64; n * n];
    for k in 0..n {
        let mut u = lift(u_star);
        u[k] = Dual64::<1>::variable(u_star[k], 0);
        let r = f(&u, &p0);
        for (i, d) in r.iter().enumerate() {
            dfdu[i * n + k] = d.eps[0];
        }
    }
    // ∂F/∂p column by column, N×M.
    let mut dfdp = vec![0.0f64; n * m];
    for k in 0..m {
        let mut pp = lift(p);
        pp[k] = Dual64::<1>::variable(p[k], 0);
        let r = f(&u0, &pp);
        for (i, d) in r.iter().enumerate() {
            dfdp[i * m + k] = d.eps[0];
        }
    }
    // ∂J/∂u and ∂J/∂p.
    let mut djdu = vec![0.0f64; n];
    for (k, slot) in djdu.iter_mut().enumerate() {
        let mut u = lift(u_star);
        u[k] = Dual64::<1>::variable(u_star[k], 0);
        *slot = j(&u, &p0).eps[0];
    }
    let mut djdp = vec![0.0f64; m];
    for (k, slot) in djdp.iter_mut().enumerate() {
        let mut pp = lift(p);
        pp[k] = Dual64::<1>::variable(p[k], 0);
        *slot = j(&u0, &pp).eps[0];
    }

    // Adjoint solve: (∂F/∂u)ᵀ λ = ∂J/∂u.
    let fact = lu(&dfdu, n)?;
    let mut lambda = djdu.clone();
    fact.solve_transpose(&mut lambda);
    // Adjoint residual for the report.
    let mut adjoint_residual = 0.0f64;
    for k in 0..n {
        let mut acc = -djdu[k];
        for (i, &li) in lambda.iter().enumerate() {
            acc = dfdu[i * n + k].mul_add(li, acc);
        }
        adjoint_residual = adjoint_residual.max(acc.abs());
    }
    // dJ/dp = ∂J/∂p − λᵀ·∂F/∂p.
    let mut grad = djdp;
    for (jj, g) in grad.iter_mut().enumerate() {
        let mut acc = 0.0f64;
        for (i, &li) in lambda.iter().enumerate() {
            acc = li.mul_add(dfdp[i * m + jj], acc);
        }
        *g -= acc;
    }
    Ok((grad, IftReport { primal_residual, adjoint_residual }))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Solve F(u, p) = 0 by Newton on f64 (test scaffolding).
    fn newton_solve<F>(f: &F, u0: &[f64], p: &[f64], iters: usize) -> Vec<f64>
    where
        F: Fn(&[Dual64<1>], &[Dual64<1>]) -> Vec<Dual64<1>>,
    {
        let n = u0.len();
        let mut u = u0.to_vec();
        for _ in 0..iters {
            let lift_p: Vec<Dual64<1>> = p.iter().map(|&v| Dual64::constant(v)).collect();
            let mut jac = vec![0.0f64; n * n];
            let mut r0 = vec![0.0f64; n];
            for k in 0..n {
                let mut ud: Vec<Dual64<1>> = u.iter().map(|&v| Dual64::constant(v)).collect();
                ud[k] = Dual64::variable(u[k], 0);
                let r = f(&ud, &lift_p);
                for (i, d) in r.iter().enumerate() {
                    jac[i * n + k] = d.eps[0];
                    if k == 0 {
                        r0[i] = d.re;
                    }
                }
            }
            let fact = lu(&jac, n).expect("Newton Jacobian nonsingular");
            let mut step = r0;
            fact.solve(&mut step);
            for (ui, si) in u.iter_mut().zip(&step) {
                *ui -= si;
            }
        }
        u
    }

    /// Nonlinear 2-system: F1 = u0³ + u1 − p0, F2 = u0 + u1³ − p1.
    fn f_sys(u: &[Dual64<1>], p: &[Dual64<1>]) -> Vec<Dual64<1>> {
        vec![u[0].powi(3) + u[1] - p[0], u[0] + u[1].powi(3) - p[1]]
    }

    /// Objective J = u0² + 2·u1 + 0.5·p0.
    fn j_obj(u: &[Dual64<1>], p: &[Dual64<1>]) -> Dual64<1> {
        u[0] * u[0] + Dual64::constant(2.0) * u[1] + Dual64::constant(0.5) * p[0]
    }

    #[test]
    fn ift_matches_finite_differences() {
        let p = [3.0, 2.0];
        let u = newton_solve(&f_sys, &[1.0, 1.0], &p, 40);
        let (grad, rep) = ift_gradient(&f_sys, &j_obj, &u, &p).unwrap();
        assert!(rep.primal_residual < 1e-12, "Newton must have converged: {rep:?}");
        assert!(rep.adjoint_residual < 1e-12, "adjoint solve must be tight: {rep:?}");
        // Central differences through the FULL pipeline (re-solve at p±h).
        let h = 1e-6;
        for k in 0..2 {
            let (mut pp, mut pm) = (p, p);
            pp[k] += h;
            pm[k] -= h;
            let up = newton_solve(&f_sys, &[1.0, 1.0], &pp, 40);
            let um = newton_solve(&f_sys, &[1.0, 1.0], &pm, 40);
            let jp = j_obj(
                &up.iter().map(|&v| Dual64::constant(v)).collect::<Vec<_>>(),
                &pp.iter().map(|&v| Dual64::constant(v)).collect::<Vec<_>>(),
            )
            .re;
            let jm = j_obj(
                &um.iter().map(|&v| Dual64::constant(v)).collect::<Vec<_>>(),
                &pm.iter().map(|&v| Dual64::constant(v)).collect::<Vec<_>>(),
            )
            .re;
            let fd = (jp - jm) / (2.0 * h);
            assert!(
                (grad[k] - fd).abs() < 1e-7 * grad[k].abs().max(1.0),
                "IFT grad[{k}] = {} vs FD {fd}",
                grad[k]
            );
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"ift\",\"verdict\":\"pass\",\"detail\":\"2-system adjoint vs re-solve FD, grad=[{:.6}, {:.6}]\"}}",
            grad[0], grad[1]
        );
    }

    #[test]
    fn ift_matches_differentiate_through_iteration() {
        // The cross-check the plan demands: adjoint-of-solution must equal
        // forward-differentiating the CONVERGED iteration (fs-ad's generic
        // Newton pattern), for each parameter direction.
        let p = [3.0, 2.0];
        let u = newton_solve(&f_sys, &[1.0, 1.0], &p, 60);
        let (grad, _) = ift_gradient(&f_sys, &j_obj, &u, &p).unwrap();
        for k in 0..2 {
            // Newton generic over Real, seeded in direction e_k of p.
            let mut pd: Vec<Dual64<1>> = p.iter().map(|&v| Dual64::constant(v)).collect();
            pd[k] = Dual64::variable(p[k], 0);
            let mut ud: Vec<Dual64<1>> =
                vec![Dual64::constant(1.0), Dual64::constant(1.0)];
            for _ in 0..60 {
                // One Newton step in dual arithmetic (2×2 solved in closed
                // form to stay generic).
                let r = f_sys(&ud, &pd);
                // Dual-valued Jacobian via nested seeding is overkill for a
                // 2×2; use analytic partials instead (exact):
                let three = Dual64::constant(3.0);
                let j00 = three * ud[0] * ud[0];
                let j01 = Dual64::constant(1.0);
                let j10 = Dual64::constant(1.0);
                let j11 = three * ud[1] * ud[1];
                let det = j00 * j11 - j01 * j10;
                let d0 = (r[0] * j11 - r[1] * j01) / det;
                let d1 = (j00 * r[1] - j10 * r[0]) / det;
                ud[0] = ud[0] - d0;
                ud[1] = ud[1] - d1;
            }
            let jd = j_obj(&ud, &pd);
            assert!(
                (grad[k] - jd.eps[0]).abs() < 1e-10 * grad[k].abs().max(1.0),
                "IFT {} vs through-iteration {} at k={k}",
                grad[k],
                jd.eps[0]
            );
        }
    }

    #[test]
    fn singular_jacobian_is_typed() {
        // F(u, p) = (u0 − u1, u0 − u1): ∂F/∂u singular everywhere.
        let f = |u: &[Dual64<1>], _p: &[Dual64<1>]| vec![u[0] - u[1], u[0] - u[1]];
        let j = |u: &[Dual64<1>], _p: &[Dual64<1>]| u[0];
        let r = ift_gradient(&f, &j, &[1.0, 1.0], &[0.0]);
        assert!(r.is_err(), "singular dF/du must surface as FactorError");
    }
}
