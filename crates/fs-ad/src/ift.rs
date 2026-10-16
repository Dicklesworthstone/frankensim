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
//! exact to floating point, O((N+M)·cost(F)).
//!
//! Bead o3ui adds the MATRIX-FREE tangent route
//! ([`ift_gradient_matrix_free`]): ∂F/∂u is only ever APPLIED (one
//! directional-dual pass per application — seed every uᵢ with εᵢ = vᵢ),
//! and the caller supplies the linear solver (fs-solver's Krylov stack
//! at L3, or anything else — fs-ad stays L1 and solver-agnostic). This
//! is the FORWARD/tangent form: one solve per parameter, the right
//! shape for few parameters and large N. The matrix-free ADJOINT form
//! (one solve total) needs Jᵀ·v products, i.e. reverse mode over
//! vectors — recorded with the FrankenTorch bridge follow-ups.

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
    let lift =
        |xs: &[f64]| -> Vec<Dual64<1>> { xs.iter().map(|&v| Dual64::<1>::constant(v)).collect() };
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
    Ok((
        grad,
        IftReport {
            primal_residual,
            adjoint_residual,
        },
    ))
}

/// Evidence attached to a matrix-free IFT gradient.
#[derive(Debug, Clone)]
pub struct MatrixFreeIftReport {
    /// ‖F(u*, p)‖∞ — the caller's honesty check, as in [`IftReport`].
    pub primal_residual: f64,
    /// Worst ‖(∂F/∂u)·wₖ + (∂F/∂p)·eₖ‖∞ over parameters — how well the
    /// caller's solver actually solved each tangent system (measured
    /// HERE with a fresh operator application, not trusted from the
    /// solver).
    pub tangent_residual: f64,
}

/// dJ/dp at u* via the TANGENT route with a caller-supplied solver:
/// for each parameter k, solve (∂F/∂u)·wₖ = −(∂F/∂p)·eₖ matrix-free,
/// then dJ/dpₖ = ∂J/∂pₖ + (∂J/∂u)·wₖ. `solve(apply, b)` must return an
/// approximate solution of A·w = b given only `apply: v ↦ A·v`.
///
/// Every Jacobian access is ONE directional dual pass — nothing N×N is
/// ever formed. Cost: M solves + M+2 objective/residual passes.
pub fn ift_gradient_matrix_free<F, J, L>(
    f: &F,
    j: &J,
    u_star: &[f64],
    p: &[f64],
    solve: &mut L,
) -> (Vec<f64>, MatrixFreeIftReport)
where
    F: Fn(&[Dual64<1>], &[Dual64<1>]) -> Vec<Dual64<1>>,
    J: Fn(&[Dual64<1>], &[Dual64<1>]) -> Dual64<1>,
    L: FnMut(&mut dyn FnMut(&[f64]) -> Vec<f64>, &[f64]) -> Vec<f64>,
{
    let n = u_star.len();
    let m = p.len();
    let lift =
        |xs: &[f64]| -> Vec<Dual64<1>> { xs.iter().map(|&v| Dual64::<1>::constant(v)).collect() };
    let u0 = lift(u_star);
    let p0 = lift(p);
    let r0 = f(&u0, &p0);
    assert_eq!(r0.len(), n, "F must return u.len() = {n} residuals");
    let primal_residual = r0.iter().fold(0.0f64, |a, d| a.max(d.re.abs()));

    // (∂F/∂u)·v in one dual pass: seed every u component along v.
    let dfdu_apply = |v: &[f64]| -> Vec<f64> {
        let u: Vec<Dual64<1>> = u_star
            .iter()
            .zip(v)
            .map(|(&ui, &vi)| Dual64 { re: ui, eps: [vi] })
            .collect();
        f(&u, &p0).iter().map(|d| d.eps[0]).collect()
    };
    // ∂J/∂u as a vector is not needed: (∂J/∂u)·wₖ is one dual pass too.
    let djdu_dot = |w: &[f64]| -> f64 {
        let u: Vec<Dual64<1>> = u_star
            .iter()
            .zip(w)
            .map(|(&ui, &wi)| Dual64 { re: ui, eps: [wi] })
            .collect();
        j(&u, &p0).eps[0]
    };

    let mut grad = vec![0.0f64; m];
    let mut tangent_residual = 0.0f64;
    for k in 0..m {
        // −(∂F/∂p)·eₖ and ∂J/∂pₖ, one seeded pass each.
        let mut pp = lift(p);
        pp[k] = Dual64::<1>::variable(p[k], 0);
        let rhs: Vec<f64> = f(&u0, &pp).iter().map(|d| -d.eps[0]).collect();
        let djdpk = j(&u0, &pp).eps[0];
        // Caller's matrix-free solve of (∂F/∂u)·w = rhs.
        let mut apply = |v: &[f64]| dfdu_apply(v);
        let w = solve(&mut apply, &rhs);
        assert_eq!(w.len(), n, "solver must return a length-{n} solution");
        // Measured solve quality (fresh application; never trusted).
        let aw = dfdu_apply(&w);
        let res = aw
            .iter()
            .zip(&rhs)
            .fold(0.0f64, |a, (&x, &b)| a.max((x - b).abs()));
        tangent_residual = tangent_residual.max(res);
        grad[k] = djdu_dot(&w) + djdpk;
    }
    (
        grad,
        MatrixFreeIftReport {
            primal_residual,
            tangent_residual,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Real;

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
        assert!(
            rep.primal_residual < 1e-12,
            "Newton must have converged: {rep:?}"
        );
        assert!(
            rep.adjoint_residual < 1e-12,
            "adjoint solve must be tight: {rep:?}"
        );
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
            let mut ud: Vec<Dual64<1>> = vec![Dual64::constant(1.0), Dual64::constant(1.0)];
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

    /// Unpreconditioned dense-storage GMRES (test scaffolding ONLY —
    /// production callers pass fs-solver's Krylov stack; fs-ad stays
    /// solver-agnostic).
    fn gmres(apply: &mut dyn FnMut(&[f64]) -> Vec<f64>, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let m = n.min(60);
        let beta = b.iter().map(|x| x * x).sum::<f64>().sqrt();
        if beta == 0.0 {
            return vec![0.0; n];
        }
        let mut v: Vec<Vec<f64>> = vec![b.iter().map(|x| x / beta).collect()];
        let mut h = vec![vec![0.0f64; m]; m + 1];
        let mut k_used = 0;
        for k in 0..m {
            let mut w = apply(&v[k]);
            for (i, vi) in v.iter().enumerate() {
                let hik = w.iter().zip(vi).map(|(a, c)| a * c).sum::<f64>();
                h[i][k] = hik;
                for (wj, vj) in w.iter_mut().zip(vi) {
                    *wj -= hik * vj;
                }
            }
            let hk1 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            h[k + 1][k] = hk1;
            k_used = k + 1;
            if hk1 < 1e-14 {
                break;
            }
            v.push(w.iter().map(|x| x / hk1).collect());
        }
        // Solve the small least-squares min ‖βe₁ − H y‖ by normal
        // equations (fixture-scale Krylov dimension).
        let kk = k_used;
        let mut ata = vec![0.0f64; kk * kk];
        let mut atb = vec![0.0f64; kk];
        for c1 in 0..kk {
            for c2 in 0..kk {
                let mut acc = 0.0;
                for row in h.iter().take(kk + 1) {
                    acc += row[c1] * row[c2];
                }
                ata[c1 * kk + c2] = acc;
            }
            atb[c1] = h[0][c1] * beta;
        }
        let fact = lu(&ata, kk).expect("GMRES normal equations nonsingular");
        let mut y = atb;
        fact.solve(&mut y);
        let mut x = vec![0.0f64; n];
        for (k, yk) in y.iter().enumerate() {
            for (xi, vi) in x.iter_mut().zip(&v[k]) {
                *xi += yk * vi;
            }
        }
        x
    }

    #[test]
    fn matrix_free_matches_dense_ift() {
        const N: usize = 50;
        // Same 2-system: the tangent route with a caller-supplied
        // Krylov solve must reproduce the dense adjoint gradient.
        let p = [3.0, 2.0];
        let u = newton_solve(&f_sys, &[1.0, 1.0], &p, 60);
        let (dense, _) = ift_gradient(&f_sys, &j_obj, &u, &p).unwrap();
        let mut solver = |apply: &mut dyn FnMut(&[f64]) -> Vec<f64>, b: &[f64]| gmres(apply, b);
        let (mf, rep) = ift_gradient_matrix_free(&f_sys, &j_obj, &u, &p, &mut solver);
        assert!(rep.primal_residual < 1e-12, "converged primal: {rep:?}");
        assert!(rep.tangent_residual < 1e-10, "tight solves: {rep:?}");
        for k in 0..2 {
            assert!(
                (mf[k] - dense[k]).abs() < 1e-9 * dense[k].abs().max(1.0),
                "matrix-free grad[{k}] = {} vs dense {}",
                mf[k],
                dense[k]
            );
        }
        // Large-N shape check: nonlinear tridiagonal system, N = 50,
        // M = 2 — the regime the tangent route exists for. Nothing N×N
        // is formed inside (by construction of the API).
        let f_big = |u: &[Dual64<1>], p: &[Dual64<1>]| -> Vec<Dual64<1>> {
            let mut r = Vec::with_capacity(N);
            for i in 0..N {
                let left = if i > 0 {
                    u[i - 1]
                } else {
                    Dual64::constant(0.0)
                };
                let right = if i + 1 < N {
                    u[i + 1]
                } else {
                    Dual64::constant(0.0)
                };
                let three = Dual64::constant(3.0);
                r.push(three * u[i] - left - right + u[i].powi(3) - p[0] - p[1] * left);
            }
            r
        };
        let j_big = |u: &[Dual64<1>], _p: &[Dual64<1>]| -> Dual64<1> {
            let mut acc = Dual64::constant(0.0);
            for &ui in u {
                acc = acc + ui * ui;
            }
            acc
        };
        let pb = [1.0, 0.3];
        let ub = newton_solve(&f_big, &vec![0.5; N], &pb, 40);
        let (dense_b, _) = ift_gradient(&f_big, &j_big, &ub, &pb).unwrap();
        let (mf_b, rep_b) = ift_gradient_matrix_free(&f_big, &j_big, &ub, &pb, &mut solver);
        assert!(rep_b.primal_residual < 1e-11, "big primal: {rep_b:?}");
        for k in 0..2 {
            assert!(
                (mf_b[k] - dense_b[k]).abs() < 1e-7 * dense_b[k].abs().max(1.0),
                "N=50 matrix-free grad[{k}] = {} vs dense {}",
                mf_b[k],
                dense_b[k]
            );
        }
        println!(
            "{{\"suite\":\"fs-ad\",\"case\":\"ift-matrix-free\",\"verdict\":\"pass\",\"detail\":\"tangent route == dense adjoint (2-system to 1e-9, N=50 tridiag to 1e-7); worst tangent residual {:.2e}\"}}",
            rep_b.tangent_residual.max(rep.tangent_residual)
        );
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
