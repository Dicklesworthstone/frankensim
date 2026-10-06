//! Right-preconditioned BiCGStab (van der Vorst 1992) with fs-sparse ILU(0).
//!
//! The recurrence residual only proposes convergence; acceptance requires the
//! recomputed true residual `||b - A x|| / ||b||` to meet the tolerance. A
//! recurrence/true mismatch restarts the method from the true residual, and
//! the shared iteration budget bounds the total work. All reductions are
//! sequential in index order, so the iterate is deterministic.

use fs_exec::CancelGate;
use fs_sparse::Csr;
use fs_sparse::precond::{Precond, ilu0};

use super::{ChtError, poll};

/// Outcome of one accepted solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct KrylovOutcome {
    pub iterations: usize,
    pub relative_residual: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        s = x.mul_add(*y, s);
    }
    s
}

fn norm(a: &[f64]) -> f64 {
    fs_math::det::sqrt(dot(a, a))
}

fn residual(a: &Csr, b: &[f64], x: &[f64], r: &mut [f64]) {
    a.spmv(x, r);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri = bi - *ri;
    }
}

/// Solve `A x = b` from the supplied initial `x`.
pub(crate) fn bicgstab_ilu0(
    system: &'static str,
    a: &Csr,
    b: &[f64],
    x: &mut [f64],
    tolerance: f64,
    max_iterations: usize,
    gate: &CancelGate,
) -> Result<KrylovOutcome, ChtError> {
    let n = b.len();
    debug_assert_eq!(a.nrows(), n);
    let b_norm = norm(b);
    if b_norm == 0.0 {
        x.fill(0.0);
        return Ok(KrylovOutcome {
            iterations: 0,
            relative_residual: 0.0,
        });
    }
    let m = ilu0(a).map_err(|e| ChtError::PreconditionerBreakdown { system, row: e.row })?;
    let mut r = vec![0.0f64; n];
    residual(a, b, x, &mut r);
    let mut rel = norm(&r) / b_norm;
    let mut iterations = 0usize;
    let mut p = vec![0.0f64; n];
    let mut v = vec![0.0f64; n];
    let mut p_hat = vec![0.0f64; n];
    let mut s_hat = vec![0.0f64; n];
    let mut t = vec![0.0f64; n];
    'restart: while rel > tolerance && iterations < max_iterations {
        let r_tilde = r.clone();
        let (mut rho, mut alpha, mut omega) = (1.0f64, 1.0f64, 1.0f64);
        p.fill(0.0);
        v.fill(0.0);
        while iterations < max_iterations {
            poll(gate)?;
            iterations += 1;
            let rho_next = dot(&r_tilde, &r);
            if rho_next == 0.0 || !rho_next.is_finite() {
                residual(a, b, x, &mut r);
                rel = norm(&r) / b_norm;
                continue 'restart;
            }
            let beta = (rho_next / rho) * (alpha / omega);
            rho = rho_next;
            for i in 0..n {
                p[i] = beta.mul_add(omega.mul_add(-v[i], p[i]), r[i]);
            }
            m.apply(&p, &mut p_hat);
            a.spmv(&p_hat, &mut v);
            let denom = dot(&r_tilde, &v);
            if denom == 0.0 || !denom.is_finite() {
                residual(a, b, x, &mut r);
                rel = norm(&r) / b_norm;
                continue 'restart;
            }
            alpha = rho / denom;
            // s overwrites r.
            for i in 0..n {
                r[i] = (-alpha).mul_add(v[i], r[i]);
            }
            if norm(&r) / b_norm <= tolerance {
                for i in 0..n {
                    x[i] = alpha.mul_add(p_hat[i], x[i]);
                }
                residual(a, b, x, &mut r);
                rel = norm(&r) / b_norm;
                continue 'restart;
            }
            m.apply(&r, &mut s_hat);
            a.spmv(&s_hat, &mut t);
            let tt = dot(&t, &t);
            if tt == 0.0 || !tt.is_finite() {
                for i in 0..n {
                    x[i] = alpha.mul_add(p_hat[i], x[i]);
                }
                residual(a, b, x, &mut r);
                rel = norm(&r) / b_norm;
                continue 'restart;
            }
            omega = dot(&t, &r) / tt;
            for i in 0..n {
                x[i] = omega.mul_add(s_hat[i], alpha.mul_add(p_hat[i], x[i]));
                r[i] = (-omega).mul_add(t[i], r[i]);
            }
            let estimate = norm(&r) / b_norm;
            if !estimate.is_finite() {
                return Err(ChtError::SolverNotConverged {
                    system,
                    iterations,
                    relative_residual: estimate,
                    tolerance,
                });
            }
            if estimate <= tolerance || omega == 0.0 {
                residual(a, b, x, &mut r);
                rel = norm(&r) / b_norm;
                continue 'restart;
            }
        }
        residual(a, b, x, &mut r);
        rel = norm(&r) / b_norm;
    }
    if rel <= tolerance {
        Ok(KrylovOutcome {
            iterations,
            relative_residual: rel,
        })
    } else {
        Err(ChtError::SolverNotConverged {
            system,
            iterations,
            relative_residual: rel,
            tolerance,
        })
    }
}
