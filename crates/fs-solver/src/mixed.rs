//! Mixed-precision Krylov (tfz.10 slice 3): f32 INNER CG under f64
//! iterative refinement — the Krylov analogue of fs-la's
//! f32-factor/f64-refine ladder, with the same evidence-object
//! philosophy. The inner solver runs entirely in f32 (storage AND
//! arithmetic: half the memory traffic — the win on bandwidth-starved
//! machines); the outer loop computes TRUE residuals in f64 and
//! corrects, so the achieved accuracy is f64-grade whenever
//! κ(A) ≪ 1/ε_f32 (the classic refinement condition; the report says
//! what was achieved rather than assuming).

use crate::op::LinearOp;
use crate::{dot, norm2};
use fs_sparse::Csr;

/// An f32 copy of an assembled operator (values truncated once;
/// structure shared conceptually with the f64 original).
pub struct CsrF32 {
    nrows: usize,
    row_ptr: Vec<usize>,
    cols: Vec<usize>,
    vals: Vec<f32>,
}

impl CsrF32 {
    /// Truncate an f64 CSR to f32 storage.
    #[must_use]
    pub fn from_csr(a: &Csr) -> CsrF32 {
        let mut row_ptr = Vec::with_capacity(a.nrows() + 1);
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        row_ptr.push(0usize);
        for r in 0..a.nrows() {
            let (rc, rv) = a.row(r);
            for (&c, &v) in rc.iter().zip(rv) {
                cols.push(c);
                #[allow(clippy::cast_possible_truncation)]
                vals.push(v as f32);
            }
            row_ptr.push(cols.len());
        }
        CsrF32 {
            nrows: a.nrows(),
            row_ptr,
            cols,
            vals,
        }
    }

    /// y = A·x in f32 arithmetic (fixed ascending order — the same
    /// determinism discipline as the f64 path).
    pub fn spmv(&self, x: &[f32], y: &mut [f32]) {
        for (r, yr) in y.iter_mut().enumerate().take(self.nrows) {
            let mut acc = 0.0f32;
            for i in self.row_ptr[r]..self.row_ptr[r + 1] {
                acc = self.vals[i].mul_add(x[self.cols[i]], acc);
            }
            *yr = acc;
        }
    }
}

/// Evidence object for a mixed-precision solve.
#[derive(Debug, Clone)]
pub struct MixedReport {
    /// Outer refinement steps taken.
    pub refine_steps: usize,
    /// Total INNER f32 CG iterations across all refinements.
    pub inner_iters: usize,
    /// Achieved relative residual (f64, true).
    pub rel_residual: f64,
    /// Target met.
    pub converged: bool,
    /// True residual after each refinement step.
    pub history: Vec<f64>,
    /// True when the inner f32 solve stalled and the driver should
    /// escalate to the plain f64 path (κ too large for f32 inner
    /// work — reported, not silently absorbed).
    pub escalate: bool,
}

/// Inner CG entirely in f32 (storage + arithmetic), fixed-order
/// reductions, run to a LOOSE inner tolerance.
fn cg_f32(a: &CsrF32, b: &[f32], tol: f32, max_iters: usize) -> (Vec<f32>, usize) {
    let n = b.len();
    let mut x = vec![0.0f32; n];
    let mut r = b.to_vec();
    let mut p = b.to_vec();
    let mut ap = vec![0.0f32; n];
    let dot32 = |a: &[f32], b: &[f32]| -> f32 {
        let mut acc = 0.0f32;
        for (x, y) in a.iter().zip(b) {
            acc = x.mul_add(*y, acc);
        }
        acc
    };
    let bnorm = dot32(b, b).sqrt().max(f32::MIN_POSITIVE);
    let mut rz = dot32(&r, &r);
    let mut iters = 0usize;
    for _ in 0..max_iters {
        if rz.sqrt() / bnorm < tol {
            break;
        }
        a.spmv(&p, &mut ap);
        let pap = dot32(&p, &ap);
        if !(pap.is_finite() && pap > 0.0) {
            break; // f32 breakdown: report upward, never mask.
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] = alpha.mul_add(p[i], x[i]);
            r[i] = alpha.mul_add(-ap[i], r[i]);
        }
        let rz_new = dot32(&r, &r);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = beta.mul_add(p[i], r[i]);
        }
        iters += 1;
    }
    (x, iters)
}

/// Solve A·x = b to f64-grade accuracy with f32 inner CG under f64
/// iterative refinement. `a64` supplies TRUE residuals; `a32` is its
/// truncation. Inner solves run to 1e-4 relative (past that, f32
/// arithmetic yields nothing).
pub fn mixed_cg_refine(
    a64: &dyn LinearOp,
    a32: &CsrF32,
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_refines: usize,
    max_inner: usize,
) -> MixedReport {
    let n = b.len();
    let bnorm = norm2(b).max(f64::MIN_POSITIVE);
    let mut history = Vec::new();
    let mut inner_total = 0usize;
    let mut ax = vec![0.0f64; n];
    let mut escalate = false;
    for step in 0..=max_refines {
        // True f64 residual.
        a64.apply(x, &mut ax);
        let r64: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
        let rel = norm2(&r64) / bnorm;
        history.push(rel);
        if rel < tol {
            return MixedReport {
                refine_steps: step,
                inner_iters: inner_total,
                rel_residual: rel,
                converged: true,
                history,
                escalate: false,
            };
        }
        if step == max_refines {
            break;
        }
        // Stall detection: refinement must contract per step.
        if history.len() >= 2 && rel > history[history.len() - 2] * 0.5 {
            escalate = true;
            break;
        }
        // Inner f32 correction solve on the SCALED residual (scaling
        // keeps the f32 dynamic range healthy).
        let rnorm = norm2(&r64).max(f64::MIN_POSITIVE);
        #[allow(clippy::cast_possible_truncation)]
        let r32: Vec<f32> = r64.iter().map(|&v| (v / rnorm) as f32).collect();
        let (c32, it) = cg_f32(a32, &r32, 1e-4, max_inner);
        inner_total += it;
        for i in 0..n {
            x[i] = rnorm.mul_add(f64::from(c32[i]), x[i]);
        }
    }
    let rel = *history.last().expect("at least one residual");
    MixedReport {
        refine_steps: history.len() - 1,
        inner_iters: inner_total,
        rel_residual: rel,
        converged: false,
        history,
        escalate,
    }
}

/// Deterministic dot re-export for battery cross-checks.
#[must_use]
pub fn dot64(a: &[f64], b: &[f64]) -> f64 {
    dot(a, b)
}
