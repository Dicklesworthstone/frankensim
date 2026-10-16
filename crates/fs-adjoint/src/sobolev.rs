//! Sobolev (H¹) gradient smoothing — the Riesz-representation step
//! (plan §8.7's "single most important practical trick in shape
//! optimization"). A raw L² gradient of a mesh-discretized functional
//! amplifies mesh noise; re-representing it in the H¹ inner product
//! ⟨g, v⟩_{H¹} = ⟨g_raw, v⟩_{L²} means solving
//! (M + α·K)·g_smooth = M·g_raw — one SPD solve through fs-solver.
//! The metric choice (α) IS a preconditioner choice; α ~ h² leaves
//! physical features intact while killing grid-frequency noise.

use fs_solver::{CgState, CsrOp};
use fs_sparse::Csr;
use fs_sparse::precond::IdentityPrecond;

/// Smooth a raw nodal gradient through the H¹ Riesz problem:
/// (M + α·K)·g = M·g_raw. Returns (g_smooth, cg_iterations).
///
/// # Panics
/// If the smoothing solve fails (SPD system at fixture scale — a
/// failure is a shape bug, not a tolerance issue).
#[must_use]
pub fn sobolev_smooth(
    mass: &Csr,
    stiffness: &Csr,
    alpha: f64,
    g_raw: &[f64],
    tol: f64,
) -> (Vec<f64>, usize) {
    let n = g_raw.len();
    assert_eq!(mass.nrows(), n, "mass dimension mismatch");
    assert_eq!(stiffness.nrows(), n, "stiffness dimension mismatch");
    // A = M + α·K (deterministic COO merge).
    let mut coo = fs_sparse::Coo::new(n, n);
    for r in 0..n {
        let (cols, vals) = mass.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            coo.push(r, c, v);
        }
        let (cols, vals) = stiffness.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            coo.push(r, c, alpha * v);
        }
    }
    let a = CsrOp::symmetric(coo.assemble());
    // RHS = M·g_raw.
    let mut b = vec![0.0f64; n];
    mass.spmv(g_raw, &mut b);
    let mut st = CgState::new(&a, &IdentityPrecond, &b);
    let rep = st.run(&a, &IdentityPrecond, tol, 10_000);
    assert!(rep.converged, "Sobolev smoothing solve failed: {rep:?}");
    (st.x, rep.iters)
}
