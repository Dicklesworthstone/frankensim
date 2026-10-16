//! MMS hooks: manufactured-solution studies THROUGH the generated
//! operator (gauntlet-g1's food). The study assembles the load from
//! the analytic forcing via the de-Rham map and mass weighting,
//! Dirichlet-pins the boundary, solves with the materialized
//! operator, and reports M0-weighted L2 errors plus observed orders
//! on the Kuhn refinement ladder.

use crate::plan::LoweredOperator;
use fs_feec::{deram0, element_geometry, kuhn_cube, mass_matrix, on_unit_cube_boundary};
use fs_sparse::precond::{IdentityPrecond, pcg};

/// Outcome of one MMS study.
#[derive(Debug, Clone)]
pub struct MmsReport {
    /// Mesh parameters of the ladder.
    pub ns: Vec<usize>,
    /// M0-weighted L2 errors per level.
    pub errors: Vec<f64>,
    /// Observed orders between consecutive levels.
    pub orders: Vec<f64>,
}

/// Run a Poisson-family MMS study through an operator BUILDER (called
/// per refinement level so the operator is genuinely regenerated from
/// the IR each time): −Δu = f with homogeneous Dirichlet data on the
/// unit cube. The operator must be symmetric positive definite on the
/// interior (PCG is the solver).
///
/// # Panics
/// If PCG fails to converge (fixture-scale meshes; a failure is an
/// operator bug, not a tolerance issue).
#[must_use]
pub fn mms_poisson_study<F>(
    ns: &[usize],
    build: F,
    u_exact: &dyn Fn([f64; 3]) -> f64,
    f_exact: &dyn Fn([f64; 3]) -> f64,
) -> MmsReport
where
    F: Fn(&fs_rep_mesh::TetComplex, &fs_feec::ElementGeometry) -> Vec<f64>,
{
    let mut errors = Vec::new();
    for &n in ns {
        let (complex, positions) = kuhn_cube(n);
        let geo = element_geometry(&complex, &positions);
        let m0 = mass_matrix(&complex, &geo, 0);
        // Load b = M0 · R0(f).
        let r0f = deram0(&positions, &|p| f_exact(p));
        let mut b = vec![0.0f64; r0f.len()];
        m0.spmv(&r0f, &mut b);
        // The generated operator, materialized by the caller.
        let a_full = build(&complex, &geo);
        let nv = positions.len();
        assert_eq!(
            a_full.len(),
            nv * nv,
            "build must return a dense row-major operator"
        );
        // Interior reduction.
        let interior: Vec<usize> = (0..nv)
            .filter(|&v| !on_unit_cube_boundary(positions[v]))
            .collect();
        let ni = interior.len();
        let mut slot = vec![usize::MAX; nv];
        for (i, &v) in interior.iter().enumerate() {
            slot[v] = i;
        }
        let mut red = fs_sparse::Coo::new(ni, ni);
        for (i, &v) in interior.iter().enumerate() {
            for (j, &w) in interior.iter().enumerate() {
                let val = a_full[v * nv + w];
                if val != 0.0 {
                    red.push(i, j, val);
                }
            }
        }
        let a = red.assemble();
        let rhs: Vec<f64> = interior.iter().map(|&v| b[v]).collect();
        let mut x = vec![0.0f64; ni];
        let report = pcg(&a, &rhs, &mut x, &IdentityPrecond, 1e-12, 20_000);
        assert!(report.converged, "MMS PCG failed at n={n}: {report:?}");
        // M0-weighted L2 error (boundary values exact).
        let mut e = vec![0.0f64; nv];
        for (i, &v) in interior.iter().enumerate() {
            e[v] = x[i] - u_exact(positions[v]);
        }
        let mut me = vec![0.0f64; nv];
        m0.spmv(&e, &mut me);
        let l2 = e.iter().zip(&me).map(|(a, b)| a * b).sum::<f64>().sqrt();
        errors.push(l2);
    }
    let orders = errors
        .windows(2)
        .zip(ns.windows(2))
        .map(|(e, n)| {
            let h_ratio = n[1] as f64 / n[0] as f64;
            (e[0] / e[1]).ln() / h_ratio.ln()
        })
        .collect();
    MmsReport {
        ns: ns.to_vec(),
        errors,
        orders,
    }
}

/// Convenience: materialize a lowered operator to a dense row-major
/// matrix (fixture scale).
#[must_use]
pub fn materialize_dense(op: &LoweredOperator<'_>) -> Vec<f64> {
    let csr = op.materialize().expect("linear operator required");
    csr.to_dense()
}
