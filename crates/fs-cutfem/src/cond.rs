//! Dense conditioning probe for the conformance batteries: the G0
//! ghost-penalty acceptance ("conditioning independent of cut
//! fraction") is verified via EIGENVALUE estimates, not solver
//! iteration counts — Jacobi rotations on the densified operator give
//! every eigenvalue of the fixture-sized systems, so λ_min/λ_max are
//! measurements, not extrapolations. Deliberately fixture-only: the
//! size gate refuses production-scale matrices (a dense O(n³) probe on
//! a real system would be a silent performance lie).

use fs_la::eigen::jacobi_eigh;
use fs_sparse::Csr;

/// Spectral summary of an SPD stiffness matrix.
#[derive(Debug, Clone, Copy)]
pub struct CondReport {
    /// Smallest eigenvalue.
    pub lambda_min: f64,
    /// Largest eigenvalue.
    pub lambda_max: f64,
    /// `lambda_max / lambda_min` (`f64::INFINITY` when `lambda_min ≤ 0`,
    /// which for a nominally SPD assembly is itself a finding).
    pub cond: f64,
}

/// Full-spectrum conditioning of a (small) symmetric matrix.
///
/// # Panics
/// If the matrix is larger than the fixture gate (4096) — this probe
/// is for conformance batteries, not production diagnostics.
#[must_use]
pub fn condition_estimate(a: &Csr) -> CondReport {
    let n = a.nrows();
    assert!(
        n <= 4096,
        "dense conditioning probe is gated to conformance fixtures (n = {n})"
    );
    let dense = a.to_dense();
    let (vals, _) = jacobi_eigh(&dense, n);
    assert!(
        vals.iter().all(|value| value.is_finite()),
        "dense conditioning probe produced a non-finite eigenvalue"
    );
    let lambda_min = vals.iter().copied().fold(f64::INFINITY, f64::min);
    let lambda_max = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let cond = if lambda_min > 0.0 {
        lambda_max / lambda_min
    } else {
        f64::INFINITY
    };
    CondReport {
        lambda_min,
        lambda_max,
        cond,
    }
}
