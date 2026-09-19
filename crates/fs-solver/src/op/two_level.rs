//! Fixed SPD additive two-level preconditioning for a symmetric fine operator.
//!
//! B = D^{-1} + P (P^T A P)^{-1} P^T. The positive diagonal action covers all
//! fine coordinates, including those absent from the coarse space. A bounded
//! coarse Cholesky solve is linear, unlike a residual-stopped inner Krylov solve.
//! This is a two-level building block, not a mesh-independent multigrid claim.

use super::LinearOp;
use fs_la::factor::{Cholesky, cholesky};
use fs_sparse::{Csr, ops::transpose, precond::Precond};
use std::ops::ControlFlow;

/// Setup limits. A single sparse apply/transpose and the small factorization
/// are not interruptible; callbacks bracket these bounded operations.
#[derive(Debug, Clone, Copy)]
pub struct TwoLevelBudget {
    /// Maximum dimension of the fine operator.
    pub max_fine_dofs: usize,
    /// Maximum scalar interpolation coefficients.
    pub max_transfer_entries: usize,
    /// Maximum coarse coordinates. Also hard-capped at 512.
    pub max_coarse_dofs: usize,
    /// One fine operator application is required per coarse coordinate.
    pub max_operator_applications: usize,
}
impl Default for TwoLevelBudget {
    fn default() -> Self {
        Self { max_fine_dofs: 250_000, max_transfer_entries: 2_000_000,
            max_coarse_dofs: 384, max_operator_applications: 384 }
    }
}

/// Setup progress remains observable after refusal or cancellation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TwoLevelWork {
    /// Completed fine applications used to form the Galerkin matrix.
    pub operator_applications: usize,
}

/// No partially prepared preconditioner is returned on an error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TwoLevelError {
    /// A shape, positivity or finite-arithmetic condition failed.
    Invalid(&'static str),
    /// A caller-provided size/work allowance was insufficient.
    Budget(&'static str),
    /// The callback requested a stop.
    Cancelled,
    /// The computed coarse action was not symmetric within roundoff allowance.
    Nonsymmetric,
    /// P^T A P could not be factored as positive definite; no diagonal shift.
    NotPositiveDefinite,
}
impl std::fmt::Display for TwoLevelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "two-level setup refused: {self:?}")
    }
}
impl std::error::Error for TwoLevelError {}

/// Prepared against an immutably borrowed fine operator. No fine matrix is
/// materialized. P and its literal transpose are retained, and only the bounded
/// Galerkin matrix is factored. Callers must not mutate A through side channels.
pub struct AdditiveTwoLevel<'a, A: LinearOp> {
    operator: &'a A,
    inverse_diagonal: Vec<f64>,
    p: Csr,
    pt: Csr,
    coarse: Cholesky,
    equilibration: Vec<f64>,
    work: TwoLevelWork,
}

impl<'a, A: LinearOp> AdditiveTwoLevel<'a, A> {
    /// Build using an explicit positive diagonal inverse and sparse full-column-
    /// rank interpolation. Remove constrained coarse coordinates before calling.
    /// The caller owns A's symmetry/SPD contract; the coarse symmetry and factor
    /// checks do not prove positivity of all fine-space modes.
    ///
    /// Setup requires exactly n_coarse fine applications, not n_fine probing.
    /// Equilibrate the coarse matrix before the existing fs-la Cholesky; accept
    /// only small rounding asymmetry, average those pairs explicitly, and refuse
    /// nonpositive pivots without regularization. Allocator metadata and peak
    /// resident memory are not certified by these dimension/entry bounds.
    pub fn new(operator: &'a A, inverse_diagonal: &[f64], p: Csr,
        budget: TwoLevelBudget, mut checkpoint: impl FnMut(TwoLevelWork) -> ControlFlow<()>)
        -> Result<Self, TwoLevelError> {
        let mut work = TwoLevelWork::default();
        poll(&mut checkpoint, work)?;
        let (n, nc) = (operator.n(), p.ncols());
        if n == 0 || nc == 0 || nc > n || p.nrows() != n || inverse_diagonal.len() != n {
            return Err(TwoLevelError::Invalid("incompatible fine/diagonal/interpolation dimensions"));
        }
        if n > budget.max_fine_dofs || nc > budget.max_coarse_dofs || nc > 512
            || p.nnz() > budget.max_transfer_entries {
            return Err(TwoLevelError::Budget("fine/coarse/transfer size"));
        }
        if nc > budget.max_operator_applications {
            return Err(TwoLevelError::Budget("Galerkin operator applications"));
        }
        if !inverse_diagonal.iter().all(|d| d.is_finite() && *d > 0.0) {
            return Err(TwoLevelError::Invalid("diagonal inverse must be finite and strictly positive"));
        }
        for i in 0..n {
            if i % 256 == 0 { poll(&mut checkpoint, work)?; }
            let (columns, values) = p.row(i);
            if !values.iter().all(|v| v.is_finite()) || columns.iter().any(|&j| j >= nc)
                || columns.windows(2).any(|w| w[0] >= w[1]) {
                return Err(TwoLevelError::Invalid("invalid interpolation entry"));
            }
        }
        let pt = transpose(&p);
        poll(&mut checkpoint, work)?;
        let mut matrix = vec![0.0; nc * nc];
        let mut basis = vec![0.0; n];
        let mut applied = vec![0.0; n];
        let mut column = vec![0.0; nc];
        for j in 0..nc {
            poll(&mut checkpoint, work)?;
            basis.fill(0.0);
            let (indices, values) = pt.row(j);
            for (&i, &value) in indices.iter().zip(values) { basis[i] = value; }
            operator.apply(&basis, &mut applied);
            work.operator_applications += 1;
            poll(&mut checkpoint, work)?;
            if !applied.iter().all(|v| v.is_finite()) {
                return Err(TwoLevelError::Invalid("nonfinite fine operator action"));
            }
            pt.spmv(&applied, &mut column);
            for i in 0..nc { matrix[i * nc + j] = column[i]; }
        }
        let mut roots = Vec::with_capacity(nc);
        for i in 0..nc {
            let d = matrix[i * nc + i];
            if !d.is_finite() || d <= 0.0 { return Err(TwoLevelError::NotPositiveDefinite); }
            roots.push(d.sqrt());
        }
        for i in 0..nc {
            poll(&mut checkpoint, work)?;
            for j in 0..=i {
                let x = (matrix[i * nc + j] / roots[i]) / roots[j];
                let y = (matrix[j * nc + i] / roots[j]) / roots[i];
                if !x.is_finite() || !y.is_finite() {
                    return Err(TwoLevelError::Invalid("coarse equilibration overflow"));
                }
                if (x - y).abs() > 1e-10 * x.abs().max(y.abs()).max(1.0) {
                    return Err(TwoLevelError::Nonsymmetric);
                }
                let value = f64::midpoint(x, y);
                matrix[i * nc + j] = value;
                matrix[j * nc + i] = value;
            }
        }
        poll(&mut checkpoint, work)?;
        let coarse = cholesky(&matrix, nc).map_err(|_| TwoLevelError::NotPositiveDefinite)?;
        let equilibration: Vec<f64> = roots.iter().map(|r| 1.0 / r).collect();
        if !equilibration.iter().all(|v| v.is_finite() && *v > 0.0) {
            return Err(TwoLevelError::Invalid("coarse scale overflow"));
        }
        poll(&mut checkpoint, work)?;
        Ok(Self { operator, inverse_diagonal: inverse_diagonal.to_vec(), p, pt,
            coarse, equilibration, work })
    }
    /// Completed setup applications (separate from outer Krylov work).
    #[must_use] pub const fn work(&self) -> TwoLevelWork { self.work }
    /// Number of directly solved coarse coordinates.
    #[must_use] pub fn coarse_dofs(&self) -> usize { self.p.ncols() }
    /// The source operator cannot change stiffness while this value is alive.
    #[must_use] pub const fn operator(&self) -> &A { self.operator }
}

impl<A: LinearOp> Precond for AdditiveTwoLevel<'_, A> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        assert_eq!(r.len(), self.operator.n());
        assert_eq!(z.len(), r.len());
        let mut rhs = vec![0.0; self.p.ncols()];
        self.pt.spmv(r, &mut rhs);
        for (value, scale) in rhs.iter_mut().zip(&self.equilibration) { *value *= scale; }
        self.coarse.solve(&mut rhs);
        for (value, scale) in rhs.iter_mut().zip(&self.equilibration) { *value *= scale; }
        self.p.spmv(&rhs, z);
        for ((value, &r), &d) in z.iter_mut().zip(r).zip(&self.inverse_diagonal) { *value += d * r; }
    }
}
fn poll(checkpoint: &mut impl FnMut(TwoLevelWork) -> ControlFlow<()>, work: TwoLevelWork)
    -> Result<(), TwoLevelError> {
    if checkpoint(work).is_break() { Err(TwoLevelError::Cancelled) } else { Ok(()) }
}
