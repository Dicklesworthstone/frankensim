//! Density-aware diagonal and geometric two-level preparation on constrained
//! adaptive elasticity. Coarse correction uses P^T A_fine P, not a separately
//! integrated coarse stiffness (cut quadrature and ghost penalties differ).
use super::*;
use fs_solver::op::two_level::{AdditiveTwoLevel, TwoLevelBudget, TwoLevelError, TwoLevelWork};
use fs_sparse::{Coo, Csr, precond::Precond};

mod reusable;
pub use reusable::{AdaptivePrepared3, AdaptiveSolveOptions3, AdaptiveSolveSpace3};

/// Exact inverse diagonal of T^T K T, prepared for one immutable density state.
/// Includes off-diagonal physical-node terms mapping to the same master and
/// the fully condensed ghost jumps. No fine-matrix or coordinate probing.
pub struct AdaptiveJacobi3<'a> {
    operator: &'a AdaptiveElasticity3,
    inverse: Vec<f64>,
    contributions: usize,
}
impl AdaptiveJacobi3<'_> {
    /// Immutable fine operator used for setup.
    #[must_use] pub fn operator(&self) -> &AdaptiveElasticity3 { self.operator }
    /// Positive inverse diagonal, including identity entries on constrained rows.
    #[must_use] pub fn inverse_diagonal(&self) -> &[f64] { &self.inverse }
    /// Number of accumulated diagonal/scatter summands admitted by the cap.
    #[must_use] pub const fn contributions(&self) -> usize { self.contributions }
}
impl Precond for AdaptiveJacobi3<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        assert_eq!(r.len(), self.inverse.len()); assert_eq!(z.len(), r.len());
        for ((z, r), d) in z.iter_mut().zip(r).zip(&self.inverse) { *z = r * d; }
    }
}
fn charge(count: &mut usize, additional: usize, cap: usize) -> Result<(), ElasticityError3> {
    *count = count.checked_add(additional)
        .filter(|n| *n <= cap).ok_or(ElasticityError3::Invalid("diagonal contribution budget exhausted"))?;
    Ok(())
}
impl AdaptiveElasticity3 {
    /// Prepare in cell/face order using actual current stiffness scales. The
    /// contribution cap bounds retained numerical accumulation work; merge-
    /// comparisons and allocator metadata are not counted as contributions.
    /// Poll once per element-node pair, face and output-node block. A cancelled
    /// or exhausted preparation does not publish a partial diagonal.
    pub fn prepare_jacobi(&self, max_contributions: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<AdaptiveJacobi3<'_>, ElasticityError3> {
        poll(&mut checkpoint)?;
        let mut diagonal = vec![0.0; self.n()]; let mut count = 0;
        for (id, cell) in self.raw.cells.iter().enumerate() {
            for a in 0..8 { for b in 0..8 {
                poll(&mut checkpoint)?;
                let (ra, rb) = (&self.rows[cell.nodes[a]], &self.rows[cell.nodes[b]]);
                let (mut i, mut j) = (0, 0);
                // Rows have sorted, unique master indices. T never couples
                // displacement components, so only equal-component K entries
                // contribute to this diagonal, but a and b need not be equal.
                while i < ra.len() && j < rb.len() {
                    if ra[i].0 < rb[j].0 { i += 1; }
                    else if rb[j].0 < ra[i].0 { j += 1; }
                    else {
                        let master = ra[i].0;
                        if !self.fixed[master] {
                            charge(&mut count, 3, max_contributions)?;
                            let factor = self.scales()[id] * ra[i].1 * rb[j].1;
                            for c in 0..3 { diagonal[3*master+c] += factor * cell.stiffness[3*a+c][3*b+c]; }
                        }
                        i += 1; j += 1;
                    }
                }
            } }
        }
        for face in &self.raw.ghosts {
            poll(&mut checkpoint)?;
            let mut jump = BTreeMap::new();
            for (&node, &derivative) in face.nodes.iter().zip(&face.jump) {
                for &(master, weight) in &self.rows[node] {
                    if !self.fixed[master] {
                        charge(&mut count, 1, max_contributions)?;
                        *jump.entry(master).or_insert(0.0) += derivative * weight;
                    }
                }
            }
            let weight = face.weight * 0.5 * (self.scales()[face.cells[0]] + self.scales()[face.cells[1]]);
            for (master, value) in jump {
                charge(&mut count, 3, max_contributions)?;
                for c in 0..3 { diagonal[3*master+c] += weight * value * value; }
            }
        }
        for (i, value) in diagonal.iter_mut().enumerate() {
            if i % 192 == 0 { poll(&mut checkpoint)?; }
            if self.fixed[i/3] { *value = 1.0; }
            if !value.is_finite() || *value <= 0.0 {
                return Err(ElasticityError3::Invalid("nonpositive or nonfinite constrained diagonal"));
            }
            *value = 1.0 / *value;
            if !value.is_finite() || *value <= 0.0 {
                return Err(ElasticityError3::Invalid("constrained inverse diagonal overflow"));
            }
        }
        poll(&mut checkpoint)?;
        Ok(AdaptiveJacobi3 { operator: self, inverse: diagonal, contributions: count })
    }
}

/// Distinguishes a failed physical preparation from generic Galerkin setup.
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptivePreconditionError3 {
    /// Geometry-owned diagonal preparation refused.
    Physics(ElasticityError3),
    /// Size/work/cancellation or coarse-factorization refusal.
    Coarse(TwoLevelError),
}
impl std::fmt::Display for AdaptivePreconditionError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "adaptive preconditioner refused: {self:?}") }
}
impl std::error::Error for AdaptivePreconditionError3 {}

impl<'a> AdaptiveTransfer3<'a> {
    // One vectorization/admission path for transient enrichment and retained
    // optimization spaces. No stiffness or numerical factor is retained here.
    fn vector_prolongation(&self, budget: TwoLevelBudget,
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Csr, AdaptivePreconditionError3> {
        let coarse_error = AdaptivePreconditionError3::Coarse;
        let mut poll_setup = || {
            if checkpoint().is_break() { Err(coarse_error(TwoLevelError::Cancelled)) } else { Ok(()) }
        };
        poll_setup()?;
        let nc = 3 * self.coarse.fixed.iter().filter(|&&fixed| !fixed).count();
        if nc == 0 { return Err(coarse_error(TwoLevelError::Invalid("empty unconstrained coarse space"))); }
        if self.fine.n() > budget.max_fine_dofs || nc > budget.max_coarse_dofs || nc > 512
            || nc > budget.max_operator_applications {
            return Err(coarse_error(TwoLevelError::Budget("adaptive coarse space/setup")));
        }
        let entries = self.rows.iter().try_fold(0usize, |n, row| n.checked_add(row.len().checked_mul(3)?))
            .ok_or_else(|| coarse_error(TwoLevelError::Budget("interpolation size overflow")))?;
        if entries > budget.max_transfer_entries {
            return Err(coarse_error(TwoLevelError::Budget("vector interpolation entries")));
        }
        let mut columns = vec![None; self.coarse.nodes.len()]; let mut next = 0;
        for (node, &fixed) in self.coarse.fixed.iter().enumerate() {
            if !fixed { columns[node] = Some(next); next += 3; }
        }
        let mut coo = Coo::new(self.fine.n(), nc);
        for (node, row) in self.rows.iter().enumerate() {
            poll_setup()?;
            for &(master, weight) in row {
                let column = columns[master].ok_or_else(|| coarse_error(TwoLevelError::Invalid("fixed transfer column")))?;
                for c in 0..3 { coo.push(3*node+c, column+c, weight); }
            }
        }
        let p = coo.assemble();
        poll_setup()?;
        Ok(p)
    }

    /// Prepare the fixed-SPD two-level action against the FINE density state.
    /// Remove fixed coarse coordinates instead of adding artificial coarse
    /// pivots. Reuse the exact Q1 transfer, its transpose, fs-solver's Galerkin
    /// builder and fs-la's direct factorization. The coarse geometry supplies
    /// a correction space; its separately integrated stiffness is NOT used.
    ///
    /// A setup-budget failure does not silently fall back or change the model.
    /// Callers can explicitly request diagonal-only preparation instead. Coarse
    /// setup applications are reported separately from subsequent Krylov work.
    pub fn prepare_two_level(&self, budget: TwoLevelBudget, max_diagonal_contributions: usize,
        mut checkpoint: impl FnMut(TwoLevelWork) -> ControlFlow<()>)
        -> Result<AdditiveTwoLevel<'a, AdaptiveElasticity3>, AdaptivePreconditionError3> {
        let p = self.vector_prolongation(budget, || checkpoint(TwoLevelWork::default()))?;
        let jacobi = self.fine.prepare_jacobi(max_diagonal_contributions,
            || checkpoint(TwoLevelWork::default())).map_err(AdaptivePreconditionError3::Physics)?;
        AdditiveTwoLevel::new(self.fine, jacobi.inverse_diagonal(), p, budget, checkpoint)
            .map_err(AdaptivePreconditionError3::Coarse)
    }
}
