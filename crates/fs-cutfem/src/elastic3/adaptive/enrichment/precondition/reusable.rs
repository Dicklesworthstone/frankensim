//! Reuse a geometric correction space while density changes between solves.
//! Only the sparse interpolation is retained: every preparation recomputes the
//! diagonal and Galerkin factor from the CURRENT fine operator. There is no
//! stale-factor reuse, geometry rebuild, or fallback to a different solver.
use super::*;
use fs_sparse::Csr;

/// Size/work limits for a reusable two-level correction space. Numerical setup
/// limits apply to each prepared density, not the whole optimization campaign.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveSolveOptions3 {
    /// Total scalar terms while constructing the geometric interpolation.
    pub max_transfer_terms: usize,
    /// Accumulation terms in each exact constrained-diagonal preparation.
    pub max_diagonal_contributions: usize,
    /// Fine/coarse/vector-transfer sizes and Galerkin setup applications.
    pub two_level: TwoLevelBudget,
}
impl Default for AdaptiveSolveOptions3 {
    fn default() -> Self {
        Self { max_transfer_terms: 2_000_000, max_diagonal_contributions: 100_000_000,
            two_level: TwoLevelBudget::default() }
    }
}

/// An adaptive operator and its inseparable geometry-owned correction space.
///
/// The operator cannot be replaced or mutably borrowed from this value. Only
/// stiffness scales can change; consequently the retained interpolation cannot
/// accidentally be reused on another mesh with the same number of unknowns.
/// The coarse operator is needed only during construction, not kept alive.
/// The existing two-level owner clones/transposes the retained sparse matrix
/// for each numerical preparation; this does not recompute geometric transfer.
pub struct AdaptiveSolveSpace3 {
    operator: AdaptiveElasticity3,
    interpolation: Option<Csr>,
    options: AdaptiveSolveOptions3,
}
impl AdaptiveSolveSpace3 {
    /// Opt into exact constrained Jacobi for repeated density evaluations.
    /// Numerical diagonal admission happens at preparation, not construction.
    #[must_use]
    pub fn jacobi(operator: AdaptiveElasticity3, max_diagonal_contributions: usize) -> Self {
        Self { operator, interpolation: None,
            options: AdaptiveSolveOptions3 { max_diagonal_contributions, ..Default::default() } }
    }

    /// Bind one geometric correction space to the owned fine operator.
    /// This constructs interpolation only: no fine stiffness applications or
    /// coarse factorization are performed until `prepare` is called.
    /// The same box/material/support/refinement admission as AdaptiveTransfer3
    /// applies. Coarse stiffness scales do not enter the correction matrix.
    pub fn two_level(operator: AdaptiveElasticity3, coarse: &AdaptiveElasticity3,
        options: AdaptiveSolveOptions3, mut checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<Self, AdaptivePreconditionError3> {
        let transfer = AdaptiveTransfer3::new(coarse, &operator, options.max_transfer_terms, &mut checkpoint)
            .map_err(AdaptivePreconditionError3::Physics)?;
        let interpolation = transfer.vector_prolongation(options.two_level, &mut checkpoint)?;
        // End the temporary transfer's borrow before moving the fine operator.
        drop(transfer);
        poll(&mut checkpoint).map_err(AdaptivePreconditionError3::Physics)?;
        Ok(Self { operator, interpolation: Some(interpolation), options })
    }

    /// Actual geometry, loads, residuals and density contractions. Read-only
    /// access permits observation without invalidating the correction space.
    #[must_use]
    pub fn elasticity(&self) -> &AdaptiveElasticity3 { &self.operator }

    /// Explicitly discard the correction space and recover its physical model.
    #[must_use]
    pub fn into_elasticity(self) -> AdaptiveElasticity3 { self.operator }

    /// Transactional density change through the original material admission.
    /// Prepared values borrow this object, so scales cannot change while any
    /// numeric preconditioner is live. No update occurs on invalid input.
    pub fn set_scales(&mut self, scales: &[f64]) -> Result<(), ElasticityError3> {
        self.operator.set_scales(scales)
    }

    /// Dimension of the retained geometric correction (zero for Jacobi).
    #[must_use]
    pub fn coarse_dofs(&self) -> usize { self.interpolation.as_ref().map_or(0, Csr::ncols) }

    /// Scalar entries retained in the vector interpolation (zero for Jacobi).
    #[must_use]
    pub fn transfer_entries(&self) -> usize { self.interpolation.as_ref().map_or(0, Csr::nnz) }

    /// Build numeric factors for the current density exactly once, then share
    /// them across its independent loads. Call again after changing density.
    /// All completed Galerkin applications are exposed to the callback even
    /// when it cancels. No partial diagonal/factor is returned on failure.
    pub fn prepare(&self, mut checkpoint: impl FnMut(TwoLevelWork) -> ControlFlow<()>)
        -> Result<AdaptivePrepared3<'_>, AdaptivePreconditionError3> {
        let diagonal = self.operator.prepare_jacobi(self.options.max_diagonal_contributions,
            || checkpoint(TwoLevelWork::default())).map_err(AdaptivePreconditionError3::Physics)?;
        match &self.interpolation {
            None => Ok(AdaptivePrepared3::Jacobi(diagonal)),
            Some(p) => AdditiveTwoLevel::new(&self.operator, diagonal.inverse_diagonal(), p.clone(),
                self.options.two_level, checkpoint)
                .map(AdaptivePrepared3::TwoLevel).map_err(AdaptivePreconditionError3::Coarse),
        }
    }
}
impl LinearOp for AdaptiveSolveSpace3 {
    fn n(&self) -> usize { self.operator.n() }
    fn apply(&self, x: &[f64], y: &mut [f64]) { self.operator.apply(x, y); }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) { self.operator.apply_transpose(x, y); }
}

/// One immutable-density numerical preparation, reusable for any RHS family.
/// Both variants retain the physical operator's lifetime; neither permits a
/// stiffness mutation while its numerical data are live.
pub enum AdaptivePrepared3<'a> {
    /// Exact diagonal of the current constrained stiffness.
    Jacobi(AdaptiveJacobi3<'a>),
    /// Fixed linear SPD additive coarse correction for the current stiffness.
    TwoLevel(AdditiveTwoLevel<'a, AdaptiveElasticity3>),
}
impl AdaptivePrepared3<'_> {
    /// Completed setup applications; diagonal-only setup uses no fine applies.
    #[must_use]
    pub fn work(&self) -> TwoLevelWork {
        match self { Self::Jacobi(_) => TwoLevelWork::default(), Self::TwoLevel(p) => p.work() }
    }
    /// The exact physical operator used to prepare this numerical action.
    #[must_use]
    pub fn operator(&self) -> &AdaptiveElasticity3 {
        match self { Self::Jacobi(p) => p.operator(), Self::TwoLevel(p) => p.operator() }
    }
}
impl Precond for AdaptivePrepared3<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        match self { Self::Jacobi(p) => p.apply(r, z), Self::TwoLevel(p) => p.apply(r, z) }
    }
}
