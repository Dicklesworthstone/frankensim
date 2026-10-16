//! Reuse geometric correction spaces while density changes between solves.
//! Sparse interpolation is retained; diagonals and Galerkin operators/factors
//! are recomputed from the CURRENT fine density. No stale-factor reuse.
use super::*;
use super::multilevel::{AdaptiveHierarchy3, AdaptiveMultilevelOptions3, AdaptiveSetupWork3};
use fs_solver::op::multilevel::SparseMultilevel;
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

/// Fine operator and its inseparable geometry-owned correction space.
/// Only stiffness scales can change; numerical preparations immutably borrow
/// the operator. Coarse geometry is needed only during construction.
pub struct AdaptiveSolveSpace3 {
    operator: AdaptiveElasticity3,
    interpolation: Option<Csr>,
    options: AdaptiveSolveOptions3,
    hierarchy: Option<AdaptiveHierarchy3>,
}
impl AdaptiveSolveSpace3 {
    /// Opt into exact constrained Jacobi for repeated density evaluations.
    #[must_use]
    pub fn jacobi(operator: AdaptiveElasticity3, max_diagonal_contributions: usize) -> Self {
        Self { operator, interpolation: None, hierarchy: None,
            options: AdaptiveSolveOptions3 { max_diagonal_contributions, ..Default::default() } }
    }
    /// Bind one geometric correction space to the owned fine operator. No
    /// stiffness applications or coarse factorization occur until preparation.
    pub fn two_level(operator: AdaptiveElasticity3, coarse: &AdaptiveElasticity3,
        options: AdaptiveSolveOptions3, mut checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<Self, AdaptivePreconditionError3> {
        let transfer = AdaptiveTransfer3::new(coarse, &operator, options.max_transfer_terms, &mut checkpoint)
            .map_err(AdaptivePreconditionError3::Physics)?;
        let interpolation = transfer.vector_prolongation(options.two_level, &mut checkpoint)?;
        drop(transfer);
        poll(&mut checkpoint).map_err(AdaptivePreconditionError3::Physics)?;
        Ok(Self { operator, interpolation: Some(interpolation), options, hierarchy: None })
    }
    /// Bind a sequence of strictly coarser geometric spaces, nearest first.
    /// Intermediate spaces may exceed 512 coordinates; only the final compact
    /// space must fit the bottom-factor cap. Every adjacent pair uses the same
    /// admitted Q1/hanging-node transfer as enrichment. Geometry construction
    /// does not perform a stiffness solve or Galerkin product.
    pub fn multilevel(operator: AdaptiveElasticity3, coarser: &[&AdaptiveElasticity3],
        options: AdaptiveMultilevelOptions3, checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<Self, AdaptivePreconditionError3> {
        let hierarchy = AdaptiveHierarchy3::new(&operator, coarser, options, checkpoint)?;
        Ok(Self { operator, interpolation: None, options: AdaptiveSolveOptions3::default(), hierarchy: Some(hierarchy) })
    }
    /// Actual retained physics; no mutable replacement of its geometry.
    #[must_use] pub fn elasticity(&self) -> &AdaptiveElasticity3 { &self.operator }
    /// Explicitly discard preparation geometry and recover the physical model.
    #[must_use] pub fn into_elasticity(self) -> AdaptiveElasticity3 { self.operator }
    /// Transactional density change through the original material admission.
    pub fn set_scales(&mut self, scales: &[f64]) -> Result<(), ElasticityError3> { self.operator.set_scales(scales) }
    /// First correction dimension; zero for Jacobi. Not necessarily bottom size.
    #[must_use] pub fn coarse_dofs(&self) -> usize {
        match &self.hierarchy {
            Some(h) => h.transfers[0].ncols(),
            None => self.interpolation.as_ref().map_or(0, Csr::ncols),
        }
    }
    /// Total retained interpolation entries; zero for Jacobi.
    #[must_use] pub fn transfer_entries(&self) -> usize {
        match &self.hierarchy {
            Some(h) => h.transfers.iter().map(Csr::nnz).sum(),
            None => self.interpolation.as_ref().map_or(0, Csr::nnz),
        }
    }
    /// Geometry-only sizes including the finest and compact correction spaces.
    #[must_use] pub fn level_sizes(&self) -> Vec<usize> {
        let mut sizes = vec![self.n()];
        if let Some(h) = &self.hierarchy { sizes.extend(h.transfers.iter().map(Csr::ncols)); }
        else if let Some(p) = &self.interpolation { sizes.push(p.ncols()); }
        sizes
    }
    /// Original fine-application progress view. Recursive sparse setup uses no
    /// fine applications; use `prepare_with_work` to observe its local products.
    /// Both entry points poll the same work boundaries and return the same action.
    pub fn prepare(&self, mut checkpoint: impl FnMut(TwoLevelWork) -> ControlFlow<()>)
        -> Result<AdaptivePrepared3<'_>, AdaptivePreconditionError3> {
        self.prepare_with_work(|w| checkpoint(TwoLevelWork { operator_applications: w.operator_applications }))
    }
    /// Prepare CURRENT density once, then reuse for its independent RHS family.
    /// Recursive setup contracts retained bulk and ghost terms into the first
    /// sparse coarse operator, without full fine assembly or coordinate probing.
    /// Failed/cancelled construction exposes spent products but no partial action.
    pub fn prepare_with_work(&self, mut checkpoint: impl FnMut(AdaptiveSetupWork3) -> ControlFlow<()>)
        -> Result<AdaptivePrepared3<'_>, AdaptivePreconditionError3> {
        let cap = self.hierarchy.as_ref().map_or(self.options.max_diagonal_contributions,
            |h| h.options.max_diagonal_contributions);
        let diagonal = self.operator.prepare_jacobi(cap,
            || checkpoint(AdaptiveSetupWork3::default())).map_err(AdaptivePreconditionError3::Physics)?;
        if let Some(hierarchy) = &self.hierarchy {
            return hierarchy.prepare(&self.operator, diagonal.inverse_diagonal(), checkpoint).map(AdaptivePrepared3::Multilevel);
        }
        match &self.interpolation {
            None => Ok(AdaptivePrepared3::Jacobi(diagonal)),
            Some(p) => AdditiveTwoLevel::new(&self.operator, diagonal.inverse_diagonal(), p.clone(),
                self.options.two_level, |w| checkpoint(AdaptiveSetupWork3 {
                    operator_applications: w.operator_applications, galerkin_products: 0,
                })).map(AdaptivePrepared3::TwoLevel).map_err(AdaptivePreconditionError3::Coarse),
        }
    }
}
impl LinearOp for AdaptiveSolveSpace3 {
    fn n(&self) -> usize { self.operator.n() }
    fn apply(&self, x: &[f64], y: &mut [f64]) { self.operator.apply(x, y); }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) { self.operator.apply_transpose(x, y); }
}

/// One immutable-density numerical action, reusable across an RHS family.
pub enum AdaptivePrepared3<'a> {
    Jacobi(AdaptiveJacobi3<'a>),
    TwoLevel(AdditiveTwoLevel<'a, AdaptiveElasticity3>),
    Multilevel(SparseMultilevel<'a, AdaptiveElasticity3>),
}
impl AdaptivePrepared3<'_> {
    /// Fine-application count; recursive assembly instead reports products in
    /// `setup_work()`. This method retains the original two-level accounting.
    #[must_use] pub fn work(&self) -> TwoLevelWork {
        TwoLevelWork { operator_applications: self.setup_work().operator_applications }
    }
    #[must_use] pub fn setup_work(&self) -> AdaptiveSetupWork3 {
        match self {
            Self::Jacobi(_) => AdaptiveSetupWork3::default(),
            Self::TwoLevel(p) => AdaptiveSetupWork3 { operator_applications: p.work().operator_applications, galerkin_products: 0 },
            Self::Multilevel(p) => AdaptiveSetupWork3 { operator_applications: 0, galerkin_products: p.work().galerkin_products },
        }
    }
    #[must_use] pub fn operator(&self) -> &AdaptiveElasticity3 {
        match self { Self::Jacobi(p) => p.operator(), Self::TwoLevel(p) => p.operator(), Self::Multilevel(p) => p.operator() }
    }
}
impl Precond for AdaptivePrepared3<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        match self { Self::Jacobi(p) => p.apply(r, z), Self::TwoLevel(p) => p.apply(r, z), Self::Multilevel(p) => p.apply(r, z) }
    }
}
