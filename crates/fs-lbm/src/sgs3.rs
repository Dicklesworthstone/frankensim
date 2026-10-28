//! Explicit, exploratory three-dimensional subgrid viscosity closures.
//!
//! Smagorinsky uses `(C_s Delta)^2 sqrt(2 S:S)`. WALE uses the traceless
//! symmetric square of the FULL velocity gradient, following Nicoud & Ducros
//! (1999), doi:10.1023/A:1009995426001, equations (1), (10) and (13). The
//! antisymmetric gradient is essential: symmetric nonequilibrium stress alone
//! cannot determine WALE. Both models vanish quadratically with filter width
//! for a fixed resolved gradient. Neither supplies unresolved fluctuations,
//! transition prediction, an energy equation, or a DNS/acceptance claim.
//!
//! No model or coefficient is a default. [`SgsUse::AcceptanceRender`] is
//! refused. The grid integration is an explicit per-step choice and never
//! changes the constructor's molecular collision law.

mod model;
pub use model::{SgsConfig, SgsError, SgsModel, SgsUse, eddy_viscosity};

/// Grid failures and closure refusals remain distinguishable to callers.
#[derive(Clone, Debug)]
pub enum SgsStepError<E> {
    Grid(crate::lattice3::GridError3<E>),
    Closure(SgsError),
}

impl<E: core::fmt::Debug> core::fmt::Display for SgsStepError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-SGS3-STEP: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SgsStepError<E> {}
impl<E> From<crate::lattice3::GridError3<E>> for SgsStepError<E> {
    fn from(error: crate::lattice3::GridError3<E>) -> Self {
        Self::Grid(error)
    }
}
impl<E> From<SgsError> for SgsStepError<E> {
    fn from(error: SgsError) -> Self {
        Self::Closure(error)
    }
}

/// A lattice whose molecular shear viscosity can be read without changing
/// its equilibrium, independent bulk rate, or higher moment relaxation.
pub trait SgsLattice3<const Q: usize>: crate::lattice3::ViscousLattice3<Q> {
    fn molecular_viscosity(collision: Self::Collision) -> Result<f64, Self::Error>;
}

impl SgsLattice3<19> for crate::lattice3::D3q19 {
    fn molecular_viscosity(collision: Self::Collision) -> Result<f64, Self::Error> {
        collision.validate()?;
        Ok(collision.kinematic_viscosity())
    }
}

impl SgsLattice3<27> for crate::d3q27::D3q27 {
    fn molecular_viscosity(collision: Self::Collision) -> Result<f64, Self::Error> {
        use crate::d3q27::CollisionModel27;
        collision.validate()?;
        match collision.model {
            CollisionModel27::Bgk { tau, .. } => Ok((tau - 0.5) / 3.0),
            CollisionModel27::CentralMoment { rates } | CollisionModel27::Cumulant { rates } => {
                rates.viscosity()
            }
        }
    }
}

/// Instantaneous, pre-step dissipation diagnostics in lattice units. These
/// are shear-only sums of `2 rho nu S_dev:S_dev` over unit-volume fluid
/// cells, not integrated heat or a turbulence-validation certificate. Bulk
/// dissipation is not inferred from the shear rate. Fixed cell order is used.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SgsReport {
    pub fluid_cells: usize,
    pub resolved_dissipation_lu: f64,
    pub sgs_dissipation_lu: f64,
    pub maximum_eddy_viscosity_lu: f64,
    pub maximum_eddy_to_molecular: f64,
}

impl SgsReport {
    /// Zero when both dissipations vanish. Infinite only when modeled
    /// dissipation is nonzero but the molecular contribution is zero.
    #[must_use]
    pub fn dissipation_ratio(self) -> f64 {
        if self.resolved_dissipation_lu == 0.0 && self.sgs_dissipation_lu == 0.0 {
            0.0
        } else {
            self.sgs_dissipation_lu / self.resolved_dissipation_lu
        }
    }
}

/// Reused velocity and collision fields: no per-cell or per-step allocation
/// in the SGS preparation sweep after this workspace has been constructed.
pub struct SgsWorkspace3<C> {
    pub(crate) velocity: Vec<[f64; 3]>,
    pub(crate) density: Vec<f64>,
    pub(crate) collisions: Vec<C>,
}

impl<C> SgsWorkspace3<C> {
    /// Prepared total-viscosity collision field for the existing per-cell
    /// stepping seam. Its contents are valid only after successful preparation.
    #[must_use]
    pub fn collisions(&self) -> &[C] {
        &self.collisions
    }
}
