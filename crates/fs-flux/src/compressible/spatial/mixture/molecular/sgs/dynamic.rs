//! Favre Germano viscosity on the original conservative mixture owner.
//!
//! The shared periodic 27-cell test filter fits the deviatoric coefficient
//! from the same-stage mixture density and momentum/density. The physical
//! filter widths are `h` and `3h`; neither a species average nor a constant
//! reference density replaces the source-admitted local density. The fit
//! supplies only `nu_t`. The parent retains the original Maxwell--Stefan
//! fluxes, adds `rho*nu_t` stress/work, source-Cp heat and common-Schmidt
//! species/partial-enthalpy transport, and checks their summed stability rate.
//!
//! This explicit exploration rung requires a fully periodic isotropic grid
//! with at least three cells per axis. It supplies no wall filter, dynamic
//! Prandtl/Schmidt model, reaction-rate closure, pressure dilatation, or SGS
//! energy equation. Negative fitted coefficients use the shared no-backscatter
//! policy; excessive positive coefficients or viscosity ratios refuse.
//! Every original SSP stage and retry recomputes the fit. Scratch is reused
//! and no fit report is published from an incomplete transport interval.

use super::*;
use fs_lbm::sgs3::dynamic::{DynamicSgsError, PeriodicDynamicSgs};

/// Dynamic momentum closure plus explicitly authored heat/species diffusivity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureDynamicSubgridConfig {
    pub transport: DynamicSubgridConfig,
    /// Positive finite common turbulent Schmidt number. Molecular species
    /// diffusivities remain with the original multicomponent source.
    pub turbulent_schmidt: f64,
}

impl MixtureDynamicSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate().map_err(|_| {
            DissipationError::Invalid("invalid dynamic mixture stress/heat controls")
        })?;
        if !positive(self.turbulent_schmidt) {
            return Err(DissipationError::Invalid(
                "positive finite turbulent Schmidt number required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Combined flow/operator payload including the nonlocal filter, staged
    /// density/velocity, and SGS viscosity/Cv arrays: 184 extra bytes per cell.
    /// Actual owned capacities are checked again after allocation.
    pub fn required_workspace_bytes_with_dynamic_sgs(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        let core = PeriodicDynamicSgs::required_bytes(flow.cell_count()).map_err(dynamic_error)?;
        let staging = flow.cell_count().checked_mul(4 * size_of::<f64>()).ok_or(
            DissipationError::Invalid("dynamic mixture workspace size overflow"),
        )?;
        Self::required_workspace_bytes_with_sgs(flow, source)?
            .checked_add(core)
            .and_then(|bytes| bytes.checked_add(staging))
            .ok_or(DissipationError::Invalid(
                "dynamic mixture workspace size overflow",
            ))
    }

    /// Prepare a dynamic SGS producer for the existing conservative transport.
    /// The byte limit covers the original flow and the complete operator.
    /// The original reacting split integrator can consume this same operator.
    pub fn with_dynamic_sgs(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        config: MixtureDynamicSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        admit_domain(flow)?;
        admit_workspace(
            Self::required_workspace_bytes_with_dynamic_sgs(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let mut result = Self::new(flow, source, maximum_workspace_bytes, &mut cancelled)?;
        let dynamic = DynamicTransport::new(&flow.grid, config)?;
        result.sgs = Some(SubgridTransport {
            closure: ViscosityClosure::Dynamic(dynamic),
            dynamic_viscosity: allocated(flow.cell_count(), 0.0)?,
            heat_capacity_v_density: allocated(flow.cell_count(), 0.0)?,
        });
        let bytes = result
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "dynamic mixture workspace size overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact dynamic controls; local-coefficient and molecular owners return None.
    pub fn dynamic_subgrid_config(&self) -> Option<MixtureDynamicSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            ViscosityClosure::Dynamic(dynamic) => Some(dynamic.config),
            ViscosityClosure::Local { .. } => None,
        }
    }
}

fn admit_domain(flow: &MixtureEuler3d) -> Result<(), DissipationError> {
    if (0..3).any(|axis| flow.boundaries.axis(axis) != AxisBoundaries::Periodic)
        || flow.grid.shape.into_iter().any(|n| n < 3)
        || flow.grid.spacing[0] != flow.grid.spacing[1]
        || flow.grid.spacing[0] != flow.grid.spacing[2]
    {
        return Err(DissipationError::Invalid(
            "dynamic mixture SGS requires periodic isotropic axes with at least three cells each",
        ));
    }
    Ok(())
}

pub(super) struct DynamicTransport {
    pub(super) config: MixtureDynamicSubgridConfig,
    filter: PeriodicDynamicSgs,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
}

impl DynamicTransport {
    fn new(
        grid: &CartesianGrid,
        config: MixtureDynamicSubgridConfig,
    ) -> Result<Self, DissipationError> {
        let core_bytes = PeriodicDynamicSgs::required_bytes(grid.len()).map_err(dynamic_error)?;
        Ok(Self {
            config,
            filter: PeriodicDynamicSgs::new(grid.shape, grid.spacing[0], core_bytes)
                .map_err(dynamic_error)?,
            density: allocated(grid.len(), 0.0)?,
            velocity: allocated(grid.len(), [0.0; 3])?,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        self.filter.workspace_bytes()
            + self.density.capacity() * size_of::<f64>()
            + self.velocity.capacity() * size_of::<[f64; 3]>()
    }

    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DissipationError> {
        for (i, p) in primitive.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            self.density[i] = p.density;
            self.velocity[i] = p.velocity;
        }
        let report = self
            .filter
            .prepare(
                self.config.transport.model,
                &self.density,
                &self.velocity,
                cancelled,
            )
            .map_err(dynamic_error)?;
        Ok(report.maximum_eddy_viscosity)
    }

    pub(super) fn viscosities(&self) -> &[f64] {
        self.filter.viscosities()
    }
}

// Preserve the original molecular error owner without introducing a dependency
// from the lower gas layer back to LBM. Cancellation keeps its typed identity.
fn dynamic_error(error: DynamicSgsError) -> DissipationError {
    match error {
        DynamicSgsError::Cancelled => DissipationError::Cancelled,
        DynamicSgsError::Invalid(reason) | DynamicSgsError::Cell { reason, .. } => {
            DissipationError::Invalid(reason)
        }
        DynamicSgsError::WorkspaceBytes { .. } => {
            DissipationError::Invalid("dynamic mixture filter workspace budget exceeded")
        }
        DynamicSgsError::Allocation => {
            DissipationError::Invalid("dynamic mixture filter allocation failed")
        }
        DynamicSgsError::Coefficient { .. } => {
            DissipationError::Invalid("dynamic mixture positive coefficient bound exceeded")
        }
    }
}
