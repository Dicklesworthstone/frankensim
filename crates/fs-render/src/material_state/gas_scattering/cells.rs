//! Borrowed, cell-centred accepted gas states on the native volume geometry.

use super::{GasGridBindingError, MaterialGasScatteringBinding};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::ScatteringGrid;
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, VolumeTransportError,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::gas_scattering::RayleighGasState;

impl MaterialGasScatteringBinding {
    /// Borrow a heterogeneous accepted gas snapshot at this binding's fixed
    /// composition, gas specs and molecular optical references.
    ///
    /// `grid` is the single geometry owner: origin and cell sizes are metres,
    /// density is kg/m³, and both buffers use the cell-centred, x-major order
    /// `i + nx * (j + ny * k)`. Each state is constant over its complete cell;
    /// there is no interpolation of density, thermodynamics or optical data.
    /// The density buffer must contain each original state's exact density.
    /// Dimensions, cell counts and physical bounds are validated before return.
    ///
    /// Supply states from `MolecularRayleighMixture::evaluate` or
    /// `MolecularRayleighMixture::bind_gas_state`. Those material operations
    /// retain the original gas evaluation and coherence checks. This adapter
    /// neither reevaluates a gas law nor infers composition from its density.
    /// Changing composition, species order, gas specs or optical source refuses,
    /// even when the resulting coefficients happen to be numerically identical.
    /// Pressure and temperature may vary freely within the material's domain.
    ///
    /// Transport uses exactly `rho[cell] * kappa(lambda)` with the original
    /// analytic inverse-fourth spectrum and depolarized phase. This is elastic
    /// molecular scattering only; local temperature does not become thermal
    /// emission. The original states and density remain borrowed without a
    /// copied field or per-cell optical allocation. Cancellation is polled at
    /// most every 1024 state checks and by the native grid admission owners.
    pub fn bind_cell_grid<'a>(
        &self,
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        states: &'a [RayleighGasState],
        budget: AbsorbingGridBudget,
    ) -> Result<MaterialGasScatteringCellGrid<'a>, GasGridBindingError> {
        checkpoint(cx)?;
        if grid.data.len() > budget.max_grid_cells || budget.max_steps_per_segment == 0 {
            return Err(GasGridBindingError::Transport(VolumeTransportError::Budget));
        }
        if states.len() != grid.data.len() {
            return Err(GasGridBindingError::StateCount);
        }
        let density_kg_m3 = grid.data;
        // Let the existing physical-grid owner validate the one geometry and
        // density buffer. No substitute grid semantics or gas law is introduced.
        let absorption = AbsorbingGrid::try_new(
            cx,
            grid,
            MassAbsorptionSpectrum::constant(0.0).map_err(GasGridBindingError::Transport)?,
            None,
            budget,
        )
        .map_err(GasGridBindingError::Transport)?;
        let model_identity = self.state.model_identity();
        let mut hash = DomainHasher::new("org.frankensim.render.material-gas-cell-grid.v1");
        hash.update(model_identity.as_bytes());
        for (index, state) in states.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            if state.model_identity() != model_identity {
                return Err(GasGridBindingError::ModelMismatch { index });
            }
            if density_kg_m3[index] != state.gas().density {
                return Err(GasGridBindingError::DensityMismatch { index });
            }
            hash.update(state.content_identity().as_bytes());
        }
        let medium = ScatteringGrid::try_new(cx, absorption, self.scattering, self.phase.clone())
            .map_err(GasGridBindingError::Transport)?;
        hash.update(medium.content_identity().as_bytes());
        Ok(MaterialGasScatteringCellGrid {
            states,
            density_kg_m3,
            medium,
            identity: hash.finalize(),
        })
    }
}

/// A heterogeneous native scattering grid borrowing the original accepted gas
/// cells and density. Geometry and traversal remain with `ScatteringGrid`.
pub struct MaterialGasScatteringCellGrid<'a> {
    states: &'a [RayleighGasState],
    density_kg_m3: &'a [f64],
    medium: ScatteringGrid<'a>,
    identity: ContentHash,
}

impl<'a> MaterialGasScatteringCellGrid<'a> {
    /// Direct input to the native path tracer and physical-grid photometry.
    #[must_use]
    pub const fn medium(&self) -> &ScatteringGrid<'a> {
        &self.medium
    }

    /// Original accepted cell states, with temperature, pressure and sources.
    /// The order is the native density grid's x-major cell order.
    #[must_use]
    pub const fn states(&self) -> &'a [RayleighGasState] {
        self.states
    }

    /// The original density buffer borrowed by the native scattering medium.
    #[must_use]
    pub const fn density_kg_m3(&self) -> &'a [f64] {
        self.density_kg_m3
    }

    /// Complete accepted state sequence, declared fixed model and native grid.
    /// State changes remain visible even when their density is unchanged.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), GasGridBindingError> {
    cx.checkpoint()
        .map_err(|_| GasGridBindingError::Transport(VolumeTransportError::Cancelled))
}

#[cfg(test)]
mod tests;
