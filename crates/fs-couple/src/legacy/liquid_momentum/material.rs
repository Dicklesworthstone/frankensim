//! Reciprocal momentum through original frozen, rheological and thermal liquids.

use super::*;
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
    MaterialFreeSurfaceMomentum3,
};

#[derive(Debug, Clone, PartialEq)]
pub enum MaterialLiquidMomentumError {
    Momentum(LiquidMomentumError),
    Fluid(MaterialFreeSurfaceImpulseError),
}
impl core::fmt::Display for MaterialLiquidMomentumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material liquid/rigid drag: {self:?}")
    }
}
impl std::error::Error for MaterialLiquidMomentumError {}
impl From<LiquidMomentumError> for MaterialLiquidMomentumError {
    fn from(error: LiquidMomentumError) -> Self {
        Self::Momentum(error)
    }
}
impl From<MaterialFreeSurfaceImpulseError> for MaterialLiquidMomentumError {
    fn from(error: MaterialFreeSurfaceImpulseError) -> Self {
        Self::Fluid(error)
    }
}

/// Exact paired drag using the liquid's ORIGINAL immutable material scales.
///
/// The grid axes coincide with the inertial world axes, with an explicit world
/// corner. Source density, dt/dx, scenario speed and retained liquid masses
/// cannot be independently substituted. The common source preserves those
/// bindings and deposits no heat; its physical loss is returned for an explicit
/// thermal owner. This is a frozen-pose source duration, not elapsed transport
/// or rigid drift. Both original states remain unchanged on refusal.
///
/// Frozen, source-rheology and enthalpy owners all retain their own original
/// native flow, viscosity/pressure samples, thermal history and chronology.
/// No flow clone, constitutive re-query or implicit heat allocation is used.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_drag<L: MaterialFreeSurfaceMomentum3 + ?Sized>(
    liquid: &mut L,
    body: &mut RigidBodyState,
    properties: MassProperties,
    indices: &[usize],
    grid_origin_world_m: Vec3,
    relaxation_rate_per_s: f64,
    source_duration_s: f64,
    controls: LiquidRigidDragControls,
    gate: &CancelGate,
) -> Result<LiquidRigidDragReport, MaterialLiquidMomentumError> {
    poll(gate)?;
    controls.admit(indices.len())?;
    let scale = liquid.momentum_scaling();
    let units = LiquidMomentumUnits3 {
        cell_size_m: scale.cell_size_m(),
        velocity_unit_m_s: scale
            .velocity_si(1.0)
            .map_err(MaterialFreeSurfaceImpulseError::from)?,
        mass_unit_kg: scale.mass_unit_kg(),
        grid_origin_world_m,
    };
    units.admit()?;
    let stage = liquid.prepare_native_impulse_stencil(
        indices,
        MaterialFreeSurfaceImpulseControls3 {
            max_cells: controls.maximum_cells,
            maximum_speed_m_s: controls.maximum_speed_m_s,
            maximum_momentum_defect_kg_m_s: controls.maximum_momentum_defect_kg_m_s,
            maximum_work_defect_j: controls.maximum_work_defect_j,
        },
        gate,
    )?;
    Ok(apply_staged(
        stage,
        body,
        properties,
        units,
        relaxation_rate_per_s,
        source_duration_s,
        controls,
        gate,
    )?)
}

/// Directional, spatially varying reciprocal drag on the ORIGINAL material
/// liquid. Native mass, physical scales and scenario speed stay source-bound.
/// The explicit symmetric substep schedule controls composition accuracy;
/// independent physical dissipation is returned without assigning heating.
/// Thermal and non-Newtonian owners participate through the same sealed
/// original-material transaction; accepted history and all clocks survive.
#[allow(clippy::too_many_arguments)]
pub fn apply_heterogeneous_rigid_drag<L: MaterialFreeSurfaceMomentum3 + ?Sized>(
    liquid: &mut L,
    body: &mut RigidBodyState,
    properties: MassProperties,
    indices: &[usize],
    laws: &[heterogeneous::PrincipalLiquidDrag3],
    grid_origin_world_m: Vec3,
    source_duration_s: f64,
    controls: heterogeneous::HeterogeneousRigidDragControls,
    gate: &CancelGate,
) -> Result<heterogeneous::HeterogeneousRigidDragReport, MaterialLiquidMomentumError> {
    poll(gate)?;
    controls.admit(indices.len(), laws, source_duration_s)?;
    let scale = liquid.momentum_scaling();
    let units = LiquidMomentumUnits3 {
        cell_size_m: scale.cell_size_m(),
        velocity_unit_m_s: scale
            .velocity_si(1.0)
            .map_err(MaterialFreeSurfaceImpulseError::from)?,
        mass_unit_kg: scale.mass_unit_kg(),
        grid_origin_world_m,
    };
    units.admit()?;
    let physical = controls.exchange;
    let stage = liquid.prepare_native_impulse_stencil(
        indices,
        MaterialFreeSurfaceImpulseControls3 {
            max_cells: physical.maximum_cells,
            maximum_speed_m_s: physical.maximum_speed_m_s,
            maximum_momentum_defect_kg_m_s: physical.maximum_momentum_defect_kg_m_s,
            maximum_work_defect_j: physical.maximum_work_defect_j,
        },
        gate,
    )?;
    Ok(heterogeneous::apply_staged(
        stage,
        body,
        properties,
        laws,
        units,
        source_duration_s,
        controls,
        gate,
    )?)
}

#[cfg(test)]
mod tests;
