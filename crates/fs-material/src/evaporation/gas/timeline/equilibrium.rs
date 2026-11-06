//! Sharp fusion, the original heat network and finite-gas exchange on one clock.
//!
//! The parent timeline owns event clipping, heat, retry, birth/exposure semantics
//! and whole-window accounting. Only phase admission and the original material
//! facade differ from its pure-liquid lane. The selected cell retains one bulk
//! enthalpy coordinate; fusion is never an extra thermal source.
//!
//! Exposed geometric area is explicitly admitted LIQUID interface area, including
//! when the body contains solid. Bulk phase fraction does not infer coverage.
//! A liquid seed must remain at every accepted heat and vapor endpoint. Complete
//! freezing, nucleation and interface handoff require an explicit event owner.
//! Geometry stays fixed, and the complete split requires timestep refinement.

use super::{
    ClosedGasParcel, DepositedCell, DepositionError, DepositionSimulation, EquilibriumEvaporation,
    EquilibriumGasHeatSolve, EvaporationLaw, Event, GasEvaporationSolve,
    GasEvaporationTimelineCandidate, GasEvaporationTimelineError, GasFilmHeatSolve,
    GasLiquidBinding, LiquidSurfaceState, PhaseInventoryStep, StepConfig, ThermochemicalGasModel,
    TimelineSurface, advance_deposition,
};
use fs_thermal::enthalpy::CellId;

/// One original condensed cell and its explicitly exposed liquid interface.
///
/// The thermal timeline must register this exact caloric identity and datum.
/// Native sources/contact/convection, gas sensible heat and vapor exchange all
/// evolve that same storage. Material properties and gas volume keep their owners.
#[derive(Clone, Copy)]
pub struct EquilibriumTimelineSurface<'a> {
    pub cell_id: CellId,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub inventory: EquilibriumEvaporation<'a>,
    pub mass_transfer_coefficient_m_s: f64,
}

impl<'a> From<EquilibriumTimelineSurface<'a>> for TimelineSurface<'a> {
    fn from(surface: EquilibriumTimelineSurface<'a>) -> Self {
        let knots = surface.inventory.curve().knots();
        Self {
            cell_id: surface.cell_id,
            law: surface.law,
            binding: surface.binding,
            enthalpy_range_j_kg: [
                knots[0].specific_enthalpy_j_kg,
                knots[knots.len() - 1].specific_enthalpy_j_kg,
            ],
            mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
            inventory: Some(surface.inventory),
        }
    }
}

/// Advance the original heat network and fusion-aware vapor exchange together.
///
/// Births and every heat substep retain the original chart, datum and positive
/// liquid inventory. The complete phase report starts at the original active
/// cell or its accepted birth, including heat-driven changes after covering.
/// Failure and cancellation publish neither thermal nor gas endpoint.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surface: EquilibriumTimelineSurface<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasEvaporationSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasEvaporationTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    advance_deposition(
        thermal,
        surface.into(),
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        None,
        checkpoint,
    )
}

/// Native network heat, reciprocal gas heat, then vapor exchange per interval.
///
/// All three stages share the same original bulk storage. The explicit bulk
/// interval also bounds native heat substeps. The sensible coefficient excludes
/// vapor-carried energy. Zero coefficient delegates exactly to the phase-only
/// timeline and does not inspect unused heat controls.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas_and_heat<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surface: EquilibriumTimelineSurface<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasEvaporationSolve,
    heat_transfer_coefficient_w_m2_k: f64,
    heat_config: EquilibriumGasHeatSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasEvaporationTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    if heat_transfer_coefficient_w_m2_k == 0.0 {
        return advance_deposition_with_gas(
            thermal,
            surface,
            gas,
            catalog,
            target_s,
            thermal_config,
            phase_config,
            checkpoint,
        );
    }
    let mut surface: TimelineSurface<'_> = surface.into();
    surface.enthalpy_range_j_kg = heat_config.bulk_enthalpy_range_j_kg;
    advance_deposition(
        thermal,
        surface,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        Some((
            heat_transfer_coefficient_w_m2_k,
            GasFilmHeatSolve {
                pure_liquid_enthalpy_range_j_kg: heat_config.bulk_enthalpy_range_j_kg,
                temperature_tolerance_k: heat_config.temperature_tolerance_k,
                maximum_iterations: heat_config.maximum_iterations,
            },
        )),
        checkpoint,
    )
}

// Summing vapor-only reports would omit network-driven melting/freezing and
// any heat after the interface becomes covered. Reconstruct the complete window.
pub(super) fn partition_report(
    original: &DepositionSimulation<'_>,
    endpoint: &DepositionSimulation<'_>,
    surface: TimelineSurface<'_>,
    transferred_mass_kg: f64,
    mass_tolerance_kg: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Option<PhaseInventoryStep>, GasEvaporationTimelineError> {
    let Some(inventory) = surface.inventory else {
        return Ok(None);
    };
    let Some(final_cell) = endpoint.cells().get(&surface.cell_id) else {
        return Ok(None);
    };
    let mut initial_cell = original.cells().get(&surface.cell_id);
    if initial_cell.is_none() {
        for (index, event) in original.pending_events().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if event.time_s > endpoint.time_s() {
                break;
            }
            if let Event::Birth(cell) = &event.event {
                if cell.id == surface.cell_id {
                    initial_cell = Some(cell);
                    break;
                }
            }
        }
    }
    let initial_cell = initial_cell.ok_or(DepositionError::InvalidInput(
        "equilibrium evaporation endpoint has no original cell or accepted birth",
    ))?;
    let bulk = |cell: &DepositedCell| LiquidSurfaceState {
        mass_kg: cell.mass_kg,
        specific_enthalpy_j_kg: cell.specific_h_j_kg,
        area_m2: 0.0,
    };
    checkpoint()?;
    Ok(Some(inventory.partition_step(
        bulk(initial_cell),
        bulk(final_cell),
        transferred_mass_kg,
        mass_tolerance_kg,
    )?))
}

#[cfg(test)]
mod tests;
