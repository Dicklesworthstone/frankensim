//! Scalar access to the shared solid/liquid sublimation inventory owner.
//!
//! The original equilibrium chart owns bulk sensible/fusion storage. The same
//! SharpMeltingCaloric used by the joint chamber owns constituent admission;
//! the original molecular law and extensive storage update own all exchange.
//! Only initially solid mass is available for emission. Pure-solid h, not
//! mixture h, determines compatibility with the sublimation latent heat.
//!
//! Exposed SOLID area is explicitly authored. Liquid evaporation, nucleation,
//! alloy mushy ranges, pressure-dependent melting and moving topology are not
//! inferred. This scalar facade refuses an all-liquid endpoint; the shared
//! chamber can report exact completion for a surrounding interface event owner.

use super::{
    ClosedGasParcel, GasSurfaceBinding, GasSurfaceCandidate, GasSurfaceError,
    GasSurfaceSolve, SolidSurfaceState, SpeciesBoundaryRelativeTransfer, SublimationError,
    SublimationLaw, ThermochemicalGasModel, admit_surface_inventory, finish_surface,
    invalid, poll, positive, species_snapshot, validate_controls,
};
use super::accommodation::{ReflectedHeatStep, ThermalAccommodation};
use super::chamber::{GasChamberSurface, mixed::{self, CondensedInventory, SharpMeltingCaloric}};
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use std::collections::BTreeMap;

#[cfg(test)]
use crate::phase::EnthalpyPhaseKnot;

/// Borrow the SAME admitted sharp-phase source as the joint chamber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumSublimation<'a> {
    source: SharpMeltingCaloric<'a>,
}

/// Actual phase masses at the same accepted bulk enthalpy and temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseInventoryState {
    pub equilibrium: EquilibriumPhaseState,
    pub mass_kg: f64,
    pub solid_mass_kg: f64,
    pub liquid_mass_kg: f64,
}

/// Conversion between solid and liquid is separate from mass sent to the gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseInventoryStep {
    pub initial: PhaseInventoryState,
    pub final_state: PhaseInventoryState,
    /// Positive liquid -> solid; negative melting. Not gas deposition.
    pub freezing_mass_kg: f64,
    /// L_f times the conversion above, already contained in bulk enthalpy.
    /// This diagnostic MUST NOT be added as another heat source.
    pub fusion_energy_j: f64,
    /// delta(m_s) + mass_into_gas - freezing_mass, in kg.
    pub solid_partition_residual_kg: f64,
}

impl<'a> EquilibriumSublimation<'a> {
    /// Admit the original sharp pure-substance fusion chart once, without
    /// constructing another interpolation table or duplicating source rules.
    pub fn new(curve: &'a EquilibriumEnthalpyPhaseCurve) -> Result<Self, GasSurfaceError> {
        Ok(Self { source: SharpMeltingCaloric::new(curve)? })
    }
    #[must_use]
    pub const fn curve(self) -> &'a EquilibriumEnthalpyPhaseCurve { self.source.curve() }
    /// Same immutable source binding for the joint/adaptive partitioned chamber.
    #[must_use]
    pub const fn source(self) -> SharpMeltingCaloric<'a> { self.source }
    #[must_use]
    pub fn fusion_latent_heat_j_kg(self) -> f64 { self.source.fusion_latent_heat_j_kg() }

    /// Source-evaluated phase partition. This scalar lane requires a remaining
    /// solid seed; liquid mass is never silently available for sublimation.
    pub fn state(self, bulk: SolidSurfaceState) -> Result<PhaseInventoryState, GasSurfaceError> {
        if !bulk.area_m2.is_finite() || bulk.area_m2 < 0.0 {
            return Err(invalid("finite nonnegative exposed solid area required"));
        }
        let state = CondensedInventory::new(self.source, bulk)?;
        if !positive(state.solid_mass_kg) {
            return Err(invalid("solid phase is absent or its inventory is unrepresentable"));
        }
        Ok(PhaseInventoryState { equilibrium: state.caloric, mass_kg: bulk.mass_kg,
            solid_mass_kg: state.solid_mass_kg, liquid_mass_kg: state.liquid_mass_kg })
    }

    /// Phase balances reconstructed from actual endpoints, not from gross rates.
    pub fn partition_step(
        self, initial: SolidSurfaceState, final_bulk: SolidSurfaceState,
        transferred_mass_kg: f64, mass_tolerance_kg: f64,
    ) -> Result<PhaseInventoryStep, GasSurfaceError> {
        if !positive(mass_tolerance_kg) || !transferred_mass_kg.is_finite() {
            return Err(invalid("finite transfer and positive partition mass tolerance required"));
        }
        let initial = self.state(initial)?;
        let final_state = self.state(final_bulk)?;
        let freezing_mass_kg = initial.liquid_mass_kg - final_state.liquid_mass_kg;
        let residual = (final_state.solid_mass_kg - initial.solid_mass_kg)
            + transferred_mass_kg - freezing_mass_kg;
        let total_residual = (final_state.mass_kg - initial.mass_kg) + transferred_mass_kg;
        let fusion_energy_j = freezing_mass_kg * self.fusion_latent_heat_j_kg();
        if !residual.is_finite() || !total_residual.is_finite() || !fusion_energy_j.is_finite()
            || residual.abs().max(total_residual.abs()) > mass_tolerance_kg {
            return Err(invalid("actual solid/liquid/vapor inventories failed their mass balance"));
        }
        Ok(PhaseInventoryStep { initial, final_state, freezing_mass_kg, fusion_energy_j,
            solid_partition_residual_kg: residual })
    }
}

/// Original gas endpoint and actual remaining solid/liquid inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumSurfaceCandidate<'curve, 'catalog> {
    inventory: EquilibriumSublimation<'curve>,
    exchange: GasSurfaceCandidate<'catalog>,
    phase: PhaseInventoryStep,
    reflection: Option<ReflectedHeatStep>,
}
impl<'curve, 'catalog> EquilibriumSurfaceCandidate<'curve, 'catalog> {
    /// Underlying solid() and solid_energy_change_j describe the WHOLE
    /// condensed body here. Recession is a vapor-transfer-equivalent thickness,
    /// not actual mixed-body geometry or interface motion from melting.
    #[must_use]
    pub const fn exchange(&self) -> &GasSurfaceCandidate<'catalog> { &self.exchange }
    #[must_use]
    pub const fn phase(&self) -> &PhaseInventoryStep { &self.phase }
    #[must_use]
    pub const fn inventory(&self) -> EquilibriumSublimation<'curve> { self.inventory }
    #[must_use]
    pub const fn reflection(&self) -> Option<&ReflectedHeatStep> { self.reflection.as_ref() }
    #[must_use]
    pub fn into_parts(self) -> (GasSurfaceCandidate<'catalog>, PhaseInventoryStep, Option<ReflectedHeatStep>) {
        (self.exchange, self.phase, self.reflection)
    }
}

/// Sublimate/deposit at an explicitly exposed solid part of a mixed inventory.
///
/// Gross emission cannot spend initially liquid mass, including material frozen
/// later in this interval. The accepted h may leave the fusion plateau into the
/// original solid branch. Reflection acts on the same declared solid area.
/// No extra fusion heat is supplied: it already resides in bulk enthalpy.
/// All intermediate preparation is private; either both owners accept or neither.
#[allow(clippy::too_many_arguments)]
pub fn advance_equilibrium_surface<'curve, 'catalog>(
    inventory: EquilibriumSublimation<'curve>, law: &dyn SublimationLaw, binding: &GasSurfaceBinding,
    bulk: SolidSurfaceState, gas: &ClosedGasParcel, catalog: &'catalog ThermochemicalGasModel,
    accommodation: Option<&ThermalAccommodation>, duration_s: f64, controls: GasSurfaceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumSurfaceCandidate<'curve, 'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_controls(controls, duration_s)?;
    let caloric = |h| inventory.curve().temperature_at_specific_enthalpy(h)
        .map_err(|_| SublimationError::Invalid("original solid/liquid caloric chart refused state"));
    let surface = GasChamberSurface { id: 0, law, binding, solid: bulk,
        temperature_from_enthalpy: &caloric };
    let gas_snapshot = species_snapshot(gas, catalog, &binding.species)?;
    let snapshot = mixed::snapshot(&surface, inventory.source, gas_snapshot, catalog, duration_s, controls)?;
    admit_surface_inventory(snapshot, duration_s, controls)?;
    let reflection = match accommodation {
        Some(law) => Some(law.evaluate(binding, gas, snapshot.temperature_k,
            bulk.area_m2, duration_s, &mut cancelled)?),
        None => None,
    };
    let prepared = mixed::prepare(&surface, inventory.source, snapshot, catalog, controls,
        reflection.as_ref().map_or(0.0, |row| row.heat_into_gas_j))?;
    let mut rows = BTreeMap::new();
    if let Some(reflection) = &reflection {
        for row in &reflection.species {
            rows.insert(row.species.clone(), SpeciesBoundaryRelativeTransfer {
                species: row.species.clone(), mass_into_gas_kg: 0.0,
                relative_energy_into_gas_j: row.heat_into_gas_j,
            });
        }
    }
    let row = rows.entry(binding.species.clone()).or_insert_with(|| SpeciesBoundaryRelativeTransfer {
        species: binding.species.clone(), mass_into_gas_kg: 0.0, relative_energy_into_gas_j: 0.0,
    });
    row.mass_into_gas_kg = snapshot.transferred_mass_kg;
    row.relative_energy_into_gas_j += prepared.phase_relative_energy_into_gas_j;
    let rows: Vec<_> = rows.into_values().collect();
    poll(&mut cancelled)?;
    let endpoint = gas.exchange_species_boundary_relative(catalog, &rows, controls.gas, &mut cancelled)?;
    if let Some(accommodation) = accommodation {
        accommodation.admit_temperatures(prepared.surface_temperature_k[1], endpoint.gas().state().state().temperature)?;
    }
    let phase = inventory.partition_step(bulk, prepared.solid, prepared.transferred_mass_kg,
        controls.gas.mass_tolerance_kg)?;
    let exchange = finish_surface(binding, bulk, gas, gas_snapshot, prepared, endpoint, controls, cancelled)?;
    Ok(EquilibriumSurfaceCandidate { inventory, exchange, phase, reflection })
}

#[cfg(test)]
mod tests;
