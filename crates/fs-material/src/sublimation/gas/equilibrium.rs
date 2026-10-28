//! Sublimation from a finite SOLID inventory on the original solid/liquid chart.
//!
//! The bulk enthalpy stores sensible and fusion energy once. Only its solid
//! mass supplies outgoing molecules; latent compatibility uses PURE SOLID h,
//! never the mixture h. The original molecular and gas-energy owners supply
//! both gross streams, reflected heat and the one accepted gas endpoint.
//!
//! This reduction assumes local solid/liquid thermal equilibrium for ONE pure
//! substance, a sharp isothermal fusion plateau and a fixed-pressure caloric
//! chart. Exposed SOLID area is explicitly authored, not inferred from volume
//! fraction. Liquid evaporation, solutes, mushy melting intervals, nucleation,
//! changing pressure/volume work and automatic surface topology are absent.
//! Complete disappearance of the solid phase refuses until an event owner can
//! change the interface law. A finite interval remains first-order explicit.

use super::{
    ClosedGasParcel, GasSpeciesSnapshot, GasSurfaceBinding, GasSurfaceCandidate,
    GasSurfaceError, GasSurfaceSolve, PreparedSurface, SolidSurfaceState,
    SpeciesBoundaryTransfer, SublimationLaw, SurfaceSnapshot, ThermochemicalGasModel,
    admit_surface_inventory, compatibility, finish_surface, invalid, phase_at, poll,
    positive, prepare_surface_with_endpoint, snapshot_surface, species_snapshot, validate_controls,
};
use super::accommodation::{ReflectedHeatStep, ThermalAccommodation};
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use std::collections::BTreeMap;

/// Borrow the original immutable chart; no second caloric table is constructed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumSublimation<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    solidus: EnthalpyPhaseKnot,
    liquidus: EnthalpyPhaseKnot,
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
    /// Admit a sharp pure-substance transition. All source knots are inspected
    /// once; queries thereafter use the original O(log knots) state evaluator.
    /// The plateau must obey h=h_s+f_liquid*(h_l-h_s), within binary64 roundoff.
    pub fn new(curve: &'a EquilibriumEnthalpyPhaseCurve) -> Result<Self, GasSurfaceError> {
        let knots = curve.knots();
        if knots.len() > crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
            return Err(invalid("equilibrium sublimation exceeds its caloric knot budget"));
        }
        let first_mixed = knots.iter().position(|k| k.liquid_mass_fraction > 0.0)
            .ok_or(invalid("solid/liquid sublimation requires an explicit fusion plateau"))?;
        if first_mixed == 0 {
            return Err(invalid("an original pure-solid endpoint is required"));
        }
        let last = knots.iter().position(|k| k.liquid_mass_fraction == 1.0)
            .ok_or(invalid("complete original fusion plateau required"))?;
        let solidus = knots[first_mixed - 1];
        let liquidus = knots[last];
        let latent = liquidus.specific_enthalpy_j_kg - solidus.specific_enthalpy_j_kg;
        if !positive(latent) || solidus.temperature_k != liquidus.temperature_k {
            return Err(invalid("sublimation requires a sharp isothermal fusion plateau"));
        }
        for (i, k) in knots.iter().enumerate() {
            if i < first_mixed {
                if k.liquid_mass_fraction != 0.0
                    || (i > 0 && k.temperature_k <= knots[i - 1].temperature_k) {
                    return Err(invalid("original solid branch must be strictly sensible"));
                }
            } else if i <= last {
                let fraction = (k.specific_enthalpy_j_kg - solidus.specific_enthalpy_j_kg) / latent;
                if k.temperature_k != solidus.temperature_k
                    || (fraction - k.liquid_mass_fraction).abs() > 128.0 * f64::EPSILON {
                    return Err(invalid("original fusion plateau violates the pure-substance lever rule"));
                }
            } else if k.liquid_mass_fraction != 1.0 || k.temperature_k <= knots[i - 1].temperature_k {
                return Err(invalid("original liquid branch must be strictly sensible"));
            }
        }
        Ok(Self { curve, solidus, liquidus })
    }

    #[must_use]
    pub const fn curve(self) -> &'a EquilibriumEnthalpyPhaseCurve { self.curve }

    #[must_use]
    pub fn fusion_latent_heat_j_kg(self) -> f64 {
        self.liquidus.specific_enthalpy_j_kg - self.solidus.specific_enthalpy_j_kg
    }

    /// Source-evaluated phase partition. No liquid is silently made available
    /// for sublimation, and an entirely liquid inventory is not a solid surface.
    pub fn state(self, bulk: SolidSurfaceState) -> Result<PhaseInventoryState, GasSurfaceError> {
        if !positive(bulk.mass_kg) || !bulk.area_m2.is_finite() || bulk.area_m2 < 0.0 {
            return Err(invalid("finite positive condensed mass and explicit solid exposure required"));
        }
        let equilibrium = self.curve.state_at_specific_enthalpy(bulk.specific_enthalpy_j_kg)
            .map_err(|_| invalid("original equilibrium chart refused bulk enthalpy"))?;
        let solid_mass_kg = bulk.mass_kg * equilibrium.solid_mass_fraction();
        let liquid_mass_kg = bulk.mass_kg * equilibrium.liquid_mass_fraction();
        if !positive(solid_mass_kg) || !liquid_mass_kg.is_finite()
            || (equilibrium.liquid_mass_fraction() > 0.0 && liquid_mass_kg == 0.0) {
            return Err(invalid("solid phase is absent or its inventory is unrepresentable"));
        }
        Ok(PhaseInventoryState { equilibrium, mass_kg: bulk.mass_kg, solid_mass_kg, liquid_mass_kg })
    }

    fn pure_solid_h(self, state: EquilibriumPhaseState) -> Result<f64, GasSurfaceError> {
        if state.is_fully_solid() {
            Ok(state.specific_enthalpy_j_kg())
        } else if state.solid_mass_fraction() > 0.0 && state.temperature_k() == self.solidus.temperature_k {
            Ok(self.solidus.specific_enthalpy_j_kg)
        } else {
            Err(invalid("no supported pure-solid caloric endpoint at the mixture temperature"))
        }
    }

    pub(super) fn snapshot<L: SublimationLaw + ?Sized>(
        self, law: &L, binding: &GasSurfaceBinding, bulk: SolidSurfaceState,
        gas: GasSpeciesSnapshot, catalog: &ThermochemicalGasModel,
        duration_s: f64, controls: GasSurfaceSolve,
    ) -> Result<SurfaceSnapshot, GasSurfaceError> {
        let initial = self.state(bulk)?;
        let pure = SolidSurfaceState {
            mass_kg: initial.solid_mass_kg,
            specific_enthalpy_j_kg: self.pure_solid_h(initial.equilibrium)?,
            area_m2: bulk.area_m2,
        };
        // Reuse the original kinetics and pure-solid/vapor compatibility gate.
        // Only storage is replaced by the complete original condensed inventory.
        let mut snapshot = snapshot_surface(law, binding, pure, gas, catalog,
            &|_| Ok(initial.equilibrium.temperature_k()), duration_s, controls)?;
        snapshot.solid = bulk;
        Ok(snapshot)
    }

    pub(super) fn prepare<L: SublimationLaw + ?Sized>(
        self, law: &L, binding: &GasSurfaceBinding, snapshot: SurfaceSnapshot,
        catalog: &ThermochemicalGasModel, controls: GasSurfaceSolve, reflected_heat_j: f64,
    ) -> Result<PreparedSurface, GasSurfaceError> {
        prepare_surface_with_endpoint(binding, snapshot, controls, reflected_heat_j, |h| {
            let state = self.state(SolidSurfaceState { specific_enthalpy_j_kg: h, ..snapshot.solid })?;
            let t = state.equilibrium.temperature_k();
            let properties = phase_at(law, t, snapshot.gas.molar_mass)?;
            let vapor_h = catalog.species_specific_enthalpy_j_kg(&binding.species, t)?;
            compatibility(vapor_h, self.pure_solid_h(state.equilibrium)?,
                binding.solid_enthalpy_offset_j_kg, properties.latent_heat_j_kg, t,
                controls.latent_compatibility_tolerance_j_kg)?;
            Ok((t, properties))
        })
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

/// The original gas endpoint with the actual remaining solid/liquid partition.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumSurfaceCandidate<'curve, 'catalog> {
    inventory: EquilibriumSublimation<'curve>,
    exchange: GasSurfaceCandidate<'catalog>,
    phase: PhaseInventoryStep,
    reflection: Option<ReflectedHeatStep>,
}
impl<'curve, 'catalog> EquilibriumSurfaceCandidate<'curve, 'catalog> {
    /// In this opt-in path, solid() and solid_energy_change_j on the underlying
    /// exchange describe the WHOLE condensed body, not just its solid phase.
    /// Its recession is a vapor-transfer-equivalent thickness, not mesh motion.
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
/// Bulk mass and h belong to the original equilibrium chart. Gross evaporation
/// cannot spend initially liquid mass, even if cooling freezes some liquid in
/// this interval. The accepted h may leave the latent plateau into the solid
/// branch; complete melting refuses rather than switching to evaporation.
/// Optional reflection acts on the SAME declared solid area. No extra fusion
/// heat is supplied: it already resides in the original caloric coordinate.
#[allow(clippy::too_many_arguments)]
pub fn advance_equilibrium_surface<'curve, 'catalog, L: SublimationLaw + ?Sized>(
    inventory: EquilibriumSublimation<'curve>, law: &L, binding: &GasSurfaceBinding,
    bulk: SolidSurfaceState, gas: &ClosedGasParcel, catalog: &'catalog ThermochemicalGasModel,
    accommodation: Option<&ThermalAccommodation>, duration_s: f64, controls: GasSurfaceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumSurfaceCandidate<'curve, 'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_controls(controls, duration_s)?;
    let gas_snapshot = species_snapshot(gas, catalog, &binding.species)?;
    let snapshot = inventory.snapshot(law, binding, bulk, gas_snapshot, catalog, duration_s, controls)?;
    admit_surface_inventory(snapshot, duration_s, controls)?;
    let reflection = match accommodation {
        Some(law) => Some(law.evaluate(binding, gas, snapshot.temperature_k,
            bulk.area_m2, duration_s, &mut cancelled)?),
        None => None,
    };
    let prepared = inventory.prepare(law, binding, snapshot, catalog, controls,
        reflection.as_ref().map_or(0.0, |row| row.heat_into_gas_j))?;
    let mut rows = BTreeMap::new();
    if let Some(reflection) = &reflection {
        for row in &reflection.species {
            rows.insert(row.species.clone(), SpeciesBoundaryTransfer {
                species: row.species.clone(), mass_into_gas_kg: 0.0,
                energy_into_gas_j: row.heat_into_gas_j,
            });
        }
    }
    let difference = snapshot.incoming_specific_energy_j_kg
        .map_or(0.0, |e| snapshot.outgoing_specific_energy_j_kg - e);
    let phase_energy = snapshot.transferred_mass_kg * snapshot.outgoing_specific_energy_j_kg
        + snapshot.incoming_mass_kg * difference;
    let row = rows.entry(binding.species.clone()).or_insert_with(|| SpeciesBoundaryTransfer {
        species: binding.species.clone(), mass_into_gas_kg: 0.0, energy_into_gas_j: 0.0,
    });
    row.mass_into_gas_kg = snapshot.transferred_mass_kg;
    row.energy_into_gas_j += phase_energy;
    let rows: Vec<_> = rows.into_values().collect();
    poll(&mut cancelled)?;
    let endpoint = gas.exchange_species_boundary(catalog, &rows, controls.gas, &mut cancelled)?;
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
