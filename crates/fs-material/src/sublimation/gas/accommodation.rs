//! Reflected-molecule heat, without replacing phase or gas caloric owners.
//!
//! The local incident gas is Maxwellian at its original temperature T_g. A
//! reflected species carries mean specific energy e(T)=h(T)-R_i*T/2. With an
//! explicitly supplied energy accommodation alpha_i and sticking probability
//! beta_i, heat into the gas is
//!
//! q_i = alpha_i (1-beta_i) p_i / sqrt(2*pi*R_i*T_g)
//!       * [h_i(T_s)-h_i(T_g)-R_i*(T_s-T_g)/2].
//!
//! beta_i is the original volatile sticking coefficient, or zero for a carrier.
//! Reflection returns the same incident mass: it is NOT another evaporation
//! flux at saturation pressure and does not withdraw carrier inventory. The
//! original sticking streams already account for their energy and latent heat.
//! Formation references cancel in the existing NASA-9 enthalpy-difference owner.
//!
//! This is an explicit frozen-rate, stationary-wall, local molecular model with
//! equilibrated internal modes and energy/angle-independent coefficients. It
//! does not solve a Knudsen layer, continuum film, enclosure distribution,
//! separate internal-mode accommodation, surface chemistry or recoil work.
//! Do not add a second gas-wall heat law for the same molecules. Coefficients
//! and their temperature interval are caller inputs, not inferred measurements.
//!
//! The energy-accommodation convention and monatomic limiting heat flux follow
//! Yamaguchi et al. (2017), equations (1)-(2),
//! <https://strathprints.strath.ac.uk/60000/>. Finite-step conservation below is
//! not a time-error or entropy certificate; refine the physical timestep.

use std::collections::BTreeMap;
use std::sync::Arc;
use fs_qty::Temperature;
use super::{
    ClosedGasParcel, GasSurfaceBinding, GasSurfaceCandidate, GasSurfaceError, GasSurfaceSolve,
    SolidSurfaceState, SpeciesBoundaryTransfer, SpeciesId, SublimationError, SublimationLaw,
    ThermochemicalGasModel, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, admit_surface_inventory,
    finish_surface, invalid, nonnegative, poll, positive, prepare_surface_with_heat,
    snapshot_surface, species_snapshot, validate_controls,
};

/// Explicit per-species thermal accommodation on one fixed surface condition.
/// This is distinct from the volatile's mass/sticking accommodation.
/// Clones share the admitted immutable table, avoiding species-string and tree
/// allocations during active-cell setup and adaptive retries. Construction
/// still copies caller inputs once; no mutable coefficient view is exposed.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalAccommodation {
    coefficients: Arc<BTreeMap<SpeciesId, f64>>,
    temperature_range_k: [f64; 2],
}
impl ThermalAccommodation {
    /// Supply exactly one coefficient, including explicit zeros, for every
    /// retained catalog species. The catalog's existing component cap bounds
    /// allocation. Reentry of depleted species therefore needs no default.
    pub fn new(
        catalog: &ThermochemicalGasModel,
        coefficients: &[(SpeciesId, f64)],
        temperature_range_k: [f64; 2],
    ) -> Result<Self, GasSurfaceError> {
        let [lo, hi] = temperature_range_k;
        if !positive(lo) || !positive(hi) || lo > hi
            || coefficients.len() != catalog.thermodynamics().components().len()
        {
            return Err(invalid("complete accommodation roster and positive temperature interval required"));
        }
        let mut ordered = BTreeMap::new();
        for (species, coefficient) in coefficients {
            if !coefficient.is_finite() || !(0.0..=1.0).contains(coefficient)
                || !catalog.thermodynamics().components().iter().any(|row| row.species() == species)
                || ordered.insert(species.clone(), *coefficient).is_some()
            {
                return Err(invalid("unique retained species and thermal accommodation in [0,1] required"));
            }
        }
        Ok(Self { coefficients: Arc::new(ordered), temperature_range_k })
    }

    /// Original explicitly supplied species coefficients; no mixture averaging.
    #[must_use]
    pub fn coefficients(&self) -> &BTreeMap<SpeciesId, f64> { self.coefficients.as_ref() }

    /// Both incident-gas and wall temperatures must remain in this interval.
    #[must_use]
    pub const fn temperature_range_k(&self) -> [f64; 2] { self.temperature_range_k }

    pub(super) fn admit_temperatures(&self, surface: f64, gas: f64) -> Result<(), GasSurfaceError> {
        let [lo, hi] = self.temperature_range_k;
        if [surface, gas].iter().any(|t| !t.is_finite() || !(lo..=hi).contains(t)) {
            return Err(invalid("gas or surface left thermal-accommodation temperature support"));
        }
        Ok(())
    }

    /// Integrated reflected heat at the SAME original gas/surface snapshot as
    /// the sticking streams. No state changes here. A positive result heats gas.
    /// Zero coefficients require no wall-temperature caloric query for that gas.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        binding: &GasSurfaceBinding,
        gas: &ClosedGasParcel,
        surface_temperature_k: f64,
        area_m2: f64,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReflectedHeatStep, GasSurfaceError> {
        poll(&mut cancelled)?;
        let gas_temperature_k = gas.state().state().temperature;
        self.admit_temperatures(surface_temperature_k, gas_temperature_k)?;
        if !nonnegative(area_m2) || !nonnegative(duration_s)
            || !binding.accommodation_coefficient.is_finite()
            || !(0.0..=1.0).contains(&binding.accommodation_coefficient)
        {
            return Err(invalid("finite nonnegative molecular exposure and valid sticking required"));
        }
        let exposure = area_m2 * duration_s;
        if !nonnegative(exposure) || (area_m2 > 0.0 && duration_s > 0.0 && exposure == 0.0) {
            return Err(invalid("reflected molecular exposure is unrepresentable"));
        }
        let model = gas.model().thermodynamics();
        let mut species = Vec::with_capacity(model.components().len());
        let mut sum = 0.0;
        let mut correction = 0.0;
        for (component, &fraction) in model.components().iter().zip(model.mole_composition().fractions()) {
            poll(&mut cancelled)?;
            let alpha = *self.coefficients.get(component.species())
                .ok_or(invalid("thermal accommodation is missing an active original gas species"))?;
            let nonsticking = if component.species() == &binding.species {
                1.0 - binding.accommodation_coefficient
            } else { 1.0 };
            let mut accommodated_incident_mass_kg = 0.0;
            let mut specific_energy_change_j_kg = 0.0;
            let mut heat_into_gas_j = 0.0;
            if alpha != 0.0 && nonsticking != 0.0 && exposure != 0.0 {
                let molar_mass = component.molar_mass().value();
                let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass;
                let denominator = fs_math::det::sqrt(2.0 * core::f64::consts::PI)
                    * fs_math::det::sqrt(r) * fs_math::det::sqrt(gas_temperature_k);
                let pressure = fraction * gas.state().state().pressure;
                if !positive(molar_mass) || !positive(r) || !positive(denominator) || !positive(pressure) {
                    return Err(invalid("original molecular incident state is unrepresentable"));
                }
                accommodated_incident_mass_kg = (pressure / denominator) * (alpha * nonsticking) * exposure;
                let enthalpy_change = component.molar_enthalpy_change(
                    Temperature::new(gas_temperature_k), Temperature::new(surface_temperature_k),
                ).map_err(|_| invalid("original reflected-species caloric endpoints refused"))?.value() / molar_mass;
                let delta_t = surface_temperature_k - gas_temperature_k;
                specific_energy_change_j_kg = enthalpy_change - 0.5 * r * delta_t;
                heat_into_gas_j = accommodated_incident_mass_kg * specific_energy_change_j_kg;
                if !positive(accommodated_incident_mass_kg)
                    || !specific_energy_change_j_kg.is_finite() || !heat_into_gas_j.is_finite()
                    || (delta_t != 0.0 && !positive(specific_energy_change_j_kg / delta_t))
                    || (specific_energy_change_j_kg != 0.0 && heat_into_gas_j == 0.0)
                {
                    return Err(invalid("reflected molecular heat is unrepresentable or nonpassive"));
                }
            }
            let (next, error) = fs_math::eft::two_sum(sum, heat_into_gas_j);
            sum = next;
            correction += error;
            if !sum.is_finite() || !correction.is_finite() {
                return Err(invalid("reflected mixture heat is unrepresentable"));
            }
            species.push(ReflectedSpeciesHeat {
                species: component.species().clone(),
                accommodated_incident_mass_kg,
                specific_energy_change_j_kg,
                heat_into_gas_j,
            });
        }
        let heat_into_gas_j = sum + correction;
        if !heat_into_gas_j.is_finite() {
            return Err(invalid("reflected mixture heat is unrepresentable"));
        }
        species.sort_unstable_by(|left, right| left.species.cmp(&right.species));
        poll(&mut cancelled)?;
        Ok(ReflectedHeatStep { species, heat_into_gas_j })
    }
}

/// One original species' reflection account. None of this mass is withdrawn.
#[derive(Clone, Debug, PartialEq)]
pub struct ReflectedSpeciesHeat {
    /// Original canonical gas species identity.
    pub species: SpeciesId,
    /// Incident encounter mass weighted by nonsticking fraction and alpha_E.
    /// Molecules return to the gas; repeated encounters are not new inventory.
    pub accommodated_incident_mass_kg: f64,
    /// e(T_surface)-e(T_gas), independent of the formation-energy reference.
    pub specific_energy_change_j_kg: f64,
    /// Signed heat delivered to gas [J]; already included in the joint exchange.
    pub heat_into_gas_j: f64,
}

/// Canonical reflected-species contributions, separate from phase mass streams.
#[derive(Clone, Debug, PartialEq)]
pub struct ReflectedHeatStep {
    /// Active species in canonical identity order, including explicit zeros.
    pub species: Vec<ReflectedSpeciesHeat>,
    /// Signed heat delivered to gas [J]; already included in the joint exchange.
    pub heat_into_gas_j: f64,
}

/// One atomic original-owner endpoint plus its separately identified reflection.
#[derive(Clone, Debug, PartialEq)]
pub struct AccommodatingSurfaceCandidate<'catalog> {
    exchange: GasSurfaceCandidate<'catalog>,
    reflection: ReflectedHeatStep,
}
impl<'catalog> AccommodatingSurfaceCandidate<'catalog> {
    /// Combined phase/reflection endpoint and total energy accounts.
    #[must_use]
    pub const fn exchange(&self) -> &GasSurfaceCandidate<'catalog> { &self.exchange }
    /// This heat is INCLUDED in exchange().step().energy_into_gas_j, not added to it.
    #[must_use]
    pub const fn reflection(&self) -> &ReflectedHeatStep { &self.reflection }
    /// Consume both the accepted endpoint and its reflection subaccount.
    #[must_use]
    pub fn into_parts(self) -> (GasSurfaceCandidate<'catalog>, ReflectedHeatStep) {
        (self.exchange, self.reflection)
    }
}

/// Sublimation/deposition and reflected heat at one jointly reconstructed endpoint.
/// No intermediate phase-only or heat-only state must be source-admissible.
/// The same original temperatures sample both terms; explicit time refinement
/// remains necessary, especially for small solid/gas heat capacities.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface_with_accommodation<'catalog, L: SublimationLaw + ?Sized>(
    law: &L,
    binding: &GasSurfaceBinding,
    solid: SolidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, SublimationError>,
    accommodation: &ThermalAccommodation,
    duration_s: f64,
    controls: GasSurfaceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AccommodatingSurfaceCandidate<'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_controls(controls, duration_s)?;
    let gas_snapshot = species_snapshot(gas, catalog, &binding.species)?;
    let snapshot = snapshot_surface(law, binding, solid, gas_snapshot, catalog,
        &temperature_from_enthalpy, duration_s, controls)?;
    admit_surface_inventory(snapshot, duration_s, controls)?;
    let reflection = accommodation.evaluate(binding, gas, snapshot.temperature_k,
        solid.area_m2, duration_s, &mut cancelled)?;
    let prepared = prepare_surface_with_heat(law, binding, snapshot, catalog,
        &temperature_from_enthalpy, controls, reflection.heat_into_gas_j)?;
    let mut rows = BTreeMap::new();
    for reflected in &reflection.species {
        rows.insert(reflected.species.clone(), SpeciesBoundaryTransfer {
            species: reflected.species.clone(), mass_into_gas_kg: 0.0,
            energy_into_gas_j: reflected.heat_into_gas_j,
        });
    }
    let difference = snapshot.incoming_specific_energy_j_kg
        .map_or(0.0, |energy| snapshot.outgoing_specific_energy_j_kg - energy);
    let phase_energy = snapshot.transferred_mass_kg * snapshot.outgoing_specific_energy_j_kg
        + snapshot.incoming_mass_kg * difference;
    let row = rows.entry(binding.species.clone()).or_insert_with(|| SpeciesBoundaryTransfer {
        species: binding.species.clone(), mass_into_gas_kg: 0.0, energy_into_gas_j: 0.0,
    });
    row.mass_into_gas_kg = snapshot.transferred_mass_kg;
    row.energy_into_gas_j += phase_energy;
    let rows: Vec<_> = rows.into_values().collect();
    let endpoint = gas.exchange_species_boundary(catalog, &rows, controls.gas, &mut cancelled)?;
    accommodation.admit_temperatures(prepared.surface_temperature_k[1], endpoint.gas().state().state().temperature)?;
    let exchange = finish_surface(binding, solid, gas, gas_snapshot, prepared, endpoint, controls, cancelled)?;
    Ok(AccommodatingSurfaceCandidate { exchange, reflection })
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod sharing_tests {
    use super::*;
    use crate::sublimation::gas::tests as fixture;

    #[test]
    fn admitted_tables_are_shared_immutable_owned_and_thread_safe() {
        fn send_sync<T: Send + Sync>() {}
        send_sync::<ThermalAccommodation>();
        let catalog = fixture::catalog(0.0);
        let carrier = SpeciesId::new("A").unwrap();
        let vapor = fixture::species();
        let mut input = vec![(carrier.clone(), 0.5), (vapor.clone(), 0.6)];
        let original = ThermalAccommodation::new(&catalog, &input, [200.0, 900.0]).unwrap();
        let independent = ThermalAccommodation::new(&catalog, &input, [200.0, 900.0]).unwrap();
        let handles: Vec<_> = (0..32).map(|_| original.clone()).collect();
        assert_eq!(original, independent);
        assert!(!std::ptr::eq(original.coefficients(), independent.coefficients()));
        for handle in &handles {
            assert!(std::ptr::eq(original.coefficients(), handle.coefficients()));
            assert_eq!(handle.temperature_range_k(), [200.0, 900.0]);
        }
        input[0].1 = 0.0;
        drop(original);
        for handle in handles {
            assert_eq!(handle.coefficients()[&carrier], 0.5);
            assert_eq!(handle.coefficients()[&vapor], 0.6);
            assert_eq!(handle, independent);
        }
    }

    #[test]
    fn shared_coefficients_preserve_original_phase_and_reflection_endpoints() {
        let catalog = fixture::catalog(0.0);
        let binding = fixture::binding(0.0);
        let chart = fixture::chart(&catalog, &binding);
        let solid = fixture::solid(400.0, chart);
        let gas = fixture::gas(&catalog, 350.0, 800.0);
        let original = ThermalAccommodation::new(&catalog, &[
            (SpeciesId::new("A").unwrap(), 0.5), (fixture::species(), 0.6),
        ], [200.0, 900.0]).unwrap();
        let shared = original.clone();
        let run = |accommodation: &ThermalAccommodation| {
            advance_surface_with_accommodation(&fixture::phase(), &binding, solid, &gas,
                &catalog, |h| fixture::temperature(h, chart), accommodation, 0.001,
                fixture::controls(), || false).unwrap()
        };
        let expected = run(&original);
        drop(original);
        let accepted = run(&shared);
        assert_eq!(accepted, expected);
        assert!(accepted.reflection().heat_into_gas_j > 0.0);
        fixture::near(accepted.exchange().step().energy_residual_j, 0.0, 1e-7);
        fixture::near(accepted.exchange().step().mass_residual_kg, 0.0, 1e-14);
    }
}
