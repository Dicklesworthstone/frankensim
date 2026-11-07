//! Closed gas/coating phase equilibrium on the original caloric sources.
//!
//! This is a thermodynamic endpoint, not a transient molecular step or a
//! nucleation/extinction clock. A positive exposed area, contact conductance and
//! molecular accommodation explicitly permit exchange. Both condensed bodies
//! and the gas reach one temperature on a caller-declared increasing branch.
//!
//! At that temperature the volatile inventory is partitioned between saturated
//! vapor and deposit; insufficient total volatile mass leaves an exactly bare,
//! undersaturated support. The original gas owner reconstructs composition,
//! pressure, transport and acoustics at most once, after the scalar root.
//!
//! The gas volume is its fixed accessible volume, excluding fixed supports.
//! An explicit dilute deposit-volume ceiling bounds omitted displacement, not
//! pressure-work error. No morphology, pore model, metastability, chemical
//! reaction or kinetic equilibration time is inferred.

use super::super::super::Total;
use super::{
    ClosedGasParcel, GasResistiveSurface, GasSurfaceError, ResistiveSurfaceState,
    ThermochemicalGasModel,
};
use crate::sublimation::gas::{
    compatibility, phase_at, poll, positive, validate_surface,
};
use crate::thermochemical_gas::energy::{
    GasEnergyError,
    parcel::transfer::species::boundary::admit_catalog,
};
use fs_thermal::sublimation::{SublimationError, SublimationProperties};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
use std::cell::Cell;
use std::collections::BTreeMap;

/// The original support-enthalpy bracket selects one increasing equilibrium
/// branch. A nonstationary solve admits both bracket endpoints through the
/// original caloric and phase sources. An already exact stationary endpoint
/// needs no root and does not probe unused bracket endpoints. Tolerances are
/// absolute.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasResistiveEquilibriumSolve {
    pub substrate_enthalpy_bracket_j_kg: [f64; 2],
    pub absolute_energy_tolerance_j: f64,
    pub absolute_mass_tolerance_kg: f64,
    pub absolute_temperature_tolerance_k: f64,
    pub absolute_pressure_tolerance_pa: f64,
    pub latent_compatibility_tolerance_j_kg: f64,
    /// Deposit volume / fixed accessible gas volume, in (0, 0.01].
    pub maximum_deposit_volume_fraction: f64,
    /// Includes the two bracket endpoints; in 2..=4096.
    pub maximum_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasResistiveEquilibriumError {
    Invalid(&'static str),
    Source(GasSurfaceError),
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        evaluations: usize,
        residual_j: f64,
    },
}
impl core::fmt::Display for GasResistiveEquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "closed coating equilibrium: {self:?}")
    }
}
impl std::error::Error for GasResistiveEquilibriumError {}
impl From<GasSurfaceError> for GasResistiveEquilibriumError {
    fn from(error: GasSurfaceError) -> Self {
        Self::Source(error)
    }
}
impl From<GasEnergyError> for GasResistiveEquilibriumError {
    fn from(error: GasEnergyError) -> Self {
        Self::Source(error.into())
    }
}
impl From<SublimationError> for GasResistiveEquilibriumError {
    fn from(error: SublimationError) -> Self {
        Self::Source(error.into())
    }
}

/// No tolerance-based coating deletion: the all-vapor branch is selected only
/// when the complete volatile inventory fits below the original saturation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasResistiveEquilibriumPhase {
    SaturatedCoating,
    AllVapor,
}

/// Original-chart endpoint balances. The initial temperature order is
/// [support, deposit, gas]. No field represents a finite-duration heat rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasResistiveEquilibriumStep {
    pub id: u64,
    pub initial_temperatures_k: [f64; 3],
    pub temperature_k: f64,
    pub phase: GasResistiveEquilibriumPhase,
    pub vapor_partial_pressure_pa: [f64; 2],
    pub saturation_pressure_pa: f64,
    /// Positive sublimation, negative deposition.
    pub transferred_mass_kg: f64,
    pub substrate_energy_change_j: f64,
    /// Includes the fixed deposit-to-gas enthalpy datum conversion.
    pub deposit_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    /// Actual caloric changes in a fixed source-relative energy chart.
    /// Inventory representation is admitted separately by the mass bounds.
    pub energy_residual_j: f64,
    /// Included in the energy residual; never an extra physical heat source.
    pub caloric_encoding_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub maximum_temperature_residual_k: f64,
    /// Actual vapor pressure minus saturation: zero for a positive deposit,
    /// nonpositive on the all-vapor branch, within the explicit pressure budget.
    pub saturation_residual_pa: f64,
    pub deposit_volume_fraction: [f64; 2],
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GasResistiveEquilibriumCandidate<'catalog> {
    state: ResistiveSurfaceState,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    step: GasResistiveEquilibriumStep,
}
impl<'catalog> GasResistiveEquilibriumCandidate<'catalog> {
    #[must_use]
    pub const fn state(&self) -> ResistiveSurfaceState {
        self.state
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    #[must_use]
    pub const fn step(&self) -> &GasResistiveEquilibriumStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ResistiveSurfaceState,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasResistiveEquilibriumStep,
    ) {
        (self.state, self.gas, self.catalog, self.step)
    }
}

fn invalid(message: &'static str) -> GasResistiveEquilibriumError {
    GasResistiveEquilibriumError::Invalid(message)
}

fn volume_fraction(
    mass_kg: f64,
    phase: SublimationProperties,
    gas_volume_m3: f64,
    maximum: f64,
) -> Result<f64, GasResistiveEquilibriumError> {
    let fraction = (mass_kg / phase.solid_density_kg_m3) / gas_volume_m3;
    if !fraction.is_finite() || fraction < 0.0 || fraction > maximum {
        return Err(invalid(
            "coating exceeds the declared fixed-gas-volume approximation",
        ));
    }
    Ok(fraction)
}

/// Use only original partial enthalpy INCREMENTS; formation terms cancel before
/// evaluating Cv storage. At most the original bounded species roster is visited.
fn gas_sensible(
    catalog: &ThermochemicalGasModel,
    before: &BTreeMap<SpeciesId, f64>,
    volatile: &SpeciesId,
    vapor_mass_kg: f64,
    actual: Option<&BTreeMap<SpeciesId, f64>>,
    from_k: f64,
    to_k: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasResistiveEquilibriumError> {
    let mut sum = Total::default();
    for source in catalog.thermodynamics().components() {
        poll(cancelled)?;
        let id = source.species();
        let mass = match actual {
            Some(inventory) => inventory.get(id).copied().unwrap_or(0.0),
            None if id == volatile => vapor_mass_kg,
            None => before.get(id).copied().unwrap_or(0.0),
        };
        if mass == 0.0 {
            continue;
        }
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / source.molar_mass().value();
        let dh = catalog.species_specific_enthalpy_change_j_kg(id, from_k, to_k)?;
        sum.add(mass * (dh - r * (to_k - from_k)))?;
    }
    Ok(sum.value()?)
}

#[derive(Clone, Copy)]
struct Trial {
    support_h: f64,
    deposit_h: f64,
    temperature_k: f64,
    temperature_residual_k: f64,
    deposit_mass_kg: f64,
    vapor_mass_kg: f64,
    dm: f64,
    phase: SublimationProperties,
    volume_fraction: f64,
    support_energy_j: f64,
    deposit_sensible_j: f64,
    encoding_energy_j: f64,
    residual_j: f64,
}

/// Find the isolated equilibrium of ONE nonvolatile support, its volatile
/// coating and one finite multicomponent original gas.
///
/// Positive exposed area, scalar contact conductance and molecular
/// accommodation are necessary exchange permissions; an authored insulating or
/// blocked interface refuses. Their magnitudes do not select a kinetic time.
/// Initial body temperatures may differ. A bare input carries its original
/// deposit caloric reference at the support temperature, with no fictitious mass.
///
/// The caller declares a continuous increasing support-enthalpy bracket.
/// Each trial derives deposit h from original vapor h(T)-latent(T)-offset and
/// admits the result through the ORIGINAL deposit inverse chart. Source failures
/// are never bypassed with an extrapolated chart. A non-increasing sampled
/// equilibrium branch refuses; this is not a global free-energy minimizer.
///
/// The conserved energy change uses original NASA enthalpy increments, original
/// latent heat and the fixed binding offset. Formation constants cancel before
/// the scalar residual; no subtraction of absolute gas/solid energy totals is
/// used for acceptance. Actual caloric encoding and species rounding retain
/// separate explicit checks. The original gas is reconstructed only after
/// convergence, retaining zero individual species but refusing total vacuum.
///
/// The endpoint does not update exposed geometry, locate nucleation/extinction
/// times, or model a moving gas volume. A later transient owner must admit its
/// own interface and geometry at this returned state.
pub fn equilibrate_resistive_chamber<'catalog>(
    surface: &GasResistiveSurface<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    controls: GasResistiveEquilibriumSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveEquilibriumCandidate<'catalog>, GasResistiveEquilibriumError> {
    poll(&mut cancelled)?;
    let [lower_h, upper_h] = controls.substrate_enthalpy_bracket_j_kg;
    if !lower_h.is_finite()
        || !upper_h.is_finite()
        || lower_h >= upper_h
        || !positive(controls.absolute_energy_tolerance_j)
        || !positive(controls.absolute_mass_tolerance_kg)
        || !positive(controls.absolute_temperature_tolerance_k)
        || !positive(controls.absolute_pressure_tolerance_pa)
        || !positive(controls.latent_compatibility_tolerance_j_kg)
        || !positive(controls.maximum_deposit_volume_fraction)
        || controls.maximum_deposit_volume_fraction > 0.01
        || !(2..=4096).contains(&controls.maximum_evaluations)
    {
        return Err(invalid("finite equilibrium bracket and positive bounded solve controls required"));
    }
    let initial = surface.state;
    validate_surface(surface.binding, initial.deposit)?;
    if !positive(initial.substrate_mass_kg)
        || !initial.substrate_specific_enthalpy_j_kg.is_finite()
        || !positive(surface.conductance_w_k)
        || !positive(initial.deposit.area_m2)
        || !positive(surface.binding.accommodation_coefficient)
    {
        return Err(invalid("equilibrium requires positive support storage and permitted thermal/phase exchange"));
    }
    admit_catalog(gas, catalog)?;
    let before = gas.species_inventory_kg()?;
    let species = &surface.binding.species;
    let source = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|source| source.species() == species)
        .ok_or(invalid("coating equilibrium requires its retained original vapor species"))?;
    let molar_mass = source.molar_mass().value();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass;
    let vapor_before = before.get(species).copied().unwrap_or(0.0);
    let total_volatile = initial.deposit.mass_kg + vapor_before;
    if !positive(total_volatile) || !positive(r) {
        return Err(invalid("positive representable total volatile inventory required"));
    }
    let support_t0 = (surface.substrate_temperature)(initial.substrate_specific_enthalpy_j_kg)?;
    let deposit_t0 = (surface.deposit_temperature)(initial.deposit.specific_enthalpy_j_kg)?;
    let gas_t0 = gas.state().state().temperature;
    if !positive(support_t0) || !positive(deposit_t0) {
        return Err(invalid("original condensed inverse returned an invalid temperature"));
    }
    if initial.deposit.mass_kg == 0.0
        && (deposit_t0 - support_t0).abs() > controls.absolute_temperature_tolerance_k
    {
        return Err(invalid("bare deposit caloric reference must match its support temperature"));
    }
    let initial_phase = phase_at(surface.law, deposit_t0, molar_mass)?;
    let initial_vapor_h = catalog.species_specific_enthalpy_j_kg(species, deposit_t0)?;
    let alignment0 = compatibility(
        initial_vapor_h,
        initial.deposit.specific_enthalpy_j_kg,
        surface.binding.solid_enthalpy_offset_j_kg,
        initial_phase.latent_heat_j_kg,
        deposit_t0,
        controls.latent_compatibility_tolerance_j_kg,
    )?;
    let initial_volume = volume_fraction(
        initial.deposit.mass_kg,
        initial_phase,
        gas.volume_m3(),
        controls.maximum_deposit_volume_fraction,
    )?;
    let vapor_h0 = catalog.species_specific_enthalpy_j_kg(species, gas_t0)?;
    // u_v(T_g0) - [h_d0 + offset], without subtracting absolute formation terms.
    let carried_relative = catalog.species_specific_enthalpy_change_j_kg(
        species, deposit_t0, gas_t0,
    )? + initial_phase.latent_heat_j_kg + alignment0 - r * gas_t0;
    if !carried_relative.is_finite() {
        return Err(invalid("equilibrium relative phase energy is unrepresentable"));
    }
    let initial_pressure = gas.species_partial_pressure_pa(species)?;
    let exact_stationary = support_t0 == deposit_t0
        && support_t0 == gas_t0
        && initial.substrate_specific_enthalpy_j_kg >= lower_h
        && initial.substrate_specific_enthalpy_j_kg <= upper_h
        && if initial.deposit.mass_kg > 0.0 {
            initial_pressure == initial_phase.saturation_pressure_pa
        } else {
            initial_pressure <= initial_phase.saturation_pressure_pa
        };
    if exact_stationary {
        poll(&mut cancelled)?;
        return Ok(GasResistiveEquilibriumCandidate {
            state: initial,
            gas: gas.clone(),
            catalog,
            step: GasResistiveEquilibriumStep {
                id: surface.id,
                initial_temperatures_k: [support_t0, deposit_t0, gas_t0],
                temperature_k: gas_t0,
                phase: if initial.deposit.mass_kg > 0.0 {
                    GasResistiveEquilibriumPhase::SaturatedCoating
                } else {
                    GasResistiveEquilibriumPhase::AllVapor
                },
                vapor_partial_pressure_pa: [initial_pressure; 2],
                saturation_pressure_pa: initial_phase.saturation_pressure_pa,
                transferred_mass_kg: 0.0,
                substrate_energy_change_j: 0.0,
                deposit_energy_change_j: 0.0,
                gas_energy_change_j: 0.0,
                energy_residual_j: 0.0,
                caloric_encoding_residual_j: 0.0,
                mass_residual_kg: 0.0,
                maximum_species_mass_residual_kg: 0.0,
                maximum_temperature_residual_k: 0.0,
                saturation_residual_pa: initial_pressure - initial_phase.saturation_pressure_pa,
                deposit_volume_fraction: [initial_volume; 2],
                evaluations: 0,
            },
        });
    }
    let evaluations = Cell::new(0usize);
    let last_residual = Cell::new(f64::INFINITY);
    let mut evaluate = |support_h: f64| -> Result<Trial, GasResistiveEquilibriumError> {
        poll(&mut cancelled)?;
        if evaluations.get() >= controls.maximum_evaluations {
            return Err(GasResistiveEquilibriumError::NotConverged {
                evaluations: evaluations.get(),
                residual_j: last_residual.get(),
            });
        }
        evaluations.set(evaluations.get() + 1);
        let t = (surface.substrate_temperature)(support_h)?;
        let phase = phase_at(surface.law, t, molar_mass)?;
        let vapor_h = catalog.species_specific_enthalpy_j_kg(species, t)?;
        let mut forward = Total::default();
        forward.add(vapor_h)?;
        forward.add(-phase.latent_heat_j_kg)?;
        forward.add(-surface.binding.solid_enthalpy_offset_j_kg)?;
        let deposit_h = forward.value()?;
        let admitted_t = (surface.deposit_temperature)(deposit_h)?;
        let temperature_residual_k = (admitted_t - t).abs();
        if !positive(admitted_t)
            || temperature_residual_k > controls.absolute_temperature_tolerance_k
        {
            return Err(invalid("original deposit caloric inverse disagrees with vapor-minus-latent enthalpy"));
        }
        let alignment = compatibility(
            vapor_h, deposit_h, surface.binding.solid_enthalpy_offset_j_kg,
            phase.latent_heat_j_kg, t, controls.latent_compatibility_tolerance_j_kg,
        )?;
        let all_vapor_pressure = (total_volatile / gas.volume_m3()) * r * t;
        if !positive(all_vapor_pressure) {
            return Err(invalid("all-vapor equilibrium pressure is unrepresentable"));
        }
        let (deposit_mass_kg, vapor_mass_kg) =
            if all_vapor_pressure <= phase.saturation_pressure_pa {
                (0.0, total_volatile)
            } else {
                let vapor = (phase.saturation_pressure_pa / (r * t)) * gas.volume_m3();
                let deposit = total_volatile - vapor;
                if !vapor.is_finite()
                    || vapor < 0.0
                    || (phase.saturation_pressure_pa > 0.0 && vapor == 0.0)
                    || !positive(deposit)
                {
                    return Err(invalid("saturated coating partition is unrepresentable"));
                }
                (deposit, vapor)
            };
        let dm = vapor_mass_kg - vapor_before;
        let deposit_change = deposit_mass_kg - initial.deposit.mass_kg;
        if !dm.is_finite()
            || (dm != 0.0 && deposit_change == 0.0)
            || (deposit_change != 0.0 && dm == 0.0)
            || (deposit_change + dm).abs() > controls.absolute_mass_tolerance_kg
        {
            return Err(invalid("equilibrium cannot represent the opposite phase mass changes"));
        }
        let fraction = volume_fraction(
            deposit_mass_kg, phase, gas.volume_m3(),
            controls.maximum_deposit_volume_fraction,
        )?;
        let source_change = catalog.species_specific_enthalpy_change_j_kg(
            species, deposit_t0, t,
        )? - (phase.latent_heat_j_kg - initial_phase.latent_heat_j_kg)
            + alignment0 - alignment;
        let encoding_energy_j = deposit_mass_kg
            * ((deposit_h - initial.deposit.specific_enthalpy_j_kg) - source_change);
        if !encoding_energy_j.is_finite() {
            return Err(invalid("deposit caloric encoding is unrepresentable"));
        }
        if encoding_energy_j.abs() > controls.absolute_energy_tolerance_j {
            return Err(GasSurfaceError::EnergyBalance {
                residual_j: encoding_energy_j,
                tolerance_j: controls.absolute_energy_tolerance_j,
            }.into());
        }
        let deposit_sensible_j = deposit_mass_kg * source_change;
        let support_energy_j = initial.substrate_mass_kg
            * (support_h - initial.substrate_specific_enthalpy_j_kg);
        let sensible = gas_sensible(
            catalog, &before, species, vapor_mass_kg, None, gas_t0, t, &mut cancelled,
        )?;
        let mut residual = Total::default();
        for term in [
            support_energy_j, deposit_sensible_j, encoding_energy_j,
            sensible, dm * carried_relative,
        ] {
            residual.add(term)?;
        }
        let residual_j = residual.value()?;
        last_residual.set(residual_j);
        Ok(Trial {
            support_h, deposit_h, temperature_k: t, temperature_residual_k,
            deposit_mass_kg, vapor_mass_kg, dm, phase, volume_fraction: fraction,
            support_energy_j, deposit_sensible_j, encoding_energy_j, residual_j,
        })
    };
    let mut lower = evaluate(lower_h)?;
    let mut upper = evaluate(upper_h)?;
    let tolerance = controls.absolute_energy_tolerance_j;
    if lower.temperature_k > upper.temperature_k
        || lower.residual_j - upper.residual_j > tolerance
    {
        return Err(invalid("equilibrium bracket is not nondecreasing"));
    }
    let accepted = if lower.residual_j.abs() <= tolerance {
        lower
    } else if upper.residual_j.abs() <= tolerance {
        upper
    } else {
        if lower.residual_j > 0.0 || upper.residual_j < 0.0 {
            return Err(GasResistiveEquilibriumError::OutsideBracket {
                lower_residual_j: lower.residual_j,
                upper_residual_j: upper.residual_j,
            });
        }
        loop {
            let h = 0.5 * lower.support_h + 0.5 * upper.support_h;
            if !h.is_finite() || h <= lower.support_h || h >= upper.support_h {
                return Err(GasResistiveEquilibriumError::NotConverged {
                    evaluations: evaluations.get(),
                    residual_j: last_residual.get(),
                });
            }
            let trial = evaluate(h)?;
            if trial.temperature_k < lower.temperature_k
                || trial.temperature_k > upper.temperature_k
                || lower.residual_j - trial.residual_j > tolerance
                || trial.residual_j - upper.residual_j > tolerance
            {
                return Err(invalid("equilibrium left the declared increasing thermodynamic branch"));
            }
            if trial.residual_j.abs() <= tolerance {
                break trial;
            }
            if trial.residual_j < 0.0 {
                lower = trial;
            } else {
                upper = trial;
            }
        }
    };
    drop(evaluate);
    poll(&mut cancelled)?;
    let mut target = before.clone();
    target.insert(species.clone(), accepted.vapor_mass_kg);
    let endpoint = gas.at_species_inventory_temperature(catalog, &target, accepted.temperature_k)?;
    poll(&mut cancelled)?;
    let actual = endpoint.species_inventory_kg()?;
    let mut maximum_species_mass_residual_kg: f64 = 0.0;
    for source in catalog.thermodynamics().components() {
        poll(&mut cancelled)?;
        let id = source.species();
        let expected = target.get(id).copied().unwrap_or(0.0);
        let stored = actual.get(id).copied().unwrap_or(0.0);
        if (expected > 0.0) != (stored > 0.0) {
            return Err(invalid("equilibrium reconstruction lost an active species"));
        }
        maximum_species_mass_residual_kg =
            maximum_species_mass_residual_kg.max((stored - expected).abs());
    }
    let mass_residual_kg = (accepted.deposit_mass_kg - initial.deposit.mass_kg)
        + (endpoint.mass_kg() - gas.mass_kg());
    let actual_vapor = actual.get(species).copied().unwrap_or(0.0);
    maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(
        ((actual_vapor - vapor_before) - accepted.dm).abs(),
    );
    if !mass_residual_kg.is_finite()
        || mass_residual_kg.abs() > controls.absolute_mass_tolerance_kg
        || maximum_species_mass_residual_kg > controls.absolute_mass_tolerance_kg
    {
        return Err(invalid("original equilibrium endpoint failed species/total mass balance"));
    }
    let temperature = endpoint.state().state().temperature;
    let maximum_temperature_residual_k = accepted.temperature_residual_k
        .max((temperature - accepted.temperature_k).abs());
    let pressure = endpoint.species_partial_pressure_pa(species)?;
    let saturation_residual_pa = pressure - accepted.phase.saturation_pressure_pa;
    if maximum_temperature_residual_k > controls.absolute_temperature_tolerance_k
        || !saturation_residual_pa.is_finite()
        || if accepted.deposit_mass_kg > 0.0 {
            saturation_residual_pa.abs() > controls.absolute_pressure_tolerance_pa
        } else {
            saturation_residual_pa > controls.absolute_pressure_tolerance_pa
        }
    {
        return Err(invalid("reconstructed equilibrium failed common-temperature or phase complementarity"));
    }
    let sensible = gas_sensible(
        catalog, &before, species, accepted.vapor_mass_kg, Some(&actual),
        gas_t0, temperature, &mut cancelled,
    )?;
    let mut residual = Total::default();
    for term in [
        accepted.support_energy_j, accepted.deposit_sensible_j, accepted.encoding_energy_j,
        sensible, accepted.dm * carried_relative,
    ] {
        residual.add(term)?;
    }
    let energy_residual_j = residual.value()?;
    if energy_residual_j.abs() > controls.absolute_energy_tolerance_j {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: energy_residual_j,
            tolerance_j: controls.absolute_energy_tolerance_j,
        }.into());
    }
    let mut deposited_energy = Total::default();
    deposited_energy.add(accepted.deposit_sensible_j)?;
    deposited_energy.add(accepted.encoding_energy_j)?;
    deposited_energy.add(-accepted.dm
        * (initial.deposit.specific_enthalpy_j_kg
            + surface.binding.solid_enthalpy_offset_j_kg))?;
    let mut gas_energy = Total::default();
    gas_energy.add(sensible)?;
    gas_energy.add(accepted.dm * (vapor_h0 - r * gas_t0))?;
    let state = ResistiveSurfaceState {
        substrate_specific_enthalpy_j_kg: accepted.support_h,
        deposit: fs_thermal::sublimation::SolidSurfaceState {
            mass_kg: accepted.deposit_mass_kg,
            specific_enthalpy_j_kg: accepted.deposit_h,
            ..initial.deposit
        },
        ..initial
    };
    let step = GasResistiveEquilibriumStep {
        id: surface.id,
        initial_temperatures_k: [support_t0, deposit_t0, gas_t0],
        temperature_k: accepted.temperature_k,
        phase: if accepted.deposit_mass_kg > 0.0 {
            GasResistiveEquilibriumPhase::SaturatedCoating
        } else {
            GasResistiveEquilibriumPhase::AllVapor
        },
        vapor_partial_pressure_pa: [initial_pressure, pressure],
        saturation_pressure_pa: accepted.phase.saturation_pressure_pa,
        transferred_mass_kg: accepted.dm,
        substrate_energy_change_j: accepted.support_energy_j,
        deposit_energy_change_j: deposited_energy.value()?,
        gas_energy_change_j: gas_energy.value()?,
        energy_residual_j,
        caloric_encoding_residual_j: accepted.encoding_energy_j,
        mass_residual_kg,
        maximum_species_mass_residual_kg,
        maximum_temperature_residual_k,
        saturation_residual_pa,
        deposit_volume_fraction: [initial_volume, accepted.volume_fraction],
        evaluations: evaluations.get(),
    };
    poll(&mut cancelled)?;
    Ok(GasResistiveEquilibriumCandidate { state, gas: endpoint, catalog, step })
}

#[cfg(test)]
mod tests;
