//! Conservative first filling of an exactly evacuated finite volume.
//!
//! An empty volume has geometry and auxiliary material conditions, but no gas
//! temperature, pressure, composition or fictitious seed mass. Positive received
//! species and source-relative internal energy construct ONE ordinary gas state
//! through the original composition and caloric owners.
//!
//! The endpoint is a thermalized, well-mixed, zero-momentum gas. Its physical
//! owner must justify that equilibration assumption and account for opposite
//! source storage, wall heat/work, and any captured kinetic energy. This local
//! initialization supplies no flow rate or collisionless startup dynamics.

use super::{
    ClosedGasParcel,
    transfer::{inventory, species::model_for_inventory},
};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergyError, GasEnergySolve, IsochoricGasParcel, poll, positive, updated_point},
    failure,
};
use fs_matdb::QueryPoint;
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{SpeciesId, mixture::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1};
use std::collections::BTreeMap;

/// A finite rigid volume containing exactly zero gas mass and internal energy.
/// No temperature or pressure is assigned to its absent gas phase.
#[derive(Clone, Debug, PartialEq)]
pub struct EvacuatedGasVolume {
    volume_m3: f64,
    context: QueryPoint,
}

/// Independent accepted-inventory and original-caloric first-fill balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvacuatedGasFillStep {
    /// Compensated received mass, rounded once to the parcel's stored mass [kg].
    pub received_mass_kg: f64,
    /// Numerical energy-chart origin, not a temperature of the original vacuum.
    pub reference_temperature_k: f64,
    /// Requested `U - sum_i m_i*u_i(T_reference)` [J].
    pub reference_relative_energy_j: f64,
    /// Actual represented species masses times original-law Delta u [J].
    pub stored_relative_energy_j: f64,
    /// Sum of represented species masses minus the received inventory [kg].
    pub mass_residual_kg: f64,
    /// Original EOS density times volume minus the received inventory [kg].
    pub density_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub iterations: usize,
}

/// A complete ordinary gas endpoint, or no published state on any refusal.
#[derive(Clone, Debug, PartialEq)]
pub struct EvacuatedGasFillCandidate {
    gas: ClosedGasParcel,
    step: EvacuatedGasFillStep,
}
impl EvacuatedGasFillCandidate {
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn step(&self) -> &EvacuatedGasFillStep {
        &self.step
    }
    pub fn into_parts(self) -> (ClosedGasParcel, EvacuatedGasFillStep) {
        (self.gas, self.step)
    }
}

fn add(total: &mut Dd, value: Dd) -> Result<(), GasEnergyError> {
    *total = *total + value;
    if !value.hi.is_finite()
        || !value.lo.is_finite()
        || !total.hi.is_finite()
        || !total.lo.is_finite()
    {
        return Err(GasEnergyError::Invalid(
            "evacuated-volume inventory or energy is unrepresentable",
        ));
    }
    Ok(())
}

impl EvacuatedGasVolume {
    /// Retain only positive geometry and the receiving volume's auxiliary axes.
    /// Even zero-valued or untyped temperature/pressure coordinates are refused.
    pub fn new(volume_m3: f64, context: QueryPoint) -> Result<Self, GasEnergyError> {
        if !positive(volume_m3)
            || context.axes().contains_key("temperature")
            || context.axes().contains_key("pressure")
        {
            return Err(GasEnergyError::Invalid(
                "an evacuated volume needs positive geometry and no gas temperature/pressure axes",
            ));
        }
        Ok(Self { volume_m3, context })
    }

    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }

    /// Construct a thermalized gas from explicit conserved first-fill data.
    ///
    /// `reference_relative_energy_j` is the received internal energy minus
    /// `sum_i m_i*u_i(reference_temperature_k)`, evaluated directly in the exact
    /// original species' caloric charts. Never obtain it by subtracting rounded
    /// absolute formation-energy totals. The reference temperature is a numerical
    /// coordinate; it is not assigned to the empty volume or a seed gas.
    ///
    /// Every named inventory must belong to the retained catalog. Exact zero
    /// species are omitted by its existing composition owner; positive traces
    /// must remain representable. The accepted total mass is stored directly,
    /// while represented species, EOS mass and original-law energy increments
    /// are independently checked against fixed extensive budgets. Original
    /// transport sources are evaluated only at the actual positive endpoint.
    ///
    /// Zero total receipt is not a gas state: its enclosing physical owner must
    /// retain this evacuated volume. This method does not mutate any input.
    #[allow(clippy::too_many_arguments)]
    pub fn fill_candidate(
        &self,
        catalog: &ThermochemicalGasModel,
        inventory_kg: &BTreeMap<SpeciesId, f64>,
        reference_temperature_k: f64,
        reference_relative_energy_j: f64,
        controls: GasEnergySolve,
        mass_tolerance_kg: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EvacuatedGasFillCandidate, GasEnergyError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if inventory_kg.is_empty()
            || inventory_kg.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1
            || !positive(reference_temperature_k)
            || !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
                .contains(&reference_temperature_k)
            || !reference_relative_energy_j.is_finite()
            || !positive(mass_tolerance_kg)
            || !positive(energy_tolerance_j)
        {
            return Err(GasEnergyError::Invalid(
                "bounded received species, a supported energy reference and positive physical budgets required",
            ));
        }
        let sources = catalog.thermodynamics().components();
        let mut received_mass = Dd::ZERO;
        for (species, &mass) in inventory_kg {
            poll(&mut cancelled)?;
            if !mass.is_finite()
                || mass < 0.0
                || sources
                    .binary_search_by(|s| s.species().cmp(species))
                    .is_err()
            {
                return Err(GasEnergyError::Invalid(
                    "received inventories must be finite nonnegative masses on original catalog species",
                ));
            }
            add(&mut received_mass, Dd::from_f64(mass))?;
        }
        let mass = received_mass.to_f64();
        let density = mass / self.volume_m3;
        let change = reference_relative_energy_j / mass;
        if !positive(mass)
            || !positive(density)
            || !change.is_finite()
            || (reference_relative_energy_j != 0.0 && change == 0.0)
        {
            return Err(GasEnergyError::Invalid(
                "first fill needs positive represented mass/density and a finite energy coordinate",
            ));
        }
        let model = model_for_inventory(catalog, inventory_kg)?;
        poll(&mut cancelled)?;
        // This is a mathematical chart for the newly received composition,
        // not a sampled vacuum state. Its pressure uses the original ideal-gas
        // constant, and final transport admission uses only the solved endpoint.
        let reference = model.thermodynamics().components()[0]
            .convention()
            .reference_pressure();
        let properties = model
            .thermodynamics()
            .evaluate_properties(Temperature::new(reference_temperature_k), reference)
            .map_err(|e| failure("evacuated fill energy reference", e))?;
        let pressure =
            density * properties.mass_specific_gas_constant().value() * reference_temperature_k;
        if !positive(pressure) {
            return Err(GasEnergyError::Invalid(
                "first-fill reference pressure is unrepresentable",
            ));
        }
        let chart = updated_point(&self.context, reference_temperature_k, pressure)?;
        // Tighten the original specific budget if the extensive J budget is
        // stronger. No step count or formation datum scales either allowance.
        let controls = GasEnergySolve {
            energy_tolerance_j_kg: controls
                .energy_tolerance_j_kg
                .min(energy_tolerance_j / mass),
            ..controls
        };
        let solved = model.reconstruct_density_internal_energy_increment(
            &chart,
            density,
            change,
            controls,
            &mut cancelled,
        )?;
        let state = solved.state().clone();
        let point = updated_point(
            &self.context,
            state.state().temperature,
            state.state().pressure,
        )?;
        let gas = ClosedGasParcel::from_isochoric(IsochoricGasParcel {
            model,
            point,
            state,
            mass_kg: mass,
            volume_m3: self.volume_m3,
        });
        let represented = inventory(gas.model(), mass)?;
        let mut represented_mass = Dd::ZERO;
        let mut stored_energy = Dd::ZERO;
        let mut maximum_species_mass_residual_kg = 0.0_f64;
        for (species, &expected) in inventory_kg {
            poll(&mut cancelled)?;
            let actual = represented.get(species).copied().unwrap_or(0.0);
            maximum_species_mass_residual_kg =
                maximum_species_mass_residual_kg.max((actual - expected).abs());
        }
        let sources = gas.model().thermodynamics().components();
        for (species, &actual) in &represented {
            poll(&mut cancelled)?;
            add(&mut represented_mass, Dd::from_f64(actual))?;
            let source = &sources[sources
                .binary_search_by(|s| s.species().cmp(species))
                .map_err(|_| {
                    GasEnergyError::Invalid("first fill lost an original species source")
                })?];
            let change = source
                .molar_internal_energy_change(
                    Temperature::new(reference_temperature_k),
                    Temperature::new(gas.state().state().temperature),
                )
                .map_err(|e| failure("evacuated fill storage increment", e))?
                .value();
            add(
                &mut stored_energy,
                Dd::from_f64(actual) * Dd::from_f64(change)
                    / Dd::from_f64(source.molar_mass().value()),
            )?;
        }
        let mass_residual_kg = (represented_mass - received_mass).to_f64();
        let density_mass_residual_kg = (Dd::from_f64(gas.state().state().density)
            * Dd::from_f64(self.volume_m3)
            - received_mass)
            .to_f64();
        let stored_relative_energy_j = stored_energy.to_f64();
        let energy_residual_j =
            (stored_energy - Dd::from_f64(reference_relative_energy_j)).to_f64();
        if [
            mass_residual_kg,
            density_mass_residual_kg,
            maximum_species_mass_residual_kg,
        ]
        .into_iter()
        .any(|e| !e.is_finite() || e.abs() > mass_tolerance_kg)
            || !stored_relative_energy_j.is_finite()
            || !energy_residual_j.is_finite()
            || energy_residual_j.abs() > energy_tolerance_j
        {
            return Err(GasEnergyError::Invalid(
                "first fill exceeded its fixed species/mass/energy budget",
            ));
        }
        let step = EvacuatedGasFillStep {
            received_mass_kg: mass,
            reference_temperature_k,
            reference_relative_energy_j,
            stored_relative_energy_j,
            mass_residual_kg,
            density_mass_residual_kg,
            maximum_species_mass_residual_kg,
            energy_residual_j,
            iterations: solved.iterations(),
        };
        poll(&mut cancelled)?;
        Ok(EvacuatedGasFillCandidate { gas, step })
    }
}

#[cfg(test)]
mod tests;
