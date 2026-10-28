//! Prescribed species and energy exchange with an external phase or boundary.
//!
//! The original gas owner reconstructs composition, temperature, pressure,
//! transport and acoustics once from the accepted extensive inventories.
//! The complete source catalog survives exact depletion of an individual
//! species; no artificial trace gas, second EOS or caloric law is introduced.
//!
//! Stream energies must use the original gas species' thermochemical datum.
//! A latent heat alone is not a vapor enthalpy. This operation supplies no
//! interface rate, solid energy law, phase equilibrium or entropy guarantee.

use super::{
    ClosedGasParcel, GasEnergyError, GasEnergySolve, MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1, SpeciesId,
    Sum, ThermochemicalGasModel, endpoint, invalid, inventory, model_for_inventory, poll, positive,
    source_error, stored_change,
};
use fs_qty::Temperature;
use std::collections::BTreeMap;

/// Integrated exchange, positive into the gas, on an explicit original species.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesBoundaryTransfer {
    pub species: SpeciesId,
    pub mass_into_gas_kg: f64,
    /// Complete signed stream energy [J], including the original enthalpy datum.
    /// For one advected stream this is signed mass times upstream partial h.
    /// Opposing streams of one species may instead supply their net mass and
    /// net energy; zero net mass does not imply zero transported energy.
    pub energy_into_gas_j: f64,
}

/// Original caloric reconstruction and absolute extensive balance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSpeciesBoundarySolve {
    pub gas: GasEnergySolve,
    pub mass_tolerance_kg: f64,
    pub energy_tolerance_j: f64,
}
impl GasSpeciesBoundarySolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.gas.validate()?;
        if !positive(self.mass_tolerance_kg) || !positive(self.energy_tolerance_j) {
            return Err(invalid(
                "positive finite boundary mass and energy tolerances required",
            ));
        }
        Ok(())
    }
}

/// Actual accepted gas inventories and energy, after original-owner inversion.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSpeciesBoundaryStep {
    /// Nonzero exchanges in canonical original SpeciesId order.
    pub species: Vec<SpeciesBoundaryTransfer>,
    pub mass_into_gas_kg: f64,
    pub energy_into_gas_j: f64,
    pub stored_energy_change_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub iterations: usize,
}

/// Proposed gas endpoint; the external phase must be accepted atomically with it.
/// Borrowing the complete catalog retains depleted sources without copying every
/// source card on each nonlinear trial. The accepted gas stores active species.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSpeciesBoundaryCandidate<'catalog> {
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    step: GasSpeciesBoundaryStep,
}
impl<'catalog> GasSpeciesBoundaryCandidate<'catalog> {
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    #[must_use]
    pub const fn step(&self) -> &GasSpeciesBoundaryStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasSpeciesBoundaryStep,
    ) {
        (self.gas, self.catalog, self.step)
    }
}

// Check source identity separately from composition: reconstructing mole
// fractions from mass may change their last rounding bit, but cannot authorize
// a different caloric law, transport claim, water binding or source receipt.
fn admit_catalog(
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
) -> Result<(), GasEnergyError> {
    for (species, transport) in gas
        .model()
        .thermodynamics()
        .components()
        .iter()
        .zip(gas.model().transport_materials())
    {
        let index = catalog
            .thermodynamics()
            .components()
            .iter()
            .position(|candidate| candidate.species() == species.species())
            .ok_or(invalid("boundary catalog is missing an active gas species"))?;
        if &catalog.thermodynamics().components()[index] != species
            || &catalog.transport_materials()[index] != transport
        {
            return Err(invalid(
                "boundary catalog substitutes an active caloric or transport source",
            ));
        }
        let original = gas
            .model()
            .source_packs()
            .iter()
            .find(|receipt| &receipt.species == species.species());
        let retained = catalog
            .source_packs()
            .iter()
            .find(|receipt| &receipt.species == species.species());
        if original != retained {
            return Err(invalid(
                "boundary catalog substitutes an active source-pack binding",
            ));
        }
    }
    let active_water = catalog.water_species().filter(|id| {
        gas.model()
            .thermodynamics()
            .components()
            .iter()
            .any(|species| species.species() == *id)
    });
    if gas.model().water_species() != active_water {
        return Err(invalid(
            "boundary catalog changes the explicit active water identity",
        ));
    }
    Ok(())
}

impl ThermochemicalGasModel {
    /// Original NASA-9 partial specific enthalpy, including its formation datum.
    ///
    /// This is an ideal-gas caloric query, not a complete material-state resolve:
    /// transport support is checked when reconstructing the accepted gas.
    /// Unsupported temperature or absent source identity refuses; no liquid,
    /// ice or phase-change enthalpy is inferred from the species name.
    pub fn species_specific_enthalpy_j_kg(
        &self,
        id: &SpeciesId,
        temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        if !positive(temperature_k) {
            return Err(invalid(
                "finite positive absolute species temperature required",
            ));
        }
        let species = self
            .thermodynamics()
            .components()
            .iter()
            .find(|species| species.species() == id)
            .ok_or(invalid(
                "species enthalpy needs its retained original source",
            ))?;
        let evaluated = species
            .evaluate(Temperature::new(temperature_k))
            .map_err(source_error)?;
        let enthalpy = evaluated.properties().h().value() / species.molar_mass().value();
        if !enthalpy.is_finite() {
            return Err(invalid("original species enthalpy is unrepresentable"));
        }
        Ok(enthalpy)
    }
}

impl ClosedGasParcel {
    /// Actual species inventory from the original mixture composition and mass.
    /// An absent active species has exactly zero mass, including after depletion.
    pub fn species_mass_kg(&self, id: &SpeciesId) -> Result<f64, GasEnergyError> {
        Ok(inventory(self.model(), self.mass_kg())?
            .get(id)
            .copied()
            .unwrap_or(0.0))
    }

    /// Partial pressure from the original resolved pressure and mole fraction.
    /// An absent active species has exactly zero partial pressure.
    pub fn species_partial_pressure_pa(&self, id: &SpeciesId) -> Result<f64, GasEnergyError> {
        let Some(index) = self
            .model()
            .thermodynamics()
            .components()
            .iter()
            .position(|species| species.species() == id)
        else {
            return Ok(0.0);
        };
        let fraction = self.model().thermodynamics().mole_composition().fractions()[index];
        let pressure = fraction * self.state().state().pressure;
        if !positive(pressure) {
            return Err(invalid(
                "active species partial pressure is unrepresentable",
            ));
        }
        Ok(pressure)
    }

    /// Reconstruct one fixed-volume gas from integrated external species streams.
    ///
    /// Retain the complete catalog between steps so a depleted species may
    /// reenter with its original sources. Other species masses, volume and
    /// non-temperature/pressure query coordinates are unchanged. Exact removal
    /// of one species is supported while total vacuum is refused.
    ///
    /// Each stream supplies its COMPLETE signed energy in the original gas
    /// datum. The external phase owner must debit that same energy and account
    /// for latent heat and flow work once. This API does not align references
    /// automatically. Net energy at zero net species mass is allowed for
    /// simultaneous counterflow. Duplicate identities refuse.
    ///
    /// No accepted input is mutated on success, cancellation or any failure.
    pub fn exchange_species_boundary<'catalog>(
        &self,
        catalog: &'catalog ThermochemicalGasModel,
        transfers: &[SpeciesBoundaryTransfer],
        controls: GasSpeciesBoundarySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasSpeciesBoundaryCandidate<'catalog>, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if transfers.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
            return Err(invalid(
                "boundary exchange exceeds the shared species bound",
            ));
        }
        admit_catalog(self, catalog)?;
        let mut ordered = BTreeMap::new();
        for row in transfers {
            poll(&mut cancelled)?;
            if !row.mass_into_gas_kg.is_finite()
                || !row.energy_into_gas_j.is_finite()
                || ordered.insert(row.species.clone(), row).is_some()
                || !catalog
                    .thermodynamics()
                    .components()
                    .iter()
                    .any(|species| species.species() == &row.species)
            {
                return Err(invalid(
                    "boundary streams need finite values and unique retained species",
                ));
            }
        }
        if ordered
            .values()
            .all(|row| row.mass_into_gas_kg == 0.0 && row.energy_into_gas_j == 0.0)
        {
            poll(&mut cancelled)?;
            return Ok(GasSpeciesBoundaryCandidate {
                gas: self.clone(),
                catalog,
                step: GasSpeciesBoundaryStep {
                    species: Vec::new(),
                    mass_into_gas_kg: 0.0,
                    energy_into_gas_j: 0.0,
                    stored_energy_change_j: 0.0,
                    mass_residual_kg: 0.0,
                    maximum_species_mass_residual_kg: 0.0,
                    energy_residual_j: 0.0,
                    iterations: 0,
                },
            });
        }
        let before = inventory(self.model(), self.mass_kg())?;
        let mut target = before.clone();
        let mut mass = Sum::default();
        let mut energy = Sum::default();
        let mut relative = Sum::default();
        let mut rows = Vec::with_capacity(ordered.len());
        for (id, row) in ordered {
            poll(&mut cancelled)?;
            if row.mass_into_gas_kg == 0.0 && row.energy_into_gas_j == 0.0 {
                continue;
            }
            let old = before.get(&id).copied().unwrap_or(0.0);
            let next = old + row.mass_into_gas_kg;
            if !next.is_finite() || next < 0.0 || (row.mass_into_gas_kg != 0.0 && next == old) {
                return Err(invalid(
                    "boundary stream exceeds or cannot change its species inventory",
                ));
            }
            target.insert(id, next);
            mass.add(row.mass_into_gas_kg)?;
            energy.add(row.energy_into_gas_j)?;
            relative.add(
                row.energy_into_gas_j
                    - row.mass_into_gas_kg * self.state().specific_internal_energy_j_kg(),
            )?;
            rows.push(row.clone());
        }
        let next_mass = self.mass_kg() + mass.value;
        if !positive(next_mass) || (mass.value != 0.0 && next_mass == self.mass_kg()) {
            return Err(invalid(
                "boundary exchange cannot produce vacuum or a vanished mass change",
            ));
        }
        let (gas, iterations) = endpoint(
            self,
            model_for_inventory(catalog, &target)?,
            next_mass,
            relative.value / next_mass,
            controls.gas,
            &mut cancelled,
        )?;
        let stored = stored_change(self, &gas);
        let mass_residual = (gas.mass_kg() - self.mass_kg()) - mass.value;
        let energy_residual = stored - energy.value;
        if ![stored, mass_residual, energy_residual]
            .iter()
            .all(|value| value.is_finite())
            || mass_residual.abs() > controls.mass_tolerance_kg
            || energy_residual.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "boundary gas endpoint failed total mass or energy balance",
            ));
        }
        let actual = inventory(gas.model(), gas.mass_kg())?;
        let mut maximum: f64 = 0.0;
        for species in catalog.thermodynamics().components() {
            poll(&mut cancelled)?;
            let id = species.species();
            let difference =
                actual.get(id).copied().unwrap_or(0.0) - target.get(id).copied().unwrap_or(0.0);
            if !difference.is_finite() {
                return Err(invalid("boundary species balance is unrepresentable"));
            }
            maximum = maximum.max(difference.abs());
        }
        if maximum > controls.mass_tolerance_kg {
            return Err(invalid(
                "boundary gas endpoint did not retain individual species inventories",
            ));
        }
        poll(&mut cancelled)?;
        Ok(GasSpeciesBoundaryCandidate {
            gas,
            catalog,
            step: GasSpeciesBoundaryStep {
                species: rows,
                mass_into_gas_kg: mass.value,
                energy_into_gas_j: energy.value,
                stored_energy_change_j: stored,
                mass_residual_kg: mass_residual,
                maximum_species_mass_residual_kg: maximum,
                energy_residual_j: energy_residual,
                iterations,
            },
        })
    }
}

#[cfg(test)]
mod tests;
