//! Simultaneous gas/film residuals for the original fully-wetting wall.
//!
//! A thermal owner supplies trial wall enthalpy and nonadvective gas heat. The
//! original Stefan law sees the trial FINAL gas, including its composition,
//! pressure and temperature. This module owns no nonlinear solver or wall heat
//! capacity. It exposes the two physical gas residuals and the opposite wall
//! heat port, then reconstructs the original gas once after convergence.
//!
//! The film complementarity is `m1 = max(0, m0 - dt*A*j_endpoint)`: no seed or
//! capacity survives dry-out. The reported exhaustion time freezes the endpoint
//! rate over the interval; physical event times still require time refinement.

use super::*;
use crate::evaporation::gas::chamber::Total;
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::admit_catalog;
use fs_qty::Temperature;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
use std::collections::BTreeMap;

/// Prepared original inventories and sources for one positive physical step.
pub struct GasWallFilmEndpoint<'a, 'catalog> {
    surface: GasWallFilmSurface<'a>,
    gas: &'a ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasEvaporationSolve,
    inventory: BTreeMap<SpeciesId, f64>,
    old_vapor_mass_kg: f64,
    gas_constant_j_kg_k: f64,
}

/// Unaccepted physical residual. A small numerical norm alone cannot publish it.
pub struct GasWallFilmEndpointTrial {
    gas: ClosedGasParcel,
    prepared: PreparedWallFilm,
    interface_liquid_enthalpy_j_kg: f64,
    heat_into_gas_j: f64,
    gas_mass_residual_kg: f64,
    maximum_mass_residual_kg: f64,
    gas_energy_residual_j: f64,
}
impl GasWallFilmEndpointTrial {
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn film(&self) -> GasWallFilmInventory {
        self.prepared.film
    }
    pub const fn step(&self) -> &GasWallFilmStep {
        &self.prepared.step
    }
    /// Actual volatile inventory change minus the endpoint phase transfer [kg].
    pub const fn gas_mass_residual_kg(&self) -> f64 {
        self.gas_mass_residual_kg
    }
    /// Maximum total, volatile, carrier and coupled film/gas mass defect [kg].
    pub const fn maximum_mass_residual_kg(&self) -> f64 {
        self.maximum_mass_residual_kg
    }
    /// Original source-relative gas storage minus vapor stream and supplied heat.
    pub const fn gas_energy_residual_j(&self) -> f64 {
        self.gas_energy_residual_j
    }
    /// Explicit nonadvective heat only; vapor-carried enthalpy is separate.
    pub const fn heat_into_gas_j(&self) -> f64 {
        self.heat_into_gas_j
    }
}

/// Joint gas/film candidate, checked again after original caloric inversion.
#[derive(Clone, Debug, PartialEq)]
pub struct GasWallFilmEndpointCandidate<'catalog> {
    candidate: GasWallFilmCandidate<'catalog>,
    heat_into_gas_j: f64,
}
impl<'catalog> GasWallFilmEndpointCandidate<'catalog> {
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.candidate.gas()
    }
    pub const fn film(&self) -> GasWallFilmInventory {
        self.candidate.film()
    }
    pub const fn step(&self) -> &GasWallFilmStep {
        self.candidate.step()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.candidate.catalog()
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        self.candidate.gas_exchange()
    }
    pub const fn heat_into_gas_j(&self) -> f64 {
        self.heat_into_gas_j
    }
    pub fn into_parts(
        self,
    ) -> (
        GasWallFilmInventory,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasWallFilmStep,
    ) {
        self.candidate.into_parts()
    }
}

impl<'a, 'catalog> GasWallFilmEndpoint<'a, 'catalog> {
    pub fn new(
        surface: GasWallFilmSurface<'a>,
        gas: &'a ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        duration_s: f64,
        controls: GasEvaporationSolve,
    ) -> Result<Self, GasEvaporationError> {
        validate_controls(controls, duration_s, surface.mass_transfer_coefficient_m_s)?;
        if !positive(duration_s) {
            return Err(invalid(
                "wall-film endpoint needs a positive physical interval",
            ));
        }
        admit_catalog(gas, catalog)?;
        // Zero exposure validates the original film and caloric binding without
        // spending donor inventory on a frozen INITIAL-rate predictor.
        prepare_wall_film(surface, gas, gas, catalog, 0.0, controls, || false)?;
        let inventory = gas.species_inventory_kg()?;
        let old_vapor_mass_kg = inventory
            .get(&surface.binding.species)
            .copied()
            .unwrap_or(0.0);
        let species = catalog
            .thermodynamics()
            .components()
            .iter()
            .find(|species| species.species() == &surface.binding.species)
            .ok_or(invalid(
                "wall-film endpoint needs its original volatile source",
            ))?;
        let gas_constant_j_kg_k =
            finite(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / species.molar_mass().value())?;
        let prepared = Self {
            surface,
            gas,
            catalog,
            duration_s,
            controls,
            inventory,
            old_vapor_mass_kg,
            gas_constant_j_kg_k,
        };
        prepared.admit_temperature(gas.state().state().temperature)?;
        Ok(prepared)
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), GasEvaporationError> {
        if !(self.controls.gas.gas.minimum_temperature_k
            ..=self.controls.gas.gas.maximum_temperature_k)
            .contains(&temperature_k)
        {
            return Err(invalid(
                "wall-film gas endpoint left its original caloric bracket",
            ));
        }
        if (temperature_k - self.gas.state().state().temperature).abs()
            > self.controls.maximum_temperature_increment_k
        {
            return Err(smaller_step(
                self.duration_s,
                "wall-film endpoint gas temperature departure",
            ));
        }
        Ok(())
    }

    /// Build a trial through the ORIGINAL gas EOS and source catalog, without a
    /// caloric root. The mass coordinate is signed gas inventory change; the
    /// residual compares it to the separately evaluated endpoint Stefan flux.
    pub fn evaluate(
        &self,
        interface_liquid_enthalpy_j_kg: f64,
        transferred_mass_coordinate_kg: f64,
        gas_temperature_k: f64,
        gas_heat_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasWallFilmEndpointTrial, GasEvaporationError> {
        poll(&mut cancelled)?;
        self.admit_temperature(gas_temperature_k)?;
        let mass = finite(self.old_vapor_mass_kg + transferred_mass_coordinate_kg)?;
        if !nonnegative(mass)
            || (transferred_mass_coordinate_kg != 0.0 && mass == self.old_vapor_mass_kg)
        {
            return Err(invalid(
                "wall-film gas coordinate cannot represent its vapor inventory",
            ));
        }
        let gas = if transferred_mass_coordinate_kg == 0.0 {
            self.gas.at_temperature(gas_temperature_k)?.0
        } else {
            let mut inventory = self.inventory.clone();
            inventory.insert(self.surface.binding.species.clone(), mass);
            self.gas.at_species_inventory_temperature(
                self.catalog,
                &inventory,
                gas_temperature_k,
            )?
        };
        self.evaluate_at_gas(interface_liquid_enthalpy_j_kg, &gas, gas_heat_j, cancelled)
    }

    /// Recheck the same residual against an actual reconstructed gas endpoint.
    /// Unchanged carrier inventories and volume remain independent gates.
    pub fn evaluate_at_gas(
        &self,
        interface_liquid_enthalpy_j_kg: f64,
        gas: &ClosedGasParcel,
        gas_heat_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasWallFilmEndpointTrial, GasEvaporationError> {
        poll(&mut cancelled)?;
        admit_catalog(gas, self.catalog)?;
        let t0 = self.gas.state().state().temperature;
        let tg = gas.state().state().temperature;
        self.admit_temperature(tg)?;
        if !gas_heat_j.is_finite() || gas.volume_m3() != self.gas.volume_m3() {
            return Err(invalid(
                "wall-film endpoint needs finite heat and unchanged gas volume",
            ));
        }
        if gas.point().axis_quantities() != self.gas.point().axis_quantities()
            || gas.point().axes().len() != self.gas.point().axes().len()
            || self.gas.point().axes().iter().any(|(axis, value)| {
                axis != "temperature"
                    && axis != "pressure"
                    && gas.point().axes().get(axis).map(|other| other.to_bits())
                        != Some(value.to_bits())
            })
        {
            return Err(invalid(
                "wall-film endpoint changes its retained gas query context",
            ));
        }
        let mut prepared = prepare_wall_film(
            GasWallFilmSurface {
                interface_liquid_enthalpy_j_kg,
                ..self.surface
            },
            self.gas,
            gas,
            self.catalog,
            self.duration_s,
            self.controls,
            &mut cancelled,
        )?;
        let dm = prepared.transfer.mass_into_gas_kg;
        let actual = gas.species_inventory_kg()?;
        let vapor_change = finite(
            actual
                .get(&self.surface.binding.species)
                .copied()
                .unwrap_or(0.0)
                - self.old_vapor_mass_kg,
        )?;
        let gas_mass_residual_kg = finite(vapor_change - dm)?;
        let total_change = finite(gas.mass_kg() - self.gas.mass_kg())?;
        let liquid_change = finite(prepared.film.mass_kg() - self.surface.film.mass_kg())?;
        let mut maximum_mass_residual_kg = gas_mass_residual_kg
            .abs()
            .max(finite(total_change - dm)?.abs())
            .max(prepared.step.liquid_mass_residual_kg.abs())
            .max(finite(liquid_change + total_change)?.abs())
            .max(finite(liquid_change + vapor_change)?.abs());
        for species in self.catalog.thermodynamics().components() {
            poll(&mut cancelled)?;
            let id = species.species();
            if id != &self.surface.binding.species {
                maximum_mass_residual_kg = maximum_mass_residual_kg.max(
                    finite(
                        actual.get(id).copied().unwrap_or(0.0)
                            - self.inventory.get(id).copied().unwrap_or(0.0),
                    )?
                    .abs(),
                );
            }
        }
        let increment = gas
            .model()
            .thermodynamics()
            .specific_internal_energy_change(Temperature::new(t0), Temperature::new(tg))
            .map_err(|error| {
                crate::thermochemical_gas::energy::GasEnergyError::Material(
                    crate::MaterialError::Parameters {
                        what: format!("wall-film original caloric increment: {error}"),
                    },
                )
            })?
            .value();
        let sensible = finite(gas.mass_kg() * increment)?;
        let pressure_flow = finite(dm * self.gas_constant_j_kg_k * t0)?;
        let relative = prepared.transfer.relative_energy_into_gas_j;
        let gas_energy_residual_j = sum([sensible, -pressure_flow, -relative, -gas_heat_j])?;
        let step = &mut prepared.step;
        step.gas_internal_energy_change_j =
            sum([step.energy_into_gas_j, -relative, -pressure_flow, sensible])?;
        step.mass_residual_kg = finite(liquid_change + vapor_change)?;
        step.energy_residual_j = sum([
            gas_energy_residual_j,
            prepared.film_energy_residual_j,
            prepared.wall_heat_roundoff_j,
        ])?;
        step.sum_absolute_energy_residual_j = sum([
            gas_energy_residual_j.abs(),
            prepared.film_energy_residual_j.abs(),
            prepared.wall_heat_roundoff_j.abs(),
        ])?;
        poll(&mut cancelled)?;
        Ok(GasWallFilmEndpointTrial {
            gas: gas.clone(),
            prepared,
            interface_liquid_enthalpy_j_kg,
            heat_into_gas_j: gas_heat_j,
            gas_mass_residual_kg,
            maximum_mass_residual_kg,
            gas_energy_residual_j,
        })
    }

    /// Gate physical kilograms and joules independently of numerical scaling.
    pub fn check_balance(
        &self,
        trial: &GasWallFilmEndpointTrial,
    ) -> Result<(), GasEvaporationError> {
        if trial.maximum_mass_residual_kg > self.controls.gas.mass_tolerance_kg {
            return Err(invalid(
                "wall-film endpoint fails its original species/total mass allowance",
            ));
        }
        if trial.step().sum_absolute_energy_residual_j > self.controls.gas.energy_tolerance_j {
            return Err(GasEvaporationError::EnergyBalance {
                residual_j: trial.step().sum_absolute_energy_residual_j,
                tolerance_j: self.controls.gas.energy_tolerance_j,
            });
        }
        Ok(())
    }

    /// Perform ONE original caloric reconstruction and ONE fresh phase residual
    /// evaluation. An enclosing numerical owner must charge that evaluation to
    /// its common budget. Failure/cancellation publishes neither gas nor film.
    pub fn finalize(
        &self,
        trial: &GasWallFilmEndpointTrial,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasWallFilmEndpointCandidate<'catalog>, GasEvaporationError> {
        poll(&mut cancelled)?;
        self.check_balance(trial)?;
        let gas = self.gas.exchange_species_boundary_relative_with_energy(
            self.catalog,
            std::slice::from_ref(&trial.prepared.transfer),
            |_| Ok(trial.heat_into_gas_j),
            self.controls.gas,
            &mut cancelled,
        )?;
        let fresh = self.evaluate_at_gas(
            trial.interface_liquid_enthalpy_j_kg,
            gas.gas(),
            trial.heat_into_gas_j,
            &mut cancelled,
        )?;
        self.check_balance(&fresh)?;
        let mut step = fresh.prepared.step;
        step.gas_internal_energy_change_j = gas.step().stored_energy_change_j;
        // The actual-endpoint residual already contains caloric reconstruction
        // error. Charge its original receipt independently in the ABSOLUTE audit
        // without counting that error twice in the signed physical balance.
        step.sum_absolute_energy_residual_j = sum([
            step.sum_absolute_energy_residual_j,
            gas.step().energy_residual_j.abs(),
        ])?;
        if step.sum_absolute_energy_residual_j > self.controls.gas.energy_tolerance_j {
            return Err(GasEvaporationError::EnergyBalance {
                residual_j: step.sum_absolute_energy_residual_j,
                tolerance_j: self.controls.gas.energy_tolerance_j,
            });
        }
        let mass_defect = sum([
            fresh.maximum_mass_residual_kg,
            gas.step().mass_residual_kg.abs(),
            gas.step().maximum_species_mass_residual_kg,
        ])?;
        if mass_defect > self.controls.gas.mass_tolerance_kg {
            return Err(invalid(
                "wall-film reconstruction exhausts its shared mass allowance",
            ));
        }
        poll(&mut cancelled)?;
        Ok(GasWallFilmEndpointCandidate {
            candidate: GasWallFilmCandidate {
                film: fresh.prepared.film,
                gas,
                step,
            },
            heat_into_gas_j: trial.heat_into_gas_j,
        })
    }
}

fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, GasEvaporationError> {
    let mut total = Total::default();
    for value in values {
        total.add(value)?;
    }
    total.value()
}

#[cfg(test)]
mod tests;
