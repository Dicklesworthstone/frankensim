//! Original implicit wall storage on the unpublished compact-volume candidates.

use super::{Conserved, Error, FrozenGas, GasNetwork1d, allocated, poll, positive};
use crate::dynamics::CaloricDomain;
use crate::multicomponent::network::wall_exchange::NetworkWallReport;
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallContact, WallControl, WallDrag, WallExchangeError, WallGasCaloric,
    WallGasState, WallThermal, exchange_caloric, exchange_caloric_with_heat_storage,
};
use fs_qty::{Pressure, Temperature, Time};

/// One compact volume's thermal boundary and optional original finite storage.
/// Distinct contacts represent distinct physical wall masses. A shared wall
/// requires its own thermal owner, rather than duplicating its stored energy.
#[derive(Clone, Debug, PartialEq)]
pub struct PlenumThermalContact<'law> {
    /// Prescribed reservoir, authored constant capacity, or original h/u state.
    pub thermal: WallThermal,
    /// Original material calorics, required for `WallThermal::Caloric` and
    /// optional for a temperature-coordinate finite wall. None retains the
    /// explicitly authored constant capacity; no Cp-to-Cv conversion is inferred.
    pub heat_storage: Option<MaterialWallCaloric<'law>>,
}

/// Whole-network balances after its accepted wave interval and thermal source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkThermalReport {
    /// Original implicit wall equation and represented gas/wall transfers.
    pub source: NetworkWallReport,
    /// Gas plus finite-wall change minus external gas/heat transport and
    /// resolved gravitational work, in the original formation-energy gauge [J].
    pub gas_wall_conservation_residual_j: f64,
    /// Gas plus finite-wall and centered-potential balance [J], when gravity
    /// is declared. External potential transport is counted exactly once.
    pub gas_wall_potential_conservation_residual_j: Option<f64>,
    /// The preceding residual in the declared absolute potential coordinate:
    /// add `reference_potential * mass_residual`, retaining the gas energy gauge.
    pub absolute_gas_wall_potential_conservation_residual_j: Option<f64>,
}

/// Prepared finite-plenum heat exchange using the original common wall kernel.
///
/// Gas amount and volume are fixed during the source. Plenums have no bulk
/// momentum, so this operator has no drag or moving-boundary work. Original
/// wall h/u and latent inventory retain their declared physical convention.
/// The source is applied after transport, over its ACTUAL accepted duration:
/// this is first-order thermal/transport splitting, not a second-order IMEX
/// method. Stiff heat exchange is backward Euler with bounded source inversion.
pub struct PlenumThermalOperator<'law> {
    gas: FrozenGas,
    volumes: Vec<f64>,
    contacts: Vec<WallContact>,
    candidate_contacts: Vec<WallContact>,
    heat_storage: Vec<Option<MaterialWallCaloric<'law>>>,
    last_report: Option<NetworkWallReport>,
}

impl<'law> PlenumThermalOperator<'law> {
    /// Bind exactly one original contact to each retained compact volume.
    pub fn new(
        network: &GasNetwork1d,
        contacts: Vec<PlenumThermalContact<'law>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if contacts.len() != network.plenums.len() {
            return Err(Error::Invalid("plenum thermal contact count"));
        }
        let count = contacts.len();
        let empty = WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        };
        let mut result = Self {
            gas: network.gas().clone(),
            volumes: allocated(count, 0.0)?,
            contacts: allocated(count, empty)?,
            candidate_contacts: allocated(count, empty)?,
            heat_storage: allocated(count, None)?,
            last_report: None,
        };
        for (i, input) in contacts.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let mut contact = WallContact {
                drag: WallDrag::None,
                thermal: input.thermal,
            };
            contact.validate()?;
            if let Some(source) = &input.heat_storage {
                contact = source.bind_contact(contact)?;
                source.temperature_k(contact.thermal)?;
            } else if matches!(contact.thermal, WallThermal::Caloric { .. }) {
                return Err(Error::Invalid(
                    "plenum caloric wall requires its original source",
                ));
            }
            result.volumes[i] = network.plenums[i].volume_m3;
            result.contacts[i] = contact;
            result.heat_storage[i] = input.heat_storage;
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Accepted contact temperatures or original h/u inventories, in plenum order.
    pub fn contacts(&self) -> &[WallContact] {
        &self.contacts
    }

    /// Retained original finite-storage owners in that same order.
    pub fn heat_storage(&self) -> &[Option<MaterialWallCaloric<'law>>] {
        &self.heat_storage
    }

    /// The accepted source interval. Failed complete network steps retain it.
    pub const fn last_report(&self) -> Option<NetworkWallReport> {
        self.last_report
    }

    pub(super) fn validate_binding(&self, network: &GasNetwork1d) -> Result<(), Error> {
        if self.gas != *network.gas()
            || self.volumes.len() != network.plenums.len()
            || self
                .volumes
                .iter()
                .zip(&network.plenums)
                .any(|(volume, plenum)| volume.to_bits() != plenum.volume_m3.to_bits())
        {
            return Err(Error::Invalid(
                "plenum thermal source or volume binding differs",
            ));
        }
        Ok(())
    }

    pub(super) fn prepare(
        &mut self,
        candidates: &mut [Conserved],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkWallReport, Error> {
        let mut sums = [0.0; 4];
        let mut compensation = [0.0; 4];
        let mut maximum_thermal_residual_j = 0.0_f64;
        let mut maximum_iterations = 0;
        for (i, state) in candidates.iter_mut().enumerate() {
            poll(cancelled)?;
            let contact = self.contacts[i];
            let source = self.heat_storage[i].as_ref();
            let hint = if let Some(source) = source {
                Some(source.temperature_k(contact.thermal)?)
            } else {
                match contact.thermal {
                    WallThermal::Reservoir { temperature_k, .. }
                    | WallThermal::Lumped { temperature_k, .. } => Some(temperature_k),
                    _ => None,
                }
            };
            let gas = FrozenCaloric {
                gas: &self.gas,
                density: state.density,
                hint,
            };
            let update = if let Some(source) = source {
                exchange_caloric_with_heat_storage(
                    &gas,
                    state.momentum,
                    state.energy,
                    self.volumes[i],
                    duration,
                    contact,
                    source,
                    None,
                    controls,
                    &mut *cancelled,
                )?
            } else {
                exchange_caloric(
                    &gas,
                    state.momentum,
                    state.energy,
                    self.volumes[i],
                    duration,
                    contact,
                    controls,
                    &mut *cancelled,
                )?
            };
            if update.momentum != 0.0 || !update.energy.is_finite() {
                return Err(Error::Invalid(
                    "plenum heat produced an inadmissible conservative state",
                ));
            }
            // Only unpublished energy changes: amount, potential and the wave
            // fluxes retain their original values. Heat causes no gravity work.
            state.energy = update.energy;
            self.candidate_contacts[i] = update.contact;
            let report = update.report;
            for (k, value) in [
                report.heat_to_gas_j,
                report.wall_energy_change_j,
                report.external_heat_to_gas_j,
                report.energy_residual_j,
            ]
            .into_iter()
            .enumerate()
            {
                let increment = value - compensation[k];
                let next = sums[k] + increment;
                if !next.is_finite() {
                    return Err(Error::Invalid(
                        "unrepresentable combined plenum heat transfer",
                    ));
                }
                compensation[k] = (next - sums[k]) - increment;
                sums[k] = next;
            }
            maximum_thermal_residual_j =
                maximum_thermal_residual_j.max(report.thermal_residual_j.abs());
            maximum_iterations = maximum_iterations.max(report.iterations);
        }
        poll(cancelled)?;
        Ok(NetworkWallReport {
            elapsed: duration,
            heat_to_gas_j: sums[0],
            wall_energy_change_j: sums[1],
            external_heat_to_gas_j: sums[2],
            energy_residual_j: sums[3],
            drag_heating_j: 0.0,
            maximum_thermal_residual_j,
            maximum_iterations,
        })
    }

    pub(super) fn commit(&mut self, report: NetworkWallReport) {
        core::mem::swap(&mut self.contacts, &mut self.candidate_contacts);
        self.last_report = Some(report);
    }
}

/// Borrow the exact frozen source; do not manufacture an equivalent species
/// bank or convert its formation-energy zero through absolute energy levels.
struct FrozenCaloric<'a> {
    gas: &'a FrozenGas,
    density: f64,
    hint: Option<f64>,
}

impl WallGasCaloric for FrozenCaloric<'_> {
    fn domain(&self) -> CaloricDomain {
        self.gas.domain()
    }

    fn primitive(
        &self,
        momentum: f64,
        energy: f64,
        mut cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallGasState, WallExchangeError> {
        poll(&mut cancelled)?;
        // A represented state exactly reproducible at the original wall
        // temperature is an admitted thermal fixed point. Avoid introducing
        // inversion tolerance into an otherwise exactly zero heat flux.
        if let Some(temperature) = self.hint.filter(|t| {
            (self.domain().minimum_temperature..=self.domain().maximum_temperature).contains(t)
        }) {
            let (thermal, _) = self.caloric(temperature)?;
            let kinetic = 0.5 * momentum * (momentum / self.density);
            if energy == thermal + kinetic {
                return Ok(WallGasState {
                    density: self.density,
                    velocity: momentum / self.density,
                    temperature,
                });
            }
        }
        let state = self.gas.decode(
            Conserved {
                density: self.density,
                momentum,
                energy,
            },
            None,
            &mut cancelled,
        )?;
        Ok(WallGasState {
            density: state.density,
            velocity: state.velocity,
            temperature: state.temperature,
        })
    }

    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        if !(self.domain().minimum_temperature..=self.domain().maximum_temperature)
            .contains(&temperature)
        {
            return Err(
                Error::CaloricPath("plenum thermal source outside declared gas interval").into(),
            );
        }
        let pressure = self.density * (self.gas.gas_constant() * temperature);
        let properties = self
            .gas
            .source()
            .evaluate_properties(Temperature::new(temperature), Pressure::new(pressure))
            .map_err(Error::from)?;
        let energy = self.density * self.gas.sensible(temperature)?;
        let cv = self.density * properties.mass_specific_cv().value();
        if !energy.is_finite() || !positive(cv) {
            return Err(Error::CaloricPath("unrepresentable original plenum calorics").into());
        }
        Ok((energy, cv))
    }
}

#[cfg(test)]
mod tests;
