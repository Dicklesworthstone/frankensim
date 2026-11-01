//! Original wall loss and implicit storage on unpublished gas candidates.

use super::{
    Conserved, Error, Euler1d, FrozenGas, GasNetwork1d, Totals, allocated, poll, positive,
};
use crate::dynamics::{CaloricDomain, Mesh1d};
use crate::multicomponent::network::wall_exchange::NetworkWallReport;
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallCellReport, WallContact, WallControl, WallDrag, WallExchangeError,
    WallGasCaloric, WallGasState, WallThermal, exchange_caloric,
    exchange_caloric_with_heat_storage,
};
use fs_qty::{Pressure, Temperature, Time};

mod branch;
use branch::BranchBinding;
pub use branch::BranchWallContact;

mod shared;
use shared::SharedWalls;
pub use shared::{SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallReport};

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

/// Prepared branch/plenum wall exchange using the original common wall kernel.
///
/// Gas amount and volume are fixed during each source. Branch drag converts
/// kinetic to sensible energy at fixed total energy; it does no extra heat or
/// gravitational work. Plenums have no bulk momentum and admit no drag. Original
/// wall h/u and latent inventory retain their declared physical convention.
/// The source is applied after transport, over its ACTUAL accepted duration:
/// this is first-order thermal/transport splitting, not a second-order IMEX
/// method. Stiff heat exchange is backward Euler with bounded source inversion.
pub struct NetworkWallOperator<'law> {
    gas: FrozenGas,
    plenum_count: usize,
    branches: Vec<Option<BranchBinding>>,
    branch_source_candidate: Vec<Totals>,
    volumes: Vec<f64>,
    contacts: Vec<WallContact>,
    candidate_contacts: Vec<WallContact>,
    heat_storage: Vec<Option<MaterialWallCaloric<'law>>>,
    shared: Option<SharedWalls<'law>>,
    last_report: Option<NetworkWallReport>,
}

/// Original plenum API name for the same prepared wall owner. Branch contacts
/// are absent unless explicitly supplied through `with_branch_contacts`.
pub type PlenumThermalOperator<'law> = NetworkWallOperator<'law>;

impl<'law> NetworkWallOperator<'law> {
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
            plenum_count: count,
            branches: allocated(network.branches.len(), None)?,
            branch_source_candidate: allocated(network.branches.len(), Totals::default())?,
            volumes: allocated(count, 0.0)?,
            contacts: allocated(count, empty)?,
            candidate_contacts: allocated(count, empty)?,
            heat_storage: allocated(count, None)?,
            shared: None,
            last_report: None,
        };
        for (i, input) in contacts.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let contact = WallContact {
                drag: WallDrag::None,
                thermal: input.thermal,
            };
            let contact = admit_contact(contact, input.heat_storage.as_ref())?;
            result.volumes[i] = network.plenums[i].volume_m3;
            result.contacts[i] = contact;
            result.heat_storage[i] = input.heat_storage;
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Accepted contact temperatures or original h/u inventories, in plenum order.
    pub fn contacts(&self) -> &[WallContact] {
        &self.contacts[..self.plenum_count]
    }

    /// Retained original finite-storage owners in that same order.
    pub fn heat_storage(&self) -> &[Option<MaterialWallCaloric<'law>>] {
        &self.heat_storage[..self.plenum_count]
    }

    /// The accepted source interval. Failed complete network steps retain it.
    pub const fn last_report(&self) -> Option<NetworkWallReport> {
        self.last_report
    }

    pub(super) fn validate_binding(&self, network: &GasNetwork1d) -> Result<(), Error> {
        if self.gas != *network.gas()
            || self.plenum_count != network.plenums.len()
            || self.branches.len() != network.branches.len()
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
        if self
            .branches
            .iter()
            .zip(&network.branches)
            .any(|(binding, branch)| {
                binding
                    .as_ref()
                    .is_some_and(|binding| binding.mesh != branch.mesh)
            })
        {
            return Err(Error::Invalid("branch wall geometry binding differs"));
        }
        Ok(())
    }

    pub(super) fn branch_source_candidates(&self) -> &[Totals] {
        &self.branch_source_candidate
    }

    pub(super) fn prepare(
        &mut self,
        branches: &mut [Euler1d],
        candidates: &mut [Conserved],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkWallReport, Error> {
        let mut total = HeatSum::default();
        self.branch_source_candidate.fill(Totals::default());
        for (i, state) in candidates.iter_mut().enumerate() {
            poll(cancelled)?;
            let report = self.prepare_cell(i, state, duration, controls, cancelled)?;
            if state.momentum != 0.0 {
                return Err(Error::Invalid(
                    "plenum heat produced an inadmissible conservative state",
                ));
            }
            total.add(report)?;
        }
        for (b, branch) in branches.iter_mut().enumerate() {
            let Some(binding) = &self.branches[b] else {
                continue;
            };
            let start = binding.start;
            let mut source = [0.0; 2];
            let mut correction = [0.0; 2];
            for (i, state) in branch.work.candidate.iter_mut().enumerate() {
                poll(cancelled)?;
                let report = self.prepare_cell(start + i, state, duration, controls, cancelled)?;
                // One branch axis only. Different support impulses are never
                // summed into a fictitious scalar network momentum.
                accumulate(
                    &mut source,
                    &mut correction,
                    [-report.wall_impulse_ns, report.heat_to_gas_j],
                )?;
                total.add(report)?;
            }
            self.branch_source_candidate[b] = Totals {
                mass: 0.0,
                momentum: source[0],
                energy: source[1],
            };
        }
        if let Some(shared) = &mut self.shared {
            total.add(shared.prepare(
                &self.gas,
                branches,
                candidates,
                &mut self.branch_source_candidate,
                duration,
                controls,
                cancelled,
            )?)?;
        }
        poll(cancelled)?;
        Ok(total.report(duration))
    }

    fn prepare_cell(
        &mut self,
        index: usize,
        state: &mut Conserved,
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, Error> {
        let contact = self.contacts[index];
        let source = self.heat_storage[index].as_ref();
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
                self.volumes[index],
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
                self.volumes[index],
                duration,
                contact,
                controls,
                &mut *cancelled,
            )?
        };
        if !update.momentum.is_finite() || !update.energy.is_finite() {
            return Err(Error::Invalid("wall source produced a non-finite state"));
        }
        // Only unpublished momentum/energy change. Fixed amount and volume
        // mean no potential change; drag heating is internal conversion of E.
        state.momentum = update.momentum;
        state.energy = update.energy;
        self.candidate_contacts[index] = update.contact;
        Ok(update.report)
    }

    pub(super) fn commit(&mut self, report: NetworkWallReport) {
        core::mem::swap(&mut self.contacts, &mut self.candidate_contacts);
        if let Some(shared) = &mut self.shared {
            shared.commit();
        }
        self.last_report = Some(report);
    }
}

fn admit_contact(
    mut contact: WallContact,
    source: Option<&MaterialWallCaloric<'_>>,
) -> Result<WallContact, Error> {
    contact.validate()?;
    if let Some(source) = source {
        contact = source.bind_contact(contact)?;
        source.temperature_k(contact.thermal)?;
    } else if matches!(contact.thermal, WallThermal::Caloric { .. }) {
        return Err(Error::Invalid("caloric wall requires its original source"));
    }
    Ok(contact)
}

#[derive(Default)]
struct HeatSum {
    values: [f64; 5],
    correction: [f64; 5],
    maximum_thermal_residual_j: f64,
    maximum_iterations: usize,
}

impl HeatSum {
    fn add(&mut self, report: WallCellReport) -> Result<(), Error> {
        accumulate(
            &mut self.values,
            &mut self.correction,
            [
                report.heat_to_gas_j,
                report.wall_energy_change_j,
                report.external_heat_to_gas_j,
                report.energy_residual_j,
                report.drag_heating_j,
            ],
        )?;
        self.maximum_thermal_residual_j = self
            .maximum_thermal_residual_j
            .max(report.thermal_residual_j.abs());
        self.maximum_iterations = self.maximum_iterations.max(report.iterations);
        Ok(())
    }

    fn report(self, duration: Time) -> NetworkWallReport {
        NetworkWallReport {
            elapsed: duration,
            heat_to_gas_j: self.values[0],
            wall_energy_change_j: self.values[1],
            external_heat_to_gas_j: self.values[2],
            energy_residual_j: self.values[3],
            drag_heating_j: self.values[4],
            maximum_thermal_residual_j: self.maximum_thermal_residual_j,
            maximum_iterations: self.maximum_iterations,
        }
    }
}

fn accumulate<const N: usize>(
    sum: &mut [f64; N],
    correction: &mut [f64; N],
    values: [f64; N],
) -> Result<(), Error> {
    for (k, value) in values.into_iter().enumerate() {
        let increment = value - correction[k];
        let next = sum[k] + increment;
        if !next.is_finite() {
            return Err(Error::Invalid("unrepresentable combined wall source"));
        }
        correction[k] = (next - sum[k]) - increment;
        sum[k] = next;
    }
    Ok(())
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
mod branch_tests;
#[cfg(test)]
mod shared_tests;
#[cfg(test)]
mod tests;
