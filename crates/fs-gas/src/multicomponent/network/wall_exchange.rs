//! Stationary radial wall sources for every branch cell and finite plenum.
//!
//! Original mixture calorics remain authoritative. Each contact is explicitly
//! authored: linear or Darcy drag and adiabatic, prescribed-temperature, or
//! finite-capacity thermal storage. The coefficients are estimated inputs, not
//! a resolved unsteady shock boundary layer or a new heat-transfer correlation.
//!
//! Work arrays and wall states are prepared once. A complete source interval
//! changes no species or flow clock and publishes all gas cells and walls
//! together. Composing owners retain candidate wall temperatures through both
//! halves and discard them when the complete wave/reaction step is retried.
//! An optional original-source transport binding refreshes developed circular
//! branch drag and conductance from each cell's current mixture at every half.

use super::MixtureNetwork1d;
use crate::multicomponent::wall_exchange::{
    WallCellReport, WallContact, WallControl, WallDrag, WallExchangeError, WallThermal,
    exchange_cell,
};
use crate::multicomponent::{Cell, Mesh1d, MixtureGas, allocated, sum};
use fs_qty::Time;

mod transport;
pub use transport::CircularLaminarWall;
use transport::CircularWallTransport;

#[cfg(feature = "material")]
mod caloric;
#[cfg(feature = "material")]
use crate::multicomponent::wall_exchange::MaterialWallCaloric;

/// Source-only exchanges over a complete interval, or both composed wall halves.
///
/// Heat is positive into gas. Stationary-wall drag retains total gas energy:
/// `drag_heating_j` is kinetic energy converted to sensible energy, not an
/// additional heat input. Gas plus finite-wall energy changes equal external
/// reservoir heat, up to the separately retained numerical residual.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkWallReport {
    /// Total duration covered by the wall operator [s].
    pub elapsed: Time,
    /// Energy transferred to gas by finite walls and reservoirs [J].
    pub heat_to_gas_j: f64,
    /// Actual represented energy change of finite walls [J].
    pub wall_energy_change_j: f64,
    /// Heat into gas from prescribed-temperature reservoirs only [J].
    pub external_heat_to_gas_j: f64,
    /// Gas kinetic energy thermalized by stationary-wall drag [J].
    pub drag_heating_j: f64,
    /// Sum of local gas-plus-wall energy changes minus reservoir heat [J].
    pub energy_residual_j: f64,
    /// Largest absolute original implicit thermal-equation residual [J].
    pub maximum_thermal_residual_j: f64,
    /// Largest scalar thermal solve iteration count in any cell.
    pub maximum_iterations: usize,
}

/// Accepted wall intervals before and after the complete transport composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitNetworkWallReport {
    /// Wall half preceding molecular transport, chemistry, or advection.
    pub first: NetworkWallReport,
    /// Wall half following molecular transport, chemistry, or advection.
    pub second: NetworkWallReport,
}

struct BranchShape {
    mesh: Mesh1d,
    start: usize,
}

/// Prepared source operator bound to original gas calorics and cell geometry.
///
/// Each branch cell and each plenum owns one distinct thermal contact. A
/// finite-capacity wall is a lumped storage attached to that cell; sharing one
/// physical wall between multiple cells requires an additional thermal owner.
/// Drag coefficients are frozen during each source interval. Plenums have no
/// bulk momentum and therefore admit only `WallDrag::None`.
/// [`Self::with_circular_laminar_transport`] derives selected branch coefficients
/// from current source-owned gas properties instead of retaining authored values.
/// With the `material` feature, `with_material_heat_storage` additionally binds
/// original nonlinear wall calorics without adding a second temperature history.
pub struct NetworkWallOperator {
    gas: MixtureGas,
    branches: Vec<BranchShape>,
    plenum_volumes: Vec<f64>,
    plenum_start: usize,
    contacts: Vec<WallContact>,
    candidate_contacts: Vec<WallContact>,
    candidate_states: Vec<[f64; 2]>,
    branch_impulses: Vec<f64>,
    candidate_impulses: Vec<f64>,
    stage_impulses: Vec<f64>,
    last_report: Option<NetworkWallReport>,
    trial_report: Option<NetworkWallReport>,
    trial_active: bool,
    transport: Option<CircularWallTransport>,
    #[cfg(feature = "material")]
    caloric: Vec<Option<MaterialWallCaloric>>,
}

impl NetworkWallOperator {
    #[allow(clippy::too_many_arguments)]
    fn exchange_for(
        &self, index: usize, state: Cell<'_>, volume_m3: f64, duration: Time,
        contact: WallContact, controls: WallControl, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<crate::multicomponent::wall_exchange::WallUpdate, WallExchangeError> {
        #[cfg(feature = "material")]
        if let Some(wall) = &self.caloric[index] {
            return crate::multicomponent::wall_exchange::exchange_cell_with_heat_storage(
                &self.gas, state, volume_m3, duration, contact, wall, controls, cancelled,
            );
        }
        #[cfg(not(feature = "material"))]
        let _ = index;
        exchange_cell(&self.gas, state, volume_m3, duration, contact, controls, cancelled)
    }

    /// Prepare one contact per original branch cell and finite plenum.
    ///
    /// Conductance is the complete exposed `h*A` in W/K for that cell, not a
    /// face conductivity. Wall storage is a complete heat capacity in J/K.
    /// An adiabatic contact explicitly opts a cell out of thermal exchange.
    pub fn new(
        network: &MixtureNetwork1d,
        branch_contacts: Vec<Vec<WallContact>>,
        plenum_contacts: Vec<WallContact>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if branch_contacts.len() != network.branches.len()
            || plenum_contacts.len() != network.plenums.len()
        {
            return Err(WallExchangeError::Invalid("wall contact topology"));
        }
        let mut branches = Vec::new();
        branches
            .try_reserve_exact(network.branches.len())
            .map_err(|_| WallExchangeError::Invalid("wall shape allocation"))?;
        let mut count = 0_usize;
        for (branch, contacts) in network.branches.iter().zip(&branch_contacts) {
            poll(&mut cancelled)?;
            if contacts.len() != branch.cell_count() {
                return Err(WallExchangeError::Invalid("wall contact cell count"));
            }
            branches.push(BranchShape {
                mesh: branch.mesh.clone(),
                start: count,
            });
            count = count
                .checked_add(branch.cell_count())
                .ok_or(WallExchangeError::Invalid("wall cell count overflow"))?;
            for contact in contacts {
                poll(&mut cancelled)?;
                contact.validate()?;
            }
        }
        let plenum_start = count;
        count = count
            .checked_add(plenum_contacts.len())
            .ok_or(WallExchangeError::Invalid("wall storage count overflow"))?;
        for contact in &plenum_contacts {
            poll(&mut cancelled)?;
            contact.validate()?;
            if contact.drag != WallDrag::None {
                return Err(WallExchangeError::Invalid(
                    "stationary plenums cannot carry a bulk drag law",
                ));
            }
        }
        let empty = WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        };
        let mut contacts = allocated(count, empty)?;
        for (shape, input) in branches.iter().zip(&branch_contacts) {
            for (destination, source) in contacts[shape.start..shape.start + input.len()]
                .chunks_mut(128)
                .zip(input.chunks(128))
            {
                poll(&mut cancelled)?;
                destination.copy_from_slice(source);
            }
        }
        contacts[plenum_start..].copy_from_slice(&plenum_contacts);
        let mut plenum_volumes = allocated(network.plenums.len(), 0.0)?;
        for (value, plenum) in plenum_volumes.iter_mut().zip(&network.plenums) {
            poll(&mut cancelled)?;
            *value = plenum.volume_m3;
        }
        let branch_count = branches.len();
        let result = Self {
            gas: network.gas().clone(),
            branches,
            plenum_volumes,
            plenum_start,
            contacts,
            candidate_contacts: allocated(count, empty)?,
            candidate_states: allocated(count, [0.0; 2])?,
            branch_impulses: allocated(branch_count, 0.0)?,
            candidate_impulses: allocated(branch_count, 0.0)?,
            stage_impulses: allocated(branch_count, 0.0)?,
            last_report: None,
            trial_report: None,
            trial_active: false,
            transport: None,
            #[cfg(feature = "material")]
            caloric: allocated(count, None)?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Original source and caloric interval used by every local exchange.
    pub fn gas(&self) -> &MixtureGas {
        &self.gas
    }

    /// Accepted wall states along one branch in its original mesh order.
    /// Bound transport branches report the coefficients used in the last accepted
    /// source half; each following half resolves them again from its actual gas.
    pub fn branch_contacts(&self, branch: usize) -> Option<&[WallContact]> {
        self.branches.get(branch).map(|shape| {
            &self.contacts[shape.start..shape.start + shape.mesh.cell_count()]
        })
    }

    /// Accepted wall states of finite plenums in retained topology order.
    pub fn plenum_contacts(&self) -> &[WallContact] {
        &self.contacts[self.plenum_start..]
    }

    /// Complete source report; failed steps leave it unchanged.
    pub const fn last_report(&self) -> Option<NetworkWallReport> {
        self.last_report
    }

    /// Support impulse per original branch axis [N s].
    ///
    /// These values are summed over both accepted wall halves when composed.
    /// Differently oriented branch axes are deliberately not combined.
    pub fn accepted_branch_impulses(&self) -> Option<&[f64]> {
        self.last_report.map(|_| self.branch_impulses.as_slice())
    }

    /// Apply one complete source interval without advancing the flow clock.
    ///
    /// Species remain bitwise unchanged. On success an old advection report is
    /// cleared because it does not describe the new physical fields. On refusal
    /// all accepted gas fields, wall temperatures, clocks and reports survive.
    pub fn advance_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: WallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkWallReport, WallExchangeError> {
        self.begin_trial(&mut cancelled)?;
        let report =
            self.advance_for_trial(network, duration, controls, &mut cancelled)?;
        // advance_for_trial already performed the final cancellation check.
        // Publication here is infallible and contains no user callback.
        self.publish_trial_report();
        Ok(report)
    }

    pub(super) fn begin_trial(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        self.trial_active = false;
        self.trial_report = None;
        for (destination, source) in self
            .candidate_contacts
            .chunks_mut(128)
            .zip(self.contacts.chunks(128))
        {
            poll(cancelled)?;
            destination.copy_from_slice(source);
        }
        for impulse in &mut self.candidate_impulses {
            poll(cancelled)?;
            *impulse = 0.0;
        }
        poll(cancelled)?;
        self.trial_active = true;
        Ok(())
    }

    pub(super) fn advance_for_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkWallReport, WallExchangeError> {
        poll(cancelled)?;
        controls.validate()?;
        if !self.trial_active
            || !duration.value().is_finite()
            || duration.value() < 0.0
            || self.gas != *network.gas()
            || self.branches.len() != network.branches.len()
            || self.plenum_volumes.len() != network.plenums.len()
        {
            return Err(WallExchangeError::Invalid(
                "wall interval, source or network ownership",
            ));
        }
        let n = self.gas.species_count();
        let width = n + 2;
        let mut balance = Balance::new();
        for (branch_index, (branch, shape)) in network
            .branches
            .iter()
            .zip(&self.branches)
            .enumerate()
        {
            poll(cancelled)?;
            if branch.gas != self.gas || branch.mesh != shape.mesh {
                return Err(WallExchangeError::Invalid("wall branch geometry changed"));
            }
            let mut impulse = 0.0;
            let mut correction = 0.0;
            for (i, values) in branch.states.chunks_exact(width).enumerate() {
                let index = shape.start + i;
                let state = Cell {
                    partial_densities: &values[..n],
                    momentum: values[n],
                    energy: values[n + 1],
                };
                let contact = match self.transport.as_mut() {
                    Some(transport) => transport.resolve(
                        &self.gas, branch_index, state, shape.mesh.volumes()[i],
                        self.candidate_contacts[index], cancelled,
                    )?,
                    None => self.candidate_contacts[index],
                };
                let update = self.exchange_for(
                    index,
                    state,
                    shape.mesh.volumes()[i],
                    duration,
                    contact,
                    controls,
                    &mut *cancelled,
                )?;
                if let Some(transport) = self.transport.as_mut() {
                    // A frozen-coefficient source step must also leave an
                    // endpoint supported by the original material/domain.
                    transport.resolve(
                        &self.gas, branch_index,
                        Cell { partial_densities: state.partial_densities,
                            momentum: update.momentum, energy: update.energy },
                        shape.mesh.volumes()[i], update.contact, cancelled,
                    )?;
                }
                self.candidate_states[index] = [update.momentum, update.energy];
                self.candidate_contacts[index] = update.contact;
                balance.add(update.report)?;
                accumulate(&mut impulse, &mut correction, update.report.wall_impulse_ns)?;
            }
            self.stage_impulses[branch_index] = impulse;
        }
        for (i, plenum) in network.plenums.iter().enumerate() {
            poll(cancelled)?;
            if plenum.volume_m3.to_bits() != self.plenum_volumes[i].to_bits() {
                return Err(WallExchangeError::Invalid("wall plenum geometry changed"));
            }
            let index = self.plenum_start + i;
            let update = self.exchange_for(
                index,
                plenum.state(),
                plenum.volume_m3,
                duration,
                self.candidate_contacts[index],
                controls,
                &mut *cancelled,
            )?;
            if update.momentum != 0.0 || update.report.wall_impulse_ns != 0.0 {
                return Err(WallExchangeError::Invalid("wall changed plenum momentum"));
            }
            self.candidate_states[index] = [0.0, update.energy];
            self.candidate_contacts[index] = update.contact;
            balance.add(update.report)?;
        }
        let report = balance.report(duration);
        let complete = match self.trial_report {
            Some(previous) => combine(previous, report)?,
            None => report,
        };
        for (total, stage) in self.candidate_impulses.iter_mut().zip(&self.stage_impulses) {
            poll(cancelled)?;
            let next = *total + *stage;
            if !next.is_finite() {
                return Err(WallExchangeError::Invalid("wall impulse accumulation"));
            }
            *total = next;
        }
        poll(cancelled)?;
        // All potentially failing work is finished. No cancellation or source
        // callback may leave a prefix of the network published.
        for (branch, shape) in network.branches.iter_mut().zip(&self.branches) {
            for (i, values) in branch.states.chunks_exact_mut(width).enumerate() {
                let updated = self.candidate_states[shape.start + i];
                values[n] = updated[0];
                values[n + 1] = updated[1];
            }
        }
        for (i, plenum) in network.plenums.iter_mut().enumerate() {
            plenum.state.energy = self.candidate_states[self.plenum_start + i][1];
        }
        network.last_report = None;
        self.trial_report = Some(complete);
        Ok(report)
    }

    pub(super) fn publish_trial_report(&mut self) {
        if let Some(report) = self.trial_report {
            core::mem::swap(&mut self.contacts, &mut self.candidate_contacts);
            core::mem::swap(&mut self.branch_impulses, &mut self.candidate_impulses);
            self.last_report = Some(report);
        }
        self.trial_active = false;
    }
}

struct Balance {
    values: [f64; 5],
    corrections: [f64; 5],
    maximum_residual: f64,
    maximum_iterations: usize,
}
impl Balance {
    fn new() -> Self {
        Self {
            values: [0.0; 5],
            corrections: [0.0; 5],
            maximum_residual: 0.0,
            maximum_iterations: 0,
        }
    }
    fn add(&mut self, report: WallCellReport) -> Result<(), WallExchangeError> {
        let terms = [
            report.heat_to_gas_j,
            report.wall_energy_change_j,
            report.external_heat_to_gas_j,
            report.drag_heating_j,
            report.energy_residual_j,
        ];
        for (i, term) in terms.into_iter().enumerate() {
            accumulate(&mut self.values[i], &mut self.corrections[i], term)?;
        }
        if !report.thermal_residual_j.is_finite() {
            return Err(WallExchangeError::Invalid("wall thermal residual"));
        }
        self.maximum_residual =
            self.maximum_residual.max(report.thermal_residual_j.abs());
        self.maximum_iterations = self.maximum_iterations.max(report.iterations);
        Ok(())
    }
    fn report(&self, elapsed: Time) -> NetworkWallReport {
        NetworkWallReport {
            elapsed,
            heat_to_gas_j: self.values[0],
            wall_energy_change_j: self.values[1],
            external_heat_to_gas_j: self.values[2],
            drag_heating_j: self.values[3],
            energy_residual_j: self.values[4],
            maximum_thermal_residual_j: self.maximum_residual,
            maximum_iterations: self.maximum_iterations,
        }
    }
}

fn combine(
    a: NetworkWallReport,
    b: NetworkWallReport,
) -> Result<NetworkWallReport, WallExchangeError> {
    let values = [
        a.elapsed.value() + b.elapsed.value(),
        sum([a.heat_to_gas_j, b.heat_to_gas_j]),
        sum([a.wall_energy_change_j, b.wall_energy_change_j]),
        sum([a.external_heat_to_gas_j, b.external_heat_to_gas_j]),
        sum([a.drag_heating_j, b.drag_heating_j]),
        sum([a.energy_residual_j, b.energy_residual_j]),
    ];
    if !values.iter().all(|x| x.is_finite()) {
        return Err(WallExchangeError::Invalid("composed wall balance"));
    }
    Ok(NetworkWallReport {
        elapsed: Time::new(values[0]),
        heat_to_gas_j: values[1],
        wall_energy_change_j: values[2],
        external_heat_to_gas_j: values[3],
        drag_heating_j: values[4],
        energy_residual_j: values[5],
        maximum_thermal_residual_j: a
            .maximum_thermal_residual_j
            .max(b.maximum_thermal_residual_j),
        maximum_iterations: a.maximum_iterations.max(b.maximum_iterations),
    })
}

fn accumulate(
    total: &mut f64,
    correction: &mut f64,
    term: f64,
) -> Result<(), WallExchangeError> {
    let y = term - *correction;
    let next = *total + y;
    if !term.is_finite() || !next.is_finite() {
        return Err(WallExchangeError::Invalid("wall extensive balance"));
    }
    *correction = (next - *total) - y;
    *total = next;
    Ok(())
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WallExchangeError> {
    if cancelled() {
        Err(WallExchangeError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
