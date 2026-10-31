//! Clock-preserving valve sweeps shared by plain and reacting shock networks.
//!
//! The original nozzle and finite-blowdown solvers remain the only valve laws.
//! All plenums are staged before any source interval publishes. Coupled owners
//! retain trial receipts across both halves and publish them only with the gas.
//! Each source half holds its supplied chamber volumes fixed; a composing owner
//! may move them between halves and separately account for mechanical work.

use super::{Orifice, ValveControl, ValveError, ValveExchange, ValveLink, MAX_VALVE_LINKS, pair, poll};
use crate::multicomponent::network::{Accumulator, MixtureNetwork1d, NetworkTotals, Plenum};
use crate::multicomponent::{Conserved, MAX_SPECIES, MixtureGas, allocated, positive};
use fs_qty::Time;

/// Root budgets for a valve source interval; no acoustic or coupling controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveSolve {
    /// Absolute endpoint mass-equation allowance [kg].
    pub absolute_mass_tolerance_kg: f64,
    /// Relative allowance against the actual exchanged mass.
    pub relative_mass_tolerance: f64,
    /// Maximum scalar iterations per link, in 1..=128.
    pub maximum_iterations: usize,
}
impl Default for ValveSolve {
    fn default() -> Self { ValveControl::default().into() }
}
impl From<ValveControl> for ValveSolve {
    fn from(value: ValveControl) -> Self {
        Self {
            absolute_mass_tolerance_kg: value.absolute_mass_tolerance_kg,
            relative_mass_tolerance: value.relative_mass_tolerance,
            maximum_iterations: value.maximum_iterations,
        }
    }
}
impl ValveSolve {
    pub(in crate::multicomponent::network) fn validate(self) -> Result<(), ValveError> {
        if !positive(self.absolute_mass_tolerance_kg)
            || !positive(self.relative_mass_tolerance)
            || self.relative_mass_tolerance > 1e-3
            || !(1..=128).contains(&self.maximum_iterations)
        { return Err(ValveError::Invalid("valve source root controls")); }
        Ok(())
    }
    fn pair_controls(self) -> ValveControl {
        ValveControl {
            absolute_mass_tolerance_kg: self.absolute_mass_tolerance_kg,
            relative_mass_tolerance: self.relative_mass_tolerance,
            maximum_iterations: self.maximum_iterations,
            ..ValveControl::default()
        }
    }
}

/// Source-only work and actual storage defects. Every valve exchange is internal.
///
/// A composed report sums the separate valve intervals' defects, not the change
/// between their endpoints: intervening chemistry/advection may change storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValveSourceReport {
    /// Sum of the accepted source durations [s]. Does not advance the gas clock.
    pub elapsed: Time,
    /// Actual complete sweeps, including exact-zero-flow sweeps.
    pub sweeps: usize,
    /// Scalar mass iterations summed over all accepted links and sweeps.
    pub iterations: usize,
    /// Accepted choked link endpoints, not rejected root trials.
    pub choked_endpoints: usize,
    /// Largest absolute endpoint mass-equation defect [kg].
    pub maximum_mass_residual_kg: f64,
    /// Actual source-only net species and referenced-energy storage defects.
    pub conservation_residual: NetworkTotals,
    /// Source-only physical-energy defect including original species references [J].
    pub absolute_energy_residual_j: f64,
}

/// The two valve half sweeps around a complete source/transport composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitNetworkValveReport {
    /// Forward-order half sweep before the inner composition.
    pub first: ValveSourceReport,
    /// Reverse-order half sweep after the inner composition.
    pub second: ValveSourceReport,
}

/// Prepared source operator over the original plenum axis and gas source.
///
/// Preparation allocates one plenum scratch bank, an accepted-volume binding and
/// two receipt banks. Every successful sweep reuses them. Changes to branch meshes
/// are irrelevant to this compact-volume operator; the exact gas source, plenum
/// count and accepted volumes must match. A moving composition adopts new volume
/// bindings only when its complete step publishes. Equal volumes do not identify
/// physical tanks: retain the original ordering, just as for branch connections.
pub struct NetworkValveOperator {
    gas: MixtureGas,
    links: Vec<ValveLink>,
    candidate: Vec<Plenum>,
    accepted_plenum_volumes: Vec<f64>,
    exchanges: Vec<ValveExchange>,
    trial_exchanges: Vec<ValveExchange>,
    last_report: Option<ValveSourceReport>,
    trial_report: Option<ValveSourceReport>,
    trial_active: bool,
}
impl NetworkValveOperator {
    /// Prepare finite-volume valve links without taking ownership of the network.
    pub fn new(
        network: &MixtureNetwork1d,
        links: Vec<ValveLink>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ValveError> {
        poll(&mut cancelled)?;
        if links.len() > MAX_VALVE_LINKS { return Err(ValveError::Invalid("valve link count")); }
        for link in &links {
            poll(&mut cancelled)?;
            if link.first == link.second || link.first >= network.plenums().len()
                || link.second >= network.plenums().len()
            { return Err(ValveError::Invalid("valve plenum topology")); }
        }
        let mut candidate = Vec::new();
        candidate.try_reserve_exact(network.plenums().len())
            .map_err(|_| ValveError::Invalid("valve plenum allocation"))?;
        let mut accepted_plenum_volumes = allocated(network.plenums().len(), 0.0)?;
        for (index, p) in network.plenums().iter().enumerate() {
            poll(&mut cancelled)?;
            accepted_plenum_volumes[index] = p.volume_m3();
            let mut partial_densities = allocated(network.gas().species_count(), 0.0)?;
            partial_densities.copy_from_slice(p.state().partial_densities);
            candidate.push(Plenum::new(p.volume_m3(), Conserved {
                partial_densities, momentum: 0.0, energy: p.state().energy,
            })?);
        }
        let empty = ValveExchange::zero(network.gas().species_count());
        let exchanges = allocated(links.len(), empty)?;
        let trial_exchanges = allocated(links.len(), empty)?;
        poll(&mut cancelled)?;
        Ok(Self { gas: network.gas().clone(), links, candidate, accepted_plenum_volumes, exchanges, trial_exchanges,
            last_report: None, trial_report: None, trial_active: false })
    }
    /// Original source bank, caloric domain and energy gauges.
    pub fn gas(&self) -> &MixtureGas { &self.gas }
    /// Original deterministic valve order and current settings.
    pub fn links(&self) -> &[ValveLink] { &self.links }
    /// Last accepted source-only report, unchanged by a failed composition.
    pub const fn last_report(&self) -> Option<ValveSourceReport> { self.last_report }
    /// Actual accepted link exchanges, summed over both halves when composed.
    pub fn accepted_exchanges(&self) -> Option<&[ValveExchange]> {
        self.last_report.map(|_| self.exchanges.as_slice())
    }
    /// Change an opening exactly between accepted steps. Historical receipts
    /// remain historical; a refused index leaves both settings and receipts intact.
    pub fn set_orifice(&mut self, index: usize, orifice: Orifice) -> Result<(), ValveError> {
        let link = self.links.get_mut(index).ok_or(ValveError::Invalid("valve index"))?;
        link.orifice = orifice;
        self.trial_active = false;
        Ok(())
    }
    /// Complete one forward source sweep atomically, without advancing any clock.
    /// Zero duration is an exact state identity. The original pair solve retains
    /// its first-order, ordered-link and compact-volume model limitations.
    pub fn advance_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: ValveSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ValveSourceReport, ValveError> {
        self.begin_trial(&mut cancelled)?;
        let report = self.advance_for_trial(network, duration, controls, false, &mut cancelled)?;
        self.publish_trial_report();
        Ok(report)
    }
    pub(in crate::multicomponent::network) fn begin_trial(
        &mut self, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ValveError> {
        self.trial_active = false;
        self.trial_report = None;
        let empty = ValveExchange::zero(self.gas.species_count());
        for chunk in self.trial_exchanges.chunks_mut(128) {
            poll(cancelled)?;
            chunk.fill(empty);
        }
        poll(cancelled)?;
        self.trial_active = true;
        Ok(())
    }
    /// Authorize motion only from the last whole-step accepted geometry.
    pub(in crate::multicomponent::network) fn admit_plenum_volumes(
        &self,
        network: &MixtureNetwork1d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ValveError> {
        if self.accepted_plenum_volumes.len() != network.plenums().len() {
            return Err(ValveError::Invalid("valve plenum count changed"));
        }
        for (&accepted, source) in self.accepted_plenum_volumes.iter().zip(&network.plenums) {
            poll(cancelled)?;
            if accepted.to_bits() != source.volume_m3().to_bits() {
                return Err(ValveError::Invalid("valve plenum volume changed"));
            }
        }
        Ok(())
    }
    /// An unsuccessful interval must be discarded by the composing owner and
    /// followed by begin_trial, never resumed from a partially computed sweep.
    pub(in crate::multicomponent::network) fn advance_for_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: ValveSolve,
        reverse: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ValveSourceReport, ValveError> {
        self.admit_plenum_volumes(network, cancelled)?;
        self.advance_for_moving_trial(network, duration, controls, reverse, cancelled)
    }
    /// The composing owner must first admit its accepted volumes, then supply
    /// only geometry reached by its prescribed motion. Every valve interval is
    /// isochoric at that stage's actual volumes: pressure work belongs to the
    /// intervening moving-transport stage, never to this internal exchange.
    pub(in crate::multicomponent::network) fn advance_for_moving_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: ValveSolve,
        reverse: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ValveSourceReport, ValveError> {
        poll(cancelled)?;
        controls.validate()?;
        if !self.trial_active || !duration.value().is_finite() || duration.value() < 0.0
            || self.gas != *network.gas() || self.candidate.len() != network.plenums().len()
        { return Err(ValveError::Invalid("valve interval, source or plenum ownership")); }
        let n = self.gas.species_count();
        let mut before = Accumulator::new(n);
        for (candidate, source) in self.candidate.iter_mut().zip(&network.plenums) {
            poll(cancelled)?;
            if !positive(source.volume_m3()) {
                return Err(ValveError::Invalid("valve trial plenum volume"));
            }
            candidate.volume_m3 = source.volume_m3();
            before.add_cell(source.state(), source.volume_m3())?;
            candidate.state.partial_densities.copy_from_slice(source.state().partial_densities);
            candidate.state.energy = source.state().energy;
        }
        let mut iterations = 0;
        let mut choked = 0;
        let mut maximum = 0.0_f64;
        if duration.value() > 0.0 {
            for j in 0..self.links.len() {
                poll(cancelled)?;
                let index = if reverse { self.links.len() - 1 - j } else { j };
                let link = self.links[index];
                let update = pair::solve(&self.gas,
                    [&self.candidate[link.first], &self.candidate[link.second]],
                    link.orifice, duration.value(), controls.pair_controls(), cancelled)?;
                self.trial_exchanges[index].add(update.exchange)?;
                iterations += update.exchange.iterations;
                choked += update.exchange.choked_endpoints;
                maximum = maximum.max(update.exchange.maximum_mass_residual_kg);
                for (side, index) in [link.first, link.second].into_iter().enumerate() {
                    self.candidate[index].state.partial_densities.copy_from_slice(&update.states[side][..n]);
                    self.candidate[index].state.energy = update.states[side][n + 1];
                }
            }
        }
        let mut after = Accumulator::new(n);
        for p in &self.candidate {
            poll(cancelled)?;
            after.add_cell(p.state(), p.volume_m3())?;
        }
        let before = before.totals();
        let after = after.totals();
        let mut residual = NetworkTotals { species_masses: [0.0; MAX_SPECIES], species_count: n,
            energy: after.energy - before.energy };
        for k in 0..n { residual.species_masses[k] = after.species_masses[k] - before.species_masses[k]; }
        let report = ValveSourceReport { elapsed: duration, sweeps: 1, iterations,
            choked_endpoints: choked, maximum_mass_residual_kg: maximum,
            absolute_energy_residual_j: residual.absolute_energy(&self.gas)?, conservation_residual: residual };
        let complete = match self.trial_report {
            None => report,
            Some(old) => combine(old, report)?,
        };
        poll(cancelled)?;
        // Everything after the last poll is infallible. No prefix of the plenum
        // bank or accepted receipt bank can become visible on refusal.
        for (destination, source) in network.plenums.iter_mut().zip(&self.candidate) {
            destination.state.partial_densities.copy_from_slice(source.state().partial_densities);
            destination.state.energy = source.state().energy;
        }
        network.last_report = None;
        self.trial_report = Some(complete);
        Ok(report)
    }
    pub(in crate::multicomponent::network) fn publish_trial_report(&mut self) {
        if self.trial_active {
            if let Some(report) = self.trial_report {
                for (accepted, candidate) in self.accepted_plenum_volumes.iter_mut().zip(&self.candidate) {
                    *accepted = candidate.volume_m3();
                }
                core::mem::swap(&mut self.exchanges, &mut self.trial_exchanges);
                self.last_report = Some(report);
            }
        }
        self.trial_active = false;
    }
}
fn combine(mut a: ValveSourceReport, b: ValveSourceReport) -> Result<ValveSourceReport, ValveError> {
    a.elapsed = Time::new(a.elapsed.value() + b.elapsed.value());
    a.sweeps = a.sweeps.checked_add(b.sweeps).ok_or(ValveError::Invalid("valve sweep count overflow"))?;
    a.iterations = a.iterations.checked_add(b.iterations).ok_or(ValveError::Invalid("valve iteration count overflow"))?;
    a.choked_endpoints = a.choked_endpoints.checked_add(b.choked_endpoints).ok_or(ValveError::Invalid("valve choke count overflow"))?;
    a.maximum_mass_residual_kg = a.maximum_mass_residual_kg.max(b.maximum_mass_residual_kg);
    a.absolute_energy_residual_j += b.absolute_energy_residual_j;
    a.conservation_residual.energy += b.conservation_residual.energy;
    for k in 0..a.conservation_residual.species_count {
        a.conservation_residual.species_masses[k] += b.conservation_residual.species_masses[k];
    }
    if !a.elapsed.value().is_finite() || !a.absolute_energy_residual_j.is_finite()
        || !a.conservation_residual.energy.is_finite()
        || !a.conservation_residual.species_masses.iter().all(|x| x.is_finite())
    { return Err(ValveError::Invalid("composed valve source balance overflow")); }
    Ok(a)
}

#[cfg(test)]
mod tests;
