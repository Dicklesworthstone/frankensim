//! Shared physical-state transactions and optional conservative source halves.

use super::super::evolve::totals;
use super::super::{Conserved, Error, allocated, poll, positive};
use super::{Accumulator, BranchEnd, MixtureNetwork1d, NetworkTotals, Plenum};
use super::valve::{NetworkValveOperator, SplitNetworkValveReport, ValveError, ValveSolve, ValveSourceReport};
use super::wall_exchange::{NetworkWallOperator, NetworkWallReport, SplitNetworkWallReport};
use super::super::wall_exchange::{WallControl, WallExchangeError};
use fs_qty::Time;

/// Optional finite-wall and valve sources for a complete coupled network step.
///
/// Construct operators from the raw network before moving it into a composing
/// owner. Each borrowed operator publishes its reports only when the whole step
/// succeeds. Openings are fixed over that accepted interval and may be changed
/// between calls. This bundle introduces no new physical law or allocation.
///
/// Source ordering is W/2, V/2, inner composition, V/2, W/2. Valve links use
/// forward then reverse retained order. Symmetric ordering does not promote the
/// first-order source solvers to a certified second-order method.
#[derive(Default)]
pub struct NetworkSources<'a, 'law> {
    /// Optional stationary drag, reservoir heat and finite-wall heat storage.
    pub wall: Option<(&'a mut NetworkWallOperator<'law>, WallControl)>,
    /// Optional original-source finite-plenum valve transport.
    pub valves: Option<(&'a mut NetworkValveOperator, ValveSolve)>,
}

pub(super) enum SourceError {
    Wall(WallExchangeError),
    Valve(ValveError),
}
impl From<WallExchangeError> for SourceError {
    fn from(error: WallExchangeError) -> Self { Self::Wall(error) }
}
impl From<ValveError> for SourceError {
    fn from(error: ValveError) -> Self { Self::Valve(error) }
}

#[derive(Clone, Copy, Default)]
pub(super) struct SourceStage {
    wall: Option<NetworkWallReport>,
    valves: Option<ValveSourceReport>,
}
impl SourceStage {
    pub(super) fn pair(self, second: Self) -> (Option<SplitNetworkWallReport>, Option<SplitNetworkValveReport>) {
        (
            self.wall.zip(second.wall).map(|(first, second)| SplitNetworkWallReport { first, second }),
            self.valves.zip(second.valves).map(|(first, second)| SplitNetworkValveReport { first, second }),
        )
    }
}
impl NetworkSources<'_, '_> {
    pub(super) fn is_empty(&self) -> bool { self.wall.is_none() && self.valves.is_none() }
    pub(super) fn validate(&self) -> Result<(), SourceError> {
        if let Some((_, controls)) = &self.wall { controls.validate()?; }
        if let Some((_, controls)) = &self.valves { controls.validate()?; }
        Ok(())
    }
    pub(super) fn admit_volume_motion(
        &self,
        network: &MixtureNetwork1d,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), SourceError> {
        // Valve sweeps retain a separate rigid-volume ownership contract.
        // Refuse before any source trial rather than silently rebind that owner.
        if self.valves.is_some() {
            return Err(
                ValveError::Invalid("moving plenums with valve sources are not supported").into(),
            );
        }
        if let Some((wall, _)) = &self.wall {
            wall.admit_plenum_volumes(network, cancelled)?;
        }
        Ok(())
    }
    pub(super) fn first_half(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceStage, SourceError> {
        self.first_half_impl(network, duration, false, cancelled)
    }
    pub(super) fn first_moving_half(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceStage, SourceError> {
        self.first_half_impl(network, duration, true, cancelled)
    }
    fn first_half_impl(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        moving: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceStage, SourceError> {
        // Reset every source receipt/storage bank before any source changes the
        // trial gas. A later failure cannot publish a prefix of the composition.
        if let Some((op, _)) = &mut self.wall { op.begin_trial(cancelled)?; }
        if let Some((op, _)) = &mut self.valves { op.begin_trial(cancelled)?; }
        let wall = self.wall_half(network, duration, moving, cancelled)?;
        let valves = self.valve_half(network, duration, false, cancelled)?;
        Ok(SourceStage { wall, valves })
    }
    pub(super) fn second_half(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceStage, SourceError> {
        self.second_half_impl(network, duration, false, cancelled)
    }
    pub(super) fn second_moving_half(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceStage, SourceError> {
        self.second_half_impl(network, duration, true, cancelled)
    }
    fn second_half_impl(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        moving: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceStage, SourceError> {
        let valves = self.valve_half(network, duration, true, cancelled)?;
        let wall = self.wall_half(network, duration, moving, cancelled)?;
        Ok(SourceStage { wall, valves })
    }
    fn wall_half(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        moving: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<NetworkWallReport>, SourceError> {
        match &mut self.wall {
            Some((op, controls)) => Ok(Some(if moving {
                op.advance_for_moving_trial(network, duration, *controls, cancelled)?
            } else {
                op.advance_for_trial(network, duration, *controls, cancelled)?
            })),
            None => Ok(None),
        }
    }
    fn valve_half(
        &mut self, network: &mut MixtureNetwork1d, duration: Time, reverse: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<ValveSourceReport>, SourceError> {
        match &mut self.valves {
            Some((op, controls)) => Ok(Some(op.advance_for_trial(network, duration, *controls, reverse, cancelled)?)),
            None => Ok(None),
        }
    }
    // Called only after the composing owner's final cancellation/admission point.
    pub(super) fn publish(&mut self) {
        if let Some((op, _)) = &mut self.wall { op.publish_trial_report(); }
        if let Some((op, _)) = &mut self.valves { op.publish_trial_report(); }
    }
}

pub(super) fn checked_totals(
    network: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NetworkTotals, Error> {
    let n = network.gas().species_count();
    let mut accumulated = Accumulator::new(n);
    for branch in &network.branches {
        let value = totals(&branch.states, branch.mesh.volumes(), n, cancelled)?;
        accumulated.add(&value.species_masses[..n], value.energy)?;
    }
    for volume in &network.plenums {
        poll(cancelled)?;
        accumulated.add_cell(volume.state(), volume.volume_m3)?;
    }
    Ok(accumulated.totals())
}

pub(super) fn duplicate(
    network: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MixtureNetwork1d, Error> {
    let mut branches = Vec::new();
    branches.try_reserve_exact(network.branches.len()).map_err(|_| Error::Allocation)?;
    for branch in &network.branches { branches.push(branch.duplicate(cancelled)?); }
    let mut plenums = Vec::new();
    plenums.try_reserve_exact(network.plenums.len()).map_err(|_| Error::Allocation)?;
    for plenum in &network.plenums {
        poll(cancelled)?;
        let mut densities = allocated(plenum.state.partial_densities.len(), 0.0)?;
        densities.copy_from_slice(&plenum.state.partial_densities);
        plenums.push(Plenum::new(plenum.volume_m3, Conserved {
            partial_densities: densities, momentum: plenum.state.momentum, energy: plenum.state.energy,
        })?);
    }
    let mut connections = Vec::new();
    connections.try_reserve_exact(network.connections.len()).map_err(|_| Error::Allocation)?;
    connections.extend_from_slice(&network.connections);
    MixtureNetwork1d::new(branches, plenums, connections, cancelled)
}

// Reset physical fields and clocks in bounded tiles without allocation.
pub(super) fn reset(
    trial: &mut MixtureNetwork1d,
    accepted: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for (trial, accepted) in trial.branches.iter_mut().zip(&accepted.branches) {
        for (destination, source) in trial.states.chunks_mut(4096).zip(accepted.states.chunks(4096)) {
            poll(cancelled)?;
            destination.copy_from_slice(source);
        }
        trial.time = accepted.time;
    }
    for (trial, accepted) in trial.plenums.iter_mut().zip(&accepted.plenums) {
        poll(cancelled)?;
        trial.volume_m3 = accepted.volume_m3;
        trial.state.partial_densities.copy_from_slice(&accepted.state.partial_densities);
        trial.state.momentum = accepted.state.momentum;
        trial.state.energy = accepted.state.energy;
    }
    trial.last_report = None;
    Ok(())
}

// Bound expensive source attempts by original acoustic scales. Actual post-source
// face speeds remain the transport owner's responsibility; shortening replays all.
pub(super) fn initial_duration(
    accepted: &MixtureNetwork1d,
    trial: &mut MixtureNetwork1d,
    maximum: f64,
    cfl: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut rate = 0.0_f64;
    for branch in &accepted.branches {
        let width = branch.gas.species_count() + 2;
        for (i, q) in branch.states.chunks_exact(width).enumerate() {
            let p = branch.gas.decode(q, None, cancelled)?;
            let area = 0.5 * branch.mesh.face_areas()[i] + 0.5 * branch.mesh.face_areas()[i + 1];
            let local = (p.velocity.abs() + p.sound_speed) * (area / branch.mesh.volumes()[i]);
            if !positive(local) { return Err(Error::Invalid("unrepresentable coupled branch CFL")); }
            rate = rate.max(local);
        }
    }
    trial.work.plenum_signal_rates.fill(0.0);
    for (i, volume) in accepted.plenums.iter().enumerate() {
        trial.work.plenum_primitive[i] = accepted.gas().primitive(volume.state(), &mut *cancelled)?;
    }
    for connection in &accepted.connections {
        poll(cancelled)?;
        let branch = &accepted.branches[connection.branch];
        let face = match connection.end {
            BranchEnd::Left => 0,
            BranchEnd::Right => branch.cell_count(),
        };
        let local = 0.5 * branch.mesh.face_areas()[face]
            * (trial.work.plenum_primitive[connection.plenum].sound_speed
                / accepted.plenums[connection.plenum].volume_m3);
        let total = trial.work.plenum_signal_rates[connection.plenum] + local;
        if !positive(local) || !positive(total) { return Err(Error::Invalid("unrepresentable coupled plenum CFL")); }
        trial.work.plenum_signal_rates[connection.plenum] = total;
        rate = rate.max(total);
    }
    Ok(maximum.min(0.8 * cfl / rate))
}
