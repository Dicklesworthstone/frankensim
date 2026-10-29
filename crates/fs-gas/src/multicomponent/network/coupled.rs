//! Shared physical-state transactions for molecular and reactive compositions.

use super::super::evolve::totals;
use super::super::{Conserved, Error, allocated, poll, positive};
use super::{Accumulator, BranchEnd, MixtureNetwork1d, NetworkTotals, Plenum};

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
    branches
        .try_reserve_exact(network.branches.len())
        .map_err(|_| Error::Allocation)?;
    for branch in &network.branches {
        branches.push(branch.duplicate(cancelled)?);
    }
    let mut plenums = Vec::new();
    plenums
        .try_reserve_exact(network.plenums.len())
        .map_err(|_| Error::Allocation)?;
    for plenum in &network.plenums {
        poll(cancelled)?;
        let mut densities = allocated(plenum.state.partial_densities.len(), 0.0)?;
        densities.copy_from_slice(&plenum.state.partial_densities);
        plenums.push(Plenum::new(
            plenum.volume_m3,
            Conserved {
                partial_densities: densities,
                momentum: plenum.state.momentum,
                energy: plenum.state.energy,
            },
        )?);
    }
    let mut connections = Vec::new();
    connections
        .try_reserve_exact(network.connections.len())
        .map_err(|_| Error::Allocation)?;
    connections.extend_from_slice(&network.connections);
    MixtureNetwork1d::new(branches, plenums, connections, cancelled)
}

// These networks were constructed from identical immutable owners and topology.
// Reset only physical fields and clocks, in bounded tiles and without allocation.
pub(super) fn reset(
    trial: &mut MixtureNetwork1d,
    accepted: &MixtureNetwork1d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for (trial, accepted) in trial.branches.iter_mut().zip(&accepted.branches) {
        for (destination, source) in trial
            .states
            .chunks_mut(4096)
            .zip(accepted.states.chunks(4096))
        {
            poll(cancelled)?;
            destination.copy_from_slice(source);
        }
        trial.time = accepted.time;
    }
    for (trial, accepted) in trial.plenums.iter_mut().zip(&accepted.plenums) {
        poll(cancelled)?;
        trial
            .state
            .partial_densities
            .copy_from_slice(&accepted.state.partial_densities);
        trial.state.momentum = accepted.state.momentum;
        trial.state.energy = accepted.state.energy;
    }
    trial.last_report = None;
    Ok(())
}

// Bound expensive source attempts by original branch and compact-volume acoustic
// scales. The actual post-source characteristic face speeds remain the transport
// owner's responsibility; a shortened step restarts the complete composition.
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
            if !positive(local) {
                return Err(Error::Invalid("unrepresentable coupled branch CFL"));
            }
            rate = rate.max(local);
        }
    }
    trial.work.plenum_signal_rates.fill(0.0);
    for (i, volume) in accepted.plenums.iter().enumerate() {
        trial.work.plenum_primitive[i] =
            accepted.gas().primitive(volume.state(), &mut *cancelled)?;
    }
    for connection in &accepted.connections {
        poll(cancelled)?;
        let branch = &accepted.branches[connection.branch];
        let face = match connection.end {
            BranchEnd::Left => 0,
            BranchEnd::Right => branch.cell_count(),
        };
        let local = 0.5
            * branch.mesh.face_areas()[face]
            * (trial.work.plenum_primitive[connection.plenum].sound_speed
                / accepted.plenums[connection.plenum].volume_m3);
        let total = trial.work.plenum_signal_rates[connection.plenum] + local;
        if !positive(local) || !positive(total) {
            return Err(Error::Invalid("unrepresentable coupled plenum CFL"));
        }
        trial.work.plenum_signal_rates[connection.plenum] = total;
        rate = rate.max(total);
    }
    Ok(maximum.min(0.8 * cfl / rate))
}
