//! One joint finite-photon history, with external and recycled energy distinct.

use super::super::{Accounts, Pending, annihilation::PreparedAnnihilation, partition};
use super::*;
use fs_thermal::enthalpy::CellId;

#[allow(clippy::too_many_arguments)]
pub(in crate::volumes::heating::luminescence::saturation) fn step_sources<const M: usize>(
    owner: &mut SaturableLuminescentHeating<M>,
    thermal: &DepositionSimulation<'_>,
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    pumped: bool,
    prepared: &PreparedRecycling,
    pairs: Option<&PreparedAnnihilation>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, OpticalHeatingError> {
    let feedback = owner
        .recycling
        .as_ref()
        .expect("bound finite photon recycling");
    let mut limits = feedback.geometry.limits;
    let laws = match pairs {
        Some(pairs) => {
            limits.annihilation = pairs.limits;
            pairs.laws.clone()
        }
        None => vec![None; owner.emitters.len()],
    };
    let incident = match &owner.beam {
        Some(beam) => beam.incident_rates(start_s, pumped, checkpoint)?,
        None => Vec::new(),
    };
    let initial: Vec<_> = owner.emitters.iter().map(|e| e.populations).collect();
    let local_pumps = if pumped
        && owner
            .emitters
            .iter()
            .any(|emitter| emitter.binding.pump.rates_s() != [0.0; 2])
    {
        owner
            .emitters
            .iter()
            .map(|emitter| emitter.binding.pump.local_photon_pump(true))
            .collect()
    } else {
        Vec::new()
    };
    let step = preserve_checkpoint(checkpoint, |cancelled| {
        prepared.solver.advance_with_local_pumps(
            &initial,
            &incident,
            &local_pumps,
            duration_s,
            &laws,
            limits,
            cancelled,
        )
    })?;
    let beam_accounts = owner
        .beam
        .as_ref()
        .map(|beam| beam.recycling_accounts(&step, duration_s, owner.energy_tolerance_j))
        .transpose()?;
    if step.cells.len() != owner.emitters.len()
        || step.emission_paths.len() != prepared.path_bins.len()
    {
        return Err(invalid(
            "finite recycling result differs from original inventory or spectrum",
        ));
    }
    let mut accounts: Accounts<2, M> = owner.accounts.clone();
    accounts.duration.add(duration_s)?;
    let mut sources = BTreeMap::<CellId, CheckedSum>::new();
    let mut populations = Vec::with_capacity(owner.emitters.len());
    let mut mean_emission_w = Vec::with_capacity(owner.emitters.len());
    for (index, emitter) in owner.emitters.iter().enumerate() {
        checkpoint()?;
        let kinetics = &step.cells[index];
        // The native owner normalizes inactive local supplies to an empty
        // report. These direct events are external energy, not recycled light.
        let local_pump = step.local_pumps.get(index).copied().unwrap_or_default();
        let mut heat = CheckedSum::default();
        let mut local = CheckedSum::default();
        local.add(kinetics.residual_j)?;
        let (singlet_light, _) = partition(&owner.singlet, kinetics.radiative_decay_energy_j[0])?;
        let (triplet_light, _) = partition(
            owner.triplet.triplet(),
            kinetics.radiative_decay_energy_j[1],
        )?;
        let mut mean = [0.0; M];
        for state in 0..2 {
            let source_heat = step.source_relaxation_heat_j[index][state];
            accounts.relaxation.add(source_heat)?;
            heat.add(source_heat)?;
            local.add(kinetics.radiative_decay_energy_j[state])?;
            local.add(-source_heat)?;
            accounts.absorbed[state].add(step.external_absorbed_pump_energy_j[index][state])?;
            accounts.absorbed[state].add(local_pump.absorbed_energy_j[state])?;
            accounts
                .prompt
                .add(step.external_pump_relaxation_heat_j[index][state])?;
            accounts.prompt.add(local_pump.relaxation_heat_j[state])?;
            accounts
                .secondary
                .add(step.recycled_pump_relaxation_heat_j[index][state])?;
            accounts.recycled.add(finite(
                step.recycled_pump_events[index][state]
                    * prepared.solver.models()[index].excitation_energies_j()[state],
            )?)?;
            heat.add(kinetics.pump_relaxation_heat_j[state])?;
            accounts.dark.add(kinetics.nonradiative_heat_j[state])?;
            heat.add(kinetics.nonradiative_heat_j[state])?;
        }
        for bin in 0..M {
            let emitted = finite(singlet_light[bin] + triplet_light[bin])?;
            accounts.emitted[bin].add(emitted)?;
            local.add(-emitted)?;
            mean[bin] = finite(emitted / duration_s)?;
        }
        accounts.annihilation.add(kinetics.annihilation_heat_j)?;
        accounts.crossing.add(kinetics.crossing_heat_j)?;
        accounts
            .reverse_crossing
            .add(kinetics.reverse_crossing_heat_absorbed_j)?;
        accounts
            .passive_pump
            .add(step.external_passive_heat_j[index])?;
        accounts
            .secondary
            .add(step.recycled_passive_heat_j[index])?;
        for value in [
            kinetics.annihilation_heat_j,
            kinetics.crossing_heat_j,
            -kinetics.reverse_crossing_heat_absorbed_j,
            step.passive_heat_j[index],
        ] {
            heat.add(value)?;
        }
        accounts.roundoff.add(local.value()?.abs())?;
        sources
            .entry(emitter.binding.thermal_cell)
            .or_default()
            .add(finite(heat.value()? / duration_s)?)?;
        populations.push(kinetics.populations);
        mean_emission_w.push((emitter.binding.grid_cell, mean));
    }
    for ((path, bins), declaration) in step
        .emission_paths
        .iter()
        .zip(&prepared.path_bins)
        .zip(prepared.solver.emission_paths())
    {
        checkpoint()?;
        for lane in 0..2 {
            if let Some(bin) = bins[lane] {
                // The ray reaches the original open boundary completely.
                accounts.escaped[bin].add(path.transmitted_energy_j[lane])?;
                accounts.roundoff.add(finite(
                    path.photon_residual[lane].abs() * declaration.path.photon_energies_j[lane],
                )?)?;
                let expected = finite(
                    step.cells[declaration.source_cell].radiative_events[declaration.source_state]
                        * declaration.photon_fractions[lane],
                )?;
                accounts.roundoff.add(finite(
                    (expected - path.incident_photons[lane]).abs()
                        * declaration.path.photon_energies_j[lane],
                )?)?;
            }
        }
    }
    for (path, declaration) in step
        .external_beams
        .iter()
        .zip(prepared.solver.external_paths())
    {
        for lane in 0..2 {
            accounts.roundoff.add(finite(
                path.photon_residual[lane].abs() * declaration.photon_energies_j[lane],
            )?)?;
        }
    }
    if accounts.roundoff.value()? > owner.energy_tolerance_j {
        return Err(invalid(
            "finite recycling photon and local energy defects exceed allowance",
        ));
    }
    let mut stored = CheckedSum::default();
    for cell in &step.cells {
        for value in cell.stored_energy_j {
            stored.add(value)?;
        }
    }
    let (_, energy) = accounts.report(stored.value()?)?;
    if energy.residual_j.abs() > owner.energy_tolerance_j {
        return Err(invalid(
            "finite recycling external pump, excitation, heat and escape do not balance",
        ));
    }
    // The native integrator subsequently validates the endpoint material state
    // before accepting these populations, heat channels and beam accounts.
    let sources = sources
        .into_iter()
        .map(|(id, value)| Ok((id, value.value()?)))
        .collect::<Result<_, OpticalHeatingError>>()?;
    let _ = thermal; // The prepared data came from this exact native state.
    checkpoint()?;
    owner.pending = Some(Pending {
        start_s,
        duration_s,
        end_s,
        populations,
        mean_emission_w,
        accounts,
        beam_accounts,
    });
    Ok(sources)
}
