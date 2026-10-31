//! Complete shared-duration molecular, valve, wall and wave transaction.
use super::*;

impl<L: TransportLaw> DissipatingNetwork1d<L> {
    pub(super) fn step_impl(
        &mut self, maximum_duration: Time, controls: DissipatingStepControl,
        mut sources: NetworkSources<'_, '_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingNetworkError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value()) || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl) || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        { return Err(DissipatingNetworkError::Invalid("molecular network duration, CFL or retry controls")); }
        sources.validate()?;
        let before = checked_totals(&self.accepted, &mut cancelled)?;
        let mut duration = initial_duration(&self.accepted, &mut self.trial,
            maximum_duration.value(), controls.transport.cfl, &mut cancelled)?;
        let mut last = DissipatingNetworkError::Invalid("no representable network increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let time = self.time().value();
            if !positive(duration) || !positive(0.5 * duration)
                || !positive(time + duration - time) || !(time + duration).is_finite()
            { return Err(DissipatingNetworkError::Invalid("molecular network increment cannot advance the clock")); }
            reset(&mut self.trial, &self.accepted, &mut cancelled)?;
            let attempt = (|| {
                let half = Time::new(0.5 * duration);
                let first_sources = sources.first_half(&mut self.trial, half, &mut cancelled)?;
                let first = self.molecular.advance_for_trial(&mut self.trial, half, controls.molecular, &mut cancelled)?;
                for (destination, source) in self.trial_molecular_ports.iter_mut().zip(self.molecular.trial_port_exchanges()) {
                    poll(&mut cancelled)?;
                    *destination = *source;
                }
                let transport = match reservoirs {
                    Some(reservoirs) => self.trial.step_with_reservoirs(
                        Time::new(duration), controls.transport, reservoirs, &mut cancelled)?,
                    None => self.trial.step(
                        Time::new(duration), controls.transport, &mut cancelled)?,
                };
                if transport.elapsed.value() != duration { return Ok((transport, None)); }
                let second = self.molecular.advance_for_trial(&mut self.trial, half, controls.molecular, &mut cancelled)?;
                for (total, second) in self.trial_molecular_ports.iter_mut().zip(self.molecular.trial_port_exchanges()) {
                    poll(&mut cancelled)?;
                    total.momentum_to_plenum += second.momentum_to_plenum;
                    total.energy_to_plenum += second.energy_to_plenum;
                    for k in 0..total.species_count {
                        total.species_masses_to_plenum[k] += second.species_masses_to_plenum[k];
                    }
                    if !total.momentum_to_plenum.is_finite() || !total.energy_to_plenum.is_finite()
                        || !total.species_masses_to_plenum[..total.species_count].iter().all(|x| x.is_finite())
                    { return Err(DissipatingNetworkError::Invalid("non-finite composed molecular port exchange")); }
                }
                let second_sources = sources.second_half(&mut self.trial, half, &mut cancelled)?;
                let (wall, valves) = first_sources.pair(second_sources);
                Ok::<_, DissipatingNetworkError>((transport,
                    Some((SplitNetworkDissipationReport { first, second }, wall, valves))))
            })();
            match attempt {
                Ok((transport, None)) => {
                    duration = transport.elapsed.value();
                    last = DissipatingNetworkError::Invalid("network transport shortened the molecular composition");
                }
                Ok((transport, Some((molecular, wall_report, valve_report)))) => {
                    let after = checked_totals(&self.trial, &mut cancelled)?;
                    let n = self.gas().species_count();
                    let mut boundary_change = NetworkTotals { species_masses: [0.0; MAX_SPECIES], species_count: n,
                        energy: sum([molecular.first.boundary_change.energy, transport.boundary_change.energy,
                            molecular.second.boundary_change.energy]) };
                    if let Some(w) = wall_report {
                        boundary_change.energy = sum([boundary_change.energy, w.first.heat_to_gas_j, w.second.heat_to_gas_j]);
                    }
                    // Valve energy and species transfers are internal. Actual
                    // storage defects remain in the whole-step residual below.
                    let mut residual = NetworkTotals { species_masses: [0.0; MAX_SPECIES], species_count: n,
                        energy: after.energy - before.energy - boundary_change.energy };
                    for k in 0..n {
                        boundary_change.species_masses[k] = sum([molecular.first.boundary_change.species_masses[k],
                            transport.boundary_change.species_masses[k], molecular.second.boundary_change.species_masses[k]]);
                        residual.species_masses[k] = after.species_masses[k] - before.species_masses[k]
                            - boundary_change.species_masses[k];
                    }
                    let mass_residual_kg = residual.mass();
                    let absolute_energy_residual_j = sum(core::iter::once(residual.energy)
                        .chain((0..n).map(|k| residual.species_masses[k] * self.gas().energy_references()[k])));
                    if !boundary_change.energy.is_finite()
                        || !boundary_change.species_masses[..n].iter().all(|x| x.is_finite())
                        || !residual.energy.is_finite() || !residual.species_masses[..n].iter().all(|x| x.is_finite())
                        || !mass_residual_kg.is_finite() || !absolute_energy_residual_j.is_finite()
                    { return Err(DissipatingNetworkError::Invalid("non-finite complete molecular network balance")); }
                    let report = DissipatingNetworkReport { elapsed: Time::new(duration), coupling_retries: retry,
                        transport, molecular, wall: wall_report, valves: valve_report, before, after, boundary_change,
                        conservation_residual: residual, mass_residual_kg, absolute_energy_residual_j };
                    poll(&mut cancelled)?;
                    sources.publish();
                    self.molecular.publish_trial_report();
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    core::mem::swap(&mut self.molecular_ports, &mut self.trial_molecular_ports);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(DissipatingNetworkError::Cancelled) => return Err(DissipatingNetworkError::Cancelled),
                Err(error) => { last = error; duration *= 0.5; }
            }
        }
        Err(DissipatingNetworkError::RetryLimit { attempts: controls.maximum_coupling_retries + 1, last: Box::new(last) })
    }
}
