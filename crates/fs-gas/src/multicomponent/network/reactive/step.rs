//! One complete transaction for chemistry, waves and any requested sources.
use super::*;

impl ReactingNetwork1d {
    pub(super) fn step_impl<L: TransportLaw>(
        &mut self, maximum_duration: Time, controls: ReactiveStepControl,
        mut operator: Option<&mut NetworkDissipationOperator<L>>,
        dissipation_controls: DissipationControl, mut sources: NetworkSources<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveNetworkReport, ReactiveNetworkError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value()) || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl) || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        { return Err(ReactiveNetworkError::Invalid("reacting network duration, CFL or retry controls")); }
        sources.validate()?;
        let before = checked_totals(&self.accepted, &mut cancelled)?;
        let mut duration = initial_duration(&self.accepted, &mut self.trial,
            maximum_duration.value(), controls.transport.cfl, &mut cancelled)?;
        let mut last = ReactiveNetworkError::Invalid("no representable reacting network increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let time = self.time().value();
            if !positive(duration) || !positive(0.5 * duration)
                || !positive(time + duration - time) || !(time + duration).is_finite()
            { return Err(ReactiveNetworkError::Invalid("reacting network increment cannot advance the clock")); }
            reset(&mut self.trial, &self.accepted, &mut cancelled)?;
            let mut chemistry = ChemistryReport::default();
            let attempt = (|| {
                let half = Time::new(0.5 * duration);
                let first_sources = sources.first_half(&mut self.trial, half, &mut cancelled)?;
                let first = match operator.as_deref_mut() {
                    Some(operator) => {
                        let report = operator.advance_for_trial(&mut self.trial, half,
                            dissipation_controls, &mut cancelled)?;
                        for (destination, source) in self.trial_molecular_ports.iter_mut()
                            .zip(operator.trial_port_exchanges())
                        {
                            poll(&mut cancelled)?;
                            *destination = *source;
                        }
                        Some(report)
                    }
                    None => None,
                };
                react_network(&mut self.trial, &mut self.reactor, &mut self.amounts,
                    half.value(), &mut chemistry, &mut cancelled)?;
                let transport = self.trial.step(Time::new(duration), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != duration {
                    return Ok((transport, false, None, None, None));
                }
                react_network(&mut self.trial, &mut self.reactor, &mut self.amounts,
                    half.value(), &mut chemistry, &mut cancelled)?;
                let second = match operator.as_deref_mut() {
                    Some(operator) => {
                        let report = operator.advance_for_trial(&mut self.trial, half,
                            dissipation_controls, &mut cancelled)?;
                        for (destination, source) in self.trial_molecular_ports.iter_mut()
                            .zip(operator.trial_port_exchanges())
                        {
                            poll(&mut cancelled)?;
                            for k in 0..destination.species_count {
                                destination.species_masses_to_plenum[k] += source.species_masses_to_plenum[k];
                            }
                            destination.momentum_to_plenum += source.momentum_to_plenum;
                            destination.energy_to_plenum += source.energy_to_plenum;
                            if !destination.species_masses_to_plenum[..destination.species_count].iter().all(|x| x.is_finite())
                                || !destination.momentum_to_plenum.is_finite() || !destination.energy_to_plenum.is_finite()
                            { return Err(ReactiveNetworkError::Invalid("unrepresentable combined molecular port exchange")); }
                        }
                        Some(report)
                    }
                    None => None,
                };
                let second_sources = sources.second_half(&mut self.trial, half, &mut cancelled)?;
                // A final valve can cool the donor or heat/change the receiver
                // after chemistry has finished. Admission must include that final
                // source, not only the states at the chemical endpoints.
                if operator.is_some() || !sources.is_empty() {
                    let chemistry_controls = self.reactor.controls();
                    for (branch, flow) in self.trial.branches.iter().enumerate() {
                        for cell in 0..flow.cell_count() {
                            admit_temperature(self.trial.gas(), flow.cell(cell).expect("admitted cell index"),
                                chemistry_controls, ChemistryLocation::Branch { branch, cell }, &mut cancelled)?;
                        }
                    }
                    for (plenum, volume) in self.trial.plenums.iter().enumerate() {
                        admit_temperature(self.trial.gas(), volume.state(), chemistry_controls,
                            ChemistryLocation::Plenum { plenum }, &mut cancelled)?;
                    }
                }
                let dissipation = first.zip(second).map(|(first, second)| SplitNetworkDissipationReport { first, second });
                let (wall, valves) = first_sources.pair(second_sources);
                Ok::<_, ReactiveNetworkError>((transport, true, dissipation, wall, valves))
            })();
            match attempt {
                Ok((transport, false, _, _, _)) => {
                    duration = transport.elapsed.value();
                    last = ReactiveNetworkError::Invalid("network transport shortened the coupled duration");
                }
                Ok((transport, true, dissipation, wall_report, valve_report)) => {
                    let after = checked_totals(&self.trial, &mut cancelled)?;
                    let n = self.gas().species_count();
                    let mut boundary_change = transport.boundary_change;
                    if let Some(d) = dissipation {
                        boundary_change.energy += d.first.boundary_change.energy + d.second.boundary_change.energy;
                        for k in 0..n {
                            boundary_change.species_masses[k] += d.first.boundary_change.species_masses[k]
                                + d.second.boundary_change.species_masses[k];
                        }
                    }
                    if let Some(w) = wall_report {
                        boundary_change.energy = sum([boundary_change.energy, w.first.heat_to_gas_j, w.second.heat_to_gas_j]);
                    }
                    // Valve transfers cancel internally. Adding their energy
                    // here would count internal donor enthalpy as external heat.
                    let mut residual = NetworkTotals { species_masses: [0.0; MAX_SPECIES], species_count: n,
                        energy: after.energy - before.energy - boundary_change.energy - chemistry.referenced_energy_change_j };
                    for k in 0..n {
                        residual.species_masses[k] = after.species_masses[k] - before.species_masses[k]
                            - boundary_change.species_masses[k] - chemistry.species_change_kg[k];
                    }
                    let mass_residual_kg = sum((0..n).map(|k| after.species_masses[k]
                        - before.species_masses[k] - boundary_change.species_masses[k]));
                    let absolute_energy_residual_j = sum(core::iter::once(after.energy - before.energy - boundary_change.energy)
                        .chain((0..n).map(|k| (after.species_masses[k] - before.species_masses[k]
                            - boundary_change.species_masses[k]) * self.gas().energy_references()[k])));
                    if !residual.energy.is_finite() || !residual.species_masses[..n].iter().all(|x| x.is_finite())
                        || !mass_residual_kg.is_finite() || !absolute_energy_residual_j.is_finite()
                        || !chemistry.entropy_change_j_k.is_finite()
                    { return Err(ReactiveNetworkError::Invalid("non-finite reacting network balance")); }
                    let report = ReactiveNetworkReport {
                        elapsed: Time::new(duration), coupling_retries: retry, transport, chemistry,
                        dissipation, wall: wall_report, valves: valve_report, before, after, boundary_change,
                        conservation_residual: residual, mass_residual_kg, absolute_energy_residual_j,
                    };
                    poll(&mut cancelled)?;
                    sources.publish();
                    if let Some(operator) = operator.as_deref_mut() {
                        operator.publish_trial_report();
                        core::mem::swap(&mut self.molecular_ports, &mut self.trial_molecular_ports);
                    }
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(ReactiveNetworkError::Cancelled) => return Err(ReactiveNetworkError::Cancelled),
                Err(error) => { last = error; duration *= 0.5; }
            }
        }
        Err(ReactiveNetworkError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1, last: Box::new(last),
        })
    }
}
