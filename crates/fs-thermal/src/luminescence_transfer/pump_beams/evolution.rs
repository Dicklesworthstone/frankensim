//! Ordered transport and bounded physical source iteration.

use super::accounts::BeamTotals;
use super::*;

struct Substep<const B: usize> {
    cells: Vec<cell::SharedTrial>,
    outgoing: Vec<[f64; B]>,
}

/// Optional attribution of actual channel events; ordinary external-only
/// callers allocate none of these additional accumulators.
pub(super) struct SourceAttribution {
    pub external_paths: usize,
    pub external_absorbed: Vec<[Sum; 2]>,
    pub recycled_absorbed: Vec<[Sum; 2]>,
    pub external_events: Vec<[Sum; 2]>,
    pub recycled_events: Vec<[Sum; 2]>,
    pub external_relaxation: Vec<[Sum; 2]>,
    pub recycled_relaxation: Vec<[Sum; 2]>,
    pub external_passive: Vec<Sum>,
    pub recycled_passive: Vec<Sum>,
}

impl SourceAttribution {
    pub fn new(cells: usize, external_paths: usize) -> Self {
        Self {
            external_paths,
            external_absorbed: vec![[Sum::default(); 2]; cells],
            recycled_absorbed: vec![[Sum::default(); 2]; cells],
            external_events: vec![[Sum::default(); 2]; cells],
            recycled_events: vec![[Sum::default(); 2]; cells],
            external_relaxation: vec![[Sum::default(); 2]; cells],
            recycled_relaxation: vec![[Sum::default(); 2]; cells],
            external_passive: vec![Sum::default(); cells],
            recycled_passive: vec![Sum::default(); cells],
        }
    }

    pub fn add_packet(
        &mut self,
        index: usize,
        beam: usize,
        packet: cell::ChannelStep,
        energy: f64,
        state_energies: [f64; 2],
    ) -> Result<(), SingletTripletError> {
        let (absorbed, events, relaxation, passive) = if beam < self.external_paths {
            (&mut self.external_absorbed[index], &mut self.external_events[index],
                &mut self.external_relaxation[index], &mut self.external_passive[index])
        } else {
            (&mut self.recycled_absorbed[index], &mut self.recycled_events[index],
                &mut self.recycled_relaxation[index], &mut self.recycled_passive[index])
        };
        for target in 0..2 {
            let count = packet.absorbed_by_state[target];
            events[target].add(count)?;
            absorbed[target].add(positive_product(count, energy, "attributed pump energy")?)?;
            if count > 0.0 {
                relaxation[target].add(positive_product(count, energy - state_energies[target], "attributed pump relaxation")?)?;
            }
        }
        passive.add(positive_product(packet.passive, energy, "attributed passive heat")?)?;
        Ok(())
    }
}

impl<const B: usize> SaturablePumpBeams<B> {
    fn channel(&self, visit: Visit, lane: usize, incident_rate: f64) -> cell::Channel {
        let path = &self.paths[visit.beam];
        let segment = &path.segments[visit.segment];
        cell::Channel {
            incident_rate,
            energy: path.photon_energies_j[lane],
            depth: segment.optical_depth_per_ground[lane],
            passive_depth: segment.passive_optical_depth[lane],
        }
    }

    fn trace(
        &self,
        means: &[f64],
        incident: &[[f64; B]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<Vec<cell::Channel>>, Vec<[f64; B]>), SingletTripletError> {
        let mut channels: Vec<_> = self.visits.iter()
            .map(|visits| Vec::with_capacity(visits.len() * B)).collect();
        let mut outgoing = incident.to_vec();
        for (beam, path) in self.paths.iter().enumerate() {
            for (ordinal, segment) in path.segments.iter().enumerate() {
                poll(cancelled)?;
                let visit = Visit { beam, segment: ordinal };
                for lane in 0..B {
                    let channel = self.channel(visit, lane, outgoing[beam][lane]);
                    outgoing[beam][lane] = cell::transmitted_rate(
                        channel, means[segment.cell_index],
                    )?;
                    channels[segment.cell_index].push(channel);
                }
            }
        }
        Ok((channels, outgoing))
    }

    #[allow(clippy::too_many_arguments)]
    fn substep(
        &self,
        initial: &[[f64; 3]],
        incident: &[[f64; B]],
        duration: f64,
        nonlinear: bool,
        limits: PumpBeamsLimits,
        evaluations: &mut usize,
        coupling_iterations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Substep<B>>, SingletTripletError> {
        // Exact dark/passive history, including cells absent from every path.
        // A geometrically cyclic network costs no coupling iterations after
        // all pumps switch off.
        if !nonlinear {
            let means: Vec<_> = initial.iter().map(|p| p[0]).collect();
            let (channels, outgoing) = self.trace(&means, incident, cancelled)?;
            let mut cells = Vec::with_capacity(self.models.len());
            for (index, &model) in self.models.iter().enumerate() {
                cells.push(cell::solve(
                    model, initial[index], &channels[index], duration,
                    limits.column, evaluations, cancelled,
                )?);
            }
            return Ok(Some(Substep { cells, outgoing }));
        }
        if let Some(order) = &self.topological {
            let mut outgoing = incident.to_vec();
            let mut trials: Vec<Option<cell::SharedTrial>> =
                (0..self.models.len()).map(|_| None).collect();
            for &index in order {
                poll(cancelled)?;
                let mut channels = Vec::with_capacity(self.visits[index].len() * B);
                for &visit in &self.visits[index] {
                    for lane in 0..B {
                        channels.push(self.channel(visit, lane, outgoing[visit.beam][lane]));
                    }
                }
                let trial = cell::solve(
                    self.models[index], initial[index], &channels,
                    duration, limits.column, evaluations, cancelled,
                )?;
                for (ordinal, visit) in self.visits[index].iter().enumerate() {
                    for lane in 0..B {
                        outgoing[visit.beam][lane] =
                            trial.channels[ordinal * B + lane].transmitted_rate;
                    }
                }
                trials[index] = Some(trial);
            }
            let cells = trials.into_iter().map(|trial| {
                trial.ok_or(SingletTripletError::Invalid("incomplete pump beam topology"))
            }).collect::<Result<_, _>>()?;
            return Ok(Some(Substep { cells, outgoing }));
        }

        let mut means: Vec<_> = initial.iter().map(|p| p[0]).collect();
        let mut previous_delta = vec![0.0_f64; means.len()];
        for _ in 0..limits.maximum_coupling_iterations {
            poll(cancelled)?;
            *coupling_iterations = coupling_iterations.checked_add(1)
                .ok_or(SingletTripletError::Unrepresentable("pump coupling iteration counter"))?;
            let (channels, outgoing) = self.trace(&means, incident, cancelled)?;
            let mut trials = Vec::with_capacity(self.models.len());
            let mut consistent = true;
            for (index, &model) in self.models.iter().enumerate() {
                let trial = cell::evaluate(
                    model, initial[index], &channels[index], duration,
                    means[index], limits.column, evaluations, cancelled,
                )?;
                consistent &= cell::consistent(
                    &trial, means[index], limits.column.kinetic.relative_tolerance,
                );
                trials.push(trial);
            }
            if consistent {
                // Every channel used exactly this same common field; its
                // independently propagated photon packet matches the one
                // shared kinetic history, within the original strict audit.
                return Ok(Some(Substep { cells: trials, outgoing }));
            }
            drop(trials);
            for (index, &model) in self.models.iter().enumerate() {
                let trial = cell::solve(
                    model, initial[index], &channels[index], duration,
                    limits.column, evaluations, cancelled,
                )?;
                let delta = trial.mean_ground - means[index];
                // Dampen an alternating source-iteration component without
                // clipping any population or altering the accepted history.
                let alternating = delta != 0.0 && previous_delta[index] != 0.0
                    && delta.is_sign_positive() != previous_delta[index].is_sign_positive();
                means[index] = if alternating {
                    means[index] + 0.5 * delta
                } else {
                    trial.mean_ground
                };
                previous_delta[index] = delta;
            }
        }
        // The caller discards this temporal grid and tries a finer one. Work
        // already spent still counts against the aggregate evaluation budget.
        Ok(None)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn uniform(
        &self,
        initial: &[[f64; 3]],
        incident: &[[f64; B]],
        duration: f64,
        substeps: usize,
        nonlinear: bool,
        limits: PumpBeamsLimits,
        evaluations: &mut usize,
        coupling_iterations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<PumpBeamsStep<B>>, SingletTripletError> {
        self.uniform_attributed(
            initial, incident, duration, substeps, nonlinear, limits,
            evaluations, coupling_iterations, None, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn uniform_attributed(
        &self,
        initial: &[[f64; 3]],
        incident: &[[f64; B]],
        duration: f64,
        substeps: usize,
        nonlinear: bool,
        limits: PumpBeamsLimits,
        evaluations: &mut usize,
        coupling_iterations: &mut usize,
        mut attribution: Option<&mut SourceAttribution>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<PumpBeamsStep<B>>, SingletTripletError> {
        let mut populations = initial.to_vec();
        let mut totals: Vec<_> = self.models.iter().map(|_| KineticTotals::new()).collect();
        let mut passive = vec![Sum::default(); self.models.len()];
        let mut beams: Vec<_> = self.paths.iter().map(|_| BeamTotals::<B>::new()).collect();
        let mut start = 0.0;
        for ordinal in 0..substeps {
            poll(cancelled)?;
            let end = duration * ((ordinal + 1) as f64 / substeps as f64);
            let h = end - start;
            if duration > 0.0 && h <= 0.0 {
                return Err(SingletTripletError::Unrepresentable("pump beam temporal resolution"));
            }
            let Some(trial) = self.substep(
                &populations, incident, h, nonlinear, limits,
                evaluations, coupling_iterations, cancelled,
            )? else {
                return Ok(None);
            };
            for (index, trial) in trial.cells.into_iter().enumerate() {
                poll(cancelled)?;
                populations[index] = trial.kinetics.populations;
                totals[index].add(&trial.kinetics)?;
                passive[index].add(trial.passive_heat_j)?;
                for (ordinal, visit) in self.visits[index].iter().enumerate() {
                    for lane in 0..B {
                        let packet = trial.channels[ordinal * B + lane];
                        beams[visit.beam].add_packet(lane, packet)?;
                        if let Some(split) = attribution.as_deref_mut() {
                            let energy = self.paths[visit.beam].photon_energies_j[lane];
                            split.add_packet(index, visit.beam, packet, energy,
                                self.models[index].excitation_energies_j())?;
                        }
                    }
                }
            }
            for (beam, rates) in beams.iter_mut().zip(trial.outgoing) {
                beam.add_transmitted(rates, h)?;
            }
            start = end;
        }
        let mut cells = Vec::with_capacity(self.models.len());
        let mut passive_heat_j = Vec::with_capacity(self.models.len());
        for (index, total) in totals.into_iter().enumerate() {
            poll(cancelled)?;
            cells.push(total.finish(
                self.models[index], initial[index], populations[index],
                limits.column.kinetic.relative_tolerance,
            )?);
            passive_heat_j.push(passive[index].value()?);
        }
        let mut summaries = Vec::with_capacity(beams.len());
        for (index, total) in beams.into_iter().enumerate() {
            poll(cancelled)?;
            summaries.push(total.finish(
                incident[index], duration, self.paths[index].photon_energies_j,
                limits.column.kinetic.relative_tolerance,
            )?);
        }
        let residual_j = self.energy_balance(
            initial, &cells, &passive_heat_j, &summaries,
            limits.column.kinetic.relative_tolerance,
        )?;
        Ok(Some(PumpBeamsStep {
            cells, passive_heat_j, beams: summaries, residual_j,
            substeps, cell_evaluations: *evaluations,
            coupling_iterations: *coupling_iterations,
            max_energy_discrepancy_j: 0.0,
        }))
    }
}
