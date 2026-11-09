//! Ordered transport and bounded physical source iteration.

use super::accounts::BeamTotals;
use super::*;

struct Substep<const B: usize> {
    cells: Vec<cell::SharedTrial>,
    outgoing: Vec<[f64; B]>,
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
