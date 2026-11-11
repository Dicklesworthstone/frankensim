//! Stationary optical histories without subtracting gross self-return traffic.
//!
//! A same-cell, same-state photon return has no net population transition.
//! Cancel that physical cycle before forming the rate residual. This permits
//! stiff steady histories whose gross event-count ulps exceed the net escaped
//! traffic. Admission checks every original state's integrated residual on
//! the original external/storage budget; the final photon and energy audit
//! still applies. This is numerical stationarity, not a certified equilibrium.

use super::super::accounts::BeamTotals;
use super::*;
use fs_math::decay::{attenuate, decayed_fraction};

impl<const B: usize> SaturablePhotonRecycling<B> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn stationary_history(
        &self,
        frozen: &SaturablePumpBeams<B>,
        initial: &[[f64; 3]],
        external: &[[f64; B]],
        local_pumps: &[LocalPhotonPump],
        duration: f64,
        limits: SaturableRecyclingLimits,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<History<B>>, SingletTripletError> {
        let count = initial.len();
        let mut flows = vec![[Sum::default(); 3]; count];
        let mut source_rates = vec![[0.0; 2]; count];
        let mut source_events = vec![[0.0; 2]; count];
        let mut partition = vec![[Sum::default(); 2]; count];
        let mut scale = Sum::default();
        let mut local_reports = Vec::with_capacity(local_pumps.len());
        for (index, model) in frozen.models.iter().enumerate() {
            poll(cancelled)?;
            let pump = local::at(local_pumps, index);
            let mut local_events = [0.0; 2];
            for state in 0..2 {
                let local_rate = positive_product(
                    initial[index][0],
                    pump.rates_s[state],
                    "stationary local photon rate",
                )?;
                transition(&mut flows[index], 0, state + 1, local_rate)?;
                local_events[state] =
                    positive_product(local_rate, duration, "stationary local photons")?;
                let n = initial[index][state + 1];
                scale.add(positive_product(
                    n,
                    model.excitation_energies_j()[state],
                    "stationary initial energy",
                )?)?;
                source_rates[index][state] = positive_product(
                    n,
                    model.radiative_rates_s()[state],
                    "stationary radiative rate",
                )?;
                source_events[index][state] = positive_product(
                    source_rates[index][state],
                    duration,
                    "stationary radiative events",
                )?;
                transition(
                    &mut flows[index],
                    state + 1,
                    0,
                    positive_product(
                        n,
                        model.nonradiative_rates_s()[state],
                        "stationary dark rate",
                    )?,
                )?;
                transition(
                    &mut flows[index],
                    state + 1,
                    2 - state,
                    positive_product(
                        n,
                        model.crossing_rates_s()[state],
                        "stationary crossing rate",
                    )?,
                )?;
                partition[index][state].add(1.0)?;
            }
            if !local_pumps.is_empty() {
                let local =
                    LocalPumpStep::from_events(local_events, pump, model.excitation_energies_j())?;
                for energy in local.absorbed_energy_j {
                    scale.add(energy)?;
                }
                local_reports.push(local);
            }
        }
        let mut incident = external.to_vec();
        let mut incident_counts = Vec::with_capacity(frozen.paths.len());
        for rates in external {
            let mut counts = [0.0; B];
            for lane in 0..B {
                counts[lane] =
                    positive_product(rates[lane], duration, "stationary external photons")?;
            }
            incident_counts.push(counts);
        }
        for emitted in &self.emission {
            let mut lanes = [0.0; B];
            let mut counts = [0.0; B];
            for lane in 0..B {
                let weight = emitted.photon_fractions[lane];
                lanes[lane] = positive_product(
                    source_rates[emitted.source_cell][emitted.source_state],
                    weight,
                    "stationary emitted lane",
                )?;
                counts[lane] = positive_product(
                    source_events[emitted.source_cell][emitted.source_state],
                    weight,
                    "stationary emitted photons",
                )?;
                partition[emitted.source_cell][emitted.source_state].add(-weight)?;
            }
            incident.push(lanes);
            incident_counts.push(counts);
        }
        let mut attribution = SourceAttribution::new(count, self.external_paths);
        let mut beams = Vec::with_capacity(incident.len());
        for (ordinal, (path, rates)) in frozen.paths.iter().zip(&incident).enumerate() {
            poll(cancelled)?;
            let source = if ordinal < self.external_paths {
                None
            } else {
                let emitted = &self.emission[ordinal - self.external_paths];
                Some((emitted.source_cell, emitted.source_state))
            };
            if source.is_none() {
                for lane in 0..B {
                    scale.add(product3(
                        rates[lane],
                        duration,
                        path.photon_energies_j[lane],
                    )?)?;
                }
            }
            let mut outgoing = *rates;
            let mut outgoing_counts = incident_counts[ordinal];
            let mut total = BeamTotals::new();
            for segment in &path.segments {
                poll(cancelled)?;
                let cell_index = segment.cell_index;
                let ground = initial[cell_index][0];
                for lane in 0..B {
                    let incoming = outgoing[lane];
                    let incoming_count = outgoing_counts[lane];
                    let depth = finite(
                        ground * segment.optical_depth_per_ground[lane].iter().sum::<f64>()
                            + segment.passive_optical_depth[lane],
                        "stationary optical depth",
                    )?;
                    let mean = if depth == 0.0 {
                        1.0
                    } else {
                        decayed_fraction(depth) / depth
                    };
                    let passive_rate =
                        product3(incoming, segment.passive_optical_depth[lane], mean)?;
                    let mut absorbed_by_state = [0.0; 2];
                    for target in 0..2 {
                        let active_depth = positive_product(
                            ground,
                            segment.optical_depth_per_ground[lane][target],
                            "stationary active depth",
                        )?;
                        let rate = product3(incoming, active_depth, mean)?;
                        absorbed_by_state[target] = product3(incoming_count, active_depth, mean)?;
                        match source {
                            Some((source_cell, source_state))
                                if source_cell == cell_index && source_state == target => {}
                            Some((source_cell, source_state)) => {
                                transition(&mut flows[source_cell], source_state + 1, 0, rate)?;
                                transition(&mut flows[cell_index], 0, target + 1, rate)?;
                            }
                            None => transition(&mut flows[cell_index], 0, target + 1, rate)?,
                        }
                    }
                    if let Some((cell, state)) = source {
                        transition(&mut flows[cell], state + 1, 0, passive_rate)?;
                    }
                    outgoing[lane] = attenuate(incoming, depth);
                    outgoing_counts[lane] = attenuate(incoming_count, depth);
                    let packet = cell::ChannelStep {
                        incident: incoming_count,
                        transmitted: outgoing_counts[lane],
                        absorbed: finite(
                            absorbed_by_state[0] + absorbed_by_state[1],
                            "stationary absorbed photons",
                        )?,
                        absorbed_by_state,
                        passive: product3(
                            incoming_count,
                            segment.passive_optical_depth[lane],
                            mean,
                        )?,
                        transmitted_rate: outgoing[lane],
                        residual: 0.0,
                    };
                    total.add_packet(lane, packet)?;
                    attribution.add_packet(
                        cell_index,
                        ordinal,
                        packet,
                        path.photon_energies_j[lane],
                        frozen.models[cell_index].excitation_energies_j(),
                    )?;
                }
            }
            if let Some((cell, state)) = source {
                for value in outgoing {
                    transition(&mut flows[cell], state + 1, 0, value)?;
                }
            }
            total.add_transmitted_photons(outgoing_counts)?;
            beams.push(total.finish_photons(
                incident_counts[ordinal],
                path.photon_energies_j,
                limits.pump.column.kinetic.relative_tolerance,
            )?);
        }
        let allowance = 0.015625 * limits.pump.column.kinetic.relative_tolerance * scale.value()?;
        for cell in 0..count {
            poll(cancelled)?;
            for state in 0..2 {
                // Retain the supplied fraction roundoff explicitly; do not
                // silently normalize a quadrature partition to unit weight.
                let unmatched = finite(
                    source_rates[cell][state] * partition[cell][state].value()?,
                    "stationary spectrum remainder",
                )?;
                transition(&mut flows[cell], state + 1, 0, unmatched)?;
            }
            for state in 0usize..3 {
                let energy = frozen.models[cell].excitation_energies_j()[state.saturating_sub(1)];
                let discrepancy = product3(flows[cell][state].value()?.abs(), duration, energy)?;
                if discrepancy > allowance {
                    return Ok(None);
                }
            }
        }
        let mut cells = Vec::with_capacity(count);
        let mut passive_heat_j = Vec::with_capacity(count);
        for (index, model) in frozen.models.iter().enumerate() {
            poll(cancelled)?;
            *evaluations =
                evaluations
                    .checked_add(1)
                    .ok_or(SingletTripletError::Unrepresentable(
                        "stationary cell evaluation counter",
                    ))?;
            cap(
                "pump column cell evaluations",
                *evaluations,
                limits.pump.column.maximum_cell_evaluations,
            )?;
            let mut cell = SaturableSingletTripletPlan::prepare(
                *model,
                [0.0; 2],
                [0.0; 2],
                0.0,
                limits.pump.column.kinetic,
                &mut *cancelled,
            )?
            .advance_with_cancel(initial[index], &mut *cancelled)?;
            for state in 0..2 {
                let mut events = attribution.external_events[index][state];
                events.add(attribution.recycled_events[index][state].value()?)?;
                if let Some(local) = local_reports.get(index) {
                    events.add(local.events[state])?;
                }
                cell.pump_events[state] = events.value()?;
                let mut absorbed = attribution.external_absorbed[index][state];
                absorbed.add(attribution.recycled_absorbed[index][state].value()?)?;
                if let Some(local) = local_reports.get(index) {
                    absorbed.add(local.absorbed_energy_j[state])?;
                }
                cell.absorbed_pump_energy_j[state] = absorbed.value()?;
                let mut heat = attribution.external_relaxation[index][state];
                heat.add(attribution.recycled_relaxation[index][state].value()?)?;
                if let Some(local) = local_reports.get(index) {
                    heat.add(local.relaxation_heat_j[state])?;
                }
                cell.pump_relaxation_heat_j[state] = heat.value()?;
                let residence =
                    positive_product(initial[index][state + 1], duration, "stationary residence")?;
                cell.radiative_events[state] = source_events[index][state];
                cell.nonradiative_events[state] = positive_product(
                    residence,
                    model.nonradiative_rates_s()[state],
                    "stationary dark events",
                )?;
                cell.crossing_events[state] = positive_product(
                    residence,
                    model.crossing_rates_s()[state],
                    "stationary crossing events",
                )?;
                cell.radiative_decay_energy_j[state] = positive_product(
                    cell.radiative_events[state],
                    model.excitation_energies_j()[state],
                    "stationary radiative energy",
                )?;
                cell.nonradiative_heat_j[state] = positive_product(
                    cell.nonradiative_events[state],
                    model.excitation_energies_j()[state],
                    "stationary dark heat",
                )?;
            }
            let gap = model.excitation_energies_j()[0] - model.excitation_energies_j()[1];
            cell.crossing_heat_j =
                positive_product(cell.crossing_events[0], gap, "stationary crossing heat")?;
            cell.reverse_crossing_heat_absorbed_j =
                positive_product(cell.crossing_events[1], gap, "stationary reverse heat")?;
            pairs::audit(
                *model,
                initial[index],
                &mut cell,
                limits.pump.column.kinetic.relative_tolerance,
            )?;
            cells.push(cell);
            let mut passive = attribution.external_passive[index];
            passive.add(attribution.recycled_passive[index].value()?)?;
            passive_heat_j.push(passive.value()?);
        }
        let beams = PumpBeamsStep {
            cells,
            passive_heat_j,
            beams,
            local_pumps: local_reports,
            residual_j: 0.0,
            substeps: 1,
            cell_evaluations: *evaluations,
            coupling_iterations: 0,
            max_energy_discrepancy_j: 0.0,
        };
        Ok(Some(History {
            beams,
            external_absorbed: attribution
                .external_absorbed
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            recycled_absorbed: attribution
                .recycled_absorbed
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            external_events: attribution
                .external_events
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            recycled_events: attribution
                .recycled_events
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            external_relaxation: attribution
                .external_relaxation
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            recycled_relaxation: attribution
                .recycled_relaxation
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            external_passive: attribution
                .external_passive
                .into_iter()
                .map(|v| v.value())
                .collect::<Result<_, _>>()?,
            recycled_passive: attribution
                .recycled_passive
                .into_iter()
                .map(|v| v.value())
                .collect::<Result<_, _>>()?,
            source_allowance_j: 0.0,
        }))
    }
}

fn transition(
    flow: &mut [Sum; 3],
    from: usize,
    to: usize,
    rate: f64,
) -> Result<(), SingletTripletError> {
    flow[from].add(-rate)?;
    flow[to].add(rate)?;
    Ok(())
}
