//! Integer radionuclides moving between well-mixed compartments.
//!
//! Each authored coefficient k [1/s] gives a movement hazard k*N in its finite
//! donor. Movement and the original nuclear clocks share one chronological
//! first-reaction race. A transferred nucleus can decay, produce a daughter,
//! move again or escape during the same observation interval. Emission remains
//! in the compartment where the original decay event occurred; movement itself
//! contributes no heat, radiation or carrier enthalpy.
//!
//! All compartments retain the exact same immutable decay model and current
//! time, but distinct random streams. Movement has its own explicit stream.
//! Only hazards whose donor population changed are resampled, at the actual
//! compensated event time. Observation endpoints never redraw clocks. Clone
//! retains the complete future trajectory, including movement clocks and counts.
//! Exact clock ties retain compartment/node order before authored movement-path
//! order. The shared finite-grid clock and branch limitations still apply.

use std::{collections::BTreeSet, sync::Arc};

pub use crate::burnup::transport::{NuclideOutflow, NuclideTransfer, NuclideTransport};

use super::production::PoissonSource;
use super::{
    Accumulator, DecayError, DecayModel, EventClock, MAX_STOCHASTIC_NUCLEI, StochasticDecayBudget,
    StochasticDecayInventory, StochasticDecayStep, StreamCheckpoint, StreamKey, before, ceiling,
    checked, defect_bound, export_channels, increment, nonnegative, sample_clock,
};
use crate::{EnergyChannels, channel_values};

pub mod production;

/// Construction bounds for the complete cohort, including inactive clocks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StochasticDecayTransportLimits {
    pub maximum_compartments: usize,
    /// Every authored transfer/outflow path times the original nuclide count.
    pub maximum_rate_coefficients: usize,
    /// Per-cell decay clocks plus every per-path, per-nuclide movement clock.
    /// The production extension also charges bound arrival and schedule slots.
    pub maximum_clocks: usize,
}

impl Default for StochasticDecayTransportLimits {
    fn default() -> Self {
        Self {
            maximum_compartments: 64,
            maximum_rate_coefficients: 65_536,
            maximum_clocks: 69_632,
        }
    }
}

/// Separate whole-cohort bounds on nuclear events, movement and clock work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StochasticDecayTransportBudget {
    /// Event count and emitted-energy roundoff allowance cover ALL compartments.
    pub decay: StochasticDecayBudget,
    /// Transfers plus exterior removals in this trial; at most 2^53.
    pub maximum_movements: u64,
    /// Clock slots scanned or rescheduled, including inactive/zero-rate slots.
    /// Original branch selection is additionally bounded by the admitted model.
    pub maximum_clock_checks: u64,
}

impl StochasticDecayTransportBudget {
    fn validate(self) -> Result<(), DecayError> {
        self.decay.validate()?;
        if self.maximum_movements > MAX_STOCHASTIC_NUCLEI {
            return Err(DecayError::Budget("stochastic movement-count ceiling"));
        }
        Ok(())
    }
}

/// Interval-local events, with nuclear emissions attributed to their actual cell.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayTransportStep {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub duration_s: f64,
    pub decay: Vec<StochasticDecayStep>,
    pub decay_events: u64,
    pub movement_events: u64,
    /// Authored path order, then original nuclide order. Circulation is counted.
    pub transfer_counts: Vec<Vec<u64>>,
    pub outflow_counts: Vec<Vec<u64>>,
    /// Outward-rounded sum of the original per-cell emission-account defects.
    pub energy_roundoff_j: f64,
    pub movement_rng_draws: u64,
    pub clock_checks: u64,
}

/// An independent physical checkpoint of the original inventories and movement.
///
/// Rates remain constant for this owner. There is no external supply, inferred
/// flow, bulk-mass evolution or transported sensible/latent heat. Construction
/// bounds the complete retained population by 2^53; outflow can only reduce it.
/// Trial cloning and reports allocate once per call. Events reuse that storage.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayTransportState {
    inventories: Vec<StochasticDecayInventory>,
    transport: Arc<NuclideTransport>,
    outgoing_paths: Arc<Vec<Vec<usize>>>,
    movement_clocks: Vec<EventClock>,
    movement_rng: StreamCheckpoint,
    transfers: Vec<Vec<u64>>,
    outflows: Vec<Vec<u64>>,
    nodes: usize,
    clock_count: usize,
    limits: StochasticDecayTransportLimits,
}

#[derive(Clone, Copy)]
enum Event {
    Change { compartment: usize },
    Decay { compartment: usize, node: usize },
    Movement { path: usize, node: usize },
    Birth { compartment: usize, node: usize },
}

impl StochasticDecayTransportState {
    /// Bind constant authored hazards without changing existing decay history.
    /// All model Arcs and current times must match. Every cell stream key and
    /// the movement key must differ, so competing clocks do not reuse a stream.
    pub fn new(
        inventories: Vec<StochasticDecayInventory>,
        transport: NuclideTransport,
        movement_key: StreamKey,
        limits: StochasticDecayTransportLimits,
    ) -> Result<Self, DecayError> {
        if inventories.is_empty() || inventories.len() > limits.maximum_compartments {
            return Err(DecayError::Budget("stochastic transport compartments"));
        }
        let model = inventories[0].model();
        let nodes = model.nuclides().len();
        let time_s = inventories[0].time_s();
        if nodes == 0 {
            return Err(DecayError::Invalid(
                "transport requires a nonempty nuclide model",
            ));
        }
        let paths = transport
            .transfers
            .len()
            .checked_add(transport.outflows.len())
            .ok_or(DecayError::Budget("stochastic transport paths"))?;
        let coefficients = paths
            .checked_mul(nodes)
            .ok_or(DecayError::Budget("stochastic transport coefficients"))?;
        let clock_count = inventories
            .len()
            .checked_mul(nodes)
            .and_then(|decay| decay.checked_add(coefficients))
            .ok_or(DecayError::Budget("stochastic transport clocks"))?;
        if coefficients > limits.maximum_rate_coefficients || clock_count > limits.maximum_clocks {
            return Err(DecayError::Budget(
                "stochastic transport coefficients or clocks",
            ));
        }
        let mut keys = BTreeSet::new();
        keys.insert((movement_key.seed, movement_key.kernel, movement_key.tile));
        for inventory in &inventories {
            if !Arc::ptr_eq(inventory.model(), model) || inventory.time_s() != time_s {
                return Err(DecayError::Invalid(
                    "transport requires the same original model and time",
                ));
            }
            let key = inventory.rng_checkpoint().key;
            if !keys.insert((key.seed, key.kernel, key.tile)) {
                return Err(DecayError::Invalid(
                    "transport and cell random stream keys must differ",
                ));
            }
        }
        population_total(&inventories)?;
        let mut outgoing_paths = vec![Vec::new(); inventories.len()];
        for (path, transfer) in transport.transfers.iter().enumerate() {
            if transfer.from_compartment >= inventories.len()
                || transfer.to_compartment >= inventories.len()
                || transfer.from_compartment == transfer.to_compartment
            {
                return Err(DecayError::Invalid(
                    "invalid directed stochastic nuclide transfer",
                ));
            }
            admit_rates(&transfer.rates_per_s, nodes)?;
            outgoing_paths[transfer.from_compartment].push(path);
        }
        for (path, outflow) in transport.outflows.iter().enumerate() {
            if outflow.compartment >= inventories.len() {
                return Err(DecayError::Invalid("invalid stochastic nuclide outflow"));
            }
            admit_rates(&outflow.rates_per_s, nodes)?;
            outgoing_paths[outflow.compartment].push(transport.transfers.len() + path);
        }
        let mut result = Self {
            inventories,
            outgoing_paths: Arc::new(outgoing_paths),
            movement_clocks: vec![EventClock::Inactive; coefficients],
            movement_rng: StreamCheckpoint::current(movement_key, 0),
            transfers: vec![vec![0; nodes]; transport.transfers.len()],
            outflows: vec![vec![0; nodes]; transport.outflows.len()],
            transport: Arc::new(transport),
            nodes,
            clock_count,
            limits,
        };
        for path in 0..paths {
            for node in 0..nodes {
                result.reschedule_movement(path, node, Accumulator::new(time_s))?;
            }
        }
        Ok(result)
    }

    pub fn time_s(&self) -> f64 {
        self.inventories[0].time_s()
    }
    pub fn model(&self) -> &Arc<DecayModel> {
        self.inventories[0].model()
    }
    pub fn inventories(&self) -> &[StochasticDecayInventory] {
        &self.inventories
    }
    pub fn transport(&self) -> &NuclideTransport {
        &self.transport
    }
    pub fn limits(&self) -> StochasticDecayTransportLimits {
        self.limits
    }
    pub fn movement_rng_checkpoint(&self) -> StreamCheckpoint {
        self.movement_rng
    }
    pub fn cumulative_transfer_counts(&self) -> &[Vec<u64>] {
        &self.transfers
    }
    pub fn cumulative_outflow_counts(&self) -> &[Vec<u64>] {
        &self.outflows
    }

    /// First representable endpoint containing the next actual cohort event.
    /// Sub-ULP events retain their chronological order within that endpoint.
    pub fn next_event_time_s(&self) -> Option<f64> {
        self.next_event().map(|(_, time)| ceiling(time))
    }

    /// Evaluate an independent candidate. Any refusal preserves all inventories,
    /// pending movement/decay clocks, random streams and cumulative accounts.
    pub fn trial_to(
        &self,
        end_time_s: f64,
        budget: StochasticDecayTransportBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, StochasticDecayTransportStep), DecayError> {
        self.trial_with_sources(end_time_s, budget, None, checkpoint)
    }

    /// One race and accounting path for closed or externally supplied cohorts.
    /// A supply caller owns cloned candidate clocks; errors publish neither side.
    fn trial_with_sources(
        &self,
        end_time_s: f64,
        budget: StochasticDecayTransportBudget,
        mut production: Option<&mut production::ProductionTrial<'_>>,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, StochasticDecayTransportStep), DecayError> {
        checkpoint()?;
        budget.validate()?;
        let start_time_s = self.time_s();
        if !nonnegative(end_time_s) || end_time_s < start_time_s {
            return Err(DecayError::Invalid(
                "transport endpoint must be finite and nondecreasing",
            ));
        }
        let initial_total = population_total(&self.inventories)?;
        let mut current_total = initial_total;
        let clock_count = self
            .clock_count
            .checked_add(
                production
                    .as_ref()
                    .map_or(0, |source| source.source_clock_count),
            )
            .ok_or(DecayError::Budget("stochastic production transport clocks"))?;
        let mut next = self.clone();
        let mut report = StochasticDecayTransportStep {
            start_time_s,
            end_time_s,
            duration_s: end_time_s - start_time_s,
            decay: self
                .inventories
                .iter()
                .map(|inventory| StochasticDecayStep {
                    start_time_s,
                    end_time_s,
                    duration_s: end_time_s - start_time_s,
                    events: 0,
                    decays: vec![0; self.nodes],
                    branch_decays: vec![0; inventory.model.branches.len()],
                    released: EnergyChannels::default(),
                    energy_roundoff_j: 0.0,
                    rng_draws: 0,
                })
                .collect(),
            decay_events: 0,
            movement_events: 0,
            transfer_counts: vec![vec![0; self.nodes]; self.transfers.len()],
            outflow_counts: vec![vec![0; self.nodes]; self.outflows.len()],
            energy_roundoff_j: 0.0,
            movement_rng_draws: 0,
            clock_checks: 0,
        };
        let mut released = vec![[Accumulator::new(0.0); 3]; self.inventories.len()];
        let mut energy_errors = vec![Accumulator::new(0.0); self.inventories.len()];
        if end_time_s > start_time_s {
            loop {
                checkpoint()?;
                charge_clocks(&mut report, clock_count, budget.maximum_clock_checks)?;
                let sources = production
                    .as_deref()
                    .map_or(&[][..], |source| &*source.sources);
                let Some((event, time)) = next.next_event_with_sources(sources) else {
                    break;
                };
                if ceiling(time) > end_time_s {
                    break;
                }
                match event {
                    Event::Change { compartment } => {
                        let context = production.as_deref_mut().ok_or(DecayError::Arithmetic)?;
                        if context.schedule_changes == context.maximum_schedule_changes {
                            return Err(DecayError::Budget(
                                "stochastic production schedule changes",
                            ));
                        }
                        // Every source slot is visited, including unchanged rates.
                        charge_clocks(&mut report, self.nodes, budget.maximum_clock_checks)?;
                        context.sources[compartment]
                            .as_mut()
                            .ok_or(DecayError::Arithmetic)?
                            .change_rates(&mut next.inventories[compartment].rng, time)?;
                        context.schedule_changes = context
                            .schedule_changes
                            .checked_add(1)
                            .ok_or(DecayError::Arithmetic)?;
                    }
                    Event::Decay { compartment, node } => {
                        if report.decay_events == budget.decay.maximum_events {
                            return Err(DecayError::Budget("stochastic cohort decay events"));
                        }
                        // Original apply_decay reschedules its parent and daughter.
                        charge_clocks(&mut report, 2, budget.maximum_clock_checks)?;
                        let daughter = next.inventories[compartment].apply_decay(
                            node,
                            time,
                            &mut report.decay[compartment],
                            &mut released[compartment],
                            &mut energy_errors[compartment],
                        )?;
                        increment(&mut report.decay_events)?;
                        next.reschedule_outgoing(compartment, node, time, &mut report, budget)?;
                        if daughter != node {
                            next.reschedule_outgoing(
                                compartment,
                                daughter,
                                time,
                                &mut report,
                                budget,
                            )?;
                        }
                    }
                    Event::Movement { path, node } => {
                        if report.movement_events == budget.maximum_movements {
                            return Err(DecayError::Budget("stochastic cohort movements"));
                        }
                        next.apply_movement(path, node, time, &mut report, budget)?;
                        if path >= next.transport.transfers.len() {
                            current_total =
                                current_total.checked_sub(1).ok_or(DecayError::Arithmetic)?;
                        }
                        increment(&mut report.movement_events)?;
                    }
                    Event::Birth { compartment, node } => {
                        let context = production.as_deref_mut().ok_or(DecayError::Arithmetic)?;
                        if context.births == context.maximum_births {
                            return Err(DecayError::Budget("stochastic production births"));
                        }
                        // The shared birth mutation reschedules decay then supply.
                        charge_clocks(&mut report, 2, budget.maximum_clock_checks)?;
                        context.sources[compartment]
                            .as_mut()
                            .ok_or(DecayError::Arithmetic)?
                            .apply_birth(
                                &mut next.inventories[compartment],
                                node,
                                time,
                                &mut current_total,
                                &mut context.produced_nuclei[compartment],
                                &mut context.births,
                            )?;
                        next.reschedule_outgoing(compartment, node, time, &mut report, budget)?;
                    }
                }
            }
        }
        for (cell, inventory) in next.inventories.iter_mut().enumerate() {
            checkpoint()?;
            let local = &mut report.decay[cell];
            local.released = export_channels(released[cell], &mut energy_errors[cell])?;
            checked(channel_values(inventory.released).total_j())?;
            local.energy_roundoff_j = defect_bound(energy_errors[cell])?;
            // Outward addition prevents a small compartment's account defect
            // disappearing against a larger one at the whole-cohort boundary.
            if local.energy_roundoff_j > 0.0 {
                report.energy_roundoff_j = checked(
                    checked(report.energy_roundoff_j + local.energy_roundoff_j)?.next_up(),
                )?;
            }
            local.rng_draws = inventory
                .rng
                .index
                .checked_sub(self.inventories[cell].rng.index)
                .ok_or(DecayError::Arithmetic)?;
            inventory.time_s = end_time_s;
        }
        if report.energy_roundoff_j > budget.decay.energy_tolerance_j {
            return Err(DecayError::Budget(
                "stochastic cohort emission energy roundoff",
            ));
        }
        let removed = report
            .outflow_counts
            .iter()
            .flatten()
            .try_fold(0_u64, |sum, &count| {
                sum.checked_add(count).ok_or(DecayError::Arithmetic)
            })?;
        let supplied = production.as_ref().map_or(0, |source| source.births);
        let final_total = population_total(&next.inventories)?;
        if initial_total
            .checked_add(supplied)
            .and_then(|total| total.checked_sub(removed))
            != Some(final_total)
            || current_total != final_total
        {
            return Err(DecayError::Arithmetic);
        }
        report.movement_rng_draws = next
            .movement_rng
            .index
            .checked_sub(self.movement_rng.index)
            .ok_or(DecayError::Arithmetic)?;
        checkpoint()?;
        Ok((next, report))
    }

    pub fn advance_to(
        &mut self,
        end_time_s: f64,
        budget: StochasticDecayTransportBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<StochasticDecayTransportStep, DecayError> {
        let (next, report) = self.trial_to(end_time_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    fn next_event(&self) -> Option<(Event, Accumulator)> {
        self.next_event_with_sources(&[])
    }

    fn next_event_with_sources(
        &self,
        sources: &[Option<PoissonSource>],
    ) -> Option<(Event, Accumulator)> {
        let mut next: Option<(Event, Accumulator)> = None;
        // All changes precede a tied nuclear/movement event or source arrival,
        // including a birth in an earlier compartment at the same actual time.
        for (compartment, source) in sources.iter().enumerate() {
            if let Some(time) = source.as_ref().and_then(PoissonSource::next_change_time)
                && next.is_none_or(|(_, earliest)| before(time, earliest))
            {
                next = Some((Event::Change { compartment }, time));
            }
        }
        for (compartment, inventory) in self.inventories.iter().enumerate() {
            if let Some((node, time)) = inventory.next_event()
                && next.is_none_or(|(_, earliest)| before(time, earliest))
            {
                next = Some((Event::Decay { compartment, node }, time));
            }
        }
        for (index, clock) in self.movement_clocks.iter().enumerate() {
            if let EventClock::Finite(time) = clock
                && next.is_none_or(|(_, earliest)| before(*time, earliest))
            {
                next = Some((
                    Event::Movement {
                        path: index / self.nodes,
                        node: index % self.nodes,
                    },
                    *time,
                ));
            }
        }
        for (compartment, source) in sources.iter().enumerate() {
            if let Some((node, time)) = source.as_ref().and_then(PoissonSource::next_birth)
                && next.is_none_or(|(_, earliest)| before(time, earliest))
            {
                next = Some((Event::Birth { compartment, node }, time));
            }
        }
        next
    }

    fn apply_movement(
        &mut self,
        path: usize,
        node: usize,
        time: Accumulator,
        report: &mut StochasticDecayTransportStep,
        budget: StochasticDecayTransportBudget,
    ) -> Result<(), DecayError> {
        let (from, to) = if path < self.transport.transfers.len() {
            let transfer = &self.transport.transfers[path];
            (transfer.from_compartment, Some(transfer.to_compartment))
        } else {
            (
                self.transport.outflows[path - self.transport.transfers.len()].compartment,
                None,
            )
        };
        self.inventories[from].populations[node] = self.inventories[from].populations[node]
            .checked_sub(1)
            .ok_or(DecayError::Arithmetic)?;
        if let Some(to) = to {
            increment(&mut self.inventories[to].populations[node])?;
            increment(&mut self.transfers[path][node])?;
            increment(&mut report.transfer_counts[path][node])?;
        } else {
            let path = path - self.transport.transfers.len();
            increment(&mut self.outflows[path][node])?;
            increment(&mut report.outflow_counts[path][node])?;
        }
        charge_clocks(
            report,
            1 + usize::from(to.is_some()),
            budget.maximum_clock_checks,
        )?;
        self.inventories[from].reschedule(node, time)?;
        self.reschedule_outgoing(from, node, time, report, budget)?;
        if let Some(to) = to {
            self.inventories[to].reschedule(node, time)?;
            self.reschedule_outgoing(to, node, time, report, budget)?;
        }
        Ok(())
    }

    fn reschedule_outgoing(
        &mut self,
        compartment: usize,
        node: usize,
        time: Accumulator,
        report: &mut StochasticDecayTransportStep,
        budget: StochasticDecayTransportBudget,
    ) -> Result<(), DecayError> {
        let count = self.outgoing_paths[compartment].len();
        charge_clocks(report, count, budget.maximum_clock_checks)?;
        for offset in 0..count {
            let path = self.outgoing_paths[compartment][offset];
            self.reschedule_movement(path, node, time)?;
        }
        Ok(())
    }

    fn reschedule_movement(
        &mut self,
        path: usize,
        node: usize,
        time: Accumulator,
    ) -> Result<(), DecayError> {
        let (from, rate) = if path < self.transport.transfers.len() {
            let transfer = &self.transport.transfers[path];
            (transfer.from_compartment, transfer.rates_per_s[node])
        } else {
            let outflow = &self.transport.outflows[path - self.transport.transfers.len()];
            (outflow.compartment, outflow.rates_per_s[node])
        };
        self.movement_clocks[path * self.nodes + node] = sample_clock(
            &mut self.movement_rng,
            rate,
            self.inventories[from].populations[node],
            time,
        )?;
        Ok(())
    }
}

fn admit_rates(rates: &[f64], nodes: usize) -> Result<(), DecayError> {
    if rates.len() != nodes || rates.iter().any(|&rate| !nonnegative(rate)) {
        return Err(DecayError::Invalid(
            "one finite nonnegative movement rate per nuclide required",
        ));
    }
    Ok(())
}

fn population_total(inventories: &[StochasticDecayInventory]) -> Result<u64, DecayError> {
    let mut total = 0_u64;
    for &count in inventories
        .iter()
        .flat_map(|inventory| inventory.populations())
    {
        total = total
            .checked_add(count)
            .ok_or(DecayError::Budget("stochastic cohort population"))?;
        if total > MAX_STOCHASTIC_NUCLEI {
            return Err(DecayError::Budget("stochastic cohort population"));
        }
    }
    Ok(total)
}

fn charge_clocks(
    report: &mut StochasticDecayTransportStep,
    count: usize,
    maximum: u64,
) -> Result<(), DecayError> {
    let count = u64::try_from(count).map_err(|_| DecayError::Budget("stochastic clock work"))?;
    report.clock_checks = report
        .clock_checks
        .checked_add(count)
        .ok_or(DecayError::Budget("stochastic clock work"))?;
    if report.clock_checks > maximum {
        return Err(DecayError::Budget("stochastic clock work"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
