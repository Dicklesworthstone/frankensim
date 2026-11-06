//! Authored movement changes on the original integer event timeline.
//!
//! The model-free transport timeline owns validation and storage of immutable
//! rate frames. This binding admits their original model, compartments, path
//! order and construction limits. Rate changes redraw only changed movement
//! hazards at their actual authored time; nuclear and Poisson clocks survive.

use std::sync::Arc;

use crate::burnup::transport::schedule::NuclideTransportTimelineError;
pub use crate::burnup::transport::schedule::{NuclideTransportChange, NuclideTransportTimeline};

use super::{
    Accumulator, DecayError, DecayModel, NuclideTransport, StochasticDecayTransportBudget,
    StochasticDecayTransportLimits, StochasticDecayTransportState, StochasticDecayTransportStep,
    charge_clocks, sample_clock,
};

impl From<NuclideTransportTimelineError> for DecayError {
    fn from(error: NuclideTransportTimelineError) -> Self {
        match error {
            NuclideTransportTimelineError::Invalid(message) => Self::Invalid(message),
            NuclideTransportTimelineError::Budget(message) => Self::Budget(message),
        }
    }
}

/// Bounded immutable movement history for one original cohort shape.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayTransportSchedule {
    model: Arc<DecayModel>,
    compartments: usize,
    limits: StochasticDecayTransportLimits,
    timeline: NuclideTransportTimeline,
}

impl StochasticDecayTransportSchedule {
    pub fn new(
        state: &StochasticDecayTransportState,
        initial: NuclideTransport,
        changes: Vec<NuclideTransportChange>,
        maximum_changes: usize,
    ) -> Result<Self, DecayError> {
        Self::new_with_checkpoint(state, initial, changes, maximum_changes, || Ok(()))
    }

    /// Admit physical dimensions once, then reuse the shared cancellable frame
    /// validation. Neither the supplied state nor any existing schedule changes.
    pub fn new_with_checkpoint(
        state: &StochasticDecayTransportState,
        initial: NuclideTransport,
        changes: Vec<NuclideTransportChange>,
        maximum_changes: usize,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<Self, DecayError> {
        checkpoint()?;
        admit_initial(state, &initial, &mut checkpoint)?;
        let timeline = NuclideTransportTimeline::new_with_checkpoint(
            initial,
            changes,
            maximum_changes,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(Self {
            model: Arc::clone(state.model()),
            compartments: state.inventories.len(),
            limits: state.limits,
            timeline,
        })
    }

    #[must_use]
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }

    #[must_use]
    pub fn timeline(&self) -> &NuclideTransportTimeline {
        &self.timeline
    }

    #[must_use]
    pub fn changes(&self) -> &[NuclideTransportChange] {
        self.timeline.changes()
    }

    pub fn transport_at(&self, time_s: f64) -> Result<&NuclideTransport, DecayError> {
        self.timeline.transport_at(time_s).map_err(Into::into)
    }

    #[must_use]
    pub fn next_after(&self, time_s: f64) -> Option<f64> {
        self.timeline.next_after(time_s)
    }

    fn validate_binding(&self, state: &StochasticDecayTransportState) -> Result<(), DecayError> {
        if !Arc::ptr_eq(&self.model, state.model())
            || self.compartments != state.inventories.len()
            || self.limits != state.limits
            || !self.timeline.initial().same_topology(&state.transport)
        {
            return Err(DecayError::Invalid(
                "movement schedule belongs to different nuclides, compartments, paths or limits",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct TransportScheduleCursor {
    pub(super) schedule: Arc<StochasticDecayTransportSchedule>,
    position: usize,
}

impl TransportScheduleCursor {
    pub(super) fn active(&self) -> &NuclideTransport {
        self.schedule
            .timeline
            .at_position(self.position)
            .expect("retained movement position was admitted")
    }
}

impl StochasticDecayTransportState {
    /// Bind or replace movement history at the current absolute time.
    ///
    /// Earlier changes are not replayed. Only rates differing from the currently
    /// active frame redraw clocks, using the original movement stream. Nuclear
    /// clocks, local energy and movement counters are retained. A first binding
    /// occupies one additional clock slot; replacing a schedule reuses that slot.
    pub fn with_transport_schedule(
        mut self,
        schedule: Arc<StochasticDecayTransportSchedule>,
    ) -> Result<Self, DecayError> {
        schedule.validate_binding(&self)?;
        let clock_count = self
            .clock_count
            .checked_add(usize::from(self.transport_schedule.is_none()))
            .ok_or(DecayError::Budget("stochastic transport schedule clock"))?;
        if clock_count > self.limits.maximum_clocks {
            return Err(DecayError::Budget("stochastic transport schedule clock"));
        }
        let position = schedule.timeline.position_at(self.time_s())?;
        let active = schedule
            .timeline
            .at_position(position)
            .ok_or(DecayError::Arithmetic)?;
        self.replace_transport_clocks(active, Accumulator::new(self.time_s()), |_| Ok(()))?;
        self.transport_schedule = Some(TransportScheduleCursor { schedule, position });
        self.clock_count = clock_count;
        Ok(self)
    }

    #[must_use]
    pub fn transport_schedule(&self) -> Option<&Arc<StochasticDecayTransportSchedule>> {
        self.transport_schedule
            .as_ref()
            .map(|bound| &bound.schedule)
    }

    pub(super) fn next_transport_change(&self) -> Option<Accumulator> {
        let bound = self.transport_schedule.as_ref()?;
        bound
            .schedule
            .timeline
            .changes()
            .get(bound.position)
            .map(|change| Accumulator::new(change.time_s))
    }

    pub(super) fn apply_transport_change(
        &mut self,
        time: Accumulator,
        report: &mut StochasticDecayTransportStep,
        budget: StochasticDecayTransportBudget,
        checkpoint: &mut impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(), DecayError> {
        let bound = self
            .transport_schedule
            .as_ref()
            .ok_or(DecayError::Arithmetic)?;
        let schedule = Arc::clone(&bound.schedule);
        let position = bound
            .position
            .checked_add(1)
            .ok_or(DecayError::Arithmetic)?;
        let active = schedule
            .timeline
            .at_position(position)
            .ok_or(DecayError::Arithmetic)?;
        let mut visited = 0_usize;
        self.replace_transport_clocks(active, time, |changed| {
            if visited % 256 == 0 {
                checkpoint()?;
            }
            visited = visited.checked_add(1).ok_or(DecayError::Arithmetic)?;
            charge_clocks(
                report,
                1 + usize::from(changed),
                budget.maximum_clock_checks,
            )
        })?;
        self.transport_schedule
            .as_mut()
            .ok_or(DecayError::Arithmetic)?
            .position = position;
        report.transport_changes = report
            .transport_changes
            .checked_add(1)
            .ok_or(DecayError::Arithmetic)?;
        Ok(())
    }

    /// The old frame remains active throughout comparison; the immutable new
    /// frame is borrowed directly, with no rate copy or per-transition allocation.
    fn replace_transport_clocks(
        &mut self,
        frame: &NuclideTransport,
        time: Accumulator,
        mut visit: impl FnMut(bool) -> Result<(), DecayError>,
    ) -> Result<(), DecayError> {
        let transfers = self.transport.transfers.len();
        let paths = frame
            .transfers
            .iter()
            .map(|path| (path.from_compartment, &path.rates_per_s))
            .chain(
                frame
                    .outflows
                    .iter()
                    .map(|path| (path.compartment, &path.rates_per_s)),
            );
        for (path, (from, rates)) in paths.enumerate() {
            for (node, &rate) in rates.iter().enumerate() {
                let old_rate = if path < transfers {
                    self.transport().transfers[path].rates_per_s[node]
                } else {
                    self.transport().outflows[path - transfers].rates_per_s[node]
                };
                let changed = old_rate != rate;
                visit(changed)?;
                if changed {
                    self.movement_clocks[path * self.nodes + node] = sample_clock(
                        &mut self.movement_rng,
                        rate,
                        self.inventories[from].populations[node],
                        time,
                    )?;
                }
            }
        }
        Ok(())
    }
}

fn admit_initial(
    state: &StochasticDecayTransportState,
    initial: &NuclideTransport,
    checkpoint: &mut impl FnMut() -> Result<(), DecayError>,
) -> Result<(), DecayError> {
    let original = &state.transport;
    if initial.transfers.len() != original.transfers.len()
        || initial.outflows.len() != original.outflows.len()
    {
        return Err(DecayError::Invalid(
            "movement schedule must retain original paths",
        ));
    }
    for (index, (path, original)) in initial
        .transfers
        .iter()
        .zip(&original.transfers)
        .enumerate()
    {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if path.from_compartment != original.from_compartment
            || path.to_compartment != original.to_compartment
            || path.rates_per_s.len() != state.nodes
        {
            return Err(DecayError::Invalid(
                "movement schedule must retain paths and nuclide shape",
            ));
        }
    }
    for (index, (path, original)) in initial.outflows.iter().zip(&original.outflows).enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if path.compartment != original.compartment || path.rates_per_s.len() != state.nodes {
            return Err(DecayError::Invalid(
                "movement schedule must retain paths and nuclide shape",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
