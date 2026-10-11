//! Exact schedule boundaries over the same transactional thermal owner.
use super::*;

impl<G: Clone> DepositionSimulation<'_, G> {
    /// Advance through the prefix of `events` at or before `end_time_s`.
    ///
    /// All supplied event headers must be ordered by time and strictly increasing
    /// identity, and none may precede the accepted state. Equal-time events are
    /// allowed. Only the returned `report.events` prefix is consumed; pass the
    /// remaining suffix on the next call. No future cell/contact is allocated or
    /// included in heat transfer before its event. A birth exactly at the horizon
    /// is included, but receives no pre-birth cooling.
    ///
    /// Every substep is shortened at the next event. The entire call, including
    /// already traversed events, is atomic on cancellation, budget, domain, or
    /// balance refusal. A clone before the call is an in-memory restart point.
    /// Resuming at the same event/substep boundaries preserves operation order;
    /// splitting a timestep at a new observation time changes the discretization.
    pub fn run_until(
        &mut self,
        cx: &Cx<'_>,
        events: &[DepositionEvent<G>],
        end_time_s: f64,
        config: DepositionStepConfig,
    ) -> Result<DepositionReport, DepositionError> {
        validate_config(config)?;
        poll(cx, 0)?;
        if !end_time_s.is_finite() || end_time_s < self.time_s {
            return Err(DepositionError::Invalid("time must be finite and nondecreasing"));
        }
        if events.len() > self.budget.max_events {
            return Err(DepositionError::Budget("schedule headers"));
        }
        let mut previous_time = self.time_s;
        let mut previous_id = self.last_event;
        for (index, event) in events.iter().enumerate() {
            poll(cx, index)?;
            if !event.time_s.is_finite() || event.time_s < previous_time
                || previous_id.is_some_and(|id| event.id <= id) {
                return Err(DepositionError::Invalid("schedule must have nondecreasing times and increasing IDs"));
            }
            previous_time = event.time_s;
            previous_id = Some(event.id);
        }
        let mut trial = self.clone();
        let before = trial.stored_enthalpy_j();
        let mut report = DepositionReport::default();
        for event in events.iter().take_while(|e| e.time_s <= end_time_s) {
            trial.march(cx, event.time_s, config, &mut report)?;
            trial.activate(cx, event, &mut report)?;
        }
        trial.march(cx, end_time_s, config, &mut report)?;
        finish(&trial, before, &mut report, config.energy_tolerance_j)?;
        poll(cx, report.steps)?;
        *self = trial;
        Ok(report)
    }
}
