//! Prescribed pump/inlet histories on the original graph thermal owners.

use super::*;

impl GraphEnthalpyTransport {
    /// Attach the native continuous flow-multiplier and inlet-enthalpy history.
    ///
    /// The original signed graph and fixed inventory remain authoritative.
    /// A common nonnegative multiplier changes all internal and external mass
    /// rates together, including pump-off intervals; it does not reverse flow,
    /// evolve composition, solve momentum or infer pump work. Inlet IDs and
    /// caloric domains are admitted by the original thermal schedule owner.
    ///
    /// Explicit steps sample at their accepted start times and clip at schedule
    /// knots. The existing dense, sparse and nonlinear backward-Euler entries
    /// sample at their common endpoint, once per interval. Split those single
    /// implicit intervals at every interior knot. Both schemes retain their
    /// original temporal order and need time refinement for varying forcing.
    pub fn with_schedule(
        mut self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        schedule: impl Into<Arc<AdvectionSchedule>>,
    ) -> Result<Self, TransientTransportError> {
        self.validate_storage(cx, thermal)?;
        let schedule = schedule.into();
        schedule.flow_at_time(&self.advection, thermal, thermal.time_s(), || {
            thermal_checkpoint(cx)
        })?;
        self.schedule = Some(schedule);
        Ok(self)
    }

    /// Shared original history, absent for the unchanged fixed-flow mode.
    pub fn schedule(&self) -> Option<&AdvectionSchedule> {
        self.schedule.as_deref()
    }

    /// Sample one original endpoint, before any Newton/Krylov evaluation. The
    /// returned binding has no schedule and therefore cannot resample it or
    /// recursively integrate it. The original public binding remains reusable.
    pub(super) fn freeze_scheduled_interval(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
    ) -> Result<Option<Self>, TransientTransportError> {
        let Some(schedule) = &self.schedule else {
            return Ok(None);
        };
        self.validate_storage(cx, thermal)?;
        if !target_s.is_finite() || target_s <= thermal.time_s() {
            return Err(TransientTransportError::Invalid(
                "scheduled implicit transport needs a future finite endpoint",
            ));
        }
        if schedule
            .next_transition_time_s(thermal.time_s())
            .is_some_and(|knot| knot < target_s)
        {
            return Err(TransientTransportError::Invalid(
                "split the native implicit transport interval at its next advection schedule knot",
            ));
        }
        let advection = schedule.flow_at_time(&self.advection, thermal, target_s, || {
            thermal_checkpoint(cx)
        })?;
        Ok(Some(Self {
            node_cells: self.node_cells.clone(),
            storage: self.storage.clone(),
            advection,
            max_connections: self.max_connections,
            reconstruction: self.reconstruction.clone(),
            schedule: None,
        }))
    }
}
