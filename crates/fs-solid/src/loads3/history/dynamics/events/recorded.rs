//! Exact acceleration event sampling into the original control recorder.
use super::*;
use crate::loads3::transient::adjoint::TransientAdjointError3;
use crate::loads3::transient::adjoint::support::trajectory::RecordedPrescribedTransient3;

fn adjoint_error(error: TransientHistoryError3) -> TransientAdjointError3 {
    match error {
        TransientHistoryError3::History(error) => error.into(),
        TransientHistoryError3::Dynamics(error) => error.into(),
    }
}

impl EventLoadHistory3<'_> {
    /// Advance a recorded acceleration history to its next exact event,
    /// requested target or step cap, using the original left/right sampler.
    ///
    /// Initialize from this schedule with `begin_transient_with_support_accelerations`,
    /// transfer `into_run`, then `record_prescribed_controls`. The supplied
    /// recorder must retain that same reference problem, material and initial
    /// source state. A different schedule is not an implicit source transfer.
    ///
    /// A finite acceleration/load event records two states at the same physical
    /// time: the accepted LEFT interval endpoint followed by its RIGHT mass
    /// equilibrium. Observations and force/control rows retain that order.
    /// Initial acceleration derivatives are in `gradient.initial.a`; later
    /// endpoint and right-limit derivatives are in `gradient.accelerations`.
    /// The caller chains these world-component cotangents to authored values.
    /// Frame parameters, knot times and adaptive decisions remain fixed.
    ///
    /// Sampling, the interval, the right-limit solve and recording publish as
    /// one transaction. A late failure retains the complete original motion and
    /// tape. Only the primal state is forked, never the recorded trajectory.
    pub fn advance_recorded_accelerations_to(
        &self,
        recorded: &mut RecordedPrescribedTransient3<'_, '_, '_, '_>,
        target_s: f64,
        maximum_step_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ScheduledStep3, TransientAdjointError3> {
        recorded.admit_next_step(cx)?;
        if !matches!(self.supports, SupportSchedule3::Acceleration(_)) {
            return Err(TransientAdjointError3::Invalid(
                "recorded acceleration events require an acceleration schedule",
            ));
        }
        let run = recorded.run();
        let start_time_s = run.state().time_s;
        let end_time_s = self
            .endpoint(start_time_s, target_s, maximum_step_s)
            .map_err(adjoint_error)?;
        let left = self
            .left_at(
                run.geometry,
                run.problem,
                end_time_s,
                run.maximum_displacement_m(),
                cx,
            )
            .map_err(adjoint_error)?;
        let right = self
            .supports
            .at(
                &self.history,
                run.geometry,
                run.problem,
                end_time_s,
                run.maximum_displacement_m(),
                cx,
            )
            .map_err(adjoint_error)?;
        let h_s = end_time_s - start_time_s;
        let (interval, jump) = if left.values() != right.values()
            || left.loads().nodal_forces_n() != right.loads().nodal_forces_n()
        {
            let event = recorded.advance_to_support_acceleration_jump(
                &left.loads().case,
                left.values(),
                &right.loads().case,
                right.values(),
                h_s,
                cx,
            )?;
            (event.interval, Some(event.jump))
        } else {
            (
                recorded.advance_with_support_accelerations(
                    &left.loads().case,
                    left.values(),
                    h_s,
                    cx,
                )?,
                None,
            )
        };
        Ok(ScheduledStep3 {
            start_time_s,
            end_time_s,
            interval,
            jump,
        })
    }
}
