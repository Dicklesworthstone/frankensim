//! Trial-local phase changes, native impulse moments and independent accounts.
use super::{
    CancelGate, JkrBranch, JkrNormalMotion, JkrNormalMotionError, JkrNormalMotionEvent,
    JkrNormalMotionState, JkrNormalMotionStep, JkrState, MAX_JKR_NORMAL_EVENTS, finite, poll,
    solve,
};

pub(super) fn run(
    law: &JkrNormalMotion<'_>,
    previous: JkrNormalMotionState,
    load: f64,
    dt: f64,
    gate: &CancelGate,
) -> Result<JkrNormalMotionStep, JkrNormalMotionError> {
    poll(gate)?;
    let before = law.response(previous)?;
    if !dt.is_finite()
        || dt <= 0.0
        || dt > law.budget.maximum_step_s
        || !load.is_finite()
        || load.abs() > law.budget.maximum_force_n
    {
        return Err(JkrNormalMotionError::Invalid(
            "finite step and held normal force inside budget",
        ));
    }
    let barrier =
        law.support
            .barrier(law.normal, load, law.budget.maximum_root_iterations, gate)?;
    let target = finite(previous.time_s + dt, "normal target time")?;
    if target <= previous.time_s {
        return Err(JkrNormalMotionError::Unresolved("normal clock precision"));
    }
    let steps = previous
        .steps
        .checked_add(1)
        .ok_or(JkrNormalMotionError::Range("normal step counter"))?;
    let mut state = previous;
    let mut elapsed = 0.0;
    let mut impulse = Sum::default();
    let mut moment = Sum::default();
    let mut work = Sum::default();
    let mut support_impulse = Sum::default();
    let mut support_moment = Sum::default();
    let mut support_work = Sum::default();
    let mut support_work_magnitude = Sum::default();
    let mut support_peak = 0.0_f64;
    let mut released = Sum::default();
    let mut travel = Sum::default();
    let mut events = [None; MAX_JKR_NORMAL_EVENTS];
    let (mut event_count, mut substeps, mut trials) = (0, 0, 0);
    while elapsed < dt {
        poll(gate)?;
        if state.contact.branch() == JkrBranch::Detached && state.contact.indentation_m() == 0.0 {
            apply_event(
                law,
                &mut state,
                elapsed,
                &mut events,
                &mut event_count,
                &mut released,
            )?;
        }
        let mut h = (dt - elapsed).min(law.budget.maximum_substep_s);
        let segment = loop {
            poll(gate)?;
            if trials == law.budget.maximum_substeps {
                return Err(JkrNormalMotionError::TrialBudget);
            }
            if h <= 0.0 || elapsed + h == elapsed {
                return Err(JkrNormalMotionError::Unresolved(
                    "normal substep time precision",
                ));
            }
            trials += 1;
            match solve::next(law, state, load, h, barrier, gate) {
                Ok(segment) => break segment,
                Err(JkrNormalMotionError::Cancelled) => {
                    return Err(JkrNormalMotionError::Cancelled);
                }
                Err(error @ JkrNormalMotionError::Range(_))
                    if state.contact.branch() == JkrBranch::Detached
                        && law.support.stiffness == 0.0 =>
                {
                    return Err(error);
                }
                Err(_) => h *= 0.5,
            }
        };
        let h = segment.duration;
        let end_time = finite(elapsed + h, "normal elapsed time")?;
        if end_time <= elapsed || end_time > dt {
            return Err(JkrNormalMotionError::Unresolved(
                "normal event time precision",
            ));
        }
        let increment = finite(
            segment.indentation - state.contact.indentation_m(),
            "normal increment",
        )?;
        let contact_impulse = finite(h * segment.force, "normal impulse")?;
        impulse.add(contact_impulse)?;
        moment.add(finite(
            contact_impulse * (elapsed + 0.5 * h),
            "normal impulse time moment",
        )?)?;
        work.add(finite(segment.force * increment, "normal smooth work")?)?;
        let elastic_impulse = finite(h * segment.support_force, "support impulse")?;
        support_impulse.add(elastic_impulse)?;
        support_moment.add(finite(
            elastic_impulse * (elapsed + 0.5 * h),
            "support impulse time moment",
        )?)?;
        let elastic_work = finite(segment.support_force * increment, "support work")?;
        support_work.add(elastic_work)?;
        support_work_magnitude.add(elastic_work.abs())?;
        support_peak = support_peak.max(segment.support_peak);
        travel.add(segment.travel)?;
        if travel.value > law.budget.maximum_accumulated_travel_m {
            return Err(JkrNormalMotionError::Range("accumulated normal travel"));
        }
        state.contact = match state.contact.branch() {
            JkrBranch::Attached => law.normal.evaluate_attached(segment.indentation)?.state(),
            JkrBranch::Detached => law.normal.detached_state(segment.indentation)?,
        };
        state.velocity_m_s = segment.velocity;
        elapsed = end_time;
        substeps += 1;
        if segment.event {
            apply_event(
                law,
                &mut state,
                elapsed,
                &mut events,
                &mut event_count,
                &mut released,
            )?;
        }
    }
    state.time_s = target;
    state.steps = steps;
    let contact = law.response(state)?;
    let increment = finite(
        state.contact.indentation_m() - previous.contact.indentation_m(),
        "normal total increment",
    )?;
    let applied = finite(load * increment, "applied normal work")?;
    let kinetic_before = finite(
        0.5 * law.mass_kg * previous.velocity_m_s * previous.velocity_m_s,
        "initial normal kinetic energy",
    )?;
    let kinetic = finite(
        0.5 * law.mass_kg * state.velocity_m_s * state.velocity_m_s,
        "normal kinetic energy",
    )?;
    let kinetic_change = finite(
        0.5 * law.mass_kg
            * (state.velocity_m_s - previous.velocity_m_s)
            * (state.velocity_m_s + previous.velocity_m_s),
        "normal kinetic change",
    )?;
    let storage_change = finite(
        contact.reversible_energy_j - before.reversible_energy_j,
        "normal storage change",
    )?;
    let support_before = law.support.energy(previous.contact.indentation_m())?;
    let support_energy = law.support.energy(state.contact.indentation_m())?;
    // The endpoint spring energy difference is evaluated in factored form,
    // independently of all substep work, so large preload does not erase it.
    let support_change = law.support.work(
        previous.contact.indentation_m(),
        state.contact.indentation_m(),
    )?;
    let kinetic_residual = finite(
        kinetic_change + work.value + support_work.value - applied,
        "normal kinetic work residual",
    )?;
    let storage_residual = finite(
        storage_change + released.value - work.value,
        "normal storage work residual",
    )?;
    let support_residual = finite(
        support_change - support_work.value,
        "support storage work residual",
    )?;
    let energy_residual = finite(
        kinetic_residual + storage_residual + support_residual,
        "normal energy residual",
    )?;
    let energy_scale = finite(
        kinetic_before
            + kinetic
            + before.reversible_energy_j.abs()
            + contact.reversible_energy_j.abs()
            + support_before
            + support_energy
            + support_work.value.abs()
            + work.value.abs()
            + applied.abs()
            + released.value,
        "normal energy scale",
    )?;
    // Channel-sized gates prevent a large, arbitrary spring preload datum from
    // masking lost contact work or kinetic work. Support uses factored endpoint
    // work and absolute path work, including cycles with zero net change.
    let contact_scale = before.reversible_energy_j.abs()
        + contact.reversible_energy_j.abs()
        + work.value.abs()
        + released.value;
    let support_scale = support_change.abs() + support_work_magnitude.value;
    let kinetic_scale =
        kinetic_before + kinetic + work.value.abs() + support_work_magnitude.value + applied.abs();
    for (residual, scale) in [
        (kinetic_residual, kinetic_scale),
        (storage_residual, contact_scale),
        (support_residual, support_scale),
        (energy_residual, energy_scale),
    ] {
        let energy_tolerance = finite(
            law.budget.energy_absolute_tolerance_j
                + law.budget.energy_relative_tolerance * finite(scale, "channel energy scale")?,
            "normal energy tolerance",
        )?;
        if residual.abs() > energy_tolerance {
            return Err(JkrNormalMotionError::EnergyBalance {
                residual_j: residual,
                tolerance_j: energy_tolerance,
            });
        }
    }
    let supplied_impulse = finite(dt * load, "applied normal impulse")?;
    let momentum_change = finite(
        law.mass_kg * (state.velocity_m_s - previous.velocity_m_s),
        "normal momentum change",
    )?;
    let reaction_impulse = finite(
        impulse.value + support_impulse.value,
        "total reaction impulse",
    )?;
    let reaction_moment = finite(
        moment.value + support_moment.value,
        "total reaction time moment",
    )?;
    let impulse_residual = finite(
        momentum_change - supplied_impulse + reaction_impulse,
        "normal impulse residual",
    )?;
    let impulse_scale = finite(
        law.mass_kg * (state.velocity_m_s.abs() + previous.velocity_m_s.abs())
            + supplied_impulse.abs()
            + impulse.value.abs()
            + support_impulse.value.abs(),
        "normal impulse scale",
    )?;
    let impulse_tolerance = finite(
        law.budget.impulse_absolute_tolerance_n_s
            + law.budget.impulse_relative_tolerance * impulse_scale,
        "normal impulse tolerance",
    )?;
    if impulse_residual.abs() > impulse_tolerance {
        return Err(JkrNormalMotionError::ImpulseBalance {
            residual_n_s: impulse_residual,
            tolerance_n_s: impulse_tolerance,
        });
    }
    // The receiving free bodies use exactly these two original force moments.
    // This independent integral checks that local root residuals did not make
    // the scalar contact coordinate drift away from their native trajectory.
    let restored_increment = finite(
        dt * previous.velocity_m_s
            + (0.5 * dt * supplied_impulse - (dt * reaction_impulse - reaction_moment))
                / law.mass_kg,
        "normal restored increment",
    )?;
    let coordinate_residual = finite(increment - restored_increment, "normal coordinate residual")?;
    let coordinate_tolerance = law.coordinate_tolerance(
        previous
            .contact
            .indentation_m()
            .abs()
            .max(state.contact.indentation_m().abs())
            .max(travel.value)
            .max(restored_increment.abs()),
    );
    if !coordinate_tolerance.is_finite() || coordinate_residual.abs() > coordinate_tolerance {
        return Err(JkrNormalMotionError::CoordinateBalance {
            residual_m: coordinate_residual,
            tolerance_m: coordinate_tolerance,
        });
    }
    let total_energy = finite(
        kinetic + contact.reversible_energy_j + support_energy,
        "normal total energy",
    )?;
    poll(gate)?;
    Ok(JkrNormalMotionStep {
        initial: previous,
        state,
        contact,
        normal_increment_m: increment,
        accumulated_travel_m: travel.value,
        normal_impulse_n_s: impulse.value,
        normal_time_moment_n_s2: moment.value,
        normal_work_j: work.value,
        support_impulse_n_s: support_impulse.value,
        support_time_moment_n_s2: support_moment.value,
        support_work_j: support_work.value,
        support_energy_j: support_energy,
        support_force_bound_n: support_peak,
        applied_work_j: applied,
        released_energy_j: released.value,
        kinetic_energy_j: kinetic,
        total_energy_j: total_energy,
        energy_residual_j: energy_residual,
        impulse_residual_n_s: impulse_residual,
        coordinate_residual_m: coordinate_residual,
        events,
        event_count,
        substeps,
        trials,
    })
}

fn apply_event(
    law: &JkrNormalMotion<'_>,
    state: &mut JkrNormalMotionState,
    time: f64,
    events: &mut [Option<JkrNormalMotionEvent>; MAX_JKR_NORMAL_EVENTS],
    count: &mut usize,
    released: &mut Sum,
) -> Result<(), JkrNormalMotionError> {
    if *count == law.budget.maximum_events {
        return Err(JkrNormalMotionError::EventBudget);
    }
    let before: JkrState = state.contact;
    let transition = law
        .normal
        .advance_monotone(before, before.indentation_m())?;
    let jump = transition.jump.ok_or(JkrNormalMotionError::Unresolved(
        "missing source branch event",
    ))?;
    state.contact = transition.response.state();
    law.response(*state)?;
    events[*count] = Some(JkrNormalMotionEvent {
        offset_s: time,
        jump,
        velocity_m_s: state.velocity_m_s,
        before,
        after: state.contact,
    });
    *count += 1;
    released.add(jump.released_energy_j)?;
    Ok(())
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), JkrNormalMotionError> {
        let adjusted = finite(value - self.correction, "normal accumulation correction")?;
        let next = finite(self.value + adjusted, "normal accumulation")?;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
        Ok(())
    }
}
