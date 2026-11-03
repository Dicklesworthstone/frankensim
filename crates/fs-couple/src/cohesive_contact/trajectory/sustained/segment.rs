//! Stop-aware finite segments through the original cohesive material response.
//!
//! A reflected trial is private numerical continuation. A shortened endpoint
//! becomes motion only after the same complete sustained-contact audit passes.
//! No body velocity, material history, force, or contact coordinate is projected.

use super::{
    AuditError, BodyLoad, NaturalMapScales, Prepared, ResponseError, StoppingControls,
    StoppingError, StoppingReport, SustainedCandidate, SustainedContact, SustainedControls,
    SustainedError, Work, add, admit, finite, mul, poll, positive, solve_actual_controlled,
    stopping,
};
use fs_exec::Cx;
use fs_solid::cohesive::dynamics::CohesiveMotion;

/// Per-segment limits include every unpublished stopping-time probe.
#[derive(Clone, Copy, Debug)]
pub struct SegmentControls {
    pub step: SustainedControls,
    pub stopping: StoppingControls,
    pub maximum_contact_solves: usize,
    pub maximum_response_evaluations: usize,
}

/// One complete admitted segment. It can end before the requested time at a
/// resolved relative sliding stop; the caller then continues from these bodies.
pub struct SegmentCandidate<'m, 'a> {
    pub step: SustainedCandidate<'m, 'a>,
    pub report: SegmentReport,
}

#[derive(Debug)]
pub struct SegmentReport {
    pub requested_duration_s: f64,
    pub accepted_duration_s: f64,
    /// Bracketed local continuation event, not exact-flow CCD or uniqueness.
    pub stopping: Option<StoppingReport>,
    /// Original contact IDs whose complete relative tangential speed stopped.
    /// Moving supports or another body can leave their world speeds nonzero.
    pub stopped_contact_ids: Vec<u64>,
    pub contact_solves: usize,
    pub response_evaluations: usize,
}

#[derive(Debug)]
pub enum SegmentError {
    Invalid(&'static str),
    Budget(&'static str),
    Step(SustainedError),
    Stopping(StoppingError<SustainedError>),
}
impl core::fmt::Display for SegmentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive stop-aware segment refused: {self:?}")
    }
}
impl std::error::Error for SegmentError {}
impl From<SustainedError> for SegmentError {
    fn from(error: SustainedError) -> Self {
        Self::Step(error)
    }
}
impl From<ResponseError> for SegmentError {
    fn from(error: ResponseError) -> Self {
        Self::Step(error.into())
    }
}

/// Advance through loaded contact up to the requested duration or an admitted
/// sliding stop. Loads are held forces; prescribed jump offsets name the
/// requested endpoint and are interpolated to every original shorter trial.
#[allow(clippy::too_many_arguments)]
pub fn segment_candidate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    maximum_duration_s: f64,
    controls: SegmentControls,
    cx: &Cx<'_>,
) -> Result<SegmentCandidate<'m, 'a>, SegmentError> {
    segment_candidate_controlled(
        motions,
        contacts,
        loads,
        scales,
        initial_forces_n,
        maximum_duration_s,
        controls,
        &mut || cx.checkpoint().is_err(),
    )
}

/// Predicate counterpart with the original body/material and audit semantics.
///
/// A full-duration reversal supplies an unpublished upper bracket. Each time
/// query re-solves all contact forces and native mechanics from the unchanged
/// accepted bodies. Safeguarded secant/bisection locates a zero of the selected
/// initial slip-direction velocity; complete tangential velocity must also
/// vanish within the declared tolerance. A transverse turn is not a stop.
///
/// The final original replay must pass all contacts' unchanged audit. If that
/// reveals another reversal, its earlier bracket is resolved before adoption.
/// Simultaneous stops are identified only from actual final velocities. No
/// claim of a unique or earliest exact continuous nonlinear-flow root is made.
#[allow(clippy::too_many_arguments)]
pub fn segment_candidate_controlled<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    maximum_duration_s: f64,
    mut controls: SegmentControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SegmentCandidate<'m, 'a>, SegmentError> {
    admit(
        motions,
        contacts,
        loads,
        scales,
        initial_forces_n,
        maximum_duration_s,
        controls.step,
        cancelled,
    )?;
    if controls.maximum_contact_solves == 0
        || controls.maximum_response_evaluations == 0
        || !positive(controls.stopping.minimum_step_s)
        || !positive(controls.stopping.time_tolerance_s)
        || !positive(controls.stopping.velocity_tolerance_m_s)
        || controls.stopping.velocity_tolerance_m_s > controls.step.audit.velocity_tolerance_m_s
        || controls.stopping.maximum_trials == 0
    {
        return Err(SegmentError::Invalid(
            "bounded stop controls no looser than the physical velocity audit required",
        ));
    }
    // This owner's initial/interpolated offset arrays and event records are
    // additional to the single native candidate already bounded by step.
    let mut extra = add(mul(motions.len(), 16)?, mul(contacts.len(), 16)?)?;
    for motion in motions {
        extra = add(extra, mul(motion.state().points().len(), 6)?)?;
    }
    controls.step.maximum_workspace_entries = controls
        .step
        .maximum_workspace_entries
        .checked_sub(extra)
        .ok_or(SegmentError::Budget("stopping offsets and event records"))?;
    let start_s = motions[0].time_s();
    let target_s = finite(start_s + maximum_duration_s)?;
    let requested_s = finite(target_s - start_s)?;
    if requested_s <= 0.0 {
        return Err(SegmentError::Invalid("unrepresentable requested duration"));
    }
    let mut initial_offsets = Vec::with_capacity(motions.len());
    for (body, motion) in motions.iter().enumerate() {
        poll(cancelled)?;
        initial_offsets.push(
            motion
                .prescribed_jumps_m_controlled(cancelled)
                .map_err(|error| ResponseError::Body { body, error })?,
        );
    }
    let mut work = Work::default();
    let mut prepared = at_duration(
        motions,
        contacts,
        loads,
        &initial_offsets,
        scales,
        initial_forces_n,
        requested_s,
        requested_s,
        controls,
        &mut work,
        cancelled,
    )?;
    let mut event = None;
    loop {
        poll(cancelled)?;
        match prepared.audit(controls.step.audit, cancelled) {
            Ok(audited) => {
                // An exact stop can coincide with the requested endpoint.
                // The original force solve already freshly replayed this time.
                if event.is_none() {
                    if let Some(contact) = prepared.inputs.iter().find(|contact| {
                        let before = contact.relative_velocity_before_m_s;
                        let after = contact.relative_velocity_after_m_s;
                        contact.force_n[0] * contact.friction_coefficient > 0.0
                            && before[1].hypot(before[2]) > controls.stopping.velocity_tolerance_m_s
                            && after[1] == 0.0
                            && after[2] == 0.0
                    }) {
                        let after = contact.relative_velocity_after_m_s;
                        event = Some(StoppingReport {
                            elapsed_s: prepared.elapsed_s,
                            bracket_s: [prepared.elapsed_s; 2],
                            velocity_m_s: [after[1], after[2]],
                            response_evaluations: 0,
                        });
                    }
                }
                let mut stopped_contact_ids = Vec::new();
                if event.is_some() {
                    for contact in &prepared.inputs {
                        poll(cancelled)?;
                        let before = contact.relative_velocity_before_m_s;
                        let after = contact.relative_velocity_after_m_s;
                        if contact.force_n[0] * contact.friction_coefficient > 0.0
                            && before[1].hypot(before[2]) > controls.stopping.velocity_tolerance_m_s
                            && after[1].hypot(after[2]) <= controls.stopping.velocity_tolerance_m_s
                        {
                            stopped_contact_ids.push(contact.id);
                        }
                    }
                    if stopped_contact_ids.is_empty() {
                        return Err(SegmentError::Invalid(
                            "stopping replay contains no complete relative-velocity stop",
                        ));
                    }
                }
                let accepted_duration_s = prepared.elapsed_s;
                poll(cancelled)?;
                return Ok(SegmentCandidate {
                    step: prepared.accept(audited, &work),
                    report: SegmentReport {
                        requested_duration_s: requested_s,
                        accepted_duration_s,
                        stopping: event,
                        stopped_contact_ids,
                        contact_solves: work.contact_solves,
                        response_evaluations: work.response_evaluations,
                    },
                });
            }
            Err(AuditError::UnresolvedReversal { contact_id }) => {
                let index = prepared
                    .inputs
                    .iter()
                    .position(|c| c.id == contact_id)
                    .ok_or(SegmentError::Invalid(
                        "reversing contact absent from original order",
                    ))?;
                let before = prepared.inputs[index].relative_velocity_before_m_s;
                let after = prepared.inputs[index].relative_velocity_after_m_s;
                let upper_s = prepared.elapsed_s;
                // The prior reflected state may not contribute to any work,
                // loss or material history. It only supplied the upper guard.
                drop(prepared);
                let mut last: Option<Prepared<'m, 'a>> = None;
                let mut response = |h: f64, poller: &mut dyn FnMut() -> bool| {
                    drop(last.take());
                    let trial = at_duration(
                        motions,
                        contacts,
                        loads,
                        &initial_offsets,
                        scales,
                        initial_forces_n,
                        requested_s,
                        h,
                        controls,
                        &mut work,
                        &mut || poller(),
                    )?;
                    let velocity = trial.inputs[index].relative_velocity_after_m_s;
                    let sampled = (trial.elapsed_s, [velocity[1], velocity[2]]);
                    last = Some(trial);
                    Ok::<_, SustainedError>(sampled)
                };
                let stopped = stopping::locate(
                    [before[1], before[2]],
                    upper_s,
                    [after[1], after[2]],
                    controls.stopping,
                    &mut response,
                    cancelled,
                )
                .map_err(SegmentError::Stopping)?;
                drop(response);
                prepared = last.ok_or(SegmentError::Invalid("missing final stopping replay"))?;
                let actual = prepared.inputs[index].relative_velocity_after_m_s;
                if prepared.elapsed_s.to_bits() != stopped.elapsed_s.to_bits()
                    || actual[1].to_bits() != stopped.velocity_m_s[0].to_bits()
                    || actual[2].to_bits() != stopped.velocity_m_s[1].to_bits()
                {
                    return Err(SegmentError::Invalid(
                        "final native time/velocity differs from stopping replay",
                    ));
                }
                event = Some(stopped);
                // All contacts and all work are checked at this actual native
                // endpoint. Another unresolved contact cannot be hidden by a
                // successful selected-contact scalar root.
            }
            Err(error) => return Err(SustainedError::Audit(error).into()),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn at_duration<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    initial_offsets: &[Vec<[f64; 3]>],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    requested_s: f64,
    h: f64,
    controls: SegmentControls,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Prepared<'m, 'a>, SustainedError> {
    poll(cancelled)?;
    if work.contact_solves >= controls.maximum_contact_solves {
        return Err(ResponseError::Budget("complete stopping contact solves").into());
    }
    let remaining = controls
        .maximum_response_evaluations
        .checked_sub(work.response_evaluations)
        .filter(|remaining| *remaining > 0)
        .ok_or(ResponseError::Budget("stopping mechanical responses"))?;
    if !positive(h) || h > requested_s {
        return Err(ResponseError::Invalid("stopping time outside requested interval").into());
    }
    let start_s = motions[0].time_s();
    let target_s = finite(start_s + h)?;
    let actual_h = finite(target_s - start_s)?;
    if actual_h <= 0.0 || actual_h > requested_s {
        return Err(ResponseError::Invalid("unrepresentable stopping query duration").into());
    }
    let fraction = finite(actual_h / requested_s)?;
    let mut offsets = Vec::with_capacity(loads.len());
    for (body, load) in loads.iter().enumerate() {
        poll(cancelled)?;
        let mut values = Vec::with_capacity(initial_offsets[body].len());
        for (before, after) in initial_offsets[body]
            .iter()
            .zip(load.endpoint_prescribed_jumps_m)
        {
            poll(cancelled)?;
            let mut value = [0.0; 3];
            for axis in 0..3 {
                value[axis] = if actual_h == requested_s {
                    after[axis]
                } else {
                    finite(before[axis] + fraction * finite(after[axis] - before[axis])?)?
                };
            }
            values.push(value);
        }
        offsets.push(values);
    }
    let shorter_loads: Vec<_> = loads
        .iter()
        .zip(&offsets)
        .map(|(load, offsets)| BodyLoad {
            mean_force_n: load.mean_force_n,
            endpoint_prescribed_jumps_m: offsets,
        })
        .collect();
    let mut step = controls.step;
    step.solve.maximum_response_evaluations =
        step.solve.maximum_response_evaluations.min(remaining);
    let prepared = solve_actual_controlled(
        motions,
        contacts,
        &shorter_loads,
        scales,
        initial_forces_n,
        actual_h,
        step,
        work,
        cancelled,
    )?;
    if prepared.target_time_s.to_bits() != target_s.to_bits() {
        return Err(ResponseError::Invalid("stopping endpoint clock changed").into());
    }
    Ok(prepared)
}

#[cfg(test)]
mod tests;
