//! Resolved stops and turns through the original cohesive material response.
//!
//! A reflected trial is private numerical continuation. A shortened endpoint
//! becomes motion only after the same complete sustained-contact audit passes.
//! No body velocity, material history, force, or contact coordinate is projected.

use super::{
    AuditError, BodyLoad, NaturalMapScales, Prepared, ResponseError, StoppingControls,
    StoppingError, StoppingReport, SustainedCandidate, SustainedContact, SustainedControls,
    SustainedError, TurningControls, TurningError, TurningReport, Work, add, admit, finite, mul,
    poll, positive, solve_actual_controlled, stopping, turning,
};
use fs_exec::Cx;
use fs_solid::cohesive::dynamics::CohesiveMotion;

/// Per-segment limits include every unpublished stopping and turning probe.
#[derive(Clone, Copy, Debug)]
pub struct SegmentControls {
    pub step: SustainedControls,
    pub stopping: StoppingControls,
    pub turning: TurningControls,
    pub maximum_contact_solves: usize,
    pub maximum_response_evaluations: usize,
}

/// One complete admitted segment. It can end before the requested time at a
/// resolved relative sliding stop or angular subdivision; the caller then
/// continues from these bodies with a newly solved force.
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
    /// Actual angular subdivision, distinct from a complete relative stop.
    pub turning: Option<TurningReport>,
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
    Turning(TurningError),
}
impl core::fmt::Display for SegmentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive resolved contact segment refused: {self:?}")
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
impl From<TurningError> for SegmentError {
    fn from(error: TurningError) -> Self {
        Self::Turning(error)
    }
}

/// Advance through loaded contact up to the requested duration or an admitted
/// sliding stop or resolved turn. Loads are held forces; prescribed jump offsets name the
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
/// Resolved turns exceeding the configured angular bound instead subdivide the
/// actual duration and recompute every force and body from the same predecessor.
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
    turning::validate(controls.turning)?;
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
    let mut turn_trigger = None;
    let mut turn_refinements = 0;
    loop {
        poll(cancelled)?;
        let angular = angular_change(&prepared, controls, cancelled)?;
        let audit = prepared.audit(controls.step.audit, cancelled);
        let subdivide = match &audit {
            Ok(_) => angular.exceeds,
            Err(AuditError::UnresolvedReversal { .. }) => angular.transverse_exceeds,
            _ => None,
        };
        if let Some(contact_id) = subdivide {
            let previous_s = prepared.elapsed_s;
            let h = turning::shorter_duration(previous_s, controls.turning, turn_refinements)?;
            let actual_s = finite(finite(start_s + h)? - start_s)?;
            if actual_s < controls.turning.minimum_step_s || actual_s >= previous_s {
                return Err(TurningError::Unrepresentable {
                    previous_s,
                    requested_s: h,
                    actual_s,
                }
                .into());
            }
            turn_trigger.get_or_insert(contact_id);
            turn_refinements = add(turn_refinements, 1)?;
            // Neither the preceding reflected candidate nor a prior stopping
            // receipt belongs to this new original-duration solve.
            drop(prepared);
            event = None;
            prepared = at_duration(
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
                cancelled,
            )?;
            continue;
        }
        match audit {
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
                        turning: turn_trigger.map(|trigger_contact_id| TurningReport {
                            trigger_contact_id,
                            elapsed_s: accepted_duration_s,
                            maximum_angle_rad: angular.maximum_angle_rad,
                            refinements: turn_refinements,
                        }),
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
                    )
                    .map_err(|error| StoppingProbeError::Failure(error.into()))?;
                    let velocity = trial.inputs[index].relative_velocity_after_m_s;
                    let sampled = (trial.elapsed_s, [velocity[1], velocity[2]]);
                    let transverse_turn = angular_change(&trial, controls, &mut || poller())
                        .map_err(StoppingProbeError::Failure)?
                        .transverse_exceeds
                        .is_some();
                    last = Some(trial);
                    if transverse_turn {
                        Err(StoppingProbeError::Turning)
                    } else {
                        Ok(sampled)
                    }
                };
                let result = stopping::locate(
                    [before[1], before[2]],
                    upper_s,
                    [after[1], after[2]],
                    controls.stopping,
                    &mut response,
                    cancelled,
                );
                drop(response);
                let stopped = match result {
                    Ok(stopped) => stopped,
                    Err(StoppingError::Response(StoppingProbeError::Turning)) => {
                        prepared = last.ok_or(SegmentError::Invalid(
                            "missing original transverse-turn response",
                        ))?;
                        event = None;
                        // The complete audit and angular bounds select the next
                        // action; this intermediate scalar probe is not a stop.
                        continue;
                    }
                    Err(error) => return Err(stopping_refusal(error)),
                };
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

struct AngularChange {
    maximum_angle_rad: f64,
    exceeds: Option<u64>,
    transverse_exceeds: Option<u64>,
}

/// Every relevant contact contributes to the admitted angular bound. A truly
/// collinear reflected endpoint still belongs to the original stopping search.
fn angular_change(
    prepared: &Prepared<'_, '_>,
    controls: SegmentControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularChange, SegmentError> {
    let mut result = AngularChange {
        maximum_angle_rad: 0.0,
        exceeds: None,
        transverse_exceeds: None,
    };
    for contact in &prepared.inputs {
        poll(cancelled)?;
        if finite(contact.force_n[0] * contact.friction_coefficient)? <= 0.0 {
            continue;
        }
        let before = contact.relative_velocity_before_m_s;
        let after = contact.relative_velocity_after_m_s;
        let Some(change) = turning::measure(
            [before[1], before[2]],
            [after[1], after[2]],
            controls.stopping.velocity_tolerance_m_s,
        )?
        else {
            continue;
        };
        result.maximum_angle_rad = result.maximum_angle_rad.max(change.angle_rad);
        if change.angle_rad > controls.turning.maximum_angle_rad {
            result.exceeds = Some(result.exceeds.map_or(contact.id, |id| id.min(contact.id)));
            if change.transverse_resolved {
                result.transverse_exceeds = Some(
                    result
                        .transverse_exceeds
                        .map_or(contact.id, |id| id.min(contact.id)),
                );
            }
        }
    }
    Ok(result)
}

enum StoppingProbeError {
    Failure(SegmentError),
    /// A complete native intermediate response resolved a sliding turn.
    Turning,
}

/// Keep ordinary stopping refusal details while the private turn marker is
/// intercepted by the surrounding segment loop, never published as mechanics.
fn stopping_refusal(error: StoppingError<StoppingProbeError>) -> SegmentError {
    let error = match error {
        StoppingError::Response(StoppingProbeError::Failure(error)) => return error,
        StoppingError::Response(StoppingProbeError::Turning) => {
            return SegmentError::Invalid("unhandled original turning endpoint");
        }
        StoppingError::Invalid(reason) => StoppingError::Invalid(reason),
        StoppingError::NonFinite(operation) => StoppingError::NonFinite(operation),
        StoppingError::Cancelled => StoppingError::Cancelled,
        StoppingError::ResponsePanicked => StoppingError::ResponsePanicked,
        StoppingError::NotSliding { initial_speed_m_s } => {
            StoppingError::NotSliding { initial_speed_m_s }
        }
        StoppingError::NotBracketed {
            upper_projection_m_s,
        } => StoppingError::NotBracketed {
            upper_projection_m_s,
        },
        StoppingError::TransverseVelocity {
            elapsed_s,
            velocity_m_s,
        } => StoppingError::TransverseVelocity {
            elapsed_s,
            velocity_m_s,
        },
        StoppingError::MinimumStep {
            bracket_s,
            minimum_step_s,
        } => StoppingError::MinimumStep {
            bracket_s,
            minimum_step_s,
        },
        StoppingError::Unrepresentable {
            bracket_s,
            requested_s,
            actual_s,
        } => StoppingError::Unrepresentable {
            bracket_s,
            requested_s,
            actual_s,
        },
        StoppingError::Budget {
            evaluations,
            maximum,
        } => StoppingError::Budget {
            evaluations,
            maximum,
        },
        StoppingError::Replay {
            expected_elapsed_s,
            actual_elapsed_s,
            velocity_m_s,
        } => StoppingError::Replay {
            expected_elapsed_s,
            actual_elapsed_s,
            velocity_m_s,
        },
    };
    SegmentError::Stopping(error)
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
