//! Time resolution of original radial momentum and regional enthalpy.
//!
//! Compare one native interval with two native half intervals, retaining only
//! the fine native states. Every radial annulus and thermal cell participates:
//! equal volume flows do not imply equal shear histories. Separate bounds on
//! sign-changing interval duration and cancelled directional volume prevent a
//! nearly zero mean flow from hiding unresolved counterflow at uniform heat.
//! This follows the existing inertial transport controller's dyadic pattern;
//! it adds neither an extrapolated physical state nor another heat integrator.

use std::cell::Cell;

use super::equipment::EquipmentHeatDestinations;
use crate::graph::rheology::transient::head::FlowHeadSource;
use fs_exec::Cx;
use fs_qty::Pressure;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

use super::{
    CarreauTransportExchange, FixedPressure, FlowInjection, ImplicitTransportError,
    TransientCarreauGraph, TransientCarreauState, TransientCarreauStep,
    TransientCarreauTransportError, TransientCarreauWorkspace, TransientRegionalCarreauControls,
    TransientRegionalCarreauReport, TransientRegionalCarreauStep,
    TransientRegionalCarreauTransport, thermal_checkpoint,
};

/// Local discretization, counterflow resolution and total work limits.
#[derive(Clone, Copy, Debug)]
pub struct RegionalCarreauResolutionControls {
    /// Maximum full comparison interval, allowing a one-ULP endpoint snap.
    /// Retained native intervals are half its actual width or smaller.
    /// Author drive discontinuities explicitly; native
    /// thermal inventory events retain their separate admission.
    pub maximum_interval_s: f64,
    /// Refinement floor; a shorter event/final interval cannot refine again.
    pub minimum_interval_s: f64,
    /// Includes discarded coarse/fine attempts and the terminal refused step.
    pub maximum_native_steps: usize,
    /// Also constrained by the original thermal inventory's max_substeps over
    /// the complete requested window, including all authored drive segments.
    pub maximum_accepted_steps: usize,
    /// All original branch solves, including discarded comparisons and every
    /// thermal residual. Each native step is clipped to the remaining limit.
    pub maximum_total_branch_solves: u64,
    /// Every original equipment query across thermal trials and discarded
    /// comparisons. Zero allows histories with no flow-dependent head sources.
    pub maximum_total_head_evaluations: u64,
    pub absolute_cell_error_j: f64,
    /// Relative to interval energy CHANGES, independent of caloric datum.
    pub relative_cell_error: f64,
    pub absolute_velocity_error_m_s: f64,
    pub relative_velocity_error: f64,
    pub absolute_flow_error_m3_s: f64,
    pub relative_flow_error: f64,
    /// Separate per-branch viscous heat, held-head work, original equipment
    /// hydraulic/input work, losses and each routed heat channel. Compare both
    /// signed transfer and total absolute half-step transfer
    /// so opposed sources or work across time cannot hide by cancellation.
    pub absolute_transfer_error_j: f64,
    pub relative_transfer_error: f64,
    /// Longest retained substep with opposite endpoint flow signs.
    pub maximum_reversal_step_s: f64,
    /// Largest min(V+, V-) of a branch's native affine flow trajectory, per
    /// retained substep. This bounds cancellation, not the net displaced volume.
    pub maximum_cancelled_volume_m3: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RegionalCarreauResolutionWork {
    pub native_steps: usize,
    /// Every drive query, including discarded comparisons and terminal source
    /// refusals. Prescribed histories query once per native attempt; thermal
    /// feedback queries once per prepared original thermal endpoint evaluation.
    /// Native attempt and endpoint limits bound both modes.
    pub drive_evaluations: usize,
    pub rejected_intervals: usize,
    /// Exact counts from completed native solves, including discarded ones.
    /// On terminal native failure its original remaining budget still bounds
    /// attempted work, but that failed owner's partial count is not inferred.
    pub completed_branch_solves: u64,
    /// Same completed-native accounting and terminal remaining-budget bound
    /// as branch solves, for original equipment source queries.
    pub completed_head_evaluations: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct RegionalCarreauResolutionEstimate {
    pub start_s: f64,
    pub end_s: f64,
    pub maximum_cell_difference_j: f64,
    pub maximum_velocity_difference_m_s: f64,
    pub maximum_flow_difference_m3_s: f64,
    pub maximum_transfer_difference_j: f64,
    pub error_ratio: f64,
    pub maximum_reversal_step_s: f64,
    pub maximum_cancelled_volume_m3: f64,
}

/// Physical native pressure, flow and optional branch-head data at one time.
/// No interpolation, source model or additional thermal stream is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct RegionalCarreauDrive {
    pub boundaries: Vec<FixedPressure>,
    pub injections: Vec<FlowInjection>,
    pub branch_heads: Option<Vec<Pressure>>,
    /// Original signed-flow head sources, one optional law per native branch.
    /// They solve inside the original radial/pressure Jacobian; accepted
    /// reports retain these exact shared sources and actual head work.
    pub branch_head_laws: Option<Vec<Option<FlowHeadSource>>>,
    /// Required recipients for positive original pump/drive routed fractions.
    /// The owner deposits those losses; exchange callbacks supply other heat.
    pub branch_heat_destinations: Option<Vec<EquipmentHeatDestinations>>,
}

/// The exact forcing used by one accepted native midpoint momentum interval.
/// Prescribed histories reuse one sample on every thermal source/derivative
/// trial. Feedback instead retains ONLY the final accepted endpoint's sample;
/// each nonlinear trial resolves its own temperatures and physical drive.
/// Samples on discarded comparison intervals carry no physical credit.
#[derive(Clone, Debug, PartialEq)]
pub struct RegionalCarreauDriveSample {
    pub start_s: f64,
    pub end_s: f64,
    pub sample_time_s: f64,
    /// None for a prescribed drive. Some(end_s) means the drive uses this
    /// accepted thermal endpoint with the separately stated midpoint time.
    pub thermal_time_s: Option<f64>,
    pub drive: RegionalCarreauDrive,
}
impl std::borrow::Borrow<RegionalCarreauDrive> for RegionalCarreauDriveSample {
    fn borrow(&self) -> &RegionalCarreauDrive {
        &self.drive
    }
}

/// All accepted native substeps, in chronological order. Their original
/// pressure/head work, viscous heat, kinetic and caloric changes remain separate.
/// Only these reports contribute to the requested physical accounting window.
#[derive(Debug)]
pub struct RegionalCarreauTrajectory<T> {
    pub steps: Vec<TransientRegionalCarreauReport<T>>,
    /// One actual midpoint drive sample per retained native report, in identical
    /// order, with explicit thermal-endpoint time for source feedback.
    pub drive_samples: Vec<RegionalCarreauDriveSample>,
    pub estimates: Vec<RegionalCarreauResolutionEstimate>,
    pub work: RegionalCarreauResolutionWork,
}

pub struct ResolvedRegionalCarreau<'m, T> {
    pub next: TransientRegionalCarreauTransport,
    pub momentum: TransientCarreauState,
    pub thermal: DepositionSimulation<'m>,
    pub report: RegionalCarreauTrajectory<T>,
}

#[derive(Debug)]
pub enum RegionalCarreauResolutionError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    WorkLimit {
        reason: &'static str,
        work: RegionalCarreauResolutionWork,
    },
    Drive {
        start_s: f64,
        end_s: f64,
        sample_time_s: f64,
        work: RegionalCarreauResolutionWork,
        source: E,
    },
    Step {
        start_s: f64,
        end_s: f64,
        work: RegionalCarreauResolutionWork,
        source: ImplicitTransportError<TransientCarreauTransportError<E>>,
    },
    Underresolved {
        estimate: RegionalCarreauResolutionEstimate,
        work: RegionalCarreauResolutionWork,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for RegionalCarreauResolutionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved regional Carreau transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RegionalCarreauResolutionError<E> {}

/// Resolve an event-free window under held pressure, flow and optional head
/// ports. The original exchange callback is pure and observes each candidate's
/// actual interval and signed mean flow, so boundary streams follow reversal.
///
/// Comparisons restart the same accepted radial AND caloric state. Only the two
/// native fine steps may advance the local trajectory; coarse/rejected states
/// carry no physical credit. Every original source, mass, rate, nonlinear and
/// conservation refusal is terminal, even during a discarded comparison. The
/// whole requested window returns synchronized replacements or leaves both
/// inputs unchanged, including after a locally accepted prefix. Workspace is
/// reusable scratch and may change on failure.
///
/// The first-order full/half difference is a local estimate, not a rigorous
/// error bound. Reversal bounds concern the native affine flow trajectory, not
/// an exact zero-crossing locator or proof that unsampled oscillations cannot
/// occur. Refine the maximum interval for time-dependent exchange and retain
/// the original spatial refinement and developed-flow applicability limits.
#[allow(clippy::too_many_arguments)]
pub fn advance_resolved<'m, T, E>(
    cx: &Cx<'_>,
    owner: &TransientRegionalCarreauTransport,
    thermal: &DepositionSimulation<'m>,
    state: &TransientCarreauState,
    workspace: &mut TransientCarreauWorkspace,
    target_s: f64,
    boundaries: &[FixedPressure],
    injections: &[FlowInjection],
    branch_heads: Option<&[Pressure]>,
    step_controls: TransientRegionalCarreauControls,
    controls: RegionalCarreauResolutionControls,
    exchange: impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &TransientCarreauStep,
    ) -> Result<CarreauTransportExchange<T>, E>,
) -> Result<ResolvedRegionalCarreau<'m, T>, RegionalCarreauResolutionError<E>> {
    let drive = RegionalCarreauDrive {
        boundaries: boundaries.to_vec(),
        injections: injections.to_vec(),
        branch_heads: branch_heads.map(<[Pressure]>::to_vec),
        branch_head_laws: None,
        branch_heat_destinations: None,
    };
    advance_resolved_with_drive(
        cx,
        owner,
        thermal,
        state,
        workspace,
        target_s,
        &[],
        step_controls,
        controls,
        |_, _| Ok(drive.clone()),
        |cx, candidate, flow, _| exchange(cx, candidate, flow),
    )
}

/// Resolve a prescribed physical drive history with exactly the same native
/// adaptive momentum/enthalpy controller as [`advance_resolved`]. `drive` is
/// queried ONCE at each proposed native interval's strict interior midpoint,
/// including the discarded coarse step and both candidate half steps. Each
/// returned forcing is held throughout that interval's thermal nonlinear solve;
/// it is never sampled at the endpoint or independently on derivative probes.
/// The exchange callback receives this exact sample alongside the thermal
/// endpoint and original flow report to author its actual incoming streams.
///
/// `discontinuities_s` must be finite and strictly increasing. It may include
/// times outside this window. Every comparison window is clipped to the next
/// authored time, and both fine steps remain on that same side of the event.
/// The next comparison starts exactly there; no callback is queried at the
/// discontinuity itself. Its comparison-width proposal restarts at the authored
/// maximum, since the previous forcing's local estimate no longer applies.
/// An endpoint exactly one floating-point successor below the next event or
/// target snaps to that authored time; native steps and work use the actual
/// resulting spans. This avoids a nonrepresentable final comparison without
/// bypassing any physical or temporal acceptance check.
/// A pressure/head jump preserves the actual velocity
/// history. No waveform, smoothing, valve topology or impulsive state reset is
/// inferred. Native thermal inventory events still require their own admission.
///
/// Time-varying injections at pressure-reservoir nodes are admitted by the
/// original momentum owner. A free rigid junction instead requires its held
/// injection to match initial, midpoint AND final flow; an incompatible change
/// refuses, even at an authored event. This adapter does not invent the pressure
/// impulse or changing-mass owner needed to change that constraint.
///
/// Only accepted fine-drive samples and native reports are returned. Every
/// native/source refusal is terminal and the complete window remains atomic.
/// The original first-order thermal coupling, local error estimates, spatial
/// resolution and maximum-interval controls still apply to smooth waveforms;
/// midpoint sampling cannot detect arbitrarily fast unsampled oscillations.
#[allow(clippy::too_many_arguments)]
pub fn advance_resolved_with_drive<'m, T, E>(
    cx: &Cx<'_>,
    owner: &TransientRegionalCarreauTransport,
    thermal: &DepositionSimulation<'m>,
    state: &TransientCarreauState,
    workspace: &mut TransientCarreauWorkspace,
    target_s: f64,
    discontinuities_s: &[f64],
    step_controls: TransientRegionalCarreauControls,
    controls: RegionalCarreauResolutionControls,
    drive: impl Fn(&Cx<'_>, f64) -> Result<RegionalCarreauDrive, E>,
    exchange: impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &TransientCarreauStep,
        &RegionalCarreauDriveSample,
    ) -> Result<CarreauTransportExchange<T>, E>,
) -> Result<ResolvedRegionalCarreau<'m, T>, RegionalCarreauResolutionError<E>> {
    advance_with_driver(
        cx,
        owner,
        thermal,
        state,
        workspace,
        target_s,
        discontinuities_s,
        step_controls,
        controls,
        DriveResolver::Prescribed(&drive),
        exchange,
    )
}

/// Resolve authored pressure/head/flow feedback in the SAME original regional
/// thermal residual and adaptive/event controller. Each source callback receives
/// the strict interval midpoint, the candidate native thermal endpoint, and the
/// original radial graph re-prepared at those local endpoint temperatures.
/// This permits temperature-controlled pressure or head laws to feed back into
/// the accepted flow and heat. The callback is pure and must validate its own
/// source support; no pump curve, equation of state or control law is inferred.
///
/// Unlike the prescribed-drive entry, EVERY prepared thermal residual and
/// derivative probe queries its own drive. Each trial restarts the same accepted
/// radial history and holds that trial's forcing during native midpoint momentum.
/// The exchange callback sees exactly that trial sample. Only the final native
/// acceptance sample is retained; its `thermal_time_s` is Some(end_s), while
/// `sample_time_s` remains the time coordinate used by the authored drive law.
/// This endpoint-temperature/midpoint-momentum coupling remains first order.
///
/// Original work, viscous heat, transport and conservation gates are unchanged:
/// temperature-dependent pressure is mechanical input, not an added caloric
/// source. Fixed reference density, fixed regional masses and zero compliance
/// remain required. This supplies no thermal expansion/EOS, valve topology,
/// controller state or flow-dependent pump-curve iteration. Free-junction
/// injections retain the original initial/midpoint/final continuity restriction.
///
/// Native attempt, thermal endpoint and branch-work limits include rejected
/// candidates. Feedback refusal is a terminal Step/Resolve/Drive error, never a
/// reason to hide the source failure by shrinking time. The complete requested
/// window publishes synchronized original owners or leaves every input intact.
#[allow(clippy::too_many_arguments)]
pub fn advance_resolved_with_thermal_drive<'m, T, E>(
    cx: &Cx<'_>,
    owner: &TransientRegionalCarreauTransport,
    thermal: &DepositionSimulation<'m>,
    state: &TransientCarreauState,
    workspace: &mut TransientCarreauWorkspace,
    target_s: f64,
    discontinuities_s: &[f64],
    step_controls: TransientRegionalCarreauControls,
    controls: RegionalCarreauResolutionControls,
    drive: impl Fn(
        &Cx<'_>,
        f64,
        &DepositionSimulation<'m>,
        &TransientCarreauGraph,
    ) -> Result<RegionalCarreauDrive, E>,
    exchange: impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &TransientCarreauStep,
        &RegionalCarreauDriveSample,
    ) -> Result<CarreauTransportExchange<T>, E>,
) -> Result<ResolvedRegionalCarreau<'m, T>, RegionalCarreauResolutionError<E>> {
    advance_with_driver(
        cx,
        owner,
        thermal,
        state,
        workspace,
        target_s,
        discontinuities_s,
        step_controls,
        controls,
        DriveResolver::Thermal(&drive),
        exchange,
    )
}

type ThermalDrive<'a, 'm, E> = dyn Fn(
        &Cx<'_>,
        f64,
        &DepositionSimulation<'m>,
        &TransientCarreauGraph,
    ) -> Result<RegionalCarreauDrive, E>
    + 'a;

enum DriveResolver<'a, 'm, E> {
    Prescribed(&'a dyn Fn(&Cx<'_>, f64) -> Result<RegionalCarreauDrive, E>),
    Thermal(&'a ThermalDrive<'a, 'm, E>),
}

#[allow(clippy::too_many_arguments)]
fn advance_with_driver<'m, T, E>(
    cx: &Cx<'_>,
    owner: &TransientRegionalCarreauTransport,
    thermal: &DepositionSimulation<'m>,
    state: &TransientCarreauState,
    workspace: &mut TransientCarreauWorkspace,
    target_s: f64,
    discontinuities_s: &[f64],
    step_controls: TransientRegionalCarreauControls,
    controls: RegionalCarreauResolutionControls,
    drive: DriveResolver<'_, 'm, E>,
    exchange: impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &TransientCarreauStep,
        &RegionalCarreauDriveSample,
    ) -> Result<CarreauTransportExchange<T>, E>,
) -> Result<ResolvedRegionalCarreau<'m, T>, RegionalCarreauResolutionError<E>> {
    checkpoint(cx)?;
    let start = thermal.time_s();
    if !target_s.is_finite()
        || target_s <= start
        || !(target_s - start).is_finite()
        || start != state.accepted_time_s()
        || !positive(controls.maximum_interval_s)
        || !positive(controls.minimum_interval_s)
        || controls.minimum_interval_s > controls.maximum_interval_s
        || !(1..=1_000_000).contains(&controls.maximum_native_steps)
        || !(2..=4096).contains(&controls.maximum_accepted_steps)
        || controls.maximum_total_branch_solves == 0
        || !tolerances(controls.absolute_cell_error_j, controls.relative_cell_error)
        || !tolerances(
            controls.absolute_velocity_error_m_s,
            controls.relative_velocity_error,
        )
        || !tolerances(
            controls.absolute_flow_error_m3_s,
            controls.relative_flow_error,
        )
        || !tolerances(
            controls.absolute_transfer_error_j,
            controls.relative_transfer_error,
        )
        || !positive(controls.maximum_reversal_step_s)
        || !nonnegative(controls.maximum_cancelled_volume_m3)
    {
        return Err(RegionalCarreauResolutionError::Invalid(
            "finite synchronized interval and bounded native resolution controls required",
        ));
    }
    for (index, &time) in discontinuities_s.iter().enumerate() {
        if index % 64 == 0 {
            checkpoint(cx)?;
        }
        if !time.is_finite() || index > 0 && time <= discontinuities_s[index - 1] {
            return Err(RegionalCarreauResolutionError::Invalid(
                "finite strictly increasing authored drive discontinuities required",
            ));
        }
    }
    let mut event = discontinuities_s.partition_point(|&time| time <= start);
    let mut next = owner.clone();
    let mut momentum = state.clone();
    let mut caloric = thermal.clone();
    let mut work = RegionalCarreauResolutionWork::default();
    let mut steps = Vec::new();
    let mut drive_samples = Vec::new();
    let mut estimates = Vec::new();
    let mut proposed = controls.maximum_interval_s;
    let maximum_accepted_steps = controls
        .maximum_accepted_steps
        .min(thermal.budget().max_substeps);
    while caloric.time_s() < target_s {
        checkpoint(cx)?;
        if steps.len() + 2 > maximum_accepted_steps {
            return Err(limit("retained native step limit", work));
        }
        let start = caloric.time_s();
        let mut crossed_event = false;
        while event < discontinuities_s.len() && discontinuities_s[event] <= start {
            event += 1;
            crossed_event = true;
        }
        if crossed_event {
            // A new forcing segment gets its own error estimate. Do not carry
            // an event-clipped floating-point width into the next segment:
            // nominally equal decimal periods can differ by one time ULP.
            proposed = controls.maximum_interval_s;
        }
        let segment_end = discontinuities_s
            .get(event)
            .copied()
            .unwrap_or(target_s)
            .min(target_s);
        let mut end = if proposed >= segment_end - start {
            segment_end
        } else {
            start + proposed
        };
        // A representable decimal cap can land one ULP below its authored
        // endpoint (for example 0.074 + 0.037 versus 0.111). Retain the actual
        // full/fine intervals; never skip the remaining physical duration.
        if end.next_up() == segment_end {
            end = segment_end;
        }
        let width = end - start;
        let middle = start + 0.5 * width;
        if !end.is_finite() || !(start < middle && middle < end) {
            return Err(RegionalCarreauResolutionError::Invalid(
                "comparison interval is not representable",
            ));
        }
        let (coarse, coarse_drive) = native_step(
            cx,
            &next,
            &caloric,
            &momentum,
            workspace,
            end,
            step_controls,
            controls,
            &mut work,
            &drive,
            &exchange,
        )?;
        let (first, first_drive) = native_step(
            cx,
            &next,
            &caloric,
            &momentum,
            workspace,
            middle,
            step_controls,
            controls,
            &mut work,
            &drive,
            &exchange,
        )?;
        let (second, second_drive) = native_step(
            cx,
            &first.next,
            &first.thermal,
            &first.momentum,
            workspace,
            end,
            step_controls,
            controls,
            &mut work,
            &drive,
            &exchange,
        )?;
        let estimate = compare(
            cx,
            &momentum,
            &caloric,
            &coarse,
            &first.report,
            &second,
            [&coarse_drive, &first_drive, &second_drive],
            controls,
        )?;
        let resolved = estimate.error_ratio <= 1.0
            && estimate.maximum_reversal_step_s <= controls.maximum_reversal_step_s
            && estimate.maximum_cancelled_volume_m3 <= controls.maximum_cancelled_volume_m3;
        if resolved {
            checkpoint(cx)?;
            next = second.next;
            momentum = second.momentum;
            caloric = second.thermal;
            steps.push(first.report);
            steps.push(second.report);
            drive_samples.push(first_drive);
            drive_samples.push(second_drive);
            estimates.push(estimate);
            proposed = if estimate.error_ratio <= 0.25 && estimate.maximum_reversal_step_s == 0.0 {
                (2.0 * width).min(controls.maximum_interval_s)
            } else {
                width
            };
        } else {
            work.rejected_intervals += 1;
            proposed = 0.5 * width;
            if proposed < controls.minimum_interval_s || start + proposed == start {
                return Err(RegionalCarreauResolutionError::Underresolved { estimate, work });
            }
        }
    }
    checkpoint(cx)?;
    Ok(ResolvedRegionalCarreau {
        next,
        momentum,
        thermal: caloric,
        report: RegionalCarreauTrajectory {
            steps,
            drive_samples,
            estimates,
            work,
        },
    })
}

#[allow(clippy::too_many_arguments)]
fn native_step<'m, T, E>(
    cx: &Cx<'_>,
    owner: &TransientRegionalCarreauTransport,
    thermal: &DepositionSimulation<'m>,
    state: &TransientCarreauState,
    workspace: &mut TransientCarreauWorkspace,
    target_s: f64,
    mut step: TransientRegionalCarreauControls,
    controls: RegionalCarreauResolutionControls,
    work: &mut RegionalCarreauResolutionWork,
    drive: &DriveResolver<'_, 'm, E>,
    exchange: &impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &TransientCarreauStep,
        &RegionalCarreauDriveSample,
    ) -> Result<CarreauTransportExchange<T>, E>,
) -> Result<
    (
        TransientRegionalCarreauStep<'m, T>,
        RegionalCarreauDriveSample,
    ),
    RegionalCarreauResolutionError<E>,
> {
    checkpoint(cx)?;
    if work.native_steps == controls.maximum_native_steps {
        return Err(limit("native step attempt limit", *work));
    }
    let remaining = controls.maximum_total_branch_solves - work.completed_branch_solves;
    if remaining == 0 {
        return Err(limit("total original branch solve limit", *work));
    }
    step.maximum_total_branch_solves = step.maximum_total_branch_solves.min(remaining);
    let head_remaining = controls.maximum_total_head_evaluations - work.completed_head_evaluations;
    step.hydraulics.head.maximum_source_evaluations = step
        .hydraulics
        .head
        .maximum_source_evaluations
        .min(head_remaining);
    let start_s = thermal.time_s();
    let sample_time_s = start_s + 0.5 * (target_s - start_s);
    if !(start_s < sample_time_s && sample_time_s < target_s) {
        return Err(RegionalCarreauResolutionError::Invalid(
            "native midpoint drive sample is not representable inside its interval",
        ));
    }
    work.native_steps += 1;
    let guess = vec![0.0; thermal.cells().len()];
    let result = match drive {
        DriveResolver::Prescribed(resolve) => {
            work.drive_evaluations = work
                .drive_evaluations
                .checked_add(1)
                .ok_or_else(|| limit("drive evaluation count limit", *work))?;
            let sample = RegionalCarreauDriveSample {
                start_s,
                end_s: target_s,
                sample_time_s,
                thermal_time_s: None,
                drive: resolve(cx, sample_time_s).map_err(|source| {
                    RegionalCarreauResolutionError::Drive {
                        start_s,
                        end_s: target_s,
                        sample_time_s,
                        work: *work,
                        source,
                    }
                })?,
            };
            checkpoint(cx)?;
            let result = owner
                .advance_with_resolved_drive(
                    cx,
                    thermal,
                    state,
                    workspace,
                    target_s,
                    step,
                    &guess,
                    |_, _, _| Ok::<_, E>(&sample.drive),
                    |cx, candidate, flow, _| exchange(cx, candidate, flow, &sample),
                )
                .map(|(result, _)| result);
            result.map(|result| (result, sample))
        }
        DriveResolver::Thermal(resolve) => {
            // The native endpoint owner admits at most 1,000,000 evaluations.
            // Count before querying so a source refusal is included as work.
            let evaluations = Cell::new(0usize);
            let result = owner.advance_with_resolved_drive(
                cx,
                thermal,
                state,
                workspace,
                target_s,
                step,
                &guess,
                |cx, candidate, source| {
                    evaluations.set(evaluations.get() + 1);
                    resolve(cx, sample_time_s, candidate, source).map(|drive| {
                        RegionalCarreauDriveSample {
                            start_s,
                            end_s: target_s,
                            sample_time_s,
                            thermal_time_s: Some(candidate.time_s()),
                            drive,
                        }
                    })
                },
                |cx, candidate, flow, sample| exchange(cx, candidate, flow, sample),
            );
            work.drive_evaluations = work
                .drive_evaluations
                .checked_add(evaluations.get())
                .ok_or_else(|| limit("drive evaluation count limit", *work))?;
            result
        }
    };
    match result {
        Ok((result, sample)) => {
            if result.report.total_branch_solves > remaining {
                return Err(RegionalCarreauResolutionError::Invalid(
                    "original branch work exceeded remaining limit",
                ));
            }
            work.completed_branch_solves += result.report.total_branch_solves;
            if result.report.total_head_evaluations > head_remaining {
                return Err(RegionalCarreauResolutionError::Invalid(
                    "original equipment work exceeded remaining trajectory limit",
                ));
            }
            work.completed_head_evaluations += result.report.total_head_evaluations;
            Ok((result, sample))
        }
        Err(source) => Err(RegionalCarreauResolutionError::Step {
            start_s: thermal.time_s(),
            end_s: target_s,
            work: *work,
            source,
        }),
    }
}

#[allow(clippy::too_many_arguments)]
fn compare<T, E>(
    cx: &Cx<'_>,
    initial_momentum: &TransientCarreauState,
    initial: &DepositionSimulation<'_>,
    coarse: &TransientRegionalCarreauStep<'_, T>,
    first: &TransientRegionalCarreauReport<T>,
    fine: &TransientRegionalCarreauStep<'_, T>,
    drives: [&RegionalCarreauDriveSample; 3],
    controls: RegionalCarreauResolutionControls,
) -> Result<RegionalCarreauResolutionEstimate, RegionalCarreauResolutionError<E>> {
    let mut estimate = RegionalCarreauResolutionEstimate {
        start_s: initial.time_s(),
        end_s: fine.thermal.time_s(),
        maximum_cell_difference_j: 0.0,
        maximum_velocity_difference_m_s: 0.0,
        maximum_flow_difference_m3_s: 0.0,
        maximum_transfer_difference_j: 0.0,
        error_ratio: 0.0,
        maximum_reversal_step_s: 0.0,
        maximum_cancelled_volume_m3: 0.0,
    };
    if coarse.thermal.time_s() != fine.thermal.time_s()
        || coarse.thermal.cells().len() != initial.cells().len()
        || fine.thermal.cells().len() != initial.cells().len()
        || coarse.momentum.branches().len() != initial_momentum.branches().len()
        || fine.momentum.branches().len() != initial_momentum.branches().len()
    {
        return Err(RegionalCarreauResolutionError::Invalid(
            "native comparison inventory changed",
        ));
    }
    for ((id, old), ((a_id, a), (b_id, b))) in initial
        .cells()
        .iter()
        .zip(coarse.thermal.cells().iter().zip(fine.thermal.cells()))
    {
        checkpoint(cx)?;
        if id != a_id
            || id != b_id
            || old.mass_kg != a.mass_kg
            || old.mass_kg != b.mass_kg
            || old.material_id != a.material_id
            || old.material_id != b.material_id
            || old.geometry != a.geometry
            || old.geometry != b.geometry
        {
            return Err(RegionalCarreauResolutionError::Invalid(
                "native comparison changed physical storage",
            ));
        }
        let difference = (old.mass_kg * (b.specific_h_j_kg - a.specific_h_j_kg)).abs();
        let change = (old.mass_kg * (b.specific_h_j_kg - old.specific_h_j_kg))
            .abs()
            .max((old.mass_kg * (a.specific_h_j_kg - old.specific_h_j_kg)).abs());
        let scale = controls.absolute_cell_error_j + controls.relative_cell_error * change;
        update(&mut estimate.error_ratio, difference, scale)?;
        estimate.maximum_cell_difference_j = estimate.maximum_cell_difference_j.max(difference);
    }
    for ((old, a), b) in initial_momentum
        .branches()
        .iter()
        .zip(coarse.momentum.branches())
        .zip(fine.momentum.branches())
    {
        if old.velocities_m_s().len() != a.velocities_m_s().len()
            || old.velocities_m_s().len() != b.velocities_m_s().len()
        {
            return Err(RegionalCarreauResolutionError::Invalid(
                "native radial coordinates changed",
            ));
        }
        for ((old, a), b) in old
            .velocities_m_s()
            .iter()
            .zip(a.velocities_m_s())
            .zip(b.velocities_m_s())
        {
            checkpoint(cx)?;
            let difference = (b - a).abs();
            let scale = controls.absolute_velocity_error_m_s
                + controls.relative_velocity_error * old.abs().max(a.abs()).max(b.abs());
            update(&mut estimate.error_ratio, difference, scale)?;
            estimate.maximum_velocity_difference_m_s =
                estimate.maximum_velocity_difference_m_s.max(difference);
        }
    }
    let transfers = [
        (&coarse.report.hydraulics, drives[0]),
        (&first.hydraulics, drives[1]),
        (&fine.report.hydraulics, drives[2]),
    ];
    for (branch, ((old, a), b)) in first
        .hydraulics
        .branches
        .iter()
        .zip(&coarse.report.hydraulics.branches)
        .zip(&fine.report.hydraulics.branches)
        .enumerate()
    {
        checkpoint(cx)?;
        transfer(
            &mut estimate,
            transfers.map(|(flow, _)| flow.branches[branch].viscous_heat_j),
            controls,
        )?;
        transfer(
            &mut estimate,
            transfers.map(|(flow, drive)| {
                let head = drive
                    .drive
                    .branch_heads
                    .as_ref()
                    .map_or(0.0, |heads| heads[branch].value());
                head * flow.branches[branch].midpoint_volume_flow_m3_s
                    * flow.branches[branch].duration_s
            }),
            controls,
        )?;
        transfer(
            &mut estimate,
            transfers.map(|(flow, _)| {
                flow.flow_heads[branch]
                    .as_ref()
                    .map_or(0.0, |head| head.work_j)
            }),
            controls,
        )?;
        let equipment = [&coarse.report, first, &fine.report].map(|report| {
            report.equipment.branches[branch].map_or([0.0; 7], |energy| {
                [
                    energy.shaft_input_j,
                    energy.electrical_input_j,
                    energy.pump_loss_j,
                    energy.drive_loss_j,
                    energy.pump_heat_j,
                    energy.drive_heat_j,
                    energy.surroundings_heat_j,
                ]
            })
        });
        for channel in 0..7 {
            transfer(
                &mut estimate,
                equipment.map(|energy| energy[channel]),
                controls,
            )?;
        }
        let difference = (b.final_volume_flow_m3_s - a.final_volume_flow_m3_s).abs();
        let scale = controls.absolute_flow_error_m3_s
            + controls.relative_flow_error
                * old
                    .initial_volume_flow_m3_s
                    .abs()
                    .max(a.final_volume_flow_m3_s.abs())
                    .max(b.final_volume_flow_m3_s.abs());
        update(&mut estimate.error_ratio, difference, scale)?;
        estimate.maximum_flow_difference_m3_s =
            estimate.maximum_flow_difference_m3_s.max(difference);
    }
    for report in [first, &fine.report] {
        for branch in &report.hydraulics.branches {
            checkpoint(cx)?;
            let (old, new) = (
                branch.initial_volume_flow_m3_s,
                branch.final_volume_flow_m3_s,
            );
            if (old < 0.0 && new > 0.0) || (old > 0.0 && new < 0.0) {
                let small = old.abs().min(new.abs());
                let fraction = small / old.abs().max(new.abs());
                let volume = (0.5 * branch.duration_s * small) * (fraction / (1.0 + fraction));
                if !nonnegative(volume) {
                    return Err(RegionalCarreauResolutionError::Invalid(
                        "unrepresentable reversal volume",
                    ));
                }
                estimate.maximum_reversal_step_s =
                    estimate.maximum_reversal_step_s.max(branch.duration_s);
                estimate.maximum_cancelled_volume_m3 =
                    estimate.maximum_cancelled_volume_m3.max(volume);
            }
        }
    }
    Ok(estimate)
}

fn transfer<E>(
    estimate: &mut RegionalCarreauResolutionEstimate,
    [coarse, first, second]: [f64; 3],
    controls: RegionalCarreauResolutionControls,
) -> Result<(), RegionalCarreauResolutionError<E>> {
    let mut residual = crate::graph::Sum::default();
    let error =
        |_| RegionalCarreauResolutionError::Invalid("unrepresentable original transfer comparison");
    for value in [coarse, -first, -second] {
        residual.add(value).map_err(error)?;
    }
    let signed = residual.value().map_err(error)?.abs();
    let fine_magnitude = first.abs() + second.abs();
    let difference = signed.max((coarse.abs() - fine_magnitude).abs());
    let scale = controls.absolute_transfer_error_j
        + controls.relative_transfer_error * coarse.abs().max(fine_magnitude);
    update(&mut estimate.error_ratio, difference, scale)?;
    estimate.maximum_transfer_difference_j = estimate.maximum_transfer_difference_j.max(difference);
    Ok(())
}

fn update<E>(
    ratio: &mut f64,
    difference: f64,
    scale: f64,
) -> Result<(), RegionalCarreauResolutionError<E>> {
    if !nonnegative(difference) || !nonnegative(scale) {
        return Err(RegionalCarreauResolutionError::Invalid(
            "unrepresentable native state comparison",
        ));
    }
    *ratio = ratio.max(if difference == 0.0 {
        0.0
    } else {
        difference / scale
    });
    Ok(())
}
fn limit<E>(
    reason: &'static str,
    work: RegionalCarreauResolutionWork,
) -> RegionalCarreauResolutionError<E> {
    RegionalCarreauResolutionError::WorkLimit { reason, work }
}
fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), RegionalCarreauResolutionError<E>> {
    thermal_checkpoint(cx).map_err(RegionalCarreauResolutionError::Thermal)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn tolerances(absolute: f64, relative: f64) -> bool {
    nonnegative(absolute)
        && nonnegative(relative)
        && relative < 1.0
        && (absolute > 0.0 || relative > 0.0)
}
