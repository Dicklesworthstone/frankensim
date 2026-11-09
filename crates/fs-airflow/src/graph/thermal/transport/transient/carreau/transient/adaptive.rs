//! Time resolution of original radial momentum and regional enthalpy.
//!
//! Compare one native interval with two native half intervals, retaining only
//! the fine native states. Every radial annulus and thermal cell participates:
//! equal volume flows do not imply equal shear histories. Separate bounds on
//! sign-changing interval duration and cancelled directional volume prevent a
//! nearly zero mean flow from hiding unresolved counterflow at uniform heat.
//! This follows the existing inertial transport controller's dyadic pattern;
//! it adds neither an extrapolated physical state nor another heat integrator.

use fs_exec::Cx;
use fs_qty::Pressure;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

use super::{
    CarreauTransportExchange, FixedPressure, FlowInjection, ImplicitTransportError,
    TransientCarreauState, TransientCarreauStep, TransientCarreauTransportError,
    TransientCarreauWorkspace, TransientRegionalCarreauControls, TransientRegionalCarreauReport,
    TransientRegionalCarreauStep, TransientRegionalCarreauTransport, thermal_checkpoint,
};

/// Local discretization, counterflow resolution and total work limits.
#[derive(Clone, Copy, Debug)]
pub struct RegionalCarreauResolutionControls {
    /// Maximum full comparison interval. Retained native intervals are half
    /// this width or smaller. Split discontinuous forcing/events beforehand.
    pub maximum_interval_s: f64,
    /// Refinement floor; a shorter clipped final interval cannot refine again.
    pub minimum_interval_s: f64,
    /// Includes discarded coarse/fine attempts and the terminal refused step.
    pub maximum_native_steps: usize,
    pub maximum_accepted_steps: usize,
    /// All original branch solves, including discarded comparisons and every
    /// thermal residual. Each native step is clipped to the remaining limit.
    pub maximum_total_branch_solves: u64,
    pub absolute_cell_error_j: f64,
    /// Relative to interval energy CHANGES, independent of caloric datum.
    pub relative_cell_error: f64,
    pub absolute_velocity_error_m_s: f64,
    pub relative_velocity_error: f64,
    pub absolute_flow_error_m3_s: f64,
    pub relative_flow_error: f64,
    /// Longest retained substep with opposite endpoint flow signs.
    pub maximum_reversal_step_s: f64,
    /// Largest min(V+, V-) of a branch's native affine flow trajectory, per
    /// retained substep. This bounds cancellation, not the net displaced volume.
    pub maximum_cancelled_volume_m3: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RegionalCarreauResolutionWork {
    pub native_steps: usize,
    pub rejected_intervals: usize,
    /// Exact counts from completed native solves, including discarded ones.
    /// On terminal native failure its original remaining budget still bounds
    /// attempted work, but that failed owner's partial count is not inferred.
    pub completed_branch_solves: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct RegionalCarreauResolutionEstimate {
    pub start_s: f64,
    pub end_s: f64,
    pub maximum_cell_difference_j: f64,
    pub maximum_velocity_difference_m_s: f64,
    pub maximum_flow_difference_m3_s: f64,
    pub error_ratio: f64,
    pub maximum_reversal_step_s: f64,
    pub maximum_cancelled_volume_m3: f64,
}

/// All accepted native substeps, in chronological order. Their original
/// pressure/head work, viscous heat, kinetic and caloric changes remain separate.
/// Only these reports contribute to the requested physical accounting window.
#[derive(Debug)]
pub struct RegionalCarreauTrajectory<T> {
    pub steps: Vec<TransientRegionalCarreauReport<T>>,
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
        || !positive(controls.maximum_reversal_step_s)
        || !nonnegative(controls.maximum_cancelled_volume_m3)
    {
        return Err(RegionalCarreauResolutionError::Invalid(
            "finite synchronized interval and bounded native resolution controls required",
        ));
    }
    let mut next = owner.clone();
    let mut momentum = state.clone();
    let mut caloric = thermal.clone();
    let mut work = RegionalCarreauResolutionWork::default();
    let mut steps = Vec::new();
    let mut estimates = Vec::new();
    let mut proposed = controls.maximum_interval_s;
    while caloric.time_s() < target_s {
        checkpoint(cx)?;
        if steps.len() + 2 > controls.maximum_accepted_steps {
            return Err(limit("retained native step limit", work));
        }
        let start = caloric.time_s();
        let end = if proposed >= target_s - start {
            target_s
        } else {
            start + proposed
        };
        let width = end - start;
        let middle = start + 0.5 * width;
        if !end.is_finite() || !(start < middle && middle < end) {
            return Err(RegionalCarreauResolutionError::Invalid(
                "comparison interval is not representable",
            ));
        }
        let coarse = native_step(
            cx,
            &next,
            &caloric,
            &momentum,
            workspace,
            end,
            boundaries,
            injections,
            branch_heads,
            step_controls,
            controls,
            &mut work,
            &exchange,
        )?;
        let first = native_step(
            cx,
            &next,
            &caloric,
            &momentum,
            workspace,
            middle,
            boundaries,
            injections,
            branch_heads,
            step_controls,
            controls,
            &mut work,
            &exchange,
        )?;
        let second = native_step(
            cx,
            &first.next,
            &first.thermal,
            &first.momentum,
            workspace,
            end,
            boundaries,
            injections,
            branch_heads,
            step_controls,
            controls,
            &mut work,
            &exchange,
        )?;
        let estimate = compare(
            cx,
            &momentum,
            &caloric,
            &coarse,
            &first.report,
            &second,
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
    boundaries: &[FixedPressure],
    injections: &[FlowInjection],
    branch_heads: Option<&[Pressure]>,
    mut step: TransientRegionalCarreauControls,
    controls: RegionalCarreauResolutionControls,
    work: &mut RegionalCarreauResolutionWork,
    exchange: &impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &TransientCarreauStep,
    ) -> Result<CarreauTransportExchange<T>, E>,
) -> Result<TransientRegionalCarreauStep<'m, T>, RegionalCarreauResolutionError<E>> {
    checkpoint(cx)?;
    if work.native_steps == controls.maximum_native_steps {
        return Err(limit("native step attempt limit", *work));
    }
    let remaining = controls.maximum_total_branch_solves - work.completed_branch_solves;
    if remaining == 0 {
        return Err(limit("total original branch solve limit", *work));
    }
    step.maximum_total_branch_solves = step.maximum_total_branch_solves.min(remaining);
    work.native_steps += 1;
    let guess = vec![0.0; thermal.cells().len()];
    let result = match branch_heads {
        Some(heads) => owner.advance_with_heads(
            cx, thermal, state, workspace, target_s, boundaries, injections, heads, step, &guess,
            exchange,
        ),
        None => owner.advance(
            cx, thermal, state, workspace, target_s, boundaries, injections, step, &guess, exchange,
        ),
    };
    match result {
        Ok(result) => {
            if result.report.total_branch_solves > remaining {
                return Err(RegionalCarreauResolutionError::Invalid(
                    "original branch work exceeded remaining limit",
                ));
            }
            work.completed_branch_solves += result.report.total_branch_solves;
            Ok(result)
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
    controls: RegionalCarreauResolutionControls,
) -> Result<RegionalCarreauResolutionEstimate, RegionalCarreauResolutionError<E>> {
    let mut estimate = RegionalCarreauResolutionEstimate {
        start_s: initial.time_s(),
        end_s: fine.thermal.time_s(),
        maximum_cell_difference_j: 0.0,
        maximum_velocity_difference_m_s: 0.0,
        maximum_flow_difference_m3_s: 0.0,
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
    for ((old, a), b) in first
        .hydraulics
        .branches
        .iter()
        .zip(&coarse.report.hydraulics.branches)
        .zip(&fine.report.hydraulics.branches)
    {
        checkpoint(cx)?;
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
