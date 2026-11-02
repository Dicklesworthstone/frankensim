//! Time resolution of the original coupled momentum/enthalpy step.
//!
//! A full native step is compared with two native half steps. Only the two
//! independently admitted half steps can be retained. No extrapolated state,
//! alternate momentum law or second thermal integrator is introduced. The
//! differences are local discretization estimates, not rigorous error bounds.
//! An explicit maximum interval remains necessary for time-dependent forcing:
//! endpoint sampling cannot detect arbitrarily fast unsampled oscillations.

use std::cell::Cell;

use fs_exec::Cx;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

use super::super::distributed::DistributedGraphInventory;
use super::{
    InertialGraph, InertialGraphHeat, InertialGraphReport, InertialThermalControls,
    InertialThermalError, InertialThermalProblem, InertialThermalReport, ResolvedInertialGraph,
    advance_inertial_graph_with_inventory, thermal_checkpoint,
};

/// Work, state-error and reversal-resolution controls for one requested window.
#[derive(Clone, Copy, Debug)]
pub struct InertialResolutionControls {
    /// Largest full-versus-two-half comparison interval, seconds. Accepted
    /// native substeps are half this interval or smaller.
    pub maximum_interval_s: f64,
    /// Refinement floor for comparison intervals. A clipped final interval may
    /// be smaller, but it cannot be refined again below this floor.
    pub minimum_interval_s: f64,
    /// All native step attempts, including discarded full and half steps.
    pub maximum_native_steps: usize,
    /// Retained native reports, including both halves of every accepted pair.
    pub maximum_accepted_steps: usize,
    /// All calls to the original resolve callback, including discarded trials
    /// and every native actual-flow refresh. The native per-step cap also applies.
    pub maximum_endpoint_evaluations: usize,
    pub absolute_cell_error_j: f64,
    /// Relative to interval ENERGY CHANGES, not formation-energy datum.
    pub relative_cell_error: f64,
    pub absolute_flow_error_m3_s: f64,
    pub relative_flow_error: f64,
    /// Longest accepted native substep whose endpoints have opposite signs.
    /// This is temporal resolution, not a zero-crossing location certificate.
    pub maximum_reversal_step_s: f64,
    /// Maximum volume in EACH direction hidden by averaging an affine signed
    /// reversal, per branch and native substep. This is min(V+, V-), not net Q dt.
    pub maximum_cancelled_volume_m3: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InertialResolutionWork {
    pub native_steps: usize,
    pub endpoint_evaluations: usize,
    pub rejected_intervals: usize,
}

/// Comparison at the same endpoint, and reversal resolution of the two native
/// half steps that were actually retained. Momentum/thermal audits remain on
/// the original reports, rather than being replaced by this estimate.
#[derive(Clone, Copy, Debug)]
pub struct InertialResolutionEstimate {
    pub start_s: f64,
    pub end_s: f64,
    pub maximum_cell_difference_j: f64,
    pub maximum_flow_difference_m3_s: f64,
    pub error_ratio: f64,
    pub maximum_reversal_step_s: f64,
    pub maximum_cancelled_volume_m3: f64,
}

#[derive(Debug)]
pub struct ResolvedInertialTrajectory<P, T> {
    /// Chronological original native reports. Sum physical accounts only over
    /// these accepted intervals; discarded comparison solves carry no credit.
    pub steps: Vec<InertialThermalReport<P, T>>,
    pub estimates: Vec<InertialResolutionEstimate>,
    pub work: InertialResolutionWork,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InertialResolutionBudget {
    NativeSteps,
    AcceptedSteps,
    EndpointEvaluations,
}

#[derive(Debug)]
pub enum InertialResolutionError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Budget {
        kind: InertialResolutionBudget,
        performed: usize,
        maximum: usize,
        work: InertialResolutionWork,
    },
    Step {
        start_s: f64,
        end_s: f64,
        work: InertialResolutionWork,
        source: InertialThermalError<E>,
    },
    Underresolved {
        estimate: InertialResolutionEstimate,
        work: InertialResolutionWork,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for InertialResolutionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved inertial convection: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for InertialResolutionError<E> {}

/// Resolve a forward, event-free window using the original coupled owners.
///
/// `inventory=None` uses original junction storage; `Some` retains the explicitly
/// authored pipe cells. Both callbacks have exactly the contract of the native
/// inertial driver and must be pure. Every source, correlation, numerical and
/// physical refusal is terminal, including during a discarded comparison solve;
/// the adapter never reinterprets an unsupported law as a timestep rejection.
/// Native scheduled-event splitting remains the caller's responsibility.
///
/// Error estimates compare EVERY active cell's energy and EVERY branch flow.
/// The coupled heat method is first order, so the full/half difference is the
/// usual first-order local estimate for the fine result. The method does not
/// claim a global enclosure, resolve arbitrary forcing frequencies or repair
/// spatial upwind diffusion. Refine the authored cells and maximum interval.
///
/// All inputs remain immutable. The entire requested window either returns
/// synchronized native owners and their accepted reports, or returns an error
/// without publishing any of its local intermediate states.
#[allow(clippy::too_many_arguments)]
pub fn advance_inertial_graph_resolved<'m, P, T, E>(
    cx: &Cx<'_>,
    problem: &InertialThermalProblem<'_>,
    inventory: Option<&DistributedGraphInventory>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    step_controls: InertialThermalControls,
    controls: InertialResolutionControls,
    resolve: impl Fn(&Cx<'_>, &DepositionSimulation<'m>, &[f64]) -> Result<ResolvedInertialGraph<P>, E>,
    exchange: impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &P,
        &InertialGraphReport,
    ) -> Result<InertialGraphHeat<T>, E>,
) -> Result<
    (
        InertialGraph,
        DepositionSimulation<'m>,
        ResolvedInertialTrajectory<P, T>,
    ),
    InertialResolutionError<E>,
> {
    checkpoint(cx)?;
    let start_s = thermal.time_s();
    if !target_s.is_finite()
        || target_s <= start_s
        || !(target_s - start_s).is_finite()
        || start_s != problem.momentum.accepted_time_s()
        || !positive(controls.minimum_interval_s)
        || !positive(controls.maximum_interval_s)
        || controls.minimum_interval_s > controls.maximum_interval_s
        || !(1..=1_000_000).contains(&controls.maximum_native_steps)
        || !(2..=4096).contains(&controls.maximum_accepted_steps)
        || !(1..=1_000_000).contains(&controls.maximum_endpoint_evaluations)
        || !tolerances(controls.absolute_cell_error_j, controls.relative_cell_error)
        || !tolerances(
            controls.absolute_flow_error_m3_s,
            controls.relative_flow_error,
        )
        || !positive(controls.maximum_reversal_step_s)
        || !nonnegative(controls.maximum_cancelled_volume_m3)
    {
        return Err(InertialResolutionError::Invalid(
            "finite synchronized forward interval and bounded resolution controls required",
        ));
    }
    let mut current_momentum = problem.momentum.clone();
    let mut current_thermal = thermal.clone();
    let mut work = InertialResolutionWork::default();
    let mut steps = Vec::new();
    let mut estimates = Vec::new();
    let mut proposed = controls.maximum_interval_s;
    while current_thermal.time_s() < target_s {
        checkpoint(cx)?;
        if steps.len() + 2 > controls.maximum_accepted_steps {
            return Err(budget(
                InertialResolutionBudget::AcceptedSteps,
                steps.len(),
                controls.maximum_accepted_steps,
                work,
            ));
        }
        let start = current_thermal.time_s();
        let remaining = target_s - start;
        let end = if proposed >= remaining {
            target_s
        } else {
            start + proposed
        };
        let width = end - start;
        let middle = start + 0.5 * width;
        if !end.is_finite() || !(start < middle && middle < end) {
            return Err(InertialResolutionError::Invalid(
                "native comparison interval is not representable",
            ));
        }
        let (coarse_momentum, coarse_thermal, coarse_report) = native_step(
            cx,
            problem,
            inventory,
            &current_momentum,
            &current_thermal,
            end,
            step_controls,
            controls,
            &mut work,
            &resolve,
            &exchange,
        )?;
        drop(coarse_report);
        let (middle_momentum, middle_thermal, first_report) = native_step(
            cx,
            problem,
            inventory,
            &current_momentum,
            &current_thermal,
            middle,
            step_controls,
            controls,
            &mut work,
            &resolve,
            &exchange,
        )?;
        let (fine_momentum, fine_thermal, second_report) = native_step(
            cx,
            problem,
            inventory,
            &middle_momentum,
            &middle_thermal,
            end,
            step_controls,
            controls,
            &mut work,
            &resolve,
            &exchange,
        )?;
        let estimate = compare(
            cx,
            &current_momentum,
            &current_thermal,
            &coarse_momentum,
            &coarse_thermal,
            &fine_momentum,
            &fine_thermal,
            [&first_report.hydraulics, &second_report.hydraulics],
            controls,
        )?;
        let resolved = estimate.error_ratio <= 1.0
            && estimate.maximum_reversal_step_s <= controls.maximum_reversal_step_s
            && estimate.maximum_cancelled_volume_m3 <= controls.maximum_cancelled_volume_m3;
        if resolved {
            checkpoint(cx)?;
            current_momentum = fine_momentum;
            current_thermal = fine_thermal;
            steps.push(first_report);
            steps.push(second_report);
            estimates.push(estimate);
            // Deterministic dyadic growth; no embedded RK controller or
            // extrapolated physics state is substituted for the native owners.
            proposed = if estimate.error_ratio <= 0.25 && estimate.maximum_reversal_step_s == 0.0 {
                (2.0 * width).min(controls.maximum_interval_s)
            } else {
                width
            };
        } else {
            work.rejected_intervals += 1;
            proposed = 0.5 * width;
            if proposed < controls.minimum_interval_s || start + proposed == start {
                return Err(InertialResolutionError::Underresolved { estimate, work });
            }
        }
    }
    checkpoint(cx)?;
    Ok((
        current_momentum,
        current_thermal,
        ResolvedInertialTrajectory {
            steps,
            estimates,
            work,
        },
    ))
}

#[allow(clippy::too_many_arguments)]
fn native_step<'m, P, T, E>(
    cx: &Cx<'_>,
    problem: &InertialThermalProblem<'_>,
    inventory: Option<&DistributedGraphInventory>,
    momentum: &InertialGraph,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    mut step: InertialThermalControls,
    controls: InertialResolutionControls,
    work: &mut InertialResolutionWork,
    resolve: &impl Fn(&Cx<'_>, &DepositionSimulation<'m>, &[f64]) -> Result<ResolvedInertialGraph<P>, E>,
    exchange: &impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &P,
        &InertialGraphReport,
    ) -> Result<InertialGraphHeat<T>, E>,
) -> Result<
    (
        InertialGraph,
        DepositionSimulation<'m>,
        InertialThermalReport<P, T>,
    ),
    InertialResolutionError<E>,
> {
    checkpoint(cx)?;
    if work.native_steps == controls.maximum_native_steps {
        return Err(budget(
            InertialResolutionBudget::NativeSteps,
            work.native_steps,
            controls.maximum_native_steps,
            *work,
        ));
    }
    let remaining = controls.maximum_endpoint_evaluations - work.endpoint_evaluations;
    if remaining == 0 {
        return Err(budget(
            InertialResolutionBudget::EndpointEvaluations,
            work.endpoint_evaluations,
            controls.maximum_endpoint_evaluations,
            *work,
        ));
    }
    step.maximum_endpoint_evaluations = step.maximum_endpoint_evaluations.min(remaining);
    work.native_steps += 1;
    let calls = Cell::new(0usize);
    let local = InertialThermalProblem {
        momentum,
        node_cells: problem.node_cells,
        maximum_connections: problem.maximum_connections,
    };
    let result = advance_inertial_graph_with_inventory(
        cx,
        &local,
        inventory,
        thermal,
        target_s,
        step,
        &vec![0.0; thermal.cells().len()],
        &momentum.volume_flows_m3_s(),
        |cx, candidate, rates| {
            calls.set(calls.get() + 1);
            resolve(cx, candidate, rates)
        },
        exchange,
    );
    // Count original callback work on failure as well as success. The native
    // per-step cap was clipped to the remaining window allowance beforehand.
    work.endpoint_evaluations += calls.get();
    match result {
        Ok((momentum, thermal, report)) => {
            if report.endpoint_evaluations != calls.get() {
                return Err(InertialResolutionError::Invalid(
                    "native endpoint accounting changed",
                ));
            }
            Ok((momentum, thermal, report))
        }
        Err(source) => Err(InertialResolutionError::Step {
            start_s: thermal.time_s(),
            end_s: target_s,
            work: *work,
            source,
        }),
    }
}

#[allow(clippy::too_many_arguments)]
fn compare<E>(
    cx: &Cx<'_>,
    initial_momentum: &InertialGraph,
    initial: &DepositionSimulation<'_>,
    coarse_momentum: &InertialGraph,
    coarse: &DepositionSimulation<'_>,
    fine_momentum: &InertialGraph,
    fine: &DepositionSimulation<'_>,
    reports: [&InertialGraphReport; 2],
    controls: InertialResolutionControls,
) -> Result<InertialResolutionEstimate, InertialResolutionError<E>> {
    let mut estimate = InertialResolutionEstimate {
        start_s: initial.time_s(),
        end_s: fine.time_s(),
        maximum_cell_difference_j: 0.0,
        maximum_flow_difference_m3_s: 0.0,
        error_ratio: 0.0,
        maximum_reversal_step_s: 0.0,
        maximum_cancelled_volume_m3: 0.0,
    };
    if coarse.time_s() != fine.time_s()
        || coarse.cells().len() != initial.cells().len()
        || fine.cells().len() != initial.cells().len()
    {
        return Err(InertialResolutionError::Invalid(
            "native comparison inventories or times differ",
        ));
    }
    for ((id, old), ((coarse_id, a), (fine_id, b))) in initial
        .cells()
        .iter()
        .zip(coarse.cells().iter().zip(fine.cells()))
    {
        checkpoint(cx)?;
        if id != coarse_id
            || id != fine_id
            || old.mass_kg != a.mass_kg
            || old.mass_kg != b.mass_kg
            || old.material_id != a.material_id
            || old.material_id != b.material_id
            || old.geometry != a.geometry
            || old.geometry != b.geometry
        {
            return Err(InertialResolutionError::Invalid(
                "native comparison changed physical storage",
            ));
        }
        let difference = (old.mass_kg * (b.specific_h_j_kg - a.specific_h_j_kg)).abs();
        let change = (old.mass_kg * (b.specific_h_j_kg - old.specific_h_j_kg))
            .abs()
            .max((old.mass_kg * (a.specific_h_j_kg - old.specific_h_j_kg)).abs());
        let scale = controls.absolute_cell_error_j + controls.relative_cell_error * change;
        if !nonnegative(difference) || !nonnegative(scale) {
            return Err(InertialResolutionError::Invalid(
                "unrepresentable native energy comparison",
            ));
        }
        estimate.maximum_cell_difference_j = estimate.maximum_cell_difference_j.max(difference);
        estimate.error_ratio = estimate.error_ratio.max(ratio(difference, scale));
    }
    let initial_rates = initial_momentum.volume_flows_m3_s();
    let coarse_rates = coarse_momentum.volume_flows_m3_s();
    let fine_rates = fine_momentum.volume_flows_m3_s();
    if initial_rates.len() != coarse_rates.len() || initial_rates.len() != fine_rates.len() {
        return Err(InertialResolutionError::Invalid(
            "native branch coordinates changed",
        ));
    }
    for ((old, a), b) in initial_rates.iter().zip(coarse_rates).zip(fine_rates) {
        checkpoint(cx)?;
        let difference = (b - a).abs();
        let scale = controls.absolute_flow_error_m3_s
            + controls.relative_flow_error * old.abs().max(a.abs()).max(b.abs());
        if !nonnegative(difference) || !nonnegative(scale) {
            return Err(InertialResolutionError::Invalid(
                "unrepresentable native flow comparison",
            ));
        }
        estimate.maximum_flow_difference_m3_s =
            estimate.maximum_flow_difference_m3_s.max(difference);
        estimate.error_ratio = estimate.error_ratio.max(ratio(difference, scale));
    }
    for report in reports {
        for (&old, &new) in report
            .initial_flows_m3_s
            .iter()
            .zip(&report.final_flows_m3_s)
        {
            checkpoint(cx)?;
            if (old < 0.0 && new > 0.0) || (old > 0.0 && new < 0.0) {
                let small = old.abs().min(new.abs());
                let large = old.abs().max(new.abs());
                let fraction = small / large;
                let volume = (0.5 * report.duration_s * small) * (fraction / (1.0 + fraction));
                if !nonnegative(volume) {
                    return Err(InertialResolutionError::Invalid(
                        "unrepresentable reversing volume",
                    ));
                }
                estimate.maximum_reversal_step_s =
                    estimate.maximum_reversal_step_s.max(report.duration_s);
                estimate.maximum_cancelled_volume_m3 =
                    estimate.maximum_cancelled_volume_m3.max(volume);
            }
        }
    }
    Ok(estimate)
}

fn budget<E>(
    kind: InertialResolutionBudget,
    performed: usize,
    maximum: usize,
    work: InertialResolutionWork,
) -> InertialResolutionError<E> {
    InertialResolutionError::Budget {
        kind,
        performed,
        maximum,
        work,
    }
}
fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), InertialResolutionError<E>> {
    thermal_checkpoint(cx).map_err(InertialResolutionError::Thermal)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn tolerances(absolute: f64, relative: f64) -> bool {
    nonnegative(absolute)
        && nonnegative(relative)
        && relative < 1.0
        && (absolute > 0.0 || relative > 0.0)
}
fn ratio(difference: f64, scale: f64) -> f64 {
    if difference == 0.0 {
        0.0
    } else {
        difference / scale
    }
}

#[cfg(test)]
mod tests;
