//! Retained graph momentum coupled to the original finite enthalpy timeline.
//!
//! One native descriptor PHS proposal enforces network continuity and advances
//! branch impulses. One native backward-Euler proposal advances original cell
//! enthalpies. Their joint residual uses the actual PHS step-average rates for
//! conservative upwind transport, including startup, shutdown and reversal.
//! Source callbacks resolve temperature-dependent losses, gravity and exchange;
//! neither accepted owner changes until both independent physical gates pass.

pub mod adaptive;

use std::cell::Cell;
use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use fs_thermal::enthalpy::advection::MassBoundary;
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};

use super::distributed::DistributedGraphInventory;
use super::{GraphEnthalpyTransport, TransientTransportError, thermal_checkpoint};
use crate::graph::inertial::{
    InertialGraph, InertialGraphControls, InertialGraphError, InertialGraphReport,
    InertialGraphTrial,
};
use crate::graph::{FixedPressure, LossGraph};

/// Accepted momentum plus an explicit node-to-inventory map. Unmapped solid
/// cells may share the original thermal residual. Every fluid node has distinct
/// fixed inventory at the momentum owner's reference density.
pub struct InertialThermalProblem<'a> {
    pub momentum: &'a InertialGraph,
    pub node_cells: &'a [CellId],
    /// Interior plus external mass connections at every actual proposal.
    pub maximum_connections: usize,
}

/// Source-owned momentum constitutive response at an actual trial enthalpy
/// field and candidate average branch rates. The retained owner admits topology,
/// branch identities, gravity geometry and reference density independently.
pub struct ResolvedInertialGraph<P> {
    pub graph: LossGraph,
    pub pressure_boundaries: Vec<FixedPressure>,
    /// Signed additional pressure head in original branch order.
    pub driving_pressure_pa: Vec<f64>,
    /// Original source owners/evaluations for the actual-flow admission callback.
    pub payload: P,
}

/// Heat exchange and external streams at the ACTUAL native average flow.
/// Internal advection is supplied by this adapter exactly once.
pub struct InertialGraphHeat<T> {
    pub mass_boundaries: Vec<MassBoundary>,
    pub heat: BoundaryStepCoefficients<T>,
}

/// Temporal mass-flux quadrature used by the original native thermal residual.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InertialAdvectionRule {
    /// Original midpoint-rate upwinding. Refine to locate interval reversals.
    SignedMean,
    /// Integrate forward/reverse swept volume of the native affine trajectory
    /// separately, with each original donor enthalpy held at the endpoint.
    AffineSwept,
}

/// Bounded joint nonlinear work and independent physical acceptance thresholds.
#[derive(Clone, Copy, Debug)]
pub struct InertialThermalControls {
    pub solver: NewtonKrylovConfig,
    /// All active thermal cells, including solid walls.
    pub maximum_cells: usize,
    pub maximum_newton_iterations: usize,
    /// Includes every derivative probe, final proposal and actual-flow source refresh.
    pub maximum_endpoint_evaluations: usize,
    pub energy_scale_j: f64,
    pub flow_scale_m3_s: f64,
    pub finite_difference_step: f64,
    pub cell_tolerance_j: f64,
    pub total_thermal_tolerance_j: f64,
    pub impulse_tolerance_pa_s: f64,
    pub mechanical_tolerance_j: f64,
    /// Constitutive consistency also needs a flow gate: impulse alone cannot
    /// constrain a candidate rate on an effectively lossless branch.
    pub average_flow_tolerance_m3_s: f64,
    pub momentum: InertialGraphControls,
}

/// Accepted sources, exchange and both native conservation accounts from the
/// same proposal. Mechanical work is not silently deposited as sensible heat.
#[derive(Debug)]
pub struct InertialThermalReport<P, T> {
    pub advection_rule: InertialAdvectionRule,
    pub hydraulics: InertialGraphReport,
    pub sources: P,
    pub heat: T,
    pub thermal: EnergyAudit,
    pub maximum_cell_residual_j: f64,
    pub endpoint_evaluations: usize,
    pub nonlinear: NewtonReport,
}

#[derive(Debug)]
pub enum InertialThermalError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Transport(TransientTransportError),
    Momentum(InertialGraphError),
    Solver(NewtonError),
    Resolve(E),
    Exchange(E),
    WorkBudget { performed: usize, maximum: usize },
    FlowResidual { value_m3_s: f64, maximum_m3_s: f64 },
    NotConverged(NewtonReport),
}
impl<E> From<DepositionError> for InertialThermalError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<TransientTransportError> for InertialThermalError<E> {
    fn from(error: TransientTransportError) -> Self {
        Self::Transport(error)
    }
}
impl<E> From<InertialGraphError> for InertialThermalError<E> {
    fn from(error: InertialGraphError) -> Self {
        Self::Momentum(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for InertialThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "inertial graph thermal coupling: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for InertialThermalError<E> {}

struct Prepared<P, T> {
    momentum: InertialGraphTrial,
    sources: P,
    heat: T,
}

/// Advance a branched fluid network and its original thermal cells atomically.
///
/// `resolve` is pure: each call compiles the actual candidate material field
/// and candidate average rates into native loss/gravity coefficients. `exchange`
/// receives the resulting ACTUAL average and old/new speed envelope, along with
/// those retained original sources. It must re-admit their material/Reynolds/
/// correlation support at these rates before returning wall exchange and any
/// explicit reservoir streams. Source refusals are terminal, including during
/// numerical derivative probes; there is no extrapolation or zero-flow fallback.
///
/// Cell-energy guesses follow native stable CellId order; flow guesses follow
/// original graph branch order. The solver closes both energy increments and
/// candidate-versus-native average rates. Final acceptance independently gates
/// physical joule, impulse and rate residuals. It resolves sources again at the
/// actual native average and rechecks the SAME momentum trajectory against that
/// law; a small rate mismatch cannot hide a sharp constitutive change. Heat is
/// reevaluated at the same endpoint enthalpies with those refreshed sources.
///
/// This combines midpoint momentum with endpoint-held backward-Euler thermal
/// coefficients, so the coupled method is first order in time. It does not
/// resolve the instant of a within-step reversal or create pipe control volumes.
/// Refine time and spatial inventory when residence times or reversal timing
/// matter. Native thermal events retain their original interval-splitting rules.
///
/// Both input owners are immutable. The returned momentum and thermal states
/// belong to the same accepted time; publish them together.
#[allow(clippy::too_many_arguments)]
pub fn advance_inertial_graph<'m, P, T, E>(
    cx: &Cx<'_>,
    problem: &InertialThermalProblem<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: InertialThermalControls,
    energy_guess_j: &[f64],
    average_flow_guess_m3_s: &[f64],
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
        InertialThermalReport<P, T>,
    ),
    InertialThermalError<E>,
> {
    advance_inertial_graph_with_inventory(
        cx,
        problem,
        None,
        InertialAdvectionRule::SignedMean,
        thermal,
        target_s,
        controls,
        energy_guess_j,
        average_flow_guess_m3_s,
        resolve,
        exchange,
    )
}

/// Advance original branch momentum with explicitly resolved pipe storage.
///
/// The inventory expands only native thermal transfers; the pressure graph and
/// its impulse coordinates remain unchanged. Each authored pipe cell keeps its
/// own enthalpy through startup, shutdown and reversal. The inventory's geometry
/// must imply the same branch inertance as the retained momentum owner. Native
/// energy and source admission are repeated at every candidate and final refresh.
/// Node cells, reference density and connection limits must match `problem`.
///
/// The time order, source callback contract and independent physical gates are
/// identical to [`advance_inertial_graph`]. Callers choose the spatial inventory;
/// the adapter never creates hidden hydraulic nodes or estimates a thermal mesh.
#[allow(clippy::too_many_arguments)]
pub fn advance_inertial_graph_distributed<'m, P, T, E>(
    cx: &Cx<'_>,
    problem: &InertialThermalProblem<'_>,
    inventory: &DistributedGraphInventory,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: InertialThermalControls,
    energy_guess_j: &[f64],
    average_flow_guess_m3_s: &[f64],
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
        InertialThermalReport<P, T>,
    ),
    InertialThermalError<E>,
> {
    advance_inertial_graph_with_inventory(
        cx,
        problem,
        Some(inventory),
        InertialAdvectionRule::SignedMean,
        thermal,
        target_s,
        controls,
        energy_guess_j,
        average_flow_guess_m3_s,
        resolve,
        exchange,
    )
}

/// Advance retained momentum with both swept directions at a flow reversal.
///
/// This calls the SAME joint native momentum/enthalpy residual and acceptance
/// gates as [`advance_inertial_graph`]. `inventory` optionally binds the same
/// original pipe cells as [`advance_inertial_graph_distributed`]. Every native
/// affine trajectory contributes its exact positive and negative swept volume
/// to separate native upwind transfers. A zero signed mean can therefore still
/// exchange heat between different original donor states.
///
/// `exchange` must supply both directional external mass streams. Use
/// [`super::affine::node_mass_streams`] on the report's initial/average/final
/// rates to obtain their required rates, and author each inlet's own original
/// enthalpy. Supplying only net rates across an external reversal refuses.
/// Existing source/correlation checks still inspect the complete flow envelope;
/// this quadrature does not authorize a law through unsupported stagnation.
///
/// This integrates a held-endpoint thermal flux, not the time ordering of
/// mixing or of boundary-layer evolution. The coupled method remains first
/// order. Additional reverse transfers consume the existing connection bound.
#[allow(clippy::too_many_arguments)]
pub fn advance_inertial_graph_affine<'m, P, T, E>(
    cx: &Cx<'_>,
    problem: &InertialThermalProblem<'_>,
    inventory: Option<&DistributedGraphInventory>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: InertialThermalControls,
    energy_guess_j: &[f64],
    average_flow_guess_m3_s: &[f64],
    resolve: impl Fn(&Cx<'_>, &DepositionSimulation<'m>, &[f64]) -> Result<ResolvedInertialGraph<P>, E>,
    exchange: impl Fn(
        &Cx<'_>,
        &DepositionSimulation<'m>,
        &P,
        &InertialGraphReport,
    ) -> Result<InertialGraphHeat<T>, E>,
) -> Result<
    (InertialGraph, DepositionSimulation<'m>, InertialThermalReport<P, T>),
    InertialThermalError<E>,
> {
    advance_inertial_graph_with_inventory(cx, problem, inventory,
        InertialAdvectionRule::AffineSwept, thermal, target_s, controls,
        energy_guess_j, average_flow_guess_m3_s, resolve, exchange)
}

#[allow(clippy::too_many_arguments)]
fn advance_inertial_graph_with_inventory<'m, P, T, E>(
    cx: &Cx<'_>,
    problem: &InertialThermalProblem<'_>,
    inventory: Option<&DistributedGraphInventory>,
    advection_rule: InertialAdvectionRule,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: InertialThermalControls,
    energy_guess_j: &[f64],
    average_flow_guess_m3_s: &[f64],
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
        InertialThermalReport<P, T>,
    ),
    InertialThermalError<E>,
> {
    thermal_checkpoint(cx)?;
    controls
        .solver
        .validate()
        .map_err(InertialThermalError::Solver)?;
    let unique: BTreeSet<_> = problem.node_cells.iter().copied().collect();
    let branches = problem.momentum.branch_count();
    if thermal.time_s() != problem.momentum.accepted_time_s()
        || !(1..=1024).contains(&controls.maximum_cells)
        || thermal.cells().len() > controls.maximum_cells
        || problem.node_cells.len() != problem.momentum.node_count()
        || unique.len() != problem.node_cells.len()
        || problem
            .node_cells
            .iter()
            .any(|id| !thermal.cells().contains_key(id))
        || branches > problem.maximum_connections
        || !(1..=1024).contains(&controls.maximum_newton_iterations)
        || !(1..=1_000_000).contains(&controls.maximum_endpoint_evaluations)
        || ![
            controls.energy_scale_j,
            controls.flow_scale_m3_s,
            controls.finite_difference_step,
        ]
        .into_iter()
        .all(positive)
        || ![
            controls.cell_tolerance_j,
            controls.total_thermal_tolerance_j,
            controls.impulse_tolerance_pa_s,
            controls.mechanical_tolerance_j,
            controls.average_flow_tolerance_m3_s,
        ]
        .into_iter()
        .all(nonnegative)
    {
        return Err(InertialThermalError::Invalid(
            "invalid synchronized timeline, inventory, scales or work limits",
        ));
    }
    if let Some(inventory) = inventory {
        if inventory.node_cells() != problem.node_cells
            || inventory.reference_density() != problem.momentum.reference_density()
            || inventory.interior_connection_count() > problem.maximum_connections
        {
            return Err(InertialThermalError::Invalid(
                "distributed inventory differs from the original node, density or connection binding",
            ));
        }
        inventory.admit_momentum(cx, problem.momentum)?;
    }
    let step = BackwardEulerStep::new(thermal, target_s, || thermal_checkpoint(cx))?;
    let count = step.cell_ids().len();
    let dimension = count
        .checked_add(branches)
        .ok_or(InertialThermalError::Invalid(
            "joint coordinate count overflow",
        ))?;
    if energy_guess_j.len() != count
        || average_flow_guess_m3_s.len() != branches
        || energy_guess_j
            .iter()
            .chain(average_flow_guess_m3_s)
            .any(|v| !v.is_finite())
    {
        return Err(InertialThermalError::Invalid(
            "complete finite thermal and branch-flow guesses required",
        ));
    }
    let evaluations = Cell::new(0usize);
    let evaluate = |x: &[f64], retained: Option<InertialGraphTrial>| {
        thermal_checkpoint(cx)?;
        if evaluations.get() >= controls.maximum_endpoint_evaluations {
            return Err(InertialThermalError::WorkBudget {
                performed: evaluations.get(),
                maximum: controls.maximum_endpoint_evaluations,
            });
        }
        evaluations.set(evaluations.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(InertialThermalError::Invalid(
                "invalid joint residual coordinates",
            ));
        }
        let energies: Vec<_> = x[..count]
            .iter()
            .map(|v| v * controls.energy_scale_j)
            .collect();
        let energies = step.resolvable_increments_j(&energies);
        let average: Vec<_> = x[count..]
            .iter()
            .map(|v| v * controls.flow_scale_m3_s)
            .collect();
        if energies.iter().chain(&average).any(|v| !v.is_finite()) {
            return Err(InertialThermalError::Invalid(
                "unrepresentable physical coordinates",
            ));
        }
        step.evaluate_with_transport(
            &energies,
            |candidate| {
                thermal_checkpoint(cx)?;
                let source_rates = retained.as_ref().map_or(average.as_slice(), |trial| {
                    trial.report().average_volume_flows_m3_s.as_slice()
                });
                let response =
                    resolve(cx, candidate, source_rates).map_err(InertialThermalError::Resolve)?;
                thermal_checkpoint(cx)?;
                let momentum = match retained {
                    Some(trial) => trial.recheck(
                        cx,
                        &response.graph,
                        &response.pressure_boundaries,
                        &response.driving_pressure_pa,
                        controls.momentum,
                    )?,
                    None => problem.momentum.propose(
                        cx,
                        &response.graph,
                        &response.pressure_boundaries,
                        &response.driving_pressure_pa,
                        &average,
                        target_s,
                        controls.momentum,
                    )?,
                };
                let coefficients = exchange(cx, candidate, &response.payload, momentum.report())
                    .map_err(InertialThermalError::Exchange)?;
                thermal_checkpoint(cx)?;
                let binding = match (inventory, advection_rule) {
                    (Some(inventory), InertialAdvectionRule::AffineSwept) => inventory.bind_affine(
                        cx, &response.graph,
                        &momentum.report().initial_flows_m3_s,
                        &momentum.report().average_volume_flows_m3_s,
                        &momentum.report().final_flows_m3_s,
                        candidate, coefficients.mass_boundaries, problem.maximum_connections,
                    )?,
                    (None, InertialAdvectionRule::AffineSwept) => GraphEnthalpyTransport::from_affine_volume_flows(
                        cx, &response.graph,
                        &momentum.report().initial_flows_m3_s,
                        &momentum.report().average_volume_flows_m3_s,
                        &momentum.report().final_flows_m3_s,
                        problem.momentum.reference_density(), problem.node_cells,
                        candidate, coefficients.mass_boundaries, problem.maximum_connections,
                    )?,
                    (Some(inventory), InertialAdvectionRule::SignedMean) => inventory.bind(
                        cx,
                        &response.graph,
                        &momentum.report().average_volume_flows_m3_s,
                        candidate,
                        coefficients.mass_boundaries,
                        problem.maximum_connections,
                    )?,
                    (None, InertialAdvectionRule::SignedMean) => GraphEnthalpyTransport::from_volume_flows(
                        cx,
                        &response.graph,
                        &momentum.report().average_volume_flows_m3_s,
                        problem.momentum.reference_density(),
                        problem.node_cells,
                        candidate,
                        coefficients.mass_boundaries,
                        problem.maximum_connections,
                    )?,
                };
                let advected = binding.endpoint_transport(cx, candidate)?;
                let BoundaryStepCoefficients { bulk, boundaries } = coefficients.heat;
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = bulk;
                for (id, watts) in advected.heat.bulk.sources_w {
                    thermal_checkpoint(cx)?;
                    let value = sources_w.entry(id).or_insert(0.0);
                    *value += watts;
                    if !value.is_finite() {
                        return Err(InertialThermalError::Invalid(
                            "combined endpoint source is nonfinite",
                        ));
                    }
                }
                Ok(TransportStepCoefficients {
                    heat: BoundaryStepCoefficients {
                        bulk: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: Prepared {
                                momentum,
                                sources: response.payload,
                                heat: payload,
                            },
                        },
                        boundaries,
                    },
                    transfers: advected.transfers,
                })
            },
            || thermal_checkpoint(cx),
        )
    };
    let difference =
        FallibleFiniteDifference::new(dimension, controls.finite_difference_step, |x: &[f64]| {
            let trial = evaluate(x, None)?;
            let mut residual: Vec<_> = trial
                .residual_j()
                .iter()
                .map(|value| value / controls.energy_scale_j)
                .collect();
            residual.extend(
                trial
                    .payload()
                    .momentum
                    .report()
                    .average_flow_residual_m3_s
                    .iter()
                    .map(|value| value / controls.flow_scale_m3_s),
            );
            Ok(residual)
        })
        .map_err(difference_error)?;
    let mut guess: Vec<_> = energy_guess_j
        .iter()
        .map(|v| v / controls.energy_scale_j)
        .collect();
    guess.extend(
        average_flow_guess_m3_s
            .iter()
            .map(|v| v / controls.flow_scale_m3_s),
    );
    let initialized = NewtonKrylovState::new(&difference, guess, controls.solver);
    if let Some(error) = difference.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = initialized.map_err(InertialThermalError::Solver)?;
    let run = solver.run_controlled(&difference, controls.maximum_newton_iterations, &mut || {
        difference.has_failed() || cx.checkpoint().is_err()
    });
    if let Some(error) = difference.take_failure() {
        return Err(difference_error(error));
    }
    thermal_checkpoint(cx)?;
    if !run.report.converged {
        return Err(InertialThermalError::NotConverged(run.report));
    }
    let final_proposal = evaluate(&solver.x, None)?;
    let retained = final_proposal.payload().momentum.clone();
    drop(final_proposal);
    // Retain the trajectory, not its candidate-flow constitutive payload. This
    // additional endpoint call is budgeted and reruns heat with fresh sources.
    let final_trial = evaluate(&solver.x, Some(retained))?;
    let flow_residual = final_trial
        .payload()
        .momentum
        .report()
        .average_flow_residual_m3_s
        .iter()
        .map(|value| value.abs())
        .fold(0.0f64, f64::max);
    if !flow_residual.is_finite() || flow_residual > controls.average_flow_tolerance_m3_s {
        return Err(InertialThermalError::FlowResidual {
            value_m3_s: flow_residual,
            maximum_m3_s: controls.average_flow_tolerance_m3_s,
        });
    }
    let maximum_cell_residual_j = final_trial.maximum_cell_residual_j();
    let (next_thermal, prepared, thermal_audit) = final_trial.accept(
        controls.cell_tolerance_j,
        controls.total_thermal_tolerance_j,
        || thermal_checkpoint(cx),
    )?;
    let (next_momentum, hydraulics) = prepared.momentum.accept(
        cx,
        controls.impulse_tolerance_pa_s,
        controls.mechanical_tolerance_j,
    )?;
    thermal_checkpoint(cx)?;
    Ok((
        next_momentum,
        next_thermal,
        InertialThermalReport {
            advection_rule,
            hydraulics,
            sources: prepared.sources,
            heat: prepared.heat,
            thermal: thermal_audit,
            maximum_cell_residual_j,
            endpoint_evaluations: evaluations.get(),
            nonlinear: run.report,
        },
    ))
}

fn difference_error<E>(
    error: FiniteDifferenceError<InertialThermalError<E>>,
) -> InertialThermalError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => InertialThermalError::Invalid(message),
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

#[cfg(test)]
mod tests;
