//! Quasi-steady hydraulic feedback on the original thermal endpoint residual.
//!
//! Each actual trial enthalpy field resolves its own hydraulic graph and wall
//! exchange. The native backward-Euler owner retains storage, caloric inversion,
//! sources, conservative advection, event boundaries and final energy acceptance.
//! Newton--Krylov solves only cell energy increments; no lagged flow or separate
//! thermal clock can be published. Momentum inertia is explicitly omitted.

use std::cell::Cell;
use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_qty::Density;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use fs_thermal::enthalpy::advection::MassBoundary;
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};

use super::{GraphEnthalpyTransport, TransientTransportError, thermal_checkpoint};
use crate::graph::{GraphSolution, LossGraph};

/// Fixed topology and conserved reference inventory. Thermal wall cells may
/// coexist with these fluid cells; every hydraulic node needs a distinct cell.
pub struct QuasiSteadyGraphProblem<'a> {
    pub topology: &'a LossGraph,
    pub node_cells: &'a [CellId],
    pub reference_density: Density,
    /// Total interior plus external mass connections at any endpoint.
    pub maximum_connections: usize,
}

/// Complete source-owned endpoint response. The callback must solve hydraulics
/// at THIS thermal state, retain the declared topology/branch identities, and
/// perform its own material, geometry and correlation applicability checks.
pub struct QuasiSteadyGraphResponse<T> {
    pub hydraulics: GraphSolution,
    /// Explicit external streams, never inferred from pressure datums.
    pub mass_boundaries: Vec<MassBoundary>,
    /// Nonadvective source/link/boundary coefficients only. Authored watts and
    /// ambient temperatures retain their native meaning; advective watts are
    /// supplied exactly once by this adapter.
    pub heat: BoundaryStepCoefficients<T>,
}

/// Numerical work, physical residuals and low-Eckert approximation limits.
#[derive(Debug, Clone, Copy)]
pub struct QuasiSteadyGraphControls {
    pub solver: NewtonKrylovConfig,
    /// Includes unmapped solid/wall cells in the same thermal residual.
    pub maximum_cells: usize,
    pub maximum_newton_iterations: usize,
    /// Includes all derivative probes and the independent final evaluation.
    pub maximum_endpoint_evaluations: usize,
    /// Positive coordinate and residual scale, joules per unknown.
    pub energy_scale_j: f64,
    /// Maximum dimensionless coordinate perturbation of a Jacobian action.
    pub finite_difference_step: f64,
    pub cell_tolerance_j: f64,
    pub total_tolerance_j: f64,
    /// Sum of dissipative hydraulic power and absolute local anomalous
    /// buoyancy work. Neither is silently deposited into sensible heat.
    pub maximum_neglected_power_w: f64,
    pub maximum_neglected_energy_j: f64,
}

/// Accepted hydraulics, heat payload and thermal audit from the SAME endpoint.
#[derive(Debug)]
pub struct QuasiSteadyGraphReport<T> {
    pub hydraulics: GraphSolution,
    pub heat: T,
    pub thermal: EnergyAudit,
    pub maximum_cell_residual_j: f64,
    pub endpoint_evaluations: usize,
    pub nonlinear: NewtonReport,
    pub neglected_power_w: f64,
    pub neglected_energy_j: f64,
}

#[derive(Debug)]
pub enum QuasiSteadyGraphError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Transport(TransientTransportError),
    Solver(NewtonError),
    Resolve(E),
    WorkBudget { performed: usize, maximum: usize },
    NeglectedWork { power_w: f64, energy_j: f64 },
    NotConverged(NewtonReport),
}
impl<E> From<DepositionError> for QuasiSteadyGraphError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<TransientTransportError> for QuasiSteadyGraphError<E> {
    fn from(error: TransientTransportError) -> Self {
        Self::Transport(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for QuasiSteadyGraphError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quasi-steady graph thermal coupling: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for QuasiSteadyGraphError<E> {}

struct Prepared<T> {
    hydraulics: GraphSolution,
    heat: T,
    neglected_power_w: f64,
}

/// Advance one algebraically coupled thermal/flow step without mutating the
/// input timeline or any accepted hydraulic history.
///
/// The pure callback is evaluated at actual native endpoint enthalpies, including
/// numerical Jacobian probes. It may query source materials and solve an updated
/// native graph, then supply flow-dependent native wall conductances. Source
/// refusals are terminal; they never become zero residuals or fallback laws.
/// Final acceptance re-evaluates hydraulics/heat at the solved energy vector
/// and gates the original per-cell and total energy residual independently.
///
/// This is the quasi-steady momentum approximation: flows may reverse but their
/// endpoint values represent the whole implicit interval. It does not retain
/// inertial momentum or locate the instant of reversal. Refine time and fluid
/// control volumes; use the inertial owner when hydraulic relaxation matters.
#[allow(clippy::too_many_arguments)]
pub fn advance_quasi_steady_graph<'m, T, E>(
    cx: &Cx<'_>,
    problem: &QuasiSteadyGraphProblem<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: QuasiSteadyGraphControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&Cx<'_>, &DepositionSimulation<'m>) -> Result<QuasiSteadyGraphResponse<T>, E>,
) -> Result<(DepositionSimulation<'m>, QuasiSteadyGraphReport<T>), QuasiSteadyGraphError<E>> {
    thermal_checkpoint(cx)?;
    controls
        .solver
        .validate()
        .map_err(QuasiSteadyGraphError::Solver)?;
    let rho = problem.reference_density.value();
    let unique: BTreeSet<_> = problem.node_cells.iter().copied().collect();
    if !(1..=1024).contains(&controls.maximum_cells)
        || thermal.cells().len() > controls.maximum_cells
        || problem.node_cells.len() != problem.topology.node_count()
        || unique.len() != problem.node_cells.len()
        || problem
            .node_cells
            .iter()
            .any(|id| !thermal.cells().contains_key(id))
        || !positive(rho)
        || problem.topology.branches().len() > problem.maximum_connections
        || !(1..=1024).contains(&controls.maximum_newton_iterations)
        || !(1..=1_000_000).contains(&controls.maximum_endpoint_evaluations)
        || !positive(controls.energy_scale_j)
        || !positive(controls.finite_difference_step)
        || ![
            controls.cell_tolerance_j,
            controls.total_tolerance_j,
            controls.maximum_neglected_power_w,
            controls.maximum_neglected_energy_j,
        ]
        .into_iter()
        .all(nonnegative)
    {
        return Err(QuasiSteadyGraphError::Invalid(
            "invalid topology, scale, residual or work budget",
        ));
    }
    let step = BackwardEulerStep::new(thermal, target_s, || thermal_checkpoint(cx))?;
    let count = step.cell_ids().len();
    if energy_guess_j.len() != count || energy_guess_j.iter().any(|v| !v.is_finite()) {
        return Err(QuasiSteadyGraphError::Invalid(
            "complete finite cell-energy guess required",
        ));
    }
    let evaluations = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        thermal_checkpoint(cx)?;
        if evaluations.get() == controls.maximum_endpoint_evaluations {
            return Err(QuasiSteadyGraphError::WorkBudget {
                performed: evaluations.get(),
                maximum: controls.maximum_endpoint_evaluations,
            });
        }
        evaluations.set(evaluations.get() + 1);
        let energies: Vec<_> = x.iter().map(|v| v * controls.energy_scale_j).collect();
        let energies = step.resolvable_increments_j(&energies);
        step.evaluate_with_transport(
            &energies,
            |candidate| {
                thermal_checkpoint(cx)?;
                let response = resolve(cx, candidate).map_err(QuasiSteadyGraphError::Resolve)?;
                thermal_checkpoint(cx)?;
                admit_topology(problem.topology, &response.hydraulics)?;
                let binding = GraphEnthalpyTransport::new(
                    cx,
                    &response.hydraulics,
                    problem.reference_density,
                    problem.node_cells,
                    candidate,
                    response.mass_boundaries,
                    problem.maximum_connections,
                )?;
                let neglected = neglected_power(cx, &response.hydraulics)?;
                let neglected_energy = neglected * step.duration_s();
                if !neglected_energy.is_finite()
                    || neglected > controls.maximum_neglected_power_w
                    || neglected_energy > controls.maximum_neglected_energy_j
                {
                    return Err(QuasiSteadyGraphError::NeglectedWork {
                        power_w: neglected,
                        energy_j: neglected_energy,
                    });
                }
                let advected = binding.endpoint_transport(cx, candidate)?;
                let BoundaryStepCoefficients { bulk, boundaries } = response.heat;
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
                        return Err(QuasiSteadyGraphError::Invalid(
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
                                hydraulics: response.hydraulics,
                                heat: payload,
                                neglected_power_w: neglected,
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
        FallibleFiniteDifference::new(count, controls.finite_difference_step, |x: &[f64]| {
            let trial = evaluate(x)?;
            Ok(trial
                .residual_j()
                .iter()
                .map(|v| v / controls.energy_scale_j)
                .collect())
        })
        .map_err(difference_error)?;
    let guess = energy_guess_j
        .iter()
        .map(|v| v / controls.energy_scale_j)
        .collect();
    let initialized = NewtonKrylovState::new(&difference, guess, controls.solver);
    if let Some(error) = difference.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = initialized.map_err(QuasiSteadyGraphError::Solver)?;
    let run = solver.run_controlled(&difference, controls.maximum_newton_iterations, &mut || {
        difference.has_failed() || cx.checkpoint().is_err()
    });
    if let Some(error) = difference.take_failure() {
        return Err(difference_error(error));
    }
    thermal_checkpoint(cx)?;
    if !run.report.converged {
        return Err(QuasiSteadyGraphError::NotConverged(run.report));
    }
    let final_trial = evaluate(&solver.x)?;
    let maximum_cell_residual_j = final_trial.maximum_cell_residual_j();
    let (next, prepared, audit) = final_trial.accept(
        controls.cell_tolerance_j,
        controls.total_tolerance_j,
        || thermal_checkpoint(cx),
    )?;
    thermal_checkpoint(cx)?;
    Ok((
        next,
        QuasiSteadyGraphReport {
            hydraulics: prepared.hydraulics,
            heat: prepared.heat,
            thermal: audit,
            maximum_cell_residual_j,
            endpoint_evaluations: evaluations.get(),
            nonlinear: run.report,
            neglected_power_w: prepared.neglected_power_w,
            neglected_energy_j: prepared.neglected_power_w * step.duration_s(),
        },
    ))
}

fn admit_topology<E>(
    graph: &LossGraph,
    flow: &GraphSolution,
) -> Result<(), QuasiSteadyGraphError<E>> {
    if flow.pressures.len() != graph.node_count()
        || flow.branches.len() != graph.branches().len()
        || flow.pressures.iter().any(|p| !p.value().is_finite())
        || flow
            .branches
            .iter()
            .zip(graph.branches())
            .any(|(a, b)| a.from != b.from || a.to != b.to || a.loss.name != b.loss.name)
    {
        return Err(QuasiSteadyGraphError::Invalid(
            "endpoint hydraulic topology or identities changed",
        ));
    }
    Ok(())
}

fn neglected_power<E>(cx: &Cx<'_>, flow: &GraphSolution) -> Result<f64, QuasiSteadyGraphError<E>> {
    let mut total = 0.0;
    if let Some(gravity) = &flow.gravity {
        if gravity.branch_materials().len() != flow.branches.len()
            || gravity.positions_m().len() != flow.pressures.len()
        {
            return Err(QuasiSteadyGraphError::Invalid(
                "gravity and hydraulic topology differ",
            ));
        }
    }
    for (index, branch) in flow.branches.iter().enumerate() {
        thermal_checkpoint(cx)?;
        let dissipation = branch.flow.value() * branch.dissipative_pressure_drop.value();
        if !nonnegative(dissipation) {
            return Err(QuasiSteadyGraphError::Invalid(
                "hydraulic loss power is invalid",
            ));
        }
        total += dissipation;
        if let Some(gravity) = &flow.gravity {
            let from = gravity.positions_m()[branch.from];
            let to = gravity.positions_m()[branch.to];
            let acceleration = gravity.acceleration_m_s2();
            let height = (0..3)
                .map(|axis| acceleration[axis] * (to[axis] - from[axis]))
                .sum::<f64>();
            let departure = gravity.branch_materials()[index].density_kg_m3()
                - gravity.reference_transport().density_kg_m3();
            total += (branch.flow.value() * departure * height).abs();
        }
        if !total.is_finite() {
            return Err(QuasiSteadyGraphError::Invalid(
                "neglected hydraulic work is unrepresentable",
            ));
        }
    }
    Ok(total)
}

fn difference_error<E>(
    error: FiniteDifferenceError<QuasiSteadyGraphError<E>>,
) -> QuasiSteadyGraphError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => QuasiSteadyGraphError::Invalid(message),
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
