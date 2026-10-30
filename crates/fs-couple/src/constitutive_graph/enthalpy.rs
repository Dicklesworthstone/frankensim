//! Material-law feedback on the original finite thermal timeline.
//!
//! Every thermal trial solves the public material equations at THAT endpoint.
//! Graph laws retain their own consistent tangents and committed history;
//! `BackwardEulerStep` retains mass, caloric/latent storage, conduction, physical
//! boundaries and events. An outer Newton solve closes this partitioned residual.
//! Its finite differences include the converged inner material response, rather
//! than freezing a current, heat source or material coefficient from the past.
//!
//! Heat is routed only from explicitly bound watt-valued graph outputs. Generic
//! graph dissipation reports may be densities and are never silently interpreted
//! as watts. The caller owns the physical meaning and completeness of these heat
//! ports and the work/storage account of the other domain. This adapter establishes
//! the original thermal balances, not passivity of arbitrary caller-authored laws.

use super::{
    MaterialEquation, MaterialGraphCandidate, MaterialGraphError, MaterialGraphResidual,
    MaterialGraphSolve, MaterialUnknown, PortKey,
};
use fs_exec::Cx;
use fs_material::graph::ConstitutiveGraph;
use fs_qty::Dims;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// The graph input's explicitly declared temperature coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermalCoordinate {
    /// Absolute thermodynamic temperature in kelvin.
    Absolute,
    /// Kelvin difference from this fixed, positive reference temperature.
    ChangeFrom(f64),
}

/// A graph input supplied by an actual finite thermal cell at each trial.
#[derive(Clone, Debug)]
pub struct GraphTemperature {
    /// Original unconnected input with temperature dimensions.
    pub input: PortKey,
    /// Existing thermal cell, with its original material and mass.
    pub cell: CellId,
    /// Absolute temperature or an explicit difference coordinate.
    pub coordinate: ThermalCoordinate,
}

/// An original graph heat output and its unique receiving thermal cell.
#[derive(Clone, Debug)]
pub struct GraphHeatPort {
    /// Original output, with power dimensions. Signed heat is permitted.
    pub output: PortKey,
    /// Existing cell. Different outputs may feed the same physical cell.
    pub cell: CellId,
}

/// Immutable graph history and equation data for one finite physical interval.
/// Every inner solve starts with the same committed history and initial guess.
pub struct GraphThermalInput<'a> {
    /// Original admitted law graph; no duplicate material registry.
    pub graph: &'a ConstitutiveGraph,
    /// Original aggregate state schema version.
    pub version: u64,
    /// Committed graph history, never mutated during any thermal trial.
    pub state: &'a [f64],
    /// Held external inputs, excluding temperature-bound and solved inputs.
    pub fixed: &'a BTreeMap<PortKey, f64>,
    /// Original square material feedback block.
    pub unknowns: &'a [MaterialUnknown],
    /// Dimensioned balances, evaluated at the same thermal endpoint.
    pub equations: &'a [MaterialEquation],
    /// Normalized initial coordinates for every inner solve.
    pub initial_coordinates: &'a [f64],
    /// Unique external inputs driven by the original caloric temperatures.
    pub temperatures: &'a [GraphTemperature],
    /// Unique seconds-valued external inputs supplied with the ORIGINAL physical
    /// interval duration. History laws therefore see one accepted thermal clock.
    pub durations: &'a [PortKey],
    /// Unique original heat outputs; no hidden conversion or duplicate deposit.
    pub heat: &'a [GraphHeatPort],
}

/// Bounded nested nonlinear work and independent physical energy admission.
#[derive(Clone, Copy, Debug)]
pub struct GraphThermalSolve {
    /// Maximum active original thermal cells.
    pub maximum_cells: usize,
    /// Maximum complete endpoint trials, including stencils and final admission.
    pub maximum_thermal_trials: usize,
    /// Total graph value/Jacobian budget across ALL inner solves.
    pub maximum_graph_evaluations: usize,
    /// Outer Newton iteration allowance.
    pub maximum_iterations: usize,
    /// Outer numerical derivative perturbation in cell-energy coordinates [J].
    pub difference_step_j: f64,
    /// Independent maximum local energy defect [J].
    pub cell_tolerance_j: f64,
    /// Independent total thermal energy defect [J].
    pub total_tolerance_j: f64,
    /// Native outer Newton/Krylov controls, acting on joules.
    pub newton: NewtonKrylovConfig,
    /// Exact-tangent inner material controls, additionally capped by total work.
    pub graph: MaterialGraphSolve,
}

/// Original physical or numerical refusal; neither history has been advanced.
#[derive(Debug)]
pub enum GraphThermalError<E> {
    /// Invalid correspondence, dimensions, work controls or finite arithmetic.
    Invalid(&'static str),
    /// Original material graph/source or inner solve failure.
    Graph(MaterialGraphError),
    /// Original material heat/contact/radiation callback refusal.
    Source(E),
    /// Original caloric, topology, event or energy admission failure.
    Thermal(DepositionError),
    /// Native outer nonlinear construction/arithmetic failure.
    Newton(NewtonError),
    /// Outer solver stopped without convergence.
    NotConverged(NewtonReport),
    /// Outer iteration was interrupted before completion.
    Interrupted(SolverRunProgress),
    /// The complete endpoint-trial allowance was consumed.
    ThermalBudget,
    /// The aggregate inner material evaluation allowance was consumed.
    GraphBudget,
    /// Cancellation before returning the joint physical candidate.
    Cancelled,
}

impl<E> From<DepositionError> for GraphThermalError<E> {
    fn from(error: DepositionError) -> Self {
        match error {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl<E> From<MaterialGraphError> for GraphThermalError<E> {
    fn from(error: MaterialGraphError) -> Self {
        match error {
            MaterialGraphError::Cancelled => Self::Cancelled,
            other => Self::Graph(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for GraphThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material graph/enthalpy coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GraphThermalError<E> {}

/// Work and source observations of the SAME accepted endpoint.
pub struct GraphThermalReport<S> {
    /// Original thermal storage, boundaries and internal-transfer accounting.
    pub thermal: EnergyAudit,
    /// Exact final material conduction/contact/radiation callback payload.
    pub source: S,
    /// Native outer solver report; not a timestep-accuracy certificate.
    pub nonlinear: NewtonReport,
    /// Largest independently checked original cell-energy residual [J].
    pub maximum_cell_residual_j: f64,
    /// Signed graph heat delivered once to the receiving cells [J].
    pub graph_heat_j: f64,
    /// Independently evaluated authored and callback external heat [J].
    pub other_external_heat_j: f64,
    /// Net numerical heat of a composed closed radiation enclosure [J]. The
    /// original thermal audit retains it; external heat above excludes it.
    pub closed_radiation_net_j: f64,
    /// All endpoint trials, including discarded stencils and final admission.
    pub thermal_trials: usize,
    /// All consumed inner value/Jacobian calls, including final graph admission.
    pub graph_evaluations: usize,
}

/// An inseparable proposed thermal endpoint and material-law history.
pub struct GraphThermalCandidate<'m, S> {
    thermal: DepositionSimulation<'m>,
    graph: MaterialGraphCandidate,
    report: GraphThermalReport<S>,
}
impl<'m, S> GraphThermalCandidate<'m, S> {
    /// Original thermal endpoint; accepted-time events remain pending.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    /// Original graph outputs, solved inputs and proposed history.
    pub const fn graph(&self) -> &MaterialGraphCandidate {
        &self.graph
    }
    /// Physical thermal account and retained material observations.
    pub const fn report(&self) -> &GraphThermalReport<S> {
        &self.report
    }
    /// Return both histories together for publication by the physical owner.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        MaterialGraphCandidate,
        GraphThermalReport<S>,
    ) {
        (self.thermal, self.graph, self.report)
    }
}

/// Close one material/thermal feedback interval on the existing thermal graph.
///
/// The callback receives the same thermal endpoint as the material laws. It
/// supplies only ADDITIONAL sources and replacement original link conductances:
/// do not repeat the bound graph heat. Existing prepared material heat networks
/// can supply conductivity, contacts and radiation here without a second thermal
/// integrator. Every callback and inner solve must be pure relative to committed
/// physical state. Sources and supports must cover all derivative stencils.
///
/// Fixed geometry/topology and backward-Euler thermal rates. Process pending
/// start events and split at interior events using the original timeline. A
/// source discontinuity needs an explicit interval cut. No phase kinetics,
/// moving mass, arbitrary graph energy certificate or timestep order beyond the
/// original first-order thermal discretization is introduced.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: GraphThermalInput<'_>,
    controls: GraphThermalSolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<StepCoefficients<S>, E>,
) -> Result<GraphThermalCandidate<'m, S>, GraphThermalError<E>> {
    let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
    checkpoint()?;
    controls
        .newton
        .validate()
        .map_err(GraphThermalError::Newton)?;
    if thermal.cells().is_empty()
        || thermal.cells().len() > controls.maximum_cells
        || initial_energy_increments_j.len() != thermal.cells().len()
        || controls.maximum_thermal_trials == 0
        || controls.maximum_graph_evaluations == 0
        || controls.maximum_iterations == 0
        || !positive(controls.difference_step_j)
        || !nonnegative(controls.cell_tolerance_j)
        || !nonnegative(controls.total_tolerance_j)
        || input.unknowns.len() != input.initial_coordinates.len()
        || input.temperatures.len() > controls.graph.maximum_terms
        || input.durations.len() > controls.graph.maximum_terms
        || input.heat.len() > controls.graph.maximum_terms
    {
        return Err(GraphThermalError::Invalid(
            "invalid correspondence sizes or work/energy controls",
        ));
    }
    let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
    let mut temperature_keys = BTreeSet::new();
    for binding in input.temperatures {
        checkpoint()?;
        if !thermal.cells().contains_key(&binding.cell)
            || !temperature_keys.insert(binding.input.clone())
            || input.fixed.contains_key(&binding.input)
            || input.unknowns.iter().any(|u| u.input == binding.input)
            || input
                .graph
                .external_input_dims(&binding.input.0, &binding.input.1)
                != Some(TEMPERATURE)
            || matches!(binding.coordinate, ThermalCoordinate::ChangeFrom(t) if !positive(t))
        {
            return Err(GraphThermalError::Invalid(
                "temperature must uniquely bind an original cell to an external temperature input",
            ));
        }
    }
    let mut heat_keys = BTreeSet::new();
    for duration in input.durations {
        checkpoint()?;
        if !temperature_keys.insert(duration.clone())
            || input.fixed.contains_key(duration)
            || input.unknowns.iter().any(|u| &u.input == duration)
            || input.graph.external_input_dims(&duration.0, &duration.1) != Some(TIME)
        {
            return Err(GraphThermalError::Invalid(
                "duration must uniquely bind an external seconds-valued input",
            ));
        }
    }
    for binding in input.heat {
        checkpoint()?;
        if !thermal.cells().contains_key(&binding.cell)
            || !heat_keys.insert(binding.output.clone())
            || input
                .graph
                .output_dims(&binding.output.0, &binding.output.1)
                != Some(POWER)
        {
            return Err(GraphThermalError::Invalid(
                "heat must uniquely bind an original watt-valued output to an active cell",
            ));
        }
    }
    let trials = Cell::new(0usize);
    let graph_evaluations = Cell::new(0usize);
    let evaluate = |energy: &[f64], final_account: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_thermal_trials {
            return Err(GraphThermalError::ThermalBudget);
        }
        trials.set(trials.get() + 1);
        step.evaluate(
            energy,
            |endpoint| {
                let remaining = controls.maximum_graph_evaluations - graph_evaluations.get();
                if remaining == 0 {
                    return Err(GraphThermalError::GraphBudget);
                }
                let mut fixed = input.fixed.clone();
                for duration in input.durations {
                    checkpoint()?;
                    fixed.insert(duration.clone(), step.duration_s());
                }
                for binding in input.temperatures {
                    checkpoint()?;
                    let t = endpoint.temperature_k(binding.cell)?;
                    let value = match binding.coordinate {
                        ThermalCoordinate::Absolute => t,
                        ThermalCoordinate::ChangeFrom(reference) => t - reference,
                    };
                    if !value.is_finite() {
                        return Err(GraphThermalError::Invalid(
                            "nonfinite thermal input coordinate",
                        ));
                    }
                    fixed.insert(binding.input.clone(), value);
                }
                let mut graph_controls = controls.graph;
                graph_controls.maximum_evaluations =
                    graph_controls.maximum_evaluations.min(remaining);
                let problem = MaterialGraphResidual::new(
                    cx,
                    input.graph,
                    input.version,
                    input.state,
                    fixed,
                    input.unknowns.to_vec(),
                    input.equations.to_vec(),
                    graph_controls,
                )?;
                let graph = problem.solve(input.initial_coordinates.to_vec());
                graph_evaluations.set(graph_evaluations.get() + problem.evaluations());
                let graph = graph?;
                checkpoint()?;
                let mut coefficients = resolve(endpoint).map_err(GraphThermalError::Source)?;
                // The independent supply account belongs to fresh final admission,
                // not every discarded Newton stencil. All trials still use the
                // original thermal rate owner in BackwardEulerStep below.
                let other_power = if final_account {
                    endpoint.external_step_power_w(
                        &coefficients.sources_w,
                        &coefficients.conductances_w_k,
                        checkpoint,
                    )?
                } else {
                    0.0
                };
                let mut graph_power = 0.0;
                for binding in input.heat {
                    checkpoint()?;
                    let power = *graph.graph.outputs.get(&binding.output).ok_or(
                        GraphThermalError::Invalid("original heat output disappeared"),
                    )?;
                    graph_power += power;
                    let target = coefficients.sources_w.entry(binding.cell).or_insert(0.0);
                    *target += power;
                    if !target.is_finite() || !graph_power.is_finite() {
                        return Err(GraphThermalError::Invalid(
                            "unrepresentable combined thermal power",
                        ));
                    }
                }
                Ok(StepCoefficients {
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    payload: (graph, coefficients.payload, graph_power, other_power),
                })
            },
            checkpoint,
        )
    };
    let numerical = |error| match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => GraphThermalError::Invalid(reason),
    };
    let problem = FallibleFiniteDifference::new(
        step.cell_ids().len(),
        controls.difference_step_j,
        |x: &[f64]| evaluate(x, false).map(|trial| trial.residual_j().to_vec()),
    )
    .map_err(numerical)?;
    let constructed = NewtonKrylovState::new(
        &problem,
        initial_energy_increments_j.to_vec(),
        controls.newton,
    );
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = constructed.map_err(GraphThermalError::Newton)?;
    let result = solver.run_controlled(&problem, controls.maximum_iterations, &mut || {
        problem.has_failed() || cx.checkpoint().is_err()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(GraphThermalError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(GraphThermalError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, true)?;
    let maximum_cell_residual_j = trial.maximum_cell_residual_j();
    let (thermal, (graph, source, graph_power, other_power), audit) = trial.accept(
        controls.cell_tolerance_j,
        controls.total_tolerance_j,
        checkpoint,
    )?;
    let graph_heat_j = step.duration_s() * graph_power;
    let other_external_heat_j = step.duration_s() * other_power;
    if !graph_heat_j.is_finite()
        || !other_external_heat_j.is_finite()
        || (audit.external_heat_j - graph_heat_j - other_external_heat_j).abs()
            > controls.total_tolerance_j
    {
        return Err(GraphThermalError::Invalid(
            "graph heat and independently evaluated external heat disagree",
        ));
    }
    checkpoint()?;
    Ok(GraphThermalCandidate {
        thermal,
        graph,
        report: GraphThermalReport {
            thermal: audit,
            source,
            nonlinear: result.report,
            maximum_cell_residual_j,
            graph_heat_j,
            other_external_heat_j,
            closed_radiation_net_j: 0.0,
            thermal_trials: trials.get(),
            graph_evaluations: graph_evaluations.get(),
        },
    })
}

impl crate::circuit_enthalpy::contacts::heat::PreparedHeatNetwork<'_, '_> {
    /// Compose original material conduction, interfaces and radiation with
    /// public constitutive-graph feedback on one finite thermal timeline.
    ///
    /// Electrical circuit bindings must be empty: this graph owns its explicit
    /// electrical equations and heat ports, so circuit losses cannot be silently
    /// dropped or duplicated. Source-integrated conductivity attachments remain
    /// active. The same local/total energy tolerance also checks the thermal
    /// residual after removing closed-enclosure radiative roundoff from supply.
    #[allow(clippy::too_many_arguments)]
    pub fn graph_step_candidate<'m>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: GraphThermalInput<'_>,
        controls: GraphThermalSolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<
        GraphThermalCandidate<'m, crate::circuit_enthalpy::contacts::heat::HeatNetworkSample>,
        GraphThermalError<crate::circuit_enthalpy::contacts::heat::HeatNetworkError>,
    > {
        self.thermal_material_network()
            .map_err(GraphThermalError::Source)?;
        let mut candidate = step_candidate(
            cx,
            thermal,
            target_s,
            input,
            controls,
            initial_energy_increments_j,
            |endpoint| {
                let source = self.resolve(cx, endpoint)?.bulk;
                Ok(StepCoefficients {
                    sources_w: source.sources_w,
                    conductances_w_k: source.conductances_w_k,
                    payload: source.source,
                })
            },
        )
        .map_err(|error| match error {
            GraphThermalError::Source(
                crate::circuit_enthalpy::contacts::heat::HeatNetworkError::Cancelled,
            ) => GraphThermalError::Cancelled,
            other => other,
        })?;
        let net = (target_s - thermal.time_s()) * candidate.report.source.radiation_net_power_w();
        let residual = candidate.report.thermal.residual_j + net;
        let external = candidate.report.other_external_heat_j - net;
        if !net.is_finite()
            || !residual.is_finite()
            || !external.is_finite()
            || residual.abs() > controls.total_tolerance_j
        {
            return Err(GraphThermalError::Invalid(
                "closed radiative roundoff cannot supply material graph heat",
            ));
        }
        candidate.report.other_external_heat_j = external;
        candidate.report.closed_radiation_net_j = net;
        cx.checkpoint().map_err(|_| GraphThermalError::Cancelled)?;
        Ok(candidate)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
