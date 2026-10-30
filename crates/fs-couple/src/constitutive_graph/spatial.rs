//! Material feedback on the original heterogeneous nodal-enthalpy field.
//!
//! The native finite-element owner retains its mesh, reference masses, phase
//! charts, conductivity, contact and radiation. Every trial evaluates the
//! original material graph at that same thermal endpoint, from the same
//! committed graph history. An outer Newton solve closes the native nodal
//! balances in joules; no second cell inventory or transport operator exists.
//!
//! Only explicitly watt-valued outputs supply heat. Their nonnegative P1
//! spatial shapes are normalized over the ORIGINAL reference volume. Generic
//! graph dissipation may be a density and is never silently deposited as watts.
//! The caller still owns the work/storage account of the graph's other domains.

use super::{
    MaterialEquation, MaterialGraphCandidate, MaterialGraphError, MaterialGraphResidual,
    MaterialGraphSolve, MaterialUnknown, PortKey, enthalpy::ThermalCoordinate,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyEndpointEvaluation, EnthalpyError, EnthalpyStepConfig,
    heterogeneous::{
        EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError,
    },
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ScalarField, ThermalInterfaces,
};
use fs_exec::Cx;
use fs_material::graph::ConstitutiveGraph;
use fs_qty::Dims;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// One graph input sampled from actual nodal temperatures, after chart inversion.
#[derive(Clone, Debug)]
pub struct GraphSpatialTemperature {
    /// Unique original external input with temperature dimensions.
    pub input: PortKey,
    /// Explicit convex weights `(vertex, weight)` in canonical mesh numbering.
    /// Vertices are distinct and weights sum to one up to roundoff. Sampling
    /// temperature, rather than averaging enthalpies from unlike charts, keeps
    /// each node's material and latent-state authority intact.
    pub vertices: Vec<(usize, f64)>,
    /// Absolute temperature or a difference from an explicit reference.
    pub coordinate: ThermalCoordinate,
}

/// A signed watt-valued output deposited through the original P1 source assembler.
#[derive(Clone, Debug)]
pub struct GraphSpatialHeatPort {
    /// Unique original graph output. An output is deposited exactly once.
    pub output: PortKey,
    /// Nonnegative dimensionless P1 shape with positive volume integral.
    /// The adapter divides by that integral to obtain W/m³ per supplied watt.
    /// This is a volumetric FE shape, not a nodal lumped-load distribution.
    pub distribution: ScalarField,
}

/// Original graph equations, committed history and explicit spatial bindings.
pub struct GraphSpatialInput<'a> {
    pub graph: &'a ConstitutiveGraph,
    pub version: u64,
    pub state: &'a [f64],
    pub fixed: &'a BTreeMap<PortKey, f64>,
    pub unknowns: &'a [MaterialUnknown],
    pub equations: &'a [MaterialEquation],
    pub initial_coordinates: &'a [f64],
    pub temperatures: &'a [GraphSpatialTemperature],
    /// Unique seconds-valued inputs receiving the original physical duration.
    pub durations: &'a [PortKey],
    pub heat: &'a [GraphSpatialHeatPort],
}

/// Bounded coupled work and independent physical acceptance controls.
#[derive(Clone, Copy, Debug)]
pub struct GraphSpatialSolve {
    pub maximum_vertices: usize,
    /// Temperature samples plus all vertex/element heat-shape preparation terms.
    pub maximum_binding_terms: usize,
    pub maximum_endpoint_trials: usize,
    /// Aggregate graph calls across every outer trial and final admission.
    pub maximum_graph_evaluations: usize,
    /// Derivative perturbation in free-node energy increments, J.
    pub difference_step_j: f64,
    /// Every original free-node energy residual must satisfy this bound, J.
    pub node_tolerance_j: f64,
    /// The native controls also govern the outer Newton solve and total energy.
    pub thermal: EnthalpyStepConfig,
    pub graph: MaterialGraphSolve,
}

/// Refusal leaves both committed histories unchanged.
#[derive(Debug)]
pub enum GraphSpatialError {
    Invalid(&'static str),
    Graph(MaterialGraphError),
    Thermal(EnthalpyError),
    Storage(HeterogeneousEnthalpyError),
    Conduction(ConductionError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    EndpointBudget,
    GraphBudget,
    /// Independent physical admission after Newton, including prescribed rows.
    Balance {
        maximum_node_residual_j: f64,
        total_residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for GraphSpatialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material graph/spatial enthalpy coupling: {self:?}")
    }
}
impl std::error::Error for GraphSpatialError {}
impl From<MaterialGraphError> for GraphSpatialError {
    fn from(error: MaterialGraphError) -> Self {
        match error {
            MaterialGraphError::Cancelled => Self::Cancelled,
            other => Self::Graph(other),
        }
    }
}
impl From<ConductionError> for GraphSpatialError {
    fn from(error: ConductionError) -> Self {
        match error {
            ConductionError::Cancelled { .. } => Self::Cancelled,
            other => Self::Conduction(other),
        }
    }
}
impl From<EnthalpyError> for GraphSpatialError {
    fn from(error: EnthalpyError) -> Self {
        match error {
            EnthalpyError::Conduction(ConductionError::Cancelled { .. }) => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl From<HeterogeneousEnthalpyError> for GraphSpatialError {
    fn from(error: HeterogeneousEnthalpyError) -> Self {
        match error {
            HeterogeneousEnthalpyError::Enthalpy(error) => error.into(),
            other => Self::Storage(other),
        }
    }
}

/// The same endpoint's native spatial balances, phases and graph heat account.
pub struct GraphSpatialReport {
    pub thermal: EnthalpyEndpointEvaluation,
    pub liquid_mass_fraction: Vec<f64>,
    pub nonlinear: NewtonReport,
    pub maximum_node_residual_j: f64,
    pub graph_heat_j: f64,
    /// Original volumetric source and physical boundary heat; excludes graph
    /// heat and any net roundoff in a closed radiative enclosure.
    pub other_external_heat_j: f64,
    /// Outward numerical net of a closed enclosure, never treated as supply.
    pub closed_radiation_net_j: f64,
    pub endpoint_trials: usize,
    pub graph_evaluations: usize,
}

/// An inseparable proposal for the nodal heat and material-law histories.
pub struct GraphSpatialCandidate {
    specific_h: Vec<f64>,
    graph: MaterialGraphCandidate,
    report: GraphSpatialReport,
}
impl GraphSpatialCandidate {
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] {
        &self.specific_h
    }
    pub const fn graph(&self) -> &MaterialGraphCandidate {
        &self.graph
    }
    pub const fn report(&self) -> &GraphSpatialReport {
        &self.report
    }
    /// Publish these histories together only after the containing physical
    /// owner has checked its additional domain work and storage accounts.
    pub fn into_parts(self) -> (Vec<f64>, MaterialGraphCandidate, GraphSpatialReport) {
        (self.specific_h, self.graph, self.report)
    }
}

struct Trial {
    h: Vec<f64>,
    phase: Vec<f64>,
    graph: MaterialGraphCandidate,
    physical: EnthalpyEndpointEvaluation,
    graph_power_w: f64,
}

/// Close one spatial material-feedback interval through the native nodal law.
///
/// `initial_energy_increments_j` contains one value per FREE vertex, in ascending
/// mesh order. Prescribed vertices use the native chart's unique bath enthalpy;
/// their reactions remain physical reservoir heat. An entirely prescribed field
/// still admits every source and proposes one physical graph-history update.
///
/// Fixed reference geometry/mass and first-order backward Euler. Native phase
/// charts, temperature-dependent conductivity, contact, convection and explicitly
/// selected radiation remain active. No mechanical-energy certificate, moving
/// mesh, phase kinetics or adjoint for the coupled iteration is implied.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    cx: &Cx<'_>,
    thermal: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    old_specific_h: &[f64],
    dt_s: f64,
    radiation: EnthalpyRadiation<'_>,
    input: GraphSpatialInput<'_>,
    controls: GraphSpatialSolve,
    initial_energy_increments_j: &[f64],
) -> Result<GraphSpatialCandidate, GraphSpatialError> {
    poll(cx)?;
    controls
        .thermal
        .newton
        .validate()
        .map_err(GraphSpatialError::Newton)?;
    let n = problem.mesh.vertex_count();
    if !std::ptr::eq(thermal.mesh(), problem.mesh)
        || n == 0
        || n > controls.maximum_vertices
        || old_specific_h.len() != n
        || !positive(dt_s)
        || controls.maximum_endpoint_trials == 0
        || controls.maximum_graph_evaluations == 0
        || controls.thermal.max_newton_iterations == 0
        || !positive(controls.difference_step_j)
        || !nonnegative(controls.node_tolerance_j)
        || !nonnegative(controls.thermal.energy_tolerance_j)
        || input.unknowns.len() != input.initial_coordinates.len()
        || input.temperatures.len() > controls.graph.maximum_terms
        || input.durations.len() > controls.graph.maximum_terms
        || input.heat.len() > controls.graph.maximum_terms
    {
        return Err(GraphSpatialError::Invalid(
            "original mesh, complete history and finite work/energy controls required",
        ));
    }
    let mut terms = input
        .heat
        .len()
        .checked_mul(n.checked_add(problem.mesh.element_count()).ok_or(
            GraphSpatialError::Invalid("spatial binding extent overflow"),
        )?)
        .ok_or(GraphSpatialError::Invalid(
            "spatial binding extent overflow",
        ))?;
    for binding in input.temperatures {
        terms = terms
            .checked_add(binding.vertices.len())
            .ok_or(GraphSpatialError::Invalid(
                "spatial binding extent overflow",
            ))?;
    }
    if terms > controls.maximum_binding_terms {
        return Err(GraphSpatialError::Invalid(
            "spatial binding work allowance exceeded",
        ));
    }
    // This original owner checks the old charts and uniquely prescribed bath
    // enthalpies before any material graph trial or thermal Newton iteration.
    let baseline =
        thermal.prescribed_endpoint_initial_guess(cx, problem.boundary, old_specific_h)?;
    let prescribed: BTreeSet<_> = problem
        .boundary
        .dirichlet()
        .iter()
        .map(|&(v, _)| v)
        .collect();
    let free: Vec<_> = (0..n).filter(|v| !prescribed.contains(v)).collect();
    if initial_energy_increments_j.len() != free.len()
        || initial_energy_increments_j.iter().any(|v| !v.is_finite())
    {
        return Err(GraphSpatialError::Invalid(
            "one finite energy increment per free vertex required",
        ));
    }
    let mut inputs = BTreeSet::new();
    let mut weight_sums = Vec::with_capacity(input.temperatures.len());
    for binding in input.temperatures {
        poll(cx)?;
        if !inputs.insert(binding.input.clone())
            || input.fixed.contains_key(&binding.input)
            || input.unknowns.iter().any(|u| u.input == binding.input)
            || input
                .graph
                .external_input_dims(&binding.input.0, &binding.input.1)
                != Some(TEMPERATURE)
            || matches!(binding.coordinate, ThermalCoordinate::ChangeFrom(t) if !positive(t))
        {
            return Err(GraphSpatialError::Invalid(
                "temperature must uniquely bind an original external input",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut sum = 0.0;
        for &(vertex, weight) in &binding.vertices {
            poll(cx)?;
            if vertex >= n || !seen.insert(vertex) || !nonnegative(weight) {
                return Err(GraphSpatialError::Invalid(
                    "temperature sampling needs distinct in-range vertices and nonnegative weights",
                ));
            }
            sum += weight;
        }
        if !positive(sum) || (sum - 1.0).abs() > 64.0 * f64::EPSILON * binding.vertices.len() as f64
        {
            return Err(GraphSpatialError::Invalid(
                "temperature sampling weights must sum to one",
            ));
        }
        weight_sums.push(sum);
    }
    for duration in input.durations {
        poll(cx)?;
        if !inputs.insert(duration.clone())
            || input.fixed.contains_key(duration)
            || input.unknowns.iter().any(|u| &u.input == duration)
            || input.graph.external_input_dims(&duration.0, &duration.1) != Some(TIME)
        {
            return Err(GraphSpatialError::Invalid(
                "duration must uniquely bind an original seconds-valued input",
            ));
        }
    }
    let base_power = integrate(cx, problem.mesh, problem.source)?;
    let mut outputs = BTreeSet::new();
    let mut distributions = Vec::with_capacity(input.heat.len());
    for binding in input.heat {
        poll(cx)?;
        if !outputs.insert(binding.output.clone())
            || input
                .graph
                .output_dims(&binding.output.0, &binding.output.1)
                != Some(POWER)
        {
            return Err(GraphSpatialError::Invalid(
                "heat must uniquely bind an original watt-valued output",
            ));
        }
        binding
            .distribution
            .validate("graph heat distribution", n)?;
        for vertex in 0..n {
            poll(cx)?;
            if !nonnegative(binding.distribution.at(vertex)) {
                return Err(GraphSpatialError::Invalid(
                    "heat distribution must be nonnegative",
                ));
            }
        }
        let volume = integrate(cx, problem.mesh, &binding.distribution)?;
        if !positive(volume) {
            return Err(GraphSpatialError::Invalid(
                "heat distribution needs positive reference-volume integral",
            ));
        }
        let density = match &binding.distribution {
            ScalarField::Uniform(value) => ScalarField::Uniform(value / volume),
            ScalarField::Nodal(values) => {
                ScalarField::Nodal(values.iter().map(|v| v / volume).collect())
            }
        };
        density.validate("normalized graph heat distribution", n)?;
        distributions.push(density);
    }
    let trials = Cell::new(0usize);
    let graph_evaluations = Cell::new(0usize);
    let evaluate = |energy: &[f64]| -> Result<Trial, GraphSpatialError> {
        poll(cx)?;
        if trials.get() >= controls.maximum_endpoint_trials {
            return Err(GraphSpatialError::EndpointBudget);
        }
        trials.set(trials.get() + 1);
        if energy.len() != free.len() || energy.iter().any(|e| !e.is_finite()) {
            return Err(GraphSpatialError::Invalid(
                "complete finite endpoint energy coordinates required",
            ));
        }
        let mut h = baseline.clone();
        for (&vertex, &delta) in free.iter().zip(energy) {
            h[vertex] += delta / thermal.reference_nodal_masses_kg()[vertex];
        }
        let mut temperatures = Vec::with_capacity(n);
        let mut phase = Vec::with_capacity(n);
        for (vertex, &value) in h.iter().enumerate() {
            poll(cx)?;
            let state = thermal
                .phase_curve_at(vertex)
                .ok_or(GraphSpatialError::Invalid(
                    "original nodal chart disappeared",
                ))?
                .state_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            temperatures.push(state.temperature_k());
            phase.push(state.liquid_mass_fraction());
        }
        let remaining = controls.maximum_graph_evaluations - graph_evaluations.get();
        if remaining == 0 {
            return Err(GraphSpatialError::GraphBudget);
        }
        let mut fixed = input.fixed.clone();
        for duration in input.durations {
            fixed.insert(duration.clone(), dt_s);
        }
        for (binding, &sum) in input.temperatures.iter().zip(&weight_sums) {
            poll(cx)?;
            let temperature: f64 = binding
                .vertices
                .iter()
                .map(|&(vertex, weight)| weight * temperatures[vertex])
                .sum::<f64>()
                / sum;
            let coordinate = match binding.coordinate {
                ThermalCoordinate::Absolute => temperature,
                ThermalCoordinate::ChangeFrom(reference) => temperature - reference,
            };
            if !coordinate.is_finite() {
                return Err(GraphSpatialError::Invalid("nonfinite sampled temperature"));
            }
            fixed.insert(binding.input.clone(), coordinate);
        }
        let mut graph_controls = controls.graph;
        graph_controls.maximum_evaluations = graph_controls.maximum_evaluations.min(remaining);
        let material = MaterialGraphResidual::new(
            cx,
            input.graph,
            input.version,
            input.state,
            fixed,
            input.unknowns.to_vec(),
            input.equations.to_vec(),
            graph_controls,
        )?;
        let graph = material.solve(input.initial_coordinates.to_vec());
        graph_evaluations.set(graph_evaluations.get() + material.evaluations());
        let graph = graph?;
        let mut source: Vec<_> = (0..n).map(|v| problem.source.at(v)).collect();
        let mut graph_power_w = 0.0;
        for (binding, density) in input.heat.iter().zip(&distributions) {
            poll(cx)?;
            let power =
                *graph
                    .graph
                    .outputs
                    .get(&binding.output)
                    .ok_or(GraphSpatialError::Invalid(
                        "original heat output disappeared",
                    ))?;
            graph_power_w += power;
            for (vertex, value) in source.iter_mut().enumerate() {
                poll(cx)?;
                *value += power * density.at(vertex);
            }
        }
        if !graph_power_w.is_finite() {
            return Err(GraphSpatialError::Invalid("nonfinite graph heat"));
        }
        let source = ScalarField::nodal("combined native volumetric heat", n, source)?;
        let physical = thermal.evaluate_endpoint_with_radiation(
            cx,
            ConductionProblem {
                source: &source,
                ..problem
            },
            interfaces,
            old_specific_h,
            &h,
            dt_s,
            controls.thermal,
            radiation,
        )?;
        poll(cx)?;
        Ok(Trial {
            h,
            phase,
            graph,
            physical,
            graph_power_w,
        })
    };
    let (coordinates, nonlinear) = if free.is_empty() {
        // Source laws and all prescribed reactions are still admitted below.
        (
            Vec::new(),
            NewtonReport {
                iterations: 0,
                residual_norm: 0.0,
                converged: true,
                history: Vec::new(),
                diagnosis: None,
            },
        )
    } else {
        let numerical = |error| match error {
            FiniteDifferenceError::Source(error) => error,
            FiniteDifferenceError::Invalid(reason) => GraphSpatialError::Invalid(reason),
        };
        let residual = FallibleFiniteDifference::new(
            free.len(),
            controls.difference_step_j,
            |energy: &[f64]| {
                evaluate(energy).map(|trial| {
                    free.iter()
                        .map(|&vertex| trial.physical.residual_j[vertex])
                        .collect()
                })
            },
        )
        .map_err(numerical)?;
        let constructed = NewtonKrylovState::new(
            &residual,
            initial_energy_increments_j.to_vec(),
            controls.thermal.newton,
        );
        if let Some(error) = residual.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = constructed.map_err(GraphSpatialError::Newton)?;
        let run = solver.run_controlled(
            &residual,
            controls.thermal.max_newton_iterations,
            &mut || residual.has_failed() || cx.checkpoint().is_err(),
        );
        if let Some(error) = residual.take_failure() {
            return Err(numerical(error));
        }
        poll(cx)?;
        if run.progress != SolverRunProgress::Complete {
            return Err(GraphSpatialError::Interrupted(run.progress));
        }
        if !run.report.converged {
            return Err(GraphSpatialError::NotConverged(run.report));
        }
        (solver.x, run.report)
    };
    let trial = evaluate(&coordinates)?;
    let maximum_node_residual_j = trial
        .physical
        .residual_j
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0, f64::max);
    let graph_heat_j = dt_s * trial.graph_power_w;
    let closed = matches!(
        radiation,
        EnthalpyRadiation::Enclosure { .. } | EnthalpyRadiation::PointwiseEnclosure(_)
    );
    let closed_radiation_net_j = if closed {
        dt_s * trial.physical.radiation_out_w
    } else {
        0.0
    };
    let external_radiation = if closed {
        0.0
    } else {
        trial.physical.radiation_out_w
    };
    let other_external_heat_j = dt_s
        * (base_power + trial.physical.dirichlet_in_w
            - trial.physical.neumann_out_w
            - trial.physical.convective_out_w
            - external_radiation);
    let total_residual_j =
        trial.physical.stored_energy_change_j - graph_heat_j - other_external_heat_j;
    let deposited_residual_j = dt_s * (trial.physical.source_w - base_power - trial.graph_power_w);
    if [
        maximum_node_residual_j,
        graph_heat_j,
        other_external_heat_j,
        closed_radiation_net_j,
        total_residual_j,
        deposited_residual_j,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(GraphSpatialError::Invalid(
            "nonfinite independent spatial heat account",
        ));
    }
    if maximum_node_residual_j > controls.node_tolerance_j
        || trial.physical.energy_residual_j.abs() > controls.thermal.energy_tolerance_j
        || total_residual_j.abs() > controls.thermal.energy_tolerance_j
        || deposited_residual_j.abs() > controls.thermal.energy_tolerance_j
    {
        return Err(GraphSpatialError::Balance {
            maximum_node_residual_j,
            total_residual_j,
        });
    }
    poll(cx)?;
    Ok(GraphSpatialCandidate {
        specific_h: trial.h,
        graph: trial.graph,
        report: GraphSpatialReport {
            thermal: trial.physical,
            liquid_mass_fraction: trial.phase,
            nonlinear,
            maximum_node_residual_j,
            graph_heat_j,
            other_external_heat_j,
            closed_radiation_net_j,
            endpoint_trials: trials.get(),
            graph_evaluations: graph_evaluations.get(),
        },
    })
}

fn integrate(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    field: &ScalarField,
) -> Result<f64, GraphSpatialError> {
    field.validate("spatial heat field", mesh.vertex_count())?;
    let mut integral = 0.0;
    for (element, vertices) in mesh.complex().tets.iter().enumerate() {
        poll(cx)?;
        let mean = vertices
            .iter()
            .map(|&vertex| field.at(vertex as usize))
            .sum::<f64>()
            * 0.25;
        integral += mesh.element_volume(element) * mean;
        if !integral.is_finite() {
            return Err(GraphSpatialError::Invalid(
                "unrepresentable P1 source integral",
            ));
        }
    }
    Ok(integral)
}
fn poll(cx: &Cx<'_>) -> Result<(), GraphSpatialError> {
    cx.checkpoint().map_err(|_| GraphSpatialError::Cancelled)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
