//! Steady incompressible networks of original Carreau-Yasuda pipes.
//!
//! Every pressure trial resolves the pipe's radial shear-dependent viscosity.
//! Junction continuity uses that integrated flow, and Newton's grounded graph
//! Laplacian uses the same profile's differential conductance. No apparent
//! viscosity, secant conductance or regularizing resistance is substituted.
//! The existing checked matrix-free CG owns each linear solve.
//!
//! All branches carry one original rheological fluid at fixed density, with a
//! prescribed uniform temperature within each pipe. Different branches may have
//! different held temperatures; no junction heat mixing or thermal evolution is
//! inferred. Pressure is a gauge relative to one common reference; orientation is only
//! bookkeeping. These are nominal, fully developed, no-slip laminar estimates.
//! Source rate ceilings, constitutive inversion and quadrature controls remain
//! the original pipe's authority. No inertia, gravity, entrance, turbulence,
//! compressibility, thermal feedback or automatic heat deposition is inferred.

use std::cell::RefCell;
use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::pipe::{
    CarreauPipeControls, CarreauPipeProfile, LaminarCarreauPipe,
};
use fs_qty::{Pressure, VolumetricFlowRate};
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, SolveReport, checked_cg};
use fs_sparse::precond::Precond;

use super::{FixedPressure, GraphError, Sum, checkpoint, finite, periodic_checkpoint};

/// Original-law geometry, temperature and prescribed-port tangents/adjoints.
pub mod sensitivity;
/// Original-law self-heating of an explicit recirculating thermal inventory.
pub mod thermal;
pub mod transient;

/// Original physical pipe and its oriented pressure ports.
#[derive(Clone, Debug, PartialEq)]
pub struct CarreauPipeBranch {
    pub name: String,
    pub from: usize,
    pub to: usize,
    pub pipe: LaminarCarreauPipe,
}

/// An ideal prescribed volume-flow port at a node, in m³/s.
/// Positive flow supplies the pipes; negative flow withdraws fluid. Repeated
/// entries at one node are summed in input order with compensated arithmetic.
/// A pressure reservoir may coexist at the node and supplies the remaining
/// exchange. A balanced pair represents an ideal flow-controlled transfer;
/// its pressure work is resolved by the network, with no equipment-efficiency
/// or head-capability model implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowInjection {
    pub node: usize,
    pub flow: VolumetricFlowRate,
}

/// Independent network, linear and constitutive accuracy/work controls.
/// Small continuity residuals do not certify the pipe quadrature's accuracy.
#[derive(Clone, Copy, Debug)]
pub struct CarreauGraphControls {
    pub absolute_flow_tolerance: VolumetricFlowRate,
    pub relative_flow_tolerance: f64,
    /// Total Newton corrections across accepted AND rejected pressure stages.
    pub maximum_newton_iterations: usize,
    /// Total attempts to advance prescribed pressures AND injections from rest
    /// to their full targets; the pressure-only entry point uses the same cap.
    pub maximum_pressure_attempts: usize,
    /// Maximum halvings of each Newton correction; at most 64.
    pub maximum_backtracks: usize,
    pub linear: CheckedCgConfig,
    pub pipe: CarreauPipeControls,
}

impl Default for CarreauGraphControls {
    fn default() -> Self {
        Self {
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
            relative_flow_tolerance: 1e-9,
            maximum_newton_iterations: 128,
            maximum_pressure_attempts: 64,
            maximum_backtracks: 32,
            linear: CheckedCgConfig {
                tolerance: 1e-11,
                max_iterations: 1024,
                max_corrections: 2,
            },
            pipe: CarreauPipeControls::default(),
        }
    }
}

impl CarreauGraphControls {
    fn validate(self) -> Result<(), GraphError> {
        if !self.absolute_flow_tolerance.value().is_finite()
            || self.absolute_flow_tolerance.value() <= 0.0
            || !self.relative_flow_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_flow_tolerance)
            || self.maximum_newton_iterations == 0
            || self.maximum_pressure_attempts == 0
            || self.maximum_backtracks > 64
            || self.linear.max_iterations == 0
            || !self.linear.tolerance.is_finite()
            || self.linear.tolerance <= 0.0
            || self.linear.tolerance >= 1.0
        {
            return Err(GraphError::InvalidInput(
                "invalid Carreau graph tolerances or work bounds",
            ));
        }
        Ok(())
    }

    fn threshold(self, magnitude: f64) -> Result<f64, GraphError> {
        finite(
            self.absolute_flow_tolerance.value() + self.relative_flow_tolerance * (0.5 * magnitude),
            "Carreau continuity tolerance",
        )
    }
}

/// No refused trial is returned as a network solution.
#[derive(Debug)]
pub enum CarreauGraphError {
    Graph(GraphError),
    /// Original branch order and unchanged material refusal.
    Pipe {
        branch: usize,
        source: MaterialError,
    },
    LinearDidNotConverge(SolveReport),
    NewtonIterationLimit {
        iterations: usize,
    },
    ResidualDidNotConverge {
        max_node_imbalance: VolumetricFlowRate,
        boundary_imbalance: VolumetricFlowRate,
    },
    /// Only the full target pressures and injections are publishable. The last
    /// refusal records a source ceiling, numerical failure or unsuccessful
    /// damped correction.
    PressureContinuationLimit {
        reached_fraction: f64,
        attempts: usize,
        last_refusal: Option<Box<CarreauGraphError>>,
    },
}

impl From<GraphError> for CarreauGraphError {
    fn from(value: GraphError) -> Self {
        Self::Graph(value)
    }
}
impl core::fmt::Display for CarreauGraphError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Carreau pipe graph: {self:?}")
    }
}
impl std::error::Error for CarreauGraphError {}

impl CarreauGraphError {
    fn terminal(&self) -> bool {
        matches!(
            self,
            Self::Graph(GraphError::Interrupted) | Self::NewtonIterationLimit { .. }
        )
    }
}

/// One accepted branch, retaining its original source and complete radial state.
#[derive(Clone, Debug, PartialEq)]
pub struct CarreauBranchFlow {
    pub name: String,
    pub from: usize,
    pub to: usize,
    /// The exact profile used by the accepted junction balance, not a replay.
    pub profile: CarreauPipeProfile,
    /// Actual differential dQ/d(delta_p), finite and positive even at rest.
    pub differential_conductance_m3_s_pa: f64,
}

/// Nominal fixed-density pressure/flow state and its signed power accounting.
#[derive(Clone, Debug)]
pub struct CarreauGraphSolution {
    pub pressures: Vec<Pressure>,
    pub branches: Vec<CarreauBranchFlow>,
    /// Positive outflow means the node supplies volume to its incident pipes.
    pub node_outflows: Vec<VolumetricFlowRate>,
    pub boundaries: Vec<FixedPressure>,
    /// Original prescribed flow ports, including signed and repeated entries.
    pub injections: Vec<FlowInjection>,
    /// Net reservoir outflows in pressure-boundary input order: actual pipe
    /// outflow minus the authored net injection at that node.
    pub reservoir_outflows: Vec<VolumetricFlowRate>,
    /// Maximum free-node |pipe outflow - prescribed injection|.
    pub max_node_imbalance: VolumetricFlowRate,
    /// Net reservoir exchange plus all prescribed injections. Fixed-node
    /// injections cancel their reservoir contribution before accumulation.
    pub boundary_imbalance: VolumetricFlowRate,
    /// Sum of prescribed pressure times NET signed reservoir outflow [W].
    pub boundary_power_w: f64,
    /// Sum of solved nodal pressure times prescribed net injection [W]. This
    /// is signed work, not viscous heat; withdrawal can return pressure work.
    pub injection_power_w: f64,
    /// Sum of original accepted pipe pressure times flow [W], nonnegative.
    pub dissipated_power_w: f64,
    /// Sum of (accepted minus requested pipe pressure) times its flow [W].
    pub constitutive_residual_power_w: f64,
    /// Sum of free-node pressure times its remaining continuity residual [W].
    pub node_residual_power_w: f64,
    /// boundary + injection + node residual + constitutive residual
    /// - dissipation [W].
    /// This is an observed floating-point defect, not an energy certificate.
    pub power_balance_residual_w: f64,
    pub newton_iterations: usize,
    /// Attempts to advance the prescribed pressure and flow forcing together.
    pub pressure_attempts: usize,
}

/// Immutable pipe topology. Parallel branches and disconnected components are
/// supported. Steady components require a pressure reservoir; the transient
/// owner can instead anchor a sealed component with its physical storage.
#[derive(Clone, Debug)]
pub struct CarreauPipeGraph {
    branches: Vec<CarreauPipeBranch>,
    adjacency: Vec<Vec<usize>>,
}

impl CarreauPipeGraph {
    /// Retain a single original rheological fluid at fixed density. Geometry,
    /// prescribed temperatures and laminar rate ceilings can differ by pipe.
    /// Each pipe remains axially developed and isothermal within its section.
    pub fn new(
        node_count: usize,
        branches: Vec<CarreauPipeBranch>,
    ) -> Result<Self, CarreauGraphError> {
        if node_count == 0 || branches.is_empty() {
            return Err(GraphError::InvalidInput("nonempty Carreau pipe graph required").into());
        }
        let mut names = BTreeSet::new();
        let mut adjacency = vec![Vec::new(); node_count];
        let reference = branches[0].pipe;
        for (index, branch) in branches.iter().enumerate() {
            let invalid = |reason| GraphError::InvalidBranch {
                branch: index,
                reason,
            };
            if branch.from >= node_count || branch.to >= node_count || branch.from == branch.to {
                return Err(invalid("endpoints must be distinct valid node indices").into());
            }
            if branch.name.trim().is_empty() || !names.insert(branch.name.clone()) {
                return Err(invalid("branch names must be nonblank and unique").into());
            }
            if branch.pipe.law() != reference.law() {
                return Err(invalid("Carreau graph requires one original rheological law").into());
            }
            adjacency[branch.from].push(index);
            adjacency[branch.to].push(index);
        }
        Ok(Self {
            branches,
            adjacency,
        })
    }

    pub fn node_count(&self) -> usize {
        self.adjacency.len()
    }
    pub fn branches(&self) -> &[CarreauPipeBranch] {
        &self.branches
    }

    /// Re-evaluate every original pipe at its prescribed local temperature.
    ///
    /// Temperatures follow original branch order. Every original thermal law,
    /// geometry and physical rate ceiling is retained; a source refusal or
    /// cancellation returns no partially updated graph. The caller owns thermal
    /// evolution and junction mixing. Each section remains isothermal, and the
    /// common fixed-density approximation must remain valid across the change.
    pub fn at_temperatures(
        &self,
        temperatures_k: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Self, CarreauGraphError> {
        checkpoint(cx)?;
        if temperatures_k.len() != self.branches.len() {
            return Err(GraphError::InvalidInput(
                "one prescribed temperature per Carreau branch required",
            )
            .into());
        }
        let mut graph = self.clone();
        for (branch, (owner, &temperature)) in
            graph.branches.iter_mut().zip(temperatures_k).enumerate()
        {
            checkpoint(cx)?;
            owner.pipe = owner
                .pipe
                .at_temperature(temperature)
                .map_err(|source| CarreauGraphError::Pipe { branch, source })?;
        }
        checkpoint(cx)?;
        Ok(graph)
    }

    /// Solve continuity using the original nonlinear radial material response.
    ///
    /// A full pressure step is attempted first. If it crosses a pipe's source
    /// support or a Newton step cannot decrease the residual, bounded pressure
    /// continuation retreats and retries. Trial failures never authorize a
    /// constitutive extension. Only a converged state at ALL original target
    /// pressures is returned, with the exact accepted profiles and signed work.
    ///
    /// Newton corrections use O(nodes + pipes) storage, checked CG and Jacobi
    /// preconditioning. Profiles and differential conductances are reused for
    /// each linear solve. Cancellation is checked before/after every original
    /// bounded pipe query and during graph/linear work. A source quadrature
    /// remains subject to its own finite controls, not a new cancellation API.
    pub fn solve(
        &self,
        boundaries: &[FixedPressure],
        controls: CarreauGraphControls,
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphSolution, CarreauGraphError> {
        self.solve_with_injections(boundaries, &[], controls, cx)
    }

    /// Solve the original pipe network with signed prescribed volume inputs.
    ///
    /// Free-node continuity is pipe outflow = net authored injection. Pressure
    /// remains unknown at those nodes, so the required source pressure and its
    /// signed work emerge from the original shear-dependent material response.
    /// A pressure boundary still anchors each connected component; its reservoir
    /// exchanges the remaining flow after any injection at the same node.
    ///
    /// Pressure differences and injections continue together from rest without
    /// changing a source ceiling. Demand magnitudes enter continuity tolerances;
    /// the original differential conductance and checked linear solve remain
    /// unchanged. Only convergence at ALL full targets publishes. Numerical
    /// storage is O(nodes + pipes), plus the retained original flow-port list.
    pub fn solve_with_injections(
        &self,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        controls: CarreauGraphControls,
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphSolution, CarreauGraphError> {
        checkpoint(cx)?;
        controls.validate()?;
        let layout = self.boundary_layout(boundaries, injections, cx)?;
        let mut pressures = layout.reference_pressure.clone();
        let mut at_rest = !layout.has_injection;
        for boundary in boundaries {
            checkpoint(cx)?;
            at_rest &= boundary.pressure.value() == layout.reference_pressure[boundary.node];
        }
        // Also admits the original constitutive controls before any retreat is
        // possible: malformed controls do not masquerade as a source ceiling.
        let rest = self.evaluate(&pressures, &layout, 0.0, controls, cx)?;
        if at_rest {
            return self.finish(pressures, boundaries, injections, &layout, rest, 0, 0, cx);
        }
        let (mut reached, mut increment) = (0.0, 1.0);
        let (mut iterations, mut attempts) = (0, 0);
        let mut last_refusal = None;
        loop {
            checkpoint(cx)?;
            if attempts == controls.maximum_pressure_attempts {
                return Err(CarreauGraphError::PressureContinuationLimit {
                    reached_fraction: reached,
                    attempts,
                    last_refusal,
                });
            }
            attempts += 1;
            let target = if increment >= 1.0 - reached {
                1.0
            } else {
                reached + increment
            };
            let mut trial = pressures.clone();
            for boundary in boundaries {
                checkpoint(cx)?;
                let reference = layout.reference_pressure[boundary.node];
                trial[boundary.node] = if target == 1.0 {
                    boundary.pressure.value()
                } else {
                    finite(
                        reference + target * (boundary.pressure.value() - reference),
                        "continued boundary pressure",
                    )?
                };
            }
            match self.solve_stage(trial, &layout, target, controls, &mut iterations, cx) {
                Ok((accepted_pressures, state)) => {
                    if target == 1.0 {
                        return self.finish(
                            accepted_pressures,
                            boundaries,
                            injections,
                            &layout,
                            state,
                            iterations,
                            attempts,
                            cx,
                        );
                    }
                    pressures = accepted_pressures;
                    reached = target;
                    increment = (2.0 * increment).min(1.0 - reached);
                    last_refusal = None;
                }
                Err(error) => {
                    if error.terminal() {
                        return Err(error);
                    }
                    increment = 0.5 * (target - reached);
                    last_refusal = Some(Box::new(error));
                    if reached + increment == reached {
                        return Err(CarreauGraphError::PressureContinuationLimit {
                            reached_fraction: reached,
                            attempts,
                            last_refusal,
                        });
                    }
                }
            }
        }
    }

    fn boundary_layout(
        &self,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        cx: &Cx<'_>,
    ) -> Result<BoundaryLayout, GraphError> {
        self.boundary_layout_with_storage(boundaries, injections, None, cx)
    }

    // A stored pressure is an INITIAL GUESS and a positive transient pressure
    // direction, never a fixed-pressure boundary. Its node stays in the same
    // Schur solve and obeys the original volume/storage equation. The steady
    // owner passes None and still requires a physical reservoir in each component.
    fn boundary_layout_with_storage(
        &self,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        initial_storage: Option<&[Option<Pressure>]>,
        cx: &Cx<'_>,
    ) -> Result<BoundaryLayout, GraphError> {
        if initial_storage.is_some_and(|values| {
            values.len() != self.node_count()
                || values
                    .iter()
                    .flatten()
                    .any(|value| !value.value().is_finite())
        }) {
            return Err(GraphError::InvalidInput("invalid stored pressure anchors"));
        }
        let mut fixed = vec![false; self.node_count()];
        let mut visited = fixed.clone();
        let mut reference_pressure = vec![0.0; self.node_count()];
        let mut queue = Vec::new();
        for boundary in boundaries {
            checkpoint(cx)?;
            if boundary.node >= self.node_count() || !boundary.pressure.value().is_finite() {
                return Err(GraphError::InvalidInput(
                    "invalid Carreau pressure boundary",
                ));
            }
            if fixed[boundary.node] {
                return Err(GraphError::InvalidInput(
                    "a node has duplicate pressure boundaries",
                ));
            }
            fixed[boundary.node] = true;
            reference_pressure[boundary.node] = boundary.pressure.value();
        }
        if let Some(storage) = initial_storage {
            for (node, pressure) in storage.iter().enumerate() {
                periodic_checkpoint(node, cx)?;
                if pressure.is_some() && fixed[node] {
                    return Err(GraphError::InvalidInput(
                        "storage node cannot be an ideal pressure reservoir",
                    ));
                }
            }
        }
        let mut injection_sums = vec![Sum::default(); self.node_count()];
        for (index, injection) in injections.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            if injection.node >= self.node_count() || !injection.flow.value().is_finite() {
                return Err(GraphError::InvalidInput(
                    "invalid Carreau nodal volume injection",
                ));
            }
            injection_sums[injection.node].add(injection.flow.value())?;
        }
        // Every component starts at its OWN reference pressure. Translating a
        // different component's gauge must never consume this fluid's rate
        // envelope or continuation budget. Boundary values are gathered before
        // replacing them by the component midpoint; each vertex/edge is visited
        // only a constant number of times.
        for root in 0..self.node_count() {
            checkpoint(cx)?;
            if visited[root] {
                continue;
            }
            queue.clear();
            queue.push(root);
            visited[root] = true;
            let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
            let mut cursor = 0;
            while cursor < queue.len() {
                checkpoint(cx)?;
                let node = queue[cursor];
                cursor += 1;
                if fixed[node] {
                    low = low.min(reference_pressure[node]);
                    high = high.max(reference_pressure[node]);
                } else if let Some(pressure) = initial_storage.and_then(|values| values[node]) {
                    low = low.min(pressure.value());
                    high = high.max(pressure.value());
                }
                for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
                    periodic_checkpoint(ordinal, cx)?;
                    let other = self.other(edge, node);
                    if !visited[other] {
                        visited[other] = true;
                        queue.push(other);
                    }
                }
            }
            if low == f64::INFINITY {
                return Err(GraphError::UnreferencedNode(root));
            }
            let midpoint = if low == high {
                low
            } else {
                0.5 * low + 0.5 * high
            };
            for (ordinal, &node) in queue.iter().enumerate() {
                periodic_checkpoint(ordinal, cx)?;
                reference_pressure[node] = midpoint;
            }
        }
        let mut free_nodes = Vec::new();
        let mut free_index = vec![None; self.node_count()];
        let mut node_injections = Vec::with_capacity(self.node_count());
        let mut has_injection = false;
        for (node, &is_fixed) in fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let injection = injection_sums[node].value()?;
            node_injections.push(injection);
            has_injection |= injection != 0.0;
            if !is_fixed {
                free_index[node] = Some(free_nodes.len());
                free_nodes.push(node);
            }
        }
        Ok(BoundaryLayout {
            fixed,
            free_nodes,
            free_index,
            reference_pressure,
            node_injections,
            has_injection,
        })
    }

    fn other(&self, edge: usize, node: usize) -> usize {
        let branch = &self.branches[edge];
        if branch.from == node {
            branch.to
        } else {
            branch.from
        }
    }

    fn evaluate(
        &self,
        pressures: &[f64],
        layout: &BoundaryLayout,
        injection_scale: f64,
        controls: CarreauGraphControls,
        cx: &Cx<'_>,
    ) -> Result<FlowState, CarreauGraphError> {
        let mut profiles = Vec::with_capacity(self.branches.len());
        let mut conductance = Vec::with_capacity(self.branches.len());
        let mut sums = vec![Sum::default(); self.node_count()];
        let mut magnitudes = sums.clone();
        for (index, branch) in self.branches.iter().enumerate() {
            checkpoint(cx)?;
            let pressure = finite(
                pressures[branch.from] - pressures[branch.to],
                "Carreau branch pressure",
            )?;
            let source = |source| CarreauGraphError::Pipe {
                branch: index,
                source,
            };
            let profile = branch
                .pipe
                .profile_at_pressure(pressure, controls.pipe)
                .map_err(source)?;
            let slope = profile.pressure_conductance_m3_s_pa().map_err(source)?;
            checkpoint(cx)?;
            let flow = profile.flow().volumetric_flow_m3_s;
            sums[branch.from].add(flow)?;
            sums[branch.to].add(-flow)?;
            magnitudes[branch.from].add(flow.abs())?;
            magnitudes[branch.to].add(flow.abs())?;
            profiles.push(profile);
            conductance.push(slope);
        }
        let mut outflows = Vec::with_capacity(self.node_count());
        let mut residuals = Vec::with_capacity(self.node_count());
        let (mut boundary_sum, mut boundary_magnitude) = (Sum::default(), Sum::default());
        let mut max_node_imbalance: f64 = 0.0;
        let mut converged = true;
        for (node, &is_fixed) in layout.fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let outflow = sums[node].value()?;
            let injection = injection_scale * layout.node_injections[node];
            let residual = finite(outflow - injection, "Carreau injection continuity residual")?;
            outflows.push(outflow);
            residuals.push(residual);
            if is_fixed {
                // The prescribed flow here cancels its own reservoir exchange.
                // Cancel it before forming convergence scales, so adding a port
                // at fixed pressure cannot change the physical pressure solve.
                boundary_sum.add(outflow)?;
                boundary_magnitude.add(outflow.abs())?;
            } else {
                boundary_sum.add(injection)?;
                boundary_magnitude.add(injection.abs())?;
                magnitudes[node].add(injection.abs())?;
                max_node_imbalance = max_node_imbalance.max(residual.abs());
                converged &= residual.abs() <= controls.threshold(magnitudes[node].value()?)?;
            }
        }
        let boundary_imbalance = boundary_sum.value()?;
        converged &= boundary_imbalance.abs() <= controls.threshold(boundary_magnitude.value()?)?;
        checkpoint(cx)?;
        Ok(FlowState {
            profiles,
            conductance,
            outflows,
            residuals,
            max_node_imbalance,
            boundary_imbalance,
            converged,
        })
    }

    fn solve_stage(
        &self,
        mut pressures: Vec<f64>,
        layout: &BoundaryLayout,
        injection_scale: f64,
        controls: CarreauGraphControls,
        iterations: &mut usize,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, FlowState), CarreauGraphError> {
        let mut state = self.evaluate(&pressures, layout, injection_scale, controls, cx)?;
        loop {
            checkpoint(cx)?;
            if state.converged {
                return Ok((pressures, state));
            }
            if *iterations == controls.maximum_newton_iterations {
                return Err(CarreauGraphError::NewtonIterationLimit {
                    iterations: *iterations,
                });
            }
            *iterations += 1;
            let correction = self.newton_correction(layout, &state, controls.linear, cx)?;
            let mut step = 1.0;
            let mut trial = pressures.clone();
            let mut accepted = None;
            let mut last_refusal = None;
            for _ in 0..=controls.maximum_backtracks {
                checkpoint(cx)?;
                for (index, &node) in layout.free_nodes.iter().enumerate() {
                    periodic_checkpoint(index, cx)?;
                    trial[node] = pressures[node] + step * correction[index];
                }
                match self.evaluate(&trial, layout, injection_scale, controls, cx) {
                    Ok(candidate) => {
                        if candidate.converged
                            || candidate.max_node_imbalance
                                < (1.0 - 1e-4 * step) * state.max_node_imbalance
                        {
                            accepted = Some(candidate);
                            break;
                        }
                    }
                    Err(error) => {
                        if error.terminal() {
                            return Err(error);
                        }
                        last_refusal = Some(error);
                    }
                }
                step *= 0.5;
            }
            let Some(candidate) = accepted else {
                return Err(last_refusal.unwrap_or_else(|| state.nonconvergence()));
            };
            pressures = trial;
            state = candidate;
        }
    }

    fn newton_correction(
        &self,
        layout: &BoundaryLayout,
        state: &FlowState,
        config: CheckedCgConfig,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, CarreauGraphError> {
        self.pressure_correction(
            layout,
            &state.residuals,
            &state.conductance,
            None,
            config,
            cx,
        )
    }

    /// The same original conductance operator serves steady pipes and the
    /// transient Schur complement. An explicit storage slope adds only a
    /// nonnegative diagonal; it never replaces a branch constitutive tangent.
    #[allow(clippy::too_many_arguments)]
    fn pressure_correction(
        &self,
        layout: &BoundaryLayout,
        residuals: &[f64],
        conductance: &[f64],
        storage_slope: Option<&[f64]>,
        config: CheckedCgConfig,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, CarreauGraphError> {
        let mut rhs = Vec::with_capacity(layout.free_nodes.len());
        let mut inverse_diagonal = Vec::with_capacity(layout.free_nodes.len());
        for (free, &node) in layout.free_nodes.iter().enumerate() {
            checkpoint(cx)?;
            rhs.push(-residuals[node]);
            let mut diagonal = Sum::default();
            if let Some(storage) = storage_slope {
                diagonal.add(storage[free])?;
            }
            for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
                periodic_checkpoint(ordinal, cx)?;
                diagonal.add(conductance[edge])?;
            }
            let inverse = finite(1.0 / diagonal.value()?, "Carreau graph diagonal inverse")?;
            if inverse <= 0.0 {
                return Err(
                    GraphError::NonFinite("Carreau graph diagonal must remain positive").into(),
                );
            }
            inverse_diagonal.push(inverse);
        }
        let operator = ConductanceLaplacian {
            graph: self,
            layout,
            conductance,
            storage_slope,
            inverse_diagonal: &inverse_diagonal,
            cx,
            failure: RefCell::new(None),
        };
        let solved = checked_cg(&operator, &operator, &rhs, config, |_| operator.check()).map_err(
            |error| match error {
                CheckedCgError::InvalidInput(message) => GraphError::InvalidInput(message).into(),
                CheckedCgError::Interrupted(error) => error,
            },
        )?;
        if !solved.report.converged_euclidean() {
            return Err(CarreauGraphError::LinearDidNotConverge(solved.report));
        }
        checkpoint(cx)?;
        Ok(solved.x)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        pressures: Vec<f64>,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        layout: &BoundaryLayout,
        state: FlowState,
        iterations: usize,
        attempts: usize,
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphSolution, CarreauGraphError> {
        let (mut boundary_power, mut node_power) = (Sum::default(), Sum::default());
        let mut injection_power = Sum::default();
        for (node, &pressure) in pressures.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let power = finite(
                pressure * state.residuals[node],
                "Carreau nodal pressure work",
            )?;
            injection_power.add(finite(
                pressure * layout.node_injections[node],
                "Carreau prescribed-flow pressure work",
            )?)?;
            if layout.fixed[node] {
                boundary_power.add(power)?;
            } else {
                node_power.add(power)?;
            }
        }
        let (mut dissipation, mut constitutive_power) = (Sum::default(), Sum::default());
        let mut branches = Vec::with_capacity(self.branches.len());
        for (index, ((branch, profile), slope)) in self
            .branches
            .iter()
            .zip(state.profiles)
            .zip(state.conductance)
            .enumerate()
        {
            periodic_checkpoint(index, cx)?;
            dissipation.add(profile.flow().dissipated_power_w)?;
            constitutive_power
                .add(profile.flow().pressure_residual_pa * profile.flow().volumetric_flow_m3_s)?;
            branches.push(CarreauBranchFlow {
                name: branch.name.clone(),
                from: branch.from,
                to: branch.to,
                profile,
                differential_conductance_m3_s_pa: slope,
            });
        }
        let boundary_power_w = boundary_power.value()?;
        let injection_power_w = injection_power.value()?;
        let node_residual_power_w = node_power.value()?;
        let dissipated_power_w = dissipation.value()?;
        let constitutive_residual_power_w = constitutive_power.value()?;
        let mut balance = Sum::default();
        for work in [
            boundary_power_w,
            injection_power_w,
            node_residual_power_w,
            constitutive_residual_power_w,
            -dissipated_power_w,
        ] {
            balance.add(work)?;
        }
        let mut retained_injections = Vec::with_capacity(injections.len());
        for (index, &injection) in injections.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            retained_injections.push(injection);
        }
        let mut reservoir_outflows = Vec::with_capacity(boundaries.len());
        for (index, boundary) in boundaries.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            reservoir_outflows.push(VolumetricFlowRate::new(state.residuals[boundary.node]));
        }
        checkpoint(cx)?;
        Ok(CarreauGraphSolution {
            pressures: pressures.into_iter().map(Pressure::new).collect(),
            branches,
            node_outflows: state
                .outflows
                .into_iter()
                .map(VolumetricFlowRate::new)
                .collect(),
            boundaries: boundaries.to_vec(),
            injections: retained_injections,
            reservoir_outflows,
            max_node_imbalance: VolumetricFlowRate::new(state.max_node_imbalance),
            boundary_imbalance: VolumetricFlowRate::new(state.boundary_imbalance),
            boundary_power_w,
            injection_power_w,
            dissipated_power_w,
            constitutive_residual_power_w,
            node_residual_power_w,
            power_balance_residual_w: balance.value()?,
            newton_iterations: iterations,
            pressure_attempts: attempts,
        })
    }
}

#[derive(Debug)]
struct BoundaryLayout {
    fixed: Vec<bool>,
    free_nodes: Vec<usize>,
    free_index: Vec<Option<usize>>,
    reference_pressure: Vec<f64>,
    node_injections: Vec<f64>,
    has_injection: bool,
}

struct FlowState {
    profiles: Vec<CarreauPipeProfile>,
    conductance: Vec<f64>,
    outflows: Vec<f64>,
    residuals: Vec<f64>,
    max_node_imbalance: f64,
    boundary_imbalance: f64,
    converged: bool,
}
impl FlowState {
    fn nonconvergence(&self) -> CarreauGraphError {
        CarreauGraphError::ResidualDidNotConverge {
            max_node_imbalance: VolumetricFlowRate::new(self.max_node_imbalance),
            boundary_imbalance: VolumetricFlowRate::new(self.boundary_imbalance),
        }
    }
}

// Original infallible LinearOp/Precond traits use the same refusal-latch pattern
// as graph::sensitivity. Checked CG polls after each bounded operator batch.
struct ConductanceLaplacian<'a, 'cx> {
    graph: &'a CarreauPipeGraph,
    layout: &'a BoundaryLayout,
    conductance: &'a [f64],
    storage_slope: Option<&'a [f64]>,
    inverse_diagonal: &'a [f64],
    cx: &'a Cx<'cx>,
    failure: RefCell<Option<CarreauGraphError>>,
}

impl ConductanceLaplacian<'_, '_> {
    fn check(&self) -> Result<(), CarreauGraphError> {
        if let Some(error) = self.failure.borrow_mut().take() {
            return Err(error);
        }
        checkpoint(self.cx).map_err(Into::into)
    }
    fn publish(&self, result: Result<(), GraphError>, output: &mut [f64]) {
        if let Err(error) = result {
            let mut failure = self.failure.borrow_mut();
            if failure.is_none() {
                *failure = Some(error.into());
            }
            output.fill(f64::NAN);
        }
    }
    fn apply_laplacian(&self, x: &[f64], y: &mut [f64]) -> Result<(), GraphError> {
        checkpoint(self.cx)?;
        for (index, &node) in self.layout.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, self.cx)?;
            let mut sum = Sum::default();
            if let Some(storage) = self.storage_slope {
                sum.add(storage[index] * x[index])?;
            }
            for (ordinal, &edge) in self.graph.adjacency[node].iter().enumerate() {
                periodic_checkpoint(ordinal, self.cx)?;
                let other = self.graph.other(edge, node);
                let neighbor = self.layout.free_index[other].map_or(0.0, |j| x[j]);
                sum.add(self.conductance[edge] * (x[index] - neighbor))?;
            }
            y[index] = sum.value()?;
        }
        Ok(())
    }
    fn apply_diagonal(&self, r: &[f64], z: &mut [f64]) -> Result<(), GraphError> {
        checkpoint(self.cx)?;
        for (index, &inverse) in self.inverse_diagonal.iter().enumerate() {
            periodic_checkpoint(index, self.cx)?;
            z[index] = finite(r[index] * inverse, "Carreau Jacobi preconditioner")?;
        }
        Ok(())
    }
}
impl LinearOp for ConductanceLaplacian<'_, '_> {
    fn n(&self) -> usize {
        self.layout.free_nodes.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.publish(self.apply_laplacian(x, y), y);
    }
    // A grounded conductance-weighted graph Laplacian: symmetric.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        LinearOp::apply(self, x, y);
    }
}
impl Precond for ConductanceLaplacian<'_, '_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.publish(self.apply_diagonal(r, z), z);
    }
}

#[cfg(test)]
mod tests;
