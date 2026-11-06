//! Coupled original-law radial Carreau pipes and linear hydraulic storage.
//!
//! Every pressure trial advances each ORIGINAL radial finite-volume momentum
//! equation from the same accepted history. At a free node the midpoint
//! equation is `Q_out - injection + C*(p_next-p_old)/dt = 0`, with
//! `p_next = 2*p_mid-p_old`. Its positive Schur Jacobian is the original
//! transient branch conductance Laplacian plus `2*C/dt`. The existing checked
//! CG owns that solve; no fitted mean-flow inertia or apparent viscosity is
//! introduced. Zero-storage junctions retain endpoint flow continuity when
//! their initial histories already satisfy the held volume inputs.
//!
//! Storage is an explicitly authored constant incremental volume compliance
//! on a declared pressure interval. Its energy is `C*p*p/2` in the same pressure
//! reference as the ports. The owner does not infer bulk modulus, nominal
//! volume, density changes, acoustic wave propagation, cavitation, temperature
//! evolution or entrance/junction losses. Radial and temporal discretization
//! errors remain; midpoint energy accounting does not certify time resolution.

use std::sync::Arc;

use fs_material::fluid_transport::rheology::pipe::transient::{
    RadialCarreauControls, RadialCarreauError, RadialCarreauPipe, RadialCarreauState,
    RadialCarreauStep, RadialCarreauWorkspace,
};

use super::{
    BoundaryLayout, CarreauGraphError, CarreauPipeGraph, CheckedCgConfig, Cx, FixedPressure,
    FlowInjection, GraphError, Pressure, Sum, VolumetricFlowRate, checkpoint, finite,
    periodic_checkpoint,
};

/// A bounded constant incremental hydraulic storage law. A caller supplies
/// `dV/dp`, not a dynamic viscosity or an inferred material compressibility.
/// Pressure bounds and all pressures use one common physical reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearHydraulicStorage {
    pub compliance_m3_pa: f64,
    pub minimum_pressure: Pressure,
    pub maximum_pressure: Pressure,
}

impl LinearHydraulicStorage {
    fn validate(self) -> Result<(), GraphError> {
        if !self.compliance_m3_pa.is_finite()
            || self.compliance_m3_pa <= 0.0
            || !self.minimum_pressure.value().is_finite()
            || !self.maximum_pressure.value().is_finite()
            || self.minimum_pressure.value() >= self.maximum_pressure.value()
        {
            return Err(GraphError::InvalidInput("invalid linear hydraulic storage"));
        }
        Ok(())
    }

    fn admits(self, pressure: f64) -> bool {
        pressure.is_finite()
            && pressure >= self.minimum_pressure.value()
            && pressure <= self.maximum_pressure.value()
    }
}

/// Independent coupled and original radial work/accuracy bounds.
#[derive(Clone, Copy, Debug)]
pub struct TransientCarreauControls {
    pub absolute_flow_tolerance: VolumetricFlowRate,
    pub relative_flow_tolerance: f64,
    pub maximum_newton_iterations: usize,
    pub maximum_backtracks: usize,
    /// All attempted original radial solves, including rejected pressure trials.
    pub maximum_branch_solves: u64,
    pub linear: CheckedCgConfig,
    pub radial: RadialCarreauControls,
}

impl Default for TransientCarreauControls {
    fn default() -> Self {
        Self {
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
            relative_flow_tolerance: 1e-9,
            maximum_newton_iterations: 64,
            maximum_backtracks: 24,
            maximum_branch_solves: 65_536,
            linear: CheckedCgConfig {
                tolerance: 1e-11,
                max_iterations: 1024,
                max_corrections: 2,
            },
            radial: RadialCarreauControls::default(),
        }
    }
}

impl TransientCarreauControls {
    fn validate(self) -> Result<(), GraphError> {
        if !self.absolute_flow_tolerance.value().is_finite()
            || self.absolute_flow_tolerance.value() <= 0.0
            || !self.relative_flow_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_flow_tolerance)
            || self.maximum_newton_iterations == 0
            || self.maximum_backtracks > 64
            || self.maximum_branch_solves == 0
            || !self.linear.tolerance.is_finite()
            || !(0.0..1.0).contains(&self.linear.tolerance)
            || self.linear.tolerance == 0.0
            || self.linear.max_iterations == 0
        {
            return Err(GraphError::InvalidInput("invalid transient graph controls"));
        }
        Ok(())
    }

    fn tolerance(self, magnitude: f64) -> Result<f64, GraphError> {
        finite(
            self.absolute_flow_tolerance.value() + self.relative_flow_tolerance * (0.5 * magnitude),
            "transient continuity tolerance",
        )
    }
}

#[derive(Debug)]
pub enum TransientCarreauError {
    Graph(CarreauGraphError),
    Radial {
        branch: usize,
        source: RadialCarreauError,
    },
    StorageDomain {
        node: usize,
        pressure_pa: f64,
    },
    InconsistentInitialFlow {
        node: usize,
        imbalance: VolumetricFlowRate,
    },
    StateMismatch,
    BranchSolveBudget,
    NewtonBudget,
    BacktrackBudget,
    Capacity,
    WorkBalance,
}

impl From<GraphError> for TransientCarreauError {
    fn from(value: GraphError) -> Self {
        Self::Graph(value.into())
    }
}
impl From<CarreauGraphError> for TransientCarreauError {
    fn from(value: CarreauGraphError) -> Self {
        Self::Graph(value)
    }
}
impl core::fmt::Display for TransientCarreauError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient Carreau graph: {self:?}")
    }
}
impl std::error::Error for TransientCarreauError {}

impl TransientCarreauError {
    fn terminal(&self) -> bool {
        matches!(
            self,
            Self::Graph(CarreauGraphError::Graph(GraphError::Interrupted))
                | Self::Radial {
                    source: RadialCarreauError::Interrupted,
                    ..
                }
                | Self::BranchSolveBudget
        )
    }
}

/// Immutable original material, topology, radial meshes and storage laws.
#[derive(Clone, Debug)]
pub struct TransientCarreauGraph {
    graph: CarreauPipeGraph,
    pipes: Vec<RadialCarreauPipe>,
    storage: Vec<Option<LinearHydraulicStorage>>,
    identity: Arc<()>,
}

/// Accepted complete history. Storage pressures exist only at storage nodes;
/// zero-storage pressures are algebraic midpoint outputs, not a hidden state.
#[derive(Clone, Debug)]
pub struct TransientCarreauState {
    branches: Vec<RadialCarreauState>,
    storage_pressures: Vec<Option<Pressure>>,
    accepted_time_s: f64,
    identity: Arc<()>,
}

impl TransientCarreauState {
    pub fn branches(&self) -> &[RadialCarreauState] {
        &self.branches
    }
    pub fn storage_pressures(&self) -> &[Option<Pressure>] {
        &self.storage_pressures
    }
    pub fn accepted_time_s(&self) -> f64 {
        self.accepted_time_s
    }
}

/// Persistent radial and graph trial buffers. The original checked CG and
/// owned returned report still allocate their graph-sized result storage.
#[derive(Debug)]
pub struct TransientCarreauWorkspace {
    radial: Vec<RadialCarreauWorkspace>,
    outputs: Vec<RadialCarreauState>,
    pressures: Vec<f64>,
    trial_pressures: Vec<f64>,
    current: Evaluation,
    trial: Evaluation,
    identity: Arc<()>,
}

#[derive(Debug)]
struct Evaluation {
    steps: Vec<Option<RadialCarreauStep>>,
    conductance: Vec<f64>,
    residuals: Vec<f64>,
    outflows: Vec<Sum>,
    final_outflows: Vec<Sum>,
    magnitudes: Vec<Sum>,
    final_magnitudes: Vec<Sum>,
    storage_pressures: Vec<Option<Pressure>>,
    storage_changes: Vec<f64>,
    max_imbalance: f64,
    max_final_rigid_imbalance: f64,
    boundary_imbalance: f64,
    converged: bool,
}

impl Evaluation {
    fn new(nodes: usize, branches: usize) -> Self {
        Self {
            steps: vec![None; branches],
            conductance: vec![0.0; branches],
            residuals: vec![0.0; nodes],
            outflows: vec![Sum::default(); nodes],
            final_outflows: vec![Sum::default(); nodes],
            magnitudes: vec![Sum::default(); nodes],
            final_magnitudes: vec![Sum::default(); nodes],
            storage_pressures: vec![None; nodes],
            storage_changes: vec![0.0; nodes],
            max_imbalance: 0.0,
            max_final_rigid_imbalance: 0.0,
            boundary_imbalance: 0.0,
            converged: false,
        }
    }

    fn clear(&mut self) {
        self.steps.fill(None);
        self.conductance.fill(0.0);
        self.residuals.fill(0.0);
        self.outflows.fill(Sum::default());
        self.final_outflows.fill(Sum::default());
        self.magnitudes.fill(Sum::default());
        self.final_magnitudes.fill(Sum::default());
        self.storage_pressures.fill(None);
        self.storage_changes.fill(0.0);
        self.max_imbalance = 0.0;
        self.max_final_rigid_imbalance = 0.0;
        self.boundary_imbalance = 0.0;
        self.converged = true;
    }
}

/// Accepted midpoint fluxes and physical work over one common interval.
#[derive(Clone, Debug)]
pub struct TransientCarreauStep {
    pub duration_s: f64,
    pub midpoint_pressures: Vec<Pressure>,
    pub branches: Vec<RadialCarreauStep>,
    pub node_outflows: Vec<VolumetricFlowRate>,
    pub storage_volume_changes_m3: Vec<f64>,
    /// Input boundary order; positive volume leaves its reservoir.
    pub reservoir_outflows: Vec<VolumetricFlowRate>,
    pub max_node_imbalance: VolumetricFlowRate,
    pub max_final_rigid_node_imbalance: VolumetricFlowRate,
    pub boundary_imbalance: VolumetricFlowRate,
    pub reservoir_work_j: f64,
    pub injection_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub storage_energy_change_j: f64,
    pub viscous_heat_j: f64,
    /// Sum of p_mid * remaining free-node continuity residual * dt.
    pub node_residual_work_j: f64,
    /// Sum of original radial pressure work - kinetic change - heat.
    pub radial_work_residual_j: f64,
    /// p_mid*delta_V - change(C*p*p/2), due only to midpoint rounding.
    pub storage_midpoint_work_residual_j: f64,
    /// External work - kinetic change - storage change - source heat.
    pub physical_work_residual_j: f64,
    /// Physical residual + node residual - radial residual - storage residual.
    pub accounting_residual_j: f64,
    /// Arithmetic consistency gate, not a physical or discretization bound.
    pub accounting_tolerance_j: f64,
    pub newton_iterations: usize,
    pub branch_solves: u64,
}

impl TransientCarreauGraph {
    /// Bind one fixed-density fluid to the original Carreau pipe topology.
    /// An absent storage law denotes a zero-volume ideal junction. Every
    /// component still requires an explicitly supplied pressure reservoir.
    pub fn new(
        graph: CarreauPipeGraph,
        density_kg_m3: f64,
        radial_cells: &[usize],
        storage: Vec<Option<LinearHydraulicStorage>>,
    ) -> Result<Self, TransientCarreauError> {
        if radial_cells.len() != graph.branches.len() || storage.len() != graph.node_count() {
            return Err(GraphError::InvalidInput("transient mesh/storage dimensions").into());
        }
        for law in storage.iter().flatten() {
            law.validate()?;
        }
        let mut pipes = Vec::new();
        pipes
            .try_reserve_exact(graph.branches.len())
            .map_err(|_| TransientCarreauError::Capacity)?;
        for (branch, (&cells, owner)) in radial_cells.iter().zip(&graph.branches).enumerate() {
            pipes.push(
                RadialCarreauPipe::new(owner.pipe, density_kg_m3, cells)
                    .map_err(|source| TransientCarreauError::Radial { branch, source })?,
            );
        }
        Ok(Self {
            graph,
            pipes,
            storage,
            identity: Arc::new(()),
        })
    }

    pub fn graph(&self) -> &CarreauPipeGraph {
        &self.graph
    }
    pub fn pipes(&self) -> &[RadialCarreauPipe] {
        &self.pipes
    }
    pub fn storage(&self) -> &[Option<LinearHydraulicStorage>] {
        &self.storage
    }

    fn admit_pressures(&self, pressures: &[Option<Pressure>]) -> Result<(), TransientCarreauError> {
        if pressures.len() != self.graph.node_count() {
            return Err(GraphError::InvalidInput("storage pressure dimensions").into());
        }
        for (node, (&law, &pressure)) in self.storage.iter().zip(pressures).enumerate() {
            match (law, pressure) {
                (None, None) => {}
                (Some(law), Some(pressure)) if law.admits(pressure.value()) => {}
                (Some(_), Some(pressure)) => {
                    return Err(TransientCarreauError::StorageDomain {
                        node,
                        pressure_pa: pressure.value(),
                    });
                }
                _ => return Err(GraphError::InvalidInput("storage pressure presence").into()),
            }
        }
        Ok(())
    }

    pub fn rest_state(
        &self,
        storage_pressures: &[Option<Pressure>],
        cx: &Cx<'_>,
    ) -> Result<TransientCarreauState, TransientCarreauError> {
        checkpoint(cx)?;
        self.admit_pressures(storage_pressures)?;
        let mut branches = Vec::with_capacity(self.pipes.len());
        for (branch, pipe) in self.pipes.iter().enumerate() {
            checkpoint(cx)?;
            branches.push(
                pipe.rest_state()
                    .map_err(|source| TransientCarreauError::Radial { branch, source })?,
            );
        }
        checkpoint(cx)?;
        Ok(TransientCarreauState {
            branches,
            storage_pressures: storage_pressures.to_vec(),
            accepted_time_s: 0.0,
            identity: self.identity.clone(),
        })
    }

    /// Import the complete original radial histories at one common time.
    /// Junction compatibility is checked against the actual boundaries and
    /// held injections at step admission, when those constraints are known.
    pub fn state_from_velocities(
        &self,
        velocities_m_s: Vec<Vec<f64>>,
        storage_pressures: &[Option<Pressure>],
        accepted_time_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientCarreauState, TransientCarreauError> {
        checkpoint(cx)?;
        self.admit_pressures(storage_pressures)?;
        if velocities_m_s.len() != self.pipes.len() {
            return Err(GraphError::InvalidInput("radial history branch dimensions").into());
        }
        let mut branches = Vec::with_capacity(self.pipes.len());
        for (branch, (pipe, velocities)) in self.pipes.iter().zip(velocities_m_s).enumerate() {
            branches.push(
                pipe.state_from_velocities(velocities, accepted_time_s, || {
                    cx.checkpoint().is_err()
                })
                .map_err(|source| TransientCarreauError::Radial { branch, source })?,
            );
        }
        checkpoint(cx)?;
        Ok(TransientCarreauState {
            branches,
            storage_pressures: storage_pressures.to_vec(),
            accepted_time_s,
            identity: self.identity.clone(),
        })
    }

    pub fn workspace(&self) -> Result<TransientCarreauWorkspace, TransientCarreauError> {
        let mut radial = Vec::with_capacity(self.pipes.len());
        let mut outputs = Vec::with_capacity(self.pipes.len());
        for (branch, pipe) in self.pipes.iter().enumerate() {
            let source = |source| TransientCarreauError::Radial { branch, source };
            radial.push(pipe.workspace().map_err(source)?);
            outputs.push(pipe.rest_state().map_err(source)?);
        }
        let nodes = self.graph.node_count();
        Ok(TransientCarreauWorkspace {
            radial,
            outputs,
            pressures: vec![0.0; nodes],
            trial_pressures: vec![0.0; nodes],
            current: Evaluation::new(nodes, self.pipes.len()),
            trial: Evaluation::new(nodes, self.pipes.len()),
            identity: self.identity.clone(),
        })
    }

    fn initial_continuity(
        &self,
        state: &TransientCarreauState,
        layout: &BoundaryLayout,
        scratch: &mut Evaluation,
        controls: TransientCarreauControls,
        cx: &Cx<'_>,
    ) -> Result<(), TransientCarreauError> {
        scratch.clear();
        for (edge, pipe) in self.pipes.iter().enumerate() {
            checkpoint(cx)?;
            let mut volume = Sum::default();
            let dr = pipe.spacing_m();
            let area_unit = (core::f64::consts::PI * dr) * dr;
            for (cell, &velocity) in state.branches[edge].velocities_m_s().iter().enumerate() {
                periodic_checkpoint(cell, cx)?;
                volume.add(area_unit * (2 * cell + 1) as f64 * velocity)?;
            }
            let q = volume.value()?;
            let branch = &self.graph.branches[edge];
            scratch.outflows[branch.from].add(q)?;
            scratch.outflows[branch.to].add(-q)?;
            scratch.magnitudes[branch.from].add(q.abs())?;
            scratch.magnitudes[branch.to].add(q.abs())?;
        }
        for &node in &layout.free_nodes {
            checkpoint(cx)?;
            if self.storage[node].is_none() {
                let imbalance = scratch.outflows[node].value()? - layout.node_injections[node];
                let magnitude =
                    scratch.magnitudes[node].value()? + layout.node_injections[node].abs();
                if imbalance.abs() > controls.tolerance(magnitude)? {
                    return Err(TransientCarreauError::InconsistentInitialFlow {
                        node,
                        imbalance: VolumetricFlowRate::new(imbalance),
                    });
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        initial: &TransientCarreauState,
        pressures: &[f64],
        layout: &BoundaryLayout,
        duration_s: f64,
        controls: TransientCarreauControls,
        radial: &mut [RadialCarreauWorkspace],
        outputs: &mut [RadialCarreauState],
        result: &mut Evaluation,
        solves: &mut u64,
        cx: &Cx<'_>,
    ) -> Result<(), TransientCarreauError> {
        result.clear();
        // Admit all actual storage endpoints before any expensive pipe trial.
        for (node, &law) in self.storage.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            if let Some(law) = law {
                let old = initial.storage_pressures[node]
                    .expect("owned storage state")
                    .value();
                let next = finite(
                    pressures[node] + (pressures[node] - old),
                    "storage endpoint",
                )?;
                if !law.admits(next) {
                    return Err(TransientCarreauError::StorageDomain {
                        node,
                        pressure_pa: next,
                    });
                }
                result.storage_pressures[node] = Some(Pressure::new(next));
                result.storage_changes[node] =
                    finite(law.compliance_m3_pa * (next - old), "stored volume change")?;
            }
        }
        for (edge, pipe) in self.pipes.iter().enumerate() {
            checkpoint(cx)?;
            if *solves == controls.maximum_branch_solves {
                return Err(TransientCarreauError::BranchSolveBudget);
            }
            *solves += 1;
            let branch = &self.graph.branches[edge];
            let pressure = finite(
                pressures[branch.from] - pressures[branch.to],
                "transient pipe pressure",
            )?;
            let step = pipe
                .step_from(
                    &initial.branches[edge],
                    &mut outputs[edge],
                    &mut radial[edge],
                    pressure,
                    duration_s,
                    controls.radial,
                    || cx.checkpoint().is_err(),
                )
                .map_err(|source| TransientCarreauError::Radial {
                    branch: edge,
                    source,
                })?;
            result.conductance[edge] = step.midpoint_pressure_conductance_m3_s_pa;
            let flow = step.midpoint_volume_flow_m3_s;
            let final_flow = step.final_volume_flow_m3_s;
            for (node, sign) in [(branch.from, 1.0), (branch.to, -1.0)] {
                result.outflows[node].add(sign * flow)?;
                result.final_outflows[node].add(sign * final_flow)?;
                result.magnitudes[node].add(flow.abs())?;
                result.final_magnitudes[node].add(final_flow.abs())?;
            }
            result.steps[edge] = Some(step);
        }
        let (mut net, mut magnitude) = (Sum::default(), Sum::default());
        for (node, &fixed) in layout.fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let flow = result.outflows[node].value()?;
            let injection = layout.node_injections[node];
            let storage_rate = finite(
                result.storage_changes[node] / duration_s,
                "storage volume rate",
            )?;
            let mut residual = Sum::default();
            residual.add(flow)?;
            residual.add(-injection)?;
            residual.add(storage_rate)?;
            result.residuals[node] = residual.value()?;
            if fixed {
                net.add(flow)?;
                magnitude.add(flow.abs())?;
            } else {
                net.add(injection)?;
                net.add(-storage_rate)?;
                magnitude.add(injection.abs())?;
                magnitude.add(storage_rate.abs())?;
                result.magnitudes[node].add(injection.abs())?;
                result.magnitudes[node].add(storage_rate.abs())?;
                result.max_imbalance = result.max_imbalance.max(result.residuals[node].abs());
                result.converged &= result.residuals[node].abs()
                    <= controls.tolerance(result.magnitudes[node].value()?)?;
                if self.storage[node].is_none() {
                    let endpoint = finite(
                        result.final_outflows[node].value()? - injection,
                        "endpoint junction continuity",
                    )?;
                    result.final_magnitudes[node].add(injection.abs())?;
                    result.max_final_rigid_imbalance =
                        result.max_final_rigid_imbalance.max(endpoint.abs());
                    result.converged &= endpoint.abs()
                        <= controls.tolerance(result.final_magnitudes[node].value()?)?;
                }
            }
        }
        result.boundary_imbalance = net.value()?;
        result.converged &=
            result.boundary_imbalance.abs() <= controls.tolerance(magnitude.value()?)?;
        checkpoint(cx)?;
        Ok(())
    }

    /// Advance all radial histories and storage pressures together. Reservoir
    /// pressures and volume injections are held over this interval. Storage
    /// nodes cannot also be ideal pressure reservoirs. A held injection at a
    /// zero-storage node requires a compatible initial radial flow; changing
    /// it discontinuously refuses instead of inventing an impulsive pressure.
    /// Failure leaves the complete accepted state and time unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        state: &mut TransientCarreauState,
        workspace: &mut TransientCarreauWorkspace,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        duration_s: f64,
        controls: TransientCarreauControls,
        cx: &Cx<'_>,
    ) -> Result<TransientCarreauStep, TransientCarreauError> {
        checkpoint(cx)?;
        controls.validate()?;
        if !Arc::ptr_eq(&self.identity, &state.identity)
            || !Arc::ptr_eq(&self.identity, &workspace.identity)
        {
            return Err(TransientCarreauError::StateMismatch);
        }
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(GraphError::InvalidInput("positive transient duration required").into());
        }
        let next_time = finite(state.accepted_time_s + duration_s, "transient graph time")?;
        if next_time <= state.accepted_time_s {
            return Err(GraphError::NonFinite("transient time increment is unresolved").into());
        }
        let layout = self.graph.boundary_layout(boundaries, injections, cx)?;
        let mut storage_slope = Vec::with_capacity(layout.free_nodes.len());
        for (node, &fixed) in layout.fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            if fixed && self.storage[node].is_some() {
                return Err(GraphError::InvalidInput(
                    "storage node cannot be an ideal pressure reservoir",
                )
                .into());
            }
            workspace.pressures[node] = state.storage_pressures[node]
                .map_or(layout.reference_pressure[node], Pressure::value);
        }
        for &node in &layout.free_nodes {
            checkpoint(cx)?;
            storage_slope.push(match self.storage[node] {
                None => 0.0,
                Some(storage) => {
                    let slope = finite(
                        2.0 * (storage.compliance_m3_pa / duration_s),
                        "storage pressure slope",
                    )?;
                    if slope <= 0.0 {
                        return Err(GraphError::NonFinite("storage slope underflows").into());
                    }
                    slope
                }
            });
        }
        for boundary in boundaries {
            workspace.pressures[boundary.node] = boundary.pressure.value();
        }
        self.initial_continuity(state, &layout, &mut workspace.current, controls, cx)?;
        let mut branch_solves = 0;
        self.evaluate(
            state,
            &workspace.pressures,
            &layout,
            duration_s,
            controls,
            &mut workspace.radial,
            &mut workspace.outputs,
            &mut workspace.current,
            &mut branch_solves,
            cx,
        )?;
        let mut iterations = 0;
        while !workspace.current.converged {
            checkpoint(cx)?;
            if iterations == controls.maximum_newton_iterations {
                return Err(TransientCarreauError::NewtonBudget);
            }
            iterations += 1;
            let correction = self.graph.pressure_correction(
                &layout,
                &workspace.current.residuals,
                &workspace.current.conductance,
                Some(&storage_slope),
                controls.linear,
                cx,
            )?;
            let mut fraction = 1.0_f64;
            let mut accepted = false;
            let mut last_refusal = None;
            for _ in 0..=controls.maximum_backtracks {
                checkpoint(cx)?;
                workspace
                    .trial_pressures
                    .copy_from_slice(&workspace.pressures);
                for (index, &node) in layout.free_nodes.iter().enumerate() {
                    periodic_checkpoint(index, cx)?;
                    workspace.trial_pressures[node] =
                        fraction.mul_add(correction[index], workspace.pressures[node]);
                }
                match self.evaluate(
                    state,
                    &workspace.trial_pressures,
                    &layout,
                    duration_s,
                    controls,
                    &mut workspace.radial,
                    &mut workspace.outputs,
                    &mut workspace.trial,
                    &mut branch_solves,
                    cx,
                ) {
                    Ok(())
                        if workspace.trial.converged
                            || workspace.trial.max_imbalance
                                < (1.0 - 1e-4 * fraction) * workspace.current.max_imbalance =>
                    {
                        accepted = true;
                        break;
                    }
                    Ok(()) => {}
                    Err(error) if error.terminal() => return Err(error),
                    Err(error) => last_refusal = Some(error),
                }
                fraction *= 0.5;
            }
            if !accepted {
                return Err(last_refusal.unwrap_or(TransientCarreauError::BacktrackBudget));
            }
            core::mem::swap(&mut workspace.current, &mut workspace.trial);
            core::mem::swap(&mut workspace.pressures, &mut workspace.trial_pressures);
        }
        let report = self.finish(
            state,
            &workspace.current,
            &workspace.pressures,
            &layout,
            boundaries,
            duration_s,
            iterations,
            branch_solves,
            cx,
        )?;
        checkpoint(cx)?;
        core::mem::swap(&mut state.branches, &mut workspace.outputs);
        core::mem::swap(
            &mut state.storage_pressures,
            &mut workspace.current.storage_pressures,
        );
        state.accepted_time_s = next_time;
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        initial: &TransientCarreauState,
        result: &Evaluation,
        pressures: &[f64],
        layout: &BoundaryLayout,
        boundaries: &[FixedPressure],
        duration_s: f64,
        newton_iterations: usize,
        branch_solves: u64,
        cx: &Cx<'_>,
    ) -> Result<TransientCarreauStep, TransientCarreauError> {
        let (mut reservoir, mut injection, mut node_work) =
            (Sum::default(), Sum::default(), Sum::default());
        let (mut stored, mut storage_defect) = (Sum::default(), Sum::default());
        let mut scale = Sum::default();
        let mut node_outflows = Vec::with_capacity(pressures.len());
        for (node, &pressure) in pressures.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let input_work = finite(
                (pressure * layout.node_injections[node]) * duration_s,
                "injection work",
            )?;
            injection.add(input_work)?;
            scale.add(input_work.abs())?;
            let residual_work = finite(
                (pressure * result.residuals[node]) * duration_s,
                "node work",
            )?;
            scale.add(residual_work.abs())?;
            if layout.fixed[node] {
                reservoir.add(residual_work)?;
            } else {
                node_work.add(residual_work)?;
            }
            if let Some(storage) = self.storage[node] {
                let old = initial.storage_pressures[node]
                    .expect("owned storage state")
                    .value();
                let new = result.storage_pressures[node]
                    .expect("accepted storage state")
                    .value();
                let change = finite(
                    (storage.compliance_m3_pa * (0.5 * old + 0.5 * new)) * (new - old),
                    "storage energy change",
                )?;
                let work = finite(
                    pressure * result.storage_changes[node],
                    "storage pressure work",
                )?;
                stored.add(change)?;
                storage_defect.add(work - change)?;
                scale.add(work.abs())?;
                scale.add(change.abs())?;
            }
            node_outflows.push(VolumetricFlowRate::new(result.outflows[node].value()?));
        }
        let (mut kinetic, mut heat, mut radial_defect) =
            (Sum::default(), Sum::default(), Sum::default());
        let mut branches = Vec::with_capacity(self.pipes.len());
        for (edge, report) in result.steps.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let report = report.expect("fully evaluated accepted branch");
            kinetic.add(report.kinetic_energy_change_j)?;
            heat.add(report.viscous_heat_j)?;
            radial_defect.add(report.work_residual_j)?;
            for value in [
                report.pressure_work_j,
                report.kinetic_energy_change_j,
                report.viscous_heat_j,
                report.work_residual_j,
            ] {
                scale.add(value.abs())?;
            }
            branches.push(report);
        }
        let reservoir_work_j = reservoir.value()?;
        let injection_work_j = injection.value()?;
        let kinetic_energy_change_j = kinetic.value()?;
        let storage_energy_change_j = stored.value()?;
        let viscous_heat_j = heat.value()?;
        let node_residual_work_j = node_work.value()?;
        let radial_work_residual_j = radial_defect.value()?;
        let storage_midpoint_work_residual_j = storage_defect.value()?;
        let mut physical = Sum::default();
        for value in [
            reservoir_work_j,
            injection_work_j,
            -kinetic_energy_change_j,
            -storage_energy_change_j,
            -viscous_heat_j,
        ] {
            physical.add(value)?;
        }
        let physical_work_residual_j = physical.value()?;
        physical.add(node_residual_work_j)?;
        physical.add(-radial_work_residual_j)?;
        physical.add(-storage_midpoint_work_residual_j)?;
        let accounting_residual_j = physical.value()?;
        let accounting_tolerance_j = finite(
            256.0 * f64::EPSILON * scale.value()?,
            "work accounting tolerance",
        )?;
        if accounting_residual_j.abs() > accounting_tolerance_j {
            return Err(TransientCarreauError::WorkBalance);
        }
        let mut reservoir_outflows = Vec::with_capacity(boundaries.len());
        for boundary in boundaries {
            checkpoint(cx)?;
            reservoir_outflows.push(VolumetricFlowRate::new(result.residuals[boundary.node]));
        }
        checkpoint(cx)?;
        Ok(TransientCarreauStep {
            duration_s,
            midpoint_pressures: pressures.iter().copied().map(Pressure::new).collect(),
            branches,
            node_outflows,
            storage_volume_changes_m3: result.storage_changes.clone(),
            reservoir_outflows,
            max_node_imbalance: VolumetricFlowRate::new(result.max_imbalance),
            max_final_rigid_node_imbalance: VolumetricFlowRate::new(
                result.max_final_rigid_imbalance,
            ),
            boundary_imbalance: VolumetricFlowRate::new(result.boundary_imbalance),
            reservoir_work_j,
            injection_work_j,
            kinetic_energy_change_j,
            storage_energy_change_j,
            viscous_heat_j,
            node_residual_work_j,
            radial_work_residual_j,
            storage_midpoint_work_residual_j,
            physical_work_residual_j,
            accounting_residual_j,
            accounting_tolerance_j,
            newton_iterations,
            branch_solves,
        })
    }
}

#[cfg(test)]
mod tests;
