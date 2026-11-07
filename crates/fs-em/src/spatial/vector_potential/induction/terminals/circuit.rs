//! Physical three-dimensional electrodes closed against the native RLC circuit.
//!
//! Every field trial advances the SAME accepted A-phi history. Its finite-step
//! terminal law is I = I0 + Y V, including diffusion, magnetic storage and
//! prescribed motion. Y is never stamped as resistor heat. The original native
//! descriptor retains all circuit storage, Kirchhoff constraints and sources.
//!
//! One positive reference resistor and one compensating native current port
//! join each non-reference electrode to its conductor's reference electrode.
//! The replayed injection j = V/Rref - I cancels that resistor exactly. With
//! the native response V = V0 + H j, only the small terminal matrix
//! [Id - H (Gref - Y)] is factored here. Reference loss is excluded from every
//! physical heat account; changing Rref cannot change the physical solution.
//!
//! Linear response matrices are cached only for the same duration, velocity,
//! fixed field/material/topology and numerical probe controls. The history
//! offsets and final native/field states are always recomputed. Probe fields
//! undergo the original material-domain and energy checks, as does final replay.
//! A source-pinned temperature refresh assembles endpoint operators once and
//! invalidates the material-dependent response. Every probe still starts from
//! the original magnetic effort and native circuit flux/charge. Geometry,
//! switching and nonlinear field evolution remain outside this owner.

mod native;
pub(in crate::spatial::vector_potential::induction) mod nonlinear;
pub use nonlinear::SpatialTerminalCircuitNonlinearControls;

use super::super::{Prepared, SpatialInductionControls, SpatialInductionStep};
use super::{SpatialTerminalInduction, SpatialTerminalInductionCandidate};
use crate::spatial::{SpatialMagneticError, SpatialMagneticMaterial};
use crate::winding::nonlinear::{NonlinearWindingError, native_seed};
use crate::winding::{
    nonlinear::resistors::GraphResistorMaterial, ohmic::WindingGraphResistorHeat,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_la::factor::{FactorError, Lu, lu};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_phs::{
    StepRecord,
    circuit::{
        Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
        assemble_circuit_with_mutual_inductance,
    },
};
use fs_solver::{NewtonError, NewtonReport, SolverRunProgress};
use std::rc::Rc;
use std::sync::Arc;

/// Native dense allocation and terminal probing are deliberately small.
pub const MAX_TERMINAL_CIRCUIT_COORDINATES: usize = 16;
pub const MAX_TERMINAL_CIRCUIT_STATES: usize = 128;
pub const MAX_TERMINAL_CIRCUIT_NODES: usize = 32;
pub const MAX_TERMINAL_CIRCUIT_BRANCHES: usize = 128;
pub const MAX_TERMINAL_CIRCUIT_SOURCE_PORTS: usize = 32;
/// Original source controls plus the private electrode-current controls.
pub const MAX_TERMINAL_CIRCUIT_PORTS: usize = 48;

/// Original circuit storage and source controls at the accepted start time.
#[derive(Clone, Debug)]
pub struct SpatialTerminalCircuitInitial {
    /// Original inductor-branch order, including mutual flux linkage [Wb].
    pub fluxes_wb: Vec<f64>,
    /// Original capacitor-branch order [C].
    pub charges_c: Vec<f64>,
    /// Original port indices, including zero-valued gaps.
    pub external_inputs: Vec<f64>,
}

/// Preparation bounds and a numerical impedance that cancels from physics.
#[derive(Clone, Copy, Debug)]
pub struct SpatialTerminalCircuitConfig {
    pub reference_resistance_ohm: f64,
    pub maximum_native_states: usize,
    /// Additional circuit vectors, terminal matrices and the cached velocity.
    /// The original field solver's memory bounds apply independently.
    pub maximum_interface_scalars: usize,
}

/// Explicit probing scales and physical final-replay allowances.
#[derive(Clone, Copy, Debug)]
pub struct SpatialTerminalCircuitControls {
    pub field: SpatialInductionControls,
    /// Positive voltage excursion for each actual native field response [V].
    pub voltage_probe_v: f64,
    /// Positive injected-current excursion for the native circuit response [A].
    pub current_probe_a: f64,
    /// A cold interval needs coordinates + 2 field trials; a cache hit needs 2.
    /// If every electrode is shorted to its conductor reference, one suffices.
    pub maximum_field_trials: usize,
    pub relative_terminal_tolerance: f64,
    pub voltage_tolerance_v: f64,
    pub current_tolerance_a: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum SpatialTerminalCircuitError {
    Invalid(&'static str),
    Field(SpatialMagneticError),
    NonlinearField(super::super::nonlinear::SpatialNonlinearInductionError),
    #[cfg(feature = "magnetic-hysteresis")]
    HystereticField(super::super::hysteretic::SpatialHystereticInductionError),
    Circuit(CircuitError),
    Initial(NonlinearWindingError),
    Factor(FactorError),
    Newton(NewtonError),
    Interrupted(SolverRunProgress),
    NotConverged(NewtonReport),
    FieldTrialBudget,
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for SpatialTerminalCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial terminal circuit: {self:?}")
    }
}
impl std::error::Error for SpatialTerminalCircuitError {}
impl From<SpatialMagneticError> for SpatialTerminalCircuitError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<CircuitError> for SpatialTerminalCircuitError {
    fn from(error: CircuitError) -> Self {
        Self::Circuit(error)
    }
}
impl From<super::super::nonlinear::SpatialNonlinearInductionError> for SpatialTerminalCircuitError {
    fn from(error: super::super::nonlinear::SpatialNonlinearInductionError) -> Self {
        use super::super::nonlinear::SpatialNonlinearInductionError;
        use crate::spatial::vector_potential::nonlinear::NonlinearVectorError;
        match error {
            SpatialNonlinearInductionError::Spatial(SpatialMagneticError::Cancelled)
            | SpatialNonlinearInductionError::Magnetic(NonlinearVectorError::Spatial(
                SpatialMagneticError::Cancelled,
            )) => Self::Cancelled,
            other => Self::NonlinearField(other),
        }
    }
}
#[cfg(feature = "magnetic-hysteresis")]
impl From<super::super::hysteretic::SpatialHystereticInductionError>
    for SpatialTerminalCircuitError
{
    fn from(error: super::super::hysteretic::SpatialHystereticInductionError) -> Self {
        use super::super::hysteretic::SpatialHystereticInductionError;
        use fs_material::field_storage::hysteresis::vector::VectorHysteresisError;
        match error {
            SpatialHystereticInductionError::Spatial(SpatialMagneticError::Cancelled)
            | SpatialHystereticInductionError::Material(VectorHysteresisError::Cancelled) => {
                Self::Cancelled
            }
            other => Self::HystereticField(other),
        }
    }
}
impl From<NewtonError> for SpatialTerminalCircuitError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl From<FactorError> for SpatialTerminalCircuitError {
    fn from(error: FactorError) -> Self {
        Self::Factor(error)
    }
}
type Result<T> = core::result::Result<T, SpatialTerminalCircuitError>;

struct Layout {
    graph: CircuitGraph,
    mutual: Vec<MutualInductance>,
    dae: CircuitDae,
    nodes: Vec<usize>,
    reference: Vec<usize>,
    coordinates: Vec<usize>,
    source_kinds: Vec<Option<bool>>,
    external_ports: usize,
    reference_ohm: f64,
}

struct Response {
    duration_s: f64,
    velocity_m_s: Vec<[f64; 3]>,
    controls: [u64; 11],
    field_s: Vec<f64>,
    native_ohm: Vec<f64>,
    matrix: Vec<f64>,
    row_scale: Vec<f64>,
    column_scale: Vec<f64>,
    factor: Lu,
}

/// Accepted field and original circuit storage; trials cannot mutate either.
#[derive(Clone)]
pub struct SpatialTerminalCircuit<'m, C: SpatialConductorMaterial = ResolvedConductorState> {
    field: SpatialTerminalInduction<'m, C>,
    layout: Rc<Layout>,
    state: Vec<f64>,
    response: Option<Rc<Response>>,
    graph_materials: Vec<GraphResistorMaterial>,
}

/// Actual first-to-second branch current and absorbed finite-step work.
#[derive(Clone, Debug)]
pub struct SpatialTerminalCircuitBranch {
    pub branch: usize,
    pub voltage_v: f64,
    pub current_a: f64,
    pub absorbed_work_j: f64,
}

pub struct SpatialTerminalCircuitStep<F = SpatialInductionStep> {
    /// Native diagnostic record, including the canceled numerical reference.
    /// Physical circuit heat is resistor_heat_j, not native.dissipated.
    pub native: StepRecord,
    pub field: F,
    pub node_voltages_v: Vec<f64>,
    /// Original electrode order; these are native circuit-node potentials.
    pub terminal_voltage_v: Vec<f64>,
    /// Positive currents enter the original conducting body.
    pub terminal_current_a: Vec<f64>,
    pub branches: Vec<SpatialTerminalCircuitBranch>,
    /// Only physical authored resistors; numerical references are excluded.
    pub graph_resistors: Vec<WindingGraphResistorHeat>,
    /// Original-branch material bindings at this accepted resistance update.
    /// Empty for a numerically authored resistor graph.
    pub graph_materials: Vec<GraphResistorMaterial>,
    pub external_source_outputs: Vec<f64>,
    /// Only original authored graph resistors, independently reconstructed.
    pub resistor_heat_j: f64,
    /// Numerical reference loss canceled by its private current port.
    pub cancelled_reference_heat_j: f64,
    pub circuit_stored_before_j: f64,
    pub circuit_stored_after_j: f64,
    /// Constitutive effort times original flux/charge increments.
    pub circuit_stored_change_j: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Stable native graph increment plus the original field's own increment.
    pub stored_change_j: f64,
    /// Physical field heat plus authored graph resistor heat, payable once.
    pub heat_j: f64,
    /// Field endpoint dissipation, excluded from every physical heat receipt.
    pub numerical_dissipation_j: f64,
    /// Original external circuit work plus impressed field work.
    pub supplied_work_j: f64,
    pub impressed_work_j: f64,
    pub parameter_work_j: f64,
    pub mechanical_work_j: f64,
    /// Original external sources only; private interface ports are excluded.
    pub external_work_j: f64,
    pub energy_residual_j: f64,
    pub relative_terminal_residual: f64,
    pub maximum_terminal_voltage_defect_v: f64,
    pub maximum_terminal_current_defect_a: f64,
    pub maximum_kirchhoff_current_defect_a: f64,
    pub field_trials: usize,
    pub reused_response: bool,
    /// Nonlinear electrode closure report; None for the direct affine path or
    /// a graph that shorts every electrode to its component reference.
    pub report: Option<NewtonReport>,
}

pub struct SpatialTerminalCircuitCandidate<'m, C: SpatialConductorMaterial = ResolvedConductorState>
{
    circuit: SpatialTerminalCircuit<'m, C>,
    step: SpatialTerminalCircuitStep,
}
impl<'m, C: SpatialConductorMaterial> SpatialTerminalCircuitCandidate<'m, C> {
    pub fn circuit(&self) -> &SpatialTerminalCircuit<'m, C> {
        &self.circuit
    }
    pub fn step(&self) -> &SpatialTerminalCircuitStep {
        &self.step
    }
    pub fn into_parts(self) -> (SpatialTerminalCircuit<'m, C>, SpatialTerminalCircuitStep) {
        (self.circuit, self.step)
    }
    /// Consistent mechanical transfer of this same current and magnetic stage.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: super::super::mechanical::SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<
        super::super::mechanical::SpatialLorentzLoads,
        super::super::mechanical::SpatialLorentzLoadError,
    > {
        let problem = self.circuit.field.problem();
        super::super::mechanical::recover(
            problem.positions_m(),
            &problem.mesh().tets,
            &self.step.field.cells,
            selected_cells,
            limits,
            cancelled,
        )
    }
}

fn invalid(message: &'static str) -> SpatialTerminalCircuitError {
    SpatialTerminalCircuitError::Invalid(message)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable terminal circuit arithmetic"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0 && value.recip().is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialTerminalCircuitError::Cancelled)
    } else {
        Ok(())
    }
}
fn close(defect: f64, tolerance: f64, message: &'static str) -> Result<()> {
    if finite(defect)?.abs() > tolerance {
        Err(SpatialTerminalCircuitError::Balance(message))
    } else {
        Ok(())
    }
}
fn signature(c: SpatialTerminalCircuitControls) -> [u64; 11] {
    let f = c.field;
    [
        f.relative_tolerance.to_bits(),
        f.absolute_residual_a.to_bits(),
        f.restart as u64,
        f.maximum_cycles as u64,
        f.maximum_krylov_scalars as u64,
        f.current_tolerance_a.to_bits(),
        f.flux_tolerance_wb.to_bits(),
        f.energy_tolerance_j.to_bits(),
        f.maximum_cell_magnetic_peclet.to_bits(),
        c.voltage_probe_v.to_bits(),
        c.current_probe_a.to_bits(),
    ]
}
fn admit_inputs(layout: &Layout, external: &[f64]) -> Result<()> {
    if external.len() != layout.external_ports
        || external
            .iter()
            .enumerate()
            .any(|(p, &v)| !v.is_finite() || (layout.source_kinds[p].is_none() && v != 0.0))
    {
        return Err(invalid(
            "complete original source inputs and zero port gaps required",
        ));
    }
    Ok(())
}

impl<'m, C: SpatialConductorMaterial> SpatialTerminalCircuit<'m, C> {
    /// Connect original electrode nodes to the original RLC/transformer graph.
    /// No graph inductor or resistor substitutes for field-owned storage/loss.
    /// Repeated node bindings represent electrical shorts between electrodes.
    /// Initial fluxes and charges are never projected by the native seed solve.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        field: SpatialTerminalInduction<'m, C>,
        graph: CircuitGraph,
        terminal_nodes: Vec<usize>,
        mutual: Vec<MutualInductance>,
        initial: SpatialTerminalCircuitInitial,
        config: SpatialTerminalCircuitConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let (layout, state) = native::prepare(
            field.terminals.references(),
            field.problem().mesh().vertex_count,
            graph,
            terminal_nodes,
            mutual,
            initial,
            config,
            cancelled,
        )?;
        Ok(Self {
            field,
            layout,
            state,
            response: None,
            graph_materials: Vec::new(),
        })
    }

    pub fn field(&self) -> &SpatialTerminalInduction<'m, C> {
        &self.field
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.layout.graph
    }
    /// Original-branch material bindings at the last resistance update.
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.graph_materials
    }

    /// Rebind every physical graph resistor to complete source-backed material
    /// and path geometry, retaining this exact field, native flux/charge and time.
    ///
    /// This explicit source reconfiguration can attach initial bindings or
    /// change resistance. Thermal consumers must separately pin the original
    /// sources and admit their temperature successors. Only dissipative graph
    /// coefficients change; no circuit storage or magnetic parameter work is
    /// created. Retained algebraic entries are the old receipt, and the next
    /// native descriptor solve determines the actual stage efforts.
    ///
    /// Identical sorted bindings retain the response cache. Rebinding clears
    /// it so changed resistance cannot reuse the previous circuit response.
    pub fn with_graph_resistors(
        &self,
        mut sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = self
            .layout
            .graph
            .branches
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Resistor { .. }))
            .count();
        if sources.len() != count || sources.len() > MAX_TERMINAL_CIRCUIT_BRANCHES {
            return Err(invalid(
                "every physical graph resistor requires one bounded source",
            ));
        }
        sources.sort_by_key(|source| source.branch);
        if sources == self.graph_materials {
            poll(cancelled)?;
            return Ok(self.clone());
        }
        let layout = native::rebind_resistors(&self.layout, &self.state, &sources, cancelled)?;
        let next = Self {
            field: self.field.clone(),
            layout,
            state: self.state.clone(),
            response: None,
            graph_materials: sources,
        };
        poll(cancelled)?;
        Ok(next)
    }
    pub fn mutual_inductance(&self) -> &[MutualInductance] {
        &self.layout.mutual
    }
    pub fn terminal_nodes(&self) -> &[usize] {
        &self.layout.nodes
    }
    pub fn native_state(&self) -> &[f64] {
        &self.state
    }
    pub fn time_s(&self) -> f64 {
        self.field.time_s()
    }

    fn native(
        &self,
        input: &[f64],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StepRecord> {
        native::advance(&self.layout, &self.state, input, dt, cancelled)
    }

    fn currents(&self, trial: &SpatialTerminalInductionCandidate<'_, C>) -> Vec<f64> {
        self.layout
            .coordinates
            .iter()
            .map(|&t| trial.step().terminals[t].incoming_current_a)
            .collect()
    }

    /// Propose one midpoint interval with original external port controls.
    /// All affine probes and both final replays share the accepted histories.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalCircuitCandidate<'m, C>> {
        self.step_prepared_candidate(
            Arc::clone(&self.field.field.prepared),
            external_inputs,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Close the original circuit against source-pinned endpoint conductivity
    /// and induced-linear permeability without resetting either history.
    ///
    /// Endpoint operators are prepared once and shared by all field trials.
    /// Every trial uses the accepted OLD magnetic source, A, impressed current,
    /// native flux and charge. Material refresh invalidates the cached response;
    /// independent source, charge, storage and parameter-work checks remain on
    /// the actual final replay. A finite caloric owner must separately account
    /// for signed magnetic parameter work and the field entropy correction.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        assembly: WeightedAssemblyLimits,
        controls: SpatialTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalCircuitCandidate<'m, C>> {
        let prepared = self
            .field
            .field
            .prepare_materials(materials, conductors, assembly, cancelled)?;
        self.step_prepared_candidate(
            prepared,
            external_inputs,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_prepared_candidate(
        &self,
        prepared: Arc<Prepared<'m, C>>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalCircuitCandidate<'m, C>> {
        poll(cancelled)?;
        admit_inputs(&self.layout, external_inputs)?;
        let k = self.layout.coordinates.len();
        if !positive(duration_s)
            || !positive(controls.voltage_probe_v)
            || !positive(controls.current_probe_a)
            || !positive(controls.relative_terminal_tolerance)
            || controls.relative_terminal_tolerance >= 1.0
            || !positive(controls.voltage_tolerance_v)
            || !positive(controls.current_tolerance_a)
            || !positive(controls.energy_tolerance_j)
            || velocity_m_s.len() != self.field.problem().mesh().vertex_count
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "positive probe scales, interval and physical replay allowances required",
            ));
        }
        if k == 0 {
            if controls.maximum_field_trials == 0 {
                return Err(SpatialTerminalCircuitError::FieldTrialBudget);
            }
            let voltage = vec![0.0; self.field.electrode_faces().len()];
            let field = self.field.step_voltage_prepared_candidate(
                Arc::clone(&prepared),
                &voltage,
                endpoint_impressed_current_a_m2,
                velocity_m_s,
                duration_s,
                controls.field,
                cancelled,
            )?;
            let native = self.native(external_inputs, duration_s, cancelled)?;
            let (field, field_step) = field.into_parts();
            let mut step = self.account(
                native,
                field_step,
                external_inputs,
                external_inputs,
                duration_s,
                controls,
                cancelled,
            )?;
            step.field_trials = 1;
            let state = step.native.x.clone();
            poll(cancelled)?;
            return Ok(SpatialTerminalCircuitCandidate {
                circuit: Self {
                    field,
                    layout: Rc::clone(&self.layout),
                    state,
                    response: None,
                    graph_materials: self.graph_materials.clone(),
                },
                step,
            });
        }
        let key = signature(controls);
        let cached = self.response.as_ref().filter(|r| {
            Arc::ptr_eq(&prepared, &self.field.field.prepared)
                && r.duration_s.to_bits() == duration_s.to_bits()
                && r.velocity_m_s == velocity_m_s
                && r.controls == key
        });
        let field_trials = if cached.is_some() { 2 } else { k + 2 };
        if field_trials > controls.maximum_field_trials {
            return Err(SpatialTerminalCircuitError::FieldTrialBudget);
        }
        let mut voltage = vec![0.0; self.field.electrode_faces().len()];
        let baseline = self.field.step_voltage_prepared_candidate(
            Arc::clone(&prepared),
            &voltage,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls.field,
            cancelled,
        )?;
        let offset = self.currents(&baseline);
        drop(baseline);
        let mut input = external_inputs.to_vec();
        input.resize(self.layout.dae.ports, 0.0);
        let open = self.native(&input, duration_s, cancelled)?;
        let open_voltage = &open.y[self.layout.external_ports..];
        let response = if let Some(cached) = cached {
            Rc::clone(cached)
        } else {
            let mut field_s = vec![0.0; k * k];
            let mut native_ohm = vec![0.0; k * k];
            for column in 0..k {
                poll(cancelled)?;
                let terminal = self.layout.coordinates[column];
                voltage[terminal] = controls.voltage_probe_v;
                let probe = self.field.step_voltage_prepared_candidate(
                    Arc::clone(&prepared),
                    &voltage,
                    endpoint_impressed_current_a_m2,
                    velocity_m_s,
                    duration_s,
                    controls.field,
                    cancelled,
                )?;
                for (row, current) in self.currents(&probe).into_iter().enumerate() {
                    field_s[row * k + column] =
                        finite((current - offset[row]) / controls.voltage_probe_v)?;
                }
                voltage[terminal] = 0.0;
                input[self.layout.external_ports + column] = controls.current_probe_a;
                let probe = self.native(&input, duration_s, cancelled)?;
                for row in 0..k {
                    native_ohm[row * k + column] = finite(
                        (probe.y[self.layout.external_ports + row] - open_voltage[row])
                            / controls.current_probe_a,
                    )?;
                }
                input[self.layout.external_ports + column] = 0.0;
            }
            let mut matrix = vec![0.0; k * k];
            for row in 0..k {
                poll(cancelled)?;
                for column in 0..k {
                    let mut value = if row == column { 1.0 } else { 0.0 };
                    value =
                        finite(value - native_ohm[row * k + column] / self.layout.reference_ohm)?;
                    for inner in 0..k {
                        value = finite(
                            native_ohm[row * k + inner].mul_add(field_s[inner * k + column], value),
                        )?;
                    }
                    matrix[row * k + column] = value;
                }
            }
            let (factor, row_scale, column_scale) = factor_terminal(&matrix, k, cancelled)?;
            Rc::new(Response {
                duration_s,
                velocity_m_s: velocity_m_s.to_vec(),
                controls: key,
                field_s,
                native_ohm,
                matrix,
                factor,
                row_scale,
                column_scale,
            })
        };
        let mut rhs = open_voltage.to_vec();
        for (row, value) in rhs.iter_mut().enumerate() {
            for (column, &current) in offset.iter().enumerate() {
                *value = finite((-response.native_ohm[row * k + column]).mul_add(current, *value))?;
            }
        }
        let (relative, relative_terminal_residual) = solve_terminal(&response, &rhs, cancelled)?;
        if relative_terminal_residual > controls.relative_terminal_tolerance {
            return Err(SpatialTerminalCircuitError::Balance(
                "original terminal Schur residual exceeds allowance",
            ));
        }
        for row in 0..k {
            let mut current = offset[row];
            for (column, &v) in relative.iter().enumerate() {
                current = finite(response.field_s[row * k + column].mul_add(v, current))?;
            }
            voltage[self.layout.coordinates[row]] = relative[row];
            input[self.layout.external_ports + row] =
                finite(relative[row] / self.layout.reference_ohm - current)?;
        }
        // Neither a cached response nor a small Schur residual admits a field:
        // replay all original constitutive, charge, Ampere and energy checks.
        let field = self.field.step_voltage_prepared_candidate(
            prepared,
            &voltage,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls.field,
            cancelled,
        )?;
        let native = self.native(&input, duration_s, cancelled)?;
        let (field, field_step) = field.into_parts();
        let mut step = self.account(
            native,
            field_step,
            &input,
            external_inputs,
            duration_s,
            controls,
            cancelled,
        )?;
        step.relative_terminal_residual = relative_terminal_residual;
        step.field_trials = field_trials;
        step.reused_response = cached.is_some();
        let state = step.native.x.clone();
        poll(cancelled)?;
        Ok(SpatialTerminalCircuitCandidate {
            circuit: Self {
                field,
                layout: Rc::clone(&self.layout),
                state,
                response: Some(response),
                graph_materials: self.graph_materials.clone(),
            },
            step,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn account(
        &self,
        native: StepRecord,
        field: SpatialInductionStep,
        input: &[f64],
        external: &[f64],
        dt: f64,
        controls: SpatialTerminalCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalCircuitStep> {
        native::account(
            &self.layout,
            &self.state,
            native,
            field,
            &self.graph_materials,
            input,
            external,
            dt,
            native::Allowances {
                voltage_v: controls.voltage_tolerance_v,
                current_a: controls.current_tolerance_a,
                energy_j: controls.energy_tolerance_j,
            },
            cancelled,
        )
    }
}

fn factor_terminal(
    matrix: &[f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Lu, Vec<f64>, Vec<f64>)> {
    let mut scaled = matrix.to_vec();
    let mut rows = Vec::with_capacity(n);
    for row in scaled.chunks_exact_mut(n) {
        poll(cancelled)?;
        let scale = row.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if !positive(scale) {
            return Err(invalid("terminal row has no representable response"));
        }
        for value in row {
            *value = finite(*value / scale)?;
        }
        rows.push(scale);
    }
    let mut columns = vec![0.0; n];
    for column in 0..n {
        poll(cancelled)?;
        let scale = (0..n)
            .map(|r| scaled[r * n + column].abs())
            .fold(0.0_f64, f64::max);
        if !positive(scale) {
            return Err(invalid("terminal column has no representable response"));
        }
        columns[column] = scale;
        for row in 0..n {
            scaled[row * n + column] = finite(scaled[row * n + column] / scale)?;
        }
    }
    let factor = lu(&scaled, n)?;
    poll(cancelled)?;
    Ok((factor, rows, columns))
}

fn solve_terminal(
    response: &Response,
    rhs: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, f64)> {
    let mut solution = rhs
        .iter()
        .zip(&response.row_scale)
        .map(|(&v, &scale)| finite(v / scale))
        .collect::<Result<Vec<_>>>()?;
    poll(cancelled)?;
    response.factor.solve(&mut solution);
    poll(cancelled)?;
    for (value, &scale) in solution.iter_mut().zip(&response.column_scale) {
        *value = finite(*value / scale)?;
    }
    let n = rhs.len();
    let mut residual = 0.0_f64;
    for row in 0..n {
        poll(cancelled)?;
        let (mut defect, mut scale) = (-rhs[row], rhs[row].abs());
        for (column, &value) in solution.iter().enumerate() {
            let coefficient = response.matrix[row * n + column];
            defect = finite(coefficient.mul_add(value, defect))?;
            scale = finite(scale + coefficient.abs() * value.abs())?;
        }
        if scale > 0.0 {
            residual = residual.max(finite(defect.abs() / scale)?);
        }
    }
    Ok((solution, residual))
}
