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
//! Geometry, material refresh, switching and nonlinear field evolution are
//! outside this fixed-configuration owner.

use super::super::{SpatialInductionControls, SpatialInductionStep};
use super::{SpatialTerminalInduction, SpatialTerminalInductionCandidate};
use crate::spatial::SpatialMagneticError;
use crate::winding::nonlinear::{NonlinearWindingError, native_seed};
use fs_la::factor::{FactorError, Lu, lu};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_phs::{
    StepRecord,
    circuit::{
        Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
        assemble_circuit_with_mutual_inductance,
    },
};
use std::rc::Rc;

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
    Circuit(CircuitError),
    Initial(NonlinearWindingError),
    Factor(FactorError),
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
}

/// Actual first-to-second branch current and absorbed finite-step work.
#[derive(Clone, Debug)]
pub struct SpatialTerminalCircuitBranch {
    pub branch: usize,
    pub voltage_v: f64,
    pub current_a: f64,
    pub absorbed_work_j: f64,
}

pub struct SpatialTerminalCircuitStep {
    /// Native diagnostic record, including the canceled numerical reference.
    /// Physical circuit heat is resistor_heat_j, not native.dissipated.
    pub native: StepRecord,
    pub field: SpatialInductionStep,
    pub node_voltages_v: Vec<f64>,
    /// Original electrode order; these are native circuit-node potentials.
    pub terminal_voltage_v: Vec<f64>,
    /// Positive currents enter the original conducting body.
    pub terminal_current_a: Vec<f64>,
    pub branches: Vec<SpatialTerminalCircuitBranch>,
    pub external_source_outputs: Vec<f64>,
    /// Only original authored graph resistors, independently reconstructed.
    pub resistor_heat_j: f64,
    /// Numerical reference loss canceled by its private current port.
    pub cancelled_reference_heat_j: f64,
    pub circuit_stored_before_j: f64,
    pub circuit_stored_after_j: f64,
    /// Constitutive effort times original flux/charge increments.
    pub circuit_stored_change_j: f64,
    /// Original external sources only; private interface ports are excluded.
    pub external_work_j: f64,
    pub energy_residual_j: f64,
    pub relative_terminal_residual: f64,
    pub maximum_terminal_voltage_defect_v: f64,
    pub maximum_terminal_current_defect_a: f64,
    pub maximum_kirchhoff_current_defect_a: f64,
    pub field_trials: usize,
    pub reused_response: bool,
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
        poll(cancelled)?;
        let count = field.electrode_faces().len();
        if terminal_nodes.len() != count || terminal_nodes.iter().any(|&n| n >= graph.node_count) {
            return Err(invalid("complete original electrode-to-node map required"));
        }
        let coordinates: Vec<_> = field
            .terminals
            .reference
            .iter()
            .enumerate()
            .filter_map(|(t, &reference)| {
                // A real wire to the same node fixes zero relative voltage.
                // It needs neither a reference resistor nor a self-loop source.
                (terminal_nodes[t] != terminal_nodes[reference]).then_some(t)
            })
            .collect();
        let k = coordinates.len();
        if k > MAX_TERMINAL_CIRCUIT_COORDINATES
            || !(2..=MAX_TERMINAL_CIRCUIT_NODES).contains(&graph.node_count)
            || graph.branches.len().saturating_add(2 * k) > MAX_TERMINAL_CIRCUIT_BRANCHES
            || graph.transformers.len() > 16
            || mutual.len() > graph.branches.len().saturating_mul(graph.branches.len()) / 2
            || !positive(config.reference_resistance_ohm)
        {
            return Err(invalid(
                "bounded native graph, complete electrode nodes and positive reference required",
            ));
        }
        let mut source_kinds = vec![None; MAX_TERMINAL_CIRCUIT_SOURCE_PORTS];
        let (mut external_ports, mut inductors, mut capacitors, mut voltage_sources) = (0, 0, 0, 0);
        for &(_, _, branch) in &graph.branches {
            poll(cancelled)?;
            let (port, kind) = match branch {
                Branch::Resistor { ohms } => {
                    if !positive(ohms) {
                        return Err(invalid("representable physical resistance required"));
                    }
                    continue;
                }
                Branch::Inductor { henries } => {
                    if !positive(henries) {
                        return Err(invalid("representable physical inductance required"));
                    }
                    inductors += 1;
                    continue;
                }
                Branch::Capacitor { farads } => {
                    if !positive(farads) {
                        return Err(invalid("representable physical capacitance required"));
                    }
                    capacitors += 1;
                    continue;
                }
                Branch::VoltageSource { port } => {
                    voltage_sources += 1;
                    (port, true)
                }
                Branch::CurrentSource { port } => (port, false),
            };
            if port >= source_kinds.len() || source_kinds[port].is_some_and(|old| old != kind) {
                return Err(invalid(
                    "source port exceeds capacity or mixes voltage/current units",
                ));
            }
            source_kinds[port] = Some(kind);
            external_ports = external_ports.max(port + 1);
        }
        let n = inductors + 2 * capacitors + graph.node_count - 1
            + voltage_sources
            + graph.transformers.len();
        let payload =
            16 * k * k + 12 * k + 8 * n + 6 * graph.branches.len() + 6 * count + 4 * external_ports;
        let payload = field
            .problem()
            .mesh()
            .vertex_count
            .checked_mul(3)
            .and_then(|v| v.checked_add(payload));
        if external_ports
            .checked_add(k)
            .is_none_or(|ports| ports > MAX_TERMINAL_CIRCUIT_PORTS)
            || n > MAX_TERMINAL_CIRCUIT_STATES
            || n > config.maximum_native_states
            || payload.is_none_or(|v| v > config.maximum_interface_scalars)
            || initial.fluxes_wb.len() != inductors
            || initial.charges_c.len() != capacitors
            || initial
                .fluxes_wb
                .iter()
                .chain(&initial.charges_c)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "bounded native states and complete original circuit storage required",
            ));
        }
        let mut lowered = graph.clone();
        for (p, &terminal) in coordinates.iter().enumerate() {
            let a = terminal_nodes[terminal];
            let b = terminal_nodes[field.terminals.reference[terminal]];
            lowered.branches.push((
                a,
                b,
                Branch::Resistor {
                    ohms: config.reference_resistance_ohm,
                },
            ));
            lowered.branches.push((
                a,
                b,
                Branch::CurrentSource {
                    port: external_ports + p,
                },
            ));
        }
        poll(cancelled)?;
        let dae = assemble_circuit_with_mutual_inductance(&lowered, &mutual)?;
        poll(cancelled)?;
        if dae.ports != external_ports + k {
            return Err(invalid(
                "native assembly changed the original/private port map",
            ));
        }
        let layout = Rc::new(Layout {
            graph,
            mutual,
            dae,
            nodes: terminal_nodes,
            coordinates,
            source_kinds,
            external_ports,
            reference_ohm: config.reference_resistance_ohm,
        });
        admit_inputs(&layout, &initial.external_inputs)?;
        let mut state = vec![0.0; layout.dae.system.state_dim()];
        for (&index, &flux) in layout.dae.flux_index.iter().zip(&initial.fluxes_wb) {
            state[index] = flux;
        }
        for (&index, &charge) in layout.dae.charge_index.iter().zip(&initial.charges_c) {
            state[index] = charge;
        }
        let mut input = initial.external_inputs;
        input.resize(layout.dae.ports, 0.0);
        let state = native_seed(&layout.dae, &state, &input)
            .map_err(SpatialTerminalCircuitError::Initial)?;
        poll(cancelled)?;
        Ok(Self {
            field,
            layout,
            state,
            response: None,
        })
    }

    pub fn field(&self) -> &SpatialTerminalInduction<'m, C> {
        &self.field
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.layout.graph
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
        poll(cancelled)?;
        // A finite-stage voltage is not a new initial condition. Preserve all
        // accepted storage and let the original descriptor solve its multipliers.
        let (step, _) = self.layout.dae.step_audited(&self.state, input, dt)?;
        poll(cancelled)?;
        Ok(step)
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
            let field = self.field.step_voltage_candidate(
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
                },
                step,
            });
        }
        let key = signature(controls);
        let cached = self.response.as_ref().filter(|r| {
            r.duration_s.to_bits() == duration_s.to_bits()
                && r.velocity_m_s == velocity_m_s
                && r.controls == key
        });
        let field_trials = if cached.is_some() { 2 } else { k + 2 };
        if field_trials > controls.maximum_field_trials {
            return Err(SpatialTerminalCircuitError::FieldTrialBudget);
        }
        let mut voltage = vec![0.0; self.field.electrode_faces().len()];
        let baseline = self.field.step_voltage_candidate(
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
                let probe = self.field.step_voltage_candidate(
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
        let field = self.field.step_voltage_candidate(
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
            },
            step,
        })
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
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
        let layout = &self.layout;
        let dae = &layout.dae;
        let midpoint: Vec<_> = self
            .state
            .iter()
            .zip(&native.x)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let mut nodes = vec![0.0; layout.graph.node_count];
        for (node, &index) in dae.node_potential_index.iter().enumerate() {
            nodes[node + 1] = finite(midpoint[index])?;
        }
        let mut maximum_voltage = 0.0_f64;
        let mut maximum_current = 0.0_f64;
        let mut reference_heat = 0.0;
        let mut interface_work = 0.0;
        for (p, &terminal) in layout.coordinates.iter().enumerate() {
            poll(cancelled)?;
            let reference = self.field.terminals.reference[terminal];
            let v = finite(nodes[layout.nodes[terminal]] - nodes[layout.nodes[reference]])?;
            let current = finite(v / layout.reference_ohm - input[layout.external_ports + p])?;
            maximum_voltage = maximum_voltage
                .max(finite(v - field.terminals[terminal].voltage_to_reference_v)?.abs());
            maximum_current = maximum_current
                .max(finite(current - field.terminals[terminal].incoming_current_a)?.abs());
            close(
                v - native.y[layout.external_ports + p],
                controls.voltage_tolerance_v,
                "native private port and original node voltage disagree",
            )?;
            reference_heat = finite(reference_heat + dt * v * (v / layout.reference_ohm))?;
            interface_work = finite(
                interface_work
                    + dt * input[layout.external_ports + p] * native.y[layout.external_ports + p],
            )?;
        }
        close(
            maximum_voltage,
            controls.voltage_tolerance_v,
            "native circuit and field electrode voltages disagree",
        )?;
        close(
            maximum_current,
            controls.current_tolerance_a,
            "native circuit and field electrode currents disagree",
        )?;
        close(
            interface_work + field.terminal_work_j - reference_heat,
            controls.energy_tolerance_j,
            "private circuit work does not cancel physical terminal work and reference loss",
        )?;

        let currents = dae.inductor_currents(&midpoint)?;
        let (mut l, mut c, mut source) = (0, 0, 0);
        let source_start =
            dae.system.differential_dim() + layout.graph.node_count - 1 + dae.charge_index.len();
        let mut kcl = vec![0.0; layout.graph.node_count];
        let mut branches = Vec::with_capacity(layout.graph.branches.len());
        let (mut resistor_heat, mut stored_change) = (0.0, 0.0);
        for (branch, &(a, b, kind)) in layout.graph.branches.iter().enumerate() {
            poll(cancelled)?;
            let v = finite(nodes[a] - nodes[b])?;
            let current = match kind {
                Branch::Resistor { ohms } => finite(v / ohms)?,
                Branch::Inductor { .. } => {
                    let index = dae.flux_index[l];
                    let current = currents[l];
                    stored_change = finite(
                        stored_change + current * finite(native.x[index] - self.state[index])?,
                    )?;
                    close(
                        finite(native.x[index] - self.state[index])? - dt * v,
                        dt * controls.voltage_tolerance_v,
                        "original inductor flux increment violates terminal voltage",
                    )?;
                    l += 1;
                    current
                }
                Branch::Capacitor { farads } => {
                    let index = dae.charge_index[c];
                    let dq = finite(native.x[index] - self.state[index])?;
                    let effort = finite(midpoint[index] / farads)?;
                    close(
                        v - effort,
                        controls.voltage_tolerance_v,
                        "original capacitor charge violates terminal voltage",
                    )?;
                    stored_change = finite(stored_change + effort * dq)?;
                    c += 1;
                    finite(dq / dt)?
                }
                Branch::VoltageSource { port } => {
                    close(
                        v - external[port],
                        controls.voltage_tolerance_v,
                        "original voltage source does not hold its authored effort",
                    )?;
                    let value = finite(midpoint[source_start + source])?;
                    source += 1;
                    value
                }
                Branch::CurrentSource { port } => -external[port],
            };
            let work = finite(dt * v * current)?;
            if matches!(kind, Branch::Resistor { .. }) {
                resistor_heat = finite(resistor_heat + work)?;
            }
            kcl[a] = finite(kcl[a] + current)?;
            kcl[b] = finite(kcl[b] - current)?;
            branches.push(SpatialTerminalCircuitBranch {
                branch,
                voltage_v: v,
                current_a: current,
                absorbed_work_j: work,
            });
        }
        let transformer_start = dae.system.state_dim() - layout.graph.transformers.len();
        for (index, transformer) in layout.graph.transformers.iter().enumerate() {
            poll(cancelled)?;
            let primary = finite(midpoint[transformer_start + index])?;
            let secondary = finite(-transformer.ratio * primary)?;
            let (a, b) = transformer.primary;
            let (c, d) = transformer.secondary;
            let vp = finite(nodes[a] - nodes[b])?;
            let vs = finite(nodes[c] - nodes[d])?;
            close(
                vp - transformer.ratio * vs,
                controls.voltage_tolerance_v,
                "original ideal transformer voltage constraint failed",
            )?;
            close(
                dt * (vp * primary + vs * secondary),
                controls.energy_tolerance_j,
                "original ideal transformer creates work",
            )?;
            for ((a, b), current) in [
                (transformer.primary, primary),
                (transformer.secondary, secondary),
            ] {
                kcl[a] = finite(kcl[a] + current)?;
                kcl[b] = finite(kcl[b] - current)?;
            }
        }
        let mut component_balance = std::collections::BTreeMap::new();
        for (terminal, transfer) in field.terminals.iter().enumerate() {
            let node = layout.nodes[terminal];
            kcl[node] = finite(kcl[node] + transfer.incoming_current_a)?;
            let balance = component_balance.entry(transfer.component).or_insert(0.0);
            *balance = finite(*balance + transfer.incoming_current_a)?;
        }
        for balance in component_balance.into_values() {
            close(
                balance,
                controls.current_tolerance_a,
                "physical conductor electrode currents do not balance",
            )?;
        }
        let maximum_kcl = kcl.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        close(
            maximum_kcl,
            controls.current_tolerance_a,
            "original graph and physical electrodes violate Kirchhoff current balance",
        )?;
        let mut external_work = 0.0;
        for (&u, &y) in external.iter().zip(&native.y) {
            external_work = finite(external_work + dt * u * y)?;
        }
        close(
            native.delta_h - stored_change,
            controls.energy_tolerance_j,
            "native and original constitutive storage changes disagree",
        )?;
        close(
            native.dissipated - resistor_heat - reference_heat,
            controls.energy_tolerance_j,
            "physical resistor and canceled reference heat do not reproduce native dissipation",
        )?;
        close(
            native.supplied - external_work - interface_work,
            controls.energy_tolerance_j,
            "original external source and private interface work disagree with native supply",
        )?;
        let energy_residual = finite(
            field.stored_after_j - field.stored_before_j
                + stored_change
                + field.heat_j
                + resistor_heat
                + field.mechanical_work_j
                - field.impressed_work_j
                - field.parameter_work_j
                - external_work,
        )?;
        close(
            energy_residual,
            controls.energy_tolerance_j,
            "coupled circuit/field storage, Ohmic heat, Lorentz work and external supply do not close",
        )?;
        let stored_before = finite(dae.system.hamiltonian(&self.state))?;
        let stored_after = finite(dae.system.hamiltonian(&native.x))?;
        let terminal_voltage_v = layout.nodes.iter().map(|&n| nodes[n]).collect();
        let terminal_current_a = field
            .terminals
            .iter()
            .map(|t| t.incoming_current_a)
            .collect();
        let external_source_outputs = native.y[..layout.external_ports].to_vec();
        Ok(SpatialTerminalCircuitStep {
            native,
            field,
            node_voltages_v: nodes,
            terminal_voltage_v,
            terminal_current_a,
            branches,
            external_source_outputs,
            resistor_heat_j: resistor_heat,
            cancelled_reference_heat_j: reference_heat,
            circuit_stored_before_j: stored_before,
            circuit_stored_after_j: stored_after,
            circuit_stored_change_j: stored_change,
            external_work_j: external_work,
            energy_residual_j: energy_residual,
            relative_terminal_residual: 0.0,
            maximum_terminal_voltage_defect_v: maximum_voltage,
            maximum_terminal_current_defect_a: maximum_current,
            maximum_kirchhoff_current_defect_a: maximum_kcl,
            field_trials: 0,
            reused_response: false,
        })
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
