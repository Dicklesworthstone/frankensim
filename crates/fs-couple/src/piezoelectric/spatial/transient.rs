//! Consistent-mass spatial piezoelectric motion and retained material memory.
//!
//! The native solid mass participates only in mechanical momentum balance;
//! electrostatics remains a charge constraint. Implicit midpoint solves one
//! coupled increment with `K + 4 M / dt^2`, the original reciprocal `G`, and
//! `-K_e`. This resolves the retained mesh's elastic waves without replacing
//! the material by a fitted actuator gain or a single oscillator.
//!
//! Support coordinates and electrode topology/types are fixed. Supports may
//! follow explicitly supplied endpoint velocities through the same midpoint
//! kinematics. Their full inertial/constitutive reactions supply boundary work.
//! Prescribed voltages and
//! charges vary linearly between accepted endpoints; force is the caller's
//! midpoint value. The accepted energy account is kinetic + elastic + electric
//! storage plus independently recovered material heat against
//! `F_mid . delta_u + phi_mid . delta_Q`. Optional source-bound Maxwell and
//! Debye memories share the same midpoint solve. Explicit acoustic guide ports
//! retain signed wave outflow separately from ceramic heat. There is no numerical
//! damping, inferred loss spectrum, thermal update, or automatic
//! time/mesh-error certificate. Small-strain geometry and material coefficients
//! stay fixed. A candidate never mutates the accepted body.

use super::{
    BoundaryMap, ElectrodeCondition, PiezoelectricElectrode, PiezoelectricElectrodeResponse,
    Result, TetPiezoelectricCell, TetPiezoelectricControls, TetPiezoelectricError,
    TetPiezoelectricModel, TetPiezoelectricSolution, equilibrated_solve, invalid, norm, poll,
    positive, prepare_boundary, validate_controls,
};
use fs_exec::Cx;
use fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric;
use fs_math::dd::Dd;
use fs_math::quadratic::{self, QuadraticError};
use std::sync::Arc;

mod energy;
use energy::Products;

/// Reciprocal guide loading and real incident acoustic waves on this solve.
pub mod acoustic;
/// Source-bound resistive loads and prescribed electrode currents on this solve.
pub mod load;
/// Source-bound intrinsic ceramic memories on the original spatial operators.
pub mod relaxation;

struct FixedBoundary {
    map: BoundaryMap,
    electrode_nodes: Vec<Vec<usize>>,
}

/// One accepted spatial state. Fields and velocities are at the SAME endpoint.
#[derive(Clone, Debug)]
pub struct TetPiezoelectricTransientState {
    time_s: f64,
    field: TetPiezoelectricSolution,
    velocity_m_s: Vec<[f64; 3]>,
    // Retain differences before adding any large absolute voltage gauge.
    relative_voltage_v: Vec<f64>,
    electrode_conditions: Vec<ElectrodeCondition>,
    nodal_free_charge_c: Vec<f64>,
    kinetic_energy_j: f64,
    material_memory_energy_j: f64,
}

impl TetPiezoelectricTransientState {
    /// Accepted time measured from initialization [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Original mesh-node displacement [m].
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.field.displacement_m
    }
    /// Original mesh-node velocity [m/s], including prescribed support motion.
    pub fn velocity_m_s(&self) -> &[[f64; 3]] {
        &self.velocity_m_s
    }
    /// Absolute nodal potentials [V]; fields and work use retained differences.
    pub fn voltage_v(&self) -> &[f64] {
        &self.field.voltage_v
    }
    /// Original cell strains, stresses, electric fields and energy.
    pub fn cells(&self) -> &[TetPiezoelectricCell] {
        &self.field.cells
    }
    /// Solved total electrode free charge and potential, in original order.
    pub fn electrodes(&self) -> &[PiezoelectricElectrodeResponse] {
        &self.field.electrodes
    }
    /// Boundary nodal charge, excluding explicitly supplied volume charge [C].
    pub fn boundary_charge_c(&self) -> &[f64] {
        &self.field.boundary_charge_c
    }
    /// Kinetic energy from the original consistent mass [J].
    pub const fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_energy_j
    }
    /// Retained mechanical and dielectric memory storage [J]. This is already
    /// included in `stored_energy_j()` and must not be added a second time.
    pub const fn material_memory_energy_j(&self) -> f64 {
        self.material_memory_energy_j
    }
    /// Positive internal elastic, electric and material memory storage [J].
    pub fn stored_energy_j(&self) -> f64 {
        self.field.stored_energy_j
    }
    /// Kinetic plus all internal material storage [J].
    pub fn total_energy_j(&self) -> f64 {
        self.kinetic_energy_j + self.field.stored_energy_j
    }
}

/// Native spatial model, immutable boundary topology, and accepted state.
/// No static support condition is imposed: consistent mass admits rigid motion.
#[derive(Clone)]
pub struct TetPiezoelectricTransient<'model, 'mesh> {
    model: &'model TetPiezoelectricModel<'mesh>,
    boundary: Arc<FixedBoundary>,
    state: TetPiezoelectricTransientState,
    intrinsic: Option<relaxation::Relaxation>,
}

/// Independent physical account for one midpoint window.
#[derive(Clone, Copy, Debug)]
pub struct TetPiezoelectricTransientStep {
    /// Accepted time interval [s].
    pub time_step_s: f64,
    /// Applied-force work plus work supplied by moving supports [J].
    pub mechanical_work_j: f64,
    /// Applied midpoint nodal forces dotted with actual displacement change [J].
    pub applied_force_work_j: f64,
    /// Full midpoint support reactions dotted with support displacement change [J].
    /// These reactions include consistent-mass inertia and material memory.
    pub support_work_j: f64,
    /// Midpoint potential dotted with actual free-charge change [J].
    pub electrical_work_j: f64,
    /// Signed work entering attached acoustic guides [J], negative when incoming
    /// sound supplies the body. Outgoing wave energy is never ceramic heat.
    pub acoustic_outward_work_j: f64,
    /// Consistent-mass kinetic energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Elastic, electric and retained material memory storage change [J].
    pub stored_energy_change_j: f64,
    /// Nonnegative intrinsic mechanical and dielectric heat in the ceramic [J].
    pub material_heat_j: f64,
    /// Nonnegative intrinsic Maxwell relaxation heat [J].
    pub mechanical_material_heat_j: f64,
    /// Nonnegative intrinsic Debye polarization heat [J].
    pub dielectric_material_heat_j: f64,
    /// Retained material memory storage before and after the step [J]. This
    /// subaccount is already included in the body's internal stored energy.
    pub material_memory_energy_j: [f64; 2],
    /// Total energy change plus material heat and signed acoustic outflow,
    /// minus mechanical and electrical supplies [J].
    pub energy_residual_j: f64,
    /// Independently reconstructed midpoint momentum residual [N].
    pub maximum_force_residual_n: f64,
    /// Largest midpoint/endpoint free-node or electrode-total residual [C].
    pub maximum_charge_residual_c: f64,
    /// Native MINRES iteration count for this window.
    pub iterations: usize,
    /// Independently recomputed equilibrated increment-system residual.
    pub relative_residual: f64,
    // Preserve the original independent storage and source contractions until
    // enclosing electrical and thermal owners finish their own work balance.
    kinetic_change: Dd,
    energy_change: Dd,
    mechanical_work: Dd,
}

impl TetPiezoelectricTransientStep {
    /// Original consistent-mass increment before diagnostic rounding [J].
    pub const fn retained_kinetic_energy_change_j(&self) -> Dd {
        self.kinetic_change
    }
    /// Kinetic plus material storage increment before diagnostic rounding [J].
    /// This is computed from the actual endpoint states, independently of work.
    pub const fn retained_energy_change_j(&self) -> Dd {
        self.energy_change
    }
    /// Applied-force and moving-support work before diagnostic rounding [J].
    pub const fn retained_mechanical_work_j(&self) -> Dd {
        self.mechanical_work
    }
}

/// Joint mechanical/electrical candidate; accepted state remains untouched.
#[derive(Clone)]
pub struct TetPiezoelectricTransientCandidate<'model, 'mesh> {
    body: TetPiezoelectricTransient<'model, 'mesh>,
    step: TetPiezoelectricTransientStep,
    midpoint_support_force_n: Vec<[f64; 3]>,
    acoustic_ports: Vec<acoustic::TransientAcousticPortResponse>,
}

impl<'model, 'mesh> TetPiezoelectricTransientCandidate<'model, 'mesh> {
    /// Candidate body after every physical audit has passed.
    pub const fn body(&self) -> &TetPiezoelectricTransient<'model, 'mesh> {
        &self.body
    }
    /// Independent work and residual report.
    pub const fn step(&self) -> &TetPiezoelectricTransientStep {
        &self.step
    }
    /// Interval-mean force exerted BY each prescribed support ON the body [N].
    /// Free coordinates are zero. This is the complete dynamic reaction,
    /// including inertia, rather than the endpoint elastic internal force.
    pub fn midpoint_support_force_n(&self) -> &[[f64; 3]] {
        &self.midpoint_support_force_n
    }
    /// Incoming/outgoing waves and physical work at each original guide port.
    pub fn acoustic_ports(&self) -> &[acoustic::TransientAcousticPortResponse] {
        &self.acoustic_ports
    }
    /// Transfer the complete accepted candidate to its owner.
    pub fn into_parts(
        self,
    ) -> (
        TetPiezoelectricTransient<'model, 'mesh>,
        TetPiezoelectricTransientStep,
    ) {
        (self.body, self.step)
    }
}

impl<'model, 'mesh> TetPiezoelectricTransient<'model, 'mesh> {
    /// Initialize prescribed displacement/velocity and solve ONLY the matching
    /// electrostatics. Initial mechanics need not be in static equilibrium.
    /// Fixed support displacements must match exactly and their velocity is zero.
    /// `Charge(0)` is an isolated, uncharged electrode, with no voltage projection.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: &'model TetPiezoelectricModel<'mesh>,
        prescribed_displacements_m: &[(usize, f64)],
        electrodes: &[PiezoelectricElectrode<'_>],
        initial_displacement_m: &[[f64; 3]],
        initial_velocity_m_s: &[[f64; 3]],
        initial_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        validate_controls(controls)?;
        let n = model.nodes_m.len();
        if initial_displacement_m.len() != n
            || initial_velocity_m_s.len() != n
            || initial_nodal_free_charge_c.len() != n
            || initial_displacement_m
                .iter()
                .chain(initial_velocity_m_s)
                .flatten()
                .chain(initial_nodal_free_charge_c)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite initial displacement, velocity and charge required",
            ));
        }
        let zero_force = vec![[0.0; 3]; n];
        let map = prepare_boundary(
            model.nodes_m,
            model.tetrahedra,
            prescribed_displacements_m,
            electrodes,
            &zero_force,
            initial_nodal_free_charge_c,
            cx,
        )?;
        for &(dof, value) in prescribed_displacements_m {
            poll(cx)?;
            if initial_displacement_m[dof / 3][dof % 3] != value
                || initial_velocity_m_s[dof / 3][dof % 3] != 0.0
            {
                return Err(invalid(
                    "initial displacement and velocity must satisfy fixed supports",
                ));
            }
        }
        let boundary = Arc::new(FixedBoundary {
            map,
            electrode_nodes: electrodes.iter().map(|e| e.nodes.to_vec()).collect(),
        });
        let conditions = electrodes.iter().map(|e| e.condition).collect::<Vec<_>>();
        let (normalized, references) = normalize_conditions(&boundary, &conditions, cx)?;
        let normalized_electrodes = electrode_views(&boundary, &normalized);
        let mut prescribed = Vec::with_capacity(3 * n);
        for (node, displacement) in initial_displacement_m.iter().enumerate() {
            poll(cx)?;
            for (axis, &value) in displacement.iter().enumerate() {
                prescribed.push((3 * node + axis, value));
            }
        }
        let field = model.solve(
            &prescribed,
            &normalized_electrodes,
            &zero_force,
            initial_nodal_free_charge_c,
            controls,
            cx,
        )?;
        let kinetic_energy_j = kinetic_energy(model, initial_velocity_m_s, cx)?;
        let state = make_state(
            field,
            initial_velocity_m_s.to_vec(),
            conditions,
            initial_nodal_free_charge_c.to_vec(),
            kinetic_energy_j,
            0.0,
            &references,
            &boundary,
            cx,
        )?;
        Ok(Self {
            model,
            boundary,
            state,
            intrinsic: None,
        })
    }

    /// Attach zero intrinsic memories to a newly initialized body. The base
    /// tensor, source identity and mesh frame must match every original cell.
    /// The base stiffness is relaxed; its permittivity is instantaneous/clamped.
    /// Memory starts at zero even when initial strain or field is nonzero.
    /// No history may be reset after time has advanced. `maximum_internal_variables`
    /// explicitly bounds all retained cell scalars before any state allocation.
    /// Frequency-limited source claims cannot authorize arbitrary transients.
    /// Nonzero thermal bias requires a separately defined memory strain and
    /// polarization reference and currently refuses.
    pub fn with_relaxation(
        mut self,
        materials: Vec<RelaxingTensorPiezoelectric>,
        maximum_internal_variables: usize,
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        validate_controls(controls)?;
        if self.model.has_thermal_bias() {
            return Err(invalid(
                "intrinsic memories require a declared thermal strain and polarization reference",
            ));
        }
        if self.intrinsic.is_some() || self.state.time_s != 0.0 {
            return Err(invalid(
                "intrinsic memories can only be attached once at initialization",
            ));
        }
        for cell in &self.state.field.cells {
            poll(cx)?;
            if norm(&cell.strain_mandel) > controls.maximum_strain_norm
                || norm(&cell.response.electric_field_v_m) > controls.maximum_electric_field_v_m
            {
                return Err(invalid(
                    "initial state leaves the requested small-signal envelope",
                ));
            }
        }
        self.intrinsic = Some(relaxation::Relaxation::new(
            self.model,
            materials,
            maximum_internal_variables,
            cx,
        )?);
        poll(cx)?;
        Ok(self)
    }

    /// Latest intrinsic heat and memory-storage accounts, in original cell
    /// order. Newly initialized memories report zeros; an unextended body is None.
    pub fn material_relaxation_steps(&self) -> Option<&[relaxation::CellRelaxationStep]> {
        self.intrinsic.as_ref().map(|memory| memory.steps())
    }

    /// Accepted state; candidate creation never changes it.
    pub const fn state(&self) -> &TetPiezoelectricTransientState {
        &self.state
    }

    /// Advance one midpoint window, retaining any attached material memories.
    /// Source types and topology must
    /// remain unchanged; the supplied voltage/charge values are ENDPOINT values.
    /// This is a fixed-reference linear dynamic solve, not an equilibrium reset.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        dt_s: f64,
        endpoint_electrode_conditions: &[ElectrodeCondition],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricTransientCandidate<'model, 'mesh>> {
        self.advance_with_loads(
            dt_s,
            None,
            endpoint_electrode_conditions,
            &[],
            &[],
            &[],
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
        .map(|(candidate, _, _)| candidate)
    }

    /// Drive the original constrained coordinates by their endpoint velocities.
    /// Supply each existing support DOF exactly once, including stationary ones;
    /// entries need not be ordered. No free coordinate may be prescribed here.
    ///
    /// The step integrates `u_next = u_old + dt*(v_old + v_next)/2` at supports
    /// and retains the complete mass/stiffness/coupling columns in the same
    /// solve. Independent reaction work pays for the resulting motion, electrical
    /// response and intrinsic material heat. This avoids differentiating an
    /// arbitrary displacement trace or introducing a separate base-force model.
    /// Ordinary fixed-support methods refuse a history with moving supports;
    /// use this method again, with zero endpoint velocity for an explicit stop.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_support_motion_candidate(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: &[(usize, f64)],
        endpoint_electrode_conditions: &[ElectrodeCondition],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricTransientCandidate<'model, 'mesh>> {
        self.advance_with_loads(
            dt_s,
            Some(endpoint_support_velocities_m_s),
            endpoint_electrode_conditions,
            &[],
            &[],
            &[],
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
        .map(|(candidate, _, _)| candidate)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_loads(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        endpoint_electrode_conditions: &[ElectrodeCondition],
        loads: &[load::ResistiveElectrodeLoad],
        reactive_loads: &[load::reactive::ReactiveElectrodeLoad],
        reactive_states: &[load::reactive::ReactiveElectrodeState],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        Vec<load::ElectrodeLoadHeat>,
        load::reactive::ReactiveUpdate,
    )> {
        self.advance_with_acoustic_loads(
            dt_s,
            endpoint_support_velocities_m_s,
            endpoint_electrode_conditions,
            loads,
            reactive_loads,
            reactive_states,
            &[],
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_acoustic_loads(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        endpoint_electrode_conditions: &[ElectrodeCondition],
        loads: &[load::ResistiveElectrodeLoad],
        reactive_loads: &[load::reactive::ReactiveElectrodeLoad],
        reactive_states: &[load::reactive::ReactiveElectrodeState],
        acoustic_ports: &[acoustic::TransientAcousticPort<'_>],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        Vec<load::ElectrodeLoadHeat>,
        load::reactive::ReactiveUpdate,
    )> {
        poll(cx)?;
        validate_controls(controls)?;
        let n = self.model.nodes_m.len();
        let nu = 3 * n;
        let inverse_dt = 1.0 / dt_s;
        let mass_coefficient = (2.0 * inverse_dt).powi(2);
        let time_s = self.state.time_s + dt_s;
        if !positive(dt_s)
            || !positive(inverse_dt)
            || !positive(mass_coefficient)
            || !time_s.is_finite()
            || time_s <= self.state.time_s
            || midpoint_nodal_forces_n.len() != n
            || endpoint_nodal_free_charge_c.len() != n
            || midpoint_nodal_forces_n
                .iter()
                .flatten()
                .chain(endpoint_nodal_free_charge_c)
                .any(|value| !value.is_finite())
            || endpoint_electrode_conditions.len() != self.state.electrode_conditions.len()
        {
            return Err(invalid(
                "representable positive timestep and complete finite source values required",
            ));
        }
        let acoustic_ports = acoustic::prepare(self.model, acoustic_ports, cx)?;
        for (old, next) in self
            .state
            .electrode_conditions
            .iter()
            .zip(endpoint_electrode_conditions)
        {
            if !matches!(
                (old, next),
                (
                    ElectrodeCondition::Voltage(_),
                    ElectrodeCondition::Voltage(_)
                ) | (ElectrodeCondition::Charge(_), ElectrodeCondition::Charge(_))
            ) {
                return Err(invalid(
                    "changing electrode source type requires a separately accounted switching event",
                ));
            }
        }
        for (id, cell) in self.state.field.cells.iter().enumerate() {
            poll(cx)?;
            if norm(&cell.strain_mandel) > controls.maximum_strain_norm
                || self.model.has_thermal_bias()
                    && norm(&self.model.elastic_strain(id, cell.strain_mandel))
                        > controls.maximum_strain_norm
                || norm(&cell.response.electric_field_v_m) > controls.maximum_electric_field_v_m
            {
                return Err(invalid(
                    "initial state leaves the requested small-signal envelope",
                ));
            }
        }
        let (mut next_conditions, references) =
            normalize_conditions(&self.boundary, endpoint_electrode_conditions, cx)?;
        let (old_conditions, _) =
            normalize_conditions(&self.boundary, &self.state.electrode_conditions, cx)?;
        let mut midpoint_conditions = old_conditions
            .iter()
            .zip(&next_conditions)
            .map(|(old, next)| match (*old, *next) {
                (ElectrodeCondition::Voltage(a), ElectrodeCondition::Voltage(b)) => {
                    ElectrodeCondition::Voltage(0.5 * a + 0.5 * b)
                }
                (ElectrodeCondition::Charge(a), ElectrodeCondition::Charge(b)) => {
                    ElectrodeCondition::Charge(0.5 * a + 0.5 * b)
                }
                _ => unreachable!("source types checked before normalization"),
            })
            .collect::<Vec<_>>();
        let mut midpoint_charge = self
            .state
            .nodal_free_charge_c
            .iter()
            .zip(endpoint_nodal_free_charge_c)
            .map(|(old, next)| 0.5 * old + 0.5 * next)
            .collect::<Vec<_>>();
        let mut boundary = self.boundary.map.clone();
        boundary.prescribed.fill(0.0);
        boundary.rhs.fill(0.0);
        boundary.voltage_reference.fill(0.0);
        prescribe_support_motion(
            self,
            dt_s,
            endpoint_support_velocities_m_s,
            &mut boundary,
            cx,
        )?;
        let momentum = mass_action(self.model, &self.state.velocity_m_s, cx)?;
        for node in 0..n {
            poll(cx)?;
            for axis in 0..3 {
                if let Some(row) = boundary.index[3 * node + axis] {
                    boundary.rhs[row] = midpoint_nodal_forces_n[node][axis]
                        - self.state.field.support_force_n[node][axis]
                        + 2.0 * inverse_dt * momentum[node][axis];
                }
            }
            if let Some(row) = boundary.index[nu + node] {
                boundary.rhs[row] -= 0.5 * endpoint_nodal_free_charge_c[node]
                    - 0.5 * self.state.nodal_free_charge_c[node];
            }
        }
        for (id, nodes) in self.boundary.electrode_nodes.iter().enumerate() {
            poll(cx)?;
            match (old_conditions[id], next_conditions[id]) {
                (ElectrodeCondition::Voltage(old), ElectrodeCondition::Voltage(next)) => {
                    for &node in nodes {
                        boundary.prescribed[nu + node] = 0.5 * next - 0.5 * old;
                    }
                }
                (ElectrodeCondition::Charge(old), ElectrodeCondition::Charge(next)) => {
                    let row = boundary.index[nu + nodes[0]].expect("admitted charge electrode");
                    boundary.rhs[row] -= 0.5 * next - 0.5 * old;
                }
                _ => unreachable!("source types checked before normalization"),
            }
        }
        let mut load_entries = load::stamp(self, loads, dt_s, &mut boundary, cx)?;
        load_entries.extend(load::reactive::stamp(
            self,
            reactive_loads,
            reactive_states,
            dt_s,
            &mut boundary,
            controls,
            cx,
        )?);
        relaxation::stamp(self, dt_s, &mut boundary, &mut load_entries, cx)?;
        acoustic::stamp(&acoustic_ports, dt_s, &mut boundary, &mut load_entries, cx)?;
        let (matrix, rhs) =
            self.model
                .reduced_system(&boundary, mass_coefficient, &load_entries, cx)?;
        let (increment, iterations, relative_residual) =
            equilibrated_solve(matrix, rhs, controls, cx)?;
        let mut full_increment = boundary.prescribed.clone();
        for (dof, index) in boundary.index.iter().enumerate() {
            if let Some(index) = index {
                full_increment[dof] = increment[*index];
            }
        }
        let mut full_midpoint = vec![0.0; 4 * n];
        let mut full_endpoint = vec![0.0; 4 * n];
        let mut midpoint_displacement = vec![[0.0; 3]; n];
        let mut displacement = vec![[0.0; 3]; n];
        let mut velocity = vec![[0.0; 3]; n];
        let mut acceleration = vec![[0.0; 3]; n];
        let mut applied_force_work = Products::default();
        for node in 0..n {
            poll(cx)?;
            for axis in 0..3 {
                let dof = 3 * node + axis;
                let delta = full_increment[dof];
                let old_u = self.state.field.displacement_m[node][axis];
                let old_v = self.state.velocity_m_s[node][axis];
                midpoint_displacement[node][axis] = old_u + delta;
                displacement[node][axis] = old_u + 2.0 * delta;
                velocity[node][axis] = 4.0 * inverse_dt * delta - old_v;
                acceleration[node][axis] = 2.0 * inverse_dt * (2.0 * inverse_dt * delta - old_v);
                full_midpoint[dof] = midpoint_displacement[node][axis];
                full_endpoint[dof] = displacement[node][axis];
                applied_force_work.add(
                    Dd::from_f64(midpoint_nodal_forces_n[node][axis]),
                    Dd::from_f64(displacement[node][axis]) - Dd::from_f64(old_u),
                );
            }
            full_midpoint[nu + node] =
                self.state.relative_voltage_v[node] + full_increment[nu + node];
            full_endpoint[nu + node] =
                self.state.relative_voltage_v[node] + 2.0 * full_increment[nu + node];
        }
        if let Some(supports) = endpoint_support_velocities_m_s {
            for &(dof, endpoint_velocity) in supports {
                poll(cx)?;
                let (node, axis) = (dof / 3, dof % 3);
                // Retain the supplied endpoint velocity without subtracting
                // nearly equal values after its kinematic elimination.
                velocity[node][axis] = endpoint_velocity;
                acceleration[node][axis] =
                    inverse_dt * (endpoint_velocity - self.state.velocity_m_s[node][axis]);
            }
        }
        for (id, nodes) in self.boundary.electrode_nodes.iter().enumerate() {
            if let ElectrodeCondition::Voltage(value) = next_conditions[id] {
                for &node in nodes {
                    full_endpoint[nu + node] = value;
                }
            }
        }
        let mut accepted_conditions = endpoint_electrode_conditions.to_vec();
        let load_heat = load::recover_load_charge(
            self,
            loads,
            dt_s,
            &full_midpoint[nu..],
            &mut accepted_conditions,
            &mut next_conditions,
            &mut midpoint_conditions,
            cx,
        )?;
        let reactive = load::reactive::recover(
            self,
            reactive_loads,
            reactive_states,
            dt_s,
            &full_midpoint[nu..],
            &mut accepted_conditions,
            &mut next_conditions,
            &mut midpoint_conditions,
            controls,
            cx,
        )?;
        let material_update =
            relaxation::advance(self, dt_s, &full_increment, &full_midpoint, controls, cx)?;
        let inertial_force = mass_action(self.model, &acceleration, cx)?;
        let mut effective_force = midpoint_nodal_forces_n
            .iter()
            .zip(&inertial_force)
            .map(|(load, inertia)| core::array::from_fn(|axis| load[axis] - inertia[axis]))
            .collect::<Vec<[f64; 3]>>();
        let mut endpoint_force = vec![[0.0; 3]; n];
        let endpoint_charge = if let Some(update) = &material_update {
            for node in 0..n {
                poll(cx)?;
                for axis in 0..3 {
                    effective_force[node][axis] -= update.midpoint.force_n[node][axis];
                    endpoint_force[node][axis] = -update.endpoint.force_n[node][axis];
                }
                midpoint_charge[node] -= update.midpoint.charge_c[node];
            }
            Some(
                endpoint_nodal_free_charge_c
                    .iter()
                    .zip(&update.endpoint.charge_c)
                    .map(|(charge, polarization)| charge - polarization)
                    .collect::<Vec<_>>(),
            )
        } else {
            None
        };
        let acoustic_responses = acoustic::recover(
            &acoustic_ports,
            dt_s,
            &full_increment,
            &self.state.field.displacement_m,
            &displacement,
            &mut effective_force,
            controls.energy_tolerance_j,
            cx,
        )?;
        let acoustic_outward_work_j = acoustic_responses
            .iter()
            .map(|port| port.net_outward_work_j)
            .sum::<f64>();
        let midpoint_electrodes = electrode_views(&self.boundary, &midpoint_conditions);
        let mut midpoint = self.model.recover(
            &full_midpoint,
            midpoint_displacement,
            full_midpoint[nu..].to_vec(),
            &boundary,
            &midpoint_electrodes,
            &effective_force,
            &midpoint_charge,
            controls,
            iterations,
            relative_residual,
            cx,
        )?;
        // All displacement coordinates are known at the endpoint. Reconstruct
        // cell fields/forces directly; no second linear solve or static reset.
        let mut endpoint_boundary = boundary;
        endpoint_boundary.index[..nu].fill(None);
        let endpoint_electrodes = electrode_views(&self.boundary, &next_conditions);
        let mut field = self.model.recover(
            &full_endpoint,
            displacement,
            full_endpoint[nu..].to_vec(),
            &endpoint_boundary,
            &endpoint_electrodes,
            &endpoint_force,
            endpoint_charge
                .as_deref()
                .unwrap_or(endpoint_nodal_free_charge_c),
            controls,
            iterations,
            relative_residual,
            cx,
        )?;
        let (mechanical_material_heat_j, dielectric_material_heat_j, material_memory_energy_j) =
            if let Some(update) = &material_update {
                update.midpoint.enrich(self.model, &mut midpoint, cx)?;
                update.endpoint.enrich(self.model, &mut field, cx)?;
                (
                    update.mechanical_heat_j,
                    update.dielectric_heat_j,
                    update.memory_energy_j,
                )
            } else {
                (0.0, 0.0, [0.0; 2])
            };
        let material_heat_j = mechanical_material_heat_j + dielectric_material_heat_j;
        let mut midpoint_support_force_n = core::mem::take(&mut midpoint.support_force_n);
        let mut support_work = Products::default();
        for node in 0..n {
            poll(cx)?;
            for axis in 0..3 {
                if self.boundary.map.index[3 * node + axis].is_none() {
                    support_work.add(
                        Dd::from_f64(midpoint_support_force_n[node][axis]),
                        Dd::from_f64(field.displacement_m[node][axis])
                            - Dd::from_f64(self.state.field.displacement_m[node][axis]),
                    );
                } else {
                    midpoint_support_force_n[node][axis] = 0.0;
                }
            }
        }
        let applied_force_work = applied_force_work.value();
        let support_work = support_work.value();
        let mechanical_work = applied_force_work + support_work;
        let applied_force_work_j = applied_force_work.to_f64();
        let support_work_j = support_work.to_f64();
        let mechanical_work_j = mechanical_work.to_f64();
        let mut electrical_work = Products::default();
        for node in 0..n {
            poll(cx)?;
            electrical_work.add(
                Dd::from_f64(full_midpoint[nu + node]),
                Dd::from_f64(endpoint_nodal_free_charge_c[node])
                    - Dd::from_f64(self.state.nodal_free_charge_c[node]),
            );
        }
        for (id, nodes) in self.boundary.electrode_nodes.iter().enumerate() {
            let delta_charge = match (old_conditions[id], next_conditions[id]) {
                (ElectrodeCondition::Voltage(_), ElectrodeCondition::Voltage(_)) => {
                    Dd::from_f64(field.electrodes[id].charge_c)
                        - Dd::from_f64(self.state.field.electrodes[id].charge_c)
                }
                (ElectrodeCondition::Charge(old), ElectrodeCondition::Charge(next)) => {
                    Dd::from_f64(next) - Dd::from_f64(old)
                }
                _ => unreachable!("source types checked before normalization"),
            };
            electrical_work.add(Dd::from_f64(full_midpoint[nu + nodes[0]]), delta_charge);
        }
        let electrical_work = electrical_work.value();
        let electrical_work_j = electrical_work.to_f64();
        let kinetic_energy_j = kinetic_energy(self.model, &velocity, cx)?;
        let kinetic_change = kinetic_change(self.model, &self.state.velocity_m_s, &velocity, cx)?;
        let stored_change = field_storage_change(self.model, &self.state.field, &field, cx)?
            + material_update.as_ref().map_or(Dd::ZERO, |update| update.memory_energy_change);
        let energy_change = kinetic_change + stored_change;
        let kinetic_energy_change_j = kinetic_change.to_f64();
        let stored_energy_change_j = stored_change.to_f64();
        let energy_residual_j = (energy_change
            + Dd::from_f64(mechanical_material_heat_j)
            + Dd::from_f64(dielectric_material_heat_j)
            + Dd::from_f64(acoustic_outward_work_j)
            - mechanical_work - electrical_work).to_f64();
        if [
            mechanical_work_j,
            applied_force_work_j,
            support_work_j,
            electrical_work_j,
            acoustic_outward_work_j,
            kinetic_energy_change_j,
            stored_energy_change_j,
            material_heat_j,
            energy_residual_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("transient work or energy is not representable"));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "kinetic and material storage, intrinsic heat and acoustic outflow do not match supplied work",
            ));
        }
        let step = TetPiezoelectricTransientStep {
            time_step_s: dt_s,
            mechanical_work_j,
            applied_force_work_j,
            support_work_j,
            electrical_work_j,
            acoustic_outward_work_j,
            kinetic_energy_change_j,
            stored_energy_change_j,
            material_heat_j,
            mechanical_material_heat_j,
            dielectric_material_heat_j,
            material_memory_energy_j,
            energy_residual_j,
            maximum_force_residual_n: midpoint.maximum_force_residual_n,
            maximum_charge_residual_c: midpoint
                .maximum_charge_residual_c
                .max(field.maximum_charge_residual_c),
            iterations,
            relative_residual,
            kinetic_change,
            energy_change,
            mechanical_work,
        };
        let mut state = make_state(
            field,
            velocity,
            accepted_conditions,
            endpoint_nodal_free_charge_c.to_vec(),
            kinetic_energy_j,
            time_s,
            &references,
            &self.boundary,
            cx,
        )?;
        state.material_memory_energy_j = material_memory_energy_j[1];
        poll(cx)?;
        Ok((
            TetPiezoelectricTransientCandidate {
                body: Self {
                    model: self.model,
                    boundary: Arc::clone(&self.boundary),
                    state,
                    intrinsic: material_update.map(|update| update.next),
                },
                step,
                midpoint_support_force_n,
                acoustic_ports: acoustic_responses,
            },
            load_heat,
            reactive,
        ))
    }
}

fn prescribe_support_motion(
    body: &TetPiezoelectricTransient<'_, '_>,
    dt_s: f64,
    endpoint: Option<&[(usize, f64)]>,
    boundary: &mut BoundaryMap,
    cx: &Cx<'_>,
) -> Result<()> {
    let nu = 3 * body.model.nodes_m.len();
    let Some(endpoint) = endpoint else {
        for dof in 0..nu {
            poll(cx)?;
            if boundary.index[dof].is_none()
                && body.state.velocity_m_s[dof / 3][dof % 3] != 0.0
            {
                return Err(invalid(
                    "moving supports require explicit endpoint support velocities",
                ));
            }
        }
        return Ok(());
    };
    if endpoint.len() != boundary.index[..nu].iter().filter(|i| i.is_none()).count() {
        return Err(invalid("complete endpoint support velocities required"));
    }
    let mut seen = vec![false; nu];
    for &(dof, velocity) in endpoint {
        poll(cx)?;
        if dof >= nu || boundary.index[dof].is_some() || seen[dof] || !velocity.is_finite() {
            return Err(invalid(
                "finite unique endpoint velocities on the original support DOFs required",
            ));
        }
        seen[dof] = true;
        let midpoint_velocity =
            0.5 * body.state.velocity_m_s[dof / 3][dof % 3] + 0.5 * velocity;
        let increment = 0.5 * dt_s * midpoint_velocity;
        if !increment.is_finite() || midpoint_velocity != 0.0 && increment == 0.0 {
            return Err(invalid("prescribed support increment is not representable"));
        }
        boundary.prescribed[dof] = increment;
    }
    Ok(())
}

fn electrode_views<'a>(
    boundary: &'a FixedBoundary,
    conditions: &[ElectrodeCondition],
) -> Vec<PiezoelectricElectrode<'a>> {
    boundary
        .electrode_nodes
        .iter()
        .zip(conditions)
        .map(|(nodes, &condition)| PiezoelectricElectrode { nodes, condition })
        .collect()
}

fn normalize_conditions(
    boundary: &FixedBoundary,
    conditions: &[ElectrodeCondition],
    cx: &Cx<'_>,
) -> Result<(Vec<ElectrodeCondition>, Vec<f64>)> {
    let mut references = Vec::with_capacity(boundary.map.voltage_reference_electrode.len());
    for &id in &boundary.map.voltage_reference_electrode {
        poll(cx)?;
        match conditions[id] {
            ElectrodeCondition::Voltage(value) if value.is_finite() => references.push(value),
            _ => {
                return Err(invalid(
                    "finite prescribed potential required for each retained gauge",
                ));
            }
        }
    }
    let mut normalized = Vec::with_capacity(conditions.len());
    for (id, condition) in conditions.iter().enumerate() {
        poll(cx)?;
        let (condition, value) = match *condition {
            ElectrodeCondition::Voltage(value) => {
                let relative = value - references[boundary.electrode_nodes[id][0]];
                (ElectrodeCondition::Voltage(relative), relative)
            }
            ElectrodeCondition::Charge(value) => (ElectrodeCondition::Charge(value), value),
        };
        if !value.is_finite() {
            return Err(invalid(
                "finite electrode charge and representable voltage differences required",
            ));
        }
        normalized.push(condition);
    }
    Ok((normalized, references))
}

fn mass_action(
    model: &TetPiezoelectricModel<'_>,
    values: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>> {
    let mut applied = vec![[0.0; 3]; values.len()];
    mass_action_into(model, values.as_flattened(), applied.as_flattened_mut(), cx)?;
    Ok(applied)
}

fn mass_action_into(
    model: &TetPiezoelectricModel<'_>,
    values: &[f64],
    applied: &mut [f64],
    cx: &Cx<'_>,
) -> Result<()> {
    for (row, value) in applied.iter_mut().enumerate() {
        if row % 3 == 0 {
            poll(cx)?;
        }
        let (cols, entries) = model.solid.mass.row(row);
        *value = cols
            .iter()
            .zip(entries)
            .map(|(&col, &entry)| entry * values[col])
            .sum();
        if !value.is_finite() {
            return Err(invalid("consistent-mass response is not representable"));
        }
    }
    Ok(())
}

fn kinetic_change(
    model: &TetPiezoelectricModel<'_>,
    before: &[[f64; 3]],
    after: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<Dd> {
    quadratic::increment(
        before.as_flattened(), after.as_flattened(),
        |values, applied, _| mass_action_into(model, values, applied, cx),
        &mut || cx.checkpoint().is_err(),
    ).map(|increment| increment.change).map_err(|error| match error {
        QuadraticError::Operator(error) => error,
        QuadraticError::Cancelled => TetPiezoelectricError::Cancelled,
        QuadraticError::Dimension => invalid("consistent-mass endpoint dimensions disagree"),
        QuadraticError::NonFinite => invalid("consistent-mass increment is not representable"),
    })
}

// Use the same original local field coordinates/tensors as material recovery.
// Intrinsic memory storage is added from its own native increment separately.
fn field_storage_change(
    model: &TetPiezoelectricModel<'_>,
    before: &TetPiezoelectricSolution,
    after: &TetPiezoelectricSolution,
    cx: &Cx<'_>,
) -> Result<Dd> {
    let mut change = Products::default();
    for (cell, (old, new)) in before.cells.iter().zip(&after.cells).enumerate() {
        poll(cx)?;
        let material = &model.materials[if model.materials.len() == 1 { 0 } else { cell }];
        let volume = Dd::from_f64(model.dielectric.volumes_m3()[cell]);
        let elastic_change = if model.has_thermal_bias() {
            energy::quadratic_change_shifted(
                material.stiffness_mandel_pa(), &old.strain_mandel, &new.strain_mandel,
                &model.thermal_strain(cell),
            )
        } else {
            energy::quadratic_change(
                material.stiffness_mandel_pa(), &old.strain_mandel, &new.strain_mandel,
            )
        };
        change.add(volume, elastic_change);
        change.add(volume, energy::quadratic_change(
            material.permittivity_f_m(),
            &old.response.electric_field_v_m, &new.response.electric_field_v_m,
        ));
    }
    Ok(change.value())
}

fn kinetic_energy(
    model: &TetPiezoelectricModel<'_>,
    velocity: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<f64> {
    let momentum = mass_action(model, velocity, cx)?;
    let energy = velocity
        .iter()
        .flatten()
        .zip(momentum.iter().flatten())
        .map(|(v, p)| 0.5 * v * p)
        .sum::<f64>();
    if !energy.is_finite()
        || energy < 0.0
        || energy == 0.0 && velocity.iter().flatten().any(|&v| v != 0.0)
    {
        return Err(invalid(
            "consistent-mass kinetic energy is not representable",
        ));
    }
    Ok(energy)
}

#[allow(clippy::too_many_arguments)]
fn make_state(
    mut field: TetPiezoelectricSolution,
    velocity_m_s: Vec<[f64; 3]>,
    electrode_conditions: Vec<ElectrodeCondition>,
    nodal_free_charge_c: Vec<f64>,
    kinetic_energy_j: f64,
    time_s: f64,
    references: &[f64],
    boundary: &FixedBoundary,
    cx: &Cx<'_>,
) -> Result<TetPiezoelectricTransientState> {
    let relative_voltage_v = field.voltage_v.clone();
    for (value, reference) in field.voltage_v.iter_mut().zip(references) {
        poll(cx)?;
        *value += reference;
        if !value.is_finite() {
            return Err(invalid("absolute nodal potential is not representable"));
        }
    }
    for (response, nodes) in field.electrodes.iter_mut().zip(&boundary.electrode_nodes) {
        response.voltage_v = field.voltage_v[nodes[0]];
    }
    if !(kinetic_energy_j + field.stored_energy_j).is_finite() {
        return Err(invalid(
            "total electromechanical energy is not representable",
        ));
    }
    Ok(TetPiezoelectricTransientState {
        time_s,
        field,
        velocity_m_s,
        relative_voltage_v,
        electrode_conditions,
        nodal_free_charge_c,
        kinetic_energy_j,
        material_memory_energy_j: 0.0,
    })
}
