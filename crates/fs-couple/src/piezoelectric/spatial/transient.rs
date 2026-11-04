//! Consistent-mass, reversible spatial piezoelectric motion.
//!
//! The native solid mass participates only in mechanical momentum balance;
//! electrostatics remains a charge constraint. Implicit midpoint solves one
//! coupled increment with `K + 4 M / dt^2`, the original reciprocal `G`, and
//! `-K_e`. This resolves the retained mesh's elastic waves without replacing
//! the material by a fitted actuator gain or a single oscillator.
//!
//! Supports and electrode topology/types are fixed. Prescribed voltages and
//! charges vary linearly between accepted endpoints; force is the caller's
//! midpoint value. The accepted energy account is kinetic + elastic + electric
//! storage against `F_mid . delta_u + phi_mid . delta_Q`. There is no numerical
//! damping, inferred physical loss, radiation, thermal update, or automatic
//! time/mesh-error certificate. Small-strain geometry and material coefficients
//! stay fixed. A candidate never mutates the accepted body.

use super::{
    BoundaryMap, ElectrodeCondition, PiezoelectricElectrode, PiezoelectricElectrodeResponse,
    Result, TetPiezoelectricCell, TetPiezoelectricControls, TetPiezoelectricError,
    TetPiezoelectricModel, TetPiezoelectricSolution, equilibrated_solve, invalid, norm, poll,
    positive, prepare_boundary, validate_controls,
};
use fs_exec::Cx;
use std::sync::Arc;

/// Source-bound resistive loads and prescribed electrode currents on this solve.
pub mod load;

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
    /// Original mesh-node velocity [m/s]. Supported coordinates are zero.
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
    /// Positive internal elastic plus electric storage [J].
    pub fn stored_energy_j(&self) -> f64 {
        self.field.stored_energy_j
    }
    /// Kinetic plus elastic plus electric energy [J].
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
}

/// Independent physical account for one midpoint window.
#[derive(Clone, Copy, Debug)]
pub struct TetPiezoelectricTransientStep {
    /// Accepted time interval [s].
    pub time_step_s: f64,
    /// Applied midpoint force dotted with actual displacement change [J].
    pub mechanical_work_j: f64,
    /// Midpoint potential dotted with actual free-charge change [J].
    pub electrical_work_j: f64,
    /// Consistent-mass kinetic energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Elastic plus electric storage change [J].
    pub stored_energy_change_j: f64,
    /// Total energy change minus the two external supplies [J].
    pub energy_residual_j: f64,
    /// Independently reconstructed midpoint momentum residual [N].
    pub maximum_force_residual_n: f64,
    /// Largest midpoint/endpoint free-node or electrode-total residual [C].
    pub maximum_charge_residual_c: f64,
    /// Native MINRES iteration count for this window.
    pub iterations: usize,
    /// Independently recomputed equilibrated increment-system residual.
    pub relative_residual: f64,
}

/// Joint mechanical/electrical candidate; accepted state remains untouched.
pub struct TetPiezoelectricTransientCandidate<'model, 'mesh> {
    body: TetPiezoelectricTransient<'model, 'mesh>,
    step: TetPiezoelectricTransientStep,
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
        })
    }

    /// Accepted state; candidate creation never changes it.
    pub const fn state(&self) -> &TetPiezoelectricTransientState {
        &self.state
    }

    /// Advance one reversible midpoint window. Source types and topology must
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
            endpoint_electrode_conditions,
            &[],
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
        .map(|(candidate, _)| candidate)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_loads(
        &self,
        dt_s: f64,
        endpoint_electrode_conditions: &[ElectrodeCondition],
        loads: &[load::ResistiveElectrodeLoad],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(
        TetPiezoelectricTransientCandidate<'model, 'mesh>,
        Vec<load::ElectrodeLoadHeat>,
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
        let midpoint_charge = self
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
        let load_entries = load::stamp(self, loads, dt_s, &mut boundary, cx)?;
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
        let mut mechanical_work_j = 0.0;
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
                mechanical_work_j +=
                    midpoint_nodal_forces_n[node][axis] * (displacement[node][axis] - old_u);
            }
            full_midpoint[nu + node] =
                self.state.relative_voltage_v[node] + full_increment[nu + node];
            full_endpoint[nu + node] =
                self.state.relative_voltage_v[node] + 2.0 * full_increment[nu + node];
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
        let inertial_force = mass_action(self.model, &acceleration, cx)?;
        let effective_force = midpoint_nodal_forces_n
            .iter()
            .zip(&inertial_force)
            .map(|(load, inertia)| core::array::from_fn(|axis| load[axis] - inertia[axis]))
            .collect::<Vec<[f64; 3]>>();
        let midpoint_electrodes = electrode_views(&self.boundary, &midpoint_conditions);
        let midpoint = self.model.recover(
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
        let field = self.model.recover(
            &full_endpoint,
            displacement,
            full_endpoint[nu..].to_vec(),
            &endpoint_boundary,
            &endpoint_electrodes,
            &vec![[0.0; 3]; n],
            endpoint_nodal_free_charge_c,
            controls,
            iterations,
            relative_residual,
            cx,
        )?;
        let mut electrical_work_j = 0.0;
        for node in 0..n {
            poll(cx)?;
            electrical_work_j += full_midpoint[nu + node]
                * (endpoint_nodal_free_charge_c[node] - self.state.nodal_free_charge_c[node]);
        }
        for (id, nodes) in self.boundary.electrode_nodes.iter().enumerate() {
            let delta_charge = match (old_conditions[id], next_conditions[id]) {
                (ElectrodeCondition::Voltage(_), ElectrodeCondition::Voltage(_)) => {
                    field.electrodes[id].charge_c - self.state.field.electrodes[id].charge_c
                }
                (ElectrodeCondition::Charge(old), ElectrodeCondition::Charge(next)) => next - old,
                _ => unreachable!("source types checked before normalization"),
            };
            electrical_work_j += full_midpoint[nu + nodes[0]] * delta_charge;
        }
        let kinetic_energy_j = kinetic_energy(self.model, &velocity, cx)?;
        let kinetic_energy_change_j = kinetic_energy_j - self.state.kinetic_energy_j;
        let stored_energy_change_j = field.stored_energy_j - self.state.field.stored_energy_j;
        let energy_residual_j = kinetic_energy_change_j + stored_energy_change_j
            - mechanical_work_j
            - electrical_work_j;
        if [
            mechanical_work_j,
            electrical_work_j,
            kinetic_energy_change_j,
            stored_energy_change_j,
            energy_residual_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("transient work or energy is not representable"));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "kinetic, elastic and electric energy change does not match supplied work",
            ));
        }
        let step = TetPiezoelectricTransientStep {
            time_step_s: dt_s,
            mechanical_work_j,
            electrical_work_j,
            kinetic_energy_change_j,
            stored_energy_change_j,
            energy_residual_j,
            maximum_force_residual_n: midpoint.maximum_force_residual_n,
            maximum_charge_residual_c: midpoint
                .maximum_charge_residual_c
                .max(field.maximum_charge_residual_c),
            iterations,
            relative_residual,
        };
        let state = make_state(
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
        poll(cx)?;
        Ok((
            TetPiezoelectricTransientCandidate {
                body: Self {
                    model: self.model,
                    boundary: Arc::clone(&self.boundary),
                    state,
                },
                step,
            },
            load_heat,
        ))
    }
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
    for (node, result) in applied.iter_mut().enumerate() {
        poll(cx)?;
        for (axis, value) in result.iter_mut().enumerate() {
            let (cols, entries) = model.solid.mass.row(3 * node + axis);
            *value = cols
                .iter()
                .zip(entries)
                .map(|(&col, &entry)| entry * values[col / 3][col % 3])
                .sum();
            if !value.is_finite() {
                return Err(invalid("consistent-mass response is not representable"));
            }
        }
    }
    Ok(applied)
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
    })
}
