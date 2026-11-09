//! Finite-amplitude membrane loss with exact passive material-channel flows.
//!
//! The original nonlinear momentum Verlet owner supplies elasticity, inertia,
//! supports and held-load work. At each endpoint geometry the material Gram
//! rows supply C(u) = sum s s^T. Every rank-one viscous flow is exponential;
//! forward/reverse sweeps and damping/elastic/damping composition are symmetric
//! and second order. Noncommuting rows retain a splitting error: this is not
//! the full linear exponential provided by viscosity::modal.
//!
//! Heat is integrated analytically along each material-channel flow. Elastic
//! integration error and floating-point loss-step error remain signed numerical
//! defects, separately bounded by the caller. No nonlinear stability, bandwidth,
//! continuum convergence or audio-render qualification follows from acceptance.

use super::{InverseMass, MembraneError, MembraneMotion, MembranePressure, kinetic, poll};
use crate::membrane::viscosity::{MembraneViscosity, MembraneViscosityError};
use fs_math::det;
use fs_solver::op::LinearOp;

/// Explicit work, time-resolution and numerical-energy limits.
#[derive(Clone, Copy, Debug)]
pub struct MembraneDampingBudget {
    /// Six original source rows per triangle, including zero rows.
    pub max_source_rows: usize,
    /// Damping inverse actions per step. Lumped inertia needs none; consistent
    /// inertia admits the worst case 4*source_rows+4 before stepping. The
    /// existing elastic step retains its separate inverse-action costs.
    pub max_mass_solves_per_step: usize,
    /// Positive finite caller-resolved elastic/splitting time limit [s].
    pub max_step_s: f64,
    /// Absolute mechanical-energy/work/heat defect allowed [J]. This controls
    /// numerical energy error, not trajectory or phase accuracy.
    pub max_energy_defect_j: f64,
}

/// Physical observations of one completely accepted nonlinear damped step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DampedMembraneStep {
    pub time_s: f64,
    pub strain_energy_j: f64,
    pub kinetic_energy_j: f64,
    pub dead_load_work_j: f64,
    /// Current-area pressure work during the original elastic drift [J].
    pub pressure_work_j: f64,
    /// Analytically integrated nonnegative material heat over the split flow.
    pub dissipated_j: f64,
    /// Mechanical-energy change + heat - held-load work [J].
    pub energy_defect_j: f64,
    /// Elastic/viscous support impulse, including consistent M_cf delta_v.
    pub support_impulse_n_s: [f64; 3],
    /// Actual damping inverse actions; excludes the original elastic step.
    pub mass_solves: usize,
}

#[derive(Debug)]
pub enum MembraneDampingError {
    Membrane(MembraneError),
    Viscosity(MembraneViscosityError),
    EnergyBudget { defect_j: f64, maximum_j: f64 },
}
impl From<MembraneError> for MembraneDampingError {
    fn from(error: MembraneError) -> Self {
        Self::Membrane(error)
    }
}
impl From<MembraneViscosityError> for MembraneDampingError {
    fn from(error: MembraneViscosityError) -> Self {
        Self::Viscosity(error)
    }
}
impl core::fmt::Display for MembraneDampingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Membrane(error) => write!(f, "damped membrane: {error}"),
            Self::Viscosity(error) => write!(f, "damped membrane: {error}"),
            Self::EnergyBudget {
                defect_j,
                maximum_j,
            } => write!(
                f,
                "membrane energy defect {defect_j} J exceeds {maximum_j} J"
            ),
        }
    }
}
impl std::error::Error for MembraneDampingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Membrane(error) => Some(error),
            Self::Viscosity(error) => Some(error),
            Self::EnergyBudget { .. } => None,
        }
    }
}
type Result<T> = core::result::Result<T, MembraneDampingError>;
type ElementRows = [[[f64; 3]; 3]; 6];
fn invalid(what: &'static str) -> MembraneDampingError {
    MembraneError::Invalid(what).into()
}

/// Resumable nonlinear damping on an existing source-bound motion owner.
///
/// Reuses the selected lumped, dense consistent or sparse consistent reference
/// inertia. Lumped channel work touches at most nine local coordinates, making
/// each sweep O(elements). Consistent channels use the original prepared inverse
/// and its nonlocal response. Storage is O(nodes+elements) plus inertia factors;
/// no global damping matrix, eigenbasis or mode-count ceiling is introduced.
///
/// An entire trial motion and loss buffers are prepared once. Warm steps
/// allocate no buffers in this owner; source constitutive evaluations retain
/// their own costs. Cloning is a complete checkpoint. Failure changes scratch
/// only; state, time and accepted heat publish after the final cancellation poll.
#[derive(Clone, Debug)]
pub struct DampedMembraneMotion<'model, 'patch> {
    viscosity: &'model MembraneViscosity<'patch>,
    motion: MembraneMotion<'patch>,
    trial: MembraneMotion<'patch>,
    budget: MembraneDampingBudget,
    free_index: Vec<usize>,
    rows: Vec<ElementRows>,
    before_momentum: Vec<f64>,
    before_velocity: Vec<f64>,
    direction: Vec<f64>,
    inverse_direction: Vec<f64>,
    element_dissipated_j: Vec<f64>,
    trial_element_dissipated_j: Vec<f64>,
    cumulative_dissipated_j: f64,
}

impl<'model, 'patch> DampedMembraneMotion<'model, 'patch> {
    /// Consume the accepted original motion without changing its physical state.
    /// The viscosity must borrow that identical patch, including its natural
    /// reference. Its declared source band remains an applicability limit;
    /// this constructor cannot infer a nonlinear trajectory's bandwidth.
    pub fn new(
        viscosity: &'model MembraneViscosity<'patch>,
        motion: MembraneMotion<'patch>,
        budget: MembraneDampingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(viscosity.patch(), motion.patch) {
            return Err(invalid(
                "motion and viscosity require the identical membrane owner",
            ));
        }
        if !budget.max_step_s.is_finite()
            || budget.max_step_s <= 0.0
            || !budget.max_energy_defect_j.is_finite()
            || budget.max_energy_defect_j < 0.0
        {
            return Err(invalid(
                "finite positive time and nonnegative energy-defect budgets required",
            ));
        }
        let elements = motion.patch.elements.len();
        let source_rows = elements
            .checked_mul(6)
            .ok_or(MembraneError::Budget("membrane damping row extent"))?;
        if source_rows > budget.max_source_rows {
            return Err(MembraneError::Budget("membrane damping source rows").into());
        }
        if matches!(motion.inverse_mass, InverseMass::Consistent { .. }) {
            let solves = source_rows
                .checked_mul(4)
                .and_then(|n| n.checked_add(4))
                .ok_or(MembraneError::Budget(
                    "membrane damping inverse-action extent",
                ))?;
            if solves > budget.max_mass_solves_per_step {
                return Err(
                    MembraneError::Budget("membrane damping consistent inverse actions").into(),
                );
            }
        }
        let mut free_index = vec![usize::MAX; motion.fixed.len()];
        for (i, &dof) in motion.free.iter().enumerate() {
            free_index[dof] = i;
        }
        let mut rows = Vec::with_capacity(elements);
        for index in 0..elements {
            poll(&mut cancelled)?;
            rows.push(viscosity.element_dissipation_rows(index, &motion.displacement_m)?);
        }
        let free = motion.free.len();
        let trial = motion.clone();
        poll(&mut cancelled)?;
        Ok(Self {
            viscosity,
            motion,
            trial,
            budget,
            free_index,
            rows,
            before_momentum: vec![0.0; free],
            before_velocity: vec![0.0; free],
            direction: vec![0.0; free],
            inverse_direction: vec![0.0; free],
            element_dissipated_j: vec![0.0; elements],
            trial_element_dissipated_j: vec![0.0; elements],
            cumulative_dissipated_j: 0.0,
        })
    }
    /// Accepted original motion, including its selected reference inertia.
    pub fn motion(&self) -> &MembraneMotion<'patch> {
        &self.motion
    }
    pub fn viscosity(&self) -> &MembraneViscosity<'patch> {
        self.viscosity
    }
    pub const fn budget(&self) -> MembraneDampingBudget {
        self.budget
    }
    /// Heat by original triangle in the last accepted step; initially all zero.
    pub fn element_dissipated_j(&self) -> &[f64] {
        &self.element_dissipated_j
    }
    /// Material heat since this wrapper was constructed, in joules.
    pub const fn cumulative_dissipated_j(&self) -> f64 {
        self.cumulative_dissipated_j
    }

    /// Symmetric damping/elastic/damping step at a resolved positive duration.
    /// Dissipative time reversal is not admitted. Stationary supports and held
    /// nodal loads retain the original motion semantics.
    pub fn step(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DampedMembraneStep> {
        self.step_with_pressure(dt_s, dead_load_n, 0.0, cancelled)
    }

    /// Compose the same material damping with held follower pressure in the
    /// original elastic map. Pressure does work only during its geometric
    /// drift, not during the two fixed-geometry damping stages. That work is
    /// reported separately and never counted as material heat. Current-area
    /// traction, material P1 samples, path work and support impulses retain
    /// [`MembraneMotion::step_with_pressure`]'s semantics and limitations.
    pub fn step_with_pressure<'p>(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DampedMembraneStep> {
        poll(&mut cancelled)?;
        let pressure = pressure.into();
        pressure.validate(self.motion.patch.elements.len(), &mut cancelled)?;
        let next_time = self.motion.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || dt_s > self.budget.max_step_s
            || 0.25 * dt_s == 0.0
            || !next_time.is_finite()
            || next_time == self.motion.time_s
            || !crate::membrane::finite_vectors(dead_load_n, self.motion.displacement_m.len())
        {
            return Err(invalid(
                "resolved positive membrane time step, finite pressure and complete loads required",
            ));
        }
        copy_state(&self.motion, &mut self.trial);
        self.trial_element_dissipated_j.fill(0.0);
        let mut mass_solves = 0;
        let first = self.damping_stage(0.5 * dt_s, &mut mass_solves, &mut cancelled)?;
        let elastic = self
            .trial
            .step_with_pressure(dt_s, dead_load_n, pressure, &mut cancelled)?;
        let last = self.damping_stage(0.5 * dt_s, &mut mass_solves, &mut cancelled)?;
        let dissipated_j = first.heat + last.heat;
        let energy_defect_j = first.defect + elastic.energy_defect_j + last.defect;
        let cumulative = self.cumulative_dissipated_j + dissipated_j;
        let support_impulse_n_s: [f64; 3] = core::array::from_fn(|axis| {
            first.support[axis] + elastic.support_impulse_n_s[axis] + last.support[axis]
        });
        if [
            dissipated_j,
            energy_defect_j,
            cumulative,
            self.trial.kinetic_energy_j,
            self.trial.strain_energy_j + self.trial.kinetic_energy_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || dissipated_j < 0.0
            || self.trial.kinetic_energy_j < 0.0
            || support_impulse_n_s.iter().any(|x| !x.is_finite())
            || self
                .trial_element_dissipated_j
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(invalid(
                "unrepresentable nonlinear damping state or work/heat balance",
            ));
        }
        if mass_solves > self.budget.max_mass_solves_per_step {
            return Err(MembraneError::Budget("membrane damping inverse actions exhausted").into());
        }
        if energy_defect_j.abs() > self.budget.max_energy_defect_j {
            return Err(MembraneDampingError::EnergyBudget {
                defect_j: energy_defect_j,
                maximum_j: self.budget.max_energy_defect_j,
            });
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.motion, &mut self.trial);
        core::mem::swap(
            &mut self.element_dissipated_j,
            &mut self.trial_element_dissipated_j,
        );
        self.cumulative_dissipated_j = cumulative;
        Ok(DampedMembraneStep {
            time_s: self.motion.time_s,
            strain_energy_j: self.motion.strain_energy_j,
            kinetic_energy_j: self.motion.kinetic_energy_j,
            dead_load_work_j: elastic.dead_load_work_j,
            pressure_work_j: elastic.pressure_work_j,
            dissipated_j,
            energy_defect_j,
            support_impulse_n_s,
            mass_solves,
        })
    }

    fn damping_stage(
        &mut self,
        duration: f64,
        mass_solves: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LossStage> {
        poll(cancelled)?;
        self.before_momentum.copy_from_slice(&self.trial.momentum);
        for (i, &dof) in self.trial.free.iter().enumerate() {
            self.before_velocity[i] = self.trial.velocity_m_s[dof / 3][dof % 3];
        }
        for index in 0..self.rows.len() {
            poll(cancelled)?;
            self.rows[index] = self
                .viscosity
                .element_dissipation_rows(index, &self.trial.displacement_m)?;
        }
        let count = self.rows.len() * 6;
        let mut result = LossStage {
            heat: 0.0,
            defect: 0.0,
            support: [0.0; 3],
        };
        // Every row advances duration/2 on both sweeps; one forward sweep
        // alone would lose second-order accuracy for noncommuting channels.
        for reverse in [false, true] {
            for offset in 0..count {
                poll(cancelled)?;
                let row = if reverse { count - 1 - offset } else { offset };
                let element = row / 6;
                let nodes = self.trial.patch.elements[element].nodes;
                let (heat, support) = decay_row(
                    &mut self.trial,
                    &self.free_index,
                    nodes,
                    self.rows[element][row % 6],
                    0.5 * duration,
                    &mut self.direction,
                    &mut self.inverse_direction,
                    mass_solves,
                )?;
                result.heat += heat;
                self.trial_element_dissipated_j[element] += heat;
                for axis in 0..3 {
                    result.support[axis] += support[axis];
                }
            }
        }
        if matches!(self.trial.inverse_mass, InverseMass::Consistent { .. }) {
            poll(cancelled)?;
            self.trial
                .inverse_mass
                .apply(&self.trial.momentum, &mut self.trial.trial_velocity);
            *mass_solves += 1;
            for (&dof, &v) in self.trial.free.iter().zip(&self.trial.trial_velocity) {
                self.trial.velocity_m_s[dof / 3][dof % 3] = v;
            }
        }
        if self.trial.inverse_mass.has_support_coupling() {
            poll(cancelled)?;
            for ((delta, new), old) in self
                .direction
                .iter_mut()
                .zip(&self.trial.momentum)
                .zip(&self.before_momentum)
            {
                *delta = new - old;
            }
            self.trial
                .inverse_mass
                .apply(&self.direction, &mut self.inverse_direction);
            *mass_solves += 1;
            let inertia = self
                .trial
                .inverse_mass
                .support_impulse(&self.inverse_direction);
            for axis in 0..3 {
                result.support[axis] += inertia[axis];
            }
        }
        let delta_kinetic: f64 = self
            .trial
            .momentum
            .iter()
            .zip(&self.before_momentum)
            .zip(&self.before_velocity)
            .zip(&self.trial.free)
            .map(|(((new, old), v0), &dof)| {
                (new - old) * (0.5 * v0 + 0.5 * self.trial.velocity_m_s[dof / 3][dof % 3])
            })
            .sum();
        self.trial.kinetic_energy_j = kinetic(
            &self.trial.momentum,
            &self.trial.inverse_mass,
            &self.trial.trial_velocity,
        );
        result.defect = delta_kinetic + result.heat;
        if !result.heat.is_finite()
            || !result.defect.is_finite()
            || !crate::membrane::finite_vectors(
                &self.trial.velocity_m_s,
                self.trial.patch.mesh.nodes.len(),
            )
        {
            return Err(invalid("unrepresentable exact membrane loss stage"));
        }
        poll(cancelled)?;
        Ok(result)
    }
}

struct LossStage {
    heat: f64,
    defect: f64,
    support: [f64; 3],
}

// Same immutable material, constraints and inertia; the elastic step
// initializes its own scratch. No allocation or partial public state changes.
fn copy_state(source: &MembraneMotion<'_>, destination: &mut MembraneMotion<'_>) {
    destination.q.copy_from_slice(&source.q);
    destination.momentum.copy_from_slice(&source.momentum);
    destination
        .displacement_m
        .copy_from_slice(&source.displacement_m);
    destination
        .velocity_m_s
        .copy_from_slice(&source.velocity_m_s);
    destination.gradient.copy_from_slice(&source.gradient);
    destination.time_s = source.time_s;
    destination.strain_energy_j = source.strain_energy_j;
    destination.kinetic_energy_j = source.kinetic_energy_j;
}

#[allow(clippy::too_many_arguments)]
fn decay_row(
    motion: &mut MembraneMotion<'_>,
    free_index: &[usize],
    nodes: [usize; 3],
    source_row: [[f64; 3]; 3],
    duration: f64,
    direction: &mut [f64],
    inverse_direction: &mut [f64],
    mass_solves: &mut usize,
) -> Result<(f64, [f64; 3])> {
    let scale = source_row
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok((0.0, [0.0; 3]));
    }
    let mut row = source_row.map(|v| v.map(|x| x / scale));
    for axis in 0..3 {
        row[0][axis] = -(row[1][axis] + row[2][axis]);
    }
    if (1..3).any(|a| (0..3).any(|axis| source_row[a][axis] != 0.0 && row[a][axis] == 0.0)) {
        return Err(invalid(
            "unrepresentable normalized membrane material channel",
        ));
    }
    // Difference form preserves exact translation nulls even at large common
    // velocities. The source row at node zero is minus the other two rows.
    let mut action = 0.0;
    for a in 1..3 {
        for axis in 0..3 {
            action += row[a][axis]
                * (motion.velocity_m_s[nodes[a]][axis] - motion.velocity_m_s[nodes[0]][axis]);
        }
    }
    if !action.is_finite() {
        return Err(invalid("unrepresentable membrane material-channel rate"));
    }
    if action == 0.0 {
        return Ok((0.0, [0.0; 3]));
    }
    let mut inverse_norm = 0.0;
    match &motion.inverse_mass {
        InverseMass::Lumped(inverse) => {
            for a in 0..3 {
                for axis in 0..3 {
                    let i = free_index[3 * nodes[a] + axis];
                    if i != usize::MAX {
                        inverse_norm += row[a][axis] * (inverse[i] * row[a][axis]);
                    }
                }
            }
        }
        InverseMass::Consistent { .. } => {
            direction.fill(0.0);
            for a in 0..3 {
                for axis in 0..3 {
                    let i = free_index[3 * nodes[a] + axis];
                    if i != usize::MAX {
                        direction[i] = row[a][axis];
                    }
                }
            }
            motion.inverse_mass.apply(direction, inverse_direction);
            *mass_solves += 1;
            for a in 0..3 {
                for axis in 0..3 {
                    let i = free_index[3 * nodes[a] + axis];
                    if i != usize::MAX {
                        inverse_norm += row[a][axis] * inverse_direction[i];
                    }
                }
            }
        }
    }
    if !inverse_norm.is_finite() || inverse_norm <= 0.0 {
        return Err(invalid(
            "material channel needs a resolved positive reference-mass norm",
        ));
    }
    let root_norm = det::sqrt(inverse_norm);
    let clock_root = (scale * root_norm) * det::sqrt(duration);
    let clock = clock_root * clock_root;
    if clock.is_nan() || clock == 0.0 {
        return Err(invalid("unrepresentable membrane material-channel clock"));
    }
    // Positive infinite clocks use their exact saturation limit; expm1
    // retains physically small transfers without endpoint subtraction.
    let lost = -det::expm1(-clock);
    let twice_lost = -det::expm1(-2.0 * clock);
    let amplitude = action / root_norm;
    let impulse = (amplitude * lost) / root_norm;
    let heat = amplitude * (amplitude * (0.5 * twice_lost));
    if !impulse.is_finite()
        || !heat.is_finite()
        || heat < 0.0
        || (heat == 0.0 && amplitude != 0.0 && twice_lost != 0.0)
    {
        return Err(invalid(
            "unrepresentable membrane material-channel heat or impulse",
        ));
    }
    let mut support = [0.0; 3];
    for a in 0..3 {
        for axis in 0..3 {
            let i = free_index[3 * nodes[a] + axis];
            let transfer = impulse * row[a][axis];
            if i == usize::MAX {
                support[axis] += transfer;
            } else {
                motion.momentum[i] -= transfer;
                if let InverseMass::Lumped(inverse) = &motion.inverse_mass {
                    motion.velocity_m_s[nodes[a]][axis] = motion.momentum[i] * inverse[i];
                }
            }
        }
    }
    if matches!(motion.inverse_mass, InverseMass::Consistent { .. }) {
        for (&dof, response) in motion.free.iter().zip(inverse_direction) {
            motion.velocity_m_s[dof / 3][dof % 3] -= impulse * *response;
        }
    }
    Ok((heat, support))
}
