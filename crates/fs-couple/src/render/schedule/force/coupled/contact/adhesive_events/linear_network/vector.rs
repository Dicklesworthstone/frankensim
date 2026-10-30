//! Set-valued two-tangent friction on the complete original linear network.
//!
//! Each contact has ONE Euclidean disk, not two independently clipped scalar
//! forces. Static and kinetic capacities coincide in this monotone rung. For
//! one backward-Euler cell, `A = I + h C + h² K`, `v1 = A⁻¹(b - h B t)`.
//! Condensation gives the strictly convex disk problem
//! `min 1/2 tᵀ H t - wᵀ t`, `|t_i| <= capacity_i`, with
//! `H = h Bᵀ A⁻¹ B`, `w = Bᵀ A⁻¹ b`. Its optimality condition is the
//! unsmoothed Coulomb graph at the endpoint velocity. This permits genuine
//! sticking and changing slip directions without an arbitrary sign choice.
//!
//! Original modal diagonals and every spring/dashpot link are preserved. The
//! native matrix factor and tangent mobility are prepared once per cell size;
//! only at-most-16 contact unknowns enter the semismooth Newton iteration.
//! All workspaces are bounded and allocated at admission. Force, disk, state,
//! cancellation and energy checks precede a successful return. The caller
//! owns accepted state; a refused trial cannot influence a deterministic retry.
//!
//! Backward Euler prevents the midpoint stop/reversal artifact but is FIRST
//! ORDER and dissipative. Its exact discrete loss
//! `1/2 |v1-v0|² + 1/2 (q1-q0)ᵀ K (q1-q0)` is a separate NUMERICAL account,
//! never friction heat or material damping. Explicit cell size and cumulative
//! numerical-loss limits control this tradeoff; refinement remains necessary.
//! The force-law gap and floating-point work estimate are not certificates or
//! a trajectory-error bound. No event times, unequal static/kinetic strengths,
//! moving frames, evolving normal loads, memory or nonlinear mechanics are
//! inferred. Endpoint link-force maxima refer to this discrete trajectory.

use super::LinearNetworkLink;
use fs_exec::CancelGate;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_math::det;

#[path = "vector/solve.rs"]
mod solve;

pub const MAX_VECTOR_NETWORK_CONTACTS: usize = 8;
const ROWS: usize = 2 * MAX_VECTOR_NETWORK_CONTACTS;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VectorFrictionRegime {
    Inactive,
    Sticking,
    Sliding,
}

/// Accuracy, physical domain and finite work limits are independent controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VectorLinearNetworkConfig {
    pub maximum_modes: usize,
    pub maximum_connections: usize,
    pub maximum_contacts: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_cells: usize,
    /// Shared across every contact solve in the whole advance call.
    pub maximum_nonlinear_iterations: usize,
    pub maximum_backtracks: usize,
    pub maximum_step_s: f64,
    pub maximum_cell_s: f64,
    /// Lower Cholesky pivot of the unit-column Gram matrix at admission.
    pub minimum_column_pivot: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub maximum_abs_modal_coordinate: f64,
    pub maximum_abs_modal_velocity: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
    pub maximum_numerical_work_uncertainty_j: f64,
    pub maximum_numerical_dissipation_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum VectorLinearNetworkError {
    Invalid(&'static str),
    Nonfinite,
    Cancelled,
    Budget(&'static str),
    LinearSolve(LuWorkspaceError),
    NotConverged,
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl From<LuWorkspaceError> for VectorLinearNetworkError {
    fn from(value: LuWorkspaceError) -> Self {
        Self::LinearSolve(value)
    }
}
impl core::fmt::Display for VectorLinearNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "vector adhesive linear network refused: {self:?}")
    }
}
impl std::error::Error for VectorLinearNetworkError {}
type Error = VectorLinearNetworkError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VectorLinearNetworkStep {
    pub contact_count: usize,
    pub regimes: [VectorFrictionRegime; MAX_VECTOR_NETWORK_CONTACTS],
    pub initial_energy_j: f64,
    pub modal_energy_j: f64,
    pub connection_energy_j: f64,
    pub external_work_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub numerical_dissipation_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub contact_friction_dissipation_j: [f64; MAX_VECTOR_NETWORK_CONTACTS],
    pub accumulated_slip_m: [f64; MAX_VECTOR_NETWORK_CONTACTS],
    pub maximum_contact_speed_m_s: [f64; MAX_VECTOR_NETWORK_CONTACTS],
    pub maximum_contact_force_n: [f64; MAX_VECTOR_NETWORK_CONTACTS],
    pub contact_impulse_n_s: [[f64; 2]; MAX_VECTOR_NETWORK_CONTACTS],
    pub contact_time_moment_n_s2: [[f64; 2]; MAX_VECTOR_NETWORK_CONTACTS],
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub cells: usize,
    pub nonlinear_iterations: usize,
}

/// Prepared complete mechanical response and bounded simultaneous contact solve.
pub struct VectorLinearNetworkFlow {
    stiffness: Vec<f64>,
    damping: Vec<f64>,
    links: Vec<LinearNetworkLink>,
    columns: Vec<Vec<f64>>,
    config: VectorLinearNetworkConfig,
    matrix: Vec<f64>,
    mechanical: LuWorkspace,
    contact_solve: LuWorkspace,
    response_columns: Vec<f64>,
    cached_h: Option<f64>,
    mobility: [f64; ROWS * ROWS],
    lipschitz: f64,
    q: Vec<f64>,
    v: Vec<f64>,
    next_v: Vec<f64>,
    free_v: Vec<f64>,
    rhs: Vec<f64>,
    link_losses: Vec<f64>,
    link_numerical_losses: Vec<f64>,
    link_forces: Vec<f64>,
    link_force_bounds: Vec<f64>,
    force: [f64; ROWS],
    trial_force: [f64; ROWS],
    free_slip: [f64; ROWS],
    residual: [f64; ROWS],
    trial_residual: [f64; ROWS],
    jacobian: [f64; ROWS * ROWS],
    direction: [f64; ROWS],
}

impl VectorLinearNetworkFlow {
    pub fn new(
        modal_stiffness: Vec<f64>,
        modal_damping: Vec<f64>,
        links: Vec<LinearNetworkLink>,
        columns: Vec<[Vec<f64>; 2]>,
        config: VectorLinearNetworkConfig,
    ) -> Result<Self, Error> {
        validate_config(config)?;
        let n = modal_stiffness.len();
        let m = 2 * columns.len();
        if n == 0
            || n > config.maximum_modes
            || links.len() > config.maximum_connections
            || columns.is_empty()
            || columns.len() > config.maximum_contacts
        {
            return Err(Error::Budget("complete network and vector contact extent"));
        }
        if modal_damping.len() != n
            || modal_stiffness
                .iter()
                .chain(&modal_damping)
                .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(Error::Invalid("passive original modal diagonals"));
        }
        for l in &links {
            if l.column.len() != n
                || l.column.iter().any(|x| !x.is_finite())
                || !l.stiffness_n_per_m.is_finite()
                || l.stiffness_n_per_m < 0.0
                || !l.damping_n_s_per_m.is_finite()
                || l.damping_n_s_per_m < 0.0
                || !l.rest_extension_m.is_finite()
            {
                return Err(Error::Invalid("passive complete original connection"));
            }
        }
        for pair in &columns {
            for b in pair {
                if b.len() != n || b.iter().any(|x| !x.is_finite()) {
                    return Err(Error::Invalid("complete two-tangent columns"));
                }
            }
        }
        let entries = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(2))
            .and_then(|x| x.checked_add(n.checked_mul(8 + 2 * m + links.len())?))
            .and_then(|x| x.checked_add(4 * links.len() + m * m + m))
            .ok_or(Error::Budget("vector workspace extent"))?;
        let bytes = entries
            .checked_mul(core::mem::size_of::<f64>())
            .and_then(|x| x.checked_add((n + m).checked_mul(core::mem::size_of::<usize>())?))
            .and_then(|x| x.checked_add(m.checked_mul(core::mem::size_of::<Vec<f64>>())?))
            .and_then(|x| {
                x.checked_add(
                    links
                        .len()
                        .checked_mul(core::mem::size_of::<LinearNetworkLink>())?,
                )
            })
            .and_then(|x| x.checked_add(core::mem::size_of::<Self>()))
            .ok_or(Error::Budget("vector workspace bytes"))?;
        if bytes > config.maximum_workspace_bytes {
            return Err(Error::Budget("vector workspace bytes"));
        }
        let columns: Vec<Vec<f64>> = columns.into_iter().flatten().collect();
        let mut gram = [0.0; ROWS * ROWS];
        let mut lengths = [0.0; ROWS];
        for i in 0..m {
            lengths[i] = det::sqrt(dot(&columns[i], &columns[i])?);
            if lengths[i] <= 0.0 {
                return Err(Error::Invalid("nonzero tangent column"));
            }
            for j in 0..=i {
                let value = finite(
                    columns[i]
                        .iter()
                        .zip(&columns[j])
                        .map(|(a, b)| (a / lengths[i]) * (b / lengths[j]))
                        .sum(),
                )?;
                gram[i * m + j] = value;
                gram[j * m + i] = value;
            }
        }
        for i in 0..m {
            for j in 0..=i {
                let mut value = gram[i * m + j];
                for k in 0..j {
                    value -= gram[i * m + k] * gram[j * m + k];
                }
                if i == j {
                    if !value.is_finite() || value <= config.minimum_column_pivot {
                        return Err(Error::Invalid("independent conditioned tangent columns"));
                    }
                    gram[i * m + j] = det::sqrt(value);
                } else {
                    gram[i * m + j] = finite(value / gram[j * m + j])?;
                }
            }
        }
        let link_count = links.len();
        Ok(Self {
            stiffness: modal_stiffness,
            damping: modal_damping,
            links,
            columns,
            config,
            matrix: vec![0.0; n * n],
            mechanical: LuWorkspace::new(n)?,
            contact_solve: LuWorkspace::new(m)?,
            response_columns: vec![0.0; n * m],
            cached_h: None,
            mobility: [0.0; ROWS * ROWS],
            lipschitz: 0.0,
            q: vec![0.0; n],
            v: vec![0.0; n],
            next_v: vec![0.0; n],
            free_v: vec![0.0; n],
            rhs: vec![0.0; n],
            link_losses: vec![0.0; link_count],
            link_numerical_losses: vec![0.0; link_count],
            link_forces: vec![0.0; link_count],
            link_force_bounds: vec![0.0; link_count],
            force: [0.0; ROWS],
            trial_force: [0.0; ROWS],
            free_slip: [0.0; ROWS],
            residual: [0.0; ROWS],
            trial_residual: [0.0; ROWS],
            jacobian: [0.0; ROWS * ROWS],
            direction: [0.0; ROWS],
        })
    }

    pub fn contact_count(&self) -> usize {
        self.columns.len() / 2
    }
    pub const fn config(&self) -> VectorLinearNetworkConfig {
        self.config
    }
    pub fn positions(&self) -> &[f64] {
        &self.q
    }
    pub fn velocities(&self) -> &[f64] {
        &self.v
    }
    pub fn connection_dissipation_j(&self) -> &[f64] {
        &self.link_losses
    }
    pub fn connection_numerical_dissipation_j(&self) -> &[f64] {
        &self.link_numerical_losses
    }
    pub fn connection_forces_n(&self) -> &[f64] {
        &self.link_forces
    }
    /// Maxima over the actual discrete cell endpoints, including the initial state.
    pub fn connection_force_bounds_n(&self) -> &[f64] {
        &self.link_force_bounds
    }

    pub fn advance(
        &mut self,
        q0: &[f64],
        v0: &[f64],
        external: &[f64],
        dt: f64,
        capacities: &[f64],
        gate: &CancelGate,
    ) -> Result<VectorLinearNetworkStep, Error> {
        poll(gate)?;
        let n = self.q.len();
        let count = self.contact_count();
        if q0.len() != n
            || v0.len() != n
            || external.len() != n
            || capacities.len() != count
            || q0.iter().chain(v0).chain(external).any(|x| !x.is_finite())
            || capacities.iter().any(|x| !x.is_finite() || *x < 0.0)
            || !dt.is_finite()
            || dt <= 0.0
            || dt > self.config.maximum_step_s
        {
            return Err(Error::Invalid("original state, capacities and step"));
        }
        let cells_f = (dt / self.config.maximum_cell_s).ceil().max(1.0);
        if cells_f > self.config.maximum_cells as f64 {
            return Err(Error::Budget("vector trajectory cells"));
        }
        let cells = cells_f as usize;
        let h = finite(dt / cells as f64)?;
        if h <= 0.0 {
            return Err(Error::Budget("representable vector cell"));
        }
        self.prepare(h, gate)?;
        self.q.copy_from_slice(q0);
        self.v.copy_from_slice(v0);
        self.check_state()?;
        self.force.fill(0.0);
        self.link_losses.fill(0.0);
        self.link_numerical_losses.fill(0.0);
        self.link_force_bounds.fill(0.0);
        let (em, el) = self.energy(q0, v0)?;
        let mut step = VectorLinearNetworkStep {
            contact_count: count,
            regimes: [VectorFrictionRegime::Inactive; MAX_VECTOR_NETWORK_CONTACTS],
            initial_energy_j: finite(em + el)?,
            modal_energy_j: 0.0,
            connection_energy_j: 0.0,
            external_work_j: 0.0,
            modal_damping_dissipation_j: 0.0,
            connection_damping_dissipation_j: 0.0,
            friction_dissipation_j: 0.0,
            numerical_dissipation_j: 0.0,
            numerical_work_uncertainty_j: 0.0,
            contact_friction_dissipation_j: [0.0; MAX_VECTOR_NETWORK_CONTACTS],
            accumulated_slip_m: [0.0; MAX_VECTOR_NETWORK_CONTACTS],
            maximum_contact_speed_m_s: [0.0; MAX_VECTOR_NETWORK_CONTACTS],
            maximum_contact_force_n: [0.0; MAX_VECTOR_NETWORK_CONTACTS],
            contact_impulse_n_s: [[0.0; 2]; MAX_VECTOR_NETWORK_CONTACTS],
            contact_time_moment_n_s2: [[0.0; 2]; MAX_VECTOR_NETWORK_CONTACTS],
            energy_residual_j: 0.0,
            energy_tolerance_j: 0.0,
            cells,
            nonlinear_iterations: 0,
        };
        for j in 0..count {
            step.maximum_contact_speed_m_s[j] = norm2([
                dot(&self.columns[2 * j], v0)?,
                dot(&self.columns[2 * j + 1], v0)?,
            ]);
        }
        self.update_link_forces()?;
        for cell in 0..cells {
            poll(gate)?;
            self.cell(external, h, capacities, cell as f64 * h, &mut step, gate)?;
            self.check_state()?;
            self.update_link_forces()?;
            if step.numerical_dissipation_j > self.config.maximum_numerical_dissipation_j {
                return Err(Error::Budget("numerical time-discretization dissipation"));
            }
        }
        let (em, el) = self.energy(&self.q, &self.v)?;
        step.modal_energy_j = em;
        step.connection_energy_j = el;
        step.connection_damping_dissipation_j = finite(self.link_losses.iter().sum())?;
        step.friction_dissipation_j = finite(step.contact_friction_dissipation_j.iter().sum())?;
        step.external_work_j = finite(
            external
                .iter()
                .zip(self.q.iter().zip(q0))
                .map(|(f, (q, q0))| f * (q - q0))
                .sum(),
        )?;
        let losses = step.modal_damping_dissipation_j
            + step.connection_damping_dissipation_j
            + step.friction_dissipation_j
            + step.numerical_dissipation_j;
        let scale = step
            .initial_energy_j
            .abs()
            .max((em + el).abs())
            .max(step.external_work_j.abs())
            .max(losses);
        step.energy_residual_j =
            finite(em + el - step.initial_energy_j - step.external_work_j + losses)?;
        step.energy_tolerance_j = finite(
            self.config.energy_absolute_tolerance_j + self.config.energy_relative_tolerance * scale,
        )?;
        if step.numerical_work_uncertainty_j > self.config.maximum_numerical_work_uncertainty_j
            || step.energy_residual_j.abs() + step.numerical_work_uncertainty_j
                > step.energy_tolerance_j
        {
            return Err(Error::EnergyAccuracy {
                residual_j: step.energy_residual_j,
                estimated_error_j: step.numerical_work_uncertainty_j,
                tolerance_j: step
                    .energy_tolerance_j
                    .min(self.config.maximum_numerical_work_uncertainty_j),
            });
        }
        poll(gate)?;
        Ok(step)
    }

    fn energy(&self, q: &[f64], v: &[f64]) -> Result<(f64, f64), Error> {
        let modal = finite(
            q.iter()
                .zip(v)
                .zip(&self.stiffness)
                .map(|((q, v), k)| 0.5 * (v * v + k * q * q))
                .sum(),
        )?;
        let mut links = 0.0;
        for l in &self.links {
            let z = dot(&l.column, q)? - l.rest_extension_m;
            links = finite(links + 0.5 * l.stiffness_n_per_m * z * z)?;
        }
        Ok((modal, links))
    }
    fn check_state(&self) -> Result<(), Error> {
        if self
            .q
            .iter()
            .any(|q| q.abs() > self.config.maximum_abs_modal_coordinate)
            || self
                .v
                .iter()
                .any(|v| v.abs() > self.config.maximum_abs_modal_velocity)
        {
            return Err(Error::Budget("physical discrete modal trajectory"));
        }
        Ok(())
    }
    fn update_link_forces(&mut self) -> Result<(), Error> {
        for (i, l) in self.links.iter().enumerate() {
            let f = finite(
                -l.stiffness_n_per_m * (dot(&l.column, &self.q)? - l.rest_extension_m)
                    - l.damping_n_s_per_m * dot(&l.column, &self.v)?,
            )?;
            self.link_forces[i] = f;
            self.link_force_bounds[i] = self.link_force_bounds[i].max(f.abs());
        }
        Ok(())
    }
}

fn validate_config(c: VectorLinearNetworkConfig) -> Result<(), Error> {
    if c.maximum_modes == 0
        || c.maximum_contacts == 0
        || c.maximum_contacts > MAX_VECTOR_NETWORK_CONTACTS
        || c.maximum_connections == usize::MAX
        || c.maximum_workspace_bytes == 0
        || c.maximum_cells == 0
        || c.maximum_nonlinear_iterations == 0
        || c.maximum_backtracks == 0
        || [
            c.maximum_step_s,
            c.maximum_cell_s,
            c.minimum_column_pivot,
            c.force_absolute_tolerance_n,
            c.velocity_absolute_tolerance_m_s,
            c.maximum_abs_modal_coordinate,
            c.maximum_abs_modal_velocity,
            c.energy_absolute_tolerance_j,
            c.maximum_numerical_work_uncertainty_j,
            c.maximum_numerical_dissipation_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || c.minimum_column_pivot >= 1.0
        || [c.force_relative_tolerance, c.energy_relative_tolerance]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0 || *x >= 1.0)
    {
        return Err(Error::Invalid(
            "positive finite vector accuracy and work budgets",
        ));
    }
    Ok(())
}
fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Nonfinite)
    }
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, Error> {
    finite(a.iter().zip(b).map(|(a, b)| a * b).sum())
}
fn norm2(x: [f64; 2]) -> f64 {
    x[0].hypot(x[1])
}
