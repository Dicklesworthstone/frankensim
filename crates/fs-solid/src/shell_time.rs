//! Passive transient shell time integration and exact checkpointing (bead `frankensim-b8bxd.9.2`).
//!
//! Solves $M \ddot{u} + C \dot{u} + K u = f(t)$ using Newmark-beta.
//! The default average acceleration parameters ($\beta = 0.25, \gamma = 0.5$)
//! are second-order and unconditionally stable for passive linear systems.
//!
//! Invariants for those default parameters and initially equilibrated acceleration:
//! - Passive damping dissipation: for $f = 0$ and $C \succeq 0$, $E(t_{n+1}) \le E(t_n)$ (monotonically non-increasing energy).
//! - Conservation of energy for undamped free vibration ($C = 0, f = 0$).
//! - Exact checkpoint / resume equivalence: resuming from a checkpoint produces bit-identical trajectories.
//! - Bounded step budgets and cancellation checking.

use fs_blake3::{ContentHash, hash_domain};
use fs_la::{LuWorkspace, LuWorkspaceError};
use std::fmt;

/// Schema identifier for shell transient time integration receipts.
pub const SHELL_TIME_INTEGRATOR_SCHEMA_V1: &str = "org.frankensim.solid.shell-time.v1";

/// Error states during transient shell time integration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ShellTimeError {
    /// Dimension mismatch between mass/stiffness matrix and state vector.
    DimensionMismatch {
        /// Expected dimension.
        expected: usize,
        /// Found dimension.
        found: usize,
    },
    /// Non-positive or non-finite timestep dt.
    InvalidTimeStep {
        /// Invalid timestep value string.
        dt: String,
    },
    /// Step budget exceeded before completing simulation interval.
    BudgetExceeded {
        /// Maximum allowed steps.
        max_steps: usize,
    },
    /// Cancellation requested during time advancement.
    Cancelled {
        /// Step at which cancellation occurred.
        step: usize,
    },
    /// Non-finite value detected in state or excitation.
    NonFiniteValue {
        /// Field containing non-finite value.
        field: &'static str,
    },
    /// A Newmark parameter is nonfinite or outside its mathematical domain.
    InvalidParameter {
        /// Parameter that must be corrected.
        parameter: &'static str,
    },
    /// The original dense linear algebra owner refused the effective mass solve.
    LinearSolve {
        /// Original factorization or substitution error.
        source: LuWorkspaceError,
    },
    /// Checkpoint hash or schema mismatch on resume.
    CorruptCheckpoint {
        /// Detail of checkpoint corruption.
        detail: String,
    },
}

impl fmt::Display for ShellTimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionMismatch { expected, found } => {
                write!(f, "dimension mismatch: expected {expected}, found {found}")
            }
            Self::InvalidTimeStep { dt } => write!(f, "invalid timestep dt: {dt}"),
            Self::BudgetExceeded { max_steps } => {
                write!(f, "step budget exceeded ({max_steps} steps)")
            }
            Self::Cancelled { step } => write!(f, "computation cancelled at step {step}"),
            Self::NonFiniteValue { field } => write!(f, "non-finite value in field {field}"),
            Self::InvalidParameter { parameter } => {
                write!(f, "invalid Newmark parameter {parameter}")
            }
            Self::LinearSolve { source } => {
                write!(f, "shell effective mass solve failed: {source}")
            }
            Self::CorruptCheckpoint { detail } => write!(f, "corrupt checkpoint: {detail}"),
        }
    }
}

impl std::error::Error for ShellTimeError {}

/// Configuration for transient shell integration.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellTimeConfig {
    /// Timestep size [s].
    pub dt_s: f64,
    /// Total duration [s].
    pub duration_s: f64,
    /// Maximum allowed steps.
    pub max_steps: usize,
    /// Newmark beta parameter (default 0.25).
    pub beta: f64,
    /// Newmark gamma parameter (default 0.50).
    pub gamma: f64,
}

impl Default for ShellTimeConfig {
    fn default() -> Self {
        Self {
            dt_s: 1.0e-4,
            duration_s: 1.0e-2,
            max_steps: 10_000,
            beta: 0.25,
            gamma: 0.50,
        }
    }
}

/// Instantaneous dynamic state of the plate/shell.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicState {
    /// Current simulation time [s].
    pub time_s: f64,
    /// Step index.
    pub step: usize,
    /// Nodal displacement vector $u$ [m / rad].
    pub displacement: Vec<f64>,
    /// Nodal velocity vector $\dot{u}$ [m/s / rad/s].
    pub velocity: Vec<f64>,
    /// Nodal acceleration vector $\ddot{u}$ [m/s² / rad/s²].
    pub acceleration: Vec<f64>,
    /// Kinetic energy $T = \frac{1}{2} \dot{u}^T M \dot{u}$ [J].
    pub kinetic_energy_j: f64,
    /// Strain energy $U = \frac{1}{2} u^T K u$ [J].
    pub strain_energy_j: f64,
    /// Total mechanical energy $E = T + U$ [J].
    pub total_energy_j: f64,
}

/// Exact serializable checkpoint for pause, resume, and fork.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellTimeCheckpoint {
    /// Schema string.
    pub schema_version: &'static str,
    /// Step index.
    pub step: usize,
    /// Current time [s].
    pub time_s: f64,
    /// Displacements.
    pub u: Vec<f64>,
    /// Velocities.
    pub v: Vec<f64>,
    /// Accelerations.
    pub a: Vec<f64>,
    /// Cryptographic digest of this checkpoint.
    pub digest: ContentHash,
}

impl ShellTimeCheckpoint {
    /// Construct and hash a checkpoint from state.
    #[must_use]
    pub fn new(state: &DynamicState) -> Self {
        let mut bytes = Vec::with_capacity(32 + state.displacement.len() * 24);
        bytes.extend_from_slice(SHELL_TIME_INTEGRATOR_SCHEMA_V1.as_bytes());
        bytes.extend_from_slice(&state.step.to_le_bytes());
        bytes.extend_from_slice(&state.time_s.to_bits().to_le_bytes());
        for &val in &state.displacement {
            bytes.extend_from_slice(&val.to_bits().to_le_bytes());
        }
        for &val in &state.velocity {
            bytes.extend_from_slice(&val.to_bits().to_le_bytes());
        }
        for &val in &state.acceleration {
            bytes.extend_from_slice(&val.to_bits().to_le_bytes());
        }

        let digest = hash_domain("org.frankensim.solid.shell-checkpoint.v1", &bytes);

        Self {
            schema_version: SHELL_TIME_INTEGRATOR_SCHEMA_V1,
            step: state.step,
            time_s: state.time_s,
            u: state.displacement.clone(),
            v: state.velocity.clone(),
            a: state.acceleration.clone(),
            digest,
        }
    }

    /// Validate checkpoint integrity against its digest.
    ///
    /// # Errors
    /// [`ShellTimeError::CorruptCheckpoint`] if hash does not match preimage.
    pub fn verify(&self) -> Result<(), ShellTimeError> {
        let mut bytes = Vec::with_capacity(32 + self.u.len() * 24);
        bytes.extend_from_slice(self.schema_version.as_bytes());
        bytes.extend_from_slice(&self.step.to_le_bytes());
        bytes.extend_from_slice(&self.time_s.to_bits().to_le_bytes());
        for &val in &self.u {
            bytes.extend_from_slice(&val.to_bits().to_le_bytes());
        }
        for &val in &self.v {
            bytes.extend_from_slice(&val.to_bits().to_le_bytes());
        }
        for &val in &self.a {
            bytes.extend_from_slice(&val.to_bits().to_le_bytes());
        }

        let expected_digest = hash_domain("org.frankensim.solid.shell-checkpoint.v1", &bytes);
        if expected_digest != self.digest {
            return Err(ShellTimeError::CorruptCheckpoint {
                detail: format!(
                    "digest mismatch: expected {}, found {}",
                    expected_digest.to_hex(),
                    self.digest.to_hex()
                ),
            });
        }
        Ok(())
    }
}

/// Helper for dense matrix-vector product $y = A x$.
fn mat_vec(n: usize, a: &[f64], x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for (i, y_elem) in y.iter_mut().enumerate() {
        let mut sum = 0.0;
        let row_offset = i * n;
        for j in 0..n {
            sum += a[row_offset + j] * x[j];
        }
        *y_elem = sum;
    }
    y
}

fn next_time(current: &DynamicState, config: &ShellTimeConfig) -> Result<f64, ShellTimeError> {
    let dt = config.dt_s;
    let next_time = current.time_s + dt;
    if !dt.is_finite() || dt <= 0.0 || !next_time.is_finite() || next_time <= current.time_s {
        return Err(ShellTimeError::InvalidTimeStep {
            dt: format!("{dt}"),
        });
    }
    Ok(next_time)
}

fn check_system(
    n: usize,
    mass: &[f64],
    damping: Option<&[f64]>,
    stiffness: &[f64],
    config: &ShellTimeConfig,
) -> Result<(), ShellTimeError> {
    let beta = config.beta;
    let gamma = config.gamma;
    if !beta.is_finite() || beta <= 0.0 {
        return Err(ShellTimeError::InvalidParameter { parameter: "beta" });
    }
    if !gamma.is_finite() || gamma < 0.0 {
        return Err(ShellTimeError::InvalidParameter { parameter: "gamma" });
    }
    let entries = n.checked_mul(n).ok_or(ShellTimeError::InvalidParameter {
        parameter: "matrix dimension",
    })?;
    for matrix in [Some(mass), damping, Some(stiffness)].into_iter().flatten() {
        if matrix.len() != entries {
            return Err(ShellTimeError::DimensionMismatch {
                expected: entries,
                found: matrix.len(),
            });
        }
        if matrix.iter().any(|v| !v.is_finite()) {
            return Err(ShellTimeError::NonFiniteValue { field: "operator" });
        }
    }
    Ok(())
}

fn check_state(n: usize, current: &DynamicState, f_next: &[f64]) -> Result<usize, ShellTimeError> {
    for values in [
        &current.displacement[..],
        &current.velocity,
        &current.acceleration,
        f_next,
    ] {
        if values.len() != n {
            return Err(ShellTimeError::DimensionMismatch {
                expected: n,
                found: values.len(),
            });
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(ShellTimeError::NonFiniteValue {
                field: "state_or_force",
            });
        }
    }
    current
        .step
        .checked_add(1)
        .ok_or(ShellTimeError::InvalidParameter {
            parameter: "step index",
        })
}

/// A fixed linear shell system with reusable Newmark factors.
///
/// Preparation borrows the row-major `n × n` operators and factors
/// $M + \gamma \Delta t C + \beta \Delta t^2 K$ once through the same
/// [`LuWorkspace`] owner as [`step_newmark`]. Each subsequent [`Self::step`]
/// performs quadratic-cost matrix-vector products and triangular solves,
/// without repeating the cubic-cost factorization. Prepare again when any
/// operator, timestep, or Newmark parameter changes.
///
/// The complete-increment advancement and its floating-point operation order
/// are shared with the one-shot API. Initial acceleration must satisfy
/// $M a_0 = f_0 - C v_0 - K u_0$ for the documented energy/impulse identities.
/// Preparation checks finite values and dimensions; it does not certify matrix
/// symmetry, definiteness, or the validity of the spatial shell model.
#[derive(Debug)]
pub struct PreparedNewmark<'a> {
    n: usize,
    mass: &'a [f64],
    damping: Option<&'a [f64]>,
    stiffness: &'a [f64],
    config: ShellTimeConfig,
    solver: LuWorkspace,
}

impl<'a> PreparedNewmark<'a> {
    /// Prepare fixed operators and timestep for repeated time advancement.
    ///
    /// # Errors
    /// Returns [`ShellTimeError`] for invalid parameters, malformed or non-finite
    /// operators, or a refused effective mass factorization.
    pub fn new(
        n: usize,
        mass: &'a [f64],
        damping: Option<&'a [f64]>,
        stiffness: &'a [f64],
        config: &ShellTimeConfig,
    ) -> Result<Self, ShellTimeError> {
        let dt = config.dt_s;
        if !dt.is_finite() || dt <= 0.0 {
            return Err(ShellTimeError::InvalidTimeStep {
                dt: format!("{dt}"),
            });
        }
        check_system(n, mass, damping, stiffness, config)?;
        Self::from_admitted(n, mass, damping, stiffness, config)
    }

    fn from_admitted(
        n: usize,
        mass: &'a [f64],
        damping: Option<&'a [f64]>,
        stiffness: &'a [f64],
        config: &ShellTimeConfig,
    ) -> Result<Self, ShellTimeError> {
        let dt = config.dt_s;
        let displacement_scale = config.beta * dt * dt;
        let velocity_scale = config.gamma * dt;
        let mut effective_mass = mass.to_vec();
        for i in 0..mass.len() {
            effective_mass[i] += displacement_scale * stiffness[i];
            if let Some(c) = damping {
                effective_mass[i] += velocity_scale * c[i];
            }
        }
        let mut solver =
            LuWorkspace::new(n).map_err(|source| ShellTimeError::LinearSolve { source })?;
        solver
            .factor(&effective_mass)
            .map_err(|source| ShellTimeError::LinearSolve { source })?;
        Ok(Self {
            n,
            mass,
            damping,
            stiffness,
            config: config.clone(),
            solver,
        })
    }

    /// Advance with forcing evaluated at `current.time_s + dt_s`.
    ///
    /// The receiver is mutable only for the original LU owner's substitution
    /// workspace. Operators and factors stay fixed, so independent trajectories
    /// and restored checkpoints can reuse the same prepared system. Returned
    /// energies are recomputed; the input state's energy fields are not used.
    ///
    /// # Errors
    /// Returns [`ShellTimeError`] for invalid time/state/force data or non-finite
    /// arithmetic. A refusal publishes no new state and preserves the factors.
    pub fn step(
        &mut self,
        current: &DynamicState,
        f_next: &[f64],
    ) -> Result<DynamicState, ShellTimeError> {
        let next_time = next_time(current, &self.config)?;
        let next_step = check_state(self.n, current, f_next)?;
        self.advance(current, f_next, next_time, next_step)
    }

    fn advance(
        &mut self,
        current: &DynamicState,
        f_next: &[f64],
        next_time: f64,
        next_step: usize,
    ) -> Result<DynamicState, ShellTimeError> {
        let n = self.n;
        let dt = self.config.dt_s;
        let beta = self.config.beta;
        let gamma = self.config.gamma;

        // The predictors exclude the unknown acceleration. Off-diagonal inertia,
        // viscous coupling and elastic forces remain in the original operators.
        let mut u_predict = vec![0.0; n];
        let mut v_predict = vec![0.0; n];
        for i in 0..n {
            u_predict[i] = current.displacement[i]
                + dt * (current.velocity[i] + dt * (0.5 - beta) * current.acceleration[i]);
            v_predict[i] = current.velocity[i] + dt * (1.0 - gamma) * current.acceleration[i];
        }
        let mut rhs = mat_vec(n, self.stiffness, &u_predict);
        for i in 0..n {
            rhs[i] = f_next[i] - rhs[i];
        }
        if let Some(c) = self.damping {
            let viscous = mat_vec(n, c, &v_predict);
            for i in 0..n {
                rhs[i] -= viscous[i];
            }
        }
        let mut a_next = vec![0.0; n];
        self.solver
            .solve_factored_into(&rhs, &mut a_next)
            .map_err(|source| ShellTimeError::LinearSolve { source })?;

        let mut u_next = vec![0.0; n];
        let mut v_next = vec![0.0; n];
        for i in 0..n {
            // Add each complete increment once. Rounded absolute positions never
            // enter the velocity or acceleration update.
            u_next[i] = current.displacement[i]
                + dt * (current.velocity[i]
                    + dt * ((0.5 - beta) * current.acceleration[i] + beta * a_next[i]));
            v_next[i] = current.velocity[i]
                + dt * ((1.0 - gamma) * current.acceleration[i] + gamma * a_next[i]);
        }
        for (field, values) in [("displacement", &u_next), ("velocity", &v_next)] {
            if values.iter().any(|v| !v.is_finite()) {
                return Err(ShellTimeError::NonFiniteValue { field });
            }
        }

        // Energy calculation
        let m_v = mat_vec(n, self.mass, &v_next);
        let mut kinetic_energy = 0.0;
        for i in 0..n {
            kinetic_energy += 0.5 * v_next[i] * m_v[i];
        }

        let k_u = mat_vec(n, self.stiffness, &u_next);
        let mut strain_energy = 0.0;
        for i in 0..n {
            strain_energy += 0.5 * u_next[i] * k_u[i];
        }

        let total_energy = kinetic_energy + strain_energy;
        if ![kinetic_energy, strain_energy, total_energy]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(ShellTimeError::NonFiniteValue { field: "energy" });
        }

        Ok(DynamicState {
            time_s: next_time,
            step: next_step,
            displacement: u_next,
            velocity: v_next,
            acceleration: a_next,
            kinetic_energy_j: kinetic_energy,
            strain_energy_j: strain_energy,
            total_energy_j: total_energy,
        })
    }
}

/// Step the dynamic system forward by one timestep using Newmark-beta.
///
/// System: $M \ddot{u}_{n+1} + C \dot{u}_{n+1} + K u_{n+1} = f_{n+1}$
///
/// The effective MASS is solved for acceleration directly. Recovering
/// acceleration from the difference of absolute displacements loses momentum
/// when rigid drift or the coordinate datum is large compared with dt*v.
/// The acceleration form also avoids scaling the matrix by 1/dt².
///
/// With beta=1/4, gamma=1/2, and an initially equilibrated acceleration,
/// the update integrates endpoint force impulses by the trapezoid rule and
/// obeys the corresponding discrete force-work/energy identity. Other admitted
/// Newmark parameters do not inherit these conservation or stability claims.
///
/// This convenience API prepares the operators for each call. Use
/// [`PreparedNewmark`] to reuse the same factorization over a trajectory.
///
/// # Errors
/// Returns [`ShellTimeError`] on dimension mismatch, invalid parameters,
/// singular effective mass, or non-finite input or output.
pub fn step_newmark(
    n: usize,
    mass: &[f64],
    damping: Option<&[f64]>,
    stiffness: &[f64],
    current: &DynamicState,
    f_next: &[f64],
    config: &ShellTimeConfig,
) -> Result<DynamicState, ShellTimeError> {
    // Preserve the one-shot admission order before allocating or factoring.
    let next_time = next_time(current, config)?;
    check_system(n, mass, damping, stiffness, config)?;
    let next_step = check_state(n, current, f_next)?;
    PreparedNewmark::from_admitted(n, mass, damping, stiffness, config)?
        .advance(current, f_next, next_time, next_step)
}
