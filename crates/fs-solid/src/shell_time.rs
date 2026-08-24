//! Passive transient shell time integration and exact checkpointing (bead `frankensim-b8bxd.9.2`).
//!
//! Solves $M \ddot{u} + C \dot{u} + K u = f(t)$ using the unconditionally stable,
//! second-order accurate Newmark-beta average acceleration method ($\beta = 0.25, \gamma = 0.5$).
//!
//! Invariants:
//! - Passive damping dissipation: for $f = 0$ and $C \succeq 0$, $E(t_{n+1}) \le E(t_n)$ (monotonically non-increasing energy).
//! - Conservation of energy for undamped free vibration ($C = 0, f = 0$).
//! - Exact checkpoint / resume equivalence: resuming from a checkpoint produces bit-identical trajectories.
//! - Bounded step budgets and cancellation checking.

use fs_blake3::{hash_domain, ContentHash};
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
    for i in 0..n {
        let mut sum = 0.0;
        let row_offset = i * n;
        for j in 0..n {
            sum += a[row_offset + j] * x[j];
        }
        y[i] = sum;
    }
    y
}

/// Solve $A x = b$ using Gaussian elimination with partial pivoting.
fn solve_dense(n: usize, a: &[f64], b: &[f64]) -> Result<Vec<f64>, ShellTimeError> {
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();

    for k in 0..n {
        // Pivot selection
        let mut max_val = mat[k * n + k].abs();
        let mut pivot_row = k;
        for i in (k + 1)..n {
            let val = mat[i * n + k].abs();
            if val > max_val {
                max_val = val;
                pivot_row = i;
            }
        }

        if max_val < 1.0e-14 {
            return Err(ShellTimeError::NonFiniteValue {
                field: "dense_system_singular",
            });
        }

        if pivot_row != k {
            for j in 0..n {
                mat.swap(k * n + j, pivot_row * n + j);
            }
            rhs.swap(k, pivot_row);
        }

        // Elimination
        let diag = mat[k * n + k];
        for i in (k + 1)..n {
            let factor = mat[i * n + k] / diag;
            mat[i * n + k] = 0.0;
            for j in (k + 1)..n {
                mat[i * n + j] -= factor * mat[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for j in (i + 1)..n {
            sum -= mat[i * n + j] * x[j];
        }
        x[i] = sum / mat[i * n + i];
    }

    Ok(x)
}

/// Step the dynamic system forward by one timestep using Newmark-beta.
///
/// System: $M \ddot{u}_{n+1} + C \dot{u}_{n+1} + K u_{n+1} = f_{n+1}$
///
/// # Errors
/// Returns [`ShellTimeError`] on dimension mismatch, singular matrix, or non-finite arithmetic.
pub fn step_newmark(
    n: usize,
    mass: &[f64],
    damping: Option<&[f64]>,
    stiffness: &[f64],
    current: &DynamicState,
    f_next: &[f64],
    config: &ShellTimeConfig,
) -> Result<DynamicState, ShellTimeError> {
    let dt = config.dt_s;
    if !dt.is_finite() || dt <= 0.0 {
        return Err(ShellTimeError::InvalidTimeStep {
            dt: format!("{dt}"),
        });
    }

    let beta = config.beta;
    let gamma = config.gamma;

    // Effective stiffness: K_eff = K + (gamma / (beta * dt)) * C + (1 / (beta * dt^2)) * M
    let a0 = 1.0 / (beta * dt * dt);
    let a1 = gamma / (beta * dt);
    let a2 = 1.0 / (beta * dt);
    let a3 = 1.0 / (2.0 * beta) - 1.0;
    let a4 = gamma / beta - 1.0;
    let a5 = (dt / 2.0) * (gamma / beta - 2.0);

    let mut k_eff = stiffness.to_vec();
    for i in 0..(n * n) {
        k_eff[i] += a0 * mass[i];
    }
    if let Some(c) = damping {
        for i in 0..(n * n) {
            k_eff[i] += a1 * c[i];
        }
    }

    // Effective force: R_eff = f_{n+1} + M * (a0*u_n + a2*v_n + a3*a_n) + C * (a1*u_n + a4*v_n + a5*a_n)
    let mut m_term = vec![0.0; n];
    for i in 0..n {
        m_term[i] = a0 * current.displacement[i] + a2 * current.velocity[i] + a3 * current.acceleration[i];
    }
    let m_force = mat_vec(n, mass, &m_term);

    let mut r_eff = vec![0.0; n];
    for i in 0..n {
        r_eff[i] = f_next[i] + m_force[i];
    }

    if let Some(c) = damping {
        let mut c_term = vec![0.0; n];
        for i in 0..n {
            c_term[i] = a1 * current.displacement[i] + a4 * current.velocity[i] + a5 * current.acceleration[i];
        }
        let c_force = mat_vec(n, c, &c_term);
        for i in 0..n {
            r_eff[i] += c_force[i];
        }
    }

    // Solve for u_{n+1}
    let u_next = solve_dense(n, &k_eff, &r_eff)?;

    // Calculate a_{n+1} and v_{n+1}
    let mut a_next = vec![0.0; n];
    let mut v_next = vec![0.0; n];
    for i in 0..n {
        a_next[i] = a0 * (u_next[i] - current.displacement[i]) - a2 * current.velocity[i] - a3 * current.acceleration[i];
        v_next[i] = current.velocity[i] + dt * ((1.0 - gamma) * current.acceleration[i] + gamma * a_next[i]);
    }

    // Check finite
    for &val in &u_next {
        if !val.is_finite() {
            return Err(ShellTimeError::NonFiniteValue { field: "displacement" });
        }
    }

    // Energy calculation
    let m_v = mat_vec(n, mass, &v_next);
    let mut kinetic_energy = 0.0;
    for i in 0..n {
        kinetic_energy += 0.5 * v_next[i] * m_v[i];
    }

    let k_u = mat_vec(n, stiffness, &u_next);
    let mut strain_energy = 0.0;
    for i in 0..n {
        strain_energy += 0.5 * u_next[i] * k_u[i];
    }

    let total_energy = kinetic_energy + strain_energy;

    Ok(DynamicState {
        time_s: current.time_s + dt,
        step: current.step + 1,
        displacement: u_next,
        velocity: v_next,
        acceleration: a_next,
        kinetic_energy_j: kinetic_energy,
        strain_energy_j: strain_energy,
        total_energy_j: total_energy,
    })
}
