//! Exact held-force linear motion of the complete admitted solid pencil.
//!
//! Every free coordinate in the original P1, TET10 or HEX assembly is retained.
//! The shared `fs-phs` quadratic-flow owner prepares one exponential for a fixed
//! clock. No modal truncation, diagonal mass, frequency retuning or independent
//! mode assumption is introduced. Warm intervals apply that same cached map.
//!
//! The state is [u,v], with H=(u^T K u+v^T M v)/2, Q=diag(K,M),
//! J=[0,M^-1;-M^-1,0] and G=[0;M^-1]. The full consistent mass is inverted once
//! through native Cholesky and admitted against the original matrix. Held force
//! work is independently f dot (u_after-u_before). It is never inferred from H.
//! Fixed supports are homogeneous and stationary; the reduced assembly cannot
//! supply constrained reactions or moving-support coupling. Force inputs and
//! outputs explicitly use `free_dofs()` order. Material damping, nonlinear
//! finite rotations and incomplete modal approximations have separate owners.

use crate::linear3::TetElasticAssembly;
use fs_exec::Cx;
use fs_la::factor::{FactorError, cholesky};
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, split_quadratic,
};
use fs_phs::{PhsError, PortHamiltonian, PreparedStepError, QuadraticStorage};
use fs_sparse::Csr;

/// Explicit complete-state capacity and numerical allowances. No defaults.
#[derive(Clone, Copy, Debug)]
pub struct ExactLinearSolidControls {
    pub maximum_free_dofs: usize,
    /// Conservative bound for retained matrices and simultaneous cold scratch.
    pub maximum_preparation_bytes: usize,
    /// The native exponential's declared scalar-product work cap.
    pub maximum_exponential_multiply_terms: u64,
    /// Admits source symmetry and max entry of M M^-1 - I; dimensionless.
    pub maximum_relative_operator_error: f64,
    /// Absolute allowance for each interval and cumulative work balance [J].
    pub maximum_work_error_j: f64,
    pub maximum_steps: u64,
}

#[derive(Debug)]
pub enum ExactLinearSolidError {
    Invalid(&'static str),
    Budget { required: usize, maximum: usize },
    Mass(FactorError),
    Operator { defect: f64, allowance: f64 },
    Propagation(PhsError),
    Energy { defect_j: f64, allowance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for ExactLinearSolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "exact linear solid refused: {self:?}")
    }
}
impl core::error::Error for ExactLinearSolidError {}
impl From<PhsError> for ExactLinearSolidError {
    fn from(error: PhsError) -> Self {
        Self::Propagation(error)
    }
}
impl From<PreparedStepError> for ExactLinearSolidError {
    fn from(error: PreparedStepError) -> Self {
        match error {
            PreparedStepError::Cancelled => Self::Cancelled,
            PreparedStepError::Solver(error) => Self::Propagation(error),
        }
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), ExactLinearSolidError> {
    cx.checkpoint()
        .map_err(|_| ExactLinearSolidError::Cancelled)
}
fn invalid(what: &'static str) -> ExactLinearSolidError {
    ExactLinearSolidError::Invalid(what)
}

/// Original physical stored/kinetic energies and signed applied work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactLinearSolidEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub supplied_work_j: f64,
    pub balance_error_j: f64,
    pub cumulative_balance_error_j: f64,
}
#[derive(Clone, Debug)]
pub struct ExactLinearSolidStep {
    pub step: u64,
    pub time_s: f64,
    /// Work-conjugate interval velocity Delta u / dt, in free-DOF order [m/s].
    pub mean_velocity_mps: Vec<f64>,
    pub energy: ExactLinearSolidEnergy,
}

/// Retained complete-coordinate motion at one explicitly fixed positive clock.
/// Preparation is dense and bounded; this owner does not claim a scalable
/// sparse exponential. Exactness concerns the admitted spatial pencil and held
/// force, not finite-element continuum accuracy or time-varying input accuracy.
pub struct ExactLinearSolid {
    flow: PreparedQuadraticFlow,
    stiffness: Vec<f64>,
    mass: Vec<f64>,
    free_dofs: Vec<usize>,
    state: Vec<f64>,
    dt_s: f64,
    controls: ExactLinearSolidControls,
    steps: u64,
    initial_energy_j: f64,
    cumulative_work_j: f64,
    work_correction_j: f64,
}

impl ExactLinearSolid {
    /// Bind the full source pencil and explicit initial free-coordinate fields.
    /// Both stationary support values and velocities must be zero in the source
    /// model. This API has no full nodal-force input that could silently discard
    /// constrained loads: each force entry belongs to exactly one free row.
    pub fn prepare(
        assembly: &TetElasticAssembly,
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: ExactLinearSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, ExactLinearSolidError> {
        poll(cx)?;
        let n = assembly.free_dofs.len();
        if n == 0
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || controls.maximum_steps == 0
            || controls.maximum_steps > (1_u64 << 53)
            || !controls.maximum_relative_operator_error.is_finite()
            || !(0.0..1.0).contains(&controls.maximum_relative_operator_error)
            || !controls.maximum_work_error_j.is_finite()
            || controls.maximum_work_error_j < 0.0
            || initial_displacement_m.len() != n
            || initial_velocity_mps.len() != n
            || initial_displacement_m
                .iter()
                .chain(initial_velocity_mps)
                .any(|x| !x.is_finite())
            || assembly.free_dofs.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(invalid(
                "complete finite initial state, sorted free map, positive clock and valid controls required",
            ));
        }
        if n > controls.maximum_free_dofs {
            return Err(ExactLinearSolidError::Budget {
                required: n,
                maximum: controls.maximum_free_dofs,
            });
        }
        let d = n
            .checked_mul(3)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| invalid("complete exponential dimension overflow"))?;
        if d > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
            return Err(ExactLinearSolidError::Budget {
                required: d,
                maximum: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
            });
        }
        // Native flow reserves 24*d^2 doubles; another 24*d^2 covers original
        // pencil, inverse, Q/J/G, eigen preparation and retained state. Cholesky
        // above its panel size additionally packs the native GEMM panels.
        let bytes = d
            .checked_mul(d)
            .and_then(|v| v.checked_mul(48 * 8))
            .and_then(|v| v.checked_add(if n > 32 { 2 * 1024 * 1024 } else { 0 }))
            .ok_or_else(|| invalid("complete exponential byte overflow"))?;
        if bytes > controls.maximum_preparation_bytes {
            return Err(ExactLinearSolidError::Budget {
                required: bytes,
                maximum: controls.maximum_preparation_bytes,
            });
        }
        let tolerance = controls.maximum_relative_operator_error;
        let stiffness = symmetric_matrix(&assembly.stiffness, n, tolerance, cx)?;
        let mass = symmetric_matrix(&assembly.mass, n, tolerance, cx)?;
        let factor = cholesky(&mass, n).map_err(ExactLinearSolidError::Mass)?;
        let mut inverse = vec![0.0; n * n];
        let mut column = vec![0.0; n];
        for j in 0..n {
            poll(cx)?;
            column.fill(0.0);
            column[j] = 1.0;
            factor.solve(&mut column);
            for i in 0..n {
                inverse[i * n + j] = column[i];
            }
        }
        // The represented inverse must be exactly symmetric for J to be skew.
        // Its full residual is tested AFTER this rounding-level symmetrization.
        for i in 0..n {
            for j in 0..i {
                let a = 0.5 * inverse[i * n + j] + 0.5 * inverse[j * n + i];
                inverse[i * n + j] = a;
                inverse[j * n + i] = a;
            }
        }
        let mut defect = 0.0_f64;
        for i in 0..n {
            poll(cx)?;
            for j in 0..n {
                let value = (0..n).fold(0.0, |s, k| mass[i * n + k].mul_add(inverse[k * n + j], s));
                let error = (value - if i == j { 1.0 } else { 0.0 }).abs();
                if !error.is_finite() {
                    return Err(invalid("unrepresentable inverse mass"));
                }
                defect = defect.max(error);
            }
        }
        if defect > tolerance {
            return Err(ExactLinearSolidError::Operator {
                defect,
                allowance: tolerance,
            });
        }
        let state_dim = 2 * n;
        let mut q = vec![0.0; state_dim * state_dim];
        let mut j = vec![0.0; state_dim * state_dim];
        let mut g = vec![0.0; state_dim * n];
        for a in 0..n {
            for b in 0..n {
                q[a * state_dim + b] = stiffness[a * n + b];
                q[(n + a) * state_dim + n + b] = mass[a * n + b];
                j[a * state_dim + n + b] = inverse[a * n + b];
                j[(n + a) * state_dim + b] = -inverse[a * n + b];
                g[(n + a) * n + b] = inverse[a * n + b];
            }
        }
        let system = PortHamiltonian::new(
            state_dim,
            n,
            j,
            vec![0.0; state_dim * state_dim],
            g,
            Box::new(QuadraticStorage::new(q, state_dim)?),
        )?;
        // Frequency matching is unused: exact propagation has no Nyquist
        // retuning/band switch. All state coordinates, including zero-K, remain.
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: state_dim,
                max_omega_dt: 3.0,
            },
        )?;
        poll(cx)?;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            dt_s,
            &[],
            &[],
            DenseExponentialBudget {
                max_dimension: d,
                max_dense_bytes: bytes / 2,
                max_multiply_terms: controls.maximum_exponential_multiply_terms,
            },
            || cx.checkpoint().is_err(),
        )?;
        let mut state = Vec::with_capacity(state_dim);
        state.extend_from_slice(initial_displacement_m);
        state.extend_from_slice(initial_velocity_mps);
        let initial_energy_j = quadratic(&stiffness, &state[..n]) + quadratic(&mass, &state[n..]);
        if !initial_energy_j.is_finite() || initial_energy_j < 0.0 {
            return Err(invalid("invalid initial physical energy"));
        }
        poll(cx)?;
        Ok(Self {
            flow,
            stiffness,
            mass,
            free_dofs: assembly.free_dofs.clone(),
            state,
            dt_s,
            controls,
            steps: 0,
            initial_energy_j,
            cumulative_work_j: 0.0,
            work_correction_j: 0.0,
        })
    }

    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    pub fn displacement_m(&self) -> &[f64] {
        &self.state[..self.free_dofs.len()]
    }
    pub fn velocity_mps(&self) -> &[f64] {
        &self.state[self.free_dofs.len()..]
    }
    pub fn time_s(&self) -> f64 {
        self.steps as f64 * self.dt_s
    }
    pub fn duration_s(&self) -> f64 {
        self.dt_s
    }
    pub fn accepted_steps(&self) -> u64 {
        self.steps
    }
    pub fn stored_energy_j(&self) -> f64 {
        quadratic(&self.stiffness, self.displacement_m())
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        quadratic(&self.mass, self.velocity_mps())
    }

    /// One complete exact held-force interval. Every force is in N, in the
    /// source free-DOF order. Cancellation/refusal leaves all accepted fields.
    pub fn advance(
        &mut self,
        held_force_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep, ExactLinearSolidError> {
        poll(cx)?;
        let n = self.free_dofs.len();
        let step = self
            .steps
            .checked_add(1)
            .filter(|v| *v <= self.controls.maximum_steps)
            .ok_or_else(|| invalid("accepted step budget exhausted"))?;
        let time_s = step as f64 * self.dt_s;
        if held_force_n.len() != n
            || held_force_n.iter().any(|f| !f.is_finite())
            || !time_s.is_finite()
            || time_s <= self.time_s()
        {
            return Err(invalid(
                "finite complete held forces and representable clock required",
            ));
        }
        let candidate = self.flow.evaluate(&self.state, held_force_n)?.state;
        poll(cx)?;
        let u = &candidate[..n];
        let v = &candidate[n..];
        let stored_before_j = self.stored_energy_j();
        let kinetic_before_j = self.kinetic_energy_j();
        let stored_after_j = quadratic(&self.stiffness, u);
        let kinetic_after_j = quadratic(&self.mass, v);
        let supplied_work_j =
            (0..n).fold(0.0, |s, i| held_force_n[i].mul_add(u[i] - self.state[i], s));
        let change = quadratic_change(&self.stiffness, self.displacement_m(), u)
            + quadratic_change(&self.mass, self.velocity_mps(), v);
        let increment = supplied_work_j - self.work_correction_j;
        let cumulative_work_j = self.cumulative_work_j + increment;
        let correction = (cumulative_work_j - self.cumulative_work_j) - increment;
        let energy = ExactLinearSolidEnergy {
            stored_before_j,
            stored_after_j,
            kinetic_before_j,
            kinetic_after_j,
            supplied_work_j,
            balance_error_j: supplied_work_j - change,
            cumulative_balance_error_j: cumulative_work_j
                - (stored_after_j + kinetic_after_j - self.initial_energy_j),
        };
        if [
            stored_before_j,
            stored_after_j,
            kinetic_before_j,
            kinetic_after_j,
            supplied_work_j,
            cumulative_work_j,
            correction,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable physical energy/work"));
        }
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > self.controls.maximum_work_error_j {
                return Err(ExactLinearSolidError::Energy {
                    defect_j,
                    allowance_j: self.controls.maximum_work_error_j,
                });
            }
        }
        let mean_velocity_mps = (0..n)
            .map(|i| (u[i] - self.state[i]) / self.dt_s)
            .collect::<Vec<_>>();
        if mean_velocity_mps.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable mean velocity"));
        }
        poll(cx)?;
        self.state = candidate;
        self.steps = step;
        self.cumulative_work_j = cumulative_work_j;
        self.work_correction_j = correction;
        Ok(ExactLinearSolidStep {
            step,
            time_s,
            mean_velocity_mps,
            energy,
        })
    }
}

fn symmetric_matrix(
    matrix: &Csr,
    n: usize,
    tolerance: f64,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, ExactLinearSolidError> {
    if matrix.nrows() != n || matrix.ncols() != n {
        return Err(invalid(
            "source matrix dimensions must match every free row",
        ));
    }
    let mut dense = vec![0.0; n * n];
    for i in 0..n {
        poll(cx)?;
        let (columns, values) = matrix.row(i);
        for (&j, &v) in columns.iter().zip(values) {
            if !v.is_finite() {
                return Err(invalid("source matrix must be finite"));
            }
            dense[i * n + j] = v;
        }
    }
    let scale = dense.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    for i in 0..n {
        for j in 0..i {
            let defect = if scale == 0.0 {
                0.0
            } else {
                (dense[i * n + j] / scale - dense[j * n + i] / scale).abs()
            };
            if defect > tolerance {
                return Err(ExactLinearSolidError::Operator {
                    defect,
                    allowance: tolerance,
                });
            }
            let value = 0.5 * dense[i * n + j] + 0.5 * dense[j * n + i];
            dense[i * n + j] = value;
            dense[j * n + i] = value;
        }
    }
    Ok(dense)
}
fn quadratic(matrix: &[f64], x: &[f64]) -> f64 {
    let n = x.len();
    0.5 * (0..n).fold(0.0, |s, i| {
        x[i].mul_add(
            (0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(x[j], s)),
            s,
        )
    })
}
fn quadratic_change(matrix: &[f64], a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    (0..n).fold(0.0, |s, i| {
        ((a[i] + b[i]) * 0.5).mul_add(
            (0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(b[j] - a[j], s)),
            s,
        )
    })
}
