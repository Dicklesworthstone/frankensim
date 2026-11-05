//! Coupled normal and vector-shear inertia on the original full-stick potential.
//!
//! For q=[delta,u1,u2], the held native mobility W=J M^-1 J^T gives
//! q''=a_external-W f. The original source's symmetric discrete gradient G
//! closes G.(q1-q0)=U1-U0, with its independently reported root-rounding defect.
//! Midpoint inertia therefore solves q1-q0=h v0+h²(a_external-W G)/2.
//! Its analytic Jacobian retains every normal/shear cross derivative. One small
//! prepared LU workspace is reused; accepted history is never a Newton iterate.
//!
//! A constant discrete effort defines a quadratic position interpolation. The
//! original path owner bounds its entire Bezier hull against radius, material
//! and sourced static-junction limits, including hidden intrastep excursions.
//! It refuses slip/fold crossing; no landing state, loss or heat is invented.
//! This is one second-order conservative discrete phase, not an exact transient
//! or a global uniqueness certificate. Explicit step limits and refinement are
//! still required. A positive symmetric part of the local inertial tangent is
//! required throughout Newton continuation; it is not a continuum stability test.

/// Original-source prescribed-temperature continuation at held coordinates.
pub mod temperature;

use super::{
    JkrMixedModeResponse, SourceJkrMixedMode,
    path::{JkrMixedModePath, JkrMixedModePathError},
};
use crate::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_math::dd::Dd;

/// Physical phase bounds, shared source-work allowance and dimensional gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeMotionBudget {
    pub maximum_step_s: f64,
    /// Bounds the whole quadratic hull's distance from its initial point [m].
    pub maximum_coordinate_increment_m: f64,
    pub maximum_speed_m_s: f64,
    /// Original source path, discrete reaction and equivalent applied force [N].
    pub maximum_force_n: f64,
    pub maximum_iterations: usize,
    /// Original endpoint/corner queries, including rejected Newton probes and
    /// final whole-path admission. Original scalar roots retain their own cap.
    pub maximum_evaluations: usize,
    pub maximum_backtracks: usize,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub velocity_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum JkrMixedModeMotionError {
    Path(JkrMixedModePathError),
    Factor(LuWorkspaceError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    Budget,
    Range(&'static str),
    /// Refine the physical duration; no arbitrary alternative root is chosen.
    InertialTangent,
    NotConverged {
        residual_m: f64,
        tolerance_m: f64,
    },
    VelocityBalance {
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<JkrMixedModePathError> for JkrMixedModeMotionError {
    fn from(error: JkrMixedModePathError) -> Self {
        Self::Path(error)
    }
}
impl From<LuWorkspaceError> for JkrMixedModeMotionError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Factor(error)
    }
}
impl core::fmt::Display for JkrMixedModeMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled full-stick JKR motion: {self:?}")
    }
}
impl std::error::Error for JkrMixedModeMotionError {}

/// Original source, shear capacity, full mobility and numerical policy token.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeMotionState {
    identity: ContentHash,
    position_m: [f64; 3],
    velocity_m_s: [f64; 3],
    time_s: f64,
    steps: u64,
}
impl JkrMixedModeMotionState {
    pub const fn position_m(self) -> [f64; 3] {
        self.position_m
    }
    pub const fn velocity_m_s(self) -> [f64; 3] {
        self.velocity_m_s
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    pub const fn steps(self) -> u64 {
        self.steps
    }
}

/// One complete admitted phase. Contact work is reversible; there is no snap
/// release or friction heat on this attached full-stick branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeMotionStep {
    initial: JkrMixedModeMotionState,
    state: JkrMixedModeMotionState,
    pub initial_contact: JkrMixedModeResponse,
    pub contact: JkrMixedModeResponse,
    pub mean_force_n: [f64; 3],
    /// Positive reactions resist positive approach/shear coordinates.
    pub force_impulse_n_s: [f64; 3],
    /// Time is local to this constant-discrete-effort phase [N s²].
    pub force_time_moment_n_s2: [f64; 3],
    pub contact_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub applied_work_j: f64,
    pub energy_residual_j: f64,
    pub source_work_residual_j: f64,
    pub coordinate_residual_m: [f64; 3],
    pub velocity_residual_m_s: [f64; 3],
    pub maximum_source_force_n: f64,
    pub iterations: usize,
    pub evaluations: usize,
}
impl JkrMixedModeMotionStep {
    pub const fn initial_state(&self) -> JkrMixedModeMotionState {
        self.initial
    }
    pub const fn state(&self) -> JkrMixedModeMotionState {
        self.state
    }
}

/// Fixed source and complete held native mobility. Setup allocates one 3x3 LU;
/// repeated proposals reuse it and allocate no trajectory or Newton storage.
pub struct JkrMixedModeMotion<'a> {
    source: &'a SourceJkrMixedMode,
    shear: &'a SourceAdhesiveShear,
    mobility: [[f64; 3]; 3],
    mass: [[f64; 3]; 3],
    budget: JkrMixedModeMotionBudget,
    identity: ContentHash,
    lu: LuWorkspace,
}

struct Trial {
    position: [f64; 3],
    path: JkrMixedModePath,
    residual: [f64; 3],
    norm: f64,
    jacobian: [f64; 9],
}

impl<'a> JkrMixedModeMotion<'a> {
    pub fn new(
        source: &'a SourceJkrMixedMode,
        shear: &'a SourceAdhesiveShear,
        mobility: [[f64; 3]; 3],
        budget: JkrMixedModeMotionBudget,
    ) -> Result<Self, JkrMixedModeMotionError> {
        if shear.normal_state_identity() != source.normal_source().normal_state().identity() {
            return Err(JkrMixedModeMotionError::StateMismatch);
        }
        validate_budget(budget)?;
        let scale = mobility
            .iter()
            .flatten()
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        if mobility.iter().flatten().any(|v| !v.is_finite())
            || (0..3).any(|i| {
                (0..i)
                    .any(|j| (mobility[i][j] - mobility[j][i]).abs() > 64.0 * f64::EPSILON * scale)
            })
            || !positive_symmetric(mobility)
        {
            return Err(JkrMixedModeMotionError::Invalid(
                "finite reciprocal positive mobility required",
            ));
        }
        let mut lu = LuWorkspace::new(3)?;
        lu.factor(&core::array::from_fn::<_, 9, _>(|i| mobility[i / 3][i % 3]))?;
        let mut mass = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut unit = [0.0; 3];
            unit[j] = 1.0;
            let mut column = [0.0; 3];
            lu.solve_factored_into(&unit, &mut column)?;
            for i in 0..3 {
                mass[i][j] = column[i];
            }
        }
        if !positive_symmetric(mass) {
            return Err(JkrMixedModeMotionError::Invalid(
                "unrepresentable native inertia",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-mixed-mode-motion.v1");
        hash.update(source.identity().as_bytes());
        hash.update(shear.resolved().card_identity().as_bytes());
        for v in mobility.into_iter().flatten().chain([
            shear.law().static_shear_strength_pa(),
            shear.law().kinetic_shear_strength_pa(),
            budget.maximum_step_s,
            budget.maximum_coordinate_increment_m,
            budget.maximum_speed_m_s,
            budget.maximum_force_n,
            budget.coordinate_absolute_tolerance_m,
            budget.coordinate_relative_tolerance,
            budget.velocity_absolute_tolerance_m_s,
            budget.velocity_relative_tolerance,
            budget.energy_absolute_tolerance_j,
            budget.energy_relative_tolerance,
        ]) {
            hash.update(&v.to_bits().to_le_bytes());
        }
        for v in [
            budget.maximum_iterations,
            budget.maximum_evaluations,
            budget.maximum_backtracks,
        ] {
            hash.update(&(v as u64).to_le_bytes());
        }
        Ok(Self {
            source,
            shear,
            mobility,
            mass,
            budget,
            identity: hash.finalize(),
            lu,
        })
    }

    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.shear
    }
    pub const fn mobility_per_kg(&self) -> [[f64; 3]; 3] {
        self.mobility
    }
    pub const fn effective_mass_matrix_kg(&self) -> [[f64; 3]; 3] {
        self.mass
    }
    pub const fn budget(&self) -> JkrMixedModeMotionBudget {
        self.budget
    }

    /// Start on an explicitly supplied attached full-stick state. This does
    /// not nucleate a bond or invent the elastic shear displacement at contact.
    pub fn initialize(
        &self,
        position_m: [f64; 3],
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<JkrMixedModeMotionState, JkrMixedModeMotionError> {
        poll(gate)?;
        self.admit_velocity(velocity_m_s)?;
        let bounds = self.source.admit_quadratic_path(
            self.shear,
            [position_m; 3],
            self.budget.maximum_evaluations,
            gate,
        )?;
        if bounds.maximum_force_n > self.budget.maximum_force_n {
            return Err(JkrMixedModeMotionError::Range("initial source force"));
        }
        if bounds.evaluations >= self.budget.maximum_evaluations {
            return Err(JkrMixedModeMotionError::Budget);
        }
        // This also rejects a marginal or numerically unresolved radial secant.
        self.source.discrete_path(position_m, position_m, gate)?;
        finite(
            0.5 * dot(velocity_m_s, mul(self.mass, velocity_m_s)),
            "initial kinetic energy",
        )?;
        Ok(JkrMixedModeMotionState {
            identity: self.identity,
            position_m,
            velocity_m_s,
            time_s: 0.0,
            steps: 0,
        })
    }

    /// Advance one held-acceleration native phase on the same accepted history.
    /// The acceleration is J M^-1 f_external, supplied by the native owner.
    /// A source, capacity, solve, budget or cancellation refusal changes no
    /// accepted token; internal LU scratch may be reused on a deterministic retry.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &mut self,
        previous: JkrMixedModeMotionState,
        external_acceleration_m_s2: [f64; 3],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<JkrMixedModeMotionStep, JkrMixedModeMotionError> {
        poll(gate)?;
        if previous.identity != self.identity {
            return Err(JkrMixedModeMotionError::StateMismatch);
        }
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > self.budget.maximum_step_s {
            return Err(JkrMixedModeMotionError::Range("phase duration"));
        }
        self.admit_velocity(previous.velocity_m_s)?;
        let end_time = finite(previous.time_s + dt_s, "motion clock")?;
        if end_time <= previous.time_s {
            return Err(JkrMixedModeMotionError::Range(
                "unrepresentable clock increment",
            ));
        }
        let steps = previous
            .steps
            .checked_add(1)
            .ok_or(JkrMixedModeMotionError::Budget)?;
        let half_dt_squared = finite(0.5 * dt_s * dt_s, "squared duration")?;
        if half_dt_squared == 0.0 {
            return Err(JkrMixedModeMotionError::Range("squared duration underflow"));
        }
        let applied_force = mul(self.mass, external_acceleration_m_s2);
        if external_acceleration_m_s2.iter().any(|v| !v.is_finite())
            || !norm(applied_force).is_finite()
            || norm(applied_force) > self.budget.maximum_force_n
        {
            return Err(JkrMixedModeMotionError::Range("equivalent applied force"));
        }
        let free_increment = core::array::from_fn(|i| {
            dt_s * previous.velocity_m_s[i] + half_dt_squared * external_acceleration_m_s2[i]
        });
        let mut evaluations = 0;
        let mut trial = self.trial(
            previous.position_m,
            previous.position_m,
            None,
            free_increment,
            half_dt_squared,
            &mut evaluations,
            gate,
        )?;
        for iteration in 0..=self.budget.maximum_iterations {
            poll(gate)?;
            let tolerance = self.coordinate_tolerance(previous.position_m, trial.position);
            if trial.norm <= tolerance {
                return self.finish(
                    previous,
                    trial,
                    applied_force,
                    external_acceleration_m_s2,
                    dt_s,
                    end_time,
                    steps,
                    iteration,
                    evaluations,
                    gate,
                );
            }
            if iteration == self.budget.maximum_iterations {
                return Err(JkrMixedModeMotionError::NotConverged {
                    residual_m: trial.norm,
                    tolerance_m: tolerance,
                });
            }
            let mut correction = [0.0; 3];
            self.lu.solve_into(
                &trial.jacobian,
                &trial.residual.map(|v| -v),
                &mut correction,
            )?;
            let mut next = None;
            let mut fraction = 1.0;
            let mut last_error = None;
            for _ in 0..=self.budget.maximum_backtracks {
                poll(gate)?;
                let position =
                    core::array::from_fn(|i| trial.position[i] + fraction * correction[i]);
                match self.trial(
                    previous.position_m,
                    position,
                    Some(trial.path.before),
                    free_increment,
                    half_dt_squared,
                    &mut evaluations,
                    gate,
                ) {
                    Ok(candidate)
                        if candidate.norm < trial.norm
                            || candidate.norm
                                <= self.coordinate_tolerance(previous.position_m, position) =>
                    {
                        next = Some(candidate);
                        break;
                    }
                    Ok(_) => {}
                    Err(JkrMixedModeMotionError::Budget) => {
                        return Err(JkrMixedModeMotionError::Budget);
                    }
                    Err(error) => {
                        poll(gate)?;
                        last_error = Some(error);
                    }
                }
                fraction *= 0.5;
            }
            trial = match next {
                Some(next) => next,
                None => {
                    return Err(last_error.unwrap_or(JkrMixedModeMotionError::NotConverged {
                        residual_m: trial.norm,
                        tolerance_m: tolerance,
                    }));
                }
            };
        }
        unreachable!("bounded Newton loop always returns")
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        from: [f64; 3],
        to: [f64; 3],
        before: Option<JkrMixedModeResponse>,
        free: [f64; 3],
        half_dt_squared: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Trial, JkrMixedModeMotionError> {
        poll(gate)?;
        let increment = sub(to, from);
        if !norm(increment).is_finite()
            || norm(increment) > self.budget.maximum_coordinate_increment_m
        {
            return Err(JkrMixedModeMotionError::Range("coordinate increment"));
        }
        let allowance = usize::from(before.is_none()) + usize::from(from != to);
        if self.budget.maximum_evaluations.saturating_sub(*evaluations) < allowance {
            return Err(JkrMixedModeMotionError::Budget);
        }
        // Retain the immutable source-owned initial response across Newton
        // probes. Only a new endpoint needs its constitutive radius solve.
        *evaluations += allowance;
        let path = match before {
            Some(before) => self.source.discrete_path_from_response(before, to, gate)?,
            None => self.source.discrete_path(from, to, gate)?,
        };
        if norm(path.mean_force_n) > self.budget.maximum_force_n {
            return Err(JkrMixedModeMotionError::Range("discrete contact force"));
        }
        let acceleration = mul(self.mobility, path.mean_force_n);
        let residual =
            core::array::from_fn(|i| increment[i] - free[i] + half_dt_squared * acceleration[i]);
        let residual_norm = finite(norm(residual), "coordinate residual")?;
        let inertial = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                self.mass[i][j] + half_dt_squared * path.end_tangent_n_per_m[i][j]
            })
        });
        if !positive_symmetric(inertial) {
            return Err(JkrMixedModeMotionError::InertialTangent);
        }
        let jacobian = core::array::from_fn(|index| {
            let i = index / 3;
            let j = index % 3;
            let column = core::array::from_fn(|k| path.end_tangent_n_per_m[k][j]);
            (if i == j { 1.0 } else { 0.0 }) + half_dt_squared * dot(self.mobility[i], column)
        });
        Ok(Trial {
            position: to,
            path,
            residual,
            norm: residual_norm,
            jacobian,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        previous: JkrMixedModeMotionState,
        trial: Trial,
        applied_force: [f64; 3],
        external: [f64; 3],
        dt: f64,
        time_s: f64,
        steps: u64,
        iterations: usize,
        mut evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrMixedModeMotionStep, JkrMixedModeMotionError> {
        let reaction = mul(self.mobility, trial.path.mean_force_n);
        let velocity =
            core::array::from_fn(|i| previous.velocity_m_s[i] + dt * (external[i] - reaction[i]));
        self.admit_velocity(velocity)?;
        let midpoint = core::array::from_fn(|i| 0.5 * previous.velocity_m_s[i] + 0.5 * velocity[i]);
        let increment = sub(trial.position, previous.position_m);
        let coordinate_residual = core::array::from_fn(|i| increment[i] - dt * midpoint[i]);
        let coordinate_tolerance = self.coordinate_tolerance(previous.position_m, trial.position);
        if norm(coordinate_residual) > coordinate_tolerance {
            return Err(JkrMixedModeMotionError::NotConverged {
                residual_m: norm(coordinate_residual),
                tolerance_m: coordinate_tolerance,
            });
        }
        let velocity_increment = sub(velocity, previous.velocity_m_s);
        let velocity_residual =
            core::array::from_fn(|i| velocity_increment[i] - dt * (external[i] - reaction[i]));
        let velocity_tolerance = self.budget.velocity_absolute_tolerance_m_s
            + self.budget.velocity_relative_tolerance
                * norm(velocity).max(norm(previous.velocity_m_s));
        if norm(velocity_residual) > velocity_tolerance {
            return Err(JkrMixedModeMotionError::VelocityBalance {
                residual_m_s: norm(velocity_residual),
                tolerance_m_s: velocity_tolerance,
            });
        }
        let control =
            core::array::from_fn(|i| previous.position_m[i] + 0.5 * dt * previous.velocity_m_s[i]);
        if norm(sub(control, previous.position_m)) > self.budget.maximum_coordinate_increment_m {
            return Err(JkrMixedModeMotionError::Range(
                "quadratic position excursion",
            ));
        }
        let bounds = self.source.admit_quadratic_path(
            self.shear,
            [previous.position_m, control, trial.position],
            self.budget.maximum_evaluations.saturating_sub(evaluations),
            gate,
        )?;
        evaluations += bounds.evaluations;
        if bounds.maximum_force_n > self.budget.maximum_force_n {
            return Err(JkrMixedModeMotionError::Range("whole-path source force"));
        }
        let kinetic_change = finite(
            dot(midpoint, mul(self.mass, velocity_increment)),
            "kinetic change",
        )?;
        let applied_work = finite(dot(applied_force, increment), "applied work")?;
        let energy_residual = finite(
            (Dd::from_f64(kinetic_change) + Dd::from_f64(trial.path.reversible_work_j)
                - Dd::from_f64(applied_work))
            .to_f64(),
            "energy residual",
        )?;
        let energy_tolerance = finite(
            self.budget.energy_absolute_tolerance_j
                + self.budget.energy_relative_tolerance
                    * (kinetic_change.abs()
                        + trial.path.reversible_work_j.abs()
                        + applied_work.abs()),
            "energy tolerance",
        )?;
        if energy_residual.abs().max(trial.path.work_residual_j.abs()) > energy_tolerance {
            return Err(JkrMixedModeMotionError::EnergyBalance {
                residual_j: energy_residual.abs().max(trial.path.work_residual_j.abs()),
                tolerance_j: energy_tolerance,
            });
        }
        let impulse = trial.path.mean_force_n.map(|force| force * dt);
        let time_moment = impulse.map(|value| value * (0.5 * dt));
        if impulse
            .into_iter()
            .chain(time_moment)
            .any(|v| !v.is_finite())
        {
            return Err(JkrMixedModeMotionError::Range("integrated contact effort"));
        }
        poll(gate)?;
        Ok(JkrMixedModeMotionStep {
            initial: previous,
            state: JkrMixedModeMotionState {
                identity: self.identity,
                position_m: trial.position,
                velocity_m_s: velocity,
                time_s,
                steps,
            },
            initial_contact: trial.path.before,
            contact: trial.path.contact,
            mean_force_n: trial.path.mean_force_n,
            force_impulse_n_s: impulse,
            force_time_moment_n_s2: time_moment,
            contact_work_j: trial.path.reversible_work_j,
            kinetic_energy_change_j: kinetic_change,
            applied_work_j: applied_work,
            energy_residual_j: energy_residual,
            source_work_residual_j: trial.path.work_residual_j,
            coordinate_residual_m: coordinate_residual,
            velocity_residual_m_s: velocity_residual,
            maximum_source_force_n: bounds.maximum_force_n,
            iterations,
            evaluations,
        })
    }

    fn coordinate_tolerance(&self, from: [f64; 3], to: [f64; 3]) -> f64 {
        self.budget.coordinate_absolute_tolerance_m
            + self.budget.coordinate_relative_tolerance
                * norm(from)
                    .max(norm(to))
                    .max(self.source.normal_source().indentation_scale_m())
    }
    fn admit_velocity(&self, velocity: [f64; 3]) -> Result<(), JkrMixedModeMotionError> {
        if velocity.iter().any(|v| !v.is_finite()) || norm(velocity) > self.budget.maximum_speed_m_s
        {
            Err(JkrMixedModeMotionError::Range("normal/shear speed"))
        } else {
            Ok(())
        }
    }
}

fn validate_budget(b: JkrMixedModeMotionBudget) -> Result<(), JkrMixedModeMotionError> {
    if [
        b.maximum_step_s,
        b.maximum_coordinate_increment_m,
        b.maximum_speed_m_s,
        b.maximum_force_n,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || v <= 0.0)
        || [
            b.coordinate_absolute_tolerance_m,
            b.coordinate_relative_tolerance,
            b.velocity_absolute_tolerance_m_s,
            b.velocity_relative_tolerance,
            b.energy_absolute_tolerance_j,
            b.energy_relative_tolerance,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
        || b.maximum_iterations == 0
        || b.maximum_iterations > 256
        || b.maximum_evaluations < 4
        || b.maximum_evaluations > 1_048_576
        || b.maximum_backtracks > 64
    {
        Err(JkrMixedModeMotionError::Invalid(
            "finite physical bounds, nonnegative gates and bounded work required",
        ))
    } else {
        Ok(())
    }
}

// Three-dimensional physical reciprocity/coercivity admission. The actual
// system is solved by the existing prepared LU, never by this screening factor.
fn positive_symmetric(matrix: [[f64; 3]; 3]) -> bool {
    let scale = matrix.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs()));
    if !scale.is_finite() || scale <= 0.0 || matrix.iter().flatten().any(|v| !v.is_finite()) {
        return false;
    }
    let a: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| 0.5 * (matrix[i][j] / scale) + 0.5 * (matrix[j][i] / scale))
    });
    if a[0][0] <= 0.0 {
        return false;
    }
    let p1 = a[1][1] - (a[1][0] / a[0][0]) * a[0][1];
    let cross = a[2][1] - (a[2][0] / a[0][0]) * a[0][1];
    let p2 = a[2][2] - (a[2][0] / a[0][0]) * a[0][2] - (cross / p1) * cross;
    p1.is_finite() && p2.is_finite() && p1 > 0.0 && p2 > 0.0
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter()
        .zip(b)
        .fold(Dd::from_f64(0.0), |sum, (x, y)| {
            sum + Dd::from_f64(x) * Dd::from_f64(y)
        })
        .to_f64()
}
fn mul(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    a.map(|row| dot(row, x))
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn finite(value: f64, field: &'static str) -> Result<f64, JkrMixedModeMotionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrMixedModeMotionError::Range(field))
    }
}
fn poll(gate: &CancelGate) -> Result<(), JkrMixedModeMotionError> {
    if gate.is_requested() {
        Err(JkrMixedModeMotionError::Cancelled)
    } else {
        Ok(())
    }
}
