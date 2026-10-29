//! Joint DMT normal adhesion and two-dimensional sourced junction shear.
//!
//! Coordinates are [signed approach, shear 1, shear 2]. The original DMT
//! potential gives the exact-work mean normal reaction. The original junction
//! law uses the analytically integrated REAL contact area on the same straight
//! approach path; attraction across a gap contributes no area. Shear traction
//! opposes the straight tangential path and its work is nonnegative friction
//! loss, separate from normal storage. There is no inferred Coulomb term,
//! elastic shear store, JKR peeling, plasticity, heat partition or bond healing.
//!
//! All entries of the original held-force compliance are retained. Sticking
//! eliminates its tangential block; sliding uses an analytic Jacobian and the
//! existing prepared LU. A sufficient small-step monotonicity screen excludes
//! unresolvable negative normal curvature and excessive normal/shear coupling.
//! A moving junction admits a zero-slip path only inside its kinetic cone.
//! Dynamic callers must reconcile its successor with the actual endpoint
//! velocity before restoring static capacity. This is a discrete path law,
//! not a certificate for unresolved intrastep reversals or physical trajectory
//! accuracy. Mechanical timestep refinement remains necessary during snap.

use super::{DmtError, DmtResponse, SourceDmtSphereContact};
use crate::interface_binding::adhesive_shear::{AdhesiveShearError, SourceAdhesiveShear};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_tribo::{TriboError, adhesive::AdhesiveContactState};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtCoupledConfig {
    pub maximum_iterations: usize,
    /// Complete constitutive path probes, including bracket/line-search probes.
    pub maximum_evaluations: usize,
    pub maximum_backtracks: usize,
    pub maximum_step_s: f64,
    pub minimum_indentation_m: f64,
    pub maximum_coordinate_increment_m: f64,
    pub maximum_slip_speed_m_s: f64,
    pub maximum_force_n: f64,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DmtCoupledRegime {
    Released,
    Sticking,
    Sliding,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DmtCoupledError {
    Normal(DmtError),
    Shear(AdhesiveShearError),
    Factor(LuWorkspaceError),
    Invalid(&'static str),
    StateMismatch,
    MechanicalStability,
    Cancelled,
    Budget,
    Range(&'static str),
    NotConverged { residual_m: f64, tolerance_m: f64 },
    ForceMismatch { residual_n: f64, tolerance_n: f64 },
}
impl core::fmt::Display for DmtCoupledError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint sourced DMT/junction trial refused: {self:?}")
    }
}
impl std::error::Error for DmtCoupledError {}
impl From<DmtError> for DmtCoupledError {
    fn from(error: DmtError) -> Self {
        Self::Normal(error)
    }
}
impl From<AdhesiveShearError> for DmtCoupledError {
    fn from(error: AdhesiveShearError) -> Self {
        Self::Shear(error)
    }
}
impl From<TriboError> for DmtCoupledError {
    fn from(error: TriboError) -> Self {
        Self::Shear(error.into())
    }
}
impl From<LuWorkspaceError> for DmtCoupledError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Factor(error)
    }
}

/// Original-source checkpoint; numerical proposals never mutate it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtCoupledState {
    identity: ContentHash,
    position_m: [f64; 3],
    regime: DmtCoupledRegime,
}
impl DmtCoupledState {
    pub const fn position_m(self) -> [f64; 3] {
        self.position_m
    }
    pub const fn regime(self) -> DmtCoupledRegime {
        self.regime
    }
}

/// Actual endpoint and analytic straight-path geometry/work, not held old area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtAdhesivePath {
    pub contact: DmtResponse,
    pub normal_work_j: f64,
    pub mean_normal_force_n: f64,
    pub mean_real_area_m2: f64,
    pub mean_repulsive_force_n: f64,
    pub static_capacity_n: f64,
    pub kinetic_capacity_n: f64,
    normal_derivative_n_per_m: f64,
    kinetic_derivative_n_per_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtCoupledStep {
    state: DmtCoupledState,
    pub path: DmtAdhesivePath,
    /// Positive efforts resist positive approach/shear coordinates.
    pub mean_force_n: [f64; 3],
    /// Regime over this path; the successor may be Released after detachment.
    pub regime: DmtCoupledRegime,
    pub friction_dissipation_j: f64,
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    pub iterations: usize,
    pub evaluations: usize,
}
impl DmtCoupledStep {
    pub const fn state(&self) -> DmtCoupledState {
        self.state
    }
}

/// Complete source-bound constitutive/mechanical solve; setup allocates only LU.
pub struct DmtCoupledLoadLine<'a> {
    normal: &'a SourceDmtSphereContact,
    shear: &'a SourceAdhesiveShear,
    compliance: [[f64; 3]; 3],
    tangent_inverse: [[f64; 2]; 2],
    schur: f64,
    config: DmtCoupledConfig,
    identity: ContentHash,
    lu: LuWorkspace,
}

impl<'a> DmtCoupledLoadLine<'a> {
    pub fn new(
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        compliance: [[f64; 3]; 3],
        config: DmtCoupledConfig,
    ) -> Result<Self, DmtCoupledError> {
        if shear.normal_state_identity() != normal.normal_state().identity() {
            return Err(DmtCoupledError::StateMismatch);
        }
        validate_config(config, normal)?;
        if compliance.iter().flatten().any(|v| !v.is_finite()) {
            return Err(DmtCoupledError::Invalid(
                "finite complete compliance required",
            ));
        }
        let scale = compliance
            .iter()
            .flatten()
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        if !positive(scale)
            || (0..3).any(|i| !positive(compliance[i][i]))
            || (0..3).any(|i| {
                (0..i).any(|j| {
                    (compliance[i][j] - compliance[j][i]).abs() > 64.0 * f64::EPSILON * scale
                })
            })
        {
            return Err(DmtCoupledError::Invalid(
                "positive reciprocal mechanical compliance required",
            ));
        }
        let mut lu = LuWorkspace::new(3)?;
        let flat: [f64; 9] = core::array::from_fn(|i| compliance[i / 3][i % 3]);
        lu.factor(&flat)?;
        let mut inverse = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut unit = [0.0; 3];
            unit[j] = 1.0;
            let mut column = [0.0; 3];
            lu.solve_factored_into(&unit, &mut column)?;
            for i in 0..3 {
                inverse[i][j] = column[i];
            }
        }
        // Analyze the symmetric part only; never alter the original matrix.
        let t01 = 0.5 * (inverse[1][2] + inverse[2][1]);
        let maximum = 0.5
            * (inverse[1][1] + inverse[2][2] + (inverse[1][1] - inverse[2][2]).hypot(2.0 * t01));
        let minimum = (inverse[1][1] * inverse[2][2] - t01 * t01) / maximum;
        let cross =
            (0.5 * (inverse[0][1] + inverse[1][0])).hypot(0.5 * (inverse[0][2] + inverse[2][0]));
        // Nbar' >= -tail/2; Abar' <= pi R/2. The symmetric cross term
        // from kinetic tau*Abar is at most tau*pi*R/4 for every slip direction.
        let coupling =
            cross + 0.25 * shear.law().kinetic_shear_strength_pa() * PI * normal.radius_m;
        let normal_margin = inverse[0][0] - 0.5 * normal.tail_stiffness_n_per_m;
        if !positive(minimum)
            || !positive(normal_margin)
            || !positive(normal_margin * minimum - coupling * coupling)
        {
            return Err(DmtCoupledError::MechanicalStability);
        }
        let d = compliance[1][1] * compliance[2][2] - compliance[1][2] * compliance[2][1];
        if !positive(d) {
            return Err(DmtCoupledError::MechanicalStability);
        }
        let tangent_inverse = [
            [compliance[2][2] / d, -compliance[1][2] / d],
            [-compliance[2][1] / d, compliance[1][1] / d],
        ];
        let cross_response = mul2(tangent_inverse, [compliance[1][0], compliance[2][0]]);
        let schur = compliance[0][0]
            - compliance[0][1] * cross_response[0]
            - compliance[0][2] * cross_response[1];
        if !positive(schur) || !positive(1.0 - 0.5 * schur * normal.tail_stiffness_n_per_m) {
            return Err(DmtCoupledError::MechanicalStability);
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.joint-dmt-junction.v1");
        hash.update(normal.identity().as_bytes());
        hash.update(shear.resolved().identity().as_bytes());
        let identity = hash.finalize();
        Ok(Self {
            normal,
            shear,
            compliance,
            tangent_inverse,
            schur,
            config,
            identity,
            lu,
        })
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.normal
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.shear
    }
    pub const fn config(&self) -> DmtCoupledConfig {
        self.config
    }
    pub const fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.compliance
    }
    pub fn initialize(&self, position_m: [f64; 3]) -> Result<DmtCoupledState, DmtCoupledError> {
        if position_m.iter().any(|v| !v.is_finite())
            || position_m[0] < self.config.minimum_indentation_m
        {
            return Err(DmtCoupledError::Range("initial coordinates"));
        }
        let contact = self.normal.evaluate(position_m[0])?;
        Ok(DmtCoupledState {
            identity: self.identity,
            position_m,
            regime: if contact.contact_area_m2 == 0.0 {
                DmtCoupledRegime::Released
            } else {
                DmtCoupledRegime::Sticking
            },
        })
    }

    /// Initialize from a declared physical tangential velocity. A moving real
    /// junction starts on the kinetic branch, so the first step cannot capture
    /// it using an unevidenced static-strength jump. Zero velocity uses the
    /// ordinary fresh-contact initialization; gap-only attraction stays released.
    pub fn initialize_with_velocity(
        &self,
        position_m: [f64; 3],
        tangential_velocity_m_s: [f64; 2],
    ) -> Result<DmtCoupledState, DmtCoupledError> {
        self.with_endpoint_velocity(self.initialize(position_m)?, tangential_velocity_m_s)
    }

    /// Retain the authenticated position and select occupancy from an
    /// independently known physical tangential velocity. A zero mean slip
    /// increment does not prove that a dynamic endpoint is at rest: the motion
    /// may reverse inside the step. Such an endpoint keeps kinetic capacity.
    /// This does not resolve that reversal or add a constitutive evaluation.
    /// A detached endpoint stays Released, regardless of its velocity.
    pub fn with_endpoint_velocity(
        &self,
        mut state: DmtCoupledState,
        tangential_velocity_m_s: [f64; 2],
    ) -> Result<DmtCoupledState, DmtCoupledError> {
        if state.identity != self.identity {
            return Err(DmtCoupledError::StateMismatch);
        }
        let speed = tangential_velocity_m_s[0].hypot(tangential_velocity_m_s[1]);
        if !speed.is_finite() || speed > self.config.maximum_slip_speed_m_s {
            return Err(DmtCoupledError::Range("physical endpoint slip speed"));
        }
        if state.regime != DmtCoupledRegime::Released {
            state.regime = if speed > 0.0 {
                DmtCoupledRegime::Sliding
            } else {
                DmtCoupledRegime::Sticking
            };
        }
        Ok(state)
    }
    /// Analytic path moments, retaining the original normal potential and shear law.
    pub fn path(&self, from_m: f64, to_m: f64) -> Result<DmtAdhesivePath, DmtCoupledError> {
        let start = self.normal.evaluate(from_m)?;
        let contact = self.normal.evaluate(to_m)?;
        let work = self.normal.reversible_work_between(from_m, to_m)?;
        let dx = checked(to_m - from_m)?;
        let normal_force = if dx == 0.0 {
            contact.normal_force_n
        } else {
            checked(work / dx)?
        };
        let (mean_x, derivative_x) = mean_positive(from_m, to_m);
        let area = checked(PI * self.normal.radius_m * mean_x)?;
        let repulsive = hertz_mean(self.normal.hertz_coefficient, from_m, to_m);
        let contact_state = AdhesiveContactState::new(area, repulsive)?;
        let static_capacity_n = self.shear.law().static_capacity_n(contact_state)?;
        // SourceAdhesiveShear selects constant junction strength with exactly
        // zero Coulomb coefficients. Its kinetic limit is therefore tau_k*A,
        // including at stopping; no fictitious positive velocity is evaluated.
        let kinetic_capacity_n = self.shear.law().kinetic_shear_strength_pa() * area;
        let (_, shifted_derivative) = mean_positive(
            from_m + self.normal.interaction_range_m,
            to_m + self.normal.interaction_range_m,
        );
        let normal_derivative = hertz_derivative(self.normal.hertz_coefficient, from_m, to_m)
            - self.normal.tail_stiffness_n_per_m * (shifted_derivative - derivative_x);
        let kinetic_derivative =
            self.shear.law().kinetic_shear_strength_pa() * PI * self.normal.radius_m * derivative_x;
        for value in [
            start.reversible_energy_j,
            work,
            normal_force,
            area,
            repulsive,
            static_capacity_n,
            kinetic_capacity_n,
            normal_derivative,
            kinetic_derivative,
        ] {
            checked(value)?;
        }
        Ok(DmtAdhesivePath {
            contact,
            normal_work_j: work,
            mean_normal_force_n: normal_force,
            mean_real_area_m2: area,
            mean_repulsive_force_n: repulsive,
            static_capacity_n,
            kinetic_capacity_n,
            normal_derivative_n_per_m: normal_derivative,
            kinetic_derivative_n_per_m: kinetic_derivative,
        })
    }

    pub fn trial(
        &mut self,
        previous: DmtCoupledState,
        free_position_m: [f64; 3],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<DmtCoupledStep, DmtCoupledError> {
        self.admit_previous(previous, dt_s, gate)?;
        if free_position_m.iter().any(|v| !v.is_finite()) {
            return Err(DmtCoupledError::Invalid("finite free coordinates required"));
        }
        let old = previous.position_m;
        let free_shear = [free_position_m[1] - old[1], free_position_m[2] - old[2]];
        let reaction_free = mul2(self.tangent_inverse, free_shear);
        let effective_free = free_position_m[0]
            - self.compliance[0][1] * reaction_free[0]
            - self.compliance[0][2] * reaction_free[1];
        let mut evaluations = 0;
        let (x, path, normal_iterations, static_admissible) =
            self.static_normal(previous, effective_free, &mut evaluations, gate)?;
        let shear_rhs = [
            free_shear[0] - self.compliance[1][0] * path.mean_normal_force_n,
            free_shear[1] - self.compliance[2][0] * path.mean_normal_force_n,
        ];
        let required = mul2(self.tangent_inverse, shear_rhs);
        let capacity = if previous.regime == DmtCoupledRegime::Sliding {
            path.kinetic_capacity_n
        } else {
            path.static_capacity_n
        };
        let magnitude = required[0].hypot(required[1]);
        if static_admissible && magnitude <= capacity {
            let force = [path.mean_normal_force_n, required[0], required[1]];
            let regime = if path.mean_real_area_m2 == 0.0 {
                DmtCoupledRegime::Released
            } else {
                DmtCoupledRegime::Sticking
            };
            return self.finish(
                previous,
                [x, old[1], old[2]],
                free_position_m,
                force,
                regime,
                path,
                dt_s,
                normal_iterations,
                evaluations,
                gate,
            );
        }
        if !positive(magnitude) {
            return Err(DmtCoupledError::Range(
                "no admitted static equilibrium or sliding direction",
            ));
        }
        // A kinetic force smaller than the required stopping reaction produces
        // a nonzero slip seed, including after an occupied-branch continuation.
        let guess = [
            path.mean_normal_force_n,
            path.kinetic_capacity_n * required[0] / magnitude,
            path.kinetic_capacity_n * required[1] / magnitude,
        ];
        let response = mul3(self.compliance, guess);
        let candidate: [f64; 3] = core::array::from_fn(|i| free_position_m[i] - response[i]);
        let static_position = [x, old[1], old[2]];
        let mut alpha = 1.0;
        let mut initial = None;
        for _ in 0..self.config.maximum_backtracks {
            let q = core::array::from_fn(|i| {
                static_position[i] + alpha * (candidate[i] - static_position[i])
            });
            match self.sliding(previous, q, free_position_m, dt_s, &mut evaluations, gate) {
                Ok(value) => {
                    initial = Some(value);
                    break;
                }
                Err(error) if domain_error(&error) => alpha *= 0.5,
                Err(error) => return Err(error),
            }
        }
        let mut current = initial.ok_or(DmtCoupledError::Range("sliding seed"))?;
        for iteration in 0..=self.config.maximum_iterations {
            let tolerance =
                self.tolerance(previous, current.position, free_position_m, current.force)?;
            let residual = norm_inf(current.residual);
            if residual <= tolerance {
                let regime = if current.path.mean_real_area_m2 == 0.0 {
                    DmtCoupledRegime::Released
                } else {
                    DmtCoupledRegime::Sliding
                };
                return self.finish(
                    previous,
                    current.position,
                    free_position_m,
                    current.force,
                    regime,
                    current.path,
                    dt_s,
                    normal_iterations + iteration,
                    evaluations,
                    gate,
                );
            }
            if iteration == self.config.maximum_iterations {
                return Err(DmtCoupledError::NotConverged {
                    residual_m: residual,
                    tolerance_m: tolerance,
                });
            }
            let mut correction = [0.0; 3];
            self.lu.solve_into(
                &current.jacobian,
                &current.residual.map(|v| -v),
                &mut correction,
            )?;
            let mut alpha = 1.0;
            let mut next = None;
            for _ in 0..self.config.maximum_backtracks {
                let q = core::array::from_fn(|i| current.position[i] + alpha * correction[i]);
                match self.sliding(previous, q, free_position_m, dt_s, &mut evaluations, gate) {
                    Ok(value) if norm_inf(value.residual) < residual => {
                        next = Some(value);
                        break;
                    }
                    Ok(_) => {}
                    Err(error) if domain_error(&error) => {}
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            current = next.ok_or(DmtCoupledError::NotConverged {
                residual_m: residual,
                tolerance_m: tolerance,
            })?;
        }
        unreachable!("bounded sliding solve returns at its iteration ceiling")
    }

    /// Reevaluate the original law on ACTUAL body motion. The supplied effort is
    /// independently checked against that path, including the occupied cone.
    pub fn advance_actual(
        &self,
        previous: DmtCoupledState,
        position: [f64; 3],
        applied: [f64; 3],
        regime: DmtCoupledRegime,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<DmtCoupledStep, DmtCoupledError> {
        self.admit_previous(previous, dt_s, gate)?;
        self.admit_motion(previous, position, dt_s)?;
        if applied
            .iter()
            .any(|force| !force.is_finite() || force.abs() > self.config.maximum_force_n)
        {
            return Err(DmtCoupledError::Range("actual applied force"));
        }
        let path = self.path(previous.position_m[0], position[0])?;
        let delta = [
            position[1] - previous.position_m[1],
            position[2] - previous.position_m[2],
        ];
        let distance = delta[0].hypot(delta[1]);
        let mut expected = applied;
        expected[0] = path.mean_normal_force_n;
        match regime {
            DmtCoupledRegime::Released => {
                if path.mean_real_area_m2 != 0.0 {
                    return Err(DmtCoupledError::Invalid("released path has real contact"));
                }
                expected[1] = 0.0;
                expected[2] = 0.0;
            }
            DmtCoupledRegime::Sticking => {
                let tolerance = self.config.coordinate_absolute_tolerance_m
                    + self.config.coordinate_relative_tolerance
                        * (position[0] - previous.position_m[0]).abs();
                if distance > tolerance {
                    return Err(DmtCoupledError::Range("actual sticking displacement"));
                }
                let cap = if previous.regime == DmtCoupledRegime::Sliding {
                    path.kinetic_capacity_n
                } else {
                    path.static_capacity_n
                };
                let magnitude = applied[1].hypot(applied[2]);
                let tol = self.force_tolerance(magnitude, cap)?;
                if magnitude > cap + tol {
                    return Err(DmtCoupledError::ForceMismatch {
                        residual_n: magnitude - cap,
                        tolerance_n: tol,
                    });
                }
            }
            DmtCoupledRegime::Sliding => {
                if !positive(distance) {
                    return Err(DmtCoupledError::Range("zero actual sliding path"));
                }
                expected[1] = path.kinetic_capacity_n * (delta[0] / distance);
                expected[2] = path.kinetic_capacity_n * (delta[1] / distance);
            }
        }
        for i in 0..3 {
            let tolerance = self.force_tolerance(applied[i], expected[i])?;
            let residual = checked(applied[i] - expected[i])?;
            if residual.abs() > tolerance {
                return Err(DmtCoupledError::ForceMismatch {
                    residual_n: residual,
                    tolerance_n: tolerance,
                });
            }
        }
        // Reconstructed free coordinate is diagnostic only; actual native
        // compatibility and total body work remain the receiver's obligations.
        let response = mul3(self.compliance, expected);
        let free = core::array::from_fn(|i| position[i] + response[i]);
        self.finish(
            previous, position, free, expected, regime, path, dt_s, 0, 1, gate,
        )
    }

    fn admit_previous(
        &self,
        previous: DmtCoupledState,
        dt: f64,
        gate: &CancelGate,
    ) -> Result<(), DmtCoupledError> {
        poll(gate)?;
        if previous.identity != self.identity {
            return Err(DmtCoupledError::StateMismatch);
        }
        if !positive(dt) || dt > self.config.maximum_step_s {
            return Err(DmtCoupledError::Range("physical timestep"));
        }
        Ok(())
    }
    fn admit_motion(
        &self,
        previous: DmtCoupledState,
        q: [f64; 3],
        dt: f64,
    ) -> Result<(), DmtCoupledError> {
        if q.iter().any(|v| !v.is_finite())
            || q[0] < self.config.minimum_indentation_m
            || q[0] > self.normal.envelope.maximum_indentation_m
            || (0..3).any(|i| {
                (q[i] - previous.position_m[i]).abs() > self.config.maximum_coordinate_increment_m
            })
        {
            return Err(DmtCoupledError::Range("physical coordinate path"));
        }
        if (q[1] - previous.position_m[1]).hypot(q[2] - previous.position_m[2]) / dt
            > self.config.maximum_slip_speed_m_s
        {
            return Err(DmtCoupledError::Range("physical slip speed"));
        }
        Ok(())
    }
    fn probe(
        &self,
        old: f64,
        x: f64,
        count: &mut usize,
        gate: &CancelGate,
    ) -> Result<DmtAdhesivePath, DmtCoupledError> {
        poll(gate)?;
        if *count >= self.config.maximum_evaluations {
            return Err(DmtCoupledError::Budget);
        }
        *count += 1;
        self.path(old, x)
    }
    fn static_normal(
        &self,
        previous: DmtCoupledState,
        free: f64,
        count: &mut usize,
        gate: &CancelGate,
    ) -> Result<(f64, DmtAdhesivePath, usize, bool), DmtCoupledError> {
        let old = previous.position_m[0];
        let mut lo = (old - self.config.maximum_coordinate_increment_m)
            .max(self.config.minimum_indentation_m);
        let mut hi = (old + self.config.maximum_coordinate_increment_m)
            .min(self.normal.envelope.maximum_indentation_m);
        let a = self.probe(old, lo, count, gate)?;
        let b = self.probe(old, hi, count, gate)?;
        // The all-stick constraint can place its root outside the source
        // envelope even when a coupled SLIDING root is inside. The appropriate
        // source boundary is then only a seed, never an accepted static state.
        // Every sliding iterate and the final actual path keep all range gates.
        if lo - free + self.schur * a.mean_normal_force_n > 0.0 {
            return Ok((lo, a, 0, false));
        }
        if hi - free + self.schur * b.mean_normal_force_n < 0.0 {
            return Ok((hi, b, 0, false));
        }
        let mut x = old;
        for iteration in 0..self.config.maximum_iterations {
            let path = self.probe(old, x, count, gate)?;
            let residual = checked(x - free + self.schur * path.mean_normal_force_n)?;
            let tolerance = self.config.coordinate_absolute_tolerance_m
                + self.config.coordinate_relative_tolerance
                    * (x - old)
                        .abs()
                        .max((free - old).abs())
                        .max((self.schur * path.mean_normal_force_n).abs());
            if residual.abs() <= tolerance {
                return Ok((x, path, iteration, true));
            }
            if residual > 0.0 {
                hi = x;
            } else {
                lo = x;
            }
            let derivative = checked(1.0 + self.schur * path.normal_derivative_n_per_m)?;
            if !positive(derivative) {
                return Err(DmtCoupledError::MechanicalStability);
            }
            let candidate = x - residual / derivative;
            let next = if candidate > lo && candidate < hi {
                candidate
            } else {
                f64::midpoint(lo, hi)
            };
            if next == x {
                return Err(DmtCoupledError::NotConverged {
                    residual_m: residual,
                    tolerance_m: tolerance,
                });
            }
            x = next;
        }
        Err(DmtCoupledError::Budget)
    }
    fn sliding(
        &self,
        previous: DmtCoupledState,
        position: [f64; 3],
        free: [f64; 3],
        dt: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Sliding, DmtCoupledError> {
        self.admit_motion(previous, position, dt)?;
        let path = self.probe(previous.position_m[0], position[0], evaluations, gate)?;
        let slip = [
            position[1] - previous.position_m[1],
            position[2] - previous.position_m[2],
        ];
        let distance = slip[0].hypot(slip[1]);
        if !positive(distance) {
            return Err(DmtCoupledError::Range("zero sliding iterate"));
        }
        let direction = slip.map(|v| v / distance);
        let force = [
            path.mean_normal_force_n,
            path.kinetic_capacity_n * direction[0],
            path.kinetic_capacity_n * direction[1],
        ];
        let response = mul3(self.compliance, force);
        let residual = core::array::from_fn(|i| position[i] - free[i] + response[i]);
        let mut derivative = [[0.0; 3]; 3];
        derivative[0][0] = path.normal_derivative_n_per_m;
        for i in 0..2 {
            derivative[i + 1][0] = path.kinetic_derivative_n_per_m * direction[i];
            for j in 0..2 {
                derivative[i + 1][j + 1] = path.kinetic_capacity_n / distance
                    * (f64::from(i == j) - direction[i] * direction[j]);
            }
        }
        let mut jacobian = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                jacobian[3 * i + j] = f64::from(i == j)
                    + (0..3)
                        .map(|k| self.compliance[i][k] * derivative[k][j])
                        .sum::<f64>();
            }
        }
        for value in residual.into_iter().chain(force).chain(jacobian) {
            checked(value)?;
        }
        Ok(Sliding {
            position,
            path,
            force,
            residual,
            jacobian,
        })
    }
    fn tolerance(
        &self,
        previous: DmtCoupledState,
        q: [f64; 3],
        free: [f64; 3],
        force: [f64; 3],
    ) -> Result<f64, DmtCoupledError> {
        let response = mul3(self.compliance, force);
        let scale = (0..3)
            .map(|i| {
                (q[i] - previous.position_m[i])
                    .abs()
                    .max((free[i] - previous.position_m[i]).abs())
                    .max(response[i].abs())
            })
            .fold(0.0_f64, f64::max);
        checked(
            self.config.coordinate_absolute_tolerance_m
                + self.config.coordinate_relative_tolerance * scale,
        )
    }
    fn force_tolerance(&self, a: f64, b: f64) -> Result<f64, DmtCoupledError> {
        checked(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance * a.abs().max(b.abs()),
        )
    }
    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        previous: DmtCoupledState,
        position: [f64; 3],
        free: [f64; 3],
        force: [f64; 3],
        regime: DmtCoupledRegime,
        path: DmtAdhesivePath,
        dt: f64,
        iterations: usize,
        evaluations: usize,
        gate: &CancelGate,
    ) -> Result<DmtCoupledStep, DmtCoupledError> {
        self.admit_motion(previous, position, dt)?;
        for f in force {
            if !f.is_finite() || f.abs() > self.config.maximum_force_n {
                return Err(DmtCoupledError::Range("physical force"));
            }
        }
        let response = mul3(self.compliance, force);
        let residual = core::array::from_fn(|i| position[i] - free[i] + response[i]);
        let tolerance = self.tolerance(previous, position, free, force)?;
        if norm_inf(residual) > tolerance {
            return Err(DmtCoupledError::NotConverged {
                residual_m: norm_inf(residual),
                tolerance_m: tolerance,
            });
        }
        let slip =
            (position[1] - previous.position_m[1]).hypot(position[2] - previous.position_m[2]);
        let loss = if regime == DmtCoupledRegime::Sliding {
            checked(path.kinetic_capacity_n * slip)?
        } else {
            0.0
        };
        let state = DmtCoupledState {
            identity: self.identity,
            position_m: position,
            regime: if path.contact.contact_area_m2 == 0.0 {
                DmtCoupledRegime::Released
            } else {
                regime
            },
        };
        poll(gate)?;
        Ok(DmtCoupledStep {
            state,
            path,
            mean_force_n: force,
            regime,
            friction_dissipation_j: loss,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            iterations,
            evaluations,
        })
    }
}

struct Sliding {
    position: [f64; 3],
    path: DmtAdhesivePath,
    force: [f64; 3],
    residual: [f64; 3],
    jacobian: [f64; 9],
}

fn mean_positive(a: f64, b: f64) -> (f64, f64) {
    if a > 0.0 && b > 0.0 {
        return (0.5 * a + 0.5 * b, 0.5);
    }
    if a <= 0.0 && b <= 0.0 {
        return (0.0, if a == 0.0 && b == 0.0 { 0.25 } else { 0.0 });
    }
    let width = (b - a).abs();
    if b > 0.0 {
        let r = b / width;
        (0.5 * b * r, r * (1.0 - 0.5 * r))
    } else {
        let r = a / width;
        (0.5 * a * r, 0.5 * r * r)
    }
}
fn hertz_mean(k: f64, a: f64, b: f64) -> f64 {
    let hi = a.max(b).max(0.0);
    let lo = a.min(b).max(0.0);
    if hi == 0.0 {
        return 0.0;
    }
    if a == b {
        return k * hi * hi.sqrt();
    }
    let r = (lo / hi).sqrt();
    let polynomial = 1.0 + r * (1.0 + r * (1.0 + r * (1.0 + r)));
    k * hi * hi.sqrt() * (0.4 * polynomial / (1.0 + r)) * ((hi - lo) / (b - a).abs())
}
fn hertz_derivative(k: f64, a: f64, b: f64) -> f64 {
    if a >= 0.0 && b >= 0.0 {
        let x = a.sqrt();
        let y = b.sqrt();
        let sum = x + y;
        if sum == 0.0 {
            return 0.0;
        }
        return 0.2 * k * (3.0 * y * y * y + 6.0 * y * y * x + 4.0 * y * x * x + 2.0 * x * x * x)
            / (sum * sum);
    }
    if a <= 0.0 && b <= 0.0 {
        return 0.0;
    }
    if b > 0.0 {
        let r = b / (b - a);
        k * b.sqrt() * r * (1.0 - 0.4 * r)
    } else {
        let r = a / (a - b);
        0.4 * k * a.sqrt() * r * r
    }
}
fn validate_config(
    c: DmtCoupledConfig,
    normal: &SourceDmtSphereContact,
) -> Result<(), DmtCoupledError> {
    if !(1..=128).contains(&c.maximum_iterations)
        || c.maximum_evaluations == 0
        || !(1..=32).contains(&c.maximum_backtracks)
        || !c.minimum_indentation_m.is_finite()
        || c.minimum_indentation_m >= normal.envelope.maximum_indentation_m
        || [
            c.maximum_step_s,
            c.maximum_coordinate_increment_m,
            c.maximum_slip_speed_m_s,
            c.maximum_force_n,
            c.coordinate_absolute_tolerance_m,
            c.force_absolute_tolerance_n,
        ]
        .into_iter()
        .any(|v| !positive(v))
        || [c.coordinate_relative_tolerance, c.force_relative_tolerance]
            .into_iter()
            .any(|v| !positive(v) || v >= 1.0)
    {
        return Err(DmtCoupledError::Invalid(
            "explicit positive physical and numerical budgets required",
        ));
    }
    Ok(())
}
fn mul2(a: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [
        a[0][0] * x[0] + a[0][1] * x[1],
        a[1][0] * x[0] + a[1][1] * x[1],
    ]
}
fn mul3(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| a[i][j] * x[j]).sum())
}
fn norm_inf(a: [f64; 3]) -> f64 {
    a.into_iter().fold(0.0, |a, b| a.max(b.abs()))
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn checked(x: f64) -> Result<f64, DmtCoupledError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DmtCoupledError::Invalid(
            "nonfinite constitutive arithmetic",
        ))
    }
}
fn poll(gate: &CancelGate) -> Result<(), DmtCoupledError> {
    if gate.is_requested() {
        Err(DmtCoupledError::Cancelled)
    } else {
        Ok(())
    }
}
fn domain_error(e: &DmtCoupledError) -> bool {
    matches!(
        e,
        DmtCoupledError::Range(_) | DmtCoupledError::Normal(DmtError::OutsideApplicability { .. })
    )
}
