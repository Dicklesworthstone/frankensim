//! Driven shear/spin under the original joint Hertz static/kinetic surface.
//!
//! In radius-scaled coordinates, `q'=a+W*w(q)` uses the complete held SPD
//! mobility. Rest is an explicit retained mode: the original static surface
//! either holds the complete drive or releases a constant-acceleration ray.
//! Moving invariant rays have exact first-zero events, followed by static
//! hold or overdriven restart. The larger static set never captures a moving
//! state merely because its speed is small.
//!
//! Curved paths reuse the original kinetic endpoint projection and hold each
//! accepted wrench over a linear velocity path. This is a bounded FIRST-ORDER
//! time discretization; its curved arrests are timestep-resolved events. A
//! capture requires strictly interior original kinetic-drive feasibility.
//! Both endpoint powers must be nonpositive. Heat is the held traction paired
//! with actual midpoint travel, not the projection's endpoint support work.
//! The original full-matrix balance is `Delta E + heat = (W^-1*a) dot travel`.
//! Geometry, load, source coefficients and drive remain fixed during a phase.

use super::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchPhase, Matrix, Vector, dot, factor_circular_response, multiply,
};
use fs_math::det;

/// Physical mode supplied by the preceding accepted phase, never a speed test.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HertzDrivenMode {
    Resting,
    Sliding,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzDrivenInput {
    pub velocity_m_s: [f64; 2],
    pub spin_rate_rad_s: f64,
    /// Original body mobility in `[u_x,u_y,a*omega]` coordinates [1/kg].
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    /// Relative acceleration from the original external body wrenches alone.
    pub driving_acceleration_m_s2: [f64; 2],
    pub driving_spin_acceleration_rad_s2: f64,
    pub dt_s: f64,
    pub mode: HertzDrivenMode,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzDrivenControls {
    /// Temporal refinement axis for curved held-traction intervals [s].
    pub maximum_step_s: f64,
    /// Refuse an unresolved curved direction/event below this step [s].
    pub minimum_step_s: f64,
    /// Original kinetic probes and projection attempts, including retries.
    pub maximum_attempts: usize,
    /// Total original Hertz support evaluations across all phase solves.
    pub maximum_evaluations: usize,
    /// Full original-matrix endpoint compatibility in scaled rate units [m/s].
    pub scaled_rate_tolerance_m_s: f64,
    /// Validates ONLY an explicitly retained resting constraint. It cannot
    /// capture a sliding input, and no actual body momentum is reset [m/s].
    pub rest_scaled_rate_tolerance_m_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct HertzDrivenReport {
    pub mode: HertzDrivenMode,
    pub velocity_after_m_s: [f64; 2],
    pub spin_rate_after_rad_s: f64,
    /// Original generalized external demand `W^-1*a_scaled`, in `[F,F,M/a]` [N].
    pub driving_scaled_wrench_n: [f64; 3],
    pub force_impulse_ns: [f64; 2],
    pub spin_angular_impulse_nms: f64,
    pub slip_increment_m: [f64; 2],
    pub spin_increment_rad: f64,
    /// Nonnegative physical held-traction work along the represented path [J].
    pub dissipated_energy_j: f64,
    /// Signed original generalized external work along the same path [J].
    pub drive_work_j: f64,
    pub kinetic_energy_change_j: f64,
    /// `kinetic_energy_change + dissipated_energy - drive_work` [J].
    pub energy_residual_j: f64,
    /// Original input plus external/contact impulses minus returned target.
    pub scaled_velocity_residual_m_s: [f64; 3],
    /// First moving-to-rest transition in this phase, if any [s].
    pub first_arrest_time_s: Option<f64>,
    /// Meaningful only when first_arrest_time_s exists. Curved endpoint
    /// capture is explicitly timestep-resolved, not an exact continuum event.
    pub arrest_time_is_exact: bool,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub evaluations: usize,
    pub projection_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HertzDrivenError {
    Contact(HertzWrenchError),
    Invalid(&'static str),
    AttemptBudget,
    EvaluationBudget,
    UnresolvedDirection,
    UnresolvedWork {
        scaled_velocity_residual_m_s: [f64; 3],
        energy_residual_j: f64,
    },
    Cancelled,
}
impl From<HertzWrenchError> for HertzDrivenError {
    fn from(error: HertzWrenchError) -> Self {
        if error == HertzWrenchError::Cancelled {
            Self::Cancelled
        } else {
            Self::Contact(error)
        }
    }
}
impl core::fmt::Display for HertzDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for HertzDrivenError {}
type Error = HertzDrivenError;

impl HertzWrenchContact {
    /// Advance one held-drive/load/geometry phase with explicit contact memory.
    ///
    /// Static hold and breakaway use the original joint surface; moving phases
    /// use only the original kinetic surface. General curved paths require
    /// refinement of maximum_step_s. All candidate state is private until the
    /// complete endpoint, work, bounded-effort and cancellation gates succeed.
    #[allow(clippy::too_many_lines)]
    pub fn advance_driven(
        &self,
        input: &HertzDrivenInput,
        config: HertzWrenchConfig,
        controls: HertzDrivenControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzDrivenReport, Error> {
        poll(&mut cancelled)?;
        config.validate()?;
        controls.validate()?;
        if !input.dt_s.is_finite() || input.dt_s <= 0.0 {
            return Err(Error::Invalid("positive finite driven duration required"));
        }
        let radius = self.parameters().contact_radius_m;
        let original = vector([
            input.velocity_m_s[0],
            input.velocity_m_s[1],
            radius * input.spin_rate_rad_s,
        ])?;
        let drive = vector([
            input.driving_acceleration_m_s2[0],
            input.driving_acceleration_m_s2[1],
            radius * input.driving_spin_acceleration_rad_s2,
        ])?;
        if !input.spin_rate_rad_s.is_finite()
            || !input.driving_spin_acceleration_rad_s2.is_finite()
            || (input.spin_rate_rad_s != 0.0 && original[2] == 0.0)
            || (input.driving_spin_acceleration_rad_s2 != 0.0 && drive[2] == 0.0)
        {
            return Err(Error::Invalid("unrepresentable scaled spin port"));
        }
        let matrix = input.scaled_mobility_per_kg;
        if (0..3).any(|i| (0..i).any(|j| matrix[i][j] != matrix[j][i])) {
            return Err(Error::Invalid("physical mobility must be symmetric"));
        }
        let factor = factor_circular_response(matrix, config)?;
        let demand = vector(
            factor
                .solve(drive)
                .ok_or(Error::Invalid("unrepresentable driving wrench"))?,
        )?;
        let mut q = match input.mode {
            HertzDrivenMode::Resting => {
                if norm(original) > controls.rest_scaled_rate_tolerance_m_s {
                    return Err(Error::Invalid(
                        "retained rest exceeds original rate residual gate",
                    ));
                }
                [0.0; 3]
            }
            HertzDrivenMode::Sliding => {
                if original == [0.0; 3] {
                    return Err(Error::Invalid("zero rate requires explicit resting mode"));
                }
                original
            }
        };
        let mut mode = input.mode;
        let mut work = Work {
            law: self,
            matrix,
            radius,
            config,
            controls,
            cancelled: &mut cancelled,
            attempts: 0,
            evaluations: 0,
            projection_iterations: 0,
        };
        let (mut impulse, mut travel) = ([0.0; 3], [0.0; 3]);
        let (mut heat, mut elapsed) = (0.0, 0.0);
        let (mut accepted, mut rejected) = (0, 0);
        let (mut first_arrest, mut arrest_exact) = (None, true);
        while elapsed < input.dt_s {
            poll(work.cancelled)?;
            let remaining = input.dt_s - elapsed;
            let (h, w, after, captured, exact) = if mode == HertzDrivenMode::Resting {
                let candidate =
                    work.project(q, drive, remaining, HertzWrenchPhase::StaticAdmissible)?;
                (
                    remaining,
                    candidate.wrench,
                    candidate.endpoint,
                    candidate.captured,
                    true,
                )
            } else {
                let w = work.kinetic(q)?;
                let acceleration = vector(add(drive, multiply(matrix, w)))?;
                if let Some((h, after, captured)) = ray(
                    q,
                    acceleration,
                    remaining,
                    controls.scaled_rate_tolerance_m_s,
                )? {
                    (h, w, after, captured, true)
                } else if w == [0.0; 3] {
                    // Genuine frictionless free motion is linear even when
                    // its direction changes; it needs no kinetic projection.
                    (
                        remaining,
                        w,
                        vector(add(q, scale(acceleration, remaining)))?,
                        false,
                        true,
                    )
                } else {
                    let mut h = remaining.min(controls.maximum_step_s);
                    loop {
                        // A representable final interval may be shorter than
                        // the refinement floor, including a last ulp left by
                        // accumulated clock arithmetic. Only its first trial
                        // is admitted; failed refinement still refuses.
                        if (h < controls.minimum_step_s && h != remaining) || elapsed + h == elapsed
                        {
                            return Err(Error::UnresolvedDirection);
                        }
                        let candidate = work.project(q, drive, h, HertzWrenchPhase::Kinetic)?;
                        let mut admissible = dot(candidate.wrench, q) <= 0.0
                            && dot(candidate.wrench, candidate.endpoint) <= 0.0;
                        if candidate.captured {
                            // Demand is enlarged, never the physical kinetic
                            // set. This conservative interior probe prevents
                            // coarse capture in the static/kinetic hysteresis
                            // band; it reuses original surface membership.
                            let strict = vector(drive.map(|x| x * (1.0 + 1024.0 * f64::EPSILON)))?;
                            let support =
                                work.project([0.0; 3], strict, h, HertzWrenchPhase::Kinetic)?;
                            admissible &= support.captured;
                        }
                        if admissible {
                            break (
                                h,
                                candidate.wrench,
                                candidate.endpoint,
                                candidate.captured,
                                false,
                            );
                        }
                        rejected += 1;
                        h *= 0.5;
                    }
                }
            };
            if !h.is_finite() || h <= 0.0 || elapsed + h == elapsed {
                return Err(Error::UnresolvedDirection);
            }
            // A rounded zero before a resolved event is not a resting state.
            // In particular q+a*h can cancel exactly one ulp before an exact
            // ray's computed stop; publishing Sliding/zero would corrupt the
            // explicit mode checkpoint on the next phase.
            if !captured && after == [0.0; 3] {
                return Err(Error::UnresolvedDirection);
            }
            let increment = vector(core::array::from_fn(|i| h * (0.5 * q[i] + 0.5 * after[i])))?;
            let loss = checked(-dot(w, increment))?;
            if loss < 0.0 || dot(w, q) > 0.0 || dot(w, after) > 0.0 {
                return Err(Error::UnresolvedDirection);
            }
            impulse = vector(add(impulse, scale(w, h)))?;
            travel = vector(add(travel, increment))?;
            heat = checked(heat + loss)?;
            elapsed = if h == remaining {
                input.dt_s
            } else {
                checked(elapsed + h)?
            };
            accepted += 1;
            if captured && mode == HertzDrivenMode::Sliding && first_arrest.is_none() {
                first_arrest = Some(elapsed);
                arrest_exact = exact;
            }
            mode = if captured {
                HertzDrivenMode::Resting
            } else {
                HertzDrivenMode::Sliding
            };
            q = after;
        }
        let recovered = vector(add(
            add(original, scale(drive, input.dt_s)),
            multiply(matrix, impulse),
        ))?;
        let residual = vector(subtract(recovered, q))?;
        let delta = factor
            .solve(subtract(recovered, original))
            .ok_or(Error::Invalid("unrepresentable kinetic change"))?;
        let kinetic = checked(dot(
            delta,
            core::array::from_fn(|i| 0.5 * original[i] + 0.5 * recovered[i]),
        ))?;
        let drive_work = checked(dot(demand, travel))?;
        let energy_residual = checked(kinetic + heat - drive_work)?;
        if norm(residual) > controls.scaled_rate_tolerance_m_s
            || energy_residual.abs() > config.energy_tolerance_j
        {
            return Err(Error::UnresolvedWork {
                scaled_velocity_residual_m_s: residual,
                energy_residual_j: energy_residual,
            });
        }
        poll(work.cancelled)?;
        Ok(HertzDrivenReport {
            mode,
            velocity_after_m_s: [q[0], q[1]],
            spin_rate_after_rad_s: checked(q[2] / radius)?,
            driving_scaled_wrench_n: demand,
            force_impulse_ns: [impulse[0], impulse[1]],
            spin_angular_impulse_nms: checked(radius * impulse[2])?,
            slip_increment_m: [travel[0], travel[1]],
            spin_increment_rad: checked(travel[2] / radius)?,
            dissipated_energy_j: heat,
            drive_work_j: drive_work,
            kinetic_energy_change_j: kinetic,
            energy_residual_j: energy_residual,
            scaled_velocity_residual_m_s: residual,
            first_arrest_time_s: first_arrest,
            arrest_time_is_exact: arrest_exact,
            accepted_steps: accepted,
            rejected_steps: rejected,
            evaluations: work.evaluations,
            projection_iterations: work.projection_iterations,
        })
    }
}

struct Candidate {
    wrench: Vector,
    endpoint: Vector,
    captured: bool,
}
struct Work<'a, C> {
    law: &'a HertzWrenchContact,
    matrix: Matrix,
    radius: f64,
    config: HertzWrenchConfig,
    controls: HertzDrivenControls,
    cancelled: &'a mut C,
    attempts: usize,
    evaluations: usize,
    projection_iterations: usize,
}
impl<C: FnMut() -> bool> Work<'_, C> {
    fn attempt(&mut self) -> Result<(), Error> {
        poll(self.cancelled)?;
        if self.attempts >= self.controls.maximum_attempts {
            return Err(Error::AttemptBudget);
        }
        if self.evaluations >= self.controls.maximum_evaluations {
            return Err(Error::EvaluationBudget);
        }
        self.attempts += 1;
        Ok(())
    }
    fn kinetic(&mut self, q: Vector) -> Result<Vector, Error> {
        self.attempt()?;
        self.evaluations += 1;
        let r = self
            .law
            .kinetic_patch
            .evaluate([q[0], q[1]], checked(q[2] / self.radius)?)
            .map_err(HertzWrenchError::Law)?;
        if det::hypot(r.force_n[0], r.force_n[1]) > self.config.maximum_force_n
            || r.spin_moment_nm.abs() > self.config.maximum_spin_moment_nm
        {
            return Err(HertzWrenchError::EffortBudget {
                force_n: det::hypot(r.force_n[0], r.force_n[1]),
                spin_moment_nm: r.spin_moment_nm.abs(),
            }
            .into());
        }
        poll(self.cancelled)?;
        vector([r.force_n[0], r.force_n[1], r.spin_moment_nm / self.radius])
    }
    fn project(
        &mut self,
        q: Vector,
        drive: Vector,
        h: f64,
        phase: HertzWrenchPhase,
    ) -> Result<Candidate, Error> {
        self.attempt()?;
        let square = checked(h * h)?;
        if square == 0.0 {
            return Err(Error::UnresolvedDirection);
        }
        let free = vector(core::array::from_fn(|i| h * (q[i] + h * drive[i])))?;
        let mut config = self.config;
        config.maximum_evaluations = config
            .maximum_evaluations
            .min(self.controls.maximum_evaluations - self.evaluations);
        config.displacement_tolerance_m = config
            .displacement_tolerance_m
            .min(h * self.controls.scaled_rate_tolerance_m_s);
        let r = self.law.advance_phase(
            &HertzWrenchInput {
                free_displacement_m: [free[0], free[1]],
                free_spin_angle_rad: checked(free[2] / self.radius)?,
                scaled_compliance_m_per_n: self.matrix.map(|row| row.map(|x| square * x)),
            },
            phase,
            config,
            &mut *self.cancelled,
        )?;
        self.evaluations = self
            .evaluations
            .checked_add(r.evaluations)
            .ok_or(Error::EvaluationBudget)?;
        self.projection_iterations = self
            .projection_iterations
            .checked_add(r.iterations)
            .ok_or(Error::AttemptBudget)?;
        let endpoint = vector([
            r.slip_displacement_m[0] / h,
            r.slip_displacement_m[1] / h,
            self.radius * (r.spin_angle_rad / h),
        ])?;
        poll(self.cancelled)?;
        Ok(Candidate {
            wrench: vector([r.force_n[0], r.force_n[1], r.spin_moment_nm / self.radius])?,
            endpoint,
            captured: r.mode == HertzWrenchMode::Sticking || endpoint == [0.0; 3],
        })
    }
}

// Original constant kinetic wrench gives an exact ray only if the COMPLETE
// driving acceleration plus mobility response preserves that direction.
fn ray(
    q: Vector,
    acceleration: Vector,
    remaining: f64,
    tolerance: f64,
) -> Result<Option<(f64, Vector, bool)>, Error> {
    let speed = norm(q);
    let magnitude = norm(acceleration);
    if magnitude == 0.0 {
        return Ok(Some((remaining, q, false)));
    }
    let direction = q.map(|v| v / speed);
    let a = acceleration.map(|v| v / magnitude);
    let parallel = dot(direction, a);
    let mismatch = norm(subtract(a, scale(direction, parallel)));
    if mismatch > 256.0 * f64::EPSILON {
        return Ok(None);
    }
    let along = checked(magnitude * parallel)?;
    let stop = if along < 0.0 {
        speed / -along
    } else {
        f64::INFINITY
    };
    let h = remaining.min(stop);
    if magnitude * mismatch * h > tolerance {
        return Ok(None);
    }
    let captured = h == stop;
    let after = if captured {
        [0.0; 3]
    } else {
        vector(add(q, scale(acceleration, h)))?
    };
    Ok(Some((h, after, captured)))
}
impl HertzDrivenControls {
    fn validate(self) -> Result<(), Error> {
        if !self.maximum_step_s.is_finite()
            || !self.minimum_step_s.is_finite()
            || self.minimum_step_s <= 0.0
            || self.maximum_step_s < self.minimum_step_s
            || self.maximum_attempts == 0
            || self.maximum_evaluations == 0
            || !self.scaled_rate_tolerance_m_s.is_finite()
            || self.scaled_rate_tolerance_m_s <= 0.0
            || !self.rest_scaled_rate_tolerance_m_s.is_finite()
            || self.rest_scaled_rate_tolerance_m_s < 0.0
        {
            return Err(Error::Invalid(
                "finite driven tolerances and bounded work required",
            ));
        }
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn checked(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid("unrepresentable driven quantity"))
    }
}
fn vector(x: Vector) -> Result<Vector, Error> {
    for value in x {
        checked(value)?;
    }
    Ok(x)
}
fn add(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| a[i] + b[i])
}
fn subtract(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: Vector, s: f64) -> Vector {
    a.map(|v| s * v)
}
fn norm(a: Vector) -> f64 {
    det::hypot(det::hypot(a[0], a[1]), a[2])
}
