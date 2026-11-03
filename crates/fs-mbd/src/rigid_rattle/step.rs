//! Bilateral point/plane motion through the original native rigid-body step.
//!
//! An initial normal point impulse is solved so the native free endpoint lies
//! on the plane; a second endpoint point impulse removes its normal speed.
//! Poses are never repaired. The construction follows the position/velocity
//! constraint ordering of Andersen, JCP 52 (1983), 24–34,
//! <https://doi.org/10.1016/0021-9991(83)90014-1>.
//!
//! The scalar solve follows local residual-decreasing corrections from zero
//! impulse, with bounded recoil and angular travel. It refuses failure of this
//! local solve instead of searching for a remote zero after a large rotation.
//! This is a smooth bilateral constraint: tensile reactions are allowed. It
//! supplies no unilateral contact, friction, impact or moving-plane law.
//! The original free map is symmetric; this composition is second order on
//! the locally connected smooth branch. Neither exact energy preservation nor
//! symplecticity of the original implicit rigid-body map is claimed.

use super::PointPlaneRattle;
use crate::{
    DynamicsError, ImpulseReceipt, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench,
};

/// Work and local-branch limits, independent of the projector's gap/speed tolerances.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPlaneStepControls {
    /// Positive cap on each initial or endpoint normal impulse [N s].
    pub maximum_impulse_ns: f64,
    /// Positive cap at most 0.25 rad on a conservative angular-travel bound.
    /// Reduce the physical timestep when this bound is exhausted.
    pub maximum_rotation_increment_rad: f64,
    /// Maximum accepted local position corrections; must be nonzero.
    pub maximum_position_iterations: u32,
    /// Maximum halvings/trials of each correction; must be nonzero.
    pub maximum_backtracks: u32,
}

/// Complete accepted motion, including both actual physical impulse events.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPlaneStepReceipt {
    /// Native momentum-form endpoint, with both constraints satisfied.
    pub state_after: RigidBodyState,
    /// Impulse at the original material point before free evolution.
    pub position_impulse: ImpulseReceipt,
    /// Original free-motion result between the two impulses.
    pub free_step: StepReceipt,
    /// Impulse at the actual endpoint material point.
    pub velocity_impulse: ImpulseReceipt,
    /// Sum of both measured midpoint impulse works [J]. This discrete work is
    /// not heat and is not projected away to manufacture energy conservation.
    pub constraint_work_j: f64,
    /// Actual endpoint signed gap [m].
    pub final_gap_m: f64,
    /// Actual endpoint signed material-point normal speed [m/s].
    pub final_normal_speed_mps: f64,
    /// Accepted local corrections; each costs at most `maximum_backtracks` free steps.
    pub position_iterations: u32,
}

/// Refusal returns no updated state, including late cancellation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PointPlaneStepError {
    /// Invalid duration, control or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The supplied initial state is not on the position and velocity manifold.
    InitialConstraint,
    /// A local impulse/rotation limit was exhausted; refine the physical step.
    LocalBranch,
    /// Bounded local corrections could not satisfy the endpoint gap or speed.
    NotConverged,
    /// Original native inertia, impulse or rigid-step refusal.
    Dynamics(DynamicsError),
    /// Caller cancelled before publication.
    Cancelled,
}

impl From<DynamicsError> for PointPlaneStepError {
    fn from(value: DynamicsError) -> Self {
        Self::Dynamics(value)
    }
}
impl core::fmt::Display for PointPlaneStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid point/plane constrained motion: {self:?}")
    }
}
impl std::error::Error for PointPlaneStepError {}

impl PointPlaneRattle {
    /// Advance an already constrained native rigid state under held forcing.
    ///
    /// The force is world-fixed and torque body-fixed, as in the original
    /// [`RigidBodyIntegrator::step`]. Each reaction is an actual impulse at the
    /// material point in its own endpoint attitude. In the absence of external
    /// symmetry-breaking loads, tangent linear momentum and the plane-normal
    /// component of total spin-plus-orbital angular momentum are preserved.
    /// No position displacement, basis repair, or momentum projection is used.
    ///
    /// The scalar iteration uses the original normal mobility as its local
    /// preconditioner and accepts only residual-decreasing backtracked steps.
    /// Finite impulse and small angular-travel limits are required. This is
    /// a local numerical branch admission, not a global root uniqueness proof.
    /// A refusal requires timestep refinement or explicitly revised budgets.
    pub fn advance(
        &self,
        state: RigidBodyState,
        integrator: RigidBodyIntegrator,
        wrench: Wrench,
        duration_s: f64,
        controls: PointPlaneStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PointPlaneStepReceipt, PointPlaneStepError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !controls.maximum_impulse_ns.is_finite()
            || controls.maximum_impulse_ns <= 0.0
            || !controls.maximum_rotation_increment_rad.is_finite()
            || controls.maximum_rotation_increment_rad <= 0.0
            || controls.maximum_rotation_increment_rad > 0.25
            || controls.maximum_position_iterations == 0
            || controls.maximum_backtracks == 0
        {
            return Err(PointPlaneStepError::Invalid(
                "positive finite local controls required",
            ));
        }
        let native = |v: fs_ga::Vec3| Vec3::new(v.x, v.y, v.z);
        let problem = Problem {
            constraint: self,
            state,
            integrator,
            wrench,
            duration_s,
            controls,
            normal: native(self.normal_world),
            arm: native(self.point_body_m),
            plane_point: native(self.plane_point_world_m),
        };
        let (gap, speed) = problem.observables(state)?;
        if gap.abs() > self.controls.position_tolerance_m
            || speed.abs() > self.controls.velocity_tolerance_mps
        {
            return Err(PointPlaneStepError::InitialConstraint);
        }
        let mobility = state
            .directional_effective_mass(self.properties, problem.arm, problem.normal)?
            .inverse_mass_kg_inverse;
        let response = duration_s * mobility;
        if !response.is_finite() || response <= 0.0 {
            return Err(PointPlaneStepError::Invalid("position impulse response"));
        }
        let mut impulse = 0.0;
        let mut trial = problem.evaluate(impulse, &mut cancelled)?;
        let mut iterations = 0;
        while trial.gap.abs() > self.controls.position_tolerance_m {
            if iterations == controls.maximum_position_iterations {
                return Err(PointPlaneStepError::NotConverged);
            }
            let mut correction = -trial.gap / response;
            if !correction.is_finite() {
                return Err(PointPlaneStepError::LocalBranch);
            }
            let mut accepted = None;
            for _ in 0..controls.maximum_backtracks {
                poll(&mut cancelled)?;
                let next_impulse = impulse + correction;
                if next_impulse == impulse {
                    break;
                }
                if next_impulse.abs() <= controls.maximum_impulse_ns {
                    match problem.evaluate(next_impulse, &mut cancelled) {
                        Ok(candidate) if candidate.gap.abs() < trial.gap.abs() => {
                            accepted = Some((next_impulse, candidate));
                            break;
                        }
                        Ok(_) | Err(PointPlaneStepError::LocalBranch) => {}
                        Err(error) => return Err(error),
                    }
                }
                correction *= 0.5;
            }
            let Some((next_impulse, candidate)) = accepted else {
                return Err(PointPlaneStepError::LocalBranch);
            };
            impulse = next_impulse;
            trial = candidate;
            iterations += 1;
        }
        poll(&mut cancelled)?;
        let free = trial.motion.state_after;
        let (_, speed) = problem.observables(free)?;
        let endpoint_mass = free
            .directional_effective_mass(self.properties, problem.arm, problem.normal)?
            .effective_mass_kg;
        let final_impulse = -speed * endpoint_mass;
        if !final_impulse.is_finite() || final_impulse.abs() > controls.maximum_impulse_ns {
            return Err(PointPlaneStepError::LocalBranch);
        }
        let velocity_impulse = free.apply_impulse_at_body_point(
            self.properties,
            problem.arm,
            problem.normal.scale(final_impulse),
        )?;
        let state_after = velocity_impulse.state_after;
        let (final_gap_m, final_normal_speed_mps) = problem.observables(state_after)?;
        if final_gap_m.abs() > self.controls.position_tolerance_m
            || final_normal_speed_mps.abs() > self.controls.velocity_tolerance_mps
        {
            return Err(PointPlaneStepError::NotConverged);
        }
        let constraint_work_j =
            trial.impulse.work.impulse_work_j + velocity_impulse.work.impulse_work_j;
        if !constraint_work_j.is_finite() {
            return Err(PointPlaneStepError::Invalid("constraint impulse work"));
        }
        poll(&mut cancelled)?;
        Ok(PointPlaneStepReceipt {
            state_after,
            position_impulse: trial.impulse,
            free_step: trial.motion,
            velocity_impulse,
            constraint_work_j,
            final_gap_m,
            final_normal_speed_mps,
            position_iterations: iterations,
        })
    }
}

struct Problem<'a> {
    constraint: &'a PointPlaneRattle,
    state: RigidBodyState,
    integrator: RigidBodyIntegrator,
    wrench: Wrench,
    duration_s: f64,
    controls: PointPlaneStepControls,
    normal: Vec3,
    arm: Vec3,
    plane_point: Vec3,
}
struct Trial {
    impulse: ImpulseReceipt,
    motion: StepReceipt,
    gap: f64,
}
impl Problem<'_> {
    fn observables(&self, state: RigidBodyState) -> Result<(f64, f64), PointPlaneStepError> {
        let point = state.point_kinematics(self.constraint.properties, self.arm)?;
        let gap = self.normal.dot(point.point_world.sub(self.plane_point));
        let speed = self.normal.dot(point.point_velocity_world);
        if !gap.is_finite() || !speed.is_finite() {
            return Err(PointPlaneStepError::Invalid("constraint gap or speed"));
        }
        Ok((gap, speed))
    }

    fn evaluate(
        &self,
        impulse_ns: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, PointPlaneStepError> {
        poll(cancelled)?;
        let impulse = self.state.apply_impulse_at_body_point(
            self.constraint.properties,
            self.arm,
            self.normal.scale(impulse_ns),
        )?;
        let moments = self.constraint.properties.principal_inertia_body();
        let minimum = moments.x.min(moments.y).min(moments.z);
        let norm = |v: Vec3| v.x.hypot(v.y).hypot(v.z);
        let angular_travel = (norm(impulse.state_after.angular_momentum_body())
            + self.duration_s * norm(self.wrench.torque_body))
            * (self.duration_s / minimum);
        if !angular_travel.is_finite()
            || angular_travel > self.controls.maximum_rotation_increment_rad
        {
            return Err(PointPlaneStepError::LocalBranch);
        }
        let motion = self.integrator.step(
            impulse.state_after,
            self.constraint.properties,
            self.wrench,
            self.duration_s,
        )?;
        poll(cancelled)?;
        let (gap, _) = self.observables(motion.state_after)?;
        Ok(Trial {
            impulse,
            motion,
            gap,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PointPlaneStepError> {
    if cancelled() {
        Err(PointPlaneStepError::Cancelled)
    } else {
        Ok(())
    }
}
