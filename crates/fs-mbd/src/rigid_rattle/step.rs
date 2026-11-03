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
//! also supports prescribed translation and rigidly rotating plane supports.
//! It supplies no unilateral contact, friction or impact law.
//! The original free map is symmetric; this composition is second order on
//! the locally connected smooth branch. Neither exact energy preservation nor
//! symplecticity of the original implicit rigid-body map is claimed.

use super::PointPlaneRattle;
use crate::{
    DynamicsError, ImpulseReceipt, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_math::dd::Dd;

#[path = "step/rotating.rs"]
mod rotating;
pub use rotating::{PlaneRigidMotion, RotatingPlaneStepReceipt};

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
    /// Actual endpoint signed material-point speed relative to the plane [m/s].
    pub final_normal_speed_mps: f64,
    /// Accepted local corrections; each costs at most `maximum_backtracks` free steps.
    pub position_iterations: u32,
}

/// Prescribed plane translation along its fixed world normal for one step.
/// Tangential translation leaves this frictionless infinite plane unchanged.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PlaneTranslation {
    /// Initial signed plane speed along its original unit normal [m/s].
    pub normal_velocity_mps: f64,
    /// Held signed normal acceleration [m/s²].
    pub normal_acceleration_mps2: f64,
}

/// One native rigid step and the prescribed support's reciprocal account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TranslatingPlaneStepReceipt {
    /// Original impulse/free-step/impulse result, with relative gap and speed.
    pub step: PointPlaneStepReceipt,
    /// Same constraint with its plane point advanced to the accepted endpoint.
    pub plane_after: PointPlaneRattle,
    /// Plane speed to use as the initial speed of the next step [m/s].
    pub normal_velocity_after_mps: f64,
    /// Represented normal plane displacement minus its prescribed value [m].
    /// Admitted against the original position tolerance, never hidden by repair.
    pub plane_displacement_residual_m: f64,
    /// Represented endpoint plane speed minus its prescribed value [m/s].
    pub plane_velocity_residual_mps: f64,
    /// Equal-and-opposite sum of the body's two reaction impulses [N s].
    /// The individual application points remain in `step`'s impulse receipts.
    pub support_impulse_world_ns: Vec3,
    /// Signed prescribed-boundary work delivered to the rigid body [J].
    /// Each impulse is paired with the support velocity at that event.
    pub boundary_work_j: f64,
    /// Actual constraint impulse work minus prescribed-boundary work [J],
    /// independently paired in the support's relative frame. This is a
    /// time-discretization defect, not frictional heat or an energy correction.
    pub constraint_work_defect_j: f64,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<PointPlaneStepReceipt, PointPlaneStepError> {
        self.advance_between(
            state,
            integrator,
            wrench,
            duration_s,
            controls,
            self.plane_point_world_m,
            self.normal_world,
            [PlaneSpeed::Normal(0.0); 2],
            cancelled,
        )
    }

    /// Advance against a translating, normally accelerating prescribed plane.
    ///
    /// The plane point in `self` is its initial point. The endpoint is
    /// `p + n*(v*dt + a*dt²/2)` and the hidden velocity constraint is
    /// `n·v_point = v + a*dt`. Both use the original fixed normal and actual
    /// full-inertia native impulse responses. The held world force/body torque
    /// have the same meaning as in [`Self::advance`].
    ///
    /// Return `plane_after` and `normal_velocity_after_mps` for continuation.
    /// Prescribed motion has an external source of momentum and work; no
    /// finite support inertia is inferred. Reactions remain bilateral. The
    /// stationary law recovers the original solver and native state exactly.
    /// No heap allocation or partial state publication occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_translating_plane(
        &self,
        state: RigidBodyState,
        integrator: RigidBodyIntegrator,
        wrench: Wrench,
        duration_s: f64,
        translation: PlaneTranslation,
        controls: PointPlaneStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TranslatingPlaneStepReceipt, PointPlaneStepError> {
        poll(&mut cancelled)?;
        let initial_speed = translation.normal_velocity_mps;
        let acceleration = translation.normal_acceleration_mps2;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !initial_speed.is_finite()
            || !acceleration.is_finite()
        {
            return Err(PointPlaneStepError::Invalid(
                "finite plane translation required",
            ));
        }
        // Keep the acceleration contribution until representability admission:
        // a small dv beside a large v must not disappear before it is checked.
        let h = Dd::from_f64(duration_s);
        let dv = Dd::from_f64(acceleration) * h;
        let prescribed_speed = Dd::from_f64(initial_speed) + dv;
        let prescribed_displacement = h * (Dd::from_f64(initial_speed) + Dd::from_f64(0.5) * dv);
        let final_speed = prescribed_speed.to_f64();
        let displacement = prescribed_displacement.to_f64();
        let velocity_residual = (Dd::from_f64(final_speed) - prescribed_speed).to_f64();
        let n = Vec3::new(
            self.normal_world.x,
            self.normal_world.y,
            self.normal_world.z,
        );
        let p = Vec3::new(
            self.plane_point_world_m.x,
            self.plane_point_world_m.y,
            self.plane_point_world_m.z,
        );
        let endpoint = p.add(n.scale(displacement));
        let represented_displacement = [n.x, n.y, n.z]
            .into_iter()
            .zip([endpoint.x, endpoint.y, endpoint.z])
            .zip([p.x, p.y, p.z])
            .fold(Dd::ZERO, |sum, ((normal, end), start)| {
                sum + Dd::from_f64(normal) * (Dd::from_f64(end) - Dd::from_f64(start))
            });
        let displacement_residual = (represented_displacement - prescribed_displacement).to_f64();
        if !final_speed.is_finite()
            || !displacement.is_finite()
            || !endpoint.is_finite()
            || !displacement_residual.is_finite()
            || !velocity_residual.is_finite()
            || displacement_residual.abs() > self.controls.position_tolerance_m
            || velocity_residual.abs() > self.controls.velocity_tolerance_mps
            || (displacement != 0.0 && n.dot(endpoint.sub(p)) == 0.0)
        {
            return Err(PointPlaneStepError::Invalid(
                "plane endpoint is unrepresentable",
            ));
        }
        let mut plane_after = *self;
        plane_after.plane_point_world_m = fs_ga::Vec3::new(endpoint.x, endpoint.y, endpoint.z);
        let step = self.advance_between(
            state,
            integrator,
            wrench,
            duration_s,
            controls,
            plane_after.plane_point_world_m,
            plane_after.normal_world,
            [
                PlaneSpeed::Normal(initial_speed),
                PlaneSpeed::Normal(final_speed),
            ],
            &mut cancelled,
        )?;
        if (Dd::from_f64(step.final_gap_m) + Dd::from_f64(displacement_residual))
            .to_f64()
            .abs()
            > self.controls.position_tolerance_m
            || (Dd::from_f64(step.final_normal_speed_mps) + Dd::from_f64(velocity_residual))
                .to_f64()
                .abs()
                > self.controls.velocity_tolerance_mps
        {
            return Err(PointPlaneStepError::NotConverged);
        }
        let mut source = Dd::ZERO;
        let mut defect = Dd::ZERO;
        for (event, speed) in [
            (step.position_impulse, Dd::from_f64(initial_speed)),
            (step.velocity_impulse, prescribed_speed),
        ] {
            let before = event.point_kinematics_before.point_velocity_world;
            let after = event.point_kinematics_after.point_velocity_world;
            for ((normal, impulse), (a, b)) in [n.x, n.y, n.z]
                .into_iter()
                .zip([
                    event.impulse_world.x,
                    event.impulse_world.y,
                    event.impulse_world.z,
                ])
                .zip(
                    [before.x, before.y, before.z]
                        .into_iter()
                        .zip([after.x, after.y, after.z]),
                )
            {
                let support = Dd::from_f64(normal) * speed;
                let j = Dd::from_f64(impulse);
                source = source + support * j;
                defect = defect
                    + ((Dd::from_f64(a) + Dd::from_f64(b)) * Dd::from_f64(0.5) - support) * j;
            }
        }
        let support_impulse_world_ns = step
            .position_impulse
            .impulse_world
            .add(step.velocity_impulse.impulse_world)
            .scale(-1.0);
        let boundary_work_j = source.to_f64();
        let constraint_work_defect_j = defect.to_f64();
        if !support_impulse_world_ns.is_finite()
            || !boundary_work_j.is_finite()
            || !constraint_work_defect_j.is_finite()
        {
            return Err(PointPlaneStepError::Invalid("prescribed support account"));
        }
        poll(&mut cancelled)?;
        Ok(TranslatingPlaneStepReceipt {
            step,
            plane_after,
            normal_velocity_after_mps: final_speed,
            plane_displacement_residual_m: displacement_residual,
            plane_velocity_residual_mps: velocity_residual,
            support_impulse_world_ns,
            boundary_work_j,
            constraint_work_defect_j,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_between(
        &self,
        state: RigidBodyState,
        integrator: RigidBodyIntegrator,
        wrench: Wrench,
        duration_s: f64,
        controls: PointPlaneStepControls,
        endpoint_plane: fs_ga::Vec3,
        endpoint_normal: fs_ga::Vec3,
        plane_speeds: [PlaneSpeed; 2],
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
            initial_normal: native(self.normal_world),
            normal: native(endpoint_normal),
            arm: native(self.point_body_m),
            plane_point: native(endpoint_plane),
            plane_speed: plane_speeds[1],
        };
        let (gap, speed) = problem.observables_at(
            state,
            native(self.plane_point_world_m),
            problem.initial_normal,
            plane_speeds[0],
        )?;
        if gap.abs() > self.controls.position_tolerance_m
            || speed.abs() > self.controls.velocity_tolerance_mps
        {
            return Err(PointPlaneStepError::InitialConstraint);
        }
        let mobility = state
            .directional_effective_mass(self.properties, problem.arm, problem.initial_normal)?
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
    initial_normal: Vec3,
    normal: Vec3,
    arm: Vec3,
    plane_point: Vec3,
    plane_speed: PlaneSpeed,
}

#[derive(Clone, Copy)]
enum PlaneSpeed {
    // Keep the established stationary/translating arithmetic unchanged.
    Normal(f64),
    Rigid {
        pivot: [Dd; 3],
        velocity: Vec3,
        angular_velocity: Vec3,
    },
}
struct Trial {
    impulse: ImpulseReceipt,
    motion: StepReceipt,
    gap: f64,
}
impl Problem<'_> {
    fn observables(&self, state: RigidBodyState) -> Result<(f64, f64), PointPlaneStepError> {
        self.observables_at(state, self.plane_point, self.normal, self.plane_speed)
    }

    fn observables_at(
        &self,
        state: RigidBodyState,
        plane_point: Vec3,
        normal: Vec3,
        plane_speed: PlaneSpeed,
    ) -> Result<(f64, f64), PointPlaneStepError> {
        let point = state.point_kinematics(self.constraint.properties, self.arm)?;
        let (gap, speed) = match plane_speed {
            PlaneSpeed::Normal(value) => (
                normal.dot(point.point_world.sub(plane_point)),
                normal.dot(point.point_velocity_world) - value,
            ),
            PlaneSpeed::Rigid { .. } => (
                rotating::gap(
                    rotating::dd_vec(normal),
                    state,
                    point.arm_world,
                    rotating::dd_vec(plane_point),
                ),
                rotating::relative_speed(
                    rotating::dd_vec(normal),
                    rotating::point_velocity(point),
                    plane_speed.velocity_at(state, point),
                ),
            ),
        };
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
            self.initial_normal.scale(impulse_ns),
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
