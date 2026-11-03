//! Rigid distance between material points on two freely moving finite bodies.
//!
//! Initial and endpoint central impulses enforce separation and range rate,
//! respectively, around the original full-inertia rigid steps. Each reaction
//! uses the native equal-and-opposite impulse owner; both bodies recoil, and
//! internal reactions preserve total world linear and angular momentum. No
//! coordinate is snapped, no body is clamped, and no stiff spring is inserted.
//!
//! This follows the primary/hidden constraint ordering of SHAKE/RATTLE; see
//! McLachlan et al., *Geometric Generalisations of SHAKE and RATTLE*,
//! <https://arxiv.org/abs/1207.3367>. The original free map is symmetric, giving
//! second-order motion on the admitted local branch. This does not assert
//! symplecticity, exact energy conservation, or a global root certificate.

use super::pair_step::Geometry;
use crate::spatial_loads::{PreparedSpatialLoads, SpatialStepReceipt};
use crate::{
    ActionReactionImpulseReceipt, DynamicsError, MassProperties, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, Vec3, Wrench, apply_equal_and_opposite_impulse_at_body_points,
};

/// Simultaneous distance constraints on a graph of independent finite bodies.
#[path = "distance/graph.rs"]
pub mod graph;

/// Original dimensioned impulse, angular-travel and nonlinear-work limits.
/// The angular cap applies to each body and to the connection direction.
pub use super::step::PointPlaneStepControls as DistanceStepControls;

/// One bilateral, strictly positive material-point separation [m].
///
/// The points are measured from each body's COM in its original body frame.
/// The law permits tension and compression. A zero-distance ball joint,
/// collision discovery, unilateral cable, friction and arbitrary joint graph
/// are outside this independent native-body pair entry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidDistanceConstraint {
    properties: [MassProperties; 2],
    points_body_m: [Vec3; 2],
    length_m: f64,
    position_tolerance_m: f64,
    velocity_tolerance_mps: f64,
}

/// Both accepted native endpoints and their physical constraint impulses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceStepReceipt<FreeStep = StepReceipt> {
    /// Final bodies in constructor order, with both constraints satisfied.
    pub states_after: [RigidBodyState; 2],
    /// Central equal-and-opposite impulses at the initial material points.
    pub position_impulse: ActionReactionImpulseReceipt,
    /// Native smooth steps before the final tangent impulses.
    pub free_steps: [FreeStep; 2],
    /// Central equal-and-opposite impulses at the endpoint material points.
    pub velocity_impulse: ActionReactionImpulseReceipt,
    /// Sum of both native pair impulse works [J]. This discrete constraint
    /// work remains explicit; it is not heat or an energy projection.
    pub constraint_work_j: f64,
    /// Actual endpoint length minus the declared distance [m].
    pub final_gap_m: f64,
    /// Actual endpoint rate of separation [m/s].
    pub final_range_rate_mps: f64,
    /// Accepted position corrections; each uses bounded backtracking.
    pub position_iterations: u32,
}

/// A refusal publishes neither body, including late cancellation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DistanceStepError {
    /// Invalid controls, collapsed geometry, or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Initial positions or velocities do not satisfy the declared constraint.
    InitialConstraint,
    /// The bounded local impulse/rotation branch was exhausted.
    LocalBranch,
    /// Position or velocity closure did not meet its physical tolerance.
    NotConverged,
    /// Original native rigid-body or pair-impulse refusal.
    Dynamics(DynamicsError),
    /// Cancellation observed before publication of the pair.
    Cancelled,
}

impl From<DynamicsError> for DistanceStepError {
    fn from(value: DynamicsError) -> Self {
        Self::Dynamics(value)
    }
}
impl core::fmt::Display for DistanceStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid distance motion: {self:?}")
    }
}
impl std::error::Error for DistanceStepError {}

impl RigidDistanceConstraint {
    /// Prepare one distance and absolute position/range-rate tolerances in SI.
    ///
    /// # Errors
    /// All points and scalars must be finite. Length and tolerances must be
    /// positive, and position tolerance must be strictly smaller than length.
    pub fn new(
        properties: [MassProperties; 2],
        points_body_m: [Vec3; 2],
        length_m: f64,
        position_tolerance_m: f64,
        velocity_tolerance_mps: f64,
    ) -> Result<Self, DistanceStepError> {
        if points_body_m.iter().any(|p| !p.is_finite())
            || !positive(length_m)
            || !positive(position_tolerance_m)
            || position_tolerance_m >= length_m
            || !positive(velocity_tolerance_mps)
        {
            return Err(DistanceStepError::Invalid(
                "distance geometry or tolerances",
            ));
        }
        Ok(Self {
            properties,
            points_body_m,
            length_m,
            position_tolerance_m,
            velocity_tolerance_mps,
        })
    }

    /// Advance both already constrained bodies under their held native wrenches.
    ///
    /// Gravity is the shared uniform field of `integrator`; each supplied
    /// force is world-fixed and each torque body-fixed. Internal impulses use
    /// the original and final connection directions respectively, retaining
    /// all products of inertia and off-center spin/orbital exchange.
    ///
    /// The position solve starts from zero impulse and uses the original
    /// pair mobility, with bounded residual-decreasing corrections. Each
    /// body and the connection direction must stay within the angular cap.
    /// This admits a local numerical branch, not a whole-path enclosure.
    /// Refine the physical timestep if local bounds are exhausted. No heap
    /// allocation or intermediate state publication occurs.
    pub fn advance(
        &self,
        states: [RigidBodyState; 2],
        integrator: RigidBodyIntegrator,
        wrenches: [Wrench; 2],
        duration_s: f64,
        controls: DistanceStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DistanceStepReceipt, DistanceStepError> {
        super::pair_step::advance(
            self, states, integrator, wrenches, duration_s, controls, cancelled,
        )
    }

    /// Advance under prepared world-held point forces and body/world couples.
    ///
    /// Every constraint trial uses the native spatial-load step, reevaluating
    /// moments at the actual endpoint attitudes. Preparation, free rotation,
    /// and work accounting remain owned by [`PreparedSpatialLoads`] and
    /// [`RigidBodyIntegrator::step_spatial`]. The supplied loads are held over
    /// this interval; they are not a state-dependent spring or drag callback.
    ///
    /// The complete native spatial receipts retain applied world linear and
    /// angular impulse, force/couple work, and their numerical defects. Internal
    /// constraint reactions are separate. The original local solve, no-pose-
    /// repair behavior, and atomic cancellation are unchanged. Torque bounds
    /// are prepared once from the loads; all trials remain allocation-free.
    pub fn advance_spatial(
        &self,
        states: [RigidBodyState; 2],
        integrator: RigidBodyIntegrator,
        loads: [PreparedSpatialLoads; 2],
        duration_s: f64,
        controls: DistanceStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DistanceStepReceipt<SpatialStepReceipt>, DistanceStepError> {
        super::pair_step::advance_spatial(
            self, states, integrator, loads, duration_s, controls, cancelled,
        )
    }

    fn geometry(&self, states: [RigidBodyState; 2]) -> Result<Geometry, DistanceStepError> {
        let a = states[0].point_kinematics(self.properties[0], self.points_body_m[0])?;
        let b = states[1].point_kinematics(self.properties[1], self.points_body_m[1])?;
        // Subtract COMs before adding their local arms: a common distant origin
        // must not round away the physical separation of two attachment points.
        let delta = states[1]
            .pose()
            .position_world()
            .sub(states[0].pose().position_world())
            .add(b.arm_world.sub(a.arm_world));
        let scale = delta.x.abs().max(delta.y.abs()).max(delta.z.abs());
        if !delta.is_finite() || !positive(scale) {
            return Err(DistanceStepError::Invalid(
                "collapsed or unrepresentable separation",
            ));
        }
        let scaled = Vec3::new(delta.x / scale, delta.y / scale, delta.z / scale);
        let norm = fs_math::det::sqrt(scaled.norm_squared());
        let normal = scaled.scale(norm.recip());
        let gap = scale * norm - self.length_m;
        let rate = normal.dot(b.point_velocity_world.sub(a.point_velocity_world));
        if !gap.is_finite() || !rate.is_finite() {
            return Err(DistanceStepError::Invalid("distance or range rate"));
        }
        Ok(Geometry {
            normal,
            arms_body_m: self.points_body_m,
            gap,
            rate,
        })
    }

    fn pair_impulse(
        &self,
        states: [RigidBodyState; 2],
        impulse: Vec3,
    ) -> Result<ActionReactionImpulseReceipt, DistanceStepError> {
        Ok(apply_equal_and_opposite_impulse_at_body_points(
            states[0],
            self.properties[0],
            self.points_body_m[0],
            states[1],
            self.properties[1],
            self.points_body_m[1],
            impulse,
        )?)
    }
}

impl super::pair_step::Constraint for RigidDistanceConstraint {
    fn properties(&self) -> [MassProperties; 2] {
        self.properties
    }
    fn tolerances(&self) -> [f64; 2] {
        [self.position_tolerance_m, self.velocity_tolerance_mps]
    }
    fn geometry(
        &self,
        states: [RigidBodyState; 2],
    ) -> Result<super::pair_step::Geometry, DistanceStepError> {
        self.geometry(states)
    }
}

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DistanceStepError> {
    if cancelled() {
        Err(DistanceStepError::Cancelled)
    } else {
        Ok(())
    }
}
