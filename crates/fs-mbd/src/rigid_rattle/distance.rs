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
pub struct DistanceStepReceipt {
    /// Final bodies in constructor order, with both constraints satisfied.
    pub states_after: [RigidBodyState; 2],
    /// Central equal-and-opposite impulses at the initial material points.
    pub position_impulse: ActionReactionImpulseReceipt,
    /// Native smooth steps before the final tangent impulses.
    pub free_steps: [StepReceipt; 2],
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DistanceStepReceipt, DistanceStepError> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !positive(controls.maximum_impulse_ns)
            || !positive(controls.maximum_rotation_increment_rad)
            || controls.maximum_rotation_increment_rad > 0.25
            || controls.maximum_position_iterations == 0
            || controls.maximum_backtracks == 0
        {
            return Err(DistanceStepError::Invalid(
                "positive finite local controls required",
            ));
        }
        let initial = self.geometry(states)?;
        if initial.gap.abs() > self.position_tolerance_m
            || initial.rate.abs() > self.velocity_tolerance_mps
        {
            return Err(DistanceStepError::InitialConstraint);
        }
        let response = duration_s * self.mobility(states, initial.normal)?;
        if !positive(response) {
            return Err(DistanceStepError::Invalid("distance impulse response"));
        }
        let problem = Problem {
            constraint: self,
            states,
            integrator,
            wrenches,
            duration_s,
            controls,
            normal: initial.normal,
        };
        let mut impulse = 0.0;
        let mut trial = problem.evaluate(impulse, &mut cancelled)?;
        let mut iterations = 0;
        while trial.geometry.gap.abs() > self.position_tolerance_m {
            if iterations == controls.maximum_position_iterations {
                return Err(DistanceStepError::NotConverged);
            }
            let mut correction = -trial.geometry.gap / response;
            if !correction.is_finite() {
                return Err(DistanceStepError::LocalBranch);
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
                        Ok(candidate)
                            if candidate.geometry.gap.abs() < trial.geometry.gap.abs() =>
                        {
                            accepted = Some((next_impulse, candidate));
                            break;
                        }
                        Ok(_) | Err(DistanceStepError::LocalBranch) => {}
                        Err(error) => return Err(error),
                    }
                }
                correction *= 0.5;
            }
            let Some((next_impulse, candidate)) = accepted else {
                return Err(DistanceStepError::LocalBranch);
            };
            impulse = next_impulse;
            trial = candidate;
            iterations += 1;
        }
        poll(&mut cancelled)?;
        let free = trial.steps.map(|s| s.state_after);
        let final_impulse = -trial.geometry.rate / self.mobility(free, trial.geometry.normal)?;
        if !final_impulse.is_finite() || final_impulse.abs() > controls.maximum_impulse_ns {
            return Err(DistanceStepError::LocalBranch);
        }
        let velocity_impulse =
            self.pair_impulse(free, trial.geometry.normal.scale(final_impulse))?;
        let states_after = [
            velocity_impulse.first.state_after,
            velocity_impulse.second.state_after,
        ];
        let final_geometry = self.geometry(states_after)?;
        if final_geometry.gap.abs() > self.position_tolerance_m
            || final_geometry.rate.abs() > self.velocity_tolerance_mps
        {
            return Err(DistanceStepError::NotConverged);
        }
        let constraint_work_j = trial.impulse.impulse_work_j + velocity_impulse.impulse_work_j;
        if !constraint_work_j.is_finite() {
            return Err(DistanceStepError::Invalid("distance constraint work"));
        }
        poll(&mut cancelled)?;
        Ok(DistanceStepReceipt {
            states_after,
            position_impulse: trial.impulse,
            free_steps: trial.steps,
            velocity_impulse,
            constraint_work_j,
            final_gap_m: final_geometry.gap,
            final_range_rate_mps: final_geometry.rate,
            position_iterations: iterations,
        })
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
        Ok(Geometry { normal, gap, rate })
    }

    fn mobility(
        &self,
        states: [RigidBodyState; 2],
        normal: Vec3,
    ) -> Result<f64, DistanceStepError> {
        let a = states[0].directional_effective_mass(
            self.properties[0],
            self.points_body_m[0],
            normal,
        )?;
        let b = states[1].directional_effective_mass(
            self.properties[1],
            self.points_body_m[1],
            normal,
        )?;
        let result = a.inverse_mass_kg_inverse + b.inverse_mass_kg_inverse;
        if !positive(result) {
            return Err(DistanceStepError::Invalid("pair mobility"));
        }
        Ok(result)
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

struct Geometry {
    normal: Vec3,
    gap: f64,
    rate: f64,
}
struct Trial {
    impulse: ActionReactionImpulseReceipt,
    steps: [StepReceipt; 2],
    geometry: Geometry,
}
struct Problem<'a> {
    constraint: &'a RigidDistanceConstraint,
    states: [RigidBodyState; 2],
    integrator: RigidBodyIntegrator,
    wrenches: [Wrench; 2],
    duration_s: f64,
    controls: DistanceStepControls,
    normal: Vec3,
}
impl Problem<'_> {
    fn evaluate(
        &self,
        impulse_ns: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, DistanceStepError> {
        poll(cancelled)?;
        let impulse = self
            .constraint
            .pair_impulse(self.states, self.normal.scale(impulse_ns))?;
        let states = [impulse.first.state_after, impulse.second.state_after];
        for i in 0..2 {
            let principal = self.constraint.properties[i].principal_inertia_body();
            let minimum = principal.x.min(principal.y).min(principal.z);
            let travel = (norm(states[i].angular_momentum_body())
                + self.duration_s * norm(self.wrenches[i].torque_body))
                * (self.duration_s / minimum);
            if !travel.is_finite() || travel > self.controls.maximum_rotation_increment_rad {
                return Err(DistanceStepError::LocalBranch);
            }
        }
        let first = self.integrator.step(
            states[0],
            self.constraint.properties[0],
            self.wrenches[0],
            self.duration_s,
        )?;
        poll(cancelled)?;
        let second = self.integrator.step(
            states[1],
            self.constraint.properties[1],
            self.wrenches[1],
            self.duration_s,
        )?;
        poll(cancelled)?;
        let geometry = self
            .constraint
            .geometry([first.state_after, second.state_after])?;
        let max_chord = 2.0 * fs_math::det::sin(0.5 * self.controls.maximum_rotation_increment_rad);
        if norm(geometry.normal.sub(self.normal)) > max_chord {
            return Err(DistanceStepError::LocalBranch);
        }
        Ok(Trial {
            impulse,
            steps: [first, second],
            geometry,
        })
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
