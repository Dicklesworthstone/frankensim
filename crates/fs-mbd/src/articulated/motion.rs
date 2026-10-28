//! Continuous free-floating articulated motion with physical base recoil.
//!
//! A second-order Lie midpoint step evolves total world linear momentum,
//! angular momentum about the system COM, COM position, base orientation, and
//! joint coordinates/rates. Native ABA supplies joint accelerations. At each
//! stage a locked composite inertia recovers the base twist at the declared
//! joint rates; the articulated Schur complement is not the appropriate map.
//! Internal actuator efforts therefore cannot create total momentum.
//!
//! This is an explicit smooth-flight lane, without contact or loop constraints.
//! It preserves the discrete spatial impulse balance, but is neither an
//! energy-exact nor a variational/symplectic integrator. Refine the duration to
//! resolve fast internal motion. Held link-frame wrenches rotate with the links.

use super::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, FreeFloatingBaseState,
    Kinematics, SpatialInertia, forward_kinematics, free_floating_forward_dynamics,
    free_floating_momentum, solve_floating_base_system, validate_external_wrenches,
    validate_forward_dynamics_header, validate_slice, validate_twist, validate_vec3,
};
use core::fmt;
use fs_ga::{GaError, Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};

/// Maximum model size admitted by one cancellable smooth-flight step.
///
/// Native ABA/FK passes are indivisible. This bound limits the work between
/// cancellation polls; the additional momentum passes poll every 64 links.
pub const MAX_FREE_FLOATING_STEP_LINKS: usize = 4096;

/// Canonical base pose/twist and compact scalar joint state.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeFloatingMotionState {
    /// Unconstrained base pose and body-coordinate twist.
    pub base: FreeFloatingBaseState,
    /// Joint coordinates in model order (radians or metres).
    pub joint_position: Vec<f64>,
    /// Joint rates in the same compact order.
    pub joint_velocity: Vec<f64>,
}

/// Declared loads held over one step; no controller or contact law is inferred.
#[derive(Debug, Clone, Copy)]
pub struct FreeFloatingMotionLoads<'a> {
    /// Internal actuator efforts, conjugate to the scalar joint rates.
    pub generalized_force: &'a [f64],
    /// Per-link wrenches about each link origin, in that link's body frame.
    /// These include both force and its moment about the declared origin.
    pub external_body_wrench: &'a [Wrench],
    /// Uniform world-frame gravitational acceleration, applied to every mass.
    pub gravity_world: Vec3,
}

/// One fully accepted trajectory step and its physical impulse/work accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeFloatingMotionStep {
    /// Updated base pose/twist and joint coordinates/rates, ready for the next step.
    pub state_after: FreeFloatingMotionState,
    /// Whole-system COM before the step, in world coordinates.
    pub center_of_mass_world_before: Vec3,
    /// Whole-system COM after the step, in world coordinates.
    pub center_of_mass_world_after: Vec3,
    /// Physical momentum/kinetic energy of the supplied state.
    pub momentum_before: ArticulatedMomentum,
    /// Independently evaluated physical momentum/kinetic energy of `state_after`.
    pub momentum_after: ArticulatedMomentum,
    /// Midpoint external force impulse, including uniform gravity, in N s.
    pub external_linear_impulse_world: Vec3,
    /// External angular impulse about the world origin, in N m s.
    /// The force moment uses the actual midpoint-stage COM origin.
    pub external_angular_impulse_world: Vec3,
    /// Midpoint actuator effort/rate work, in joules.
    pub actuator_work: f64,
    /// Midpoint body-wrench work plus exact uniform-gravity displacement work.
    pub external_work: f64,
    /// Actual kinetic change minus actuator/external work, in joules.
    /// This is numerical integration error, not a physical dissipation channel.
    pub kinetic_work_residual: f64,
}

/// Atomic refusal of a continuous articulated step.
#[derive(Debug, Clone, PartialEq)]
pub enum FreeFloatingMotionError {
    /// Duration must be finite and strictly positive.
    InvalidDuration,
    /// The native non-interruptible tree passes exceed this lane's work bound.
    TooManyLinks {
        /// Supplied number of links.
        got: usize,
        /// Maximum admitted number of links.
        maximum: usize,
    },
    /// Cancellation was observed before publishing any updated state.
    Cancelled,
    /// Physical model, state, limits, or numerical admission refused.
    Dynamics(ArticulatedError),
}

impl fmt::Display for FreeFloatingMotionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDuration => {
                f.write_str("free-flight duration must be finite and positive")
            }
            Self::TooManyLinks { got, maximum } => {
                write!(f, "free-flight model has {got} links; maximum is {maximum}")
            }
            Self::Cancelled => f.write_str("free-flight step cancelled before state publication"),
            Self::Dynamics(source) => write!(f, "free-flight dynamics refused: {source}"),
        }
    }
}

impl std::error::Error for FreeFloatingMotionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(source) => Some(source),
            _ => None,
        }
    }
}

impl From<ArticulatedError> for FreeFloatingMotionError {
    fn from(value: ArticulatedError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<GaError> for FreeFloatingMotionError {
    fn from(value: GaError) -> Self {
        Self::Dynamics(ArticulatedError::Geometry(value))
    }
}

/// Advance a free-floating articulated trajectory by one second-order step.
///
/// Inputs are immutable. All limits, intermediate numerical admissions and the
/// final cancellation poll must succeed before a state is returned. Cancellation
/// is checked between native tree passes and every 64 links in momentum passes.
/// Work and storage are O(links); only fixed-size 6x6 inertias are solved.
///
/// The midpoint COM is `Xh = X0 + dt P0/(2 M)`, while the accepted COM uses
/// `X1 = X0 + dt (P0 + P1)/(2 M)`. Hence the orbital angular change is exactly
/// `dt Xh cross Fh` in real arithmetic. Advancing central angular momentum by
/// `dt tau_COM,h` consequently retains the same world-origin angular impulse
/// as the applied forces, even while joints move the base relative to the COM.
///
/// # Errors
/// Refuses invalid duration, an oversized or non-fixed-root model, invalid
/// states/loads, joint position/rate/effort limit violations at any stage,
/// nonphysical or ill-conditioned inertia, nonfinite results, or cancellation.
pub fn step_free_floating(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FreeFloatingMotionStep, FreeFloatingMotionError> {
    poll(&mut cancelled)?;
    if !duration.is_finite() || duration <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidDuration);
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        });
    }
    model.free_floating_complexity()?;
    validate_forward_dynamics_header(
        model,
        state.base.twist_body,
        &state.joint_position,
        &state.joint_velocity,
        loads.generalized_force,
        loads.gravity_world,
        loads.external_body_wrench,
    )?;
    validate_external_wrenches(loads.external_body_wrench)?;
    let initial_map = MomentumMap::new(model, state, &mut cancelled)?;
    let initial = initial_map.coordinates(state.base)?;
    poll(&mut cancelled)?;
    let momentum_before = free_floating_momentum(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(&mut cancelled)?;
    let initial_dynamics = free_floating_forward_dynamics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        loads.generalized_force,
        loads.gravity_world,
        loads.external_body_wrench,
    )?;
    let initial_load = initial_map.loads(state.base, loads, &mut cancelled)?;
    let half = 0.5 * duration;
    let midpoint_coordinates = MomentumCoordinates {
        center_world: initial.center_world + initial.linear_world.scale(half / initial_map.mass()),
        linear_world: initial.linear_world + initial_load.force_world.scale(half),
        angular_com_world: initial.angular_com_world + initial_load.torque_com_world.scale(half),
    };
    let mut midpoint = FreeFloatingMotionState {
        base: state.base,
        joint_position: advance(&state.joint_position, &state.joint_velocity, half),
        joint_velocity: advance(
            &state.joint_velocity,
            &initial_dynamics.generalized_acceleration,
            half,
        ),
    };
    let midpoint_rotation = rotate(
        state.base.world_from_base.rotation(),
        state.base.twist_body.angular,
        half,
    )?;
    let midpoint_map = MomentumMap::new(model, &midpoint, &mut cancelled)?;
    midpoint.base = midpoint_map.reconstruct(midpoint_rotation, midpoint_coordinates)?;
    poll(&mut cancelled)?;
    let midpoint_dynamics = free_floating_forward_dynamics(
        model,
        midpoint.base,
        &midpoint.joint_position,
        &midpoint.joint_velocity,
        loads.generalized_force,
        loads.gravity_world,
        loads.external_body_wrench,
    )?;
    let midpoint_load = midpoint_map.loads(midpoint.base, loads, &mut cancelled)?;
    let linear_impulse = midpoint_load.force_world.scale(duration);
    let angular_impulse = (midpoint_load.torque_com_world
        + midpoint_coordinates
            .center_world
            .cross(midpoint_load.force_world))
    .scale(duration);
    let final_coordinates = MomentumCoordinates {
        center_world: initial.center_world
            + (initial.linear_world + linear_impulse.scale(0.5))
                .scale(duration / initial_map.mass()),
        linear_world: initial.linear_world + linear_impulse,
        angular_com_world: initial.angular_com_world
            + midpoint_load.torque_com_world.scale(duration),
    };
    let mut state_after = FreeFloatingMotionState {
        base: state.base,
        joint_position: advance(&state.joint_position, &midpoint.joint_velocity, duration),
        joint_velocity: advance(
            &state.joint_velocity,
            &midpoint_dynamics.generalized_acceleration,
            duration,
        ),
    };
    let final_rotation = rotate(
        state.base.world_from_base.rotation(),
        midpoint.base.twist_body.angular,
        duration,
    )?;
    let final_map = MomentumMap::new(model, &state_after, &mut cancelled)?;
    state_after.base = final_map.reconstruct(final_rotation, final_coordinates)?;
    poll(&mut cancelled)?;
    let momentum_after = free_floating_momentum(
        model,
        state_after.base,
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    let actuator_work = duration
        * loads
            .generalized_force
            .iter()
            .zip(&midpoint.joint_velocity)
            .map(|(force, rate)| force * rate)
            .sum::<f64>();
    let gravity_work = initial_map.mass()
        * loads
            .gravity_world
            .dot(final_coordinates.center_world - initial.center_world);
    let external_work = duration * midpoint_load.body_wrench_power + gravity_work;
    let kinetic_work_residual = momentum_after.kinetic_energy
        - momentum_before.kinetic_energy
        - actuator_work
        - external_work;
    validate_vec3(linear_impulse, "free_floating_motion.linear_impulse")?;
    validate_vec3(angular_impulse, "free_floating_motion.angular_impulse")?;
    validate_slice(
        &[actuator_work, external_work, kinetic_work_residual],
        "free_floating_motion.work",
    )?;
    poll(&mut cancelled)?;
    Ok(FreeFloatingMotionStep {
        state_after,
        center_of_mass_world_before: initial.center_world,
        center_of_mass_world_after: final_coordinates.center_world,
        momentum_before,
        momentum_after,
        external_linear_impulse_world: linear_impulse,
        external_angular_impulse_world: angular_impulse,
        actuator_work,
        external_work,
        kinetic_work_residual,
    })
}

#[derive(Clone, Copy)]
struct MomentumCoordinates {
    center_world: Vec3,
    linear_world: Vec3,
    angular_com_world: Vec3,
}

/// Configuration-dependent momentum map at fixed relative joint rates.
struct MomentumMap {
    center_base: Vec3,
    locked: SpatialInertia,
    shape_linear_base: Vec3,
    shape_angular_com_base: Vec3,
    shape: Kinematics,
}

impl MomentumMap {
    fn new(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FreeFloatingMotionError> {
        poll(cancelled)?;
        // With a stationary identity base these are base-from-link poses and
        // link twists due solely to the declared relative joint rates.
        let shape = forward_kinematics(
            model,
            BaseState::stationary(Se3::identity()),
            &state.joint_position,
            &state.joint_velocity,
        )?;
        let mut mass = 0.0;
        let mut center_base = Vec3::default();
        for (i, link) in model.links.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let next_mass = mass + link.inertia.mass();
            let center = shape.world_from_link[i].transform_point(link.inertia.center_of_mass())?;
            center_base =
                center_base + (center - center_base).scale(link.inertia.mass() / next_mass);
            mass = next_mass;
        }
        validate_vec3(center_base, "free_floating_motion.center_base")?;
        let mut inertia = Mat3 { m: [0.0; 9] };
        let mut shape_linear_base = Vec3::default();
        let mut shape_angular_com_base = Vec3::default();
        for (i, link) in model.links.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let pose = shape.world_from_link[i];
            let rotation = pose.rotation();
            let r = pose.transform_point(link.inertia.center_of_mass())? - center_base;
            let r_array = [r.x, r.y, r.z];
            let r_squared = r.dot(r);
            let matrix = rotation.matrix();
            let rotated_inertia = matrix
                .compose(link.inertia.inertia_com())
                .compose(matrix.transpose());
            for row in 0..3 {
                for column in 0..3 {
                    let diagonal = if row == column { r_squared } else { 0.0 };
                    inertia.m[row * 3 + column] += rotated_inertia.m[row * 3 + column]
                        + link.inertia.mass() * (diagonal - r_array[row] * r_array[column]);
                }
            }
            let velocity = shape.body_twist[i];
            let linear = rotation
                .rotate(velocity.linear + velocity.angular.cross(link.inertia.center_of_mass()))?
                .scale(link.inertia.mass());
            shape_linear_base = shape_linear_base + linear;
            shape_angular_com_base = shape_angular_com_base
                + rotation.rotate(link.inertia.inertia_com().apply(velocity.angular))?
                + r.cross(linear);
        }
        validate_vec3(shape_linear_base, "free_floating_motion.shape_linear")?;
        validate_vec3(shape_angular_com_base, "free_floating_motion.shape_angular")?;
        let locked = SpatialInertia::new(mass, Vec3::default(), inertia)?;
        Ok(Self {
            center_base,
            locked,
            shape_linear_base,
            shape_angular_com_base,
            shape,
        })
    }

    fn mass(&self) -> f64 {
        self.locked.mass()
    }

    fn coordinates(
        &self,
        base: FreeFloatingBaseState,
    ) -> Result<MomentumCoordinates, FreeFloatingMotionError> {
        let rotation = base.world_from_base.rotation();
        let linear_base = (base.twist_body.linear
            + base.twist_body.angular.cross(self.center_base))
        .scale(self.mass())
            + self.shape_linear_base;
        let angular_com_base =
            self.locked.inertia_com().apply(base.twist_body.angular) + self.shape_angular_com_base;
        let result = MomentumCoordinates {
            center_world: base.world_from_base.transform_point(self.center_base)?,
            linear_world: rotation.rotate(linear_base)?,
            angular_com_world: rotation.rotate(angular_com_base)?,
        };
        validate_coordinates(result)?;
        Ok(result)
    }

    fn reconstruct(
        &self,
        rotation: So3,
        coordinates: MomentumCoordinates,
    ) -> Result<FreeFloatingBaseState, FreeFloatingMotionError> {
        validate_coordinates(coordinates)?;
        let inverse = rotation.inverse();
        let rhs = Wrench::new(
            inverse.rotate(coordinates.angular_com_world)? - self.shape_angular_com_base,
            inverse.rotate(coordinates.linear_world)? - self.shape_linear_base,
        );
        // Fixed qdot: this is the LOCKED composite inertia, not the ABA Schur
        // inertia that allows the joints to accelerate under a base wrench.
        let solved = solve_floating_base_system(&self.locked.matrix()?, rhs)?;
        let twist = Twist::new(
            solved.angular,
            solved.linear - solved.angular.cross(self.center_base),
        );
        validate_twist(twist, "free_floating_motion.reconstructed_twist")?;
        let position = coordinates.center_world - rotation.rotate(self.center_base)?;
        Ok(FreeFloatingBaseState::new(
            Se3::from_parts(rotation, position)?,
            twist,
        ))
    }

    fn loads(
        &self,
        base: FreeFloatingBaseState,
        loads: FreeFloatingMotionLoads<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StageLoad, FreeFloatingMotionError> {
        let rotation = base.world_from_base.rotation();
        let mut force_world = loads.gravity_world.scale(self.mass());
        let mut torque_com_world = Vec3::default();
        let mut body_wrench_power = 0.0;
        for (i, wrench) in loads.external_body_wrench.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let pose = self.shape.world_from_link[i];
            let force = rotation.rotate(pose.rotation().rotate(wrench.force)?)?;
            // Translate by the local COM arm before rotation, avoiding loss of
            // lever-arm precision when the entire system is far from the origin.
            let arm = rotation.rotate(pose.translation() - self.center_base)?;
            force_world = force_world + force;
            torque_com_world = torque_com_world
                + rotation.rotate(pose.rotation().rotate(wrench.torque)?)?
                + arm.cross(force);
            let velocity = pose
                .inverse()?
                .adjoint()
                .apply_twist(base.twist_body)
                .plus(self.shape.body_twist[i]);
            body_wrench_power += wrench.pairing(velocity);
        }
        validate_vec3(force_world, "free_floating_motion.force_world")?;
        validate_vec3(torque_com_world, "free_floating_motion.torque_com_world")?;
        validate_slice(
            &[body_wrench_power],
            "free_floating_motion.body_wrench_power",
        )?;
        Ok(StageLoad {
            force_world,
            torque_com_world,
            body_wrench_power,
        })
    }
}

struct StageLoad {
    force_world: Vec3,
    torque_com_world: Vec3,
    body_wrench_power: f64,
}

fn advance(value: &[f64], rate: &[f64], duration: f64) -> Vec<f64> {
    value
        .iter()
        .zip(rate)
        .map(|(value, rate)| value + duration * rate)
        .collect()
}

fn rotate(rotation: So3, angular: Vec3, duration: f64) -> Result<So3, GaError> {
    rotation.compose(So3::exp(So3Tangent::new(angular.scale(duration)))?)
}

fn validate_coordinates(value: MomentumCoordinates) -> Result<(), ArticulatedError> {
    validate_vec3(value.center_world, "free_floating_motion.center_world")?;
    validate_vec3(value.linear_world, "free_floating_motion.linear_world")?;
    validate_vec3(
        value.angular_com_world,
        "free_floating_motion.angular_com_world",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FreeFloatingMotionError> {
    if cancelled() {
        Err(FreeFloatingMotionError::Cancelled)
    } else {
        Ok(())
    }
}
