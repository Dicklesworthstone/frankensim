//! Continuous articulated motion with free or prescribed base boundaries.
//!
//! The free-floating entries below retain physical base recoil. The
//! [`prescribed`] module advances fixed or externally driven supports through
//! the same native articulated dynamics and recovers their force and work.
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

/// Continuous motion with a fixed or externally driven support boundary.
pub mod prescribed;

/// World-fixed loads and finite axial connections on original material points.
pub mod spatial_loads;

/// Momentum-preserving changes of full link inertia at fixed configuration.
pub mod reconfigure;

/// Momentum-conserving capture and release of rigid link payloads.
pub mod inventory;

use super::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, FreeFloatingBaseState,
    Kinematics, SpatialInertia, forward_kinematics, free_floating_forward_dynamics,
    free_floating_momentum, momentum_from_kinematics, solve_floating_base_system,
    validate_external_wrenches, validate_forward_dynamics_header, validate_slice, validate_twist,
    validate_vec3,
};
use core::fmt;
use fs_ga::{GaError, Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_math::dd::Dd;
use std::borrow::Cow;

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
    /// Midpoint body-wrench work plus uniform-gravity work along the computed
    /// COM increment, retaining displacements below the absolute pose's ulp.
    pub external_work: f64,
    /// Midpoint work of internal spatial connections, in joules; zero in the
    /// held-wrench API. Their storage and loss are owned by the connection lane.
    pub internal_body_work: f64,
    /// Actual kinetic change minus actuator/external/internal work, in joules.
    /// This is numerical integration error, not a physical dissipation channel.
    pub kinetic_work_residual: f64,
}

/// Atomic refusal of a continuous articulated step.
#[derive(Debug, Clone, PartialEq)]
pub enum FreeFloatingMotionError {
    /// Duration must be finite and strictly positive.
    InvalidDuration,
    /// A spatial load or connection violates its declared physical admission.
    InvalidLoad(&'static str),
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
            Self::InvalidLoad(what) => write!(f, "free-flight load refused: {what}"),
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
    Ok(step_evaluated(
        model,
        state,
        loads,
        duration,
        &mut cancelled,
        |map, base, cancelled| {
            Ok(EvaluatedStage {
                body_wrench: Cow::Borrowed(loads.external_body_wrench),
                load: map.loads(base, loads, cancelled)?,
            })
        },
        None,
    )?
    .step)
}

/// Advance with internal actuator efforts evaluated at both actual motion stages.
///
/// `effort` is called in start/midpoint order with elapsed times `0` and
/// `duration / 2`, the complete stage state, and one entry per original joint
/// coordinate. Before each call the output slice is reset to
/// `loads.generalized_force`; the callback may add to or replace those held
/// efforts. Revolute/helical efforts are torques in N m and prismatic efforts
/// are forces in N. An absolute clock can be captured by the callback.
///
/// Torque-speed curves, torsional spring/damper laws and smooth controllers
/// therefore use the actual stage rates and configuration. The original ABA
/// supplies equal/opposite base recoil, and the same midpoint effort supplies
/// `actuator_work`. Internal effort is never added to external spatial impulse.
/// A single O(joints) scratch buffer is reused across the two evaluations.
///
/// Second-order accuracy requires a smooth load law and a timestep resolving
/// its response; this is an explicit scheme, without stiff stability, contact
/// switching or a coupled electrical-state solve. The callback owns any
/// constitutive storage or dissipation interpretation of its effort. Native
/// inputs remain immutable on refusal, but callback-owned side effects cannot
/// be rolled back. Cancellation is checked before and after each callback;
/// callers own cancellation within a long-running callback.
///
/// # Errors
/// The original motion admissions, including finite efforts and native effort
/// limits at both stages, cancellation, or an error returned by the callback.
pub fn step_free_floating_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FreeFloatingMotionStep, FreeFloatingMotionError> {
    Ok(step_evaluated(
        model,
        state,
        loads,
        duration,
        &mut cancelled,
        |map, base, cancelled| {
            Ok(EvaluatedStage {
                body_wrench: Cow::Borrowed(loads.external_body_wrench),
                load: map.loads(base, loads, cancelled)?,
            })
        },
        Some(&mut effort),
    )?
    .step)
}

type EffortLaw<'a> = dyn FnMut(f64, &FreeFloatingMotionState, &mut [f64]) -> Result<(), FreeFloatingMotionError>
    + 'a;

fn step_evaluated<'a, C: FnMut() -> bool>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'a>,
    duration: f64,
    cancelled: &mut C,
    mut evaluate: impl FnMut(
        &MomentumMap,
        FreeFloatingBaseState,
        &mut C,
    ) -> Result<EvaluatedStage<'a>, FreeFloatingMotionError>,
    mut effort: Option<&mut EffortLaw<'_>>,
) -> Result<EvaluatedStep, FreeFloatingMotionError> {
    poll(cancelled)?;
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
    let initial_map = MomentumMap::new(model, state, cancelled)?;
    let initial = initial_map.coordinates(state.base)?;
    poll(cancelled)?;
    let kinematics_before = forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    let momentum_before = momentum_from_kinematics(model, &kinematics_before)?;
    poll(cancelled)?;
    let initial_stage = evaluate(&initial_map, state.base, cancelled)?;
    let mut stage_effort = Cow::Borrowed(loads.generalized_force);
    if let Some(law) = effort.as_deref_mut() {
        poll(cancelled)?;
        law(0.0, state, stage_effort.to_mut())?;
        poll(cancelled)?;
    }
    let initial_dynamics = free_floating_forward_dynamics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &stage_effort,
        loads.gravity_world,
        &initial_stage.body_wrench,
    )?;
    let initial_load = initial_stage.load;
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
    let midpoint_map = MomentumMap::new(model, &midpoint, cancelled)?;
    midpoint.base = midpoint_map.reconstruct(midpoint_rotation, midpoint_coordinates)?;
    poll(cancelled)?;
    let midpoint_stage = evaluate(&midpoint_map, midpoint.base, cancelled)?;
    if let Some(law) = effort.as_deref_mut() {
        poll(cancelled)?;
        stage_effort
            .to_mut()
            .copy_from_slice(loads.generalized_force);
        law(half, &midpoint, stage_effort.to_mut())?;
        poll(cancelled)?;
    }
    let midpoint_dynamics = free_floating_forward_dynamics(
        model,
        midpoint.base,
        &midpoint.joint_position,
        &midpoint.joint_velocity,
        &stage_effort,
        loads.gravity_world,
        &midpoint_stage.body_wrench,
    )?;
    let midpoint_load = midpoint_stage.load;
    let linear_impulse = midpoint_load.force_world.scale(duration);
    let angular_impulse = (midpoint_load.torque_com_world
        + midpoint_coordinates
            .center_world
            .cross(midpoint_load.force_world))
    .scale(duration);
    let center_displacement =
        (initial.linear_world + linear_impulse.scale(0.5)).scale(duration / initial_map.mass());
    let final_coordinates = MomentumCoordinates {
        center_world: initial.center_world + center_displacement,
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
    let final_map = MomentumMap::new(model, &state_after, cancelled)?;
    state_after.base = final_map.reconstruct(final_rotation, final_coordinates)?;
    poll(cancelled)?;
    let kinematics_after = forward_kinematics(
        model,
        BaseState::prescribed(
            state_after.base.world_from_base,
            state_after.base.twist_body,
            Twist::zero(),
        ),
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    let momentum_after = momentum_from_kinematics(model, &kinematics_after)?;
    let kinetic_change =
        kinetic_increment(model, &kinematics_before, &kinematics_after, cancelled)?;
    let actuator_work_dd = Dd::from_f64(duration)
        * sum_products(
            stage_effort.iter().copied(),
            midpoint.joint_velocity.iter().copied(),
        );
    let gravity_work = Dd::from_f64(initial_map.mass())
        * sum_products(
            components(loads.gravity_world),
            components(center_displacement),
        );
    let external_work_dd =
        Dd::from_f64(duration) * Dd::from_f64(midpoint_load.body_wrench_power) + gravity_work;
    let internal_body_work_dd =
        Dd::from_f64(duration) * Dd::from_f64(midpoint_load.internal_body_power);
    let actuator_work = actuator_work_dd.to_f64();
    let external_work = external_work_dd.to_f64();
    let internal_body_work = internal_body_work_dd.to_f64();
    let dissipated_energy = duration * midpoint_load.dissipation_power;
    let kinetic_work_residual =
        (kinetic_change - actuator_work_dd - external_work_dd - internal_body_work_dd).to_f64();
    validate_vec3(linear_impulse, "free_floating_motion.linear_impulse")?;
    validate_vec3(angular_impulse, "free_floating_motion.angular_impulse")?;
    validate_slice(
        &[
            actuator_work,
            external_work,
            internal_body_work,
            dissipated_energy,
            kinetic_work_residual,
        ],
        "free_floating_motion.work",
    )?;
    poll(cancelled)?;
    Ok(EvaluatedStep {
        stored_energy_before: initial_load.stored_energy,
        dissipated_energy,
        step: FreeFloatingMotionStep {
            state_after,
            center_of_mass_world_before: initial.center_world,
            center_of_mass_world_after: final_coordinates.center_world,
            momentum_before,
            momentum_after,
            external_linear_impulse_world: linear_impulse,
            external_angular_impulse_world: angular_impulse,
            actuator_work,
            external_work,
            internal_body_work,
            kinetic_work_residual,
        },
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
    shape_linear_base: VectorSum,
    shape_angular_com_base: VectorSum,
    shape: Kinematics,
}

impl MomentumMap {
    fn new(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FreeFloatingMotionError> {
        Self::from_joint_state(
            model,
            &state.joint_position,
            &state.joint_velocity,
            cancelled,
        )
    }

    // Shared material geometry for either free or prescribed base motion.
    // No floating-base reconstruction or fixed-root admission occurs here.
    fn from_joint_state(
        model: &ArticulatedModel,
        joint_position: &[f64],
        joint_velocity: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FreeFloatingMotionError> {
        poll(cancelled)?;
        // With a stationary identity base these are base-from-link poses and
        // link twists due solely to the declared relative joint rates.
        let shape = forward_kinematics(
            model,
            BaseState::stationary(Se3::identity()),
            joint_position,
            joint_velocity,
        )?;
        let mut mass = Dd::ZERO;
        let mut center = VectorSum::default();
        for (i, link) in model.links.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            mass = mass + Dd::from_f64(link.inertia.mass());
            let point = shape.world_from_link[i].transform_point(link.inertia.center_of_mass())?;
            let weight = Dd::from_f64(link.inertia.mass()) / mass;
            for (mean, value) in center.0.iter_mut().zip(components(point)) {
                *mean = *mean + (Dd::from_f64(value) - *mean) * weight;
            }
        }
        let center_base = center.finish();
        validate_vec3(center_base, "free_floating_motion.center_base")?;
        let mut inertia = [Dd::ZERO; 9];
        let mut shape_linear_base = VectorSum::default();
        let mut shape_angular_com_base = VectorSum::default();
        for (i, link) in model.links.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let pose = shape.world_from_link[i];
            let rotation = pose.rotation();
            let r = pose.transform_point(link.inertia.center_of_mass())? - center_base;
            let r_array = [r.x, r.y, r.z];
            let matrix = rotation.matrix();
            let rotated_inertia = matrix
                .compose(link.inertia.inertia_com())
                .compose(matrix.transpose());
            for row in 0..3 {
                for column in 0..3 {
                    // The diagonal is the sum of the two perpendicular
                    // squared arms. Subtracting r[row]^2 from r.dot(r) loses
                    // transverse orbital inertia on elongated assemblies.
                    let orbital = if row == column {
                        let a = r_array[(row + 1) % 3];
                        let b = r_array[(row + 2) % 3];
                        a * a + b * b
                    } else {
                        -r_array[row] * r_array[column]
                    };
                    let entry = &mut inertia[row * 3 + column];
                    *entry = *entry
                        + Dd::from_f64(rotated_inertia.m[row * 3 + column])
                        + Dd::from_f64(link.inertia.mass()) * Dd::from_f64(orbital);
                }
            }
            let velocity = shape.body_twist[i];
            let linear = rotation
                .rotate(velocity.linear + velocity.angular.cross(link.inertia.center_of_mass()))?
                .scale(link.inertia.mass());
            shape_linear_base.add(linear);
            shape_angular_com_base
                .add(rotation.rotate(link.inertia.inertia_com().apply(velocity.angular))?);
            shape_angular_com_base.add_cross(r, linear);
        }
        validate_vec3(
            shape_linear_base.finish(),
            "free_floating_motion.shape_linear",
        )?;
        validate_vec3(
            shape_angular_com_base.finish(),
            "free_floating_motion.shape_angular",
        )?;
        let locked = SpatialInertia::new(
            mass.to_f64(),
            Vec3::default(),
            Mat3 {
                m: inertia.map(Dd::to_f64),
            },
        )?;
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
        // Preserve the residual of a moving base and counter-moving shape.
        // Rounding shape momentum before the combination can erase real recoil.
        let mut velocity = VectorSum::from(base.twist_body.linear);
        velocity.add_cross(base.twist_body.angular, self.center_base);
        let mut linear_base = self.shape_linear_base;
        for (sum, value) in linear_base.0.iter_mut().zip(velocity.0) {
            *sum = *sum + value * Dd::from_f64(self.mass());
        }
        let mut angular_com_base = self.shape_angular_com_base;
        angular_com_base.add_matrix_product(self.locked.inertia_com(), base.twist_body.angular);
        let result = MomentumCoordinates {
            center_world: base.world_from_base.transform_point(self.center_base)?,
            linear_world: rotation.rotate(linear_base.finish())?,
            angular_com_world: rotation.rotate(angular_com_base.finish())?,
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
            self.shape_angular_com_base
                .subtract_from(inverse.rotate(coordinates.angular_com_world)?),
            self.shape_linear_base
                .subtract_from(inverse.rotate(coordinates.linear_world)?),
        );
        // Fixed qdot: this is the LOCKED composite inertia, not the ABA Schur
        // inertia that allows the joints to accelerate under a base wrench.
        let solved = solve_floating_base_system(&self.locked.matrix()?, rhs)?;
        // At the system COM the locked inertia is block diagonal. Recover
        // translation by the physical scalar mass, without the rotational
        // block's scaling and Cholesky roundoff contaminating a large boost.
        let linear = Vec3::new(
            rhs.force.x / self.mass(),
            rhs.force.y / self.mass(),
            rhs.force.z / self.mass(),
        );
        let twist = Twist::new(
            solved.angular,
            linear - solved.angular.cross(self.center_base),
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
        Ok(self.loads_accumulated(base, loads, cancelled)?.finish())
    }

    fn loads_accumulated(
        &self,
        base: FreeFloatingBaseState,
        loads: FreeFloatingMotionLoads<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StageLoadAccumulator, FreeFloatingMotionError> {
        let rotation = base.world_from_base.rotation();
        let mut accumulated = StageLoadAccumulator::default();
        accumulated.add_force(loads.gravity_world.scale(self.mass()));
        for (i, wrench) in loads.external_body_wrench.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let pose = self.shape.world_from_link[i];
            let force = rotation.rotate(pose.rotation().rotate(wrench.force)?)?;
            // Translate by the local COM arm before rotation, avoiding loss of
            // lever-arm precision when the entire system is far from the origin.
            let arm = rotation.rotate(pose.translation() - self.center_base)?;
            accumulated.add_force(force);
            accumulated.add_torque(rotation.rotate(pose.rotation().rotate(wrench.torque)?)?);
            accumulated.add_torque_cross(arm, force);
            let velocity = pose
                .inverse()?
                .adjoint()
                .apply_twist(base.twist_body)
                .plus(self.shape.body_twist[i]);
            accumulated.add_external_pairing(*wrench, velocity);
        }
        let total = accumulated.finish();
        validate_vec3(total.force_world, "free_floating_motion.force_world")?;
        validate_vec3(
            total.torque_com_world,
            "free_floating_motion.torque_com_world",
        )?;
        validate_slice(
            &[total.body_wrench_power],
            "free_floating_motion.body_wrench_power",
        )?;
        Ok(accumulated)
    }
}

/// Fixed-size compensated reductions retain physically small resultant momentum
/// until all opposing contributions have arrived. No per-link heap is required.
#[derive(Clone, Copy, Default)]
struct VectorSum([Dd; 3]);

impl From<Vec3> for VectorSum {
    fn from(value: Vec3) -> Self {
        Self(components(value).map(Dd::from_f64))
    }
}

impl VectorSum {
    fn add(&mut self, value: Vec3) {
        for (sum, value) in self.0.iter_mut().zip(components(value)) {
            *sum = *sum + Dd::from_f64(value);
        }
    }

    fn add_cross(&mut self, arm: Vec3, force: Vec3) {
        let arm = components(arm).map(Dd::from_f64);
        let force = components(force).map(Dd::from_f64);
        for axis in 0..3 {
            let a = (axis + 1) % 3;
            let b = (axis + 2) % 3;
            self.0[axis] = self.0[axis] + arm[a] * force[b] - arm[b] * force[a];
        }
    }

    fn add_matrix_product(&mut self, matrix: Mat3, value: Vec3) {
        let value = components(value).map(Dd::from_f64);
        for row in 0..3 {
            for (column, value) in value.iter().enumerate() {
                self.0[row] = self.0[row] + Dd::from_f64(matrix.m[3 * row + column]) * *value;
            }
        }
    }

    fn subtract_from(self, value: Vec3) -> Vec3 {
        let mut difference = Self::from(value);
        for (component, value) in difference.0.iter_mut().zip(self.0) {
            *component = *component - value;
        }
        difference.finish()
    }

    fn finish(self) -> Vec3 {
        Vec3::new(self.0[0].to_f64(), self.0[1].to_f64(), self.0[2].to_f64())
    }
}

/// Kept through held, world-point and attachment evaluation before rounding.
/// The finalized StageLoad remains the same small native stage interface.
#[derive(Clone, Copy, Default)]
struct StageLoadAccumulator {
    force_world: VectorSum,
    torque_com_world: VectorSum,
    body_wrench_power: Dd,
    internal_body_power: Dd,
    dissipation_power: Dd,
    stored_energy: Dd,
}

impl StageLoadAccumulator {
    fn add_force(&mut self, force: Vec3) {
        self.force_world.add(force);
    }

    fn add_torque(&mut self, torque: Vec3) {
        self.torque_com_world.add(torque);
    }

    fn add_torque_cross(&mut self, arm: Vec3, force: Vec3) {
        self.torque_com_world.add_cross(arm, force);
    }

    fn add_external_pairing(&mut self, wrench: Wrench, velocity: Twist) {
        self.body_wrench_power =
            self.body_wrench_power + sum_products(wrench.to_array(), velocity.to_array());
    }

    fn add_internal_power(&mut self, power: f64) {
        self.internal_body_power = self.internal_body_power + Dd::from_f64(power);
    }

    fn add_dissipation(&mut self, power: f64) {
        self.dissipation_power = self.dissipation_power + Dd::from_f64(power);
    }

    fn add_storage(&mut self, energy: f64) {
        self.stored_energy = self.stored_energy + Dd::from_f64(energy);
    }

    fn finish(self) -> StageLoad {
        StageLoad {
            force_world: self.force_world.finish(),
            torque_com_world: self.torque_com_world.finish(),
            body_wrench_power: self.body_wrench_power.to_f64(),
            internal_body_power: self.internal_body_power.to_f64(),
            dissipation_power: self.dissipation_power.to_f64(),
            stored_energy: self.stored_energy.to_f64(),
        }
    }
}

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn sum_products(a: impl IntoIterator<Item = f64>, b: impl IntoIterator<Item = f64>) -> Dd {
    a.into_iter().zip(b).fold(Dd::ZERO, |sum, (a, b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

/// Polarization avoids subtracting two large total kinetic energies. COM
/// translation uses world axes; central spin uses the fixed material inertia.
fn kinetic_link_increment(
    inertia: SpatialInertia,
    pose_before: Se3,
    before: Twist,
    pose_after: Se3,
    after: Twist,
) -> Result<Dd, GaError> {
    let center = inertia.center_of_mass();
    let v0 = pose_before
        .rotation()
        .rotate(before.linear + before.angular.cross(center))?;
    let v1 = pose_after
        .rotation()
        .rotate(after.linear + after.angular.cross(center))?;
    let mut spin_change = VectorSum::default();
    spin_change.add_matrix_product(inertia.inertia_com(), after.angular - before.angular);
    let mut energy = Dd::ZERO;
    for ((old, new), delta) in components(v0)
        .into_iter()
        .zip(components(v1))
        .zip(components(v1 - v0))
    {
        let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
        energy = energy + Dd::from_f64(inertia.mass()) * Dd::from_f64(delta) * mean;
    }
    for ((old, new), impulse) in components(before.angular)
        .into_iter()
        .zip(components(after.angular))
        .zip(spin_change.0)
    {
        let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
        energy = energy + impulse * mean;
    }
    Ok(energy)
}

fn kinetic_increment(
    model: &ArticulatedModel,
    before: &Kinematics,
    after: &Kinematics,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, FreeFloatingMotionError> {
    let mut energy = Dd::ZERO;
    for (i, link) in model.links.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        energy = energy
            + kinetic_link_increment(
                link.inertia,
                before.world_from_link[i],
                before.body_twist[i],
                after.world_from_link[i],
                after.body_twist[i],
            )?;
    }
    Ok(energy)
}

#[derive(Clone, Copy)]
struct StageLoad {
    force_world: Vec3,
    torque_com_world: Vec3,
    body_wrench_power: f64,
    internal_body_power: f64,
    dissipation_power: f64,
    stored_energy: f64,
}

struct EvaluatedStage<'a> {
    body_wrench: Cow<'a, [Wrench]>,
    load: StageLoad,
}

struct EvaluatedStep {
    step: FreeFloatingMotionStep,
    stored_energy_before: f64,
    dissipated_energy: f64,
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
