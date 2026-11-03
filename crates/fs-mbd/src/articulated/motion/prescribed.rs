//! Continuous articulated motion with a prescribed support trajectory.
//!
//! Native ABA supplies the two explicit midpoint stages; native FK and RNEA
//! retain the original link geometry, inertia and full support reaction. The
//! support is exerted on the tree, including root-actuator transmission. Its
//! work pairs with base motion, while actuator work pairs with relative joint
//! rates. No heavy surrogate base or force-duration impulse is introduced.
//!
//! The caller supplies smooth, mutually consistent base pose/twist/acceleration
//! samples. Their derivative consistency is not certified here. This explicit
//! second-order lane must resolve fast modes; it is not stiffly stable,
//! symplectic, energy-exact, or a contact/loop-constraint integrator. Reported
//! work defects remain numerical diagnostics, never heat or repaired state.

use super::super::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, Kinematics,
    forward_dynamics, forward_kinematics, inverse_dynamics, momentum_from_kinematics,
    validate_external_wrenches, validate_forward_dynamics_header, validate_slice, validate_twist,
    validate_wrench,
};
use super::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, MAX_FREE_FLOATING_STEP_LINKS, components,
    kinetic_link_increment, sum_products,
};
use fs_ga::{GaError, Se3, Vec3, Wrench};
use fs_math::dd::Dd;
use std::borrow::Cow;

/// Complete accepted state, including the support boundary for continuation.
#[derive(Clone, Debug, PartialEq)]
pub struct PrescribedBaseMotionState {
    /// Finite physical time in seconds.
    pub time_s: f64,
    /// Prescribed support pose, body twist and non-gravity spatial acceleration.
    pub base: BaseState,
    /// Original compact joint coordinates, in radians or metres.
    pub joint_position: Vec<f64>,
    /// Original compact joint rates, in radians/s or metres/s.
    pub joint_velocity: Vec<f64>,
}

/// Remaining samples of the prescribed base path over one motion interval.
///
/// The initial sample is the accepted state's base; it is never replaced.
/// Both samples must describe the same smooth path as that initial sample.
/// BaseState::acceleration_body uses the original Featherstone convention:
/// its linear part is R^T xddot - omega_body cross v_body, and its angular
/// part is R^T alpha_world. Gravity is supplied separately in the loads.
/// Derivative consistency between these samples is the caller's responsibility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedBasePath {
    /// Base boundary at state.time_s + duration / 2.
    pub midpoint: BaseState,
    /// Boundary at state.time_s + duration, retained for continuation.
    pub endpoint: BaseState,
}

impl PrescribedBasePath {
    /// Samples for a support held at one pose throughout the interval.
    ///
    /// The supplied state's base must also be stationary at this same pose.
    #[must_use]
    pub fn stationary(world_from_base: Se3) -> Self {
        let base = BaseState::stationary(world_from_base);
        Self {
            midpoint: base,
            endpoint: base,
        }
    }
}

/// Accepted endpoint and physical work accounts of the original body tree.
#[derive(Clone, Debug, PartialEq)]
pub struct PrescribedBaseMotionStep {
    /// Accepted state ready for the next smooth step or native impulse.
    pub state_after: PrescribedBaseMotionState,
    /// Native whole-tree momentum and kinetic storage before the step.
    pub momentum_before: ArticulatedMomentum,
    /// Native whole-tree momentum and kinetic storage at the endpoint.
    pub momentum_after: ArticulatedMomentum,
    /// Midpoint external body/world force impulse including uniform gravity, N s.
    /// Internal actuators, spatial connections and the support are excluded.
    pub external_linear_impulse_world: Vec3,
    /// The same source angular impulse about the world origin, in N m s.
    pub external_angular_impulse_world: Vec3,
    /// Midpoint native RNEA support force ON the tree, integrated in N s.
    pub support_linear_impulse_world: Vec3,
    /// Support angular impulse about the world origin, in N m s.
    /// Original root arms and central couples remain separate until reduction.
    pub support_angular_impulse_world: Vec3,
    /// Actual per-link linear-momentum increments, reduced before rounding, N s.
    /// This is not a difference of rounded whole-tree endpoint totals.
    pub linear_momentum_change_world: Vec3,
    /// Actual central-spin plus orbital increments about the world origin, N m s.
    pub angular_momentum_change_world: Vec3,
    /// Actual change minus independent midpoint external and support impulses.
    /// This measures numerical/path error; it is not repaired or assigned to heat.
    pub linear_momentum_residual_world: Vec3,
    /// Angular counterpart of linear_momentum_residual_world about the same origin.
    pub angular_momentum_residual_world: Vec3,
    /// Kinetic change from polarized link COM and central-spin increments, J.
    /// Absolute endpoint energy totals are retained separately as diagnostics.
    pub kinetic_energy_change: f64,
    /// Midpoint support wrench ON the tree about the base origin, in base axes.
    /// Includes both force and couple, including root-joint actuator transmission.
    pub support_wrench_body: Wrench,
    /// The same midpoint wrench about the world origin, in world axes.
    pub support_wrench_world: Wrench,
    /// Midpoint support wrench paired with the prescribed base twist, W.
    pub support_power: f64,
    /// Signed support work from midpoint quadrature, J.
    pub support_work: f64,
    /// Midpoint effort times relative joint rate, integrated over the step, J.
    pub actuator_work: f64,
    /// Midpoint external body-wrench work plus uniform-gravity work, J.
    /// Internal finite connections and the support are excluded.
    pub external_work: f64,
    /// Uniform gravity work from original relative geometry and accepted base
    /// translation, J. This is already included in external_work.
    pub gravity_work: f64,
    /// Midpoint internal spatial-connection work, J; zero for held loads alone.
    pub internal_body_work: f64,
    /// Kinetic change minus actuator, support, external and internal work, J.
    /// A diagnostic of integration/path error, not dissipation or a bound.
    pub kinetic_work_residual: f64,
}

/// Refusal publishes no partial trajectory or support state.
#[derive(Clone, Debug, PartialEq)]
pub enum PrescribedBaseMotionError {
    /// Invalid duration/clock, unrepresentable update, load or stage data.
    Invalid(&'static str),
    /// Bounds the original indivisible native tree passes between polls.
    TooManyLinks {
        /// Number of supplied links.
        got: usize,
        /// Maximum supported links.
        maximum: usize,
    },
    /// Cancellation observed before returning any accepted endpoint.
    Cancelled,
    /// Original model, configuration, rate, effort or geometry refusal.
    Dynamics(ArticulatedError),
}

impl core::fmt::Display for PrescribedBaseMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prescribed-base motion refused: {self:?}")
    }
}
impl std::error::Error for PrescribedBaseMotionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(source) => Some(source),
            _ => None,
        }
    }
}
impl From<ArticulatedError> for PrescribedBaseMotionError {
    fn from(value: ArticulatedError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<GaError> for PrescribedBaseMotionError {
    fn from(value: GaError) -> Self {
        Self::Dynamics(ArticulatedError::Geometry(value))
    }
}
impl From<FreeFloatingMotionError> for PrescribedBaseMotionError {
    fn from(value: FreeFloatingMotionError) -> Self {
        match value {
            FreeFloatingMotionError::InvalidDuration => Self::Invalid("duration"),
            FreeFloatingMotionError::InvalidLoad(what) => Self::Invalid(what),
            FreeFloatingMotionError::TooManyLinks { got, maximum } => {
                Self::TooManyLinks { got, maximum }
            }
            FreeFloatingMotionError::Cancelled => Self::Cancelled,
            FreeFloatingMotionError::Dynamics(source) => Self::Dynamics(source),
        }
    }
}

/// Advance the original articulated tree with a fixed or driven support.
///
/// The held load carrier is shared with free-flight mechanics: joint efforts
/// remain relative actuator efforts, link wrenches rotate with their original
/// frames, and gravity is a world acceleration. All root joint kinds supported
/// by native prescribed-base ABA are admitted, including a scalar root joint.
///
/// Inputs are immutable. Native configuration/rate/effort limits apply at both
/// dynamics stages and to the returned configuration/rates. Cancellation is
/// checked between native passes and every 64 links in additional reductions.
/// Work/storage are O(links); no generalized dense mass matrix is constructed.
///
/// Use the native prescribed-base impulse API for a zero-duration velocity
/// jump, then continue with its returned rates and the unchanged configuration.
/// The base path here is smooth; discontinuities are not resolved as forces.
pub fn step_prescribed_base(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBaseMotionStep, PrescribedBaseMotionError> {
    Ok(step_evaluated(
        model,
        state,
        path,
        loads,
        duration,
        &mut cancelled,
        |stage, kinematics, cancelled| held_stage(model, stage, loads, kinematics, cancelled),
        None,
    )?
    .step)
}

/// Advance with joint efforts evaluated at the two actual motion stages.
///
/// The callback is called in start/midpoint order with elapsed times 0 and
/// duration/2 and the complete stage state, including physical time and base
/// motion. Its output is reset to held efforts before each call. The original
/// effort limits apply to each result. The same midpoint result supplies work.
/// Constitutive storage/dissipation behind an effort law remains caller-owned.
///
/// The law must be smooth on the resolved interval. Cancellation is checked
/// before/after calls; callers own polling inside a long callback. Refusal
/// preserves all borrowed physical inputs, but cannot undo callback side effects.
pub fn step_prescribed_base_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBaseMotionStep, PrescribedBaseMotionError> {
    Ok(step_evaluated(
        model,
        state,
        path,
        loads,
        duration,
        &mut cancelled,
        |stage, kinematics, cancelled| held_stage(model, stage, loads, kinematics, cancelled),
        Some(&mut effort),
    )?
    .step)
}

pub(super) type EffortLaw<'a> = dyn FnMut(f64, &PrescribedBaseMotionState, &mut [f64]) -> Result<(), PrescribedBaseMotionError>
    + 'a;

pub(super) struct PrescribedStage<'a> {
    pub(super) body_wrench: Cow<'a, [Wrench]>,
    /// Original external sources before internal forces enter body_wrench.
    pub(super) external: SpatialSum,
    pub(super) external_power: f64,
    pub(super) internal_power: f64,
    pub(super) stored_energy: f64,
    pub(super) dissipation_power: f64,
}

pub(super) struct EvaluatedPrescribedStep {
    pub(super) step: PrescribedBaseMotionStep,
    pub(super) stored_energy_before: f64,
    pub(super) dissipated_energy: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn step_evaluated<'a, C: FnMut() -> bool>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: FreeFloatingMotionLoads<'a>,
    duration: f64,
    cancelled: &mut C,
    mut evaluate: impl FnMut(
        &PrescribedBaseMotionState,
        &Kinematics,
        &mut C,
    ) -> Result<PrescribedStage<'a>, PrescribedBaseMotionError>,
    mut effort: Option<&mut EffortLaw<'_>>,
) -> Result<EvaluatedPrescribedStep, PrescribedBaseMotionError> {
    poll(cancelled)?;
    let half = 0.5 * duration;
    let midpoint_time = state.time_s + half;
    let endpoint_time = state.time_s + duration;
    if !duration.is_finite()
        || duration <= 0.0
        || half <= 0.0
        || !state.time_s.is_finite()
        || !endpoint_time.is_finite()
        || !(midpoint_time > state.time_s && endpoint_time > midpoint_time)
    {
        return Err(PrescribedBaseMotionError::Invalid(
            "duration or representable clock",
        ));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(PrescribedBaseMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        });
    }
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
    for base in [state.base, path.midpoint, path.endpoint] {
        validate_twist(base.twist_body, "prescribed_motion.base_twist")?;
        validate_twist(
            base.acceleration_body,
            "prescribed_motion.base_acceleration",
        )?;
    }
    let before = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(cancelled)?;
    let momentum_before = momentum_from_kinematics(model, &before)?;
    poll(cancelled)?;
    let initial_load = evaluate(state, &before, cancelled)?;
    admit_stage(&initial_load)?;
    let mut stage_effort = Cow::Borrowed(loads.generalized_force);
    if let Some(law) = effort.as_deref_mut() {
        poll(cancelled)?;
        law(0.0, state, stage_effort.to_mut())?;
        poll(cancelled)?;
    }
    let initial_dynamics = forward_dynamics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &stage_effort,
        loads.gravity_world,
        &initial_load.body_wrench,
    )?;
    poll(cancelled)?;
    let midpoint = PrescribedBaseMotionState {
        time_s: midpoint_time,
        base: path.midpoint,
        joint_position: advance(&state.joint_position, &state.joint_velocity, half)?,
        joint_velocity: advance(
            &state.joint_velocity,
            &initial_dynamics.generalized_acceleration,
            half,
        )?,
    };
    let middle = forward_kinematics(
        model,
        midpoint.base,
        &midpoint.joint_position,
        &midpoint.joint_velocity,
    )?;
    poll(cancelled)?;
    let midpoint_load = evaluate(&midpoint, &middle, cancelled)?;
    admit_stage(&midpoint_load)?;
    if let Some(law) = effort.as_deref_mut() {
        poll(cancelled)?;
        stage_effort
            .to_mut()
            .copy_from_slice(loads.generalized_force);
        law(half, &midpoint, stage_effort.to_mut())?;
        poll(cancelled)?;
    }
    let midpoint_dynamics = forward_dynamics(
        model,
        midpoint.base,
        &midpoint.joint_position,
        &midpoint.joint_velocity,
        &stage_effort,
        loads.gravity_world,
        &midpoint_load.body_wrench,
    )?;
    poll(cancelled)?;
    let inverse = inverse_dynamics(
        model,
        midpoint.base,
        &midpoint.joint_position,
        &midpoint.joint_velocity,
        &midpoint_dynamics.generalized_acceleration,
        loads.gravity_world,
        &midpoint_load.body_wrench,
    )?;
    poll(cancelled)?;
    // Use the original root transform directly: subtracting two world poses
    // could erase its arm under a large common translation. The root wrench
    // already includes actuator transmission; adding it again double-counts work.
    let root = &model.links()[0];
    let root_coordinate = model
        .dof_index(0)
        .map_or(0.0, |i| midpoint.joint_position[i]);
    let base_from_root = root
        .parent_from_child_zero()
        .compose(root.joint().motion(root_coordinate)?)?;
    let support_wrench_body = inverse.body_wrench[0].transform_by(&base_from_root);
    let support_wrench_world = support_wrench_body.transform_by(&midpoint.base.world_from_base);
    validate_wrench(support_wrench_body, "prescribed_motion.support_body")?;
    validate_wrench(support_wrench_world, "prescribed_motion.support_world")?;
    // Integrate the actual midpoint reaction, independently of endpoint momentum.
    // Do not round a small root couple into a large world-origin force moment.
    let mut support = SpatialSum::default();
    support.add_world_wrench(
        middle.world_from_link[0]
            .rotation()
            .rotate(inverse.body_wrench[0].force)?,
        middle.world_from_link[0]
            .rotation()
            .rotate(inverse.body_wrench[0].torque)?,
        &[
            midpoint.base.world_from_base.translation(),
            midpoint
                .base
                .world_from_base
                .rotation()
                .rotate(base_from_root.translation())?,
        ],
    );
    let support = support.scaled(duration);
    let external = midpoint_load.external.scaled(duration);

    let state_after = PrescribedBaseMotionState {
        time_s: endpoint_time,
        base: path.endpoint,
        joint_position: advance(&state.joint_position, &midpoint.joint_velocity, duration)?,
        joint_velocity: advance(
            &state.joint_velocity,
            &midpoint_dynamics.generalized_acceleration,
            duration,
        )?,
    };
    let after = forward_kinematics(
        model,
        state_after.base,
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    poll(cancelled)?;
    let momentum_after = momentum_from_kinematics(model, &after)?;
    poll(cancelled)?;
    let (kinetic_change, gravity_work, momentum_change) = increments(
        model,
        &before,
        &after,
        state,
        &state_after,
        loads.gravity_world,
        cancelled,
    )?;
    let momentum_residual = momentum_change.subtract(external).subtract(support);
    let external = external.checked_wrench("prescribed_motion.external_impulse")?;
    let support = support.checked_wrench("prescribed_motion.support_impulse")?;
    let momentum_change = momentum_change.checked_wrench("prescribed_motion.momentum_change")?;
    let momentum_residual =
        momentum_residual.checked_wrench("prescribed_motion.momentum_residual")?;
    let actuator_power = sum_products(
        stage_effort.iter().copied(),
        midpoint.joint_velocity.iter().copied(),
    );
    let support_power = sum_products(
        support_wrench_body.to_array(),
        midpoint.base.twist_body.to_array(),
    );
    let dt = Dd::from_f64(duration);
    let actuator_work = dt * actuator_power;
    let support_work = dt * support_power;
    let external_work = dt * Dd::from_f64(midpoint_load.external_power) + gravity_work;
    let internal_body_work = dt * Dd::from_f64(midpoint_load.internal_power);
    let dissipated_energy = dt * Dd::from_f64(midpoint_load.dissipation_power);
    let residual =
        kinetic_change - actuator_work - support_work - external_work - internal_body_work;
    validate_slice(
        &[
            kinetic_change.to_f64(),
            gravity_work.to_f64(),
            actuator_work.to_f64(),
            support_power.to_f64(),
            support_work.to_f64(),
            external_work.to_f64(),
            internal_body_work.to_f64(),
            dissipated_energy.to_f64(),
            residual.to_f64(),
        ],
        "prescribed_motion.work",
    )?;
    poll(cancelled)?;
    Ok(EvaluatedPrescribedStep {
        stored_energy_before: initial_load.stored_energy,
        dissipated_energy: dissipated_energy.to_f64(),
        step: PrescribedBaseMotionStep {
            state_after,
            momentum_before,
            momentum_after,
            external_linear_impulse_world: external.force,
            external_angular_impulse_world: external.torque,
            support_linear_impulse_world: support.force,
            support_angular_impulse_world: support.torque,
            linear_momentum_change_world: momentum_change.force,
            angular_momentum_change_world: momentum_change.torque,
            linear_momentum_residual_world: momentum_residual.force,
            angular_momentum_residual_world: momentum_residual.torque,
            kinetic_energy_change: kinetic_change.to_f64(),
            support_wrench_body,
            support_wrench_world,
            support_power: support_power.to_f64(),
            support_work: support_work.to_f64(),
            actuator_work: actuator_work.to_f64(),
            external_work: external_work.to_f64(),
            gravity_work: gravity_work.to_f64(),
            internal_body_work: internal_body_work.to_f64(),
            kinetic_work_residual: residual.to_f64(),
        },
    })
}

fn held_stage<'a>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    loads: FreeFloatingMotionLoads<'a>,
    kinematics: &Kinematics,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PrescribedStage<'a>, PrescribedBaseMotionError> {
    let mut power = Dd::ZERO;
    for (i, (wrench, velocity)) in loads
        .external_body_wrench
        .iter()
        .zip(&kinematics.body_twist)
        .enumerate()
    {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        power = power + sum_products(wrench.to_array(), velocity.to_array());
    }
    // Native relative geometry retains local arms under a distant world origin.
    let relative = if state.base.world_from_base == Se3::identity() {
        Cow::Borrowed(kinematics)
    } else {
        poll(cancelled)?;
        let relative = forward_kinematics(
            model,
            BaseState::stationary(Se3::identity()),
            &state.joint_position,
            &state.joint_velocity,
        )?;
        poll(cancelled)?;
        Cow::Owned(relative)
    };
    let external = external_sources(
        model,
        state.base.world_from_base,
        kinematics,
        &relative,
        loads,
        cancelled,
    )?;
    Ok(PrescribedStage {
        body_wrench: Cow::Borrowed(loads.external_body_wrench),
        external,
        external_power: power.to_f64(),
        internal_power: 0.0,
        stored_energy: 0.0,
        dissipation_power: 0.0,
    })
}

fn admit_stage(stage: &PrescribedStage<'_>) -> Result<(), PrescribedBaseMotionError> {
    stage
        .external
        .checked_wrench("prescribed_motion.stage_sources")?;
    validate_slice(
        &[
            stage.external_power,
            stage.internal_power,
            stage.stored_energy,
            stage.dissipation_power,
        ],
        "prescribed_motion.stage_energy",
    )?;
    if stage.dissipation_power < 0.0 {
        return Err(PrescribedBaseMotionError::Invalid(
            "negative physical dissipation",
        ));
    }
    Ok(())
}

fn advance(
    values: &[f64],
    rates: &[f64],
    duration: f64,
) -> Result<Vec<f64>, PrescribedBaseMotionError> {
    values
        .iter()
        .zip(rates)
        .map(|(value, rate)| {
            let increment = duration * rate;
            let endpoint = value + increment;
            #[allow(clippy::float_cmp)]
            let swallowed = increment != 0.0 && endpoint == *value;
            if !increment.is_finite() || !endpoint.is_finite() || swallowed {
                return Err(PrescribedBaseMotionError::Invalid(
                    "unrepresentable joint update",
                ));
            }
            Ok(endpoint)
        })
        .collect()
}

// Spatial reductions keep local couples and each orbital contribution separate
// until after cancellation. A Wrench's f64 torque alone cannot retain those tails.
#[derive(Clone, Copy, Default)]
pub(super) struct SpatialSum {
    linear: [Dd; 3],
    angular: [Dd; 3],
}

impl SpatialSum {
    fn add_force(&mut self, force: [Dd; 3], position_parts: &[Vec3]) {
        add_components(&mut self.linear, force);
        for arm in position_parts {
            add_components(
                &mut self.angular,
                cross_dd(components(*arm).map(Dd::from_f64), force),
            );
        }
    }

    pub(super) fn add_world_wrench(&mut self, force: Vec3, couple: Vec3, position_parts: &[Vec3]) {
        self.add_force(components(force).map(Dd::from_f64), position_parts);
        add_components(&mut self.angular, components(couple).map(Dd::from_f64));
    }

    fn scaled(self, duration: f64) -> Self {
        let dt = Dd::from_f64(duration);
        Self {
            linear: self.linear.map(|value| dt * value),
            angular: self.angular.map(|value| dt * value),
        }
    }

    fn subtract(self, other: Self) -> Self {
        Self {
            linear: core::array::from_fn(|i| self.linear[i] - other.linear[i]),
            angular: core::array::from_fn(|i| self.angular[i] - other.angular[i]),
        }
    }

    fn checked_wrench(self, label: &'static str) -> Result<Wrench, PrescribedBaseMotionError> {
        for value in self.linear.iter().chain(&self.angular) {
            validate_slice(&[value.hi, value.lo], label)?;
        }
        let vector = |v: [Dd; 3]| Vec3::new(v[0].to_f64(), v[1].to_f64(), v[2].to_f64());
        Ok(Wrench::new(vector(self.angular), vector(self.linear)))
    }
}

fn add_components(sum: &mut [Dd; 3], value: [Dd; 3]) {
    for (sum, value) in sum.iter_mut().zip(value) {
        *sum = *sum + value;
    }
}

fn cross_dd(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Original held external loads and gravity, before any internal additions.
/// The caller already owns native world and base-relative stage geometry.
pub(super) fn external_sources(
    model: &ArticulatedModel,
    world_from_base: Se3,
    world: &Kinematics,
    relative: &Kinematics,
    loads: FreeFloatingMotionLoads<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialSum, PrescribedBaseMotionError> {
    let mut sources = SpatialSum::default();
    for (i, (link, wrench)) in model
        .links()
        .iter()
        .zip(loads.external_body_wrench)
        .enumerate()
    {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let rotation = world.world_from_link[i].rotation();
        let origin = world_from_base.translation();
        let arm = world_from_base
            .rotation()
            .rotate(relative.world_from_link[i].translation())?;
        sources.add_world_wrench(
            rotation.rotate(wrench.force)?,
            rotation.rotate(wrench.torque)?,
            &[origin, arm],
        );
        let inertia = link.inertia();
        let gravity_force =
            components(loads.gravity_world).map(|g| Dd::from_f64(inertia.mass()) * Dd::from_f64(g));
        sources.add_force(
            gravity_force,
            &[origin, arm, rotation.rotate(inertia.center_of_mass())?],
        );
    }
    Ok(sources)
}

fn increments(
    model: &ArticulatedModel,
    before: &Kinematics,
    after: &Kinematics,
    state: &PrescribedBaseMotionState,
    state_after: &PrescribedBaseMotionState,
    gravity: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Dd, Dd, SpatialSum), PrescribedBaseMotionError> {
    let mut kinetic = Dd::ZERO;
    let mut gravitational = Dd::ZERO;
    let mut momentum = SpatialSum::default();
    // Native relative FK keeps original local arms and their changes separate
    // from the base translation, including for orbital momentum without gravity.
    poll(cancelled)?;
    let first = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(cancelled)?;
    let last = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    poll(cancelled)?;
    for (i, link) in model.links().iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        let inertia = link.inertia();
        let center = inertia.center_of_mass();
        let old = before.body_twist[i];
        let new = after.body_twist[i];
        let pose0 = before.world_from_link[i];
        let pose1 = after.world_from_link[i];
        kinetic = kinetic + kinetic_link_increment(inertia, pose0, old, pose1, new)?;
        // Translational polarization in world axes prevents changing body frames
        // from manufacturing work from a large common inertial translation.
        let v0 = pose0
            .rotation()
            .rotate(old.linear + old.angular.cross(center))?;
        let v1 = pose1
            .rotation()
            .rotate(new.linear + new.angular.cross(center))?;
        let old_parts = [
            state.base.world_from_base.translation(),
            state
                .base
                .world_from_base
                .rotation()
                .rotate(first.world_from_link[i].translation())?,
            pose0.rotation().rotate(center)?,
        ];
        let new_parts = [
            state_after.base.world_from_base.translation(),
            state_after
                .base
                .world_from_base
                .rotation()
                .rotate(last.world_from_link[i].translation())?,
            pose1.rotation().rotate(center)?,
        ];
        let mut displacement = [Dd::ZERO; 3];
        let mut mean_position = [Dd::ZERO; 3];
        for (old_part, new_part) in old_parts.into_iter().zip(new_parts) {
            for ((axis, a), b) in components(old_part)
                .into_iter()
                .enumerate()
                .zip(components(new_part))
            {
                displacement[axis] = displacement[axis] + Dd::from_f64(b) - Dd::from_f64(a);
                mean_position[axis] =
                    mean_position[axis] + Dd::from_f64(0.5 * a) + Dd::from_f64(0.5 * b);
            }
        }
        let old_v = components(v0);
        let new_v = components(v1);
        let mass = Dd::from_f64(inertia.mass());
        let change = core::array::from_fn(|axis| {
            mass * (Dd::from_f64(new_v[axis]) - Dd::from_f64(old_v[axis]))
        });
        let mean_momentum = core::array::from_fn(|axis| {
            mass * (Dd::from_f64(0.5 * old_v[axis]) + Dd::from_f64(0.5 * new_v[axis]))
        });
        add_components(&mut momentum.linear, change);
        // Polarization: delta(r cross p) = mean(r) cross delta(p)
        //                                    + delta(r) cross mean(p).
        add_components(&mut momentum.angular, cross_dd(mean_position, change));
        add_components(&mut momentum.angular, cross_dd(displacement, mean_momentum));
        let spin0 = components(
            pose0
                .rotation()
                .rotate(inertia.inertia_com().apply(old.angular))?,
        );
        let spin1 = components(
            pose1
                .rotation()
                .rotate(inertia.inertia_com().apply(new.angular))?,
        );
        add_components(
            &mut momentum.angular,
            core::array::from_fn(|axis| Dd::from_f64(spin1[axis]) - Dd::from_f64(spin0[axis])),
        );
        for (g, dx) in components(gravity).into_iter().zip(displacement) {
            gravitational = gravitational + mass * Dd::from_f64(g) * dx;
        }
    }
    Ok((kinetic, gravitational, momentum))
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PrescribedBaseMotionError> {
    if cancelled() {
        Err(PrescribedBaseMotionError::Cancelled)
    } else {
        Ok(())
    }
}
