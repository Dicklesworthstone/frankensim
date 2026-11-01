//! Live full-inertia changes for a closed, freely moving articulated inventory.
//!
//! At unchanged configuration, preserve every base and joint momentum:
//! `M_next(q) u_after = M_previous(q) u_before`. The original articulated-body
//! inverse-mass recursion supplies the new rates in linear work and storage;
//! no dense generalized matrix or physical external impulse is introduced.
//!
//! Only centroidal inertia tensors may change. Every link retains its mass,
//! COM, body frame, joint and attachment, and its place in the tree. This is a
//! prescribed parameter change, not a model of the internal shape motion that
//! produced it. The caller owns that closed-inventory interpretation and its
//! energy source. Mass exchange, COM/attachment motion, shape velocities,
//! contacts, constraints and fluid storage require their own physical laws.

/// Continuous prescribed inertia histories through the native smooth dynamics.
#[path = "reconfigure/step.rs"]
pub mod step;

/// Constant affine joint and gear constraints during tensor reconfiguration.
#[path = "reconfigure/constraints.rs"]
pub mod constraints;

use super::{FreeFloatingMotionState, MAX_FREE_FLOATING_STEP_LINKS};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, FreeFloatingBaseState,
    JointKind, Kinematics, ensure_finite_scalar, forward_kinematics, free_floating_impulse_delta,
    inverse_dynamics, momentum_from_kinematics, validate_slice, validate_vec3,
};
use core::fmt;
use fs_ga::{GaError, Twist, Vec3, Wrench};

/// Physical momentum acceptance and an absolute parameter-work budget.
#[derive(Debug, Clone, Copy)]
pub struct ReconfigurationControls<'a> {
    /// Relative componentwise momentum tolerance, finite and in `[0, 1e-8]`.
    pub relative_momentum_tolerance: f64,
    /// Positive absolute tolerance for base and world angular momenta, N m s.
    pub angular_momentum_tolerance_nms: f64,
    /// Positive absolute tolerance for base and world linear momenta, N s.
    pub linear_momentum_tolerance_ns: f64,
    /// Positive absolute tolerance per joint, in compact model order. Units
    /// are N m s for revolute/helical joints and N s for prismatic joints.
    pub joint_momentum_tolerances: &'a [f64],
    /// Finite nonnegative cap on `abs(parameter_work_j)`. Both input work and
    /// released mechanical energy must fit the declared budget.
    pub maximum_parameter_work_j: f64,
}

/// One accepted tensor change with independently reconstructed momenta.
#[derive(Debug, Clone, PartialEq)]
pub struct InertiaReconfiguration {
    /// Original pose and coordinates, with the momentum-preserving new rates.
    pub state_after: FreeFloatingMotionState,
    /// Physical world momentum and lab-frame kinetic energy before the change.
    pub momentum_before: ArticulatedMomentum,
    /// Physical world momentum and lab-frame kinetic energy of the returned state.
    pub momentum_after: ArticulatedMomentum,
    /// Native `[base angular, base linear, joint]` conjugate momenta before.
    pub generalized_momentum_before: Vec<f64>,
    /// The same mass action independently evaluated on the returned state.
    pub generalized_momentum_after: Vec<f64>,
    /// Signed mechanical kinetic change in joules; positive work enters the
    /// moving tree. Evaluated in a common translating frame with the actual
    /// linear-momentum residual correction, avoiding cancellation of unchanged
    /// bulk kinetic energy. No internal deformation energy is inferred.
    pub parameter_work_j: f64,
    /// Largest componentwise residual across all 6+n generalized momenta and
    /// independently accumulated world angular/linear momenta. At most one.
    pub maximum_normalized_momentum_residual: f64,
}

/// Atomic refusal: no model, state, or work budget is mutated.
#[derive(Debug, Clone, PartialEq)]
pub enum ReconfigurationError {
    /// Incompatible inventories/kinematics or invalid acceptance controls.
    Invalid(&'static str),
    /// Native geometry, physical inertia, state, limit, or inverse-mass refusal.
    Dynamics(ArticulatedError),
    /// Cancellation observed before publishing a result.
    Cancelled,
    /// Independently reconstructed momentum exceeds its declared tolerance.
    MomentumResidual {
        /// Maximum componentwise residual divided by its allowed error.
        maximum_normalized: f64,
    },
    /// Signed parameter work exceeds the absolute admitted budget.
    ParameterWork {
        /// Actual signed mechanical work, joules.
        work_j: f64,
        /// Supplied absolute budget, joules.
        maximum_j: f64,
    },
    /// Indivisible native tree passes exceed the cancellable lane's work bound.
    TooManyLinks {
        /// Supplied link count.
        got: usize,
        /// Maximum admitted link count.
        maximum: usize,
    },
}

impl fmt::Display for ReconfigurationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "inertia reconfiguration refused: {what}"),
            Self::Dynamics(source) => write!(f, "inertia reconfiguration dynamics: {source}"),
            Self::Cancelled => f.write_str("inertia reconfiguration cancelled before publication"),
            Self::MomentumResidual { maximum_normalized } => write!(
                f,
                "inertia reconfiguration normalized momentum residual {maximum_normalized} exceeds one"
            ),
            Self::ParameterWork { work_j, maximum_j } => write!(
                f,
                "inertia reconfiguration work {work_j} J exceeds absolute budget {maximum_j} J"
            ),
            Self::TooManyLinks { got, maximum } => write!(
                f,
                "inertia reconfiguration has {got} links; maximum is {maximum}"
            ),
        }
    }
}

impl std::error::Error for ReconfigurationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(source) => Some(source),
            _ => None,
        }
    }
}

impl From<ArticulatedError> for ReconfigurationError {
    fn from(value: ArticulatedError) -> Self {
        Self::Dynamics(value)
    }
}

impl From<GaError> for ReconfigurationError {
    fn from(value: GaError) -> Self {
        Self::Dynamics(ArticulatedError::Geometry(value))
    }
}

/// Transfer a live free-floating state to changed full centroidal inertias.
///
/// Models must have identical links, transforms, joints and limits, with
/// bit-identical per-link mass and COM. Only the validated centroidal tensor
/// may change. The whole 6+n momentum vector is retained, including the joint
/// momenta responsible for base recoil. Native state and speed limits apply
/// to the returned rates; no rate clipping or artificial impulse is allowed.
///
/// An identical model returns an exact clone of the supplied state after
/// admission, with zero work. All other outcomes are atomic refusals or a
/// complete new state. Cancellation is polled between native tree passes and
/// every 64 links in additional work; native passes admit at most 4096 links.
pub fn transfer_inertia(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    controls: ReconfigurationControls<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InertiaReconfiguration, ReconfigurationError> {
    poll(&mut cancelled)?;
    validate_compatible_models(previous, next)?;
    validate_controls(previous, controls)?;
    let before = state_kinematics(previous, state)?;
    poll(&mut cancelled)?;
    let momentum_before = momentum_from_kinematics(previous, &before)?;
    poll(&mut cancelled)?;
    let generalized_momentum_before = generalized_momentum(previous, state)?;
    poll(&mut cancelled)?;
    if previous == next {
        return Ok(InertiaReconfiguration {
            state_after: state.clone(),
            momentum_before,
            momentum_after: momentum_before,
            generalized_momentum_after: generalized_momentum_before.clone(),
            generalized_momentum_before,
            parameter_work_j: 0.0,
            maximum_normalized_momentum_residual: 0.0,
        });
    }

    let (state_after, boost_world) =
        raw_momentum_candidate(previous, next, state, momentum_before, &mut cancelled)?;
    let after = state_kinematics(next, &state_after)?;
    poll(&mut cancelled)?;
    let momentum_after = momentum_from_kinematics(next, &after)?;
    poll(&mut cancelled)?;
    let generalized_momentum_after = generalized_momentum(next, &state_after)?;
    poll(&mut cancelled)?;
    let maximum_normalized_momentum_residual = audit_momentum(
        &generalized_momentum_before,
        &generalized_momentum_after,
        momentum_before,
        momentum_after,
        controls,
    )?;
    let parameter_work_j =
        parameter_work(previous, next, &before, &after, boost_world, &mut cancelled)?;
    if parameter_work_j.abs() > controls.maximum_parameter_work_j {
        return Err(ReconfigurationError::ParameterWork {
            work_j: parameter_work_j,
            maximum_j: controls.maximum_parameter_work_j,
        });
    }
    poll(&mut cancelled)?;
    Ok(InertiaReconfiguration {
        state_after,
        momentum_before,
        momentum_after,
        generalized_momentum_before,
        generalized_momentum_after,
        parameter_work_j,
        maximum_normalized_momentum_residual,
    })
}

// Internal candidate only: the constrained owner must add its simultaneous
// reactions before applying endpoint speed, momentum and work admission.
fn raw_momentum_candidate(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    momentum_before: ArticulatedMomentum,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(FreeFloatingMotionState, Vec3), ReconfigurationError> {
    // A uniform translation has the same generalized momentum in both models:
    // only the centroidal tensors differ. Remove it before the native solve
    // so a large Galilean boost cannot swamp the changing rotational momentum.
    let mut total_mass = 0.0;
    for (index, link) in previous.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        total_mass += link.inertia().mass();
        ensure_finite_scalar(total_mass, "reconfigure.total_mass", index)?;
    }
    let linear = momentum_before.linear_momentum_world;
    let boost_world = Vec3::new(
        linear.x / total_mass,
        linear.y / total_mass,
        linear.z / total_mass,
    );
    validate_vec3(boost_world, "reconfigure.boost_world")?;
    let boost_body = state
        .base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(boost_world)?;
    let centered_before = forward_kinematics(
        previous,
        BaseState::prescribed(
            state.base.world_from_base,
            Twist::new(
                state.base.twist_body.angular,
                state.base.twist_body.linear - boost_body,
            ),
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(cancelled)?;
    let mut link_momenta = Vec::with_capacity(previous.link_count());
    for (index, link) in previous.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        link_momenta.push(link.inertia().momentum(centered_before.body_twist[index])?);
    }
    // Per-link momenta are individual contributions, never accumulated RNEA
    // subtree wrenches. The result is the new velocity, not a velocity increment.
    let zero_joint_rhs = vec![0.0; previous.dof_count()];
    let solved =
        free_floating_impulse_delta(next, &state.joint_position, &zero_joint_rhs, &link_momenta)?;
    poll(cancelled)?;
    Ok((
        FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                state.base.world_from_base,
                Twist::new(solved.base.angular, solved.base.linear + boost_body),
            ),
            joint_position: state.joint_position.clone(),
            joint_velocity: solved.joint,
        },
        boost_world,
    ))
}

fn validate_compatible_models(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
) -> Result<(), ReconfigurationError> {
    for model in [previous, next] {
        if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
            return Err(ReconfigurationError::TooManyLinks {
                got: model.link_count(),
                maximum: MAX_FREE_FLOATING_STEP_LINKS,
            });
        }
        if model.links()[0].joint().kind() != JointKind::Fixed {
            return Err(ArticulatedError::FloatingBaseRootJointNotFixed.into());
        }
    }
    if previous.link_count() != next.link_count() || previous.dof_count() != next.dof_count() {
        return Err(ReconfigurationError::Invalid("tree dimensions changed"));
    }
    for (old, new) in previous.links().iter().zip(next.links()) {
        if old.name() != new.name()
            || old.parent() != new.parent()
            || old.parent_from_child_zero() != new.parent_from_child_zero()
            || old.joint() != new.joint()
        {
            return Err(ReconfigurationError::Invalid(
                "link geometry or joint changed",
            ));
        }
        let old_center = old.inertia().center_of_mass();
        let new_center = new.inertia().center_of_mass();
        if old.inertia().mass().to_bits() != new.inertia().mass().to_bits()
            || old_center.x.to_bits() != new_center.x.to_bits()
            || old_center.y.to_bits() != new_center.y.to_bits()
            || old_center.z.to_bits() != new_center.z.to_bits()
        {
            return Err(ReconfigurationError::Invalid("link mass or COM changed"));
        }
    }
    Ok(())
}

pub(in crate::articulated::motion) fn validate_controls(
    model: &ArticulatedModel,
    controls: ReconfigurationControls<'_>,
) -> Result<(), ReconfigurationError> {
    if !controls.relative_momentum_tolerance.is_finite()
        || !(0.0..=1e-8).contains(&controls.relative_momentum_tolerance)
    {
        return Err(ReconfigurationError::Invalid("relative momentum tolerance"));
    }
    if !controls.angular_momentum_tolerance_nms.is_finite()
        || controls.angular_momentum_tolerance_nms <= 0.0
        || !controls.linear_momentum_tolerance_ns.is_finite()
        || controls.linear_momentum_tolerance_ns <= 0.0
        || controls.joint_momentum_tolerances.len() != model.dof_count()
        || controls
            .joint_momentum_tolerances
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(ReconfigurationError::Invalid(
            "absolute momentum tolerances",
        ));
    }
    if !controls.maximum_parameter_work_j.is_finite() || controls.maximum_parameter_work_j < 0.0 {
        return Err(ReconfigurationError::Invalid(
            "absolute parameter-work budget",
        ));
    }
    Ok(())
}

pub(in crate::articulated::motion) fn state_kinematics(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
) -> Result<Kinematics, ReconfigurationError> {
    Ok(forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )?)
}

fn generalized_momentum(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
) -> Result<Vec<f64>, ReconfigurationError> {
    // RNEA with zero velocity and actual rates supplied as acceleration is an
    // independent M(q)u evaluation, without Coriolis or gravity contributions.
    let zero_rates = vec![0.0; model.dof_count()];
    let zero_wrenches = vec![Wrench::default(); model.link_count()];
    let evaluated = inverse_dynamics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            Twist::zero(),
            state.base.twist_body,
        ),
        &state.joint_position,
        &zero_rates,
        &state.joint_velocity,
        Vec3::default(),
        &zero_wrenches,
    )?;
    // Root link and floating base need not share an origin or orientation.
    let base = evaluated.body_wrench[0].transform_by(&model.links()[0].parent_from_child_zero());
    let mut result = Vec::with_capacity(6 + model.dof_count());
    result.extend_from_slice(&base.to_array());
    result.extend_from_slice(&evaluated.generalized_force);
    validate_slice(&result, "reconfigure.generalized_momentum")?;
    Ok(result)
}

pub(in crate::articulated::motion) fn audit_momentum(
    before: &[f64],
    after: &[f64],
    world_before: ArticulatedMomentum,
    world_after: ArticulatedMomentum,
    controls: ReconfigurationControls<'_>,
) -> Result<f64, ReconfigurationError> {
    let mut maximum = 0.0_f64;
    for (index, (&old, &new)) in before.iter().zip(after).enumerate() {
        let absolute = match index {
            0..=2 => controls.angular_momentum_tolerance_nms,
            3..=5 => controls.linear_momentum_tolerance_ns,
            _ => controls.joint_momentum_tolerances[index - 6],
        };
        maximum = maximum.max(normalized_residual(old, new, absolute, controls)?);
    }
    for (old, new, absolute) in [
        (
            world_before.angular_momentum_world,
            world_after.angular_momentum_world,
            controls.angular_momentum_tolerance_nms,
        ),
        (
            world_before.linear_momentum_world,
            world_after.linear_momentum_world,
            controls.linear_momentum_tolerance_ns,
        ),
    ] {
        for (old, new) in [(old.x, new.x), (old.y, new.y), (old.z, new.z)] {
            maximum = maximum.max(normalized_residual(old, new, absolute, controls)?);
        }
    }
    if maximum > 1.0 {
        return Err(ReconfigurationError::MomentumResidual {
            maximum_normalized: maximum,
        });
    }
    Ok(maximum)
}

fn normalized_residual(
    before: f64,
    after: f64,
    absolute: f64,
    controls: ReconfigurationControls<'_>,
) -> Result<f64, ReconfigurationError> {
    let allowed = absolute + controls.relative_momentum_tolerance * before.abs().max(after.abs());
    if !allowed.is_finite() {
        return Err(ReconfigurationError::Invalid("momentum tolerance overflow"));
    }
    Ok((after - before).abs() / allowed)
}

fn parameter_work(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
    before: &Kinematics,
    after: &Kinematics,
    boost_world: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ReconfigurationError> {
    let mut centered_work = 0.0;
    let mut actual_linear_residual = Vec3::default();
    for (index, (old, new)) in previous.links().iter().zip(next.links()).enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let rotation = before.world_from_link[index].rotation();
        let local_boost = Twist::new(Vec3::default(), rotation.inverse().rotate(boost_world)?);
        let old_velocity = before.body_twist[index];
        let new_velocity = after.body_twist[index];
        centered_work += new
            .inertia()
            .kinetic_energy(new_velocity.minus(local_boost))?
            - old
                .inertia()
                .kinetic_energy(old_velocity.minus(local_boost))?;
        // Linear momentum depends only on unchanged mass and COM. Evaluate its
        // actual difference directly, before a large common boost is summed.
        let delta = old
            .inertia()
            .momentum(new_velocity.minus(old_velocity))?
            .force;
        actual_linear_residual = actual_linear_residual + rotation.rotate(delta)?;
        ensure_finite_scalar(centered_work, "reconfigure.centered_work", index)?;
        validate_vec3(actual_linear_residual, "reconfigure.linear_residual")?;
    }
    // K_v = K - v.P + (m/2)|v|², hence ΔK = ΔK_v + v.ΔP for fixed mass.
    // Retain the finite returned state's residual; do not declare it zero.
    let work = centered_work + boost_world.dot(actual_linear_residual);
    ensure_finite_scalar(work, "reconfigure.parameter_work", 0)?;
    Ok(work)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReconfigurationError> {
    if cancelled() {
        Err(ReconfigurationError::Cancelled)
    } else {
        Ok(())
    }
}
