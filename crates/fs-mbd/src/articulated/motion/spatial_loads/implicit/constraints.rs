//! Constant affine joint locks and gears in free-floating implicit motion.
//!
//! At each actual trial the original free ABA supplies a_free. The shared
//! native row owner solves (A M_free^-1 A^T) lambda = -A a_free and adds its
//! joint response. M_free includes the carrier's finite inertia and recoil.
//! The original momentum-coordinate equations reconstruct that carrier's
//! actual midpoint and endpoint twist from the constrained joint motion.
//!
//! Spatial loads, finite bushings and smooth actuator/external-force feedback
//! are evaluated once at the same trial. Ideal reactions are internal: no reaction enters
//! external spatial impulse and the carrier retains its finite inertia. Native
//! inverse-mass probes form the bounded row mobility; the loaded dynamics and
//! force laws are not evaluated a second time. Caller actuator work and signed
//! joint-reaction work
//! remain separate from connection storage and damping loss.
//!
//! Rows A q = b are constant and tangent rates satisfy A qdot = 0. Initial,
//! accepted midpoint and final residuals must meet the original row tolerances.
//! Numerical derivative probes may leave the affine manifold; no projection or
//! clipping alters their equations. This does not activate a constraint, model
//! backlash or varying ratios, or replace the nonlinear material-point loops.

use super::{ImplicitSpatialConfig, ImplicitSpatialError, ImplicitSpatialStep, step_with_forcing};
use crate::articulated::ArticulatedModel;
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
    validate_condition_floor,
};
use crate::articulated::motion::spatial_loads::SpatialMotionLoads;
use crate::articulated::motion::spatial_loads::bushing::coupled::CoupledSpatialBushing;
use crate::articulated::motion::spatial_loads::bushing::{BushingLaw, SpatialBushing};
use crate::articulated::motion::{
    EffortLaw, ForcingLaw, FreeFloatingMotionError, FreeFloatingMotionState,
    MAX_FREE_FLOATING_STEP_LINKS, poll, sum_products,
};
use fs_ga::{Twist, Vec3, Wrench};
use fs_math::dd::Dd;

/// Row regularity and signed ideal-reaction work admission.
#[derive(Debug, Clone, Copy)]
pub struct FreeConstraintControls {
    /// Native reciprocal condition floor in [4096*EPSILON, 1].
    pub minimum_reciprocal_condition: f64,
    /// Nonnegative maximum absolute midpoint joint-reaction work, in joules.
    /// Exact tangent reactions do no work; finite defects remain measurable.
    pub maximum_reaction_work_j: f64,
}

/// Original native reactions at the independently rechecked accepted midpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeConstraintReaction {
    /// Row multipliers in caller order and scaling (N m or N).
    pub multipliers: Vec<f64>,
    /// Original compact joint forces/torques A^T lambda (N or N m).
    pub joint_effort: Vec<f64>,
    /// Native carrier acceleration response to these internal joint reactions,
    /// expressed in the actual base axes (angular rad/s^2, linear m/s^2).
    pub base_acceleration_delta: Twist,
    /// Actual A qddot after the native response, in row units per second squared.
    pub acceleration_residuals: Vec<f64>,
    /// Actual joint-reaction effort paired with native relative joint rates, W.
    pub power_w: f64,
    /// Reciprocal condition of the normalized row mobility; one for empty rows.
    pub reciprocal_condition: f64,
}

/// One free implicit step with independently admitted affine/tangent rows.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedFreeSpatialStep {
    /// Original motion and spatial storage/loss. Caller actuator_work excludes
    /// reactions, whose signed work is retained as constraint_work.
    pub implicit: ImplicitSpatialStep,
    pub midpoint_reaction: FreeConstraintReaction,
    /// A q - b in original row order and units at the three admitted states.
    pub position_residuals_before: Vec<f64>,
    pub position_residuals_midpoint: Vec<f64>,
    pub position_residuals_after: Vec<f64>,
    /// A qdot in original row order and units per second.
    pub velocity_residuals_before: Vec<f64>,
    pub velocity_residuals_midpoint: Vec<f64>,
    pub velocity_residuals_after: Vec<f64>,
}

/// Simultaneous environmental forcing and internal reactions, with the actual
/// finite translational solve defect retained in physical world coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedForcedFreeSpatialStep {
    /// Original constrained state, accepted reactions, residuals and work.
    pub step: ConstrainedFreeSpatialStep,
    /// Trial minus momentum-derived midpoint COM velocity, in world m/s.
    /// Its torque-origin contribution is retained by the original finalizer.
    pub com_velocity_residual_world: Vec3,
}

/// Advance world/axial loads and constant ideal joint rows with carrier recoil.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_spatial(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    step_free_floating_bushings(
        model,
        state,
        loads,
        &[],
        rows,
        duration,
        config,
        controls,
        cancelled,
    )
}

/// Jointly solve smooth caller feedback, original spatial loads and ideal rows.
///
/// Feedback follows the parent repeatable-trial contract, once per residual
/// evaluation. Its physical effort limits apply before passive row reactions.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_spatial_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    effort: impl FnMut(f64, &FreeFloatingMotionState, &mut [f64]) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    step_free_floating_bushings_with_effort(
        model,
        state,
        loads,
        &[],
        rows,
        duration,
        config,
        controls,
        effort,
        cancelled,
    )
}

/// Resolve environmental forcing, internal actuation and ideal joint rows
/// together at the actual implicit midpoint.
///
/// The callback receives elapsed time h/2 and the full trial state, including
/// COM velocity and carrier recoil. Link-body wrenches about the original link
/// origins and relative joint efforts are reset to the held values every call.
/// The native loaded ABA runs once; the existing free inverse-mass row solve
/// then adds internal reactions. These have separate signed constraint work,
/// and never enter external spatial impulse or the caller's actuator limits.
///
/// The accepted force/reaction samples supply all accounts without reevaluating
/// the callback. Initial, midpoint and final row tolerances and the independent
/// reaction-work allowance still apply. No projection repairs a loose solve.
/// Empty rows exactly reproduce the original unconstrained forcing lane.
///
/// Callbacks must be smooth and repeatable, including at off-manifold numerical
/// derivative probes, and must not commit history during trials. The parent's
/// linear midpoint stability, finite residual, timestep, limit, cancellation
/// and atomic-refusal contracts apply; callback side effects remain its own.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_spatial_with_forcing(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    forcing: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedForcedFreeSpatialStep, ImplicitSpatialError> {
    step_free_floating_bushings_with_forcing(
        model,
        state,
        loads,
        &[],
        rows,
        duration,
        config,
        controls,
        forcing,
        cancelled,
    )
}

/// Solve original finite diagonal bushings with ideal joint rows and recoil.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    advance(
        model, state, loads, bushings, rows, duration, config, controls, None, cancelled,
    )
}

/// Solve diagonal bushings, feedback and ideal rows at the same actual stage.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_bushings_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    advance(
        model,
        state,
        loads,
        bushings,
        rows,
        duration,
        config,
        controls,
        Some(&mut effort),
        cancelled,
    )
}

/// Resolve the same external forcing and ideal rows with finite diagonal
/// bushings. Attachment loads, recoil and callback loads use one trial state;
/// connection storage/loss and internal-reaction work remain distinct.
/// See [`step_free_floating_spatial_with_forcing`] for callback and refusal rules.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_bushings_with_forcing(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    mut forcing: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedForcedFreeSpatialStep, ImplicitSpatialError> {
    advance_with_forcing(
        model,
        state,
        loads,
        bushings,
        rows,
        duration,
        config,
        controls,
        None,
        Some(&mut forcing),
        cancelled,
    )
}

/// Solve original full passive 6x6 K/C bushings with ideal rows and carrier recoil.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_coupled_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    advance(
        model, state, loads, bushings, rows, duration, config, controls, None, cancelled,
    )
}

/// Solve full passive coupling, smooth feedback and ideal rows simultaneously.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_coupled_bushings_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    advance(
        model,
        state,
        loads,
        bushings,
        rows,
        duration,
        config,
        controls,
        Some(&mut effort),
        cancelled,
    )
}

/// Resolve full passive 6x6 K/C coupling, environmental force, internal effort
/// and ideal joint rows through the same native midpoint and reaction owners.
/// See [`step_free_floating_spatial_with_forcing`] for callback and refusal rules.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_coupled_bushings_with_forcing(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    mut forcing: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedForcedFreeSpatialStep, ImplicitSpatialError> {
    advance_with_forcing(
        model,
        state,
        loads,
        bushings,
        rows,
        duration,
        config,
        controls,
        None,
        Some(&mut forcing),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    effort: Option<&mut EffortLaw<'_>>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
    Ok(advance_with_forcing(
        model, state, loads, bushings, rows, duration, config, controls, effort, None, cancelled,
    )?
    .step)
}

#[allow(clippy::too_many_arguments)]
fn advance_with_forcing<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: FreeConstraintControls,
    effort: Option<&mut EffortLaw<'_>>,
    forcing: Option<&mut ForcingLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedForcedFreeSpatialStep, ImplicitSpatialError> {
    poll(&mut cancelled)?;
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    if state.joint_position.len() != model.dof_count()
        || state.joint_velocity.len() != model.dof_count()
    {
        return Err(ConstraintError::Invalid("free constraint state dimensions").into());
    }
    validate_condition_floor(controls.minimum_reciprocal_condition)?;
    if !controls.maximum_reaction_work_j.is_finite() || controls.maximum_reaction_work_j < 0.0 {
        return Err(ConstraintError::Invalid("absolute reaction-work budget").into());
    }
    let initial = admit_constraints(
        rows,
        &state.joint_position,
        &state.joint_velocity,
        &mut cancelled,
    )?;
    let mut accepted = FreeConstraintReaction {
        multipliers: Vec::new(),
        joint_effort: vec![0.0; model.dof_count()],
        base_acceleration_delta: Twist::zero(),
        acceleration_residuals: Vec::new(),
        power_w: 0.0,
        reciprocal_condition: 1.0,
    };
    let forced = if rows.is_empty() {
        step_with_forcing(
            model,
            state,
            loads,
            bushings,
            duration,
            config,
            effort,
            None,
            forcing,
            &mut cancelled,
        )?
    } else {
        // The parent independently repeats the final residual. Retaining only
        // that last reaction associates this receipt with the accepted midpoint.
        let mut reaction = |stage: &FreeFloatingMotionState,
                            acceleration: &mut [f64],
                            stop: &mut dyn FnMut() -> bool|
         -> Result<Dd, ImplicitSpatialError> {
            let mut cancelled = || stop();
            let mut prepared = PreparedJointConstraints::free_trial(
                model,
                &stage.joint_position,
                rows,
                controls.minimum_reciprocal_condition,
                &mut cancelled,
            )?;
            let mut rhs = apply_rows(rows, acceleration, &mut cancelled)?;
            for value in &mut rhs {
                *value = -*value;
            }
            let response = prepared.solve_response(&rhs, &mut cancelled)?;
            for (a, delta) in acceleration.iter_mut().zip(&response.joint_delta) {
                *a += delta;
            }
            let residuals = apply_rows(rows, acceleration, &mut cancelled)?;
            let power = sum_products(
                response.joint_effort.iter().copied(),
                stage.joint_velocity.iter().copied(),
            );
            accepted = FreeConstraintReaction {
                multipliers: response.multipliers,
                joint_effort: response.joint_effort,
                base_acceleration_delta: response.base_delta,
                acceleration_residuals: residuals,
                power_w: power.to_f64(),
                reciprocal_condition: prepared.reciprocal_condition,
            };
            Ok(power)
        };
        step_with_forcing(
            model,
            state,
            loads,
            bushings,
            duration,
            config,
            effort,
            Some(&mut reaction),
            forcing,
            &mut cancelled,
        )?
    };
    let implicit = forced.step;
    let middle = admit_constraints(
        rows,
        &implicit.midpoint_state.joint_position,
        &implicit.midpoint_state.joint_velocity,
        &mut cancelled,
    )?;
    let endpoint = &implicit.spatial.motion.state_after;
    let final_rows = admit_constraints(
        rows,
        &endpoint.joint_position,
        &endpoint.joint_velocity,
        &mut cancelled,
    )?;
    let work = implicit.spatial.motion.constraint_work;
    if work.abs() > controls.maximum_reaction_work_j {
        return Err(ConstraintError::ReactionWork {
            work_j: work,
            maximum_j: controls.maximum_reaction_work_j,
        }
        .into());
    }
    poll(&mut cancelled)?;
    Ok(ConstrainedForcedFreeSpatialStep {
        step: ConstrainedFreeSpatialStep {
            implicit,
            midpoint_reaction: accepted,
            position_residuals_before: initial.position,
            position_residuals_midpoint: middle.position,
            position_residuals_after: final_rows.position,
            velocity_residuals_before: initial.velocity,
            velocity_residuals_midpoint: middle.velocity,
            velocity_residuals_after: final_rows.velocity,
        },
        com_velocity_residual_world: forced.com_velocity_residual_world,
    })
}
