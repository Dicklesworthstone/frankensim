//! Constant affine joint locks and gears in free-floating implicit motion.
//!
//! At each actual trial the original free ABA supplies a_free. The shared
//! native row owner solves (A M_free^-1 A^T) lambda = -A a_free and adds its
//! joint response. M_free includes the carrier's finite inertia and recoil.
//! The original momentum-coordinate equations reconstruct that carrier's
//! actual midpoint and endpoint twist from the constrained joint motion.
//!
//! Spatial loads, finite bushings and smooth actuator feedback are evaluated
//! once at the same trial. Ideal reactions are internal: no reaction enters
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

use super::{
    ImplicitSpatialConfig, ImplicitSpatialError, ImplicitSpatialStep,
    step_with_bushings_effort_and_reaction,
};
use crate::articulated::ArticulatedModel;
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
    validate_condition_floor,
};
use crate::articulated::motion::spatial_loads::SpatialMotionLoads;
use crate::articulated::motion::spatial_loads::bushing::coupled::CoupledSpatialBushing;
use crate::articulated::motion::spatial_loads::bushing::{BushingLaw, SpatialBushing};
use crate::articulated::motion::{
    EffortLaw, FreeFloatingMotionError, FreeFloatingMotionState, MAX_FREE_FLOATING_STEP_LINKS,
    poll, sum_products,
};
use fs_ga::Twist;
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
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedFreeSpatialStep, ImplicitSpatialError> {
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
    let implicit = if rows.is_empty() {
        super::step_with_bushings_and_effort(
            model,
            state,
            loads,
            bushings,
            duration,
            config,
            effort,
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
        step_with_bushings_effort_and_reaction(
            model,
            state,
            loads,
            bushings,
            duration,
            config,
            effort,
            Some(&mut reaction),
            &mut cancelled,
        )?
    };
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
    Ok(ConstrainedFreeSpatialStep {
        implicit,
        midpoint_reaction: accepted,
        position_residuals_before: initial.position,
        position_residuals_midpoint: middle.position,
        position_residuals_after: final_rows.position,
        velocity_residuals_before: initial.velocity,
        velocity_residuals_midpoint: middle.velocity,
        velocity_residuals_after: final_rows.velocity,
    })
}
