//! Ideal constant affine joint constraints in supported implicit motion.
//!
//! Each row is A q = b, with A qdot = 0. At every actual nonlinear trial the
//! original prescribed-base ABA acceleration is augmented by the native
//! inverse-mass response to A^T lambda, where
//! (A M^-1 A^T) lambda = -A a_free. The existing normalized row factorization
//! owns regularity and the solve. Only this bounded row-space matrix is dense.
//!
//! World loads, finite connections and smooth actuator feedback are evaluated
//! once at that same trial by the original implicit owner. Passive reactions
//! bypass actuator effort limits while genuine caller efforts retain them.
//! RNEA sees the complete resulting acceleration, preserving support force,
//! torque, external momentum exchange and signed support work. Joint-reaction
//! work has its own native account; it never becomes actuation or damping loss.
//!
//! Initial, accepted midpoint and final position/rate residuals must satisfy the
//! original row tolerances. Numerical probes can lie off the affine manifold;
//! clipping them would invalidate the nonlinear derivative. No state projection,
//! constraint activation, varying ratio, backlash or nonlinear point-loop model
//! is inferred. The underlying midpoint phase/geometry limitations still apply.

use super::{
    ImplicitPrescribedSpatialError, ImplicitPrescribedSpatialStep, ImplicitSpatialConfig,
    step_with_bushings_effort_and_reaction,
};
use crate::articulated::ArticulatedModel;
use crate::articulated::motion::prescribed::{
    EffortLaw, PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath, poll,
};
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
    validate_condition_floor,
};
use crate::articulated::motion::spatial_loads::SpatialMotionLoads;
use crate::articulated::motion::spatial_loads::bushing::coupled::CoupledSpatialBushing;
use crate::articulated::motion::spatial_loads::bushing::{BushingLaw, SpatialBushing};
use crate::articulated::motion::{
    FreeFloatingMotionError, MAX_FREE_FLOATING_STEP_LINKS, sum_products,
};
use fs_math::dd::Dd;

/// Row regularity and signed ideal-reaction work admission.
#[derive(Debug, Clone, Copy)]
pub struct SupportedConstraintControls {
    /// Native reciprocal condition floor in [4096*EPSILON, 1].
    pub minimum_reciprocal_condition: f64,
    /// Nonnegative maximum absolute midpoint joint-reaction work, in joules.
    /// Exact tangent reactions do no work; finite residuals remain measurable.
    pub maximum_reaction_work_j: f64,
}

/// Native ideal reactions at the independently rechecked accepted midpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct SupportedConstraintReaction {
    /// Row multipliers, in caller row order and scaling (N m or N).
    pub multipliers: Vec<f64>,
    /// Original compact joint forces/torques A^T lambda (N or N m).
    pub joint_effort: Vec<f64>,
    /// Actual A qddot after the native response, in row units per second squared.
    pub acceleration_residuals: Vec<f64>,
    /// Actual reaction effort paired with native joint rates, W.
    pub power_w: f64,
    /// Reciprocal condition of the normalized row mobility; one for empty rows.
    pub reciprocal_condition: f64,
}

/// One supported implicit step with independently admitted affine/tangent rows.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedPrescribedSpatialStep {
    /// Native motion and spatial storage/loss. Its actuator_work belongs only
    /// to the caller; constraint_work contains the separate ideal-reaction work.
    pub implicit: ImplicitPrescribedSpatialStep,
    pub midpoint_reaction: SupportedConstraintReaction,
    /// A q - b in original row order and units at the three admitted states.
    pub position_residuals_before: Vec<f64>,
    pub position_residuals_midpoint: Vec<f64>,
    pub position_residuals_after: Vec<f64>,
    /// A qdot in original row order and units per second.
    pub velocity_residuals_before: Vec<f64>,
    pub velocity_residuals_midpoint: Vec<f64>,
    pub velocity_residuals_after: Vec<f64>,
}

/// Advance world/axial loads with constant ideal joint locks or gear relations.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_spatial(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_prescribed_base_bushings(
        model,
        state,
        path,
        loads,
        &[],
        rows,
        duration,
        config,
        controls,
        cancelled,
    )
}

/// Solve smooth caller feedback together with the original spatial loads and rows.
///
/// Feedback follows the parent implicit repeatable-trial contract and is called
/// once per residual evaluation. Its native effort limits apply before reaction.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_spatial_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_prescribed_base_bushings_with_effort(
        model,
        state,
        path,
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

/// Advance finite diagonal bushings on the same constrained implicit stage.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    advance(
        model, state, path, loads, bushings, rows, duration, config, controls, None, cancelled,
    )
}

/// Advance the original diagonal bushings, actuator feedback and ideal rows.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    advance(
        model,
        state,
        path,
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

/// Advance full passive 6x6 bushings and ideal rows with original native geometry.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_coupled_bushings(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    advance(
        model, state, path, loads, bushings, rows, duration, config, controls, None, cancelled,
    )
}

/// Solve full passive K/C coupling, smooth feedback and ideal joint rows together.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_coupled_bushings_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[CoupledSpatialBushing],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    advance(
        model,
        state,
        path,
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
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    rows: &[JointConstraint],
    duration: f64,
    config: ImplicitSpatialConfig,
    controls: SupportedConstraintControls,
    effort: Option<&mut EffortLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
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
        return Err(ConstraintError::Invalid("supported constraint state dimensions").into());
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
    let mut accepted = SupportedConstraintReaction {
        multipliers: Vec::new(),
        joint_effort: vec![0.0; model.dof_count()],
        acceleration_residuals: Vec::new(),
        power_w: 0.0,
        reciprocal_condition: 1.0,
    };
    let implicit = if rows.is_empty() {
        super::step_with_bushings_and_effort(
            model,
            state,
            path,
            loads,
            bushings,
            duration,
            config,
            effort,
            &mut cancelled,
        )?
    } else {
        // Only the last evaluation is retained: the parent independently
        // re-evaluates and admits that final midpoint before endpoint publication.
        let mut reaction = |stage: &PrescribedBaseMotionState,
                            acceleration: &mut [f64],
                            stop: &mut dyn FnMut() -> bool|
         -> Result<Dd, ImplicitPrescribedSpatialError> {
            let mut cancelled = || stop();
            let mut prepared = PreparedJointConstraints::prescribed_trial(
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
            accepted = SupportedConstraintReaction {
                multipliers: response.multipliers,
                joint_effort: response.joint_effort,
                acceleration_residuals: residuals,
                power_w: power.to_f64(),
                reciprocal_condition: prepared.reciprocal_condition,
            };
            Ok(power)
        };
        step_with_bushings_effort_and_reaction(
            model,
            state,
            path,
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
    Ok(ConstrainedPrescribedSpatialStep {
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
