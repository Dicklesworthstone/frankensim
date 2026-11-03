//! Prescribed inertia histories with constant affine joint constraints.
//!
//! Original native dynamics advances the motion. At each actual force stage,
//! joint reactions enforce `A qddot = 0`; both tensor transfers retain `A qdot = 0`.

use super::super::{
    ReconfigurationControls, ReconfigurationError, validate_compatible_models,
    validate_controls as validate_reconfiguration_controls,
};
use super::{
    ConstrainedInertiaReconfiguration, ConstraintControls, ConstraintError, JointConstraint,
    PreparedJointConstraints, admit_constraints, apply_rows, transfer_inertia_constrained,
    validate_controls,
};
use crate::articulated::ArticulatedModel;
use crate::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating_with_effort_and_reaction, sum_products,
};
use core::fmt;
use fs_ga::Vec3;
use fs_math::dd::Dd;
use std::cell::RefCell;

/// Internal reactions evaluated at one actual native motion stage.
#[derive(Debug)]
pub struct ConstraintReactionStage {
    /// Native elapsed time, zero or half the requested duration [s].
    pub elapsed_s: f64,
    /// Reaction multipliers in the original row order and scaling.
    pub multipliers: Vec<f64>,
    /// `A^T multipliers`, in the original compact joint order [N or N m].
    pub joint_effort: Vec<f64>,
    /// Actual reaction effort dotted with the stage joint rates [W].
    /// Finite admitted tangency residuals are retained, never declared zero.
    pub power_w: f64,
    /// Caller effort dotted with stage rates BEFORE reactions are added [W].
    pub caller_actuator_power_w: f64,
    /// Condition estimate of the row-equilibrated mobility; one for no rows.
    pub reciprocal_condition: f64,
}

/// One fully accepted constrained changing-inertia interval.
#[derive(Debug)]
pub struct ConstrainedInertiaMotionStep {
    pub first_transfer: ConstrainedInertiaReconfiguration,
    /// Original native receipt with separate caller actuator and constraint work.
    pub motion: FreeFloatingMotionStep,
    pub final_transfer: ConstrainedInertiaReconfiguration,
    /// Start and midpoint reactions, in native evaluation order.
    pub smooth_reactions: Vec<ConstraintReactionStage>,
    /// Native midpoint quadrature of the caller's effort alone [J].
    pub caller_actuator_work_j: f64,
    /// Native midpoint quadrature of the smooth constraint reaction [J].
    pub smooth_reaction_work_j: f64,
    pub net_parameter_work_j: f64,
    /// Sum of absolute parameter work over BOTH transfers [J].
    pub absolute_parameter_work_j: f64,
    /// Both transfer reactions plus the smooth reaction work [J].
    pub net_reaction_work_j: f64,
    /// Sum of absolute work from both transfer reactions and smooth motion [J].
    pub absolute_reaction_work_j: f64,
    /// Actual endpoint kinetic change minus all supplied physical work [J].
    /// Numerical error is not heat. Large common boosts can obscure small
    /// changes in this lab-frame subtraction; transfer work remains conditioned.
    pub kinetic_work_residual: f64,
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
}

impl ConstrainedInertiaMotionStep {
    #[must_use]
    pub fn state_after(&self) -> &FreeFloatingMotionState {
        &self.final_transfer.state_after
    }

    #[must_use]
    pub fn into_state_after(self) -> FreeFloatingMotionState {
        self.final_transfer.state_after
    }
}

/// Atomic refusal with the original constraint or native motion cause retained.
#[derive(Debug)]
pub enum ConstrainedInertiaMotionError {
    Constraint(ConstraintError),
    Motion(FreeFloatingMotionError),
    ParameterWorkBudget { requested: f64, maximum: f64 },
    ReactionWorkBudget { requested: f64, maximum: f64 },
    InvalidAccounting(&'static str),
}

impl fmt::Display for ConstrainedInertiaMotionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Constraint(error) => write!(f, "constrained inertia motion refused: {error}"),
            Self::Motion(error) => write!(f, "native constrained motion refused: {error}"),
            Self::ParameterWorkBudget { requested, maximum } => write!(
                f,
                "gross parameter work {requested} J exceeds interval budget {maximum} J",
            ),
            Self::ReactionWorkBudget { requested, maximum } => write!(
                f,
                "gross constraint work {requested} J exceeds interval budget {maximum} J",
            ),
            Self::InvalidAccounting(what) => {
                write!(f, "constrained motion accounting refused: {what}")
            }
        }
    }
}
impl std::error::Error for ConstrainedInertiaMotionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Constraint(error) => Some(error),
            Self::Motion(error) => Some(error),
            _ => None,
        }
    }
}
impl From<ConstraintError> for ConstrainedInertiaMotionError {
    fn from(value: ConstraintError) -> Self {
        Self::Constraint(value)
    }
}
impl From<ReconfigurationError> for ConstrainedInertiaMotionError {
    fn from(value: ReconfigurationError) -> Self {
        Self::Constraint(ConstraintError::Reconfiguration(value))
    }
}
impl From<FreeFloatingMotionError> for ConstrainedInertiaMotionError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value)
    }
}

/// Advance prescribed inertia through ideal constant affine gear/lock constraints.
///
/// Each row is `A q = b`, with time-independent coefficients and offset. The
/// supplied state must already satisfy position and velocity tolerances. Native
/// stage reactions enforce `A qddot = 0`, so affine constraints remain compatible
/// with the original midpoint position/rate updates. No position correction,
/// nonlinear loop closure, rolling/contact law or constraint activation is inferred.
///
/// Models retain per-link mass, COM, topology, attachments, joint definitions and
/// limits. Their tensors sample elapsed times zero, half-duration and duration.
/// The initial and final constrained transfers surround one original native step
/// at midpoint inertia. The row mobility is rebuilt at each ACTUAL stage
/// configuration. One original ABA evaluates genuine caller efforts; its joint
/// acceleration receives the native ideal reaction response. Motor effort limits
/// apply to the caller alone. Passive reactions never consume a motor rating.
///
/// Smooth prescribed tensors, O(duration²)-accurate midpoint samples, smooth
/// resolved loads and adequate mechanical-mode resolution give second-order
/// refinement. The native explicit midpoint step is not energy-exact or
/// time-reversible. Parameter work and finite reaction work remain separate;
/// each complete-interval budget counts absolute work before signs can cancel.
///
/// # Errors
/// Invalid models, rows, duration or controls; constraint/state/effort admission;
/// native dynamics or row solve failure; work exhaustion; nonfinite accounting;
/// or cancellation. All inputs remain unchanged and no partial state is returned.
#[allow(clippy::too_many_arguments)]
pub fn advance(
    previous: &ArticulatedModel,
    midpoint: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    constraints: &[JointConstraint],
    controls: ConstraintControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedInertiaMotionStep, ConstrainedInertiaMotionError> {
    advance_with_effort(
        previous,
        midpoint,
        next,
        state,
        loads,
        duration,
        constraints,
        controls,
        |_, _, _| Ok(()),
        cancelled,
    )
}

/// Advance with the original native stage-dependent actuator callback.
///
/// Calls occur at elapsed times zero and half-duration. The first sees the
/// VIRTUAL constrained state transferred to midpoint inertia, not the supplied
/// old-inertia rates; the second sees the actual native predicted midpoint.
/// Each effort slice begins with `loads.generalized_force`. After the caller
/// returns, joint reactions are solved from those actual applied efforts.
/// Transfers introduce no callback or elapsed time. The callback owns its
/// constitutive storage/dissipation and side effects, which cannot be rolled back.
///
/// # Errors
/// The refusals of [`advance`] or the original callback error.
#[allow(clippy::too_many_arguments)]
pub fn advance_with_effort(
    previous: &ArticulatedModel,
    midpoint: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    constraints: &[JointConstraint],
    controls: ConstraintControls<'_>,
    effort: impl FnMut(f64, &FreeFloatingMotionState, &mut [f64]) -> Result<(), FreeFloatingMotionError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedInertiaMotionStep, ConstrainedInertiaMotionError> {
    poll(&mut cancelled)?;
    if !duration.is_finite() || duration <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidDuration.into());
    }
    validate_controls(controls)?;
    validate_reconfiguration_controls(previous, controls.reconfiguration)?;
    validate_compatible_models(previous, midpoint)?;
    poll(&mut cancelled)?;
    validate_compatible_models(midpoint, next)?;
    poll(&mut cancelled)?;
    let first_transfer = transfer_inertia_constrained(
        previous,
        midpoint,
        state,
        constraints,
        controls,
        &mut cancelled,
    )?;
    let (motion, smooth_reactions) = constrained_motion(
        midpoint,
        &first_transfer.state_after,
        loads,
        duration,
        constraints,
        controls.minimum_reciprocal_condition,
        effort,
        &mut cancelled,
    )?;
    let smooth_reaction_work_j = finite(motion.constraint_work, "smooth reaction work")?;
    let caller_actuator_work_j = finite(motion.actuator_work, "caller actuator work")?;
    let used_reaction = finite(
        first_transfer.reaction_work_j.abs() + smooth_reaction_work_j.abs(),
        "partial gross reaction work",
    )?;
    check_reaction_work(used_reaction, controls.maximum_reaction_work_j)?;
    poll(&mut cancelled)?;
    let final_transfer = transfer_inertia_constrained(
        midpoint,
        next,
        &motion.state_after,
        constraints,
        ConstraintControls {
            reconfiguration: ReconfigurationControls {
                maximum_parameter_work_j: controls.reconfiguration.maximum_parameter_work_j
                    - first_transfer.parameter_work_j.abs(),
                ..controls.reconfiguration
            },
            maximum_reaction_work_j: controls.maximum_reaction_work_j - used_reaction,
            ..controls
        },
        &mut cancelled,
    )?;
    let net_parameter_work_j = finite(
        first_transfer.parameter_work_j + final_transfer.parameter_work_j,
        "net parameter work",
    )?;
    let absolute_parameter_work_j = finite(
        first_transfer.parameter_work_j.abs() + final_transfer.parameter_work_j.abs(),
        "gross parameter work",
    )?;
    if absolute_parameter_work_j > controls.reconfiguration.maximum_parameter_work_j {
        return Err(ConstrainedInertiaMotionError::ParameterWorkBudget {
            requested: absolute_parameter_work_j,
            maximum: controls.reconfiguration.maximum_parameter_work_j,
        });
    }
    let net_reaction_work_j = finite(
        first_transfer.reaction_work_j + smooth_reaction_work_j + final_transfer.reaction_work_j,
        "net reaction work",
    )?;
    let absolute_reaction_work_j = finite(
        used_reaction + final_transfer.reaction_work_j.abs(),
        "gross reaction work",
    )?;
    check_reaction_work(absolute_reaction_work_j, controls.maximum_reaction_work_j)?;
    let before = &first_transfer.momentum_before;
    let after = &final_transfer.momentum_after;
    let kinetic_work_residual = finite(
        after.kinetic_energy
            - before.kinetic_energy
            - net_parameter_work_j
            - net_reaction_work_j
            - caller_actuator_work_j
            - motion.external_work
            - motion.internal_body_work,
        "endpoint kinetic work residual",
    )?;
    let linear_momentum_residual = after.linear_momentum_world
        - before.linear_momentum_world
        - motion.external_linear_impulse_world;
    let angular_momentum_residual = after.angular_momentum_world
        - before.angular_momentum_world
        - motion.external_angular_impulse_world;
    for v in [linear_momentum_residual, angular_momentum_residual] {
        for value in [v.x, v.y, v.z] {
            finite(value, "endpoint momentum residual")?;
        }
    }
    poll(&mut cancelled)?;
    Ok(ConstrainedInertiaMotionStep {
        first_transfer,
        motion,
        final_transfer,
        smooth_reactions,
        caller_actuator_work_j,
        smooth_reaction_work_j,
        net_parameter_work_j,
        absolute_parameter_work_j,
        net_reaction_work_j,
        absolute_reaction_work_j,
        kinetic_work_residual,
        linear_momentum_residual,
        angular_momentum_residual,
    })
}

#[allow(clippy::too_many_arguments)]
fn constrained_motion(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    rows: &[JointConstraint],
    minimum_reciprocal_condition: f64,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(FreeFloatingMotionStep, Vec<ConstraintReactionStage>), ConstrainedInertiaMotionError> {
    let cancellation = RefCell::new(cancelled);
    let failure = RefCell::new(None);
    let mut reactions = Vec::with_capacity(2);
    let mut caller = |elapsed, stage: &FreeFloatingMotionState, efforts: &mut [f64]| {
        let mut checkpoint = || (*cancellation.borrow_mut())();
        if let Err(error) = admit_constraints(
            rows,
            &stage.joint_position,
            &stage.joint_velocity,
            &mut checkpoint,
        ) {
            *failure.borrow_mut() = Some(error);
            return Err(FreeFloatingMotionError::InvalidLoad(
                "constraint stage admission",
            ));
        }
        effort(elapsed, stage, efforts)?;
        if !rows.is_empty() && checkpoint() {
            *failure.borrow_mut() = Some(ConstraintError::Cancelled);
            return Err(FreeFloatingMotionError::Cancelled);
        }
        Ok(())
    };
    let mut reaction = |elapsed,
                        stage: &FreeFloatingMotionState,
                        efforts: &[f64],
                        acceleration: &mut [f64],
                        stop: &mut dyn FnMut() -> bool| {
        let mut checkpoint = || stop();
        match reaction_stage(
            model,
            stage,
            elapsed,
            efforts,
            acceleration,
            rows,
            minimum_reciprocal_condition,
            &mut checkpoint,
        ) {
            Ok((reaction, power)) => {
                reactions.push(reaction);
                Ok(power)
            }
            Err(error) => {
                // The native hook has a narrower error type. Preserve the
                // original typed cause without publishing a partial step.
                *failure.borrow_mut() = Some(error);
                Err(FreeFloatingMotionError::InvalidLoad(
                    "constraint reaction solve",
                ))
            }
        }
    };
    let motion = step_free_floating_with_effort_and_reaction(
        model,
        state,
        loads,
        duration,
        &mut caller,
        &mut reaction,
        &mut || (*cancellation.borrow_mut())(),
    );
    if let Some(error) = failure.into_inner() {
        return Err(error.into());
    }
    Ok((motion?, reactions))
}

#[allow(clippy::too_many_arguments)]
fn reaction_stage(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    elapsed_s: f64,
    efforts: &[f64],
    acceleration: &mut [f64],
    rows: &[JointConstraint],
    minimum_reciprocal_condition: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ConstraintReactionStage, Dd), ConstraintError> {
    let caller_actuator_power_w = sum_products(
        efforts.iter().copied(),
        state.joint_velocity.iter().copied(),
    )
    .to_f64();
    if !caller_actuator_power_w.is_finite() {
        return Err(ConstraintError::Invalid("caller actuator power overflow"));
    }
    if rows.is_empty() {
        return Ok((
            ConstraintReactionStage {
                elapsed_s,
                multipliers: Vec::new(),
                joint_effort: vec![0.0; efforts.len()],
                power_w: 0.0,
                caller_actuator_power_w,
                reciprocal_condition: 1.0,
            },
            Dd::ZERO,
        ));
    }
    let rhs = apply_rows(rows, acceleration, cancelled)?
        .into_iter()
        .map(|value| -value)
        .collect::<Vec<_>>();
    let mut prepared = PreparedJointConstraints::new(
        model,
        &state.joint_position,
        rows,
        minimum_reciprocal_condition,
        cancelled,
    )?;
    let response = prepared.solve_response(&rhs, cancelled)?;
    let power = sum_products(
        response.joint_effort.iter().copied(),
        state.joint_velocity.iter().copied(),
    );
    let power_w = power.to_f64();
    if !power_w.is_finite() {
        return Err(ConstraintError::Invalid(
            "constraint reaction power overflow",
        ));
    }
    // Do not pass ideal reactions back through the motor effort validator or
    // repeat ABA. The original response includes carrier recoil, recovered by
    // the native momentum map at the resulting joint rates.
    for (value, delta) in acceleration.iter_mut().zip(&response.joint_delta) {
        *value += delta;
    }
    Ok((
        ConstraintReactionStage {
            elapsed_s,
            multipliers: response.multipliers,
            joint_effort: response.joint_effort,
            power_w,
            caller_actuator_power_w,
            reciprocal_condition: prepared.reciprocal_condition,
        },
        power,
    ))
}

fn finite(value: f64, what: &'static str) -> Result<f64, ConstrainedInertiaMotionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ConstrainedInertiaMotionError::InvalidAccounting(what))
    }
}
fn check_reaction_work(requested: f64, maximum: f64) -> Result<(), ConstrainedInertiaMotionError> {
    if requested > maximum {
        Err(ConstrainedInertiaMotionError::ReactionWorkBudget { requested, maximum })
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConstrainedInertiaMotionError> {
    if cancelled() {
        Err(FreeFloatingMotionError::Cancelled.into())
    } else {
        Ok(())
    }
}
