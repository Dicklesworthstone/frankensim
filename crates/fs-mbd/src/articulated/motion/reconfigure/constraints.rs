//! Fixed joint and gear manifolds during prescribed full-tensor changes.
//!
//! Constant affine rows `A q = b` retain their tangent motion while the native
//! articulated inertia changes. Internal reactions solve
//! `M_next u_next = M_previous u_previous + [0; A^T lambda]` and `A qdot_next=0`.
//! Only the bounded row-space mobility is dense; every mass action is native ABA.
//! This does not represent contact, moving gear ratios, backlash, compliance,
//! point/weld geometry or additional internal shape momentum.

/// Continuous native motion with stage-resolved ideal joint reactions.
#[path = "constraints/step.rs"]
pub mod step;

use super::{
    InertiaReconfiguration, ReconfigurationControls, ReconfigurationError, audit_momentum,
    generalized_momentum, parameter_work, raw_momentum_candidate, state_kinematics,
    transfer_inertia, validate_compatible_models,
};
use crate::articulated::motion::{FreeFloatingMotionState, MAX_FREE_FLOATING_STEP_LINKS};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, free_floating_impulse_delta,
    momentum_from_kinematics, prescribed_base_impulse_delta, validate_slice, validate_twist,
};
use core::fmt;
use fs_ga::{Twist, Wrench};
use fs_la::{LuWorkspace, LuWorkspaceError};

/// Maximum independently reacting constant joint rows in one operation.
pub const MAX_JOINT_CONSTRAINTS: usize = 32;

/// Physical coordinate represented by a complete constraint row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintCoordinate {
    /// Row position in radians, rate in rad/s, multiplier in N m s for an event.
    Angle,
    /// Row position in metres, rate in m/s, multiplier in N s for an event.
    Length,
}

/// One fixed ideal relation `coefficients dot q = position`.
///
/// Entries follow compact scalar-joint order. Each coefficient carries the
/// row-coordinate unit divided by that joint's coordinate unit: revolute and
/// helical coordinates are radians, prismatic coordinates are metres. A locked
/// joint uses one nonzero coefficient; gear ratios use several. Rescaling a
/// row also rescales its position and both absolute tolerances.
#[derive(Debug, Clone, PartialEq)]
pub struct JointConstraint {
    /// Unit of the complete row, including its residuals and multiplier.
    pub coordinate: ConstraintCoordinate,
    /// Finite, nonzero row with exactly one entry per scalar joint.
    pub coefficients: Vec<f64>,
    /// Constant right-hand side in the declared row-coordinate unit.
    pub position: f64,
    /// Positive absolute allowed position residual, in the declared row unit.
    pub position_tolerance: f64,
    /// Positive absolute allowed tangent-rate residual, in row units/second.
    pub velocity_tolerance: f64,
}

/// Momentum, physical work and independent-row admission for one tensor change.
#[derive(Debug, Clone, Copy)]
pub struct ConstraintControls<'a> {
    /// Original momentum tolerances and absolute PARAMETER-work allowance.
    pub reconfiguration: ReconfigurationControls<'a>,
    /// Finite nonnegative cap on absolute actual ideal-reaction midpoint work.
    /// Exact tangent motion has zero work; admitted finite residuals remain here.
    pub maximum_reaction_work_j: f64,
    /// Floor for the reciprocal infinity-norm condition number of the
    /// diagonal-normalized row mobility, in `[4096*EPSILON, 1]`.
    pub minimum_reciprocal_condition: f64,
}

/// One accepted tensor update and simultaneous ideal joint reactions.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedInertiaReconfiguration {
    /// Unchanged pose/coordinates with the complete constrained endpoint rates.
    pub state_after: FreeFloatingMotionState,
    /// Actual world momentum and lab-frame energy before the event.
    pub momentum_before: ArticulatedMomentum,
    /// Actual world momentum and lab-frame energy after the event.
    pub momentum_after: ArticulatedMomentum,
    /// Native `[base angular, base linear, joint]` momenta before the event.
    pub generalized_momentum_before: Vec<f64>,
    /// Independently reconstructed momenta after the event.
    pub generalized_momentum_after: Vec<f64>,
    /// Full `[0; A^T lambda]`, with six leading zero base components.
    pub generalized_constraint_impulse: Vec<f64>,
    /// Row impulses in input order. Units are N m s for Angle and N s for Length.
    pub multipliers: Vec<f64>,
    /// Signed conditioned kinetic change less actual reaction midpoint work, J.
    pub parameter_work_j: f64,
    /// `sum lambda_i * (A_i v_before + A_i v_after)/2`, J. No heat is inferred.
    pub reaction_work_j: f64,
    /// Conditioned total mechanical kinetic change, J; the sum of both works.
    pub kinetic_energy_change_j: f64,
    /// Maximum normalized reaction-adjusted generalized or world momentum error.
    pub maximum_normalized_momentum_residual: f64,
    /// Actual `A q - b` in input row order; positions are unchanged.
    pub position_residuals: Vec<f64>,
    /// Actual initial `A qdot`, in each row's coordinate units/second.
    pub velocity_residuals_before: Vec<f64>,
    /// Actual final `A qdot`, in each row's coordinate units/second.
    pub velocity_residuals_after: Vec<f64>,
    /// Numerical reciprocal condition of the normalized mobility; one if empty.
    pub reciprocal_condition: f64,
}

/// Atomic refusal without state, model, geometry or budget mutation.
#[derive(Debug, Clone, PartialEq)]
pub enum ConstraintError {
    /// Original model, state, momentum or parameter-work admission refused.
    Reconfiguration(ReconfigurationError),
    /// Original reusable LU owner refused factorization or substitution.
    Solve(LuWorkspaceError),
    /// Invalid row, controls, dimensions or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The row count exceeds the bounded independent constraint solve.
    TooManyRows {
        /// Requested rows.
        got: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Initial or returned coordinates do not satisfy the declared manifold.
    Position {
        /// Input row index.
        row: usize,
        /// Unit of the residual.
        coordinate: ConstraintCoordinate,
        /// Actual signed residual.
        residual: f64,
        /// Declared absolute tolerance.
        maximum: f64,
    },
    /// Initial or returned motion does not satisfy the tangent-rate allowance.
    Velocity {
        /// Input row index.
        row: usize,
        /// Coordinate whose time derivative supplies the residual unit.
        coordinate: ConstraintCoordinate,
        /// Actual signed row rate.
        residual: f64,
        /// Declared absolute tolerance.
        maximum: f64,
    },
    /// Dependent or numerically unresolved row mobility is never regularized.
    IllConditioned {
        /// Numerical reciprocal condition from retained-factor unit-RHS solves.
        reciprocal_condition: f64,
        /// Caller-required minimum.
        minimum: f64,
    },
    /// Actual ideal-reaction work exceeds its separate absolute allowance.
    ReactionWork {
        /// Signed reaction work in joules.
        work_j: f64,
        /// Declared absolute cap in joules.
        maximum_j: f64,
    },
    /// Cancellation was observed before a result could be published.
    Cancelled,
}

impl fmt::Display for ConstraintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "joint-constrained inertia refused: {self:?}")
    }
}
impl std::error::Error for ConstraintError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reconfiguration(source) => Some(source),
            Self::Solve(source) => Some(source),
            _ => None,
        }
    }
}
impl From<ReconfigurationError> for ConstraintError {
    fn from(error: ReconfigurationError) -> Self {
        match error {
            ReconfigurationError::Cancelled => Self::Cancelled,
            other => Self::Reconfiguration(other),
        }
    }
}
impl From<ArticulatedError> for ConstraintError {
    fn from(error: ArticulatedError) -> Self {
        ReconfigurationError::from(error).into()
    }
}
impl From<LuWorkspaceError> for ConstraintError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Solve(error)
    }
}

/// Change full link tensors while respecting constant ideal joint/gear rows.
///
/// Models retain the original inventory, COM, attachments and joint definitions.
/// The supplied configuration AND rates must already satisfy every constraint.
/// Original native ABA supplies the unconstrained momentum candidate and each
/// row response. One normalized row-space solve supplies all internal reactions;
/// only the complete endpoint undergoes native speed and work admission.
///
/// Generalized momentum is audited against `[0; A^T lambda]`, and actual world
/// linear/angular momentum is independently conserved. Parameter work uses the
/// original common-translation conditioning, less actual reaction midpoint work.
/// Both signed work channels have independent absolute caps. No rejected
/// intermediate velocity is published, clipped or charged as impact heat.
///
/// At most 32 independent rows and 4096 links are admitted. With fixed row count,
/// native mass actions and stored velocity columns scale linearly in tree size;
/// only the bounded row-space matrix is factored. Diagonal normalization removes
/// row-unit scaling, and retained-factor unit-RHS solves measure its numerical
/// infinity-norm condition. Singular/unresolved rows refuse without regularization.
/// Empty rows reproduce the original transfer; unchanged models retain state bits
/// after row/mobility admission. This event advances no configuration or time.
///
/// # Errors
/// Invalid models, rows, initial/final state, conditioning, momentum/work bounds,
/// native solve failure or cancellation. Every input remains unchanged.
#[allow(clippy::too_many_lines)]
pub fn transfer_inertia_constrained(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    constraints: &[JointConstraint],
    controls: ConstraintControls<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedInertiaReconfiguration, ConstraintError> {
    poll(&mut cancelled)?;
    validate_controls(controls)?;
    if constraints.is_empty() {
        let result = transfer_inertia(
            previous,
            next,
            state,
            controls.reconfiguration,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        return Ok(unconstrained_receipt(result, Vec::new(), Vec::new(), 1.0));
    }
    validate_compatible_models(previous, next)?;
    super::validate_controls(previous, controls.reconfiguration)?;
    let initial = admit_constraints(
        constraints,
        &state.joint_position,
        &state.joint_velocity,
        &mut cancelled,
    )?;
    let before = state_kinematics(previous, state)?;
    poll(&mut cancelled)?;
    let mut mobility = PreparedJointConstraints::new(
        next,
        &state.joint_position,
        constraints,
        controls.minimum_reciprocal_condition,
        &mut cancelled,
    )?;
    if previous == next {
        let result = transfer_inertia(
            previous,
            next,
            state,
            controls.reconfiguration,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        return Ok(unconstrained_receipt(
            result,
            initial.position,
            initial.velocity,
            mobility.reciprocal_condition,
        ));
    }
    let momentum_before = momentum_from_kinematics(previous, &before)?;
    poll(&mut cancelled)?;
    let generalized_momentum_before = generalized_momentum(previous, state)?;
    poll(&mut cancelled)?;
    let (mut state_after, boost) =
        raw_momentum_candidate(previous, next, state, momentum_before, &mut cancelled)?;
    let free_rates = apply_rows(constraints, &state_after.joint_velocity, &mut cancelled)?;
    let rhs: Vec<_> = free_rates.iter().map(|rate| -rate).collect();
    let response = mobility.solve_response(&rhs, &mut cancelled)?;
    state_after.base.twist_body = state_after.base.twist_body.plus(response.base_delta);
    for (index, (rate, delta)) in state_after
        .joint_velocity
        .iter_mut()
        .zip(&response.joint_delta)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        *rate += delta;
    }
    let endpoint = admit_constraints(
        constraints,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut cancelled,
    )?;
    let after = state_kinematics(next, &state_after)?;
    poll(&mut cancelled)?;
    let momentum_after = momentum_from_kinematics(next, &after)?;
    poll(&mut cancelled)?;
    let generalized_momentum_after = generalized_momentum(next, &state_after)?;
    poll(&mut cancelled)?;
    let mut generalized_constraint_impulse = vec![0.0; 6];
    generalized_constraint_impulse.extend_from_slice(&response.joint_effort);
    let expected: Vec<_> = generalized_momentum_before
        .iter()
        .zip(&generalized_constraint_impulse)
        .map(|(momentum, reaction)| momentum + reaction)
        .collect();
    validate_slice(&expected, "constrained_inertia.expected_momentum")?;
    let maximum_normalized_momentum_residual = audit_momentum(
        &expected,
        &generalized_momentum_after,
        momentum_before,
        momentum_after,
        controls.reconfiguration,
    )?;
    let kinetic_energy_change_j =
        parameter_work(previous, next, &before, &after, boost, &mut cancelled)?;
    let mut reaction_work_j = 0.0;
    for ((&impulse, &old), &new) in response
        .multipliers
        .iter()
        .zip(&initial.velocity)
        .zip(&endpoint.velocity)
    {
        reaction_work_j = finite(
            impulse.mul_add(0.5 * old + 0.5 * new, reaction_work_j),
            "reaction midpoint work",
        )?;
    }
    if reaction_work_j.abs() > controls.maximum_reaction_work_j {
        return Err(ConstraintError::ReactionWork {
            work_j: reaction_work_j,
            maximum_j: controls.maximum_reaction_work_j,
        });
    }
    let parameter_work_j = finite(kinetic_energy_change_j - reaction_work_j, "parameter work")?;
    if parameter_work_j.abs() > controls.reconfiguration.maximum_parameter_work_j {
        return Err(ReconfigurationError::ParameterWork {
            work_j: parameter_work_j,
            maximum_j: controls.reconfiguration.maximum_parameter_work_j,
        }
        .into());
    }
    poll(&mut cancelled)?;
    Ok(ConstrainedInertiaReconfiguration {
        state_after,
        momentum_before,
        momentum_after,
        generalized_momentum_before,
        generalized_momentum_after,
        generalized_constraint_impulse,
        multipliers: response.multipliers,
        parameter_work_j,
        reaction_work_j,
        kinetic_energy_change_j,
        maximum_normalized_momentum_residual,
        position_residuals: initial.position,
        velocity_residuals_before: initial.velocity,
        velocity_residuals_after: endpoint.velocity,
        reciprocal_condition: mobility.reciprocal_condition,
    })
}

fn unconstrained_receipt(
    result: InertiaReconfiguration,
    position: Vec<f64>,
    velocity: Vec<f64>,
    reciprocal_condition: f64,
) -> ConstrainedInertiaReconfiguration {
    ConstrainedInertiaReconfiguration {
        generalized_constraint_impulse: vec![0.0; result.generalized_momentum_before.len()],
        multipliers: vec![0.0; velocity.len()],
        state_after: result.state_after,
        momentum_before: result.momentum_before,
        momentum_after: result.momentum_after,
        generalized_momentum_before: result.generalized_momentum_before,
        generalized_momentum_after: result.generalized_momentum_after,
        parameter_work_j: result.parameter_work_j,
        reaction_work_j: 0.0,
        kinetic_energy_change_j: result.parameter_work_j,
        maximum_normalized_momentum_residual: result.maximum_normalized_momentum_residual,
        position_residuals: position,
        velocity_residuals_before: velocity.clone(),
        velocity_residuals_after: velocity,
        reciprocal_condition,
    }
}

pub(in crate::articulated::motion) fn validate_controls(
    controls: ConstraintControls<'_>,
) -> Result<(), ConstraintError> {
    if !controls.maximum_reaction_work_j.is_finite() || controls.maximum_reaction_work_j < 0.0 {
        return Err(ConstraintError::Invalid("absolute reaction-work budget"));
    }
    validate_condition_floor(controls.minimum_reciprocal_condition)
}

pub(in crate::articulated::motion) fn validate_condition_floor(
    minimum: f64,
) -> Result<(), ConstraintError> {
    if !minimum.is_finite() || !(4096.0 * f64::EPSILON..=1.0).contains(&minimum) {
        return Err(ConstraintError::Invalid("reciprocal-condition floor"));
    }
    Ok(())
}

fn validate_rows(
    rows: &[JointConstraint],
    joints: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ConstraintError> {
    if rows.len() > MAX_JOINT_CONSTRAINTS {
        return Err(ConstraintError::TooManyRows {
            got: rows.len(),
            maximum: MAX_JOINT_CONSTRAINTS,
        });
    }
    for row in rows {
        poll(cancelled)?;
        if row.coefficients.len() != joints
            || !row.position.is_finite()
            || !row.position_tolerance.is_finite()
            || row.position_tolerance <= 0.0
            || !row.velocity_tolerance.is_finite()
            || row.velocity_tolerance <= 0.0
        {
            return Err(ConstraintError::Invalid(
                "joint row dimensions or tolerances",
            ));
        }
        let mut nonzero = false;
        for (index, &value) in row.coefficients.iter().enumerate() {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            finite(value, "joint row coefficient")?;
            nonzero |= value != 0.0;
        }
        if !nonzero {
            return Err(ConstraintError::Invalid("zero joint constraint row"));
        }
    }
    Ok(())
}

pub(in crate::articulated::motion) fn apply_rows(
    rows: &[JointConstraint],
    values: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, ConstraintError> {
    let mut output = Vec::with_capacity(rows.len());
    for row in rows {
        poll(cancelled)?;
        if row.coefficients.len() != values.len() {
            return Err(ConstraintError::Invalid("joint row input dimension"));
        }
        let mut value = 0.0;
        for (index, (&coefficient, &entry)) in row.coefficients.iter().zip(values).enumerate() {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            value = finite(coefficient.mul_add(entry, value), "joint row action")?;
        }
        output.push(value);
    }
    Ok(output)
}

fn admit_positions(
    rows: &[JointConstraint],
    position: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, ConstraintError> {
    let mut residuals = apply_rows(rows, position, cancelled)?;
    for (index, (residual, row)) in residuals.iter_mut().zip(rows).enumerate() {
        *residual = finite(*residual - row.position, "joint position residual")?;
        if residual.abs() > row.position_tolerance {
            return Err(ConstraintError::Position {
                row: index,
                coordinate: row.coordinate,
                residual: *residual,
                maximum: row.position_tolerance,
            });
        }
    }
    Ok(residuals)
}

pub(in crate::articulated::motion) struct ConstraintResiduals {
    pub(in crate::articulated::motion) position: Vec<f64>,
    pub(in crate::articulated::motion) velocity: Vec<f64>,
}

pub(in crate::articulated::motion) fn admit_constraints(
    rows: &[JointConstraint],
    position: &[f64],
    velocity: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ConstraintResiduals, ConstraintError> {
    poll(cancelled)?;
    if position.len() != velocity.len() {
        return Err(ConstraintError::Invalid("joint position/rate dimensions"));
    }
    validate_rows(rows, position.len(), cancelled)?;
    let position = admit_positions(rows, position, cancelled)?;
    let velocity = apply_rows(rows, velocity, cancelled)?;
    for (index, (&residual, row)) in velocity.iter().zip(rows).enumerate() {
        if residual.abs() > row.velocity_tolerance {
            return Err(ConstraintError::Velocity {
                row: index,
                coordinate: row.coordinate,
                residual,
                maximum: row.velocity_tolerance,
            });
        }
    }
    poll(cancelled)?;
    Ok(ConstraintResiduals { position, velocity })
}

#[derive(Clone, Copy)]
enum ConstraintBoundary {
    Free,
    FreeTrial,
    PrescribedTrial,
}

struct VelocityColumn {
    base: Twist,
    joint: Vec<f64>,
}

/// Cached original mass responses and one normalized small row factorization.
pub(in crate::articulated::motion) struct PreparedJointConstraints<'a> {
    rows: &'a [JointConstraint],
    joints: usize,
    columns: Vec<VelocityColumn>,
    scale: Vec<f64>,
    factor: LuWorkspace,
    pub(in crate::articulated::motion) reciprocal_condition: f64,
}

pub(in crate::articulated::motion) struct ConstraintResponse {
    pub(in crate::articulated::motion) multipliers: Vec<f64>,
    pub(in crate::articulated::motion) joint_effort: Vec<f64>,
    pub(in crate::articulated::motion) base_delta: Twist,
    pub(in crate::articulated::motion) joint_delta: Vec<f64>,
}

impl<'a> PreparedJointConstraints<'a> {
    #[allow(clippy::too_many_lines)]
    pub(in crate::articulated::motion) fn new(
        model: &ArticulatedModel,
        joint_position: &[f64],
        rows: &'a [JointConstraint],
        minimum_reciprocal_condition: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConstraintError> {
        Self::for_boundary(
            model,
            joint_position,
            rows,
            minimum_reciprocal_condition,
            ConstraintBoundary::Free,
            cancelled,
        )
    }

    /// Free-carrier row mobility at an actual nonlinear trial configuration.
    /// Trial probes retain native model/configuration admission but may lie off
    /// the affine manifold; the motion owner admits actual accepted row states.
    pub(in crate::articulated::motion) fn free_trial(
        model: &ArticulatedModel,
        joint_position: &[f64],
        rows: &'a [JointConstraint],
        minimum_reciprocal_condition: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConstraintError> {
        Self::for_boundary(
            model,
            joint_position,
            rows,
            minimum_reciprocal_condition,
            ConstraintBoundary::FreeTrial,
            cancelled,
        )
    }

    /// Constant rows remain meaningful away from their affine manifold during
    /// nonlinear probes. Actual initial and accepted states are admitted by the
    /// constrained motion owner. Native joint configuration limits still apply.
    pub(in crate::articulated::motion) fn prescribed_trial(
        model: &ArticulatedModel,
        joint_position: &[f64],
        rows: &'a [JointConstraint],
        minimum_reciprocal_condition: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConstraintError> {
        Self::for_boundary(
            model,
            joint_position,
            rows,
            minimum_reciprocal_condition,
            ConstraintBoundary::PrescribedTrial,
            cancelled,
        )
    }

    #[allow(clippy::too_many_lines)]
    fn for_boundary(
        model: &ArticulatedModel,
        joint_position: &[f64],
        rows: &'a [JointConstraint],
        minimum_reciprocal_condition: f64,
        boundary: ConstraintBoundary,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConstraintError> {
        poll(cancelled)?;
        validate_condition_floor(minimum_reciprocal_condition)?;
        match boundary {
            ConstraintBoundary::Free | ConstraintBoundary::FreeTrial => {
                validate_compatible_models(model, model)?;
                model.validate_configuration(joint_position)?;
            }
            ConstraintBoundary::PrescribedTrial => {
                if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
                    return Err(ReconfigurationError::TooManyLinks {
                        got: model.link_count(),
                        maximum: MAX_FREE_FLOATING_STEP_LINKS,
                    }
                    .into());
                }
                model.validate_configuration(joint_position)?;
            }
        }
        let joints = model.dof_count();
        if joint_position.len() != joints {
            return Err(ConstraintError::Invalid(
                "prepared joint position dimension",
            ));
        }
        validate_rows(rows, joints, cancelled)?;
        if matches!(boundary, ConstraintBoundary::Free) {
            admit_positions(rows, joint_position, cancelled)?;
        }
        let n = rows.len();
        let mut factor = LuWorkspace::new(n)?;
        if n == 0 {
            factor.factor(&[])?;
            return Ok(Self {
                rows,
                joints,
                columns: Vec::new(),
                scale: Vec::new(),
                factor,
                reciprocal_condition: 1.0,
            });
        }
        let mut columns = Vec::with_capacity(n);
        let mut mobility = vec![0.0; n * n];
        let zero_body_impulses = vec![Wrench::default(); model.link_count()];
        for (column, constraint) in rows.iter().enumerate() {
            poll(cancelled)?;
            let delta = match boundary {
                ConstraintBoundary::Free | ConstraintBoundary::FreeTrial => {
                    free_floating_impulse_delta(
                        model,
                        joint_position,
                        &constraint.coefficients,
                        &zero_body_impulses,
                    )?
                }
                ConstraintBoundary::PrescribedTrial => {
                    prescribed_base_impulse_delta(
                        model,
                        joint_position,
                        &constraint.coefficients,
                        &zero_body_impulses,
                    )?
                    .velocity
                }
            };
            poll(cancelled)?;
            let action = apply_rows(rows, &delta.joint, cancelled)?;
            for (row, value) in action.into_iter().enumerate() {
                mobility[row * n + column] = value;
            }
            columns.push(VelocityColumn {
                base: delta.base,
                joint: delta.joint,
            });
        }
        let mut scale = Vec::with_capacity(n);
        for row in 0..n {
            let diagonal = mobility[row * n + row];
            if !diagonal.is_finite() || diagonal <= 0.0 {
                return Err(ConstraintError::Invalid("nonpositive joint row mobility"));
            }
            scale.push(diagonal.sqrt());
        }
        let mut normalized = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            normalized[i * n + i] = 1.0;
            for j in 0..i {
                let a = finite(
                    (mobility[i * n + j] / scale[i]) / scale[j],
                    "normalized mobility",
                )?;
                let b = finite(
                    (mobility[j * n + i] / scale[j]) / scale[i],
                    "normalized mobility",
                )?;
                if (a - b).abs() > 4096.0 * f64::EPSILON * a.abs().max(b.abs()).max(1.0) {
                    return Err(ConstraintError::Invalid("joint mobility reciprocity"));
                }
                let entry = 0.5 * a + 0.5 * b;
                normalized[i * n + j] = entry;
                normalized[j * n + i] = entry;
            }
        }
        poll(cancelled)?;
        factor.factor(&normalized)?;
        poll(cancelled)?;
        // Exact column solves of this small factor give an independent matrix
        // infinity-norm condition diagnostic, rather than trusting pivot sizes.
        let norm = normalized
            .chunks_exact(n)
            .map(|row| row.iter().map(|value| value.abs()).sum::<f64>())
            .fold(0.0_f64, f64::max);
        let mut inverse_sums = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        let mut solution = vec![0.0; n];
        for column in 0..n {
            poll(cancelled)?;
            rhs[column] = 1.0;
            factor.solve_factored_into(&rhs, &mut solution)?;
            rhs[column] = 0.0;
            for (total, entry) in inverse_sums.iter_mut().zip(&solution) {
                *total = finite(*total + entry.abs(), "inverse mobility norm")?;
            }
        }
        let inverse_norm = inverse_sums.into_iter().fold(0.0_f64, f64::max);
        let reciprocal_condition = finite((1.0 / norm) / inverse_norm, "reciprocal condition")?;
        if reciprocal_condition < minimum_reciprocal_condition {
            return Err(ConstraintError::IllConditioned {
                reciprocal_condition,
                minimum: minimum_reciprocal_condition,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            rows,
            joints,
            columns,
            scale,
            factor,
            reciprocal_condition,
        })
    }

    pub(in crate::articulated::motion) fn solve_response(
        &mut self,
        rhs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConstraintResponse, ConstraintError> {
        poll(cancelled)?;
        if rhs.len() != self.rows.len() {
            return Err(ConstraintError::Invalid(
                "constraint right-hand side dimension",
            ));
        }
        let scaled: Vec<_> = rhs
            .iter()
            .zip(&self.scale)
            .map(|(&value, &scale)| finite(value / scale, "scaled constraint right-hand side"))
            .collect::<Result<_, _>>()?;
        let mut multipliers = vec![0.0; self.rows.len()];
        self.factor.solve_factored_into(&scaled, &mut multipliers)?;
        for (value, scale) in multipliers.iter_mut().zip(&self.scale) {
            *value = finite(*value / scale, "constraint multiplier")?;
        }
        let mut base_delta = Twist::zero();
        let mut joint_delta = vec![0.0; self.joints];
        let mut joint_effort = vec![0.0; self.joints];
        for ((row, column), &multiplier) in self.rows.iter().zip(&self.columns).zip(&multipliers) {
            poll(cancelled)?;
            base_delta = base_delta.plus(column.base.scale(multiplier));
            for index in 0..self.joints {
                if index % 64 == 0 {
                    poll(cancelled)?;
                }
                joint_effort[index] =
                    row.coefficients[index].mul_add(multiplier, joint_effort[index]);
                joint_delta[index] = column.joint[index].mul_add(multiplier, joint_delta[index]);
            }
        }
        validate_twist(base_delta, "constraint.base_response")?;
        validate_slice(&joint_delta, "constraint.joint_response")?;
        validate_slice(&joint_effort, "constraint.joint_effort")?;
        poll(cancelled)?;
        Ok(ConstraintResponse {
            multipliers,
            joint_effort,
            base_delta,
            joint_delta,
        })
    }
}

fn finite(value: f64, field: &'static str) -> Result<f64, ConstraintError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ConstraintError::Invalid(field))
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConstraintError> {
    if cancelled() {
        Err(ConstraintError::Cancelled)
    } else {
        Ok(())
    }
}
