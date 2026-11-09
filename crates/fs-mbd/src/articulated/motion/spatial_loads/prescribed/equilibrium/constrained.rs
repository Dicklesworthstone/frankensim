//! Loaded equilibrium on the original constant joint-lock/gear manifold.
//!
//! Solve the native equations R(q)-u-A^T lambda=0 and A q=b together. Ideal
//! reactions never pass through actuator limits, penalty springs or invented
//! dissipation. The same RNEA/load owner and Newton--Krylov dispatch serve the
//! unconstrained problem. Only the existing bounded row-mobility admission is
//! dense; the augmented tangent remains matrix free.

use super::*;
use crate::articulated::motion::reconfigure::constraints::{
    JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
};
use fs_math::dd::Dd;

/// Explicit units for the extra unknowns and equations.
#[derive(Clone, Copy, Debug)]
pub struct JointEquilibriumControls<'a> {
    /// Positive multiplier scales: N m for Angle rows, N for Length rows.
    pub reaction_scales: &'a [f64],
    /// Positive position-residual scales in each row's declared coordinate unit.
    /// Native row.position_tolerance is ALSO enforced at the accepted root.
    pub position_scales: &'a [f64],
    /// Original normalized row-mobility condition floor, in [4096*EPSILON,1].
    pub minimum_reciprocal_condition: f64,
}

/// Independently balanced original state and separately accounted ideal forces.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstrainedArticulatedEquilibrium {
    pub state: PrescribedBaseMotionState,
    /// Original RNEA holding efforts, before subtracting any applied forces.
    pub required_effort: Vec<f64>,
    /// Caller actuator efforts, independently admitted against actuator limits.
    pub applied_effort: Vec<f64>,
    /// Row multipliers in caller order and scaling: N m or N by row coordinate.
    pub multipliers: Vec<f64>,
    /// A^T lambda in original joint effort units. This is not actuator effort.
    pub reaction_effort: Vec<f64>,
    /// Required minus applied minus ideal reaction, independently recomputed.
    pub joint_effort_residual: Vec<f64>,
    pub scaled_force_residual_norm: f64,
    /// A q-b in original row units, admitted against each original tolerance.
    pub position_residuals: Vec<f64>,
    /// Original normalized row-mobility condition at the accepted configuration.
    pub reciprocal_condition: f64,
    pub support_wrench_body: Wrench,
    pub support_wrench_world: Wrench,
    pub attachment_energy_j: f64,
    /// Report for the complete force/position residual, not a stability proof.
    pub nonlinear: NewtonReport,
    pub residual_evaluations: usize,
}

/// Solve supported static equilibrium with fixed ideal locks and gear ratios.
///
/// A seed may be off the affine manifold; original joint limits still apply to
/// every trial. No projection is used, and success requires independent force
/// AND native position admission. The result initializes the existing supported
/// constrained motion owner. Contact activation, backlash, varying ratios and
/// nonlinear loop closure are outside this constant-row model. Small-signal
/// analysis must retain the same constraints; the unconstrained response owner
/// cannot substitute for that constrained tangent.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    rows: &[JointConstraint],
    config: EquilibriumConfig<'_>,
    controls: JointEquilibriumControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedEquilibrium, EquilibriumError> {
    dispatch(
        model,
        seed,
        loads,
        bushings,
        rows,
        config,
        controls,
        None,
        FrameMode::Stationary,
        cancelled,
    )
}

/// The same repeatable smooth actuator callback as native supported motion.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_with_effort(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    rows: &[JointConstraint],
    config: EquilibriumConfig<'_>,
    controls: JointEquilibriumControls<'_>,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedEquilibrium, EquilibriumError> {
    dispatch(
        model,
        seed,
        loads,
        bushings,
        rows,
        config,
        controls,
        Some(&mut effort),
        FrameMode::Stationary,
        cancelled,
    )
}

/// Uniformly rotating relative equilibrium with the original symmetry admission.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_rotating(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    rows: &[JointConstraint],
    config: EquilibriumConfig<'_>,
    controls: JointEquilibriumControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedEquilibrium, EquilibriumError> {
    dispatch(
        model,
        seed,
        loads,
        bushings,
        rows,
        config,
        controls,
        None,
        FrameMode::UniformRotation,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn dispatch(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    rows: &[JointConstraint],
    config: EquilibriumConfig<'_>,
    controls: JointEquilibriumControls<'_>,
    effort: Option<&mut EffortLaw<'_>>,
    frame: FrameMode,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedEquilibrium, EquilibriumError> {
    poll(&mut cancelled)?;
    if [controls.reaction_scales, controls.position_scales]
        .iter()
        .any(|values| {
            values.len() != rows.len() || values.iter().any(|v| !v.is_finite() || *v <= 0.0)
        })
    {
        return Err(EquilibriumError::Invalid(
            "complete positive constraint scales required",
        ));
    }
    crate::articulated::motion::reconfigure::constraints::validate_condition_floor(
        controls.minimum_reciprocal_condition,
    )?;
    if matches!(frame, FrameMode::UniformRotation) {
        rotating::admit_frame(seed, loads, &mut cancelled)?;
    }
    // Empty rows retain the original solver's exact state, reports and work.
    if rows.is_empty() {
        let result = super::dispatch(
            model, seed, loads, bushings, config, effort, frame, cancelled,
        )?;
        return Ok(ConstrainedArticulatedEquilibrium {
            reaction_effort: vec![0.0; model.dof_count()],
            multipliers: Vec::new(),
            position_residuals: Vec::new(),
            reciprocal_condition: 1.0,
            state: result.state,
            required_effort: result.required_effort,
            applied_effort: result.applied_effort,
            joint_effort_residual: result.joint_effort_residual,
            scaled_force_residual_norm: result.scaled_residual_norm,
            support_wrench_body: result.support_wrench_body,
            support_wrench_world: result.support_wrench_world,
            attachment_energy_j: result.attachment_energy_j,
            nonlinear: result.nonlinear,
            residual_evaluations: result.residual_evaluations,
        });
    }
    match bushings {
        EquilibriumBushings::Diagonal(b) => solve(
            model, seed, loads, b, rows, config, controls, effort, frame, cancelled,
        ),
        EquilibriumBushings::Coupled(b) => solve(
            model, seed, loads, b, rows, config, controls, effort, frame, cancelled,
        ),
    }
}

struct Trial {
    native: Evaluation,
    multipliers: Vec<f64>,
    reaction: Vec<f64>,
    balance: Vec<f64>,
    scaled: Vec<f64>,
}

fn finite(value: f64) -> Result<f64, EquilibriumError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EquilibriumError::Invalid(
            "unrepresentable constrained equilibrium",
        ))
    }
}
fn scaled(value: f64, scale: f64) -> Result<f64, EquilibriumError> {
    let result = finite(value / scale)?;
    if value != 0.0 && result == 0.0 {
        return Err(EquilibriumError::Invalid(
            "constraint residual scaling erased a balance",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn solve<B: BushingLaw>(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    rows: &[JointConstraint],
    config: EquilibriumConfig<'_>,
    controls: JointEquilibriumControls<'_>,
    effort: Option<&mut EffortLaw<'_>>,
    frame: FrameMode,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedEquilibrium, EquilibriumError> {
    admit_inputs(model, seed, loads, bushings, config, frame, &mut cancelled)?;
    // Original row owner admits dimensions, coefficients, bounded inventory and
    // independence. Off-manifold trial coordinates are deliberately permitted.
    PreparedJointConstraints::prescribed_trial(
        model,
        &seed.joint_position,
        rows,
        controls.minimum_reciprocal_condition,
        &mut cancelled,
    )?;
    let n = model.dof_count();
    let dimension = n + rows.len();
    let evaluations = Cell::new(0usize);
    let cancel = RefCell::new(cancelled);
    let effort = RefCell::new(effort);
    let evaluate = |x: &[f64]| -> Result<Trial, EquilibriumError> {
        let mut cancel = cancel.borrow_mut();
        poll(&mut *cancel)?;
        let used = evaluations.get();
        if used >= config.maximum_evaluations {
            return Err(EquilibriumError::EvaluationBudget {
                used,
                maximum: config.maximum_evaluations,
            });
        }
        evaluations.set(used + 1);
        if x.len() != dimension {
            return Err(EquilibriumError::Invalid("constrained trial dimension"));
        }
        let state = PrescribedBaseMotionState {
            time_s: seed.time_s,
            base: seed.base,
            joint_position: x[..n]
                .iter()
                .zip(config.coordinate_scales)
                .zip(&seed.joint_position)
                .map(|((&v, &scale), &initial)| v.mul_add(scale, initial))
                .collect(),
            joint_velocity: vec![0.0; n],
        };
        let native = evaluate_native(
            model,
            state,
            loads,
            bushings,
            effort.borrow_mut().as_deref_mut(),
            &mut *cancel,
        )?;
        let multipliers = x[n..]
            .iter()
            .zip(controls.reaction_scales)
            .map(|(&v, &scale)| {
                let value = finite(v * scale)?;
                if v != 0.0 && value == 0.0 {
                    return Err(EquilibriumError::Invalid(
                        "reaction scale erased a multiplier",
                    ));
                }
                Ok(value)
            })
            .collect::<Result<Vec<_>, EquilibriumError>>()?;
        let mut reaction = Vec::with_capacity(n);
        let mut balance = Vec::with_capacity(n);
        let mut residual = Vec::with_capacity(dimension);
        for joint in 0..n {
            if joint % 64 == 0 {
                poll(&mut *cancel)?;
            }
            // Match the original constrained-motion owner's row order and FMA
            // arithmetic. Admit the ACTUAL representable reaction that will be
            // published/applied, not an unrounded hidden value.
            let mut sum = 0.0;
            for (row, &multiplier) in rows.iter().zip(&multipliers) {
                sum = row.coefficients[joint].mul_add(multiplier, sum);
            }
            reaction.push(finite(sum)?);
            let value = finite(
                (Dd::from_f64(native.required[joint])
                    - Dd::from_f64(native.applied[joint])
                    - Dd::from_f64(sum))
                .to_f64(),
            )?;
            balance.push(value);
            residual.push(scaled(value, config.effort_scales[joint])?);
        }
        let position = apply_rows(rows, &native.state.joint_position, &mut *cancel)?;
        for ((value, row), &scale) in position.iter().zip(rows).zip(controls.position_scales) {
            residual.push(scaled(finite(*value - row.position)?, scale)?);
        }
        Ok(Trial {
            native,
            multipliers,
            reaction,
            balance,
            scaled: residual,
        })
    };
    let (solution, mut report) = solve_scaled(
        dimension,
        config,
        |x| evaluate(x).map(|v| v.scaled),
        || (*cancel.borrow_mut())(),
    )?;
    let value = evaluate(&solution)?;
    let residual = fs_solver::norm2(&value.scaled);
    report.residual_norm = residual;
    if !residual.is_finite() || residual > config.residual_tolerance {
        report.converged = false;
        return Err(EquilibriumError::Unresolved(report));
    }
    let mut cancel = cancel.borrow_mut();
    let row_admission = admit_constraints(
        rows,
        &value.native.state.joint_position,
        &value.native.state.joint_velocity,
        &mut *cancel,
    )?;
    let row_mobility = PreparedJointConstraints::prescribed_trial(
        model,
        &value.native.state.joint_position,
        rows,
        controls.minimum_reciprocal_condition,
        &mut *cancel,
    )?;
    let (support_wrench_body, support_wrench_world, _) =
        support_wrenches(model, &value.native.state, value.native.root_wrench)?;
    poll(&mut *cancel)?;
    Ok(ConstrainedArticulatedEquilibrium {
        scaled_force_residual_norm: fs_solver::norm2(&value.scaled[..n]),
        state: value.native.state,
        required_effort: value.native.required,
        applied_effort: value.native.applied,
        multipliers: value.multipliers,
        reaction_effort: value.reaction,
        joint_effort_residual: value.balance,
        position_residuals: row_admission.position,
        reciprocal_condition: row_mobility.reciprocal_condition,
        support_wrench_body,
        support_wrench_world,
        attachment_energy_j: value.native.energy,
        nonlinear: report,
        residual_evaluations: evaluations.get(),
    })
}
