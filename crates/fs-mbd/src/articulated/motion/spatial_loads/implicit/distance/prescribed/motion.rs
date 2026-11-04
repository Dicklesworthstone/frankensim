//! Continuous material/joint closures on the actual prescribed support path.
//!
//! The original loop position solve and hidden-constraint impulse surround the
//! native supported midpoint. Original ABA mass actions, RNEA support recovery,
//! physical load owners and actual impulse events supply every reaction/account.

use super::*;
use crate::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBasePath, admit_stage, finish_step, prepare_step,
};
use crate::articulated::motion::spatial_loads::prescribed::{
    PrescribedSpatialMotionStep, evaluate_shared, finish,
    implicit::{ForcingLaw, MidpointEquations as SupportedMidpointEquations},
};

/// Accepted closed-mechanism motion and independently computed support exchange.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedBilateralLoopStep {
    /// Canonical continuation state, after the final physical tangency event.
    pub state_after: PrescribedBaseMotionState,
    /// Actual independently verified midpoint; the forcing law was last sampled here.
    pub midpoint_state: PrescribedBaseMotionState,
    /// Native smooth supported motion, including RNEA reaction and finite storage/loss.
    /// Its endpoint precedes the final constraint impulse.
    pub midpoint: PrescribedSpatialMotionStep,
    /// Actual summed initial reaction, including its independent support impulse/work.
    pub initial_reaction: PrescribedBaseImpulseResponse,
    /// Native final tangent reaction at the prescribed endpoint boundary.
    pub final_reaction: PrescribedBaseImpulseResponse,
    pub rows: Vec<BilateralRow>,
    /// Flattened material-row coefficients, in N s, before/after smooth motion.
    pub initial_loop_impulses_ns: Vec<f64>,
    pub final_loop_impulses_ns: Vec<f64>,
    /// Affine joint-row coefficients in original row units, before/after motion.
    pub initial_joint_multipliers: Vec<f64>,
    pub final_joint_multipliers: Vec<f64>,
    pub position_residual_before_m: Vec<f64>,
    pub position_residual_after_m: Vec<f64>,
    pub velocity_residual_before_mps: Vec<f64>,
    pub velocity_residual_after_mps: Vec<f64>,
    pub point_residual_after: Vec<PointClosureResidual>,
    pub joint_position_residuals_after: Vec<f64>,
    pub joint_velocity_residuals_after: Vec<f64>,
    pub initial_scaled_pivot: f64,
    pub final_scaled_pivot: f64,
    /// Original Newton report with independent actual-impulse verification.
    pub nonlinear: NewtonReport,
    pub residual_evaluations: usize,
    /// Actual per-link increments across both events and smooth motion, J.
    pub kinetic_energy_change: f64,
    /// Both actual applied constraint-impulse works, excluding support work, J.
    pub constraint_impulse_work: f64,
    /// Native smooth support work plus both native support-event works, J.
    pub support_work: f64,
    /// Native support reaction ON the tree about the world origin, [N m s, N s].
    pub support_impulse_world: Wrench,
    /// Actual material-row resultant, including finite-gap couples and rounding.
    /// These are constraint events, not environmental or support forcing.
    pub constraint_impulse_world: Wrench,
    /// Sum of polarized native increments; not subtraction of large endpoint totals.
    pub momentum_change_world: Wrench,
    /// Actual increment minus environmental, support and constraint impulses.
    /// Includes the native smooth path/quadrature defect; never repairs momentum.
    pub momentum_balance_residual_world: Wrench,
    /// Kinetic + stored-energy change + loss minus external, actuator, support
    /// and constraint work. Numerical path/time error remains visible, never heat.
    pub mechanical_work_residual: f64,
}

/// Advance distance/spherical closures and optional affine joint rows on a
/// fixed or driven support, retaining original diagonal passive bushings.
///
/// The caller supplies mutually consistent support pose, body twist and spatial
/// acceleration at the actual initial/midpoint/endpoint times. The root may be
/// mobile. Representable endpoint coordinates and the initial reaction impulse
/// solve position closure simultaneously. Every trial derives its midpoint rate
/// from (q1-q0)/h, including ABA, forcing and independent final verification.
/// This configuration chart preserves exact stationary coordinates without
/// clipping a solved rate; the original physical endpoint admission still applies.
/// A native endpoint impulse imposes
/// tangency without changing position. Both events and smooth RNEA independently
/// expose support force/couple and signed work. No finite-mass surrogate or
/// momentum-residual-derived support is used.
///
/// At most 32 independent scalar rows and 4096 links are admitted. Shared row
/// normalization and Newton--Krylov own the bounded solve; each trial adds one
/// native O(links) inverse-mass action and FK to the supported midpoint. Only the
/// row mobility is dense. Initial/endpoint position and physical tangency are
/// checked; midpoint/probe states need not be on the manifold. This is not a
/// branch enclosure, constraint activation, symplecticity or exact nonlinear
/// momentum/energy claim. Refine nonlinear geometry and the prescribed path.
#[allow(clippy::too_many_arguments)]
pub fn step_with_bushings(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    bushings: &[SpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBilateralLoopStep, BilateralLoopError> {
    step_evaluated(
        model, state, path, loads, loops, joint_rows, bushings, duration, config, None, cancelled,
    )
}

/// Resolve environmental link wrenches and internal effort with the complete
/// closed mechanism. The callback sees elapsed h/2 and the actual absolute clock,
/// support path and trial motion after initial reaction. Both output slices reset
/// from held loads. Smooth, repeatable off-manifold probes and staged history are
/// required. The independently rechecked accepted sample owns dynamics and work;
/// finalization never resamples it. All borrowed inputs remain unchanged on error.
#[allow(clippy::too_many_arguments)]
pub fn step_with_bushings_and_forcing(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    bushings: &[SpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    mut forcing: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBilateralLoopStep, BilateralLoopError> {
    step_evaluated(
        model,
        state,
        path,
        loads,
        loops,
        joint_rows,
        bushings,
        duration,
        config,
        Some(&mut forcing),
        cancelled,
    )
}

/// The same supported closed mechanism with full passive six-coordinate K/C.
#[allow(clippy::too_many_arguments)]
pub fn step_with_coupled_bushings_and_forcing(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    mut forcing: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBilateralLoopStep, BilateralLoopError> {
    step_evaluated(
        model,
        state,
        path,
        loads,
        loops,
        joint_rows,
        bushings,
        duration,
        config,
        Some(&mut forcing),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn step_evaluated<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    bushings: &[B],
    duration: f64,
    config: BilateralLoopConfig,
    forcing: Option<&mut ForcingLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedBilateralLoopStep, BilateralLoopError> {
    poll(&mut cancelled)?;
    validate_config(config.nonlinear)?;
    validate_bilateral_controls(config)?;
    if loops.is_empty() && joint_rows.is_empty()
        || loops.len() > MAX_BILATERAL_ROWS
        || joint_rows.len() > MAX_BILATERAL_ROWS
        || loops.iter().map(|c| c.rows()).sum::<usize>() + joint_rows.len() > MAX_BILATERAL_ROWS
    {
        return Err(BilateralLoopError::Invalid(
            "independent bilateral row count must be in 1..=32",
        ));
    }
    let initial = prepare_step(model, state, path, loads.held, duration, &mut cancelled)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    admit_loops(model, loops)?;
    admit_constraints(
        joint_rows,
        &state.joint_position,
        &state.joint_velocity,
        &mut cancelled,
    )?;
    let initial_load = evaluate_shared(
        model,
        state,
        &initial.kinematics,
        loads,
        bushings,
        &mut cancelled,
    )?;
    admit_stage(&initial_load)?;
    let before = shape_geometry(model, state, loops)?;
    check_position(&before, config, "initial")?;
    check_velocity(&before, config, "initial")?;
    let initial_mobility = prescribed_mobility(
        model,
        &state.joint_position,
        &before,
        joint_rows,
        config,
        "initial",
        &mut cancelled,
    )?;
    let n = model.dof_count();
    let loop_count = before.gaps.len();
    let count = loop_count + joint_rows.len();
    let mut scales = vec![config.nonlinear.angular_velocity_scale; n];
    for (link, part) in model.links().iter().enumerate() {
        if part.joint().kind() == JointKind::Prismatic
            && let Some(dof) = model.dof_index(link)
        {
            scales[dof] = config.nonlinear.linear_velocity_scale;
        }
    }
    let (impulse_scales, gap_scales, target) = constraint_scales(
        loop_count,
        joint_rows,
        &scales,
        &initial_mobility,
        duration,
        config,
    )?;
    let mut guess: Vec<_> = state
        .joint_velocity
        .iter()
        .zip(&scales)
        .map(|(v, s)| v / s)
        .collect();
    guess.resize(n + count, 0.0);
    let forcing_values = if forcing.is_some() {
        (
            loads.held.external_body_wrench.to_vec(),
            loads.held.generalized_force.to_vec(),
        )
    } else {
        (Vec::new(), Vec::new())
    };
    let forcing = RefCell::new((forcing, forcing_values));
    let cancelled = RefCell::new(cancelled);
    let evaluations = Cell::new(0usize);
    let charge = || -> Result<(), BilateralLoopError> {
        let used = evaluations.get();
        if used >= config.nonlinear.maximum_evaluations {
            return Err(ImplicitSpatialError::EvaluationBudget {
                used,
                maximum: config.nonlinear.maximum_evaluations,
            }
            .into());
        }
        evaluations.set(used + 1);
        Ok(())
    };
    let residual = |x: &[f64]| -> Result<Vec<f64>, BilateralLoopError> {
        let mut cancel = cancelled.borrow_mut();
        poll(&mut *cancel)?;
        charge()?;
        let (endpoint_position, midpoint_coordinates) =
            endpoint_chart(state, &x[..n], &scales, duration)?;
        let multipliers: Vec<_> = x[n..]
            .iter()
            .zip(&impulse_scales)
            .map(|(x, s)| x * s)
            .collect();
        let (body, joint) =
            assemble_impulses(model, &before, joint_rows, &multipliers, &mut *cancel)?;
        // One native mass action for the summed trial impulse. Unit columns and
        // nonlinear probes are not independently admitted as physical endpoints.
        let delta = prescribed_base_impulse_delta(model, &state.joint_position, &joint, &body)?;
        poll(&mut *cancel)?;
        let kicked = PrescribedBaseMotionState {
            joint_velocity: state
                .joint_velocity
                .iter()
                .zip(delta.velocity.joint)
                .map(|(old, jump)| old + jump)
                .collect(),
            ..state.clone()
        };
        let equations = SupportedMidpointEquations {
            model,
            state: &kicked,
            path,
            loads,
            bushings,
            duration,
            midpoint_time: initial.midpoint_time,
            scales: &scales,
        };
        let midpoint = {
            let mut forcing = forcing.borrow_mut();
            if let (Some(law), (wrenches, efforts)) = &mut *forcing {
                equations.evaluate_with_forcing(
                    &midpoint_coordinates,
                    &mut **law,
                    wrenches,
                    efforts,
                    None,
                    &mut *cancel,
                )?
            } else {
                equations.evaluate(&midpoint_coordinates, None, None, &mut *cancel)?
            }
        };
        endpoint_residual(
            model,
            state,
            path,
            loops,
            joint_rows,
            duration,
            &endpoint_position,
            midpoint.residual,
            &gap_scales,
            &mut *cancel,
        )
    };
    let (solution, mut report) =
        solve_loop_equations(guess, target, config.nonlinear, residual, &cancelled)?;
    let initial_impulse: Vec<_> = solution[n..]
        .iter()
        .zip(&impulse_scales)
        .map(|(x, s)| x * s)
        .collect();
    let mut cancel = cancelled.borrow_mut();
    let initial_reaction = apply_reaction(
        model,
        state,
        &before,
        joint_rows,
        &initial_impulse,
        &mut *cancel,
    )?;
    let kicked = PrescribedBaseMotionState {
        base: initial_reaction.base_after,
        joint_velocity: initial_reaction.joint_velocity_after.clone(),
        ..state.clone()
    };
    // Independently verify the ACTUAL native event before any endpoint publication.
    charge()?;
    let kicked_initial = prepare_step(model, &kicked, path, loads.held, duration, &mut *cancel)?;
    let equations = SupportedMidpointEquations {
        model,
        state: &kicked,
        path,
        loads,
        bushings,
        duration,
        midpoint_time: initial.midpoint_time,
        scales: &scales,
    };
    let (endpoint_position, midpoint_coordinates) =
        endpoint_chart(state, &solution[..n], &scales, duration)?;
    let actual_midpoint = {
        let mut forcing = forcing.borrow_mut();
        if let (Some(law), (wrenches, efforts)) = &mut *forcing {
            equations.evaluate_with_forcing(
                &midpoint_coordinates,
                &mut **law,
                wrenches,
                efforts,
                None,
                &mut *cancel,
            )?
        } else {
            equations.evaluate(&midpoint_coordinates, None, None, &mut *cancel)?
        }
    };
    let actual_residual = endpoint_residual(
        model,
        state,
        path,
        loops,
        joint_rows,
        duration,
        &endpoint_position,
        actual_midpoint.residual,
        &gap_scales,
        &mut *cancel,
    )?;
    report.residual_norm = fs_solver::norm2(&actual_residual);
    if !report.residual_norm.is_finite() || report.residual_norm > target {
        report.converged = false;
        return Err(ImplicitSpatialError::Unresolved(report).into());
    }
    let midpoint_state = actual_midpoint.stage.state.clone();
    let evaluated = finish_step(
        model,
        &kicked,
        path,
        loads.held,
        duration,
        kicked_initial,
        initial_load.stored_energy,
        actual_midpoint.stage,
        &mut *cancel,
    )?;
    let midpoint = finish(model, evaluated, loads, bushings, &mut *cancel)?;
    let endpoint = &midpoint.motion.state_after;
    let end = shape_geometry(model, endpoint, loops)?;
    check_position(&end, config, "endpoint")?;
    let final_mobility = prescribed_mobility(
        model,
        &endpoint.joint_position,
        &end,
        joint_rows,
        config,
        "endpoint",
        &mut *cancel,
    )?;
    let mut end_rates = end.rates.clone();
    end_rates.extend(apply_rows(
        joint_rows,
        &endpoint.joint_velocity,
        &mut *cancel,
    )?);
    let final_impulse = final_mobility.cancel_rates(&end_rates)?;
    let final_reaction = apply_reaction(
        model,
        endpoint,
        &end,
        joint_rows,
        &final_impulse,
        &mut *cancel,
    )?;
    let state_after = PrescribedBaseMotionState {
        base: final_reaction.base_after,
        joint_velocity: final_reaction.joint_velocity_after.clone(),
        ..endpoint.clone()
    };
    let actual = shape_geometry(model, &state_after, loops)?;
    check_position(&actual, config, "endpoint")?;
    check_velocity(&actual, config, "endpoint")?;
    let joint_after = admit_constraints(
        joint_rows,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut *cancel,
    )?;
    let smooth = &midpoint.motion;
    let constraint_impulse_work = sum(&[
        initial_reaction.applied_impulse_work,
        final_reaction.applied_impulse_work,
    ]);
    let support_work = sum(&[
        initial_reaction.support_work,
        smooth.support_work,
        final_reaction.support_work,
    ]);
    let kinetic_energy_change = sum(&[
        initial_reaction.kinetic_energy_change,
        smooth.kinetic_energy_change,
        final_reaction.kinetic_energy_change,
    ]);
    let mechanical_work_residual = sum(&[
        initial_reaction.kinetic_work_residual,
        midpoint.mechanical_work_residual,
        final_reaction.kinetic_work_residual,
    ]);
    let support_impulse_world = sum_wrenches(&[
        initial_reaction.support_impulse_world,
        Wrench::new(
            smooth.support_angular_impulse_world,
            smooth.support_linear_impulse_world,
        ),
        final_reaction.support_impulse_world,
    ])?;
    let constraint_impulse_world = sum_wrenches(&[
        initial_reaction.external_impulse_world,
        final_reaction.external_impulse_world,
    ])?;
    let momentum_change_world = sum_wrenches(&[
        Wrench::new(
            initial_reaction.angular_momentum_change_world,
            initial_reaction.linear_momentum_change_world,
        ),
        Wrench::new(
            smooth.angular_momentum_change_world,
            smooth.linear_momentum_change_world,
        ),
        Wrench::new(
            final_reaction.angular_momentum_change_world,
            final_reaction.linear_momentum_change_world,
        ),
    ])?;
    let momentum_balance_residual_world = sum_wrenches(&[
        Wrench::new(
            initial_reaction.angular_momentum_residual_world,
            initial_reaction.linear_momentum_residual_world,
        ),
        Wrench::new(
            smooth.angular_momentum_residual_world,
            smooth.linear_momentum_residual_world,
        ),
        Wrench::new(
            final_reaction.angular_momentum_residual_world,
            final_reaction.linear_momentum_residual_world,
        ),
    ])?;
    validate_slice(
        &[
            constraint_impulse_work,
            support_work,
            kinetic_energy_change,
            mechanical_work_residual,
        ],
        "supported_loop.work",
    )?;
    poll(&mut *cancel)?;
    Ok(PrescribedBilateralLoopStep {
        state_after,
        midpoint_state,
        midpoint,
        initial_reaction,
        final_reaction,
        rows: before.rows,
        initial_loop_impulses_ns: initial_impulse[..loop_count].to_vec(),
        final_loop_impulses_ns: final_impulse[..loop_count].to_vec(),
        initial_joint_multipliers: initial_impulse[loop_count..].to_vec(),
        final_joint_multipliers: final_impulse[loop_count..].to_vec(),
        position_residual_before_m: before.gaps,
        position_residual_after_m: actual.gaps,
        velocity_residual_before_mps: before.rates,
        velocity_residual_after_mps: actual.rates,
        point_residual_after: actual.points.iter().map(|p| p.residual).collect(),
        joint_position_residuals_after: joint_after.position,
        joint_velocity_residuals_after: joint_after.velocity,
        initial_scaled_pivot: initial_mobility.minimum_pivot,
        final_scaled_pivot: final_mobility.minimum_pivot,
        nonlinear: report,
        residual_evaluations: evaluations.get(),
        kinetic_energy_change,
        constraint_impulse_work,
        support_work,
        support_impulse_world,
        constraint_impulse_world,
        momentum_change_world,
        momentum_balance_residual_world,
        mechanical_work_residual,
    })
}

fn shape_geometry(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    loops: &[BilateralLoop],
) -> Result<Geometry, BilateralLoopError> {
    let shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    geometry(
        &shape,
        FreeFloatingBaseState::new(state.base.world_from_base, state.base.twist_body),
        loops,
    )
}

// The nonlinear variables are scaled endpoint displacements, not rates later
// projected onto a constraint. Form the representable configuration first and
// derive the secant rate before evaluating ANY dynamics or load. In particular,
// q1 == q0 has exactly zero rate even at a nonzero locked coordinate. This also
// keeps the derivative probes in the same floating chart as the accepted trial.
fn endpoint_chart(
    state: &PrescribedBaseMotionState,
    displacement: &[f64],
    scales: &[f64],
    duration: f64,
) -> Result<(Vec<f64>, Vec<f64>), BilateralLoopError> {
    let mut position = Vec::with_capacity(scales.len());
    let mut midpoint_coordinates = Vec::with_capacity(scales.len());
    for ((q, x), scale) in state.joint_position.iter().zip(displacement).zip(scales) {
        let endpoint = q + duration * (x * scale);
        position.push(endpoint);
        midpoint_coordinates.push(((endpoint - q) / duration) / scale);
    }
    validate_slice(&position, "supported_loop.endpoint_chart")?;
    validate_slice(&midpoint_coordinates, "supported_loop.midpoint_chart")?;
    Ok((position, midpoint_coordinates))
}

#[allow(clippy::too_many_arguments)]
fn endpoint_residual(
    model: &ArticulatedModel,
    initial: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    duration: f64,
    position: &[f64],
    mut residual: Vec<f64>,
    gap_scales: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, BilateralLoopError> {
    poll(cancelled)?;
    let endpoint = PrescribedBaseMotionState {
        time_s: initial.time_s + duration,
        base: path.endpoint,
        joint_position: position.to_vec(),
        joint_velocity: vec![0.0; position.len()],
    };
    let mut gaps = shape_geometry(model, &endpoint, loops)?.gaps;
    gaps.extend(joint_gaps(joint_rows, &endpoint.joint_position, cancelled)?);
    residual.extend(gaps.iter().zip(gap_scales).map(|(g, s)| g / s));
    validate_slice(&residual, "supported_loop.equations")?;
    poll(cancelled)?;
    Ok(residual)
}

fn apply_reaction(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    geometry: &Geometry,
    joint_rows: &[JointConstraint],
    multipliers: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PrescribedBaseImpulseResponse, BilateralLoopError> {
    poll(cancelled)?;
    let (body, joint) = assemble_impulses(model, geometry, joint_rows, multipliers, cancelled)?;
    let impulse = apply_prescribed_base_impulses(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &joint,
        &body,
    )?;
    poll(cancelled)?;
    Ok(impulse)
}

fn sum(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(Dd::ZERO, |s, x| s + Dd::from_f64(*x))
        .to_f64()
}

fn sum_wrenches(values: &[Wrench]) -> Result<Wrench, BilateralLoopError> {
    let mut sum = [Dd::ZERO; 6];
    for value in values {
        for (s, x) in sum.iter_mut().zip(value.to_array()) {
            *s = *s + Dd::from_f64(x);
        }
    }
    let s = sum.map(Dd::to_f64);
    validate_slice(&s, "supported_loop.momentum")?;
    Ok(Wrench::new(
        Vec3::new(s[0], s[1], s[2]),
        Vec3::new(s[3], s[4], s[5]),
    ))
}
