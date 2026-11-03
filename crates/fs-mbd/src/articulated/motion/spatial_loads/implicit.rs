//! Implicit finite-connection dynamics in the original tree's momentum coordinates.
//!
//! Only the midpoint joint rates and base angular rate are nonlinear unknowns.
//! All external spatial impulses are eliminated through the original locked
//! inertia map. Thus internal attachment forces never become external momentum
//! sources, including when many stiff connections share links.
//!
//! This is a second-order Lie midpoint scheme, not generalized-alpha, a rigid
//! constraint solver, or a variational integrator. Its linear oscillator limit
//! is the A-stable implicit midpoint rule, whose stability function tends to -1
//! in the stiff limit: stability does not resolve an oscillation or remove it.
//! Nonlinear geometry and finite solve error have measured work defects; no
//! unconditional energy/passivity or exact trajectory claim is made.
//! The linear rule and its A-stable, non-L-stable limit are also documented by
//! the original MOOSE implementation:
//! <https://mooseframework.inl.gov/docs/doxygen/moose/classImplicitMidpoint.html>.

/// Regular rigid distance and spherical/ball loops with position and velocity closure.
pub mod distance;

/// Constant ideal joint locks and gear relations with native carrier recoil.
pub mod constraints;

use super::super::{
    ArticulatedModel, BaseState, EffortLaw, EvaluatedStep, FreeFloatingBaseState,
    FreeFloatingMotionError, FreeFloatingMotionState, FreeFloatingMotionStep,
    MAX_FREE_FLOATING_STEP_LINKS, MomentumCoordinates, MomentumFrame, MomentumMap, Se3, StageLoad,
    Twist, Vec3, advance, forward_kinematics, free_floating_forward_dynamics,
    free_floating_momentum, momentum_from_kinematics, physical_increments, poll, rotate,
    sum_products, validate_external_wrenches, validate_forward_dynamics_header, validate_slice,
    validate_vec3,
};
use super::{SpatialMotionLoads, SpatialMotionStep, bushing, evaluate, finish_with_bushings};
use crate::articulated::motion::reconfigure::constraints::ConstraintError;
use crate::articulated::{ArticulatedError, JointKind};
use fs_ga::GaError;
use fs_math::dd::Dd;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::nonlinear::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use std::cell::{Cell, RefCell};

/// Maximum nonlinear residual evaluations, including differencing and final verification.
pub const MAX_IMPLICIT_SPATIAL_EVALUATIONS: usize = 16384;

/// Bounded nonlinear solve and explicit physical rate scales.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplicitSpatialConfig {
    /// Positive angular-rate scale in rad/s for base, revolute and helical rates.
    pub angular_velocity_scale: f64,
    /// Positive rate scale in m/s for prismatic coordinates.
    pub linear_velocity_scale: f64,
    /// Positive absolute 2-norm tolerance on the scaled midpoint equations.
    /// No relative initial-residual tolerance can weaken this acceptance test.
    pub residual_tolerance: f64,
    /// Positive maximum-component stencil size in scaled rate coordinates.
    /// This is a numerical Jacobian parameter, not a trajectory tolerance.
    pub difference_step: f64,
    /// Newton attempts in 1..=64, including rejected attempts.
    pub maximum_iterations: usize,
    /// Residual evaluations in 2..=16384. The accepted final stage is charged.
    pub maximum_evaluations: usize,
    /// Krylov vectors in 1..=64, reduced to the actual nonlinear dimension.
    /// At most eight restart cycles are permitted per Newton attempt.
    pub krylov_restart: usize,
}

impl Default for ImplicitSpatialConfig {
    fn default() -> Self {
        Self {
            angular_velocity_scale: 1.0,
            linear_velocity_scale: 1.0,
            residual_tolerance: 1.0e-10,
            difference_step: 1.0e-5,
            maximum_iterations: 24,
            maximum_evaluations: 2048,
            krylov_restart: 24,
        }
    }
}

/// Accepted finite attachment motion and the actual simultaneous solve.
#[derive(Debug, Clone, PartialEq)]
pub struct ImplicitSpatialStep {
    /// Exact accepted stage used by the force and work quadratures. Material
    /// consumers can check physical applicability without reconstructing it.
    pub midpoint_state: FreeFloatingMotionState,
    /// Original state, spatial impulses, attachment storage, loss and work defects.
    pub spatial: SpatialMotionStep,
    /// Native Newton--Krylov report, with an independently repeated final residual.
    pub nonlinear: NewtonReport,
    /// All actual nonlinear residual evaluations, including Jacobian stencils.
    pub residual_evaluations: usize,
    /// Actual final scaled equation 2-norm; at most `residual_tolerance`.
    pub scaled_residual_norm: f64,
    /// Midpoint-rate equation errors in each original joint's rad/s or m/s.
    pub joint_velocity_residual: Vec<f64>,
    /// Reconstructed minus trial midpoint angular rate in body coordinates, rad/s.
    pub angular_velocity_residual: Vec3,
}

/// No updated state is published on any refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum ImplicitSpatialError {
    /// The caller's numerical controls are malformed or exceed hard bounds.
    Invalid(&'static str),
    /// Original model, state, load, geometric, limit or cancellation refusal.
    Motion(FreeFloatingMotionError),
    /// Original affine-row, native mobility or ideal-reaction admission refused.
    Constraint(ConstraintError),
    /// The bounded residual work was exhausted, including on final verification.
    EvaluationBudget {
        /// Residual evaluations already consumed.
        used: usize,
        /// Caller-admitted maximum.
        maximum: usize,
    },
    /// The generic numerical Jacobian owner could not represent its stencil.
    Difference(&'static str),
    /// The generic nonlinear owner refused construction.
    Solver(NewtonError),
    /// The generic solver stopped with a callback fault or dimension mismatch.
    SolverProgress(fs_solver::nonlinear::SolverRunProgress),
    /// The simultaneous midpoint equations did not meet the actual tolerance.
    Unresolved(NewtonReport),
}

impl core::fmt::Display for ImplicitSpatialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit articulated connections: {self:?}")
    }
}
impl std::error::Error for ImplicitSpatialError {}
impl From<FreeFloatingMotionError> for ImplicitSpatialError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value)
    }
}
impl From<ConstraintError> for ImplicitSpatialError {
    fn from(value: ConstraintError) -> Self {
        Self::Constraint(value)
    }
}
impl From<ArticulatedError> for ImplicitSpatialError {
    fn from(value: ArticulatedError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<GaError> for ImplicitSpatialError {
    fn from(value: GaError) -> Self {
        Self::Motion(value.into())
    }
}

/// Advance a complete free-floating tree with simultaneously implicit spatial loads.
///
/// At a solution, `qh=q0+h*vh/2`, `vh=v0+h*ah/2`, and
/// `Rh=R0 exp(h*omega_h/2)`, with `omega_h` recovered from the same stage's
/// total world momenta. The original ABA evaluates `ah` with all connections,
/// world loads and body wrenches at this stage. Endpoints use the same midpoint
/// rates; original link inertias reconstruct the final base twist and pose.
///
/// Native Newton--Krylov and its fallible numerical-Jacobian adapter own the
/// solve. Each residual costs O(links+loads); extra Krylov storage is
/// O(restart*(joints+3)), with no dense generalized mass or Jacobian matrix.
/// Initial and final checks retain all native limits. A forbidden trial or
/// derivative stencil refuses rather than clipping a joint or changing a force.
/// Every evaluation and generic-solver boundary polls cancellation.
///
/// Smooth positive-length finite connections are required. Resolve changing
/// geometry and stiff oscillation phase by timestep refinement, even if a coarse
/// linear step is stable. This does not impose loop closure, contact or rigid
/// attachment constraints, identify a unique nonlinear branch, or provide an
/// adjoint. The supplied state and all loads remain unchanged on every error.
/// Only initial, trial-midpoint and endpoint geometry is evaluated; this does
/// not establish continuous-path admissibility or perform collision detection.
pub fn step_free_floating_spatial_implicit(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    step_with_bushings::<bushing::SpatialBushing>(
        model,
        state,
        loads,
        &[],
        duration,
        config,
        cancelled,
    )
}

/// Advance spatial loads and state-dependent internal joint efforts in one
/// simultaneous implicit Lie-midpoint solve.
///
/// Every residual and numerical-Jacobian probe evaluates `effort` at elapsed
/// time `duration / 2` and the complete reconstructed trial midpoint state.
/// Its output is reset to `loads.held.generalized_force` before every call;
/// the callback may add to or replace those efforts. Revolute/helical efforts
/// are torques in N m and prismatic efforts are forces in N. Native ABA gives
/// the equal/opposite carrier reaction. Internal effort is never included in
/// external spatial impulse, and the independently verified accepted midpoint
/// effort supplies the returned actuator-work quadrature.
///
/// This makes smooth torque-speed, joint spring/damper and feedback laws part
/// of the implicit equations, rather than holding their starting torque over
/// the step. The linear implicit-midpoint limit is A-stable but not L-stable:
/// a coarse stable step can retain unresolved oscillation or alternate a
/// stiff decaying mode. Refine to resolve response and geometry. An arbitrary
/// effort law is not assumed passive; its storage and dissipation are owned
/// by the caller, separately from the spatial connections' accounts.
///
/// The law must be repeatable for identical inputs throughout the solve and
/// must not advance accepted constitutive or controller history during trial
/// evaluations. Callback order/count are solver-dependent, including probes
/// and final verification. Native inputs remain immutable on any refusal;
/// caller-owned callback side effects cannot be rolled back. Cancellation is
/// polled before and after each call; long-running callbacks own their checks.
///
/// # Errors
/// Original implicit/motion/load refusals, invalid or over-limit trial effort,
/// cancellation, or a callback error. A forbidden trial is refused, not clipped.
pub fn step_free_floating_spatial_implicit_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    config: ImplicitSpatialConfig,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    step_with_bushings_and_effort::<bushing::SpatialBushing>(
        model,
        state,
        loads,
        &[],
        duration,
        config,
        Some(&mut effort),
        cancelled,
    )
}

pub(super) fn step_with_bushings<B: bushing::BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    step_with_bushings_and_effort(
        model, state, loads, bushings, duration, config, None, cancelled,
    )
}

pub(super) fn step_with_bushings_and_effort<B: bushing::BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    config: ImplicitSpatialConfig,
    effort: Option<&mut EffortLaw<'_>>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    step_with_bushings_effort_and_reaction(
        model, state, loads, bushings, duration, config, effort, None, cancelled,
    )
}

type ReactionLaw<'a> = dyn FnMut(
        &FreeFloatingMotionState,
        &mut [f64],
        &mut dyn FnMut() -> bool,
    ) -> Result<Dd, ImplicitSpatialError>
    + 'a;

#[allow(clippy::too_many_arguments)]
fn step_with_bushings_effort_and_reaction<B: bushing::BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    config: ImplicitSpatialConfig,
    effort: Option<&mut EffortLaw<'_>>,
    reaction: Option<&mut ReactionLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, ImplicitSpatialError> {
    poll(&mut cancelled)?;
    if !duration.is_finite() || duration <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidDuration.into());
    }
    validate_config(config)?;
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    model.free_floating_complexity()?;
    validate_forward_dynamics_header(
        model,
        state.base.twist_body,
        &state.joint_position,
        &state.joint_velocity,
        loads.held.generalized_force,
        loads.held.gravity_world,
        loads.held.external_body_wrench,
    )?;
    validate_external_wrenches(loads.held.external_body_wrench)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    let initial_map = MomentumMap::new(model, state, &mut cancelled)?;
    let initial = initial_map.coordinates(state.base)?;
    let initial_load =
        bushing::evaluate(&initial_map, state.base, loads, bushings, &mut cancelled)?.load;
    let momentum_before = free_floating_momentum(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    let mut scales = vec![config.angular_velocity_scale; model.dof_count() + 3];
    for (link, part) in model.links().iter().enumerate() {
        if part.joint.kind() == JointKind::Prismatic
            && let Some(dof) = model.dof_index(link)
        {
            scales[dof] = config.linear_velocity_scale;
        }
    }
    let n = model.dof_count();
    let mut guess = state.joint_velocity.clone();
    let omega = state.base.twist_body.angular;
    guess.extend([omega.x, omega.y, omega.z]);
    for (value, scale) in guess.iter_mut().zip(&scales) {
        *value /= scale;
    }
    let gravity_only_external = loads.world_points.is_empty()
        && loads
            .held
            .external_body_wrench
            .iter()
            .all(|wrench| wrench.force == Vec3::default() && wrench.torque == Vec3::default());
    let equations = MidpointEquations {
        model,
        state,
        loads,
        bushings,
        duration,
        initial,
        mass: initial_map.mass(),
        scales: &scales,
        gravity_only_external,
    };
    let evaluations = Cell::new(0usize);
    let cancelled = RefCell::new(cancelled);
    // One O(joints) output buffer is reused across all nonlinear and stencil
    // evaluations. The held-effort entry points need no extra allocation.
    let effort_values = if effort.is_some() {
        loads.held.generalized_force.to_vec()
    } else {
        Vec::new()
    };
    let effort = RefCell::new((effort, effort_values));
    let reaction = RefCell::new(reaction);
    let evaluate_stage = |x: &[f64]| {
        let mut cancel = cancelled.borrow_mut();
        poll(&mut *cancel)?;
        let used = evaluations.get();
        if used >= config.maximum_evaluations {
            return Err(ImplicitSpatialError::EvaluationBudget {
                used,
                maximum: config.maximum_evaluations,
            });
        }
        evaluations.set(used + 1);
        let mut effort = effort.borrow_mut();
        let (law, values) = &mut *effort;
        let stage_effort = law.as_deref_mut().map(|law| (law, values.as_mut_slice()));
        let mut reaction = reaction.borrow_mut();
        equations.evaluate_with_effort_and_reaction(
            x,
            stage_effort,
            reaction.as_deref_mut(),
            &mut *cancel,
        )
    };
    let problem = FallibleFiniteDifference::new(n + 3, config.difference_step, |x: &[f64]| {
        evaluate_stage(x).map(|stage| stage.residual)
    })
    .map_err(difference_error)?;
    let solver_config = NewtonKrylovConfig {
        absolute_tolerance: config.residual_tolerance,
        relative_tolerance: 0.0,
        linear_restart: config.krylov_restart.min(n + 3),
        max_linear_cycles: 8,
        // Numerical directional differences have their own rounding floor.
        // Do not demand an over-solved inner action; every outer residual still
        // must independently meet the caller's strict absolute tolerance.
        forcing_minimum: 1.0e-3,
        forcing_maximum: 0.1,
        ..NewtonKrylovConfig::default()
    };
    let construction = NewtonKrylovState::new(&problem, guess, solver_config);
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = construction.map_err(ImplicitSpatialError::Solver)?;
    let result = solver.run_controlled(&problem, config.maximum_iterations, &mut || {
        problem.has_failed() || (*cancelled.borrow_mut())()
    });
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    match result.progress {
        fs_solver::nonlinear::SolverRunProgress::Complete => {}
        fs_solver::nonlinear::SolverRunProgress::Paused => {
            return Err(FreeFloatingMotionError::Cancelled.into());
        }
        other => return Err(ImplicitSpatialError::SolverProgress(other)),
    }
    let mut report = result.report;
    if !report.converged {
        return Err(ImplicitSpatialError::Unresolved(report));
    }
    let midpoint = evaluate_stage(&solver.x)?;
    let residual_norm = fs_solver::norm2(&midpoint.residual);
    report.residual_norm = residual_norm;
    if !residual_norm.is_finite() || residual_norm > config.residual_tolerance {
        report.converged = false;
        return Err(ImplicitSpatialError::Unresolved(report));
    }
    let joint_velocity_residual = midpoint.residual[..n]
        .iter()
        .zip(&scales)
        .map(|(r, s)| r * s)
        .collect();
    let angular_velocity_residual = Vec3::new(
        -midpoint.residual[n],
        -midpoint.residual[n + 1],
        -midpoint.residual[n + 2],
    )
    .scale(config.angular_velocity_scale);
    let mut cancel = cancelled.borrow_mut();
    let midpoint_state = midpoint.state.clone();
    let spatial = finish_midpoint(
        model,
        state,
        loads,
        duration,
        MidpointInitial {
            map: initial_map,
            coordinates: initial,
            load: initial_load,
            momentum: momentum_before,
        },
        midpoint,
        bushings,
        &mut *cancel,
    )?;
    Ok(ImplicitSpatialStep {
        midpoint_state,
        spatial,
        nonlinear: report,
        residual_evaluations: evaluations.get(),
        scaled_residual_norm: residual_norm,
        joint_velocity_residual,
        angular_velocity_residual,
    })
}

struct MidpointInitial {
    map: MomentumMap,
    coordinates: MomentumCoordinates,
    load: StageLoad,
    momentum: crate::articulated::ArticulatedMomentum,
}

fn finish_midpoint<B: bushing::BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    initial: MidpointInitial,
    midpoint: MidpointEvaluation,
    bushings: &[B],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialMotionStep, ImplicitSpatialError> {
    let MidpointInitial {
        map: initial_map,
        coordinates: initial,
        load: initial_load,
        momentum: momentum_before,
    } = initial;
    let load = midpoint.load;
    let linear_impulse = load.force_world.scale(duration);
    let external_impulse = load.external_impulse(
        midpoint.state.base.world_from_base,
        midpoint.center_base,
        duration,
    )?;
    let center_displacement =
        (initial.linear_world + linear_impulse.scale(0.5)).scale(duration / initial_map.mass());
    let final_coordinates = MomentumCoordinates {
        center_world: initial.center_world + center_displacement,
        linear_world: initial.linear_world + linear_impulse,
        angular_com_world: initial.angular_com_world + load.torque_com_world.scale(duration),
    };
    let mut state_after = FreeFloatingMotionState {
        base: state.base,
        joint_position: advance(
            &state.joint_position,
            &midpoint.state.joint_velocity,
            duration,
        ),
        joint_velocity: midpoint
            .state
            .joint_velocity
            .iter()
            .zip(&state.joint_velocity)
            .map(|(mid, old)| 2.0 * mid - old)
            .collect(),
    };
    let final_rotation = rotate(
        state.base.world_from_base.rotation(),
        midpoint.trial_angular,
        duration,
    )?;
    let final_map = MomentumMap::new(model, &state_after, cancelled)?;
    state_after.base = final_map.reconstruct(final_rotation, final_coordinates)?;
    poll(cancelled)?;
    let before = forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(cancelled)?;
    let after = forward_kinematics(
        model,
        BaseState::prescribed(
            state_after.base.world_from_base,
            state_after.base.twist_body,
            Twist::zero(),
        ),
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    let momentum_after = momentum_from_kinematics(model, &after)?;
    let (kinetic_change, _, momentum_change) = physical_increments(
        model,
        MomentumFrame {
            world: &before,
            relative: &initial_map.shape,
            world_from_base: state.base.world_from_base,
        },
        MomentumFrame {
            world: &after,
            relative: &final_map.shape,
            world_from_base: state_after.base.world_from_base,
        },
        Vec3::default(),
        cancelled,
    )?;
    let momentum_residual = momentum_change.subtract(external_impulse);
    let external_impulse = external_impulse.checked_wrench("implicit_spatial.impulse")?;
    let momentum_change = momentum_change.checked_wrench("implicit_spatial.change")?;
    let momentum_residual = momentum_residual.checked_wrench("implicit_spatial.residual")?;
    let dt = Dd::from_f64(duration);
    let actuator_work_dd = dt * midpoint.actuator_power;
    let constraint_work_dd = dt * midpoint.constraint_power;
    let gravity = loads.held.gravity_world;
    let gravity_work_dd = Dd::from_f64(initial_map.mass())
        * sum_products(
            [gravity.x, gravity.y, gravity.z],
            [
                center_displacement.x,
                center_displacement.y,
                center_displacement.z,
            ],
        );
    let external_work_dd = dt * Dd::from_f64(load.body_wrench_power) + gravity_work_dd;
    let internal_body_work_dd = dt * Dd::from_f64(load.internal_body_power);
    let actuator_work = actuator_work_dd.to_f64();
    let constraint_work = constraint_work_dd.to_f64();
    let external_work = external_work_dd.to_f64();
    let internal_body_work = internal_body_work_dd.to_f64();
    let kinetic_work_residual = (kinetic_change
        - actuator_work_dd
        - constraint_work_dd
        - external_work_dd
        - internal_body_work_dd)
        .to_f64();
    validate_vec3(linear_impulse, "implicit_spatial.linear_impulse")?;
    validate_slice(
        &[
            kinetic_change.to_f64(),
            actuator_work,
            constraint_work,
            external_work,
            internal_body_work,
            kinetic_work_residual,
            duration * load.dissipation_power,
        ],
        "implicit_spatial.work",
    )?;
    let spatial = finish_with_bushings(
        model,
        EvaluatedStep {
            stored_energy_before: initial_load.stored_energy,
            dissipated_energy: duration * load.dissipation_power,
            step: FreeFloatingMotionStep {
                state_after,
                center_of_mass_world_before: initial.center_world,
                center_of_mass_world_after: final_coordinates.center_world,
                momentum_before,
                momentum_after,
                kinetic_energy_change: kinetic_change.to_f64(),
                linear_momentum_change_world: momentum_change.force,
                angular_momentum_change_world: momentum_change.torque,
                external_linear_impulse_world: external_impulse.force,
                external_angular_impulse_world: external_impulse.torque,
                linear_momentum_residual_world: momentum_residual.force,
                angular_momentum_residual_world: momentum_residual.torque,
                actuator_work,
                constraint_work,
                external_work,
                internal_body_work,
                kinetic_work_residual,
            },
        },
        loads,
        bushings,
        cancelled,
    )?;
    Ok(spatial)
}

struct MidpointEquations<'a, B: bushing::BushingLaw> {
    model: &'a ArticulatedModel,
    state: &'a FreeFloatingMotionState,
    loads: SpatialMotionLoads<'a>,
    bushings: &'a [B],
    duration: f64,
    initial: MomentumCoordinates,
    mass: f64,
    scales: &'a [f64],
    gravity_only_external: bool,
}

struct MidpointEvaluation {
    state: FreeFloatingMotionState,
    center_base: Vec3,
    trial_angular: Vec3,
    load: StageLoad,
    actuator_power: Dd,
    constraint_power: Dd,
    residual: Vec<f64>,
}

impl<B: bushing::BushingLaw> MidpointEquations<'_, B> {
    fn evaluate(
        &self,
        x: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MidpointEvaluation, ImplicitSpatialError> {
        self.evaluate_with_effort(x, None, cancelled)
    }

    fn evaluate_with_effort(
        &self,
        x: &[f64],
        effort: Option<(&mut EffortLaw<'_>, &mut [f64])>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MidpointEvaluation, ImplicitSpatialError> {
        self.evaluate_with_effort_and_reaction(x, effort, None, cancelled)
    }

    fn evaluate_with_effort_and_reaction(
        &self,
        x: &[f64],
        effort: Option<(&mut EffortLaw<'_>, &mut [f64])>,
        reaction: Option<&mut ReactionLaw<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MidpointEvaluation, ImplicitSpatialError> {
        let n = self.model.dof_count();
        let half = 0.5 * self.duration;
        let velocity: Vec<_> = x[..n].iter().zip(self.scales).map(|(x, s)| x * s).collect();
        let omega = Vec3::new(
            x[n] * self.scales[n],
            x[n + 1] * self.scales[n + 1],
            x[n + 2] * self.scales[n + 2],
        );
        let rotation = rotate(self.state.base.world_from_base.rotation(), omega, half)?;
        let mut state = FreeFloatingMotionState {
            base: self.state.base,
            joint_position: advance(&self.state.joint_position, &velocity, half),
            joint_velocity: velocity,
        };
        let map = MomentumMap::new(self.model, &state, cancelled)?;
        // External resultants depend on configuration, not the trial velocity or
        // base origin. Eliminate their midpoint momentum equations exactly before
        // computing physical work or any connection's actual strain rate.
        let (force_world, torque_com_world) = if self.gravity_only_external {
            (
                self.loads.held.gravity_world.scale(self.mass),
                Vec3::default(),
            )
        } else {
            let probe =
                FreeFloatingBaseState::stationary(Se3::from_parts(rotation, Vec3::default())?);
            let external = evaluate(
                &map,
                probe,
                SpatialMotionLoads {
                    connections: &[],
                    ..self.loads
                },
                cancelled,
            )?
            .load;
            (external.force_world, external.torque_com_world)
        };
        let linear_world = self.initial.linear_world + force_world.scale(half);
        let coordinates = MomentumCoordinates {
            center_world: self.initial.center_world + linear_world.scale(half / self.mass),
            linear_world,
            angular_com_world: self.initial.angular_com_world + torque_com_world.scale(half),
        };
        state.base = map.reconstruct(rotation, coordinates)?;
        let stage = bushing::evaluate(&map, state.base, self.loads, self.bushings, cancelled)?;
        poll(cancelled)?;
        let generalized_force = if let Some((law, values)) = effort {
            values.copy_from_slice(self.loads.held.generalized_force);
            law(half, &state, values)?;
            poll(cancelled)?;
            &*values
        } else {
            self.loads.held.generalized_force
        };
        let mut dynamics = free_floating_forward_dynamics(
            self.model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
            generalized_force,
            self.loads.held.gravity_world,
            &stage.body_wrench,
        )?;
        let constraint_power = if let Some(law) = reaction {
            poll(cancelled)?;
            let power = law(&state, &mut dynamics.generalized_acceleration, cancelled)?;
            poll(cancelled)?;
            validate_slice(
                &dynamics.generalized_acceleration,
                "implicit_spatial.constrained_acceleration",
            )?;
            power
        } else {
            Dd::ZERO
        };
        let mut residual: Vec<_> = state
            .joint_velocity
            .iter()
            .zip(&self.state.joint_velocity)
            .zip(&dynamics.generalized_acceleration)
            .zip(self.scales)
            .map(|(((mid, old), a), s)| (mid - old - half * a) / s)
            .collect();
        let angular = (omega - state.base.twist_body.angular).scale(1.0 / self.scales[n]);
        residual.extend([angular.x, angular.y, angular.z]);
        validate_slice(&residual, "implicit_spatial.equations")?;
        let actuator_power = sum_products(
            generalized_force.iter().copied(),
            state.joint_velocity.iter().copied(),
        );
        validate_slice(
            &[actuator_power.to_f64(), constraint_power.to_f64()],
            "implicit_spatial.stage_power",
        )?;
        poll(cancelled)?;
        Ok(MidpointEvaluation {
            state,
            center_base: map.center_base,
            trial_angular: omega,
            load: stage.load,
            actuator_power,
            constraint_power,
            residual,
        })
    }
}

fn validate_config(config: ImplicitSpatialConfig) -> Result<(), ImplicitSpatialError> {
    if let Some(reason) = config_error(config) {
        Err(ImplicitSpatialError::Invalid(reason))
    } else {
        Ok(())
    }
}

/// Shared hard bounds for free and prescribed-base nonlinear midpoint solves.
pub(super) fn config_error(config: ImplicitSpatialConfig) -> Option<&'static str> {
    if [
        config.angular_velocity_scale,
        config.linear_velocity_scale,
        config.residual_tolerance,
        config.difference_step,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Some("finite positive rate scales and tolerances required");
    }
    if config.maximum_iterations == 0
        || config.maximum_iterations > 64
        || config.maximum_evaluations < 2
        || config.maximum_evaluations > MAX_IMPLICIT_SPATIAL_EVALUATIONS
        || config.krylov_restart == 0
        || config.krylov_restart > 64
    {
        return Some("nonlinear work controls exceed admitted bounds");
    }
    None
}

fn difference_error(error: FiniteDifferenceError<ImplicitSpatialError>) -> ImplicitSpatialError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => ImplicitSpatialError::Difference(reason),
    }
}
