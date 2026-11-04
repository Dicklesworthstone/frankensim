//! Implicit midpoint dynamics for fixed and prescribed moving supports.
//!
//! The nonlinear unknowns are the original scalar joint rates. Native ABA
//! supplies their acceleration at the same trial stage as every world load,
//! finite connection and feedback effort. The shared supported-motion owner
//! then recovers native RNEA support reactions, impulses and work.
//!
//! For a resolved linear oscillator this is the A-stable implicit midpoint
//! rule. It is not L-stable: very stiff unresolved oscillations are not filtered.
//! Resolve phase and nonlinear geometry by timestep refinement. The supplied
//! base samples must describe one smooth physical path; no derivative
//! consistency, continuous-path admission, rigid closure or contact is inferred.

/// Constant affine joint locks and gear relations on the same implicit stage.
pub mod constraints;

use super::super::SpatialMotionLoads;
use super::super::bushing::{self, BushingLaw, SpatialBushing};
pub use super::super::implicit::ImplicitSpatialConfig;
use super::super::implicit::config_error;
use super::{PrescribedSpatialMotionStep, evaluate_shared, finish};
use crate::articulated::motion::prescribed::{
    EffortLaw, PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
    PrescribedMidpoint, PrescribedStage, admit_stage, finish_step, poll, prepare_step,
};
use crate::articulated::motion::reconfigure::constraints::ConstraintError;
use crate::articulated::motion::{FreeFloatingMotionError, advance, sum_products};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, JointKind, Kinematics, forward_dynamics,
    forward_kinematics, validate_slice,
};
use fs_ga::Wrench;
use fs_math::dd::Dd;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::nonlinear::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use std::borrow::Cow;
use std::cell::{Cell, RefCell};

pub(in crate::articulated::motion::spatial_loads) type ForcingLaw<'a> = dyn FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>
    + 'a;

/// An accepted supported trajectory and its actual implicit midpoint solve.
#[derive(Debug, Clone, PartialEq)]
pub struct ImplicitPrescribedSpatialStep {
    /// Exact accepted stage used for dynamics, loads, effort and work quadrature.
    pub midpoint_state: PrescribedBaseMotionState,
    /// Original supported motion, support impulses/work and connection accounts.
    pub spatial: PrescribedSpatialMotionStep,
    /// Native Newton--Krylov report with an independently repeated final residual.
    /// An all-fixed tree has no unknowns and completes with zero Newton attempts.
    pub nonlinear: NewtonReport,
    /// All residual evaluations, including numerical stencils and final verification.
    pub residual_evaluations: usize,
    /// Actual final scaled 2-norm, no greater than the admitted absolute tolerance.
    pub scaled_residual_norm: f64,
    /// Midpoint joint-rate equation errors in their original rad/s or m/s units.
    pub joint_velocity_residual: Vec<f64>,
}

/// A refusal publishes no endpoint or partially updated support state.
#[derive(Debug, Clone, PartialEq)]
pub enum ImplicitPrescribedSpatialError {
    /// Malformed or out-of-bounds numerical controls.
    Invalid(&'static str),
    /// Original supported dynamics, load, geometry, clock or cancellation refusal.
    Motion(PrescribedBaseMotionError),
    /// Original affine joint-row admission, mobility, reaction or work refusal.
    Constraint(ConstraintError),
    /// Residual work was exhausted, including on independent final verification.
    EvaluationBudget {
        /// Residual evaluations already consumed.
        used: usize,
        /// Caller-admitted maximum.
        maximum: usize,
    },
    /// The native numerical-Jacobian owner refused its stencil.
    Difference(&'static str),
    /// The native nonlinear owner refused construction.
    Solver(NewtonError),
    /// The generic solver stopped with a callback fault or dimension mismatch.
    SolverProgress(fs_solver::nonlinear::SolverRunProgress),
    /// Actual simultaneous midpoint equations did not meet the required tolerance.
    Unresolved(NewtonReport),
}

impl core::fmt::Display for ImplicitPrescribedSpatialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit prescribed-base connections: {self:?}")
    }
}
impl std::error::Error for ImplicitPrescribedSpatialError {}
impl From<PrescribedBaseMotionError> for ImplicitPrescribedSpatialError {
    fn from(value: PrescribedBaseMotionError) -> Self {
        Self::Motion(value)
    }
}
impl From<ConstraintError> for ImplicitPrescribedSpatialError {
    fn from(value: ConstraintError) -> Self {
        Self::Constraint(value)
    }
}
impl From<ArticulatedError> for ImplicitPrescribedSpatialError {
    fn from(value: ArticulatedError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<FreeFloatingMotionError> for ImplicitPrescribedSpatialError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value.into())
    }
}

/// Advance a supported mechanism with simultaneously implicit spatial loads.
///
/// The solved equations are qh=q0+h*vh/2 and vh=v0+h*ah/2. Native prescribed-base
/// ABA evaluates ah using the caller's midpoint base pose, twist and acceleration.
/// Endpoints use q1=q0+h*vh and v1=v0+h*ah. Original link inertias, scalar root
/// joints, held/world forces and axial laws are retained without a surrogate base.
///
/// Native Newton--Krylov and fallible directional finite differences own the solve.
/// Each residual costs O(links+loads); Krylov storage is O(restart*joints), with
/// no dense mass or Jacobian matrix. Every evaluation consumes the bounded work
/// budget and polls cancellation. Forbidden trial or endpoint geometry refuses
/// rather than clipping coordinates. Inputs remain unchanged on every error.
///
/// Initial, trial-midpoint and endpoint attachment geometry are admitted. These
/// sampled checks do not prove the continuous path stays in an admissible chart.
/// No unique nonlinear branch or convergence is certified. Exact energy/momentum
/// conservation is not promised; the shared receipt exposes actual physical work
/// and momentum defects.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_spatial_implicit(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_prescribed_base_bushings_implicit(
        model,
        state,
        path,
        loads,
        &[],
        duration,
        config,
        cancelled,
    )
}

/// Solve smooth joint feedback together with all spatial loads at one midpoint.
///
/// Every callback receives elapsed time h/2 and the actual trial stage with its
/// absolute clock. Its effort buffer is reset to held generalized force before
/// each invocation. Newton, derivative stencils and final verification may call
/// the law repeatedly in nonmonotone trial order. The law must be deterministic
/// and smooth on admitted trials; callback-owned effects cannot be rolled back.
///
/// Native effort/rate/configuration limits and cancellation apply to every trial.
/// A changing or failed law cannot publish an unchecked final endpoint.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_spatial_implicit_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    config: ImplicitSpatialConfig,
    effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_prescribed_base_bushings_implicit_with_effort(
        model,
        state,
        path,
        loads,
        &[],
        duration,
        config,
        effort,
        cancelled,
    )
}

/// Solve finite six-axis bushings and axial connections on a prescribed support.
///
/// The original bushing owner supplies finite translation, principal-log rotation,
/// exact SO(3) pullback, elastic storage and damping. Its forces share one native
/// implicit stage, support reaction and work account with every other load.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_implicit(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_with_bushings_and_effort(
        model, state, path, loads, bushings, duration, config, None, cancelled,
    )
}

/// Solve the same finite bushings with jointly implicit smooth effort feedback.
///
/// Trial callbacks follow the repeated-evaluation contract of the spatial effort
/// entry. Diagonal and coupled bushing adapters share this same nonlinear owner.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_implicit_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_with_bushings_and_effort(
        model,
        state,
        path,
        loads,
        bushings,
        duration,
        config,
        Some(&mut effort),
        cancelled,
    )
}

/// Resolve environmental link wrenches and internal effort at the same actual
/// implicit stage as the prescribed support, spatial loads and finite bushings.
///
/// The callback receives elapsed time h/2; the state carries the absolute time.
/// Both output slices are reset from held loads at every numerical trial and
/// independent final verification. Wrenches are about original link origins in
/// their body axes. The last accepted sample supplies ABA, native support RNEA,
/// momentum and work; it is not sampled again during finalization. Forces are
/// finite and subject to native admission; no support reaction is inferred from
/// a balance residual. Use a smooth repeatable law, stage callback-owned history
/// until success, and resolve nonlinear motion and any load discontinuities.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_implicit_with_forcing(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    config: ImplicitSpatialConfig,
    mut forcing: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [Wrench],
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_with_bushings_forcing_and_reaction(
        model,
        state,
        path,
        loads,
        bushings,
        duration,
        config,
        None,
        Some(&mut forcing),
        None,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(in crate::articulated::motion::spatial_loads) fn step_with_bushings_and_effort<
    B: BushingLaw,
>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    config: ImplicitSpatialConfig,
    effort: Option<&mut EffortLaw<'_>>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_with_bushings_effort_and_reaction(
        model, state, path, loads, bushings, duration, config, effort, None, cancelled,
    )
}

/// Additional native joint acceleration and its separately owned physical power.
/// Caller efforts and their native effort admission are never replaced or widened.
pub(in crate::articulated::motion::spatial_loads) type ReactionLaw<'a> = dyn FnMut(
        &PrescribedBaseMotionState,
        &mut [f64],
        &mut dyn FnMut() -> bool,
    ) -> Result<Dd, ImplicitPrescribedSpatialError>
    + 'a;

#[allow(clippy::too_many_arguments)]
fn step_with_bushings_effort_and_reaction<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    config: ImplicitSpatialConfig,
    effort: Option<&mut EffortLaw<'_>>,
    reaction: Option<&mut ReactionLaw<'_>>,
    cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    step_with_bushings_forcing_and_reaction(
        model, state, path, loads, bushings, duration, config, effort, None, reaction, cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn step_with_bushings_forcing_and_reaction<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    config: ImplicitSpatialConfig,
    effort: Option<&mut EffortLaw<'_>>,
    forcing: Option<&mut ForcingLaw<'_>>,
    reaction: Option<&mut ReactionLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ImplicitPrescribedSpatialStep, ImplicitPrescribedSpatialError> {
    poll(&mut cancelled)?;
    if let Some(reason) = config_error(config) {
        return Err(ImplicitPrescribedSpatialError::Invalid(reason));
    }
    let initial = prepare_step(model, state, path, loads.held, duration, &mut cancelled)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    let initial_load = evaluate_shared(
        model,
        state,
        &initial.kinematics,
        loads,
        bushings,
        &mut cancelled,
    )?;
    admit_stage(&initial_load)?;
    let n = model.dof_count();
    let mut scales = vec![config.angular_velocity_scale; n];
    for (link, part) in model.links().iter().enumerate() {
        if part.joint().kind() == JointKind::Prismatic
            && let Some(dof) = model.dof_index(link)
        {
            scales[dof] = config.linear_velocity_scale;
        }
    }
    let guess: Vec<_> = state
        .joint_velocity
        .iter()
        .zip(&scales)
        .map(|(rate, scale)| rate / scale)
        .collect();
    let equations = MidpointEquations {
        model,
        state,
        path,
        loads,
        bushings,
        duration,
        midpoint_time: initial.midpoint_time,
        scales: &scales,
    };
    let evaluations = Cell::new(0usize);
    let cancelled = RefCell::new(cancelled);
    // Reuse one effort buffer across every trial, numerical stencil and recheck.
    let effort_values = if effort.is_some() {
        loads.held.generalized_force.to_vec()
    } else {
        Vec::new()
    };
    let effort = RefCell::new((effort, effort_values));
    let forcing_values = if forcing.is_some() {
        (
            loads.held.external_body_wrench.to_vec(),
            loads.held.generalized_force.to_vec(),
        )
    } else {
        (Vec::new(), Vec::new())
    };
    let forcing = RefCell::new((forcing, forcing_values));
    let reaction = RefCell::new(reaction);
    let evaluate_stage = |x: &[f64]| {
        let mut cancel = cancelled.borrow_mut();
        poll(&mut *cancel)?;
        let used = evaluations.get();
        if used >= config.maximum_evaluations {
            return Err(ImplicitPrescribedSpatialError::EvaluationBudget {
                used,
                maximum: config.maximum_evaluations,
            });
        }
        evaluations.set(used + 1);
        let mut effort = effort.borrow_mut();
        let (law, values) = &mut *effort;
        let stage_effort = law.as_deref_mut().map(|law| (law, values.as_mut_slice()));
        let mut reaction = reaction.borrow_mut();
        let mut forcing = forcing.borrow_mut();
        if let (Some(law), (wrenches, efforts)) = &mut *forcing {
            equations.evaluate_with_forcing(
                x,
                &mut **law,
                wrenches,
                efforts,
                reaction.as_deref_mut(),
                &mut *cancel,
            )
        } else {
            equations.evaluate(x, stage_effort, reaction.as_deref_mut(), &mut *cancel)
        }
    };
    let (solution, mut report) = if n == 0 {
        // All-fixed supported trees still evaluate the actual loads and native
        // dynamics twice. There is no zero-dimensional Jacobian to manufacture.
        evaluate_stage(&guess)?;
        (
            guess,
            NewtonReport {
                iterations: 0,
                residual_norm: 0.0,
                converged: true,
                history: Vec::new(),
                diagnosis: None,
            },
        )
    } else {
        let problem = FallibleFiniteDifference::new(n, config.difference_step, |x: &[f64]| {
            evaluate_stage(x).map(|stage| stage.residual)
        })
        .map_err(difference_error)?;
        let solver_config = NewtonKrylovConfig {
            absolute_tolerance: config.residual_tolerance,
            relative_tolerance: 0.0,
            linear_restart: config.krylov_restart.min(n),
            max_linear_cycles: 8,
            forcing_minimum: 1.0e-3,
            forcing_maximum: 0.1,
            ..NewtonKrylovConfig::default()
        };
        let construction = NewtonKrylovState::new(&problem, guess, solver_config);
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = construction.map_err(ImplicitPrescribedSpatialError::Solver)?;
        let result = solver.run_controlled(&problem, config.maximum_iterations, &mut || {
            problem.has_failed() || (*cancelled.borrow_mut())()
        });
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        match result.progress {
            fs_solver::nonlinear::SolverRunProgress::Complete => {}
            fs_solver::nonlinear::SolverRunProgress::Paused => {
                return Err(PrescribedBaseMotionError::Cancelled.into());
            }
            other => return Err(ImplicitPrescribedSpatialError::SolverProgress(other)),
        }
        if !result.report.converged {
            return Err(ImplicitPrescribedSpatialError::Unresolved(result.report));
        }
        (solver.x, result.report)
    };
    let midpoint = evaluate_stage(&solution)?;
    let residual_norm = fs_solver::norm2(&midpoint.residual);
    report.residual_norm = residual_norm;
    if !residual_norm.is_finite() || residual_norm > config.residual_tolerance {
        report.converged = false;
        return Err(ImplicitPrescribedSpatialError::Unresolved(report));
    }
    let joint_velocity_residual = midpoint
        .residual
        .iter()
        .zip(&scales)
        .map(|(error, scale)| error * scale)
        .collect();
    let midpoint_state = midpoint.stage.state.clone();
    let mut cancel = cancelled.borrow_mut();
    let evaluated = finish_step(
        model,
        state,
        path,
        loads.held,
        duration,
        initial,
        initial_load.stored_energy,
        midpoint.stage,
        &mut *cancel,
    )?;
    let spatial = finish(model, evaluated, loads, bushings, &mut *cancel)?;
    Ok(ImplicitPrescribedSpatialStep {
        midpoint_state,
        spatial,
        nonlinear: report,
        residual_evaluations: evaluations.get(),
        scaled_residual_norm: residual_norm,
        joint_velocity_residual,
    })
}

pub(in crate::articulated::motion::spatial_loads) struct MidpointEquations<'a, B: BushingLaw> {
    pub(in crate::articulated::motion::spatial_loads) model: &'a ArticulatedModel,
    pub(in crate::articulated::motion::spatial_loads) state: &'a PrescribedBaseMotionState,
    pub(in crate::articulated::motion::spatial_loads) path: PrescribedBasePath,
    pub(in crate::articulated::motion::spatial_loads) loads: SpatialMotionLoads<'a>,
    pub(in crate::articulated::motion::spatial_loads) bushings: &'a [B],
    pub(in crate::articulated::motion::spatial_loads) duration: f64,
    pub(in crate::articulated::motion::spatial_loads) midpoint_time: f64,
    pub(in crate::articulated::motion::spatial_loads) scales: &'a [f64],
}

pub(in crate::articulated::motion::spatial_loads) struct MidpointEvaluation<'a> {
    pub(in crate::articulated::motion::spatial_loads) stage: PrescribedMidpoint<'a>,
    pub(in crate::articulated::motion::spatial_loads) residual: Vec<f64>,
}

impl<'a, B: BushingLaw> MidpointEquations<'a, B> {
    fn trial_state(
        &self,
        x: &[f64],
    ) -> Result<PrescribedBaseMotionState, ImplicitPrescribedSpatialError> {
        if x.len() != self.scales.len() {
            return Err(ImplicitPrescribedSpatialError::Invalid(
                "trial rate dimension",
            ));
        }
        let velocity: Vec<_> = x
            .iter()
            .zip(self.scales)
            .map(|(x, scale)| x * scale)
            .collect();
        Ok(PrescribedBaseMotionState {
            time_s: self.midpoint_time,
            base: self.path.midpoint,
            // Numerical off-manifold probes near an exact lock may have a
            // sub-ULP coordinate increment. They remain finite/FK-admitted chart
            // samples, not published updates. finish_step still uses the native
            // strict physical endpoint advance; no rate is clipped or projected.
            joint_position: advance(&self.state.joint_position, &velocity, 0.5 * self.duration),
            joint_velocity: velocity,
        })
    }

    pub(in crate::articulated::motion::spatial_loads) fn evaluate(
        &self,
        x: &[f64],
        effort: Option<(&mut EffortLaw<'_>, &mut [f64])>,
        reaction: Option<&mut ReactionLaw<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MidpointEvaluation<'a>, ImplicitPrescribedSpatialError> {
        let half = 0.5 * self.duration;
        let state = self.trial_state(x)?;
        poll(cancelled)?;
        let kinematics = forward_kinematics(
            self.model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
        )?;
        poll(cancelled)?;
        let load = evaluate_shared(
            self.model,
            &state,
            &kinematics,
            self.loads,
            self.bushings,
            cancelled,
        )?;
        admit_stage(&load)?;
        let generalized_force = if let Some((law, values)) = effort {
            values.copy_from_slice(self.loads.held.generalized_force);
            poll(cancelled)?;
            law(half, &state, values)?;
            poll(cancelled)?;
            &*values
        } else {
            self.loads.held.generalized_force
        };
        self.evaluate_dynamics(
            state,
            kinematics,
            load,
            generalized_force,
            reaction,
            cancelled,
        )
    }

    pub(in crate::articulated::motion::spatial_loads) fn evaluate_with_forcing(
        &self,
        x: &[f64],
        forcing: &mut ForcingLaw<'_>,
        wrenches: &mut [Wrench],
        efforts: &mut [f64],
        reaction: Option<&mut ReactionLaw<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MidpointEvaluation<'a>, ImplicitPrescribedSpatialError> {
        let state = self.trial_state(x)?;
        poll(cancelled)?;
        let kinematics = forward_kinematics(
            self.model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
        )?;
        wrenches.copy_from_slice(self.loads.held.external_body_wrench);
        efforts.copy_from_slice(self.loads.held.generalized_force);
        poll(cancelled)?;
        forcing(0.5 * self.duration, &state, wrenches, efforts)?;
        poll(cancelled)?;
        let loads = SpatialMotionLoads {
            held: crate::articulated::motion::FreeFloatingMotionLoads {
                external_body_wrench: wrenches,
                generalized_force: efforts,
                ..self.loads.held
            },
            ..self.loads
        };
        let evaluated = evaluate_shared(
            self.model,
            &state,
            &kinematics,
            loads,
            self.bushings,
            cancelled,
        )?;
        // The accepted load must outlive the reusable callback buffers. Already
        // owned spatial/bushing loads move directly; held-only loads copy once.
        let load = PrescribedStage {
            body_wrench: Cow::Owned(evaluated.body_wrench.into_owned()),
            external: evaluated.external,
            external_power: evaluated.external_power,
            internal_power: evaluated.internal_power,
            stored_energy: evaluated.stored_energy,
            dissipation_power: evaluated.dissipation_power,
        };
        admit_stage(&load)?;
        self.evaluate_dynamics(state, kinematics, load, efforts, reaction, cancelled)
    }

    fn evaluate_dynamics(
        &self,
        state: PrescribedBaseMotionState,
        kinematics: Kinematics,
        load: PrescribedStage<'a>,
        generalized_force: &[f64],
        reaction: Option<&mut ReactionLaw<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MidpointEvaluation<'a>, ImplicitPrescribedSpatialError> {
        let half = 0.5 * self.duration;
        let mut dynamics = forward_dynamics(
            self.model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
            generalized_force,
            self.loads.held.gravity_world,
            &load.body_wrench,
        )?;
        poll(cancelled)?;
        let constraint_power = if let Some(law) = reaction {
            let power = law(&state, &mut dynamics.generalized_acceleration, cancelled)?;
            poll(cancelled)?;
            validate_slice(
                &dynamics.generalized_acceleration,
                "implicit_prescribed.constrained_acceleration",
            )?;
            power
        } else {
            Dd::ZERO
        };
        let residual: Vec<_> = state
            .joint_velocity
            .iter()
            .zip(&self.state.joint_velocity)
            .zip(&dynamics.generalized_acceleration)
            .zip(self.scales)
            .map(|(((mid, old), acceleration), scale)| (mid - old - half * acceleration) / scale)
            .collect();
        validate_slice(&residual, "implicit_prescribed.equations")?;
        let actuator_power = sum_products(
            generalized_force.iter().copied(),
            state.joint_velocity.iter().copied(),
        );
        validate_slice(
            &[actuator_power.to_f64(), constraint_power.to_f64()],
            "implicit_prescribed.stage_power",
        )?;
        poll(cancelled)?;
        Ok(MidpointEvaluation {
            stage: PrescribedMidpoint {
                state,
                kinematics,
                load,
                generalized_acceleration: dynamics.generalized_acceleration,
                actuator_power,
                constraint_power,
            },
            residual,
        })
    }
}

fn difference_error(
    error: FiniteDifferenceError<ImplicitPrescribedSpatialError>,
) -> ImplicitPrescribedSpatialError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => {
            ImplicitPrescribedSpatialError::Difference(reason)
        }
    }
}
