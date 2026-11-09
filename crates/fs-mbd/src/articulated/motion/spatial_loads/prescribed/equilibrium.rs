//! Static and uniformly rotating equilibrium of the original supported tree.
//!
//! The residual is the native zero-joint-acceleration RNEA holding effort minus the actual
//! actuator effort. World point loads, finite axial links and diagonal/full
//! bushings use the same geometry and constitutive evaluator as time evolution.
//! Newton--Krylov and numerical directional derivatives remain in fs-solver.
//! Finding a root does not establish its stability, uniqueness or reachability.

/// Relative equilibrium under a genuinely uniform supported rotation.
pub mod rotating;

use super::super::bushing::coupled::CoupledSpatialBushing;
use super::super::bushing::{self, BushingLaw, SpatialBushing};
use super::super::implicit::MAX_IMPLICIT_SPATIAL_EVALUATIONS;
use super::{SpatialMotionLoads, evaluate_shared};
use crate::articulated::motion::prescribed::{
    EffortLaw, PrescribedBaseMotionError, PrescribedBaseMotionState, admit_stage, poll,
    support_wrenches,
};
use crate::articulated::motion::{FreeFloatingMotionError, MAX_FREE_FLOATING_STEP_LINKS};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, forward_kinematics, inverse_dynamics,
    validate_external_wrenches, validate_forward_dynamics_header,
};
use fs_ga::{Twist, Wrench};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

/// Original finite-attachment inventory; coefficients are never reinterpreted.
#[derive(Clone, Copy, Debug)]
pub enum EquilibriumBushings<'a> {
    /// Native diagonal translation and principal-log rotation laws; may be empty.
    Diagonal(&'a [SpatialBushing]),
    /// Native passive six-coordinate laws including reciprocal cross couplings.
    Coupled(&'a [CoupledSpatialBushing]),
}

/// Explicit physical scaling and bounded numerical work.
#[derive(Clone, Copy, Debug)]
pub struct EquilibriumConfig<'a> {
    /// Positive position-increment magnitudes, one per original joint: metres for
    /// prismatic coordinates and radians for revolute/helical coordinates.
    pub coordinate_scales: &'a [f64],
    /// Positive effort magnitudes in the same order: N or N m respectively.
    /// The accepted absolute physical residual in row i cannot exceed
    /// `residual_tolerance * effort_scales[i]`.
    pub effort_scales: &'a [f64],
    /// Positive absolute 2-norm limit on the scaled physical force balances.
    /// Relative initial-residual convergence is disabled.
    pub residual_tolerance: f64,
    /// Positive maximum-component difference stencil in scaled coordinates.
    /// A numerical Jacobian parameter, not a displacement accuracy certificate.
    pub difference_step: f64,
    /// Newton attempts in 1..=64, including rejected attempts.
    pub maximum_iterations: usize,
    /// Residual calls in 2..=16384, including stencils and final verification.
    pub maximum_evaluations: usize,
    /// Krylov vectors in 1..=64, reduced to the actual joint count.
    pub krylov_restart: usize,
}

/// A force-balanced native state ready for the existing supported motion owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedEquilibrium {
    /// Same support pose, motion and clock as the seed, solved coordinates and
    /// zero joint rates. The stationary entry points also require a still base;
    /// the rotating entry point explicitly retains its driven base rotation.
    /// This is explicit initialization, not a time-integrated trajectory.
    pub state: PrescribedBaseMotionState,
    /// Native zero-joint-acceleration holding efforts, N or N m by joint kind.
    pub required_effort: Vec<f64>,
    /// Actual held/callback actuator efforts at the independently checked root.
    pub applied_effort: Vec<f64>,
    /// Required minus applied effort in the original physical units.
    pub joint_effort_residual: Vec<f64>,
    /// Actual 2-norm of those residuals divided by their explicit effort scales.
    pub scaled_residual_norm: f64,
    /// Full support wrench ON the tree, about the base origin in base axes.
    pub support_wrench_body: Wrench,
    /// The same support wrench about the world origin, in world axes.
    pub support_wrench_world: Wrench,
    /// Original axial/bushing storage at the root, J. No relaxation work, loss
    /// or heat is inferred from displacement between the seed and this root.
    pub attachment_energy_j: f64,
    /// Original solver report with an independently repeated final residual.
    pub nonlinear: NewtonReport,
    /// All actual residual calls, including numerical differences and recheck.
    pub residual_evaluations: usize,
}

/// Refusal never modifies the seed or publishes an unbalanced candidate.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumError {
    /// Invalid physical scales, seed or bounded controls.
    Invalid(&'static str),
    /// Native geometry, load, effort, material-chart or cancellation refusal.
    Motion(PrescribedBaseMotionError),
    /// Residual work was exhausted, possibly during final verification.
    EvaluationBudget {
        /// Residual evaluations already consumed.
        used: usize,
        /// Caller-admitted maximum.
        maximum: usize,
    },
    /// Native finite-difference stencil refusal.
    Difference(&'static str),
    /// Native nonlinear construction refusal.
    Solver(NewtonError),
    /// Solver callback or dimension failure; no candidate is accepted.
    SolverProgress(SolverRunProgress),
    /// Required force balance was not achieved. Singular loaded mechanisms
    /// and branches without a nearby equilibrium may return this refusal.
    Unresolved(NewtonReport),
}
impl core::fmt::Display for EquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "supported articulated equilibrium: {self:?}")
    }
}
impl std::error::Error for EquilibriumError {}
impl From<PrescribedBaseMotionError> for EquilibriumError {
    fn from(value: PrescribedBaseMotionError) -> Self {
        Self::Motion(value)
    }
}
impl From<ArticulatedError> for EquilibriumError {
    fn from(value: ArticulatedError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<FreeFloatingMotionError> for EquilibriumError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value.into())
    }
}

/// Solve stationary force balance without fictitious damping or time evolution.
///
/// The seed must have zero base/joint velocity and zero base acceleration. Its
/// coordinates select a local branch and its clock is retained exactly. Native
/// joint limits apply to every trial; hard-stop/contact reactions are not added.
/// Prescribed motion, Coulomb sticking and closed kinematic loops require their
/// dedicated owners. An already balanced mechanism can have a neutral mode;
/// neither stability nor uniqueness is asserted by successful force balance.
///
/// Each residual has O(links+loads) work and O(links) temporary storage. Krylov
/// storage is O(restart*joints), with no dense mass or tangent matrix. Bounded
/// native passes and source evaluations poll cancellation. A rejected source
/// or difference stencil stops the solve rather than clipping its geometry.
pub fn equilibrate(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    config: EquilibriumConfig<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedEquilibrium, EquilibriumError> {
    dispatch(
        model,
        seed,
        loads,
        bushings,
        config,
        None,
        FrameMode::Stationary,
        cancelled,
    )
}

/// Use the same smooth effort callback as supported time evolution.
///
/// The elapsed-time argument is always zero and the absolute clock is the
/// seed's clock. Every call receives the original held efforts as its starting
/// buffer. Calls include rejected trials, numerical stencils and final recheck;
/// the callback must be repeatable and must not commit constitutive history.
/// Constitutive storage behind a joint-effort law remains caller-owned.
/// Long callbacks own their internal cancellation checks; callback side effects
/// cannot be rolled back even though no native candidate is published on error.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_with_effort(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    config: EquilibriumConfig<'_>,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedEquilibrium, EquilibriumError> {
    dispatch(
        model,
        seed,
        loads,
        bushings,
        config,
        Some(&mut effort),
        FrameMode::Stationary,
        cancelled,
    )
}

#[derive(Clone, Copy)]
enum FrameMode {
    Stationary,
    UniformRotation,
}

#[allow(clippy::too_many_arguments)]
fn dispatch(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    config: EquilibriumConfig<'_>,
    effort: Option<&mut EffortLaw<'_>>,
    frame: FrameMode,
    cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedEquilibrium, EquilibriumError> {
    match bushings {
        EquilibriumBushings::Diagonal(values) => {
            solve(model, seed, loads, values, config, effort, frame, cancelled)
        }
        EquilibriumBushings::Coupled(values) => {
            solve(model, seed, loads, values, config, effort, frame, cancelled)
        }
    }
}

struct Evaluation {
    state: PrescribedBaseMotionState,
    required: Vec<f64>,
    applied: Vec<f64>,
    physical_residual: Vec<f64>,
    scaled_residual: Vec<f64>,
    root_wrench: Wrench,
    energy: f64,
}

#[allow(clippy::too_many_arguments)]
fn solve<B: BushingLaw>(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    config: EquilibriumConfig<'_>,
    effort: Option<&mut EffortLaw<'_>>,
    frame: FrameMode,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedEquilibrium, EquilibriumError> {
    poll(&mut cancelled)?;
    let n = model.dof_count();
    if !seed.time_s.is_finite()
        || match frame {
            FrameMode::Stationary => seed.base.twist_body != Twist::zero(),
            FrameMode::UniformRotation => seed.base.twist_body.linear != Twist::zero().linear,
        }
        || seed.base.acceleration_body != Twist::zero()
        || seed.joint_velocity.iter().any(|v| *v != 0.0)
    {
        return Err(EquilibriumError::Invalid(
            "zero joint rates, admitted support motion and finite clock required",
        ));
    }
    if config.coordinate_scales.len() != n
        || config.effort_scales.len() != n
        || config
            .coordinate_scales
            .iter()
            .chain(config.effort_scales)
            .any(|v| !v.is_finite() || *v <= 0.0)
        || [config.residual_tolerance, config.difference_step]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(EquilibriumError::Invalid(
            "complete positive physical scales and tolerances required",
        ));
    }
    if config.maximum_iterations == 0
        || config.maximum_iterations > 64
        || config.maximum_evaluations < 2
        || config.maximum_evaluations > MAX_IMPLICIT_SPATIAL_EVALUATIONS
        || config.krylov_restart == 0
        || config.krylov_restart > 64
    {
        return Err(EquilibriumError::Invalid(
            "bounded nonlinear controls required",
        ));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(PrescribedBaseMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    validate_forward_dynamics_header(
        model,
        seed.base.twist_body,
        &seed.joint_position,
        &seed.joint_velocity,
        loads.held.generalized_force,
        loads.held.gravity_world,
        loads.held.external_body_wrench,
    )?;
    validate_external_wrenches(loads.held.external_body_wrench)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    // Solve increments about the original seed. A zero increment preserves the
    // supplied coordinates exactly, even with nonbinary physical scales.
    let guess = vec![0.0; n];
    let evaluations = Cell::new(0usize);
    let cancel = RefCell::new(cancelled);
    let effort = RefCell::new(effort);
    let evaluate = |x: &[f64]| -> Result<Evaluation, EquilibriumError> {
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
        if x.len() != n {
            return Err(EquilibriumError::Invalid("trial coordinate dimension"));
        }
        let state = PrescribedBaseMotionState {
            time_s: seed.time_s,
            base: seed.base,
            joint_position: x
                .iter()
                .zip(config.coordinate_scales)
                .zip(&seed.joint_position)
                .map(|((&increment, &scale), &initial)| increment.mul_add(scale, initial))
                .collect(),
            joint_velocity: vec![0.0; n],
        };
        let kinematics = forward_kinematics(
            model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
        )?;
        poll(&mut *cancel)?;
        let stage = evaluate_shared(model, &state, &kinematics, loads, bushings, &mut *cancel)?;
        admit_stage(&stage)?;
        let mut applied = loads.held.generalized_force.to_vec();
        if let Some(law) = effort.borrow_mut().as_deref_mut() {
            law(0.0, &state, &mut applied)?;
            poll(&mut *cancel)?;
        }
        model.validate_applied_effort(&applied)?;
        let inverse = inverse_dynamics(
            model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
            &state.joint_velocity,
            loads.held.gravity_world,
            &stage.body_wrench,
        )?;
        poll(&mut *cancel)?;
        let physical_residual: Vec<_> = inverse
            .generalized_force
            .iter()
            .zip(&applied)
            .map(|(a, b)| a - b)
            .collect();
        let scaled_residual: Vec<_> = physical_residual
            .iter()
            .zip(config.effort_scales)
            .map(|(a, b)| a / b)
            .collect();
        if physical_residual
            .iter()
            .chain(&scaled_residual)
            .any(|v| !v.is_finite())
        {
            return Err(EquilibriumError::Invalid("unrepresentable force balance"));
        }
        Ok(Evaluation {
            state,
            required: inverse.generalized_force,
            applied,
            physical_residual,
            scaled_residual,
            root_wrench: inverse.body_wrench[0],
            energy: stage.stored_energy,
        })
    };
    let (solution, mut report) = if n == 0 {
        evaluate(&guess)?;
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
            evaluate(x).map(|value| value.scaled_residual)
        })
        .map_err(difference_error)?;
        let construction = NewtonKrylovState::new(
            &problem,
            guess,
            NewtonKrylovConfig {
                absolute_tolerance: config.residual_tolerance,
                relative_tolerance: 0.0,
                linear_restart: config.krylov_restart.min(n),
                max_linear_cycles: 8,
                forcing_minimum: 1.0e-3,
                forcing_maximum: 0.1,
                ..NewtonKrylovConfig::default()
            },
        );
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = construction.map_err(EquilibriumError::Solver)?;
        let result = solver.run_controlled(&problem, config.maximum_iterations, &mut || {
            problem.has_failed() || (*cancel.borrow_mut())()
        });
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        match result.progress {
            SolverRunProgress::Complete => {}
            SolverRunProgress::Paused => return Err(PrescribedBaseMotionError::Cancelled.into()),
            other => return Err(EquilibriumError::SolverProgress(other)),
        }
        if !result.report.converged {
            return Err(EquilibriumError::Unresolved(result.report));
        }
        (solver.x, result.report)
    };
    let final_state = evaluate(&solution)?;
    let residual_norm = fs_solver::norm2(&final_state.scaled_residual);
    report.residual_norm = residual_norm;
    if !residual_norm.is_finite() || residual_norm > config.residual_tolerance {
        report.converged = false;
        return Err(EquilibriumError::Unresolved(report));
    }
    let (support_wrench_body, support_wrench_world, _) =
        support_wrenches(model, &final_state.state, final_state.root_wrench)?;
    poll(&mut *cancel.borrow_mut())?;
    Ok(ArticulatedEquilibrium {
        state: final_state.state,
        required_effort: final_state.required,
        applied_effort: final_state.applied,
        joint_effort_residual: final_state.physical_residual,
        scaled_residual_norm: residual_norm,
        support_wrench_body,
        support_wrench_world,
        attachment_energy_j: final_state.energy,
        nonlinear: report,
        residual_evaluations: evaluations.get(),
    })
}

fn difference_error(error: FiniteDifferenceError<EquilibriumError>) -> EquilibriumError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => EquilibriumError::Difference(reason),
    }
}
