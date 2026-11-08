//! Coupled floating-pose initialization on the original reservoir load owner.

use super::pose::HydrostaticPoseAxes;
use super::{
    DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticReservoir, finite, poll,
};
use crate::principal::PrincipalFrame;
use fs_mbd::{RigidBodyState, Vec3};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonStallDiagnosis, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

/// Hard bound on original reservoir load queries, including difference stencils.
pub const MAX_FLOATING_EVALUATIONS: usize = 8192;
/// Hard bound on native Newton attempts, including rejected globalization steps.
pub const MAX_FLOATING_NEWTON_ITERATIONS: usize = 128;

/// Physical residual allowances, a local pose domain and explicit numerical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FloatingEquilibriumControls {
    /// Direction projected perpendicular to gravity to define world roll.
    /// Pitch is `up cross roll`; a zero or gravity-parallel direction refuses.
    pub roll_direction_world: Vec3,
    /// Positive absolute norm allowance for actual buoyancy plus body weight [N].
    pub maximum_force_defect_n: f64,
    /// Positive absolute norm allowance for the actual buoyancy torque [N m].
    pub maximum_torque_defect_nm: f64,
    /// Positive translation scale [m] for the dimensionless nonlinear coordinates.
    /// Near the reference, the heave perturbation is `length_scale_m * x[0]`.
    pub length_scale_m: f64,
    /// Positive maximum-component difference step in the nonlinear coordinates.
    /// Near the reference this means `length_scale_m * difference_step` metres
    /// in heave, or `difference_step` radians in tilt. This is a numerical
    /// Jacobian scale, not an analytic tangent or a solution-error bound.
    pub difference_step: f64,
    /// Positive maximum absolute COM heave displacement from the reference [m].
    /// All load queries, including stencils and rejected trials, remain bounded.
    pub maximum_heave_offset_m: f64,
    /// Maximum norm of the horizontal world rotation vector, in `(0, pi)` [rad].
    /// This is a local search about the reference; it does not search all poses.
    pub maximum_tilt_radians: f64,
    /// Native Newton attempts, in `1..=MAX_FLOATING_NEWTON_ITERATIONS`.
    pub maximum_newton_iterations: usize,
    /// Actual reservoir load calls, in `1..=MAX_FLOATING_EVALUATIONS`.
    /// Repeated identical coordinates reuse a two-entry immutable load cache.
    pub maximum_evaluations: usize,
}

/// A zero-momentum floating pose satisfying both physical residual allowances.
/// Force/torque balance alone does not imply stability, uniqueness or damping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FloatingEquilibrium {
    /// Candidate on the original principal frame, with both momenta zero.
    pub state: RigidBodyState,
    /// Actual original load, including any solved finite-tank interfaces.
    pub load: HydrostaticLoad,
    /// COM displacement from the reference, positive opposite gravity [m].
    pub heave_offset_m: f64,
    /// Horizontal world rotation vector left-composed onto the reference [rad].
    pub tilt_rotation_vector_world_radians: Vec3,
    /// Actual buoyancy plus weight, without a compensating external force [N].
    pub net_force_world_n: Vec3,
    /// Actual torque about the COM in the world frame [N m]. No holding torque
    /// is applied. The original principal-frame torque is in `load.wrench()`.
    pub torque_world_nm: Vec3,
    /// Actual reservoir queries, including those used by numerical differences.
    pub evaluations: usize,
    /// Completed native Newton attempts; zero if the reference already balances.
    pub iterations: usize,
}

/// A failed initialization never returns an approximate candidate or changes history.
#[derive(Clone, Debug, PartialEq)]
pub enum FloatingEquilibriumError {
    /// Original geometry, reservoir, native pose, cancellation or input refusal.
    Physics(HydrostaticError),
    /// The native finite-difference adapter refused an unrepresentable stencil.
    Difference(&'static str),
    /// Native nonlinear construction failure.
    Solver(NewtonError),
    /// Native callback interruption other than an ordinary cancellation request.
    Interrupted(SolverRunProgress),
    /// The explicit original-load budget ended, including inside a stencil.
    EvaluationBudget {
        evaluations: usize,
        maximum_evaluations: usize,
    },
    /// Native iteration/globalization stopped without a resolved local root.
    /// This is not a proof that no equilibrium exists in another pose domain.
    NotConverged {
        diagnosis: Option<NewtonStallDiagnosis>,
        scaled_residual_norm: f64,
        evaluations: usize,
        iterations: usize,
    },
    /// Full physical force or torque failed after the three projected residuals
    /// converged. In particular, an inconsistent reservoir cannot hide a lateral
    /// force or vertical torque in the neutral coordinates.
    ResidualDefect {
        force_defect_n: f64,
        torque_defect_nm: f64,
        evaluations: usize,
    },
}
impl core::fmt::Display for FloatingEquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hydrostatic floating equilibrium: {self:?}")
    }
}
impl std::error::Error for FloatingEquilibriumError {}
impl From<HydrostaticError> for FloatingEquilibriumError {
    fn from(error: HydrostaticError) -> Self {
        Self::Physics(error)
    }
}

fn numerical(error: FiniteDifferenceError<FloatingEquilibriumError>) -> FloatingEquilibriumError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => FloatingEquilibriumError::Difference(reason),
    }
}

#[derive(Clone, Copy)]
struct Trial {
    coordinates: [f64; 3],
    state: RigidBodyState,
    load: HydrostaticLoad,
    heave_offset_m: f64,
    tilt_world: Vec3,
    force_world: Vec3,
    torque_world: Vec3,
}

impl FloatingEquilibriumControls {
    fn validate(self) -> Result<(), HydrostaticError> {
        if [
            self.maximum_force_defect_n,
            self.maximum_torque_defect_nm,
            self.length_scale_m,
            self.difference_step,
            self.maximum_heave_offset_m,
            self.maximum_tilt_radians,
            self.length_scale_m / self.maximum_heave_offset_m,
        ]
        .into_iter()
        .any(|value| !value.is_finite() || value <= 0.0)
            || self.maximum_tilt_radians >= core::f64::consts::PI
            || !(1..=MAX_FLOATING_NEWTON_ITERATIONS).contains(&self.maximum_newton_iterations)
            || !(1..=MAX_FLOATING_EVALUATIONS).contains(&self.maximum_evaluations)
        {
            return Err(HydrostaticError::Invalid(
                "floating equilibrium requires positive finite physical/difference scales, a local pose domain and bounded work",
            ));
        }
        Ok(())
    }

    // Saturation keeps line-search and directional-difference queries inside
    // the same physical domain. The Jacobian still belongs to the original
    // fallible finite-difference owner; no invalid physical trial is retried.
    fn pose_coordinates(self, x: [f64; 3]) -> (f64, [f64; 2]) {
        let heave = self.maximum_heave_offset_m
            * (x[0] * (self.length_scale_m / self.maximum_heave_offset_m)).tanh();
        let scale = x[1].abs().max(x[2].abs());
        let tilt = if scale == 0.0 {
            [0.0; 2]
        } else {
            let direction = [x[1] / scale, x[2] / scale];
            let norm = direction[0].hypot(direction[1]);
            let angle =
                self.maximum_tilt_radians * ((scale / self.maximum_tilt_radians) * norm).tanh();
            direction.map(|component| angle * (component / norm))
        };
        (heave, tilt)
    }
}

/// Initialize a locally torque-free floating pose through the native nonlinear solver.
///
/// The three unknowns are COM heave and two horizontal-world tilt coordinates.
/// Lateral COM position stays fixed. The tilt is one exponential relative to
/// the reference orientation, with no additional rotation about gravity. This
/// chooses a yaw-neutral local chart, not an external yaw constraint or torque.
/// Each pose has zero momentum and is evaluated by the SAME `load_frame` used
/// for dynamics. Material sources remain frozen; finite tanks re-solve their
/// conserved inventories at every new pose.
///
/// Heave is `H*tanh(x[0]*L/H)`. The horizontal rotation has the direction of
/// `(x[1],x[2])` and magnitude `Theta*tanh(norm(x[1..])/Theta)`. Thus even
/// rejected Newton trials and numerical-difference stencils stay within the
/// caller's heave/tilt envelope. Tank contact, overflow, unresolved inner roots
/// and any other original reservoir refusal still terminate the entire solve;
/// they are never replaced with a fictitious residual for backtracking.
///
/// The original `FallibleFiniteDifference` and `NewtonKrylovState` own all
/// Jacobian actions, Krylov solves and Armijo globalization. The residual is
/// `[F_up/F_tolerance, torque_roll/T_tolerance, torque_pitch/T_tolerance]`.
/// Native absolute tolerance is 0.5, relative tolerance zero; its three-column
/// FGMRES uses at most eight cycles per attempt and a `1e-4` forcing floor to
/// avoid demanding analytic-Jacobian precision from numerical differences.
/// Final admission independently checks FULL world force and torque norms.
/// These are residual checks, not a pose-error enclosure or stability proof.
///
/// This initializes a new state, not a dissipative step or a conservative live
/// replacement. A balanced unstable attitude or a neutral family is admissible.
/// No-root, ill-conditioned, unresolved or out-of-budget solves return an error
/// without publishing an approximate root. Cancellation is latched through both
/// native solver callbacks and the original reservoir's per-cell checkpoints.
/// Native solver scratch is bounded by dimension three and the attempt cap;
/// two exact-coordinate load entries avoid repeated base-pose clipping.
pub fn solve_floating_equilibrium<R: HydrostaticReservoir + ?Sized, C: FnMut() -> bool>(
    reservoir: &R,
    mesh: &DisplacementMesh,
    body: &PrincipalFrame,
    reference: RigidBodyState,
    controls: FloatingEquilibriumControls,
    mut cancelled: C,
) -> Result<FloatingEquilibrium, FloatingEquilibriumError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let gravity = reservoir.gravity();
    let axes = HydrostaticPoseAxes::new(gravity, controls.roll_direction_world)?;
    let weight = gravity.acceleration_world().scale(body.properties().mass());
    if !weight.is_finite() {
        return Err(HydrostaticError::Invalid("body weight is unrepresentable").into());
    }

    let cancel = RefCell::new(cancelled);
    let cancellation_seen = Cell::new(false);
    let stop = || {
        if cancellation_seen.get() {
            true
        } else {
            let requested = (*cancel.borrow_mut())();
            cancellation_seen.set(requested);
            requested
        }
    };
    let evaluations = Cell::new(0usize);
    let cache: RefCell<[Option<Trial>; 2]> = RefCell::new([None; 2]);
    let evaluate = |coordinates: &[f64]| -> Result<Trial, FloatingEquilibriumError> {
        poll(&mut || stop())?;
        let x: [f64; 3] = coordinates.try_into().map_err(|_| {
            HydrostaticError::Invalid("floating pose needs three nonlinear coordinates")
        })?;
        if x.iter().any(|value| !value.is_finite()) {
            return Err(HydrostaticError::Invalid("nonfinite floating pose coordinates").into());
        }
        {
            let mut entries = cache.borrow_mut();
            for i in 0..2 {
                if let Some(trial) = entries[i]
                    && trial
                        .coordinates
                        .iter()
                        .zip(x)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                {
                    entries.swap(0, i);
                    return Ok(trial);
                }
            }
        }
        if evaluations.get() == controls.maximum_evaluations {
            return Err(FloatingEquilibriumError::EvaluationBudget {
                evaluations: evaluations.get(),
                maximum_evaluations: controls.maximum_evaluations,
            });
        }
        let (heave_offset_m, tilt) = controls.pose_coordinates(x);
        let state = axes.displaced_state(reference, heave_offset_m, tilt)?;
        evaluations.set(evaluations.get() + 1);
        let load = reservoir.load_frame(mesh, body, state, || stop())?;
        let force_world = load.wrench().force_world.add(weight);
        let torque_world = state
            .pose()
            .orientation()
            .rotate_body_to_world(load.wrench().torque_body);
        if !force_world.is_finite() || !torque_world.is_finite() {
            return Err(HydrostaticError::Invalid("floating residual is unrepresentable").into());
        }
        poll(&mut || stop())?;
        let trial = Trial {
            coordinates: x,
            state,
            load,
            heave_offset_m,
            tilt_world: axes.roll.scale(tilt[0]).add(axes.pitch.scale(tilt[1])),
            force_world,
            torque_world,
        };
        let mut entries = cache.borrow_mut();
        entries[1] = entries[0];
        entries[0] = Some(trial);
        Ok(trial)
    };
    let problem = FallibleFiniteDifference::new(3, controls.difference_step, |x: &[f64]| {
        let trial = evaluate(x)?;
        Ok::<_, FloatingEquilibriumError>(vec![
            finite(trial.force_world.dot(axes.up) / controls.maximum_force_defect_n)?,
            finite(trial.torque_world.dot(axes.roll) / controls.maximum_torque_defect_nm)?,
            finite(trial.torque_world.dot(axes.pitch) / controls.maximum_torque_defect_nm)?,
        ])
    })
    .map_err(numerical)?;
    let config = NewtonKrylovConfig {
        absolute_tolerance: 0.5,
        relative_tolerance: 0.0,
        linear_restart: 3,
        max_linear_cycles: 8,
        forcing_minimum: 1e-4,
        forcing_maximum: 0.25,
        ..NewtonKrylovConfig::default()
    };
    let construction = NewtonKrylovState::new(&problem, vec![0.0; 3], config);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction.map_err(FloatingEquilibriumError::Solver)?;
    let outcome = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        stop() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    poll(&mut || stop())?;
    if outcome.progress != SolverRunProgress::Complete {
        return Err(FloatingEquilibriumError::Interrupted(outcome.progress));
    }
    if !outcome.report.converged {
        return Err(FloatingEquilibriumError::NotConverged {
            diagnosis: outcome.report.diagnosis,
            scaled_residual_norm: outcome.report.residual_norm,
            evaluations: evaluations.get(),
            iterations: outcome.report.iterations,
        });
    }
    let trial = evaluate(&solver.x)?;
    let force_defect_n = finite(
        trial
            .force_world
            .x
            .hypot(trial.force_world.y)
            .hypot(trial.force_world.z),
    )?;
    let torque_defect_nm = finite(
        trial
            .torque_world
            .x
            .hypot(trial.torque_world.y)
            .hypot(trial.torque_world.z),
    )?;
    if force_defect_n > controls.maximum_force_defect_n
        || torque_defect_nm > controls.maximum_torque_defect_nm
    {
        return Err(FloatingEquilibriumError::ResidualDefect {
            force_defect_n,
            torque_defect_nm,
            evaluations: evaluations.get(),
        });
    }
    poll(&mut || stop())?;
    Ok(FloatingEquilibrium {
        state: trial.state,
        load: trial.load,
        heave_offset_m: trial.heave_offset_m,
        tilt_rotation_vector_world_radians: trial.tilt_world,
        net_force_world_n: trial.force_world,
        torque_world_nm: trial.torque_world,
        evaluations: evaluations.get(),
        iterations: outcome.report.iterations,
    })
}
