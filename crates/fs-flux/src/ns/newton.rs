//! Steady and fully implicit BDF1/BDF2/SDIRK2 flow with a consistent
//! constitutive Jacobian and line search.
//!
//! Every accepted endpoint solves the original nonlinear momentum equations.
//! Damping changes the iteration, never the viscosity law or physical timestep.

use std::cell::RefCell;

use super::boundary::{AdmittedBoundary, StokesBoundary, StokesBoundaryError};
use super::linearization::{FluxLinearizationError, FluxViscosityTangent, rheology_residual, sample};
use super::{FluxParams, FluxSolution, FluxSystem};

mod sdirk2;
pub use sdirk2::{FluxSdirk2Error, FluxSdirk2Stage, FluxSdirk2Step};

/// Explicit nonlinear work and acceptance controls.
#[derive(Clone, Copy, Debug)]
pub struct FluxNewtonControls {
    /// Maximum Jacobian solves, excluding an already converged initial state.
    pub maximum_iterations: usize,
    /// Maximum halvings per line search, from zero through 52.
    pub maximum_backtracks: usize,
    /// Both original and row-equilibrated endpoint residuals must meet this.
    pub residual_tolerance: f64,
}

impl Default for FluxNewtonControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 40,
            maximum_backtracks: 16,
            residual_tolerance: 1e-9,
        }
    }
}

/// Accepted endpoint and actual nonlinear work, with the previous state intact.
pub struct FluxNewtonStep {
    /// Endpoint field; residual fields measure the reassembled nonlinear
    /// equations, and `iters` counts the last linear Newton-correction solve.
    pub solution: FluxSolution,
    /// Completed Newton-correction solves.
    pub iterations: usize,
    /// Total rejected line-search candidates.
    pub backtracks: usize,
    /// Last relative velocity update; zero for an initially converged state.
    pub relative_update: f64,
    /// Larger of the endpoint original and row-equilibrated relative defects.
    pub nonlinear_residual: f64,
}

/// No endpoint is returned after invalid data, failed linear work or stagnation.
#[derive(Clone, Debug, PartialEq)]
pub enum FluxNewtonError {
    /// Invalid physical/numerical controls or unrepresentable iteration state.
    Invalid(&'static str),
    /// Existing boundary topology, value or net-flux admission refused.
    Boundary(StokesBoundaryError),
    /// BDF2 history, step ratio or mass-history arithmetic refused.
    History(StokesBoundaryError),
    /// The exact constitutive linearization refused its state or input.
    Linearization(FluxLinearizationError),
    /// Singular or unconverged Newton correction.
    LinearSolve {
        iteration: usize,
        error: StokesBoundaryError,
    },
    /// No finite, admissible candidate reduced the actual nonlinear residual.
    LineSearch {
        iteration: usize,
        attempted_steps: usize,
        residual: f64,
    },
    /// The explicit Newton budget ended without the required endpoint residual.
    NonConvergence { iterations: usize, residual: f64 },
}

impl core::fmt::Display for FluxNewtonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "Newton flow: {reason}"),
            Self::Boundary(error) => write!(f, "Newton flow boundary: {error}"),
            Self::History(error) => write!(f, "Newton flow history: {error}"),
            Self::Linearization(error) => write!(f, "Newton flow: {error}"),
            Self::LinearSolve { iteration, error } => {
                write!(f, "Newton flow correction {iteration}: {error}")
            }
            Self::LineSearch {
                iteration,
                attempted_steps,
                residual,
            } => write!(
                f,
                "Newton flow: iteration {iteration} found no descent in {attempted_steps} trials (residual {residual})"
            ),
            Self::NonConvergence {
                iterations,
                residual,
            } => write!(
                f,
                "Newton flow: exhausted {iterations} iterations (residual {residual})"
            ),
        }
    }
}
impl std::error::Error for FluxNewtonError {}

struct Trial {
    residual: Vec<f64>,
    original: f64,
    scaled: f64,
}

fn weighted_norm(residual: &[f64], scales: &[f64]) -> f64 {
    residual
        .iter()
        .zip(scales)
        .fold(0.0, |value, (r, s)| value.hypot(r / s))
}

impl FluxSystem<'_> {
    /// Solve steady generalized-Newtonian flow with bounded damped Newton
    /// iteration and the same constitutive Jacobian and residual acceptance as
    /// [`Self::bdf1_rheology_newton_step`]. `initial_velocity` is an immutable
    /// velocity-moment guess, not a previous physical state or time history.
    /// No mass term or artificial timestep is added to the steady equations.
    ///
    /// Set `convection` to true for steady Navier–Stokes, or false for creeping
    /// flow. Both cases differentiate the full symmetric viscous stress,
    /// viscosity-dependent interface weights and boundary loads. The law and
    /// all fixed external fields must be pure throughout the solve.
    ///
    /// # Errors
    /// Refuses invalid data, incompatible boundary flux, failed Jacobian/linear
    /// solves, exhausted line searches and unmet nonlinear convergence.
    #[allow(clippy::too_many_arguments)]
    pub fn steady_rheology_newton(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        convection: bool,
        controls: FluxNewtonControls,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        self.rheology_newton_solve_checked(
            params,
            law,
            force,
            boundary,
            initial_velocity,
            convection,
            0.0,
            None,
            controls,
            &|| true,
        )
    }

    /// Advance one fully implicit generalized-Newtonian BDF1 step with bounded
    /// damped Newton iteration. Velocity moments are m/s, time seconds, forces
    /// acceleration, pressure kinematic, and density one fixed reference value.
    /// The law supplies viscosity and its analytic strain-rate derivative.
    /// All callbacks must be pure; external material fields and boundary values
    /// belong to the endpoint and stay fixed throughout the solve.
    ///
    /// The existing exact constitutive/convection Jacobian supplies corrections.
    /// Backtracking tests the actual nonlinear residual using a row-weighted
    /// Armijo criterion (coefficient 1e-4); row weights are held fixed within
    /// each search. Trial assemblies do not build Jacobians or factor matrices.
    /// An invalid trial constitutive state is rejected and damping is reduced.
    /// Neither viscosity nor the timestep is clipped or silently changed.
    ///
    /// Acceptance requires both endpoint residual norms, independently assembled
    /// at the returned velocity, to meet the nonlinear tolerance. A small update
    /// alone never establishes convergence. The existing boundary owner checks
    /// global flux compatibility before pressure pinning can omit a continuity
    /// row. The old velocity is immutable; all rejected work remains private.
    /// Upwind flux switches use the Jacobian's semismooth active branch. Newton
    /// convergence is conditional; the iteration and backtrack budgets are hard
    /// limits, not a promise that arbitrary rheology or forcing will converge.
    ///
    /// # Errors
    /// Refuses invalid data, incompatible boundary flux, failed Jacobian/linear
    /// solves, exhausted line searches and unmet nonlinear convergence.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn bdf1_rheology_newton_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        self.bdf1_rheology_newton_step_checked(
            params,
            law,
            force,
            boundary,
            previous_velocity,
            timestep_s,
            controls,
            &|| true,
        )
    }

    /// Advance equal-step, fully implicit generalized-Newtonian BDF2 with
    /// consistent damped Newton iteration. `older_velocity` belongs one
    /// timestep before `previous_velocity`; `None` is exactly BDF1 startup.
    /// All constitutive, force and boundary callbacks belong to the new time.
    ///
    /// The derivative is `(3 u_new - 4 u_previous + u_older) / (2 dt)`.
    /// Only the mass history is combined; the Newton initial iterate remains
    /// the actual previous velocity, and viscosity and convection are solved
    /// together at the returned endpoint. No law is frozen or extrapolated.
    ///
    /// # Errors
    /// Uses the history, nonlinear and source admission of
    /// [`Self::bdf2_rheology_newton_variable_step`].
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_rheology_newton_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        self.bdf2_rheology_newton_variable_step(
            params,
            law,
            force,
            boundary,
            previous_velocity,
            older_velocity.map(|older| (older, timestep_s)),
            timestep_s,
            controls,
        )
    }

    /// Fully implicit BDF2 with unequal intervals and the exact constitutive
    /// Jacobian. `older_history` pairs `u_(n-1)` with its interval to `u_n`;
    /// `timestep_s` is the interval from `u_n` to the requested endpoint.
    /// `None` calls the original BDF1 path, including its acceptance criteria.
    ///
    /// For `r = dt / previous_dt`, the mass derivative has coefficients
    /// `[(1+2r)/(1+r), -(1+r), r*r/(1+r)] / dt`. The existing BDF2 history
    /// owner forms the normalized mass load once, outside Newton iteration.
    /// It admits `0 < r <= 2`, below the BDF2 zero-stability ratio limit;
    /// neither steps nor history values are silently clipped. The starting
    /// iterate is the actual `u_n`, not the mass-history extrapolation.
    ///
    /// Smooth resolved flows retain second-order temporal accuracy when the
    /// histories, endpoint data and nonlinear tolerances resolve that error.
    /// BDF2 is not a guarantee of monotone single-state kinetic energy, and
    /// endpoint viscous heating remains an instantaneous physical rate.
    /// Both original and equilibrated nonlinear endpoint residuals must pass.
    ///
    /// # Errors
    /// Refuses invalid or unrepresentable histories/intervals, a step ratio
    /// outside `(0, 2]`, and the same boundary, constitutive and nonlinear
    /// failures as [`Self::bdf1_rheology_newton_step`]. Histories are immutable.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_rheology_newton_variable_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        self.bdf2_rheology_newton_variable_step_checked(
            params,
            law,
            force,
            boundary,
            previous_velocity,
            older_history,
            timestep_s,
            controls,
            &|| true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf2_rheology_newton_variable_step_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        if older_history.is_none() {
            return self.bdf1_rheology_newton_step_checked(
                params,
                law,
                force,
                boundary,
                previous_velocity,
                timestep_s,
                controls,
                source_is_valid,
            );
        }
        let (history, effective_timestep) = self
            .mixed_bdf2_history(previous_velocity, older_history, timestep_s)
            .map_err(FluxNewtonError::History)?;
        self.rheology_newton_solve_checked(
            params,
            law,
            force,
            boundary,
            previous_velocity,
            true,
            effective_timestep.recip(),
            Some(history.as_deref().unwrap_or(previous_velocity)),
            controls,
            source_is_valid,
        )
    }

    // A material source keeps its original error outside the scalar callback.
    // Finish only the assembly already in progress, then stop before solving,
    // accepting a trial, or publishing an endpoint using fallback coefficients.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(crate) fn bdf1_rheology_newton_step_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        if !timestep_s.is_finite() || timestep_s <= 0.0 || !timestep_s.recip().is_finite() {
            return Err(FluxNewtonError::Invalid(
                "finite positive timestep with representable reciprocal required",
            ));
        }
        self.rheology_newton_solve_checked(
            params,
            law,
            force,
            boundary,
            previous_velocity,
            true,
            timestep_s.recip(),
            Some(previous_velocity),
            controls,
            source_is_valid,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(crate) fn rheology_newton_solve_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        convection: bool,
        inverse_timestep: f64,
        previous_velocity: Option<&[f64]>,
        controls: FluxNewtonControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        self.rheology_newton_solve_admitted_checked(
            params,
            law,
            force,
            boundary,
            initial_velocity,
            convection,
            inverse_timestep,
            previous_velocity,
            controls,
            source_is_valid,
            None,
        )
    }

    // One bounded Newton/line-search owner for the original velocity-boundary
    // problem and admitted mixed Stokes flow. Mixed boundaries are sampled
    // once by their caller and reused by every Jacobian and trial residual.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn rheology_newton_solve_admitted_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        initial_velocity: &[f64],
        convection: bool,
        inverse_timestep: f64,
        previous_velocity: Option<&[f64]>,
        controls: FluxNewtonControls,
        source_is_valid: &dyn Fn() -> bool,
        admitted: Option<&AdmittedBoundary>,
    ) -> Result<FluxNewtonStep, FluxNewtonError> {
        if self.mesh.tris.is_empty()
            || initial_velocity.len() != self.n_u
            || initial_velocity.iter().any(|v| !v.is_finite())
            || !params.tol.is_finite()
            || params.tol <= 0.0
            || params.max_iters == 0
            || !params.sigma.is_finite()
            || params.sigma <= 0.0
            || controls.maximum_iterations == 0
            || controls.maximum_backtracks > 52
            || !controls.residual_tolerance.is_finite()
            || controls.residual_tolerance <= 0.0
        {
            return Err(FluxNewtonError::Invalid(
                "finite initial velocity, positive tolerances/penalty, nonzero solve budgets and 0..=52 backtracks required",
            ));
        }
        if admitted.is_none() {
            let boundaries: Vec<_> = self
                .mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(_, edge)| edge.tris.1 == usize::MAX)
                .map(|(edge, _)| StokesBoundary::Velocity {
                    edge,
                    value: boundary,
                })
                .collect();
            self.admit_stokes_boundaries(&boundaries)
                .map_err(FluxNewtonError::Boundary)?;
        }
        let mut state = vec![0.0; self.n];
        state[..self.n_u].copy_from_slice(initial_velocity);
        let mut backtracks = 0;
        let mut last_update = 0.0;
        let mut last_linear_iterations = 0;
        for iteration in 0..=controls.maximum_iterations {
            let linearization = self.linearize_rheology_admitted(
                params,
                law,
                force,
                boundary,
                &state,
                convection,
                inverse_timestep,
                previous_velocity,
                admitted,
            );
            if !source_is_valid() {
                return Err(FluxNewtonError::Invalid("constitutive source refused"));
            }
            let linearization = linearization.map_err(FluxNewtonError::Linearization)?;
            let defect = linearization
                .relative_residual
                .max(linearization.scaled_relative_residual);
            if defect.is_finite() && defect <= controls.residual_tolerance {
                return Ok(FluxNewtonStep {
                    solution: FluxSolution {
                        x: state,
                        n_u: self.n_u,
                        iters: last_linear_iterations,
                        rel_residual: linearization.relative_residual,
                        scaled_rel_residual: linearization.scaled_relative_residual,
                    },
                    iterations: iteration,
                    backtracks,
                    relative_update: last_update,
                    nonlinear_residual: defect,
                });
            }
            if iteration == controls.maximum_iterations {
                return Err(FluxNewtonError::NonConvergence {
                    iterations: iteration,
                    residual: defect,
                });
            }
            let scales: Vec<_> = (0..self.n)
                .map(|row| {
                    linearization
                        .jacobian
                        .row(row)
                        .1
                        .iter()
                        .fold(0.0_f64, |v, coefficient| v.max(coefficient.abs()))
                        .max(f64::MIN_POSITIVE)
                })
                .collect();
            let merit = weighted_norm(&linearization.residual, &scales);
            if !merit.is_finite() {
                return Err(FluxNewtonError::Invalid("nonfinite weighted residual"));
            }
            let rhs: Vec<_> = linearization.residual.iter().map(|r| -r).collect();
            let correction = self
                .solve_linear_result(&linearization.jacobian, &rhs, params)
                .map_err(|error| FluxNewtonError::LinearSolve {
                    iteration: iteration + 1,
                    error,
                })?;
            if !correction.converged(params.tol) || correction.x.iter().any(|v| !v.is_finite()) {
                return Err(FluxNewtonError::LinearSolve {
                    iteration: iteration + 1,
                    error: StokesBoundaryError::NotConverged {
                        iterations: correction.iters,
                        relative_residual: correction.rel_residual,
                        scaled_relative_residual: correction.scaled_rel_residual,
                    },
                });
            }
            last_linear_iterations = correction.iters;
            let mut alpha = 1.0;
            let mut accepted = None;
            for attempt in 0..=controls.maximum_backtracks {
                let candidate: Vec<_> = state
                    .iter()
                    .zip(&correction.x)
                    .map(|(x, step)| x + alpha * step)
                    .collect();
                let trial = self.newton_trial(
                    params,
                    law,
                    force,
                    boundary,
                    &candidate,
                    convection,
                    inverse_timestep,
                    previous_velocity,
                    admitted,
                );
                if !source_is_valid() {
                    return Err(FluxNewtonError::Invalid("constitutive source refused"));
                }
                if let Ok(trial) = trial {
                    let trial_defect = trial.original.max(trial.scaled);
                    let trial_merit = weighted_norm(&trial.residual, &scales);
                    if trial_merit <= (1.0 - 1e-4 * alpha) * merit
                        || (trial_defect.is_finite() && trial_defect <= controls.residual_tolerance)
                    {
                        accepted = Some(candidate);
                        break;
                    }
                }
                backtracks += 1;
                if attempt < controls.maximum_backtracks {
                    alpha *= 0.5;
                }
            }
            let next = accepted.ok_or(FluxNewtonError::LineSearch {
                iteration: iteration + 1,
                attempted_steps: controls.maximum_backtracks + 1,
                residual: defect,
            })?;
            let mut difference = 0.0_f64;
            let mut magnitude = 0.0_f64;
            let mut previous_magnitude = 0.0_f64;
            for i in 0..self.n_u {
                difference = difference.hypot(next[i] - state[i]);
                magnitude = magnitude.hypot(next[i]);
                previous_magnitude = previous_magnitude.hypot(state[i]);
            }
            last_update = difference / magnitude.max(previous_magnitude).max(f64::MIN_POSITIVE);
            state = next;
        }
        unreachable!("bounded Newton loop accepts or refuses its final state")
    }

    #[allow(clippy::too_many_arguments)]
    fn newton_trial(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundary: &dyn Fn([f64; 2]) -> [f64; 2],
        state: &[f64],
        convection: bool,
        inverse_timestep: f64,
        previous: Option<&[f64]>,
        admitted: Option<&AdmittedBoundary>,
    ) -> Result<Trial, FluxLinearizationError> {
        if state.iter().any(|v| !v.is_finite()) {
            return Err(FluxLinearizationError::Invalid("nonfinite Newton trial"));
        }
        let rates = self.shear_rates(state);
        if rates.iter().any(|v| !v.is_finite()) {
            return Err(FluxLinearizationError::Invalid("nonfinite trial strain"));
        }
        let failure = RefCell::new(None);
        let field = |cell, point| {
            if failure.borrow().is_some() {
                return 1.0;
            }
            match sample(law, cell, point, rates[cell]) {
                Ok(value) => value.nu,
                Err(error) => {
                    *failure.borrow_mut() = Some(error);
                    1.0
                }
            }
        };
        let (operator, rhs) = self.assemble_operator(
            params,
            Some(&field),
            force,
            boundary,
            convection.then_some(&state[..self.n_u]),
            inverse_timestep,
            previous,
            None,
            admitted,
            None,
            None,
        );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let (residual, original, scaled) =
            rheology_residual(self, &operator, &rhs, state, admitted.is_some());
        if residual.iter().any(|v| !v.is_finite()) {
            return Err(FluxLinearizationError::Invalid("nonfinite trial residual"));
        }
        Ok(Trial {
            residual,
            original,
            scaled,
        })
    }
}
