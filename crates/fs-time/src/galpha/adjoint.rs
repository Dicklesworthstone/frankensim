//! Discrete first-order generalized-alpha pullbacks for
//! `M(p) rate + r(t, u, p) = forcing`.
//!
//! The forward solve is the production Jansen--Whiting--Hulbert step. Its
//! converged residual is differentiated by the implicit function theorem;
//! Newton, globalization and Krylov iterations are not differentiated. Time,
//! step size, spectral radius and solver policies are fixed. Accuracy depends
//! on primal/adjoint residuals and on the caller's matched derivative actions.

use super::{
    FirstOrderProblem, FirstOrderState, ImplicitStepTelemetry, OperatorFirstOrderGeneralizedAlpha,
    TimeSolveError, first_order_poll,
};
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolveReport, StallDiagnosis};
use std::cell::RefCell;

/// One immutable parameter point with explicitly transposed derivatives.
/// Implementations must be pure during a pullback, overwrite every output,
/// and bound their own callback work. Neither mass nor tangent is assumed
/// symmetric. Mass may depend on parameters, but not on state or time, as in
/// [`FirstOrderProblem`].
pub trait FirstOrderVjp: FirstOrderProblem {
    /// Number of model parameters, which may be zero.
    fn parameter_count(&self) -> usize;
    /// Overwrite `output` with `M(p)^T seed`.
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String>;
    /// Overwrite `output` with `r_u(t,u,p)^T seed` at fixed time and parameters.
    fn tangent_transpose_apply(
        &self,
        t: f64,
        u: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String>;
    /// Overwrite `parameters` with
    /// `(d[M(p) rate + r(t,u,p)]/dp)^T seed`, holding `t`, `u`, `rate` fixed.
    /// Parameter-dependent mass is part of this derivative. An independent
    /// model must explicitly write zeros. The external forcing is excluded;
    /// the returned forcing cotangent lets its owner apply its own chain rule.
    fn residual_parameter_vjp(
        &self,
        t: f64,
        u: &[f64],
        rate: &[f64],
        seed: &[f64],
        parameters: &mut [f64],
    ) -> Result<(), String>;
}

/// Independent controls for the transposed effective-system solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FirstOrderAdjointConfig {
    /// Positive FGMRES restart length.
    pub restart: usize,
    /// Positive maximum number of restart cycles.
    pub max_cycles: usize,
    /// Finite positive true-relative-residual tolerance.
    pub tolerance: f64,
}

impl Default for FirstOrderAdjointConfig {
    fn default() -> Self {
        Self {
            restart: 24,
            max_cycles: 16,
            tolerance: 1.0e-11,
        }
    }
}

/// A refused pullback returns no partial gradient and never mutates its input.
#[derive(Debug, Clone)]
pub enum FirstOrderAdjointError {
    /// Inconsistent dimensions, nonfinite inputs or invalid solver controls.
    InvalidInput(&'static str),
    /// The admitted numerical scratch ceiling is too small.
    WorkspaceLimit {
        /// Conservative number of required scalar components.
        required: usize,
        /// Caller-supplied maximum number of scalar components.
        limit: usize,
    },
    /// The forward production step refused or cancellation was requested.
    Step(TimeSolveError),
    /// The transposed system did not converge within its configured work.
    NotConverged(SolveReport),
    /// A model derivative callback explicitly refused its evaluation.
    Derivative(String),
    /// A derivative callback left an entry unwritten or nonfinite.
    NonFiniteDerivative,
    /// Finite input values produced an unrepresentable intermediate/result.
    NonFiniteAccumulation,
}

impl From<TimeSolveError> for FirstOrderAdjointError {
    fn from(error: TimeSolveError) -> Self {
        Self::Step(error)
    }
}

impl std::fmt::Display for FirstOrderAdjointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "first-order generalized-alpha adjoint failed: {self:?}")
    }
}

impl std::error::Error for FirstOrderAdjointError {}

/// One endpoint plus the pullback of a caller-supplied `(state, rate)` seed.
#[derive(Debug, Clone)]
pub struct FirstOrderStepGradient {
    /// State at the end of the production step.
    pub value: Vec<f64>,
    /// Rate at the end of the production step.
    pub rate: Vec<f64>,
    /// Cotangent of the initial state, holding initial rate fixed.
    pub initial: Vec<f64>,
    /// Cotangent of the initial rate, holding initial state fixed.
    pub initial_rate: Vec<f64>,
    /// Model parameter cotangent, holding both initial vectors/forcing fixed.
    pub parameters: Vec<f64>,
    /// Cotangent of the supplied forcing at `t_n + alpha_f*h`.
    pub forcing: Vec<f64>,
    /// Complete forward Newton/Krylov convergence report.
    pub primal: ImplicitStepTelemetry,
    /// True-residual report for the transposed effective-system solve.
    pub adjoint: SolveReport,
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

/// Applies `[alpha_m/(gamma*h) M + alpha_f r_u]^T`. Fallible model
/// derivatives cannot be expressed by `LinearOp`; retain their first error and
/// poison the result so FGMRES stops without publishing a candidate gradient.
struct EffectiveTranspose<'a, M: ?Sized> {
    model: &'a M,
    t: f64,
    u: &'a [f64],
    mass_scale: f64,
    tangent_scale: f64,
    error: RefCell<Option<FirstOrderAdjointError>>,
}

impl<M: FirstOrderVjp + ?Sized> LinearOp for EffectiveTranspose<'_, M> {
    fn n(&self) -> usize {
        self.u.len()
    }

    fn apply(&self, seed: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.error.borrow().is_some() {
            return;
        }
        let mut mass = vec![f64::NAN; self.n()];
        let result = self
            .model
            .mass_transpose_apply(seed, &mut mass)
            .and_then(|()| {
                self.model
                    .tangent_transpose_apply(self.t, self.u, seed, output)
            });
        if let Err(error) = result {
            *self.error.borrow_mut() = Some(FirstOrderAdjointError::Derivative(error));
            output.fill(f64::NAN);
            return;
        }
        if !finite(&mass) || !finite(output) {
            *self.error.borrow_mut() = Some(FirstOrderAdjointError::NonFiniteDerivative);
            output.fill(f64::NAN);
            return;
        }
        for (value, mass_value) in output.iter_mut().zip(mass) {
            *value = self
                .mass_scale
                .mul_add(mass_value, self.tangent_scale * *value);
        }
        if !finite(output) {
            *self.error.borrow_mut() = Some(FirstOrderAdjointError::NonFiniteAccumulation);
        }
    }
}

impl OperatorFirstOrderGeneralizedAlpha {
    /// Conservative live scalar-storage ceiling for `step_vjp`, including
    /// both Krylov workspaces, temporary forward state, Newton/linear history,
    /// derivative scratch and returned gradients. Existing input history is
    /// borrowed, never cloned. Allocator metadata and callback-owned storage
    /// are outside this bound. `None` means checked size arithmetic overflow.
    #[must_use]
    pub fn adjoint_workspace_components(
        &self,
        parameters: usize,
        config: FirstOrderAdjointConfig,
    ) -> Option<usize> {
        let primal_restart = self.solve.newton.linear_restart;
        let adjoint_restart = config.restart;
        let restarts = primal_restart.checked_add(adjoint_restart)?;
        let vectors = restarts
            .checked_mul(2)?
            .checked_add(40)?
            .checked_mul(self.n)?;
        let primal_hessenberg = primal_restart.checked_add(1)?.checked_mul(primal_restart)?;
        let adjoint_hessenberg = adjoint_restart
            .checked_add(1)?
            .checked_mul(adjoint_restart)?;
        vectors
            .checked_add(primal_hessenberg)?
            .checked_add(adjoint_hessenberg)?
            .checked_add(restarts.checked_mul(8)?)?
            .checked_add(self.solve.max_newton_iterations.checked_mul(32)?)?
            .checked_add(self.solve.newton.max_linear_cycles.checked_mul(6)?)?
            .checked_add(config.max_cycles.checked_mul(4)?)?
            .checked_add(parameters.checked_mul(2)?)?
            .checked_add(64)
    }

    /// Pull back a production first-order generalized-alpha step.
    ///
    /// The adjoint preconditioner must approximate the transposed effective
    /// system; it is never implicitly assumed symmetric. The primal uses
    /// the model's forward preconditioner hook (identity by default). Endpoint
    /// state AND rate seeds are needed when chaining multiple steps. If the
    /// initial rate was computed from the initial state or model parameters,
    /// its returned cotangent must also enter that initialization's chain rule.
    /// Direct objective terms likewise belong to the caller.
    ///
    /// Cancellation is polled at forward Newton-attempt boundaries, between
    /// adjoint restart cycles, around final derivative calls and before return.
    /// A single callback, Newton attempt or Krylov cycle is not preemptible.
    /// No adaptive-controller, event, timestep or spectral-radius derivative
    /// and no interval enclosure of the gradient is claimed.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_vjp<M, P, Cancel>(
        &self,
        initial: &FirstOrderState,
        model: &M,
        forcing: &[f64],
        terminal: (&[f64], &[f64]),
        adjoint_preconditioner: &P,
        config: FirstOrderAdjointConfig,
        max_workspace_components: usize,
        cancelled: &mut Cancel,
    ) -> Result<FirstOrderStepGradient, FirstOrderAdjointError>
    where
        M: FirstOrderVjp + ?Sized,
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        first_order_poll(cancelled)?;
        if model.dimension() != self.n
            || [
                initial.u.as_slice(),
                initial.rate.as_slice(),
                forcing,
                terminal.0,
                terminal.1,
            ]
            .iter()
            .any(|values| values.len() != self.n || !finite(values))
        {
            return Err(FirstOrderAdjointError::InvalidInput(
                "finite dimension-matched initial state, rate, forcing and seeds required",
            ));
        }
        if config.restart == 0
            || config.max_cycles == 0
            || !config.tolerance.is_finite()
            || config.tolerance <= 0.0
        {
            return Err(FirstOrderAdjointError::InvalidInput(
                "positive finite adjoint tolerance and positive work limits required",
            ));
        }
        let p = model.parameter_count();
        let required = self.adjoint_workspace_components(p, config).ok_or(
            FirstOrderAdjointError::InvalidInput("workspace dimension overflow"),
        )?;
        if required > max_workspace_components {
            return Err(FirstOrderAdjointError::WorkspaceLimit {
                required,
                limit: max_workspace_components,
            });
        }
        let a = 1.0 / (self.gamma * self.h);
        let b = (1.0 - self.gamma) / self.gamma;
        if !a.is_finite() {
            return Err(FirstOrderAdjointError::InvalidInput(
                "finite reciprocal of gamma*h required",
            ));
        }
        // Deliberately avoid cloning the caller's arbitrarily long history.
        let mut next = FirstOrderState::new(initial.t, &initial.u, &initial.rate);
        next.steps = initial.steps;
        let primal = self.step_controlled(&mut next, model, forcing, cancelled)?;
        if !finite(&next.u) || !finite(&next.rate) {
            return Err(FirstOrderAdjointError::NonFiniteAccumulation);
        }
        let t_eval = self.h.mul_add(self.alpha_f, initial.t);
        let mut u_eval = vec![0.0; self.n];
        let mut rate_eval = vec![0.0; self.n];
        let mut rhs = vec![0.0; self.n];
        for i in 0..self.n {
            u_eval[i] = self
                .alpha_f
                .mul_add(next.u[i], (1.0 - self.alpha_f) * initial.u[i]);
            rate_eval[i] = self
                .alpha_m
                .mul_add(next.rate[i], (1.0 - self.alpha_m) * initial.rate[i]);
            rhs[i] = a.mul_add(terminal.1[i], terminal.0[i]);
        }
        if !finite(&u_eval) || !finite(&rate_eval) || !finite(&rhs) {
            return Err(FirstOrderAdjointError::NonFiniteAccumulation);
        }
        let operator = EffectiveTranspose {
            model,
            t: t_eval,
            u: &u_eval,
            mass_scale: self.alpha_m / (self.gamma * self.h),
            tangent_scale: self.alpha_f,
            error: RefCell::new(None),
        };
        let mut solve = FgmresState::new(&rhs, config.restart);
        let mut report = None;
        for cycle in 0..config.max_cycles {
            first_order_poll(cancelled)?;
            let candidate = solve.run(&operator, adjoint_preconditioner, &rhs, config.tolerance, 1);
            first_order_poll(cancelled)?;
            if let Some(error) = operator.error.borrow_mut().take() {
                return Err(error);
            }
            if candidate.converged {
                report = Some(candidate);
                break;
            }
            if candidate.diagnosis == Some(StallDiagnosis::Breakdown)
                || cycle + 1 == config.max_cycles
            {
                return Err(FirstOrderAdjointError::NotConverged(candidate));
            }
        }
        let adjoint = report.expect("positive cycle budget either converges or returns an error");
        let lambda = solve.x;
        let mut mass_bar = vec![f64::NAN; self.n];
        let mut internal_bar = vec![f64::NAN; self.n];
        let mut parameters = vec![f64::NAN; p];
        first_order_poll(cancelled)?;
        let result = model.mass_transpose_apply(&lambda, &mut mass_bar);
        first_order_poll(cancelled)?;
        result.map_err(FirstOrderAdjointError::Derivative)?;
        let result = model.tangent_transpose_apply(t_eval, &u_eval, &lambda, &mut internal_bar);
        first_order_poll(cancelled)?;
        result.map_err(FirstOrderAdjointError::Derivative)?;
        let result =
            model.residual_parameter_vjp(t_eval, &u_eval, &rate_eval, &lambda, &mut parameters);
        first_order_poll(cancelled)?;
        result.map_err(FirstOrderAdjointError::Derivative)?;
        if !finite(&mass_bar) || !finite(&internal_bar) || !finite(&parameters) {
            return Err(FirstOrderAdjointError::NonFiniteDerivative);
        }
        // R_u0 = -alpha_m*a*M + (1-alpha_f)*r_u and
        // R_rate0 = (1-alpha_m/gamma)*M. The direct endpoint-rate
        // derivatives contribute -a*bar_rate1 and -b*bar_rate1 respectively.
        let mut initial_bar = vec![0.0; self.n];
        let mut initial_rate_bar = vec![0.0; self.n];
        for i in 0..self.n {
            initial_bar[i] = -a * terminal.1[i] + self.alpha_m * a * mass_bar[i]
                - (1.0 - self.alpha_f) * internal_bar[i];
            initial_rate_bar[i] =
                -b * terminal.1[i] - (1.0 - self.alpha_m / self.gamma) * mass_bar[i];
        }
        for value in &mut parameters {
            *value = -*value;
        }
        if !finite(&initial_bar) || !finite(&initial_rate_bar) || !finite(&lambda) {
            return Err(FirstOrderAdjointError::NonFiniteAccumulation);
        }
        first_order_poll(cancelled)?;
        Ok(FirstOrderStepGradient {
            value: next.u,
            rate: next.rate,
            initial: initial_bar,
            initial_rate: initial_rate_bar,
            parameters,
            forcing: lambda,
            primal,
            adjoint,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galpha::ImplicitSolveConfig;

    struct ScalarModel {
        derivative_mode: u8,
    }

    impl FirstOrderProblem for ScalarModel {
        fn dimension(&self) -> usize {
            1
        }
        fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
            output[0] = 1.5 * input[0];
        }
        fn internal_force(&self, t: f64, u: &[f64], output: &mut [f64]) {
            output[0] = u[0] + 0.2 * u[0] * u[0] + t;
        }
        fn tangent_apply(&self, _t: f64, u: &[f64], input: &[f64], output: &mut [f64]) {
            output[0] = (1.0 + 0.4 * u[0]) * input[0];
        }
    }

    impl FirstOrderVjp for ScalarModel {
        fn parameter_count(&self) -> usize {
            0
        }
        fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
            if self.derivative_mode == 1 {
                return Err("mass derivative unavailable".to_owned());
            }
            if self.derivative_mode != 2 {
                output[0] = 1.5 * seed[0];
            }
            Ok(())
        }
        fn tangent_transpose_apply(
            &self,
            _t: f64,
            u: &[f64],
            seed: &[f64],
            output: &mut [f64],
        ) -> Result<(), String> {
            output[0] = (1.0 + 0.4 * u[0]) * seed[0];
            Ok(())
        }
        fn residual_parameter_vjp(
            &self,
            _t: f64,
            _u: &[f64],
            _rate: &[f64],
            _seed: &[f64],
            parameters: &mut [f64],
        ) -> Result<(), String> {
            if self.derivative_mode == 3 {
                return Err("parameter derivative unavailable".to_owned());
            }
            parameters.fill(0.0);
            Ok(())
        }
    }

    struct Identity;
    impl FlexiblePreconditioner for Identity {
        fn apply(&self, _iteration: usize, input: &[f64], output: &mut [f64]) {
            output.copy_from_slice(input);
        }
    }

    struct Zero;
    impl FlexiblePreconditioner for Zero {
        fn apply(&self, _iteration: usize, _input: &[f64], output: &mut [f64]) {
            output.fill(0.0);
        }
    }

    // G4: every reachable cancellation boundary refuses atomically, including
    // a previously accepted history, and the unmodified input can be retried.
    #[test]
    fn cancellation_preserves_first_order_forward_and_reverse_state() {
        let model = ScalarModel { derivative_mode: 0 };
        let method =
            OperatorFirstOrderGeneralizedAlpha::new(1, 0.1, 0.35, ImplicitSolveConfig::default());
        let config = FirstOrderAdjointConfig::default();
        let budget = method.adjoint_workspace_components(0, config).unwrap();
        let mut initial = FirstOrderState::new(0.0, &[0.8], &[-0.4]);
        method.step(&mut initial, &model, &[0.2]).unwrap();
        let preserved = initial.clone();

        let mut forward = initial.clone();
        let mut forward_polls = 0;
        method
            .step_controlled(&mut forward, &model, &[0.2], &mut || {
                forward_polls += 1;
                false
            })
            .unwrap();
        let mut ordinary = initial.clone();
        method.step(&mut ordinary, &model, &[0.2]).unwrap();
        assert_eq!(forward, ordinary);
        for stop in 1..=forward_polls {
            let mut attempted = initial.clone();
            let mut polls = 0;
            let result = method.step_controlled(&mut attempted, &model, &[0.2], &mut || {
                polls += 1;
                polls == stop
            });
            assert_eq!(
                result,
                Err(TimeSolveError::Cancelled),
                "forward poll {stop}"
            );
            assert_eq!(attempted, preserved);
        }

        let mut reverse_polls = 0;
        let reference = method
            .step_vjp(
                &initial,
                &model,
                &[0.2],
                (&[0.7], &[-0.25]),
                &Identity,
                config,
                budget,
                &mut || {
                    reverse_polls += 1;
                    false
                },
            )
            .unwrap();
        assert_eq!(reference.value, ordinary.u);
        assert_eq!(reference.rate, ordinary.rate);
        assert_eq!(reference.parameters, [] as [f64; 0]);
        for stop in 1..=reverse_polls {
            let mut polls = 0;
            let result = method.step_vjp(
                &initial,
                &model,
                &[0.2],
                (&[0.7], &[-0.25]),
                &Identity,
                config,
                budget,
                &mut || {
                    polls += 1;
                    polls == stop
                },
            );
            assert!(
                matches!(
                    result,
                    Err(FirstOrderAdjointError::Step(TimeSolveError::Cancelled))
                ),
                "reverse poll {stop}: {result:?}"
            );
            assert_eq!(initial, preserved);
        }
        let retry = method
            .step_vjp(
                &initial,
                &model,
                &[0.2],
                (&[0.7], &[-0.25]),
                &Identity,
                config,
                budget,
                &mut || false,
            )
            .unwrap();
        assert_eq!(format!("{reference:?}"), format!("{retry:?}"));
    }

    // G0/G4: resource admission, malformed seeds, failed/missing derivatives,
    // exhausted primal/adjoint solves and zero-seed success remain distinct.
    #[test]
    #[allow(clippy::too_many_lines)]
    fn first_order_adjoint_refuses_incomplete_derivatives_and_exhausted_budgets() {
        let method =
            OperatorFirstOrderGeneralizedAlpha::new(1, 0.1, 0.35, ImplicitSolveConfig::default());
        let config = FirstOrderAdjointConfig::default();
        let budget = method.adjoint_workspace_components(0, config).unwrap();
        let model = ScalarModel { derivative_mode: 0 };
        let initial = FirstOrderState::new(0.0, &[0.8], &[-0.4]);
        let result = method.step_vjp(
            &initial,
            &model,
            &[0.2],
            (&[0.7], &[-0.25]),
            &Identity,
            config,
            budget - 1,
            &mut || false,
        );
        assert!(
            matches!(result, Err(FirstOrderAdjointError::WorkspaceLimit { required, .. }) if required == budget)
        );
        assert!(
            method
                .adjoint_workspace_components(usize::MAX, config)
                .is_none()
        );
        let result = method.step_vjp(
            &initial,
            &model,
            &[0.2],
            (&[f64::NAN], &[-0.25]),
            &Identity,
            config,
            budget,
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(FirstOrderAdjointError::InvalidInput(_))
        ));
        let result = method.step_vjp(
            &initial,
            &model,
            &[0.2],
            (&[0.7], &[-0.25]),
            &Identity,
            FirstOrderAdjointConfig {
                restart: 0,
                ..config
            },
            budget,
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(FirstOrderAdjointError::InvalidInput(_))
        ));
        for mode in 1..=3 {
            let result = method.step_vjp(
                &initial,
                &ScalarModel {
                    derivative_mode: mode,
                },
                &[0.2],
                (&[0.7], &[-0.25]),
                &Identity,
                config,
                budget,
                &mut || false,
            );
            if mode == 2 {
                assert!(matches!(
                    result,
                    Err(FirstOrderAdjointError::NonFiniteDerivative)
                ));
            } else {
                assert!(matches!(result, Err(FirstOrderAdjointError::Derivative(_))));
            }
        }
        let result = method.step_vjp(
            &initial,
            &model,
            &[0.2],
            (&[0.7], &[-0.25]),
            &Zero,
            config,
            budget,
            &mut || false,
        );
        assert!(
            matches!(result, Err(FirstOrderAdjointError::NotConverged(report)) if !report.converged)
        );
        let limited = OperatorFirstOrderGeneralizedAlpha::new(
            1,
            0.1,
            0.35,
            ImplicitSolveConfig {
                max_newton_iterations: 1,
                ..ImplicitSolveConfig::default()
            },
        );
        let result = limited.step_vjp(
            &initial,
            &model,
            &[0.2],
            (&[0.7], &[-0.25]),
            &Identity,
            config,
            budget,
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(FirstOrderAdjointError::Step(TimeSolveError::NotConverged(
                _
            )))
        ));
        let zero = method
            .step_vjp(
                &initial,
                &model,
                &[0.2],
                (&[0.0], &[0.0]),
                &Identity,
                config,
                budget,
                &mut || false,
            )
            .unwrap();
        assert_eq!(zero.initial, [0.0]);
        assert_eq!(zero.initial_rate, [0.0]);
        assert_eq!(zero.forcing, [0.0]);
        assert!(zero.adjoint.converged);
    }
}
