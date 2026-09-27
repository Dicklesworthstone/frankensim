//! Structural generalized-alpha pullbacks for `M(p) a + C(p) v + r(q,p) = f`.
//!
//! Differentiate the converged Chung--Hulbert residual and the Newmark
//! correctors, including all displacement, velocity and acceleration seeds.
//! The production Newton step is reused; Newton/globalization/Krylov iteration
//! paths are not differentiated. Accuracy is limited by the reported primal
//! and adjoint residuals and the caller's matched derivative actions.

use super::{
    ImplicitStepTelemetry, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    TimeSolveError, structural_poll,
};
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolveReport, StallDiagnosis};
use std::cell::RefCell;

/// A fixed structural parameter point with explicit transposed actions.
/// Implementations must be pure during a call and overwrite every output.
/// Mass and damping may depend on parameters, but not on state/time as in
/// [`SecondOrderProblem`]. No symmetry is assumed for M, C, or the tangent.
pub trait SecondOrderVjp: SecondOrderProblem {
    /// Number of model parameters, which may be zero.
    fn parameter_count(&self) -> usize;
    /// Overwrite `output` with `M(p)^T seed`.
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String>;
    /// Overwrite `output` with `C(p)^T seed`.
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String>;
    /// Overwrite `output` with `r_q(q,p)^T seed`.
    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String>;
    /// Overwrite `parameters` with
    /// `(d[M(p) a + C(p) v + r(q,p)]/dp)^T seed`, holding q/v/a fixed.
    /// Both parameter-dependent mass and damping are part of this derivative.
    /// An independent model explicitly writes zeros. External forcing is
    /// excluded; its owner can chain the separately returned forcing seed.
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        parameters: &mut [f64],
    ) -> Result<(), String>;
}

/// Independent bounded controls for the transposed effective-system solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondOrderAdjointConfig {
    /// Positive FGMRES restart length.
    pub restart: usize,
    /// Positive maximum number of restart cycles.
    pub max_cycles: usize,
    /// Finite positive true-relative-residual tolerance.
    pub tolerance: f64,
}

impl Default for SecondOrderAdjointConfig {
    fn default() -> Self {
        Self {
            restart: 24,
            max_cycles: 16,
            tolerance: 1.0e-11,
        }
    }
}

/// A refused structural pullback never returns partial gradients.
#[derive(Debug, Clone)]
pub enum SecondOrderAdjointError {
    /// Inconsistent dimensions, nonfinite data or invalid controls.
    InvalidInput(&'static str),
    /// Numerical scratch requirements exceed the admitted ceiling.
    WorkspaceLimit {
        /// Conservative required number of scalar components.
        required: usize,
        /// Supplied scalar-component ceiling.
        limit: usize,
    },
    /// The production step refused, or cancellation was requested.
    Step(TimeSolveError),
    /// Transposed FGMRES did not meet its true-residual target.
    NotConverged(SolveReport),
    /// A derivative callback refused its evaluation.
    Derivative(String),
    /// A derivative output was left unwritten or nonfinite.
    NonFiniteDerivative,
    /// Finite inputs produced unrepresentable arithmetic.
    NonFiniteAccumulation,
}

impl From<TimeSolveError> for SecondOrderAdjointError {
    fn from(error: TimeSolveError) -> Self {
        Self::Step(error)
    }
}
impl std::fmt::Display for SecondOrderAdjointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "structural generalized-alpha adjoint failed: {self:?}")
    }
}
impl std::error::Error for SecondOrderAdjointError {}

/// Endpoint and pullback of the supplied terminal `(q, v, a)` cotangent.
#[derive(Debug, Clone)]
pub struct SecondOrderStepGradient {
    /// Endpoint displacement from the production step.
    pub q: Vec<f64>,
    /// Endpoint velocity from the production Newmark corrector.
    pub v: Vec<f64>,
    /// Endpoint acceleration from the production Newmark corrector.
    pub a: Vec<f64>,
    /// Initial displacement cotangent, holding initial v/a fixed.
    pub initial_q: Vec<f64>,
    /// Initial velocity cotangent, holding initial q/a fixed.
    pub initial_v: Vec<f64>,
    /// Initial acceleration cotangent, holding initial q/v fixed.
    pub initial_a: Vec<f64>,
    /// Model parameter cotangent, holding initial state and forcing fixed.
    pub parameters: Vec<f64>,
    /// Cotangent of supplied forcing at `t_n + (1-alpha_f)*h`.
    pub forcing: Vec<f64>,
    /// Complete production Newton/Krylov report.
    pub primal: ImplicitStepTelemetry,
    /// Transposed effective-system true-residual report.
    pub adjoint: SolveReport,
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}

struct StructuralTranspose<'a, M: ?Sized> {
    model: &'a M,
    q: &'a [f64],
    mass_scale: f64,
    damping_scale: f64,
    tangent_scale: f64,
    error: RefCell<Option<SecondOrderAdjointError>>,
}

impl<M: SecondOrderVjp + ?Sized> LinearOp for StructuralTranspose<'_, M> {
    fn n(&self) -> usize {
        self.q.len()
    }

    fn apply(&self, seed: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.error.borrow().is_some() {
            return;
        }
        let mut mass = vec![f64::NAN; self.n()];
        let mut damping = vec![f64::NAN; self.n()];
        let result = self
            .model
            .mass_transpose_apply(seed, &mut mass)
            .and_then(|()| self.model.damping_transpose_apply(seed, &mut damping))
            .and_then(|()| self.model.tangent_transpose_apply(self.q, seed, output));
        if let Err(error) = result {
            *self.error.borrow_mut() = Some(SecondOrderAdjointError::Derivative(error));
            output.fill(f64::NAN);
            return;
        }
        if !finite(&mass) || !finite(&damping) || !finite(output) {
            *self.error.borrow_mut() = Some(SecondOrderAdjointError::NonFiniteDerivative);
            output.fill(f64::NAN);
            return;
        }
        for i in 0..self.n() {
            output[i] = self.mass_scale.mul_add(
                mass[i],
                self.damping_scale
                    .mul_add(damping[i], self.tangent_scale * output[i]),
            );
        }
        if !finite(output) {
            *self.error.borrow_mut() = Some(SecondOrderAdjointError::NonFiniteAccumulation);
        }
    }
}

impl OperatorGeneralizedAlpha {
    /// Conservative numerical scalar-storage bound, including primal/adjoint
    /// Krylov work, temporary forward state, Newton/linear histories, derivative
    /// scratch and returned gradients. The caller's existing state/history is
    /// borrowed. Allocator metadata and callback-owned storage are excluded.
    /// `None` indicates checked size arithmetic overflow.
    #[must_use]
    pub fn adjoint_workspace_components(
        &self,
        parameters: usize,
        config: SecondOrderAdjointConfig,
    ) -> Option<usize> {
        let primal_restart = self.solve.newton.linear_restart;
        let restarts = primal_restart.checked_add(config.restart)?;
        restarts
            .checked_mul(2)?
            .checked_add(52)?
            .checked_mul(self.n)?
            .checked_add(primal_restart.checked_add(1)?.checked_mul(primal_restart)?)?
            .checked_add(config.restart.checked_add(1)?.checked_mul(config.restart)?)?
            .checked_add(restarts.checked_mul(8)?)?
            .checked_add(self.solve.max_newton_iterations.checked_mul(32)?)?
            .checked_add(self.solve.newton.max_linear_cycles.checked_mul(6)?)?
            .checked_add(config.max_cycles.checked_mul(4)?)?
            .checked_add(parameters.checked_mul(2)?)?
            .checked_add(64)
    }

    /// Pull back the production structural generalized-alpha step.
    ///
    /// The supplied adjoint preconditioner suits the TRANSPOSED effective
    /// matrix; no symmetry is assumed. The primal uses the existing Newton
    /// solver's identity inner preconditioner. All three endpoint seeds matter
    /// when chaining a trajectory. If initial acceleration was computed from
    /// q/v/parameters, chain its cotangent through that consistency solve;
    /// direct objective and forcing-model terms likewise belong to the caller.
    ///
    /// Time, timestep, spectral radius and solver policies are held fixed.
    /// Cancellation is checked at primal Newton-attempt boundaries, adjoint
    /// restart boundaries, around final derivatives and before return. Work
    /// inside one callback/attempt/cycle is not preempted. Failed solves or
    /// derivatives and cancellation return no partial gradient. No controller,
    /// event or interval-gradient guarantee is claimed.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_vjp<M, P, Cancel>(
        &self,
        initial: &SecondOrderState,
        model: &M,
        forcing: &[f64],
        terminal: (&[f64], &[f64], &[f64]),
        adjoint_preconditioner: &P,
        config: SecondOrderAdjointConfig,
        max_workspace_components: usize,
        cancelled: &mut Cancel,
    ) -> Result<SecondOrderStepGradient, SecondOrderAdjointError>
    where
        M: SecondOrderVjp + ?Sized,
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        structural_poll(cancelled)?;
        if model.dimension() != self.n
            || [
                initial.q.as_slice(),
                initial.v.as_slice(),
                initial.a.as_slice(),
                forcing,
                terminal.0,
                terminal.1,
                terminal.2,
            ]
            .iter()
            .any(|values| values.len() != self.n || !finite(values))
        {
            return Err(SecondOrderAdjointError::InvalidInput(
                "finite dimension-matched q/v/a, forcing and terminal seeds required",
            ));
        }
        if config.restart == 0
            || config.max_cycles == 0
            || !config.tolerance.is_finite()
            || config.tolerance <= 0.0
        {
            return Err(SecondOrderAdjointError::InvalidInput(
                "positive finite adjoint tolerance and positive work budgets required",
            ));
        }
        let p = model.parameter_count();
        let required = self.adjoint_workspace_components(p, config).ok_or(
            SecondOrderAdjointError::InvalidInput("workspace dimension overflow"),
        )?;
        if required > max_workspace_components {
            return Err(SecondOrderAdjointError::WorkspaceLimit {
                required,
                limit: max_workspace_components,
            });
        }
        let acceleration_q = 1.0 / (self.beta * self.h * self.h);
        let acceleration_v = 1.0 / (self.beta * self.h);
        let acceleration_a = 0.5 / self.beta - 1.0;
        let velocity_a = self.gamma * self.h;
        if !finite(&[acceleration_q, acceleration_v, acceleration_a, velocity_a]) {
            return Err(SecondOrderAdjointError::InvalidInput(
                "finite Newmark derivative coefficients required",
            ));
        }
        // Reuse the actual nonlinear residual and correctors. Copy only the
        // three live state vectors, never the caller's trajectory history.
        let mut next = SecondOrderState::new(initial.t, &initial.q, &initial.v, &initial.a);
        next.steps = initial.steps;
        let primal = self.step_controlled(&mut next, model, forcing, cancelled)?;
        if !finite(&next.q) || !finite(&next.v) || !finite(&next.a) {
            return Err(SecondOrderAdjointError::NonFiniteAccumulation);
        }
        let mut q_eval = vec![0.0; self.n];
        let mut v_eval = vec![0.0; self.n];
        let mut a_eval = vec![0.0; self.n];
        let mut rhs = vec![0.0; self.n];
        for i in 0..self.n {
            q_eval[i] = (1.0 - self.alpha_f).mul_add(next.q[i], self.alpha_f * initial.q[i]);
            v_eval[i] = (1.0 - self.alpha_f).mul_add(next.v[i], self.alpha_f * initial.v[i]);
            a_eval[i] = (1.0 - self.alpha_m).mul_add(next.a[i], self.alpha_m * initial.a[i]);
            rhs[i] = acceleration_q.mul_add(
                velocity_a.mul_add(terminal.1[i], terminal.2[i]),
                terminal.0[i],
            );
        }
        if !finite(&q_eval) || !finite(&v_eval) || !finite(&a_eval) || !finite(&rhs) {
            return Err(SecondOrderAdjointError::NonFiniteAccumulation);
        }
        let operator = StructuralTranspose {
            model,
            q: &q_eval,
            mass_scale: (1.0 - self.alpha_m) / (self.beta * self.h * self.h),
            damping_scale: (1.0 - self.alpha_f) * self.gamma / (self.beta * self.h),
            tangent_scale: 1.0 - self.alpha_f,
            error: RefCell::new(None),
        };
        let mut solve = FgmresState::new(&rhs, config.restart);
        let mut report = None;
        for cycle in 0..config.max_cycles {
            structural_poll(cancelled)?;
            let candidate = solve.run(&operator, adjoint_preconditioner, &rhs, config.tolerance, 1);
            structural_poll(cancelled)?;
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
                return Err(SecondOrderAdjointError::NotConverged(candidate));
            }
        }
        let adjoint = report.expect("positive cycle budget either converges or refuses");
        let lambda = solve.x;
        let mut mass_bar = vec![f64::NAN; self.n];
        let mut damping_bar = vec![f64::NAN; self.n];
        let mut internal_bar = vec![f64::NAN; self.n];
        let mut parameters = vec![f64::NAN; p];
        structural_poll(cancelled)?;
        let result = model.mass_transpose_apply(&lambda, &mut mass_bar);
        structural_poll(cancelled)?;
        result.map_err(SecondOrderAdjointError::Derivative)?;
        let result = model.damping_transpose_apply(&lambda, &mut damping_bar);
        structural_poll(cancelled)?;
        result.map_err(SecondOrderAdjointError::Derivative)?;
        let result = model.tangent_transpose_apply(&q_eval, &lambda, &mut internal_bar);
        structural_poll(cancelled)?;
        result.map_err(SecondOrderAdjointError::Derivative)?;
        let result =
            model.residual_parameter_vjp(&q_eval, &v_eval, &a_eval, &lambda, &mut parameters);
        structural_poll(cancelled)?;
        result.map_err(SecondOrderAdjointError::Derivative)?;
        if !finite(&mass_bar)
            || !finite(&damping_bar)
            || !finite(&internal_bar)
            || !finite(&parameters)
        {
            return Err(SecondOrderAdjointError::NonFiniteDerivative);
        }
        let mut initial_q = vec![0.0; self.n];
        let mut initial_v = vec![0.0; self.n];
        let mut initial_a = vec![0.0; self.n];
        for i in 0..self.n {
            // Reverse the residual's intermediate-state interpolation, then
            // v1=v0+h[(1-gamma)a0+gamma*a1] and
            // a1=A(q1-q0)-B*v0-D*a0. This preserves all three initial seeds.
            let v_bar = terminal.1[i] - (1.0 - self.alpha_f) * damping_bar[i];
            let a_bar = terminal.2[i] - (1.0 - self.alpha_m) * mass_bar[i] + velocity_a * v_bar;
            initial_q[i] = -self.alpha_f * internal_bar[i] - acceleration_q * a_bar;
            initial_v[i] = -self.alpha_f * damping_bar[i] + v_bar - acceleration_v * a_bar;
            initial_a[i] = -self.alpha_m * mass_bar[i] + self.h * (1.0 - self.gamma) * v_bar
                - acceleration_a * a_bar;
        }
        for value in &mut parameters {
            *value = -*value;
        }
        if !finite(&initial_q) || !finite(&initial_v) || !finite(&initial_a) || !finite(&lambda) {
            return Err(SecondOrderAdjointError::NonFiniteAccumulation);
        }
        structural_poll(cancelled)?;
        Ok(SecondOrderStepGradient {
            q: next.q,
            v: next.v,
            a: next.a,
            initial_q,
            initial_v,
            initial_a,
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

    struct Oscillator(u8);
    impl SecondOrderProblem for Oscillator {
        fn dimension(&self) -> usize {
            1
        }
        fn mass_apply(&self, x: &[f64], out: &mut [f64]) {
            out[0] = 1.4 * x[0];
        }
        fn damping_apply(&self, x: &[f64], out: &mut [f64]) {
            out[0] = 0.25 * x[0];
        }
        fn internal_force(&self, q: &[f64], out: &mut [f64]) {
            out[0] = 2.0 * q[0] + 0.3 * q[0] * q[0] * q[0];
        }
        fn tangent_apply(&self, q: &[f64], x: &[f64], out: &mut [f64]) {
            out[0] = (2.0 + 0.9 * q[0] * q[0]) * x[0];
        }
    }
    impl SecondOrderVjp for Oscillator {
        fn parameter_count(&self) -> usize {
            1
        }
        fn mass_transpose_apply(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
            if self.0 == 1 {
                return Err("mass transpose unavailable".to_owned());
            }
            self.mass_apply(seed, out);
            Ok(())
        }
        fn damping_transpose_apply(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
            if self.0 != 2 {
                self.damping_apply(seed, out);
            }
            Ok(())
        }
        fn tangent_transpose_apply(
            &self,
            q: &[f64],
            seed: &[f64],
            out: &mut [f64],
        ) -> Result<(), String> {
            self.tangent_apply(q, seed, out);
            if self.0 == 3 {
                out[0] = f64::INFINITY;
            }
            Ok(())
        }
        fn residual_parameter_vjp(
            &self,
            q: &[f64],
            _v: &[f64],
            _a: &[f64],
            seed: &[f64],
            out: &mut [f64],
        ) -> Result<(), String> {
            if self.0 == 4 {
                return Err("stiffness derivative unavailable".to_owned());
            }
            if self.0 != 5 {
                out[0] = q[0] * seed[0];
            }
            Ok(())
        }
    }

    struct Preconditioner(f64);
    impl FlexiblePreconditioner for Preconditioner {
        fn apply(&self, _iteration: usize, input: &[f64], out: &mut [f64]) {
            out[0] = self.0 * input[0];
        }
    }
    fn config() -> SecondOrderAdjointConfig {
        SecondOrderAdjointConfig {
            restart: 1,
            max_cycles: 3,
            tolerance: 1.0e-12,
        }
    }
    fn pull<Cancel: FnMut() -> bool>(
        method: &OperatorGeneralizedAlpha,
        initial: &SecondOrderState,
        model: &Oscillator,
        preconditioner: f64,
        cancelled: &mut Cancel,
    ) -> Result<SecondOrderStepGradient, SecondOrderAdjointError> {
        method.step_vjp(
            initial,
            model,
            &[0.1],
            (&[0.2], &[0.3], &[-0.15]),
            &Preconditioner(preconditioner),
            config(),
            method.adjoint_workspace_components(1, config()).unwrap(),
            cancelled,
        )
    }

    // G4: cancellation at every reached boundary preserves a nonempty forward
    // trajectory; the immutable reverse input and all three vectors are retryable.
    #[test]
    fn structural_cancellation_is_atomic_for_forward_and_adjoint() {
        let method = OperatorGeneralizedAlpha::new(1, 0.15, 0.35, ImplicitSolveConfig::default());
        let model = Oscillator(0);
        let mut initial = SecondOrderState::new(0.0, &[0.7], &[0.2], &[-0.3]);
        method.step(&mut initial, &model, &[0.1]).unwrap();
        let original = initial.clone();
        let mut ordinary = initial.clone();
        method.step(&mut ordinary, &model, &[0.1]).unwrap();
        let mut forward_polls = 0;
        let mut positive = initial.clone();
        method
            .step_controlled(&mut positive, &model, &[0.1], &mut || {
                forward_polls += 1;
                false
            })
            .unwrap();
        assert_eq!(positive, ordinary);
        for stop in 1..=forward_polls {
            let mut attempted = initial.clone();
            let mut count = 0;
            let result = method.step_controlled(&mut attempted, &model, &[0.1], &mut || {
                count += 1;
                count == stop
            });
            assert_eq!(
                result,
                Err(TimeSolveError::Cancelled),
                "forward boundary {stop}"
            );
            assert_eq!(attempted, original);
        }
        let mut reverse_polls = 0;
        let gradient = pull(&method, &initial, &model, 1.0, &mut || {
            reverse_polls += 1;
            false
        })
        .unwrap();
        assert_eq!(gradient.q, ordinary.q);
        assert_eq!(gradient.v, ordinary.v);
        assert_eq!(gradient.a, ordinary.a);
        for stop in 1..=reverse_polls {
            let mut count = 0;
            let result = pull(&method, &initial, &model, 1.0, &mut || {
                count += 1;
                count == stop
            });
            assert!(
                matches!(
                    result,
                    Err(SecondOrderAdjointError::Step(TimeSolveError::Cancelled))
                ),
                "reverse boundary {stop}: {result:?}"
            );
            assert_eq!(initial, original);
        }
        let retry = pull(&method, &initial, &model, 1.0, &mut || false).unwrap();
        assert_eq!(format!("{gradient:?}"), format!("{retry:?}"));
    }

    // G0/G4: missing mass/damping/tangent/parameter derivatives, exhausted
    // solves, invalid seeds and resource admission cannot publish gradients.
    #[test]
    #[allow(clippy::too_many_lines)]
    fn structural_adjoint_refusals_and_zero_seed() {
        let method = OperatorGeneralizedAlpha::new(1, 0.15, 0.35, ImplicitSolveConfig::default());
        let initial = SecondOrderState::new(0.0, &[0.7], &[0.2], &[-0.3]);
        let budget = method.adjoint_workspace_components(1, config()).unwrap();
        let result = method.step_vjp(
            &initial,
            &Oscillator(0),
            &[0.1],
            (&[0.2], &[0.3], &[-0.15]),
            &Preconditioner(1.0),
            config(),
            budget - 1,
            &mut || false,
        );
        assert!(
            matches!(result, Err(SecondOrderAdjointError::WorkspaceLimit { required, .. }) if required == budget)
        );
        assert!(
            method
                .adjoint_workspace_components(usize::MAX, config())
                .is_none()
        );
        for mode in 1..=5 {
            let result = pull(&method, &initial, &Oscillator(mode), 1.0, &mut || false);
            if mode == 1 || mode == 4 {
                assert!(matches!(
                    result,
                    Err(SecondOrderAdjointError::Derivative(_))
                ));
            } else {
                assert!(matches!(
                    result,
                    Err(SecondOrderAdjointError::NonFiniteDerivative)
                ));
            }
        }
        let result = pull(&method, &initial, &Oscillator(0), 0.0, &mut || false);
        assert!(
            matches!(result, Err(SecondOrderAdjointError::NotConverged(report)) if !report.converged)
        );
        let limited = OperatorGeneralizedAlpha::new(
            1,
            0.15,
            0.35,
            ImplicitSolveConfig {
                max_newton_iterations: 1,
                ..ImplicitSolveConfig::default()
            },
        );
        let result = pull(&limited, &initial, &Oscillator(0), 1.0, &mut || false);
        assert!(matches!(
            result,
            Err(SecondOrderAdjointError::Step(TimeSolveError::NotConverged(
                _
            )))
        ));
        let result = method.step_vjp(
            &initial,
            &Oscillator(0),
            &[0.1],
            (&[0.2], &[0.3], &[f64::NAN]),
            &Preconditioner(1.0),
            config(),
            budget,
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(SecondOrderAdjointError::InvalidInput(_))
        ));
        let result = method.step_vjp(
            &initial,
            &Oscillator(0),
            &[0.1],
            (&[0.2], &[0.3], &[-0.15]),
            &Preconditioner(1.0),
            SecondOrderAdjointConfig {
                max_cycles: 0,
                ..config()
            },
            budget,
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(SecondOrderAdjointError::InvalidInput(_))
        ));
        let zero = method
            .step_vjp(
                &initial,
                &Oscillator(0),
                &[0.1],
                (&[0.0], &[0.0], &[0.0]),
                &Preconditioner(1.0),
                config(),
                budget,
                &mut || false,
            )
            .unwrap();
        assert_eq!(zero.initial_q, [0.0]);
        assert_eq!(zero.initial_v, [0.0]);
        assert_eq!(zero.initial_a, [0.0]);
        assert_eq!(zero.parameters, [0.0]);
        assert_eq!(zero.forcing, [0.0]);
        assert!(zero.adjoint.converged);
    }
}
