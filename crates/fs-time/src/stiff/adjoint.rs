//! Discrete ARS(2,2,2) pullbacks for autonomous `u' = L(p) u + N(u,p)`.
//!
//! Differentiate each converged stage equation by the implicit function
//! theorem, never through FGMRES iterations. The two transposed solves use
//! the same shifted operator as the primal. Accuracy is limited by primal
//! and adjoint residuals and the supplied derivative actions. Time, step size,
//! solver policy and initial-condition parameterization are held fixed.

use super::{ImexSolveError, ImexStage, OperatorImex2, ShiftedLinearOp, imex_poll};
use fs_solver::{FlexiblePreconditioner, LinearOp, SolveReport};

pub mod trajectory;

/// A fixed parameter point with matched primal and transposed derivatives.
/// Implementations must be pure throughout one call, overwrite every output,
/// and implement the actual `L^T` action for a nonsymmetric `LinearOp`.
/// Callback work must be bounded; cancellation is polled around callbacks.
pub trait ImexVjp: LinearOp {
    fn parameter_count(&self) -> usize;
    fn nonlinear(&self, state: &[f64], output: &mut [f64]);
    /// `(N_u^T seed, N_p^T seed)` at fixed state and parameter point.
    fn nonlinear_vjp(
        &self,
        state: &[f64],
        seed: &[f64],
        state_bar: &mut [f64],
        parameter_bar: &mut [f64],
    ) -> Result<(), String>;
    /// `(d[L(p) state]/dp)^T seed`, holding `state` fixed.
    /// An explicitly parameter-independent L must write zeros here.
    fn linear_parameter_vjp(
        &self,
        state: &[f64],
        seed: &[f64],
        parameter_bar: &mut [f64],
    ) -> Result<(), String>;
}

#[derive(Debug, Clone)]
pub enum ImexAdjointError {
    InvalidInput(&'static str),
    WorkspaceLimit { required: usize, limit: usize },
    Step(ImexSolveError),
    Derivative(String),
    NonFiniteDerivative,
    NonFiniteAccumulation,
}

impl From<ImexSolveError> for ImexAdjointError {
    fn from(error: ImexSolveError) -> Self {
        Self::Step(error)
    }
}
impl std::fmt::Display for ImexAdjointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "IMEX discrete adjoint failed: {self:?}")
    }
}
impl std::error::Error for ImexAdjointError {}

/// Endpoint and pullback of a caller-supplied endpoint cotangent.
#[derive(Debug, Clone)]
pub struct ImexStepGradient {
    pub value: Vec<f64>,
    pub initial: Vec<f64>,
    pub parameters: Vec<f64>,
    /// Stage-one and stage-two primal true-residual reports.
    pub primal: [SolveReport; 2],
    /// Stage-one and stage-two transposed true-residual reports.
    pub adjoint: [SolveReport; 2],
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

fn accumulate(total: &mut [f64], contribution: &[f64]) -> Result<(), ImexAdjointError> {
    if !finite(contribution) {
        return Err(ImexAdjointError::NonFiniteDerivative);
    }
    for (sum, value) in total.iter_mut().zip(contribution) {
        *sum += value;
    }
    if !finite(total) {
        return Err(ImexAdjointError::NonFiniteAccumulation);
    }
    Ok(())
}

struct Transpose<'a, L>(&'a L);
impl<L: LinearOp> LinearOp for Transpose<'_, L> {
    fn n(&self) -> usize {
        self.0.n()
    }
    fn apply(&self, x: &[f64], out: &mut [f64]) {
        self.0.apply_transpose(x, out);
    }
    fn apply_transpose(&self, x: &[f64], out: &mut [f64]) {
        self.0.apply(x, out);
    }
}

impl OperatorImex2 {
    /// Conservative live scalar-storage bound for `step_vjp`, including the
    /// primal stages, FGMRES basis/Hessenberg storage, residual histories and
    /// gradient outputs. Excludes allocator metadata and callback-owned memory.
    /// Uses the configured restart length, even when it exceeds the dimension.
    pub fn adjoint_workspace_components(&self, parameters: usize) -> Option<usize> {
        let m = self.solve.restart;
        let vectors = m.checked_mul(2)?.checked_add(24)?.checked_mul(self.n)?;
        let krylov = m
            .checked_add(1)?
            .checked_mul(m)?
            .checked_add(m.checked_mul(8)?)?;
        vectors
            .checked_add(krylov)?
            .checked_add(self.solve.max_cycles.checked_mul(8)?)?
            .checked_add(parameters.checked_mul(2)?)?
            .checked_add(32)
    }

    /// Pull back one production IMEX step. Both preconditioners are explicit:
    /// the adjoint preconditioner must suit `I - gamma*h*L^T`; a nonsymmetric
    /// primal preconditioner is not silently treated as its own transpose.
    ///
    /// Reuses the forward stage implementation. Four bounded linear solves
    /// (two primal, two transposed) must converge before a gradient is returned.
    /// Cancellation or callback/solve failure returns no partial gradient and
    /// does not mutate the input. Polling granularity is one restart cycle or
    /// one callback, as in `step_controlled`. Chain these pullbacks in reverse
    /// time for terminal objectives and sum their parameter contributions;
    /// direct objective/initial-condition derivatives belong to the caller.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_vjp<M, P, Q, Cancel>(
        &self,
        initial: &[f64],
        model: &M,
        primal_preconditioner: &P,
        adjoint_preconditioner: &Q,
        terminal: &[f64],
        max_workspace_components: usize,
        cancelled: &mut Cancel,
    ) -> Result<ImexStepGradient, ImexAdjointError>
    where
        M: ImexVjp,
        P: FlexiblePreconditioner,
        Q: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        imex_poll(cancelled)?;
        if model.n() != self.n
            || initial.len() != self.n
            || terminal.len() != self.n
            || !finite(initial)
            || !finite(terminal)
        {
            return Err(ImexAdjointError::InvalidInput(
                "finite dimension-matched state and seed required",
            ));
        }
        let p = model.parameter_count();
        let required =
            self.adjoint_workspace_components(p)
                .ok_or(ImexAdjointError::InvalidInput(
                    "workspace dimension overflow",
                ))?;
        if required > max_workspace_components {
            return Err(ImexAdjointError::WorkspaceLimit {
                required,
                limit: max_workspace_components,
            });
        }
        let stages = self.stages(
            initial,
            model,
            primal_preconditioner,
            &|state, out| model.nonlinear(state, out),
            cancelled,
        )?;
        self.reverse_stages(
            initial,
            model,
            adjoint_preconditioner,
            terminal,
            p,
            stages,
            cancelled,
        )
    }

    // Reuse the already replayed and checked primal stages in trajectory
    // sweeps; no second forward solve is needed for a leaf pullback.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn reverse_stages<M, Q, Cancel>(
        &self,
        initial: &[f64],
        model: &M,
        adjoint_preconditioner: &Q,
        terminal: &[f64],
        p: usize,
        stages: super::ImexStages,
        cancelled: &mut Cancel,
    ) -> Result<ImexStepGradient, ImexAdjointError>
    where
        M: ImexVjp,
        Q: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        let shifted = ShiftedLinearOp {
            linear: model,
            shift: -self.gamma * self.h,
        };
        let transposed = Transpose(&shifted);
        let (lambda_two, report_two) = self.solve_stage(
            &transposed,
            adjoint_preconditioner,
            terminal,
            ImexStage::Two,
            cancelled,
        )?;
        let delta = 1.0 - 1.0 / (2.0 * self.gamma);
        let mut parameters = vec![0.0; p];
        let mut parameter_scratch = vec![f64::NAN; p];
        let mut one_bar = vec![f64::NAN; self.n];
        let mut seed = lambda_two
            .iter()
            .map(|v| self.h * (1.0 - delta) * v)
            .collect::<Vec<_>>();
        if !finite(&seed) {
            return Err(ImexAdjointError::NonFiniteAccumulation);
        }
        imex_poll(cancelled)?;
        model
            .nonlinear_vjp(&stages.one, &seed, &mut one_bar, &mut parameter_scratch)
            .map_err(ImexAdjointError::Derivative)?;
        imex_poll(cancelled)?;
        if !finite(&one_bar) {
            return Err(ImexAdjointError::NonFiniteDerivative);
        }
        accumulate(&mut parameters, &parameter_scratch)?;
        let mut state_scratch = vec![f64::NAN; self.n];
        model.apply_transpose(&lambda_two, &mut state_scratch);
        imex_poll(cancelled)?;
        if !finite(&state_scratch) {
            return Err(ImexAdjointError::NonFiniteDerivative);
        }
        for i in 0..self.n {
            one_bar[i] += self.h * (1.0 - self.gamma) * state_scratch[i];
        }
        if !finite(&one_bar) {
            return Err(ImexAdjointError::NonFiniteAccumulation);
        }
        let (lambda_one, report_one) = self.solve_stage(
            &transposed,
            adjoint_preconditioner,
            &one_bar,
            ImexStage::One,
            cancelled,
        )?;
        // The original N(u) is shared by both stages, so its two seeds add.
        for i in 0..self.n {
            seed[i] = self.h * (delta * lambda_two[i] + self.gamma * lambda_one[i]);
        }
        if !finite(&seed) {
            return Err(ImexAdjointError::NonFiniteAccumulation);
        }
        state_scratch.fill(f64::NAN);
        parameter_scratch.fill(f64::NAN);
        imex_poll(cancelled)?;
        model
            .nonlinear_vjp(initial, &seed, &mut state_scratch, &mut parameter_scratch)
            .map_err(ImexAdjointError::Derivative)?;
        imex_poll(cancelled)?;
        if !finite(&state_scratch) {
            return Err(ImexAdjointError::NonFiniteDerivative);
        }
        accumulate(&mut parameters, &parameter_scratch)?;
        let mut initial_bar = state_scratch;
        for i in 0..self.n {
            initial_bar[i] += lambda_two[i] + lambda_one[i];
        }
        if !finite(&initial_bar) {
            return Err(ImexAdjointError::NonFiniteAccumulation);
        }
        // Parameter dependence of both implicit matrices and the explicit L*u1
        // term is essential: omitting either changes the discretized gradient.
        for (state, lambda, weight) in [
            (&stages.next, &lambda_two, self.gamma),
            (&stages.one, &lambda_two, 1.0 - self.gamma),
            (&stages.one, &lambda_one, self.gamma),
        ] {
            for i in 0..self.n {
                seed[i] = self.h * weight * lambda[i];
            }
            if !finite(&seed) {
                return Err(ImexAdjointError::NonFiniteAccumulation);
            }
            parameter_scratch.fill(f64::NAN);
            imex_poll(cancelled)?;
            model
                .linear_parameter_vjp(state, &seed, &mut parameter_scratch)
                .map_err(ImexAdjointError::Derivative)?;
            imex_poll(cancelled)?;
            accumulate(&mut parameters, &parameter_scratch)?;
        }
        imex_poll(cancelled)?;
        Ok(ImexStepGradient {
            value: stages.next,
            initial: initial_bar,
            parameters,
            primal: [stages.report_one, stages.report_two],
            adjoint: [report_one, report_two],
        })
    }
}
