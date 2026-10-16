//! Consistent initial rates/accelerations and their matrix-free pullbacks.
//!
//! Generalized-alpha carries a rate or acceleration as an independent state
//! variable. When that variable comes from the initial physical residual, its
//! dependence on displacement/state, velocity, forcing and model parameters
//! must enter an inverse problem's chain rule. These routines solve that mass
//! system and its transpose with the existing bounded FGMRES implementation.
//! A nonsingular mass model and matched pure derivatives are caller obligations;
//! solver residuals do not establish uniqueness or physical model validity.

use super::adjoint::FirstOrderVjp;
use super::second_order_adjoint::SecondOrderVjp;
use super::{FirstOrderProblem, SecondOrderProblem};
use fs_solver::{FgmresState, FlexiblePreconditioner, LinearOp, SolveReport, StallDiagnosis};
use std::cell::RefCell;

/// Bounded controls used independently by the primal and transposed mass solves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InitialSolveConfig {
    /// Positive FGMRES restart length.
    pub restart: usize,
    /// Maximum restart cycles per mass solve.
    pub max_cycles: usize,
    /// Positive finite true-relative-residual tolerance.
    pub tolerance: f64,
}
impl Default for InitialSolveConfig {
    fn default() -> Self {
        Self {
            restart: 24,
            max_cycles: 16,
            tolerance: 1e-11,
        }
    }
}

/// A refused initialization never mutates inputs or returns a partial gradient.
#[derive(Debug, Clone)]
pub enum InitialSolveError {
    InvalidInput(&'static str),
    WorkspaceLimit { required: usize, limit: usize },
    NotConverged(SolveReport),
    Derivative(String),
    NonFiniteOutput,
    NonFiniteAccumulation,
    Cancelled,
}
impl std::fmt::Display for InitialSolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "consistent implicit initialization failed: {self:?}")
    }
}
impl std::error::Error for InitialSolveError {}

/// Consistent rate/acceleration and its recomputed mass-system residual report.
#[derive(Debug, Clone)]
pub struct InitialSolveResult {
    pub value: Vec<f64>,
    pub report: SolveReport,
}

/// Pullback of the consistent first-order rate with state and forcing independent.
#[derive(Debug, Clone)]
pub struct FirstOrderInitialGradient {
    pub rate: Vec<f64>,
    pub initial: Vec<f64>,
    pub parameters: Vec<f64>,
    pub forcing: Vec<f64>,
    pub primal: SolveReport,
    pub adjoint: SolveReport,
}

/// Pullback of the consistent structural acceleration with q/v/forcing independent.
#[derive(Debug, Clone)]
pub struct SecondOrderInitialGradient {
    pub acceleration: Vec<f64>,
    pub initial_q: Vec<f64>,
    pub initial_v: Vec<f64>,
    pub parameters: Vec<f64>,
    pub forcing: Vec<f64>,
    pub primal: SolveReport,
    pub adjoint: SolveReport,
}

/// Conservative live scalar-storage ceiling for either initialization pullback.
/// Includes Krylov basis/Hessenberg storage, histories and returned vectors;
/// excludes caller inputs, callback-owned storage and allocator metadata.
/// `None` indicates checked dimension/byte-capacity arithmetic overflow.
#[must_use]
pub fn initial_workspace_components(
    n: usize,
    p: usize,
    config: InitialSolveConfig,
) -> Option<usize> {
    let m = config.restart;
    let required = m
        .checked_mul(2)?
        .checked_add(24)?
        .checked_mul(n)?
        .checked_add(m.checked_add(1)?.checked_mul(m)?)?
        .checked_add(m.checked_mul(8)?)?
        .checked_add(config.max_cycles.checked_mul(8)?)?
        .checked_add(p.checked_mul(2)?)?
        .checked_add(32)?;
    (required <= isize::MAX as usize / std::mem::size_of::<f64>()).then_some(required)
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|x| x.is_finite())
}
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), InitialSolveError> {
    if cancelled() {
        Err(InitialSolveError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(
    n: usize,
    p: usize,
    vectors: &[&[f64]],
    config: InitialSolveConfig,
    cap: usize,
) -> Result<(), InitialSolveError> {
    if n == 0
        || config.restart == 0
        || config.max_cycles == 0
        || !config.tolerance.is_finite()
        || config.tolerance <= 0.0
    {
        return Err(InitialSolveError::InvalidInput(
            "nonempty state and positive finite solver controls required",
        ));
    }
    let required = initial_workspace_components(n, p, config).ok_or(
        InitialSolveError::InvalidInput("workspace dimension overflow"),
    )?;
    if required > cap {
        return Err(InitialSolveError::WorkspaceLimit {
            required,
            limit: cap,
        });
    }
    if vectors.iter().any(|v| v.len() != n || !finite(v)) {
        return Err(InitialSolveError::InvalidInput(
            "finite dimension-matched vectors required",
        ));
    }
    Ok(())
}

struct MassAction<'a, F> {
    n: usize,
    action: &'a F,
    error: RefCell<Option<InitialSolveError>>,
}
impl<F: Fn(&[f64], &mut [f64]) -> Result<(), String>> LinearOp for MassAction<'_, F> {
    fn n(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if self.error.borrow().is_some() {
            return;
        }
        if let Err(error) = (self.action)(x, out) {
            *self.error.borrow_mut() = Some(InitialSolveError::Derivative(error));
            out.fill(f64::NAN);
        } else if !finite(out) {
            *self.error.borrow_mut() = Some(InitialSolveError::NonFiniteOutput);
        }
    }
}

fn solve<F, P, C>(
    rhs: &[f64],
    action: F,
    preconditioner: &P,
    config: InitialSolveConfig,
    cancelled: &mut C,
) -> Result<InitialSolveResult, InitialSolveError>
where
    F: Fn(&[f64], &mut [f64]) -> Result<(), String>,
    P: FlexiblePreconditioner,
    C: FnMut() -> bool,
{
    if !finite(rhs) {
        return Err(InitialSolveError::NonFiniteAccumulation);
    }
    let op = MassAction {
        n: rhs.len(),
        action: &action,
        error: RefCell::new(None),
    };
    let mut state = FgmresState::new(rhs, config.restart);
    for cycle in 0..config.max_cycles {
        poll(cancelled)?;
        let report = state.run(&op, preconditioner, rhs, config.tolerance, 1);
        poll(cancelled)?;
        if let Some(error) = op.error.borrow_mut().take() {
            return Err(error);
        }
        if report.converged {
            if !finite(&state.x) {
                return Err(InitialSolveError::NonFiniteAccumulation);
            }
            return Ok(InitialSolveResult {
                value: state.x,
                report,
            });
        }
        if report.diagnosis == Some(StallDiagnosis::Breakdown) || cycle + 1 == config.max_cycles {
            return Err(InitialSolveError::NotConverged(report));
        }
    }
    unreachable!("admission requires a positive cycle budget")
}

/// Solve `M rate = forcing - r(time,initial)` before first-order stepping.
/// The preconditioner approximates M. Cancellation is polled around callbacks
/// and between restart cycles; one callback or cycle must finish before polling.
#[allow(clippy::too_many_arguments)]
pub fn first_order_rate<M, P, C>(
    model: &M,
    time: f64,
    initial: &[f64],
    forcing: &[f64],
    preconditioner: &P,
    config: InitialSolveConfig,
    max_workspace_components: usize,
    cancelled: &mut C,
) -> Result<InitialSolveResult, InitialSolveError>
where
    M: FirstOrderProblem + ?Sized,
    P: FlexiblePreconditioner,
    C: FnMut() -> bool,
{
    poll(cancelled)?;
    admit(
        model.dimension(),
        0,
        &[initial, forcing],
        config,
        max_workspace_components,
    )?;
    if !time.is_finite() {
        return Err(InitialSolveError::InvalidInput(
            "finite initial time required",
        ));
    }
    let mut rhs = vec![f64::NAN; initial.len()];
    model.internal_force(time, initial, &mut rhs);
    poll(cancelled)?;
    if !finite(&rhs) {
        return Err(InitialSolveError::NonFiniteOutput);
    }
    for (value, load) in rhs.iter_mut().zip(forcing) {
        *value = *load - *value;
    }
    solve(
        &rhs,
        |x, out| {
            model.mass_apply(x, out);
            Ok(())
        },
        preconditioner,
        config,
        cancelled,
    )
}

/// Solve `M acceleration = forcing - C velocity - r(displacement)` before
/// structural stepping. The force is supplied at the initial physical time,
/// rather than at the generalized-alpha step's intermediate load time.
#[allow(clippy::too_many_arguments)]
pub fn second_order_acceleration<M, P, C>(
    model: &M,
    displacement: &[f64],
    velocity: &[f64],
    forcing: &[f64],
    preconditioner: &P,
    config: InitialSolveConfig,
    max_workspace_components: usize,
    cancelled: &mut C,
) -> Result<InitialSolveResult, InitialSolveError>
where
    M: SecondOrderProblem + ?Sized,
    P: FlexiblePreconditioner,
    C: FnMut() -> bool,
{
    poll(cancelled)?;
    admit(
        model.dimension(),
        0,
        &[displacement, velocity, forcing],
        config,
        max_workspace_components,
    )?;
    let mut rhs = vec![f64::NAN; displacement.len()];
    let mut damping = vec![f64::NAN; displacement.len()];
    model.internal_force(displacement, &mut rhs);
    poll(cancelled)?;
    model.damping_apply(velocity, &mut damping);
    poll(cancelled)?;
    if !finite(&rhs) || !finite(&damping) {
        return Err(InitialSolveError::NonFiniteOutput);
    }
    for ((value, load), loss) in rhs.iter_mut().zip(forcing).zip(damping) {
        *value = *load - loss - *value;
    }
    solve(
        &rhs,
        |x, out| {
            model.mass_apply(x, out);
            Ok(())
        },
        preconditioner,
        config,
        cancelled,
    )
}

/// Pull back a seed on the consistent initial rate. Time is fixed. Add the
/// returned initial-state cotangent to any direct state seed; chain forcing
/// and state parameterizations separately. Model parameter partials include
/// both M(p) and r(t,u,p). The adjoint preconditioner must suit M^T.
#[allow(clippy::too_many_arguments)]
pub fn first_order_rate_vjp<M, P, Q, C>(
    model: &M,
    time: f64,
    initial: &[f64],
    forcing: &[f64],
    seed: &[f64],
    primal_preconditioner: &P,
    adjoint_preconditioner: &Q,
    config: InitialSolveConfig,
    max_workspace_components: usize,
    cancelled: &mut C,
) -> Result<FirstOrderInitialGradient, InitialSolveError>
where
    M: FirstOrderVjp + ?Sized,
    P: FlexiblePreconditioner,
    Q: FlexiblePreconditioner,
    C: FnMut() -> bool,
{
    poll(cancelled)?;
    let p = model.parameter_count();
    admit(
        model.dimension(),
        p,
        &[initial, forcing, seed],
        config,
        max_workspace_components,
    )?;
    let primal = first_order_rate(
        model,
        time,
        initial,
        forcing,
        primal_preconditioner,
        config,
        max_workspace_components,
        cancelled,
    )?;
    let adjoint = solve(
        seed,
        |x, out| model.mass_transpose_apply(x, out),
        adjoint_preconditioner,
        config,
        cancelled,
    )?;
    let mut initial_bar = vec![f64::NAN; initial.len()];
    let mut parameters = vec![f64::NAN; p];
    poll(cancelled)?;
    let result = model.tangent_transpose_apply(time, initial, &adjoint.value, &mut initial_bar);
    poll(cancelled)?;
    result.map_err(InitialSolveError::Derivative)?;
    let result = model.residual_parameter_vjp(
        time,
        initial,
        &primal.value,
        &adjoint.value,
        &mut parameters,
    );
    poll(cancelled)?;
    result.map_err(InitialSolveError::Derivative)?;
    if !finite(&initial_bar) || !finite(&parameters) {
        return Err(InitialSolveError::NonFiniteOutput);
    }
    for value in initial_bar.iter_mut().chain(&mut parameters) {
        *value = -*value;
    }
    poll(cancelled)?;
    Ok(FirstOrderInitialGradient {
        rate: primal.value,
        initial: initial_bar,
        parameters,
        forcing: adjoint.value,
        primal: primal.report,
        adjoint: adjoint.report,
    })
}

/// Pull back a seed on consistent initial acceleration. The returned q/v
/// cotangents include the internal-force and damping terms. Model parameter
/// derivatives include mass, damping and internal force; forcing receives its
/// own cotangent for the caller's load-parameter chain rule. All inputs remain
/// unchanged on failure. The adjoint preconditioner must approximate M^T.
#[allow(clippy::too_many_arguments)]
pub fn second_order_acceleration_vjp<M, P, Q, C>(
    model: &M,
    displacement: &[f64],
    velocity: &[f64],
    forcing: &[f64],
    seed: &[f64],
    primal_preconditioner: &P,
    adjoint_preconditioner: &Q,
    config: InitialSolveConfig,
    max_workspace_components: usize,
    cancelled: &mut C,
) -> Result<SecondOrderInitialGradient, InitialSolveError>
where
    M: SecondOrderVjp + ?Sized,
    P: FlexiblePreconditioner,
    Q: FlexiblePreconditioner,
    C: FnMut() -> bool,
{
    poll(cancelled)?;
    let p = model.parameter_count();
    admit(
        model.dimension(),
        p,
        &[displacement, velocity, forcing, seed],
        config,
        max_workspace_components,
    )?;
    let primal = second_order_acceleration(
        model,
        displacement,
        velocity,
        forcing,
        primal_preconditioner,
        config,
        max_workspace_components,
        cancelled,
    )?;
    let adjoint = solve(
        seed,
        |x, out| model.mass_transpose_apply(x, out),
        adjoint_preconditioner,
        config,
        cancelled,
    )?;
    let mut initial_q = vec![f64::NAN; displacement.len()];
    let mut initial_v = vec![f64::NAN; displacement.len()];
    let mut parameters = vec![f64::NAN; p];
    poll(cancelled)?;
    let result = model.tangent_transpose_apply(displacement, &adjoint.value, &mut initial_q);
    poll(cancelled)?;
    result.map_err(InitialSolveError::Derivative)?;
    let result = model.damping_transpose_apply(&adjoint.value, &mut initial_v);
    poll(cancelled)?;
    result.map_err(InitialSolveError::Derivative)?;
    let result = model.residual_parameter_vjp(
        displacement,
        velocity,
        &primal.value,
        &adjoint.value,
        &mut parameters,
    );
    poll(cancelled)?;
    result.map_err(InitialSolveError::Derivative)?;
    if !finite(&initial_q) || !finite(&initial_v) || !finite(&parameters) {
        return Err(InitialSolveError::NonFiniteOutput);
    }
    for value in initial_q
        .iter_mut()
        .chain(&mut initial_v)
        .chain(&mut parameters)
    {
        *value = -*value;
    }
    poll(cancelled)?;
    Ok(SecondOrderInitialGradient {
        acceleration: primal.value,
        initial_q,
        initial_v,
        parameters,
        forcing: adjoint.value,
        primal: primal.report,
        adjoint: adjoint.report,
    })
}
