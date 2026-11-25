//! Fallible matrix-free Newton--Krylov solver that carries `Result` end-to-end.
//!
//! Unlike [`crate::NonlinearProblem`], whose callbacks are infallible and forced
//! callers to smuggle component refusals through `NaN` and side latches,
//! [`FallibleNonlinearProblem`] propagates the owner's concrete error `E`
//! directly through residual evaluation, Jacobian-vector actions, and
//! preconditioning without ever encoding a refusal as `NaN`.
//!
//! Inner products and Givens rotations use the crate's deterministic reductions
//! (`crate::dot`, `crate::norm2`, `fs_math::det::hypot`), and cooperative
//! cancellation is polled through [`fs_exec::Cx`] at every residual, Jacobian,
//! preconditioner, Arnoldi-column, and line-search boundary.

#![forbid(unsafe_code)]

use core::fmt;

use fs_exec::Cx;

use crate::krylov::StallDiagnosis;
use crate::nonlinear::{
    Globalization, GlobalizationDecision, NewtonIteration, NewtonKrylovConfig, NewtonReport,
    NewtonStallDiagnosis,
};
use crate::{dot, norm2};

/// Fallible nonlinear residual and Jacobian-vector action at an exact state.
pub trait FallibleNonlinearProblem {
    /// Owner refusal type carried verbatim without `NaN` encoding.
    type Error;

    /// State and residual vector dimension.
    fn dimension(&self) -> usize;

    /// Overwrite `residual` with `F(x)`.
    ///
    /// # Errors
    /// Propagates any owner or domain refusal directly.
    fn residual(&self, x: &[f64], residual: &mut [f64], cx: &Cx<'_>) -> Result<(), Self::Error>;

    /// Overwrite `output` with `J(x) * direction`.
    ///
    /// # Errors
    /// Propagates any owner or tangent refusal directly.
    fn jacobian_apply(
        &self,
        x: &[f64],
        direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;

    /// Optional right preconditioner action at `(x, outer_iteration, inner_iteration)`.
    /// Defaults to identity.
    ///
    /// # Errors
    /// Propagates any preconditioner refusal directly.
    fn preconditioner_apply(
        &self,
        _x: &[f64],
        _outer_iteration: usize,
        _inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Self::Error> {
        output.copy_from_slice(residual);
        Ok(())
    }
}

/// Failure returned by [`solve_fallible_newton_krylov`].
#[derive(Debug)]
pub enum FallibleNewtonError<E> {
    /// Invalid solver configuration.
    InvalidConfig(&'static str),
    /// State length did not match `problem.dimension()`.
    Dimension {
        /// Dimension declared by the problem.
        problem: usize,
        /// Length of the supplied state vector.
        state: usize,
    },
    /// Cooperative cancellation was requested through `Cx`.
    Cancelled,
    /// The problem refused a residual, Jacobian, or preconditioner evaluation.
    Problem(E),
    /// A state coordinate was non-finite.
    NonFiniteState {
        /// Outer Newton iteration index.
        iteration: usize,
        /// First non-finite coordinate.
        index: usize,
        /// IEEE-754 bit pattern.
        bits: u64,
    },
    /// A problem callback returned a non-finite residual coordinate despite `Ok(())`.
    NonFiniteResidual {
        /// Outer Newton iteration index.
        iteration: usize,
        /// First non-finite coordinate.
        index: usize,
        /// IEEE-754 bit pattern.
        bits: u64,
    },
    /// A Jacobian or preconditioner callback returned non-finite output despite `Ok(())`.
    NonFiniteOperator {
        /// Outer Newton iteration index.
        iteration: usize,
        /// Stage where non-finiteness was observed.
        stage: &'static str,
    },
}

impl<E: fmt::Display> fmt::Display for FallibleNewtonError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig(msg) => write!(f, "invalid Newton-Krylov configuration: {msg}"),
            Self::Dimension { problem, state } => {
                write!(
                    f,
                    "dimension mismatch: problem has {problem}, state has {state}"
                )
            }
            Self::Cancelled => write!(f, "fallible Newton-Krylov cancelled"),
            Self::Problem(err) => write!(f, "nonlinear problem refused: {err}"),
            Self::NonFiniteState {
                iteration,
                index,
                bits,
            } => write!(
                f,
                "non-finite state at iteration {iteration}, index {index} (bits {bits:#018x})"
            ),
            Self::NonFiniteResidual {
                iteration,
                index,
                bits,
            } => write!(
                f,
                "non-finite residual at iteration {iteration}, index {index} (bits {bits:#018x})"
            ),
            Self::NonFiniteOperator { iteration, stage } => {
                write!(f, "non-finite {stage} output at iteration {iteration}")
            }
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for FallibleNewtonError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Problem(err) => Some(err),
            _ => None,
        }
    }
}

/// Completed fallible Newton--Krylov solve (either converged or diagnosed stall).
#[derive(Debug, Clone, PartialEq)]
pub struct FallibleNewtonSolution {
    /// Final state iterate.
    pub x: Vec<f64>,
    /// Residual vector at `x`.
    pub residual: Vec<f64>,
    /// Iteration history, final residual norm, and convergence/stall verdict.
    pub report: NewtonReport,
    /// Total residual evaluations performed.
    pub residual_evaluations: u64,
    /// Total Jacobian-vector evaluations performed.
    pub jacobian_evaluations: u64,
}

fn check_cancel<E>(cx: &Cx<'_>) -> Result<(), FallibleNewtonError<E>> {
    cx.checkpoint().map_err(|_| FallibleNewtonError::Cancelled)
}

#[allow(clippy::too_many_lines)]
fn validate_config<E>(config: NewtonKrylovConfig) -> Result<(), FallibleNewtonError<E>> {
    if !config.absolute_tolerance.is_finite() || config.absolute_tolerance < 0.0 {
        return Err(FallibleNewtonError::InvalidConfig(
            "absolute_tolerance must be finite and nonnegative",
        ));
    }
    if !config.relative_tolerance.is_finite()
        || config.relative_tolerance < 0.0
        || config.relative_tolerance > 1.0
    {
        return Err(FallibleNewtonError::InvalidConfig(
            "relative_tolerance must be in [0, 1]",
        ));
    }
    if config.linear_restart == 0 || config.max_linear_cycles == 0 {
        return Err(FallibleNewtonError::InvalidConfig(
            "linear_restart and max_linear_cycles must be positive",
        ));
    }
    if !config.forcing_minimum.is_finite()
        || !config.forcing_maximum.is_finite()
        || config.forcing_minimum <= 0.0
        || config.forcing_maximum < config.forcing_minimum
        || config.forcing_maximum >= 1.0
    {
        return Err(FallibleNewtonError::InvalidConfig(
            "forcing bounds must satisfy 0 < minimum <= maximum < 1",
        ));
    }
    if !config.forcing_gamma.is_finite()
        || config.forcing_gamma <= 0.0
        || config.forcing_gamma > 1.0
    {
        return Err(FallibleNewtonError::InvalidConfig(
            "forcing_gamma must be in (0, 1]",
        ));
    }
    if !config.forcing_exponent.is_finite()
        || config.forcing_exponent < 1.0
        || config.forcing_exponent > 2.0
    {
        return Err(FallibleNewtonError::InvalidConfig(
            "forcing_exponent must be in [1, 2]",
        ));
    }
    match config.globalization {
        Globalization::LineSearch(line) => {
            if !line.armijo.is_finite() || line.armijo <= 0.0 || line.armijo >= 1.0 {
                return Err(FallibleNewtonError::InvalidConfig(
                    "armijo must be in (0, 1)",
                ));
            }
            if !line.contraction.is_finite() || line.contraction <= 0.0 || line.contraction >= 1.0 {
                return Err(FallibleNewtonError::InvalidConfig(
                    "contraction must be in (0, 1)",
                ));
            }
            if !line.minimum_step.is_finite() || line.minimum_step <= 0.0 || line.minimum_step > 1.0
            {
                return Err(FallibleNewtonError::InvalidConfig(
                    "minimum_step must be in (0, 1]",
                ));
            }
        }
        Globalization::TrustRegion(trust) => {
            if !trust.initial_radius.is_finite()
                || !trust.minimum_radius.is_finite()
                || !trust.maximum_radius.is_finite()
                || trust.minimum_radius <= 0.0
                || trust.initial_radius < trust.minimum_radius
                || trust.maximum_radius < trust.initial_radius
            {
                return Err(FallibleNewtonError::InvalidConfig(
                    "trust radii must satisfy 0 < min <= initial <= max",
                ));
            }
            if !trust.acceptance_ratio.is_finite()
                || !trust.expansion_ratio.is_finite()
                || trust.acceptance_ratio <= 0.0
                || trust.expansion_ratio <= trust.acceptance_ratio
                || trust.expansion_ratio >= 1.0
            {
                return Err(FallibleNewtonError::InvalidConfig(
                    "trust ratios must satisfy 0 < acceptance < expansion < 1",
                ));
            }
            if !trust.shrink.is_finite()
                || trust.shrink <= 0.0
                || trust.shrink >= 1.0
                || !trust.expansion.is_finite()
                || trust.expansion <= 1.0
            {
                return Err(FallibleNewtonError::InvalidConfig(
                    "trust shrink must be in (0, 1) and expansion > 1",
                ));
            }
        }
    }
    Ok(())
}

fn eval_residual<P: FallibleNonlinearProblem>(
    problem: &P,
    x: &[f64],
    iteration: usize,
    cx: &Cx<'_>,
    res_evals: &mut u64,
) -> Result<(Vec<f64>, f64), FallibleNewtonError<P::Error>> {
    check_cancel(cx)?;
    for (index, &val) in x.iter().enumerate() {
        if !val.is_finite() {
            return Err(FallibleNewtonError::NonFiniteState {
                iteration,
                index,
                bits: val.to_bits(),
            });
        }
    }
    let mut residual = vec![0.0; problem.dimension()];
    *res_evals = res_evals.saturating_add(1);
    problem
        .residual(x, &mut residual, cx)
        .map_err(FallibleNewtonError::Problem)?;
    check_cancel(cx)?;
    for (index, &val) in residual.iter().enumerate() {
        if !val.is_finite() {
            return Err(FallibleNewtonError::NonFiniteResidual {
                iteration,
                index,
                bits: val.to_bits(),
            });
        }
    }
    let norm = norm2(&residual);
    if !norm.is_finite() {
        return Err(FallibleNewtonError::NonFiniteResidual {
            iteration,
            index: 0,
            bits: norm.to_bits(),
        });
    }
    Ok((residual, norm))
}

struct LinearOutcome {
    direction: Vec<f64>,
    iters: usize,
    rel_residual: f64,
    converged: bool,
    diagnosis: Option<StallDiagnosis>,
}

/// Deterministic right-preconditioned FGMRES directly over [`FallibleNonlinearProblem`].
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_linear_fgmres<P: FallibleNonlinearProblem>(
    problem: &P,
    x_point: &[f64],
    rhs: &[f64],
    outer_iteration: usize,
    restart: usize,
    max_cycles: usize,
    tolerance: f64,
    cx: &Cx<'_>,
    jac_evals: &mut u64,
) -> Result<LinearOutcome, FallibleNewtonError<P::Error>> {
    let n = problem.dimension();
    let mut x = vec![0.0; n];
    let bnorm = norm2(rhs);
    if !bnorm.is_finite() {
        return Err(FallibleNewtonError::NonFiniteOperator {
            iteration: outer_iteration,
            stage: "rhs norm",
        });
    }
    if bnorm == 0.0 {
        return Ok(LinearOutcome {
            direction: x,
            iters: 0,
            rel_residual: 0.0,
            converged: true,
            diagnosis: None,
        });
    }
    let mut applied = vec![0.0; n];
    let mut iters = 0usize;
    let mut rel = 1.0_f64;
    let mut broken = false;

    for _cycle in 0..max_cycles {
        check_cancel(cx)?;
        if iters == 0 {
            applied.fill(0.0);
        } else {
            applied.fill(0.0);
            *jac_evals = jac_evals.saturating_add(1);
            problem
                .jacobian_apply(x_point, &x, &mut applied, cx)
                .map_err(FallibleNewtonError::Problem)?;
            check_cancel(cx)?;
            if applied.iter().any(|v| !v.is_finite()) {
                return Err(FallibleNewtonError::NonFiniteOperator {
                    iteration: outer_iteration,
                    stage: "jacobian_apply",
                });
            }
        }
        let mut residual: Vec<f64> = rhs.iter().zip(&applied).map(|(r, a)| r - a).collect();
        let beta = norm2(&residual);
        rel = beta / bnorm;
        if rel <= tolerance {
            break;
        }
        if !beta.is_finite() || beta == 0.0 {
            broken = true;
            break;
        }
        for value in &mut residual {
            *value /= beta;
        }

        let m = restart.min(n.max(1));
        let mut basis = vec![residual];
        let mut preconditioned_basis: Vec<Vec<f64>> = Vec::with_capacity(m);
        let mut h = vec![0.0; (m + 1) * m];
        let mut cs = vec![0.0_f64; m];
        let mut sn = vec![0.0_f64; m];
        let mut g = vec![0.0; m + 1];
        g[0] = beta;
        let mut columns = 0usize;

        for column in 0..m {
            check_cancel(cx)?;
            let mut z = vec![0.0; n];
            problem
                .preconditioner_apply(x_point, outer_iteration, iters, &basis[column], &mut z, cx)
                .map_err(FallibleNewtonError::Problem)?;
            check_cancel(cx)?;
            if z.iter().any(|v| !v.is_finite()) {
                return Err(FallibleNewtonError::NonFiniteOperator {
                    iteration: outer_iteration,
                    stage: "preconditioner_apply",
                });
            }

            let mut w = vec![0.0; n];
            *jac_evals = jac_evals.saturating_add(1);
            problem
                .jacobian_apply(x_point, &z, &mut w, cx)
                .map_err(FallibleNewtonError::Problem)?;
            check_cancel(cx)?;
            if w.iter().any(|v| !v.is_finite()) {
                return Err(FallibleNewtonError::NonFiniteOperator {
                    iteration: outer_iteration,
                    stage: "jacobian_apply",
                });
            }

            // Twice-orthogonalized Modified Gram-Schmidt for numerical stability.
            for _reorth in 0..2 {
                for (row, vector) in basis.iter().enumerate() {
                    check_cancel(cx)?;
                    let coeff = dot(vector, &w);
                    h[row * m + column] += coeff;
                    for (val, b_val) in w.iter_mut().zip(vector) {
                        *val = coeff.mul_add(-b_val, *val);
                    }
                }
            }
            let next_norm = norm2(&w);
            h[(column + 1) * m + column] = next_norm;
            for row in 0..column {
                let upper = h[row * m + column];
                let lower = h[(row + 1) * m + column];
                h[row * m + column] = cs[row].mul_add(upper, sn[row] * lower);
                h[(row + 1) * m + column] = (-sn[row]).mul_add(upper, cs[row] * lower);
            }
            let diagonal = h[column * m + column];
            let denom = fs_math::det::hypot(diagonal, next_norm);
            if !denom.is_finite() || denom == 0.0 {
                broken = true;
                break;
            }
            cs[column] = diagonal / denom;
            sn[column] = next_norm / denom;
            h[column * m + column] = denom;
            h[(column + 1) * m + column] = 0.0;
            g[column + 1] = -sn[column] * g[column];
            g[column] *= cs[column];
            preconditioned_basis.push(z);
            columns = column + 1;
            iters += 1;
            if next_norm == 0.0 || g[column + 1].abs() / bnorm <= tolerance {
                break;
            }
            for val in &mut w {
                *val /= next_norm;
            }
            basis.push(w);
        }
        if broken || columns == 0 {
            break;
        }
        let mut coeffs = vec![0.0; columns];
        for row in (0..columns).rev() {
            let mut val = g[row];
            for col in (row + 1)..columns {
                val = h[row * m + col].mul_add(-coeffs[col], val);
            }
            let diag = h[row * m + row];
            if !diag.is_finite() || diag == 0.0 {
                broken = true;
                break;
            }
            coeffs[row] = val / diag;
        }
        if broken {
            break;
        }
        for (coeff, vec_z) in coeffs.iter().zip(&preconditioned_basis) {
            for (xi, zi) in x.iter_mut().zip(vec_z) {
                *xi = coeff.mul_add(*zi, *xi);
            }
        }
        applied.fill(0.0);
        *jac_evals = jac_evals.saturating_add(1);
        problem
            .jacobian_apply(x_point, &x, &mut applied, cx)
            .map_err(FallibleNewtonError::Problem)?;
        check_cancel(cx)?;
        if applied.iter().any(|v| !v.is_finite()) {
            return Err(FallibleNewtonError::NonFiniteOperator {
                iteration: outer_iteration,
                stage: "jacobian_apply",
            });
        }
        let true_res: Vec<f64> = rhs.iter().zip(&applied).map(|(r, a)| r - a).collect();
        rel = norm2(&true_res) / bnorm;
        if rel <= tolerance {
            break;
        }
    }

    let converged = !broken && rel.is_finite() && rel <= tolerance;
    let diagnosis = if converged {
        None
    } else if broken || !rel.is_finite() {
        Some(StallDiagnosis::Breakdown)
    } else {
        Some(StallDiagnosis::BudgetExhausted)
    };
    Ok(LinearOutcome {
        direction: x,
        iters,
        rel_residual: rel,
        converged,
        diagnosis,
    })
}

/// Solve `F(x) = 0` with matrix-free Newton--Krylov and direct typed error propagation.
///
/// # Errors
/// Returns [`FallibleNewtonError::Problem`] with the owner's exact refusal (never encoded
/// as `NaN`), [`FallibleNewtonError::Cancelled`] on cooperative `Cx` cancellation, or a
/// configuration/dimension/non-finiteness error.
#[allow(clippy::too_many_lines)]
pub fn solve_fallible_newton_krylov<P: FallibleNonlinearProblem>(
    problem: &P,
    mut x: Vec<f64>,
    config: NewtonKrylovConfig,
    max_iterations: usize,
    cx: &Cx<'_>,
) -> Result<FallibleNewtonSolution, FallibleNewtonError<P::Error>> {
    validate_config(config)?;
    let n = problem.dimension();
    if x.len() != n {
        return Err(FallibleNewtonError::Dimension {
            problem: n,
            state: x.len(),
        });
    }
    let mut residual_evaluations = 0_u64;
    let mut jacobian_evaluations = 0_u64;
    let (mut residual, mut residual_norm) =
        eval_residual(problem, &x, 0, cx, &mut residual_evaluations)?;
    let initial_norm = residual_norm.max(f64::MIN_POSITIVE);
    let target = config
        .absolute_tolerance
        .max(config.relative_tolerance * initial_norm);
    let mut previous_norm = residual_norm;
    let mut trust_radius = match config.globalization {
        Globalization::LineSearch(_) => 0.0,
        Globalization::TrustRegion(trust) => trust.initial_radius,
    };
    let mut history = Vec::new();
    let mut iterations = 0usize;
    let mut terminal: Option<NewtonStallDiagnosis> = None;

    for _ in 0..max_iterations {
        check_cancel(cx)?;
        if residual_norm <= target || terminal.is_some() {
            break;
        }
        let residual_before = residual_norm;
        let forcing = if iterations == 0 {
            config.forcing_maximum
        } else {
            let ratio = residual_before / previous_norm.max(f64::MIN_POSITIVE);
            (config.forcing_gamma * fs_math::det::pow(ratio, config.forcing_exponent))
                .clamp(config.forcing_minimum, config.forcing_maximum)
        };
        let rhs: Vec<f64> = residual.iter().map(|v| -*v).collect();
        let linear = solve_linear_fgmres(
            problem,
            &x,
            &rhs,
            iterations,
            config.linear_restart,
            config.max_linear_cycles,
            forcing,
            cx,
            &mut jacobian_evaluations,
        )?;
        if !linear.converged {
            terminal = Some(NewtonStallDiagnosis::LinearSolveFailed(
                linear.diagnosis.unwrap_or(StallDiagnosis::BudgetExhausted),
            ));
            break;
        }
        let direction_norm = norm2(&linear.direction);
        if !direction_norm.is_finite() {
            terminal = Some(NewtonStallDiagnosis::NonFinite);
            break;
        }

        let iteration = iterations;
        let accepted = match config.globalization {
            Globalization::LineSearch(line) => {
                let mut alpha = 1.0_f64;
                let mut last_alpha = alpha;
                let mut found = None;
                while alpha >= line.minimum_step {
                    check_cancel(cx)?;
                    last_alpha = alpha;
                    let trial_x: Vec<f64> = x
                        .iter()
                        .zip(&linear.direction)
                        .map(|(&xi, &di)| alpha.mul_add(di, xi))
                        .collect();
                    let (trial_res, trial_norm) = eval_residual(
                        problem,
                        &trial_x,
                        iteration + 1,
                        cx,
                        &mut residual_evaluations,
                    )?;
                    if trial_norm <= (1.0 - line.armijo * alpha) * residual_before
                        || trial_norm <= target
                    {
                        found = Some((
                            trial_x,
                            trial_res,
                            trial_norm,
                            alpha,
                            None,
                            GlobalizationDecision::LineSearchAccepted,
                        ));
                        break;
                    }
                    alpha *= line.contraction;
                }
                if found.is_none() {
                    history.push(NewtonIteration {
                        iteration,
                        residual_before,
                        residual_after: residual_before,
                        forcing,
                        linear_relative_residual: linear.rel_residual,
                        linear_iterations: linear.iters,
                        newton_step_norm: direction_norm,
                        step_length: last_alpha,
                        reduction_ratio: None,
                        decision: GlobalizationDecision::LineSearchRejected,
                    });
                    iterations += 1;
                    terminal = Some(NewtonStallDiagnosis::LineSearchRejected);
                }
                found
            }
            Globalization::TrustRegion(trust) => {
                let alpha = if direction_norm > trust_radius {
                    trust_radius / direction_norm
                } else {
                    1.0
                };
                let trial_x: Vec<f64> = x
                    .iter()
                    .zip(&linear.direction)
                    .map(|(&xi, &di)| alpha.mul_add(di, xi))
                    .collect();
                let (trial_res, trial_norm) = eval_residual(
                    problem,
                    &trial_x,
                    iteration + 1,
                    cx,
                    &mut residual_evaluations,
                )?;
                let mut j_step = vec![0.0; n];
                let step_vec: Vec<f64> = linear.direction.iter().map(|&d| alpha * d).collect();
                jacobian_evaluations = jacobian_evaluations.saturating_add(1);
                problem
                    .jacobian_apply(&x, &step_vec, &mut j_step, cx)
                    .map_err(FallibleNewtonError::Problem)?;
                check_cancel(cx)?;
                if j_step.iter().any(|v| !v.is_finite()) {
                    return Err(FallibleNewtonError::NonFiniteOperator {
                        iteration,
                        stage: "jacobian_apply",
                    });
                }
                let model: Vec<f64> = residual.iter().zip(&j_step).map(|(r, js)| r + js).collect();
                let predicted = residual_before - norm2(&model);
                let actual = residual_before - trial_norm;
                let ratio = (predicted > 0.0 && predicted.is_finite() && actual.is_finite())
                    .then_some(actual / predicted);
                if trial_norm <= target
                    || (actual > 0.0 && ratio.is_some_and(|r| r >= trust.acceptance_ratio))
                {
                    if ratio.is_some_and(|r| r >= trust.expansion_ratio) && alpha < 1.0 {
                        trust_radius = (trust_radius * trust.expansion).min(trust.maximum_radius);
                    }
                    Some((
                        trial_x,
                        trial_res,
                        trial_norm,
                        alpha,
                        ratio,
                        GlobalizationDecision::TrustRegionAccepted,
                    ))
                } else {
                    trust_radius *= trust.shrink;
                    history.push(NewtonIteration {
                        iteration,
                        residual_before,
                        residual_after: residual_before,
                        forcing,
                        linear_relative_residual: linear.rel_residual,
                        linear_iterations: linear.iters,
                        newton_step_norm: direction_norm,
                        step_length: alpha,
                        reduction_ratio: ratio,
                        decision: GlobalizationDecision::TrustRegionRejected,
                    });
                    iterations += 1;
                    if trust_radius < trust.minimum_radius {
                        terminal = Some(NewtonStallDiagnosis::TrustRegionRadiusExhausted);
                    }
                    None
                }
            }
        };

        if let Some((trial_x, trial_res, trial_norm, step_length, ratio, decision)) = accepted {
            previous_norm = residual_before;
            x = trial_x;
            residual = trial_res;
            residual_norm = trial_norm;
            history.push(NewtonIteration {
                iteration,
                residual_before,
                residual_after: trial_norm,
                forcing,
                linear_relative_residual: linear.rel_residual,
                linear_iterations: linear.iters,
                newton_step_norm: direction_norm,
                step_length,
                reduction_ratio: ratio,
                decision,
            });
            iterations += 1;
        }
    }

    let converged = residual_norm <= target;
    let diagnosis = if converged {
        None
    } else {
        Some(terminal.unwrap_or(NewtonStallDiagnosis::BudgetExhausted))
    };
    Ok(FallibleNewtonSolution {
        x,
        residual,
        report: NewtonReport {
            iterations,
            residual_norm,
            converged,
            history,
            diagnosis,
        },
        residual_evaluations,
        jacobian_evaluations,
    })
}

#[cfg(test)]
mod tests {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

    use super::*;

    fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    #[derive(Debug, PartialEq, Eq)]
    struct DomainRefusal(&'static str);

    struct CoupledQuadratic {
        refuse_negative_x0: bool,
    }

    impl FallibleNonlinearProblem for CoupledQuadratic {
        type Error = DomainRefusal;

        fn dimension(&self) -> usize {
            2
        }

        fn residual(
            &self,
            x: &[f64],
            residual: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), Self::Error> {
            let _ = cx;
            if self.refuse_negative_x0 && x[0] < 0.0 {
                return Err(DomainRefusal("negative temperature"));
            }
            residual[0] = x[0] * x[0] + x[1] - 3.0;
            residual[1] = x[0] + x[1] * x[1] - 3.0;
            Ok(())
        }

        fn jacobian_apply(
            &self,
            x: &[f64],
            direction: &[f64],
            output: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), Self::Error> {
            let _ = cx;
            output[0] = 2.0 * x[0] * direction[0] + direction[1];
            output[1] = direction[0] + 2.0 * x[1] * direction[1];
            Ok(())
        }
    }

    #[test]
    fn fallible_newton_krylov_converges_and_propagates_owner_errors_without_nan() {
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let problem = CoupledQuadratic {
                refuse_negative_x0: true,
            };
            let cfg = NewtonKrylovConfig {
                absolute_tolerance: 1e-13,
                relative_tolerance: 1e-13,
                linear_restart: 4,
                max_linear_cycles: 4,
                ..NewtonKrylovConfig::default()
            };
            let sol = solve_fallible_newton_krylov(&problem, vec![2.5, 0.5], cfg, 25, cx).unwrap();
            assert!(sol.report.converged, "report: {:?}", sol.report);
            // Positive root: x0 = x1 = (-1 + sqrt(13)) / 2
            let expected = (-1.0 + 13.0_f64.sqrt()) * 0.5;
            assert!((sol.x[0] - expected).abs() < 1e-12);
            assert!((sol.x[1] - expected).abs() < 1e-12);

            // Owner refusal travels directly as Err(Problem(DomainRefusal)), never NaN!
            let err = solve_fallible_newton_krylov(&problem, vec![-0.5, 1.0], cfg, 25, cx)
                .expect_err("domain refusal must propagate");
            assert!(matches!(
                err,
                FallibleNewtonError::Problem(DomainRefusal("negative temperature"))
            ));
        });

        // Pre-cancelled Cx immediately returns Cancelled without evaluating problem.
        gate.request();
        with_cx(&gate, |cx| {
            let problem = CoupledQuadratic {
                refuse_negative_x0: false,
            };
            let err = solve_fallible_newton_krylov(
                &problem,
                vec![2.0, 2.0],
                NewtonKrylovConfig::default(),
                25,
                cx,
            )
            .expect_err("cancelled context must abort");
            assert!(matches!(err, FallibleNewtonError::Cancelled));
        });
    }
}
