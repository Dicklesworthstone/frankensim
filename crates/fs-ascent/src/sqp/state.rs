//! Fallible, resumable small-dense SQP using the existing QP and BFGS kernels.

use super::{KktResidual, SqpReport, bfgs_update, solve_inequality_qp, violation};
use fs_exec::Cx;

/// One coherent value/derivative evaluation at a single decision point.
/// Jacobians are row-major, with one row per constraint. Inequality residuals
/// retain their original sign: `ci <= 0` is feasible, not clipped to zero.
#[derive(Debug, Clone, PartialEq)]
pub struct SqpSample {
    /// Original objective, not a penalized merit value.
    pub f: f64,
    /// Objective gradient in decision coordinates.
    pub gradient: Vec<f64>,
    /// Equality residuals.
    pub ce: Vec<f64>,
    /// Inequality residuals.
    pub ci: Vec<f64>,
    /// Equality Jacobian, `ce.len() * gradient.len()` entries.
    pub je: Vec<f64>,
    /// Inequality Jacobian, `ci.len() * gradient.len()` entries.
    pub ji: Vec<f64>,
}

/// Refusal preserves the caller's error without panic or a fake derivative.
#[derive(Debug, Clone, PartialEq)]
pub enum SqpError<E> {
    /// Original callback failure. Its attempted evaluation remains counted.
    Evaluation(E),
    /// An invalid option, start, sample, or arithmetic result.
    Invalid(&'static str),
    /// A callback changed the declared dimensions or returned malformed data.
    Shape {
        /// Payload being checked.
        field: &'static str,
        /// Required scalar count.
        expected: usize,
        /// Supplied scalar count.
        actual: usize,
    },
    /// Cancellation leaves accepted numerical state intact, but spent work is
    /// not refunded. An unfinished line search restarts on continuation.
    Cancelled,
}

impl<E: core::fmt::Display> core::fmt::Display for SqpError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Evaluation(error) => write!(f, "SQP evaluation failed: {error}"),
            Self::Invalid(what) => write!(f, "invalid SQP input or arithmetic: {what}"),
            Self::Shape { field, expected, actual } => {
                write!(f, "SQP {field} requires {expected} entries, received {actual}")
            }
            Self::Cancelled => write!(f, "SQP cancelled with accepted state retained"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for SqpError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Evaluation(error) => Some(error), _ => None }
    }
}

/// Why a segment ended. A small step alone never establishes convergence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SqpStop {
    /// All four KKT residuals satisfy the requested tolerance.
    Converged,
    /// The caller's additional accepted-step allowance was reached.
    IterationLimit,
    /// The cumulative callback ceiling was reached.
    EvaluationLimit,
    /// No usable linearized QP or merit-decreasing trial was found.
    /// This is not a nonlinear infeasibility certificate.
    Stalled,
}

/// Segment outcome and a certificate evaluated from retained derivatives.
#[derive(Debug, Clone)]
pub struct SqpRunReport {
    /// Stop attribution, independently of the certificate.
    pub stop: SqpStop,
    /// `iters` counts accepted steps; `evals` counts all sample attempts.
    pub solution: SqpReport,
}

fn poll<E>(cx: Option<&Cx<'_>>) -> Result<(), SqpError<E>> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| SqpError::Cancelled)?; }
    Ok(())
}

fn shape<E>(field: &'static str, expected: usize, actual: usize) -> Result<(), SqpError<E>> {
    if expected != actual { return Err(SqpError::Shape { field, expected, actual }); }
    Ok(())
}

impl SqpSample {
    fn validate<E>(&self, n: usize, ne: usize, ni: usize) -> Result<(), SqpError<E>> {
        shape("gradient", n, self.gradient.len())?;
        shape("equalities", ne, self.ce.len())?;
        shape("inequalities", ni, self.ci.len())?;
        shape("equality Jacobian", ne.checked_mul(n).ok_or(SqpError::Invalid("Jacobian size overflow"))?, self.je.len())?;
        shape("inequality Jacobian", ni.checked_mul(n).ok_or(SqpError::Invalid("Jacobian size overflow"))?, self.ji.len())?;
        if !self.f.is_finite() || self.gradient.iter().chain(&self.ce).chain(&self.ci)
            .chain(&self.je).chain(&self.ji).any(|x| !x.is_finite()) {
            return Err(SqpError::Invalid("sample entries must be finite; reject unavailable trials explicitly"));
        }
        Ok(())
    }

    fn lagrangian_gradient(&self, lambda: &[f64], nu: &[f64]) -> Option<Vec<f64>> {
        let n = self.gradient.len();
        let mut g = self.gradient.clone();
        for (rows, weights) in [(&self.je, lambda), (&self.ji, nu)] {
            for (row, &weight) in rows.chunks_exact(n).zip(weights) {
                for (g, &a) in g.iter_mut().zip(row) { *g += a * weight; }
            }
        }
        g.iter().all(|v| v.is_finite()).then_some(g)
    }

    fn certificate(&self, lambda: &[f64], nu: &[f64]) -> Option<KktResidual> {
        let g = self.lagrangian_gradient(lambda, nu)?;
        let mut complementarity = 0.0_f64;
        for (&c, &v) in self.ci.iter().zip(nu) {
            let product = (c * v).abs();
            if !product.is_finite() { return None; }
            complementarity = complementarity.max(product);
        }
        Some(KktResidual {
            stationarity: g.iter().map(|v| v.abs()).fold(0.0, f64::max),
            feasibility: self.ce.iter().map(|v| v.abs())
                .chain(self.ci.iter().map(|v| v.max(0.0))).fold(0.0, f64::max),
            dual_feasibility: nu.iter().map(|v| (-v).max(0.0)).fold(0.0, f64::max),
            complementarity,
        })
    }
}

/// An accepted SQP checkpoint, including its full derivative sample and BFGS
/// model. Clone retains all numerical state; no initial evaluation is repeated.
/// Public access is read-only so cached derivatives cannot be accidentally
/// paired with a different point. Callbacks must keep the same problem meaning.
#[derive(Debug, Clone)]
pub struct SqpState {
    x: Vec<f64>,
    sample: SqpSample,
    hessian: Vec<f64>,
    lambda: Vec<f64>,
    nu: Vec<f64>,
    penalty: f64,
    iterations: usize,
    evaluations: usize,
    history: Vec<f64>,
    rejected_trials: usize,
}

impl SqpState {
    /// Evaluate the initial point once. `max_kkt_dimension` bounds decision plus
    /// all constraint dimensions before the solver allocates dense matrices.
    /// Callback-owned allocations are not covered by this bound.
    ///
    /// `Ok(None)` explicitly marks an unavailable domain trial; it is an error
    /// at initialization. `Err(E)` always propagates. Successful samples must
    /// have finite, dimensionally consistent values and derivatives.
    pub fn try_new<E>(
        x0: &[f64],
        max_kkt_dimension: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<Option<SqpSample>, E>,
        cx: Option<&Cx<'_>>,
    ) -> Result<Self, SqpError<E>> {
        poll(cx)?;
        if x0.is_empty() || x0.len() > max_kkt_dimension || x0.iter().any(|x| !x.is_finite()) {
            return Err(SqpError::Invalid("start must be finite, nonempty and inside the dimension cap"));
        }
        let sample = evaluate(x0).map_err(SqpError::Evaluation)?
            .ok_or(SqpError::Invalid("initial point is outside the evaluation domain"))?;
        poll(cx)?;
        let (n, ne, ni) = (x0.len(), sample.ce.len(), sample.ci.len());
        let dim = n.checked_add(ne).and_then(|v| v.checked_add(ni))
            .ok_or(SqpError::Invalid("KKT dimension overflow"))?;
        if dim > max_kkt_dimension || dim.checked_mul(dim).is_none() {
            return Err(SqpError::Invalid("KKT dimension exceeds the dense solve cap"));
        }
        sample.validate(n, ne, ni)?;
        let mut hessian = vec![0.0; n * n];
        for i in 0..n { hessian[i * n + i] = 1.0; }
        let history = vec![sample.f];
        Ok(Self {
            x: x0.to_vec(), sample, hessian, lambda: vec![0.0; ne], nu: vec![0.0; ni],
            penalty: 10.0, iterations: 0, evaluations: 1, history, rejected_trials: 0,
        })
    }

    /// Accepted decision point.
    #[must_use]
    pub fn point(&self) -> &[f64] { &self.x }
    /// Cached value and derivative sample at the accepted point.
    #[must_use]
    pub fn sample(&self) -> &SqpSample { &self.sample }
    /// Total sample attempts, including failures and rejected trials.
    #[must_use]
    pub fn evaluations(&self) -> usize { self.evaluations }
    /// Total accepted steps.
    #[must_use]
    pub fn iterations(&self) -> usize { self.iterations }
    /// Objective values at the initial and each accepted point.
    #[must_use]
    pub fn history(&self) -> &[f64] { &self.history }
    /// Domain-unavailable or merit-rejected callback results.
    #[must_use]
    pub fn rejected_trials(&self) -> usize { self.rejected_trials }

    fn report<E>(&self, stop: SqpStop, tol: f64) -> Result<SqpRunReport, SqpError<E>> {
        let kkt = self.sample.certificate(&self.lambda, &self.nu)
            .ok_or(SqpError::Invalid("non-finite KKT certificate"))?;
        let converged = kkt.within_tolerance(tol);
        Ok(SqpRunReport { stop, solution: SqpReport {
            x: self.x.clone(), f: self.sample.f, kkt, lambda: self.lambda.clone(), nu: self.nu.clone(),
            iters: self.iterations, evals: self.evaluations, converged,
        } })
    }

    /// Continue for an additional accepted-step allowance. `max_evals` is a
    /// cumulative hard ceiling, including the constructor evaluation. It is
    /// checked before every trial, never after an overshoot. A final-budget
    /// trial can still be accepted. A lower limit cannot undo already spent work.
    ///
    /// Only work counters change on evaluation failure or cancellation. Partial
    /// searches are discarded; resuming may repeat probes, not accepted steps.
    /// Polls bracket the bounded dense QP/BFGS phases and each callback. These
    /// dense kernels are not internally interruptible; no latency bound is claimed.
    pub fn try_run<E>(
        &mut self,
        evaluate: &mut impl FnMut(&[f64]) -> Result<Option<SqpSample>, E>,
        tol: f64,
        additional_iters: usize,
        max_evals: usize,
        cx: Option<&Cx<'_>>,
    ) -> Result<SqpRunReport, SqpError<E>> {
        if !tol.is_finite() || tol <= 0.0 { return Err(SqpError::Invalid("tolerance must be finite and positive")); }
        let n = self.x.len();
        let mut completed = 0usize;
        loop {
            poll(cx)?;
            if self.evaluations >= max_evals { return self.report(SqpStop::EvaluationLimit, tol); }
            let current = self.report(SqpStop::IterationLimit, tol)?;
            if current.solution.converged { return self.report(SqpStop::Converged, tol); }
            if completed == additional_iters { return Ok(current); }
            let Some(qp) = solve_inequality_qp(&self.hessian, &self.sample.gradient,
                &self.sample.je, &self.sample.ce, &self.sample.ji, &self.sample.ci) else {
                poll(cx)?;
                return self.report(SqpStop::Stalled, tol);
            };
            poll(cx)?;
            let Some(kkt) = self.sample.certificate(&qp.lambda, &qp.nu) else {
                return self.report(SqpStop::Stalled, tol);
            };
            let dnorm = qp.d.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if dnorm < tol && kkt.within_tolerance(tol) {
                self.lambda = qp.lambda;
                self.nu = qp.nu;
                return self.report(SqpStop::Converged, tol);
            }
            let dual_norm = qp.lambda.iter().chain(&qp.nu).map(|v| v.abs()).fold(0.0, f64::max);
            let penalty = self.penalty.max(1.1 * dual_norm);
            let v0 = violation(&self.sample.ce, &self.sample.ci);
            let gd: f64 = self.sample.gradient.iter().zip(&qp.d).map(|(g, d)| g * d).sum();
            let slope = gd - penalty * v0;
            if !penalty.is_finite() || !slope.is_finite() || slope >= 0.0 {
                return self.report(SqpStop::Stalled, tol);
            }
            let mut alpha = 1.0_f64;
            let mut accepted = false;
            for _ in 0..40 {
                poll(cx)?;
                if self.evaluations >= max_evals { return self.report(SqpStop::EvaluationLimit, tol); }
                let xt: Vec<f64> = self.x.iter().zip(&qp.d).map(|(x, d)| alpha.mul_add(*d, *x)).collect();
                if xt == self.x { return self.report(SqpStop::Stalled, tol); }
                if xt.iter().any(|v| !v.is_finite()) { alpha *= 0.5; continue; }
                self.evaluations += 1;
                let trial = evaluate(&xt).map_err(SqpError::Evaluation)?;
                poll(cx)?;
                if let Some(trial) = trial {
                    trial.validate(n, self.lambda.len(), self.nu.len())?;
                    let change = (trial.f - self.sample.f) + penalty * (violation(&trial.ce, &trial.ci) - v0);
                    if change.is_finite() && change < 0.0 && change <= 1e-4 * alpha * slope {
                        let old_g = self.sample.lagrangian_gradient(&qp.lambda, &qp.nu)
                            .ok_or(SqpError::Invalid("non-finite Lagrangian gradient"))?;
                        let new_g = trial.lagrangian_gradient(&qp.lambda, &qp.nu)
                            .ok_or(SqpError::Invalid("non-finite trial Lagrangian gradient"))?;
                        if trial.certificate(&qp.lambda, &qp.nu).is_none() {
                            return Err(SqpError::Invalid("non-finite trial KKT certificate"));
                        }
                        let s: Vec<f64> = xt.iter().zip(&self.x).map(|(a, b)| a - b).collect();
                        let y: Vec<f64> = new_g.iter().zip(&old_g).map(|(a, b)| a - b).collect();
                        let mut next_b = self.hessian.clone();
                        if s.iter().chain(&y).all(|v| v.is_finite()) {
                            bfgs_update(&mut next_b, n, &s, &y);
                        }
                        // A failed curvature update cannot poison an otherwise
                        // valid accepted iterate; retain the previous model.
                        if next_b.iter().any(|v| !v.is_finite()) { next_b.clone_from(&self.hessian); }
                        poll(cx)?;
                        self.x = xt;
                        self.sample = trial;
                        self.hessian = next_b;
                        self.lambda = qp.lambda;
                        self.nu = qp.nu;
                        self.penalty = penalty;
                        self.iterations += 1;
                        self.history.push(self.sample.f);
                        completed += 1;
                        accepted = true;
                        break;
                    }
                }
                self.rejected_trials += 1;
                alpha *= 0.5;
            }
            if !accepted { return self.report(SqpStop::Stalled, tol); }
        }
    }
}
