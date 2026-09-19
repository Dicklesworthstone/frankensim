//! Fallible, hard-budget execution over the existing L-BFGS checkpoint.

use super::{LbfgsReport, LbfgsState, inf_norm};
use crate::stop::{StopObservation, StopReason, StopRule};
use crate::wolfe::try_strong_wolfe_with_budget;
use std::collections::VecDeque;

/// A refused input or failed objective evaluation. The accepted checkpoint
/// remains usable after a run error; `evals` includes failed callback attempts.
#[derive(Debug, Clone, PartialEq)]
pub enum LbfgsError<E> {
    /// A dimension-independent solver precondition was violated.
    InvalidInput(&'static str),
    /// The callback did not return one derivative per decision coordinate.
    GradientLength {
        /// Decision dimension.
        expected: usize,
        /// Returned gradient length.
        actual: usize,
    },
    /// A required finite scalar or vector component was not finite.
    NonFinite {
        /// Input or computation that failed.
        what: &'static str,
        /// Vector component, or None for a scalar.
        component: Option<usize>,
        /// Exact offending value.
        bits: u64,
    },
    /// The callback's original typed refusal, without panic conversion.
    Evaluation(E),
}

impl<E: core::fmt::Display> core::fmt::Display for LbfgsError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(what) => write!(f, "L-BFGS: {what}"),
            Self::GradientLength { expected, actual } => {
                write!(f, "L-BFGS needs {expected} gradient components, received {actual}")
            }
            Self::NonFinite { what, component, bits } => {
                write!(f, "L-BFGS non-finite {what} at {component:?}: {bits:#018x}")
            }
            Self::Evaluation(error) => write!(f, "L-BFGS objective: {error}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for LbfgsError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Evaluation(error) => Some(error),
            _ => None,
        }
    }
}

fn finite<E>(what: &'static str, values: &[f64]) -> Result<(), LbfgsError<E>> {
    for (component, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(LbfgsError::NonFinite {
                what, component: Some(component), bits: value.to_bits(),
            });
        }
    }
    Ok(())
}

fn checked_fg<E>(f: f64, g: &[f64], n: usize, trial: bool) -> Result<(), LbfgsError<E>> {
    if !f.is_finite() && !(trial && f == f64::INFINITY) {
        return Err(LbfgsError::NonFinite { what: "objective", component: None, bits: f.to_bits() });
    }
    if g.len() != n {
        return Err(LbfgsError::GradientLength { expected: n, actual: g.len() });
    }
    finite("gradient", g)
}

// Budgets are resource ceilings even under an All whose other predicates are
// false. Preserve the normal stop algebra for the non-resource predicates.
fn budget(rule: &StopRule) -> usize {
    match rule {
        StopRule::Budget(maximum) => *maximum,
        StopRule::Any(rules) | StopRule::All(rules) => {
            rules.iter().map(budget).min().unwrap_or(usize::MAX)
        }
        _ => usize::MAX,
    }
}

impl LbfgsState {
    /// Fallible initialization; validates point and memory before one callback.
    /// A successful checkpoint starts with `evals == 1`. The supplied callback
    /// is attempted once even if it returns an error; no partial state escapes.
    /// A positive problem evaluation limit therefore always admits this call.
    pub fn try_new<E>(
        x0: &[f64],
        memory: usize,
        fg: &mut dyn FnMut(&[f64]) -> Result<(f64, Vec<f64>), E>,
    ) -> Result<Self, LbfgsError<E>> {
        if x0.is_empty() || memory == 0 {
            return Err(LbfgsError::InvalidInput("point and memory must be nonempty"));
        }
        finite("initial point", x0)?;
        let (f, g) = fg(x0).map_err(LbfgsError::Evaluation)?;
        checked_fg(f, &g, x0.len(), false)?;
        Ok(Self {
            x: x0.to_vec(), f, g, pairs: VecDeque::new(), memory,
            iters: 0, evals: 1, history: vec![f],
        })
    }

    fn fallible_report(&self, reason: StopReason) -> LbfgsReport {
        LbfgsReport {
            reason, grad_norm: inf_norm(&self.g), f: self.f,
            iters: self.iters, evals: self.evals,
        }
    }

    fn fallible_stop(&self, rule: &StopRule, ceiling: usize) -> Option<StopReason> {
        if self.evals >= ceiling {
            return Some(StopReason::Budget);
        }
        rule.check(&StopObservation {
            grad_norm: inf_norm(&self.g), objective: self.f,
            evals: self.evals, history: &self.history,
        })
    }

    /// Run with fallible callbacks and an absolute cumulative evaluation ceiling.
    /// `usize::MAX` explicitly opts out of that ceiling; every Budget leaf in
    /// `rule` is also a hard ceiling, including leaves nested under All.
    ///
    /// The same two-loop recursion and Wolfe probes as the legacy run are used.
    /// Only accepted steps change x/f/g/history/curvature. Interrupted searches
    /// leave these intact but retain their spent evaluations. A subsequent call
    /// restarts that search; it does not refund trials or claim in-flight replay.
    /// Split runs at accepted-iteration boundaries retain the complete state.
    ///
    /// A +infinity trial with a finite gradient is the existing rejected-domain
    /// convention. Other malformed observations return a typed refusal. Callback
    /// panics and allocation aborts are not caught. Legacy `new`/`run` are unchanged.
    pub fn try_run<E>(
        &mut self,
        fg: &mut dyn FnMut(&[f64]) -> Result<(f64, Vec<f64>), E>,
        rule: &StopRule,
        max_iters: usize,
        max_evals: usize,
    ) -> Result<LbfgsReport, LbfgsError<E>> {
        if self.x.is_empty() || self.memory == 0 {
            return Err(LbfgsError::InvalidInput("point and memory must be nonempty"));
        }
        finite("accepted point", &self.x)?;
        checked_fg(self.f, &self.g, self.x.len(), false)?;
        let ceiling = max_evals.min(budget(rule));
        for _ in 0..max_iters {
            if let Some(reason) = self.fallible_stop(rule, ceiling) {
                return Ok(self.fallible_report(reason));
            }
            let mut d = self.direction();
            let mut slope: f64 = d.iter().zip(&self.g).map(|(di, gi)| di * gi).sum();
            if !slope.is_finite() || slope >= 0.0 {
                d = self.g.iter().map(|g| -g).collect();
                slope = -self.g.iter().map(|g| g * g).sum::<f64>();
            }
            if !slope.is_finite() || slope >= 0.0 {
                let reason = if self.g.iter().all(|g| *g == 0.0) {
                    StopReason::GradNorm
                } else {
                    StopReason::Stall
                };
                return Ok(self.fallible_report(reason));
            }
            let mut spent = 0usize;
            let mut last = None;
            let outcome = {
                let mut phi = |alpha: f64| -> Result<(f64, f64), LbfgsError<E>> {
                    let xt: Vec<f64> = self.x.iter().zip(&d)
                        .map(|(xi, di)| alpha.mul_add(*di, *xi)).collect();
                    finite("trial point", &xt)?;
                    // Count before invoking: failed/cancelled evaluations cost work.
                    spent += 1;
                    let (f, g) = fg(&xt).map_err(LbfgsError::Evaluation)?;
                    checked_fg(f, &g, xt.len(), true)?;
                    if f == f64::INFINITY {
                        return Ok((f64::INFINITY, 0.0));
                    }
                    let derivative: f64 = g.iter().zip(&d).map(|(gi, di)| gi * di).sum();
                    if !derivative.is_finite() {
                        return Err(LbfgsError::NonFinite {
                            what: "directional derivative", component: None, bits: derivative.to_bits(),
                        });
                    }
                    last = Some((f, g, xt));
                    Ok((f, derivative))
                };
                try_strong_wolfe_with_budget(
                    &mut phi, self.f, slope, 1.0, 1e-4, 0.9, ceiling - self.evals,
                )
            };
            self.evals += spent; // the pre-call ceiling makes this addition safe
            let outcome = outcome?;
            if !outcome.success {
                return Ok(self.fallible_report(if self.evals == ceiling {
                    StopReason::Budget
                } else {
                    StopReason::Stall
                }));
            }
            let (f, g, x) = last.expect("a successful Wolfe search has a finite last probe");
            if x == self.x {
                return Ok(self.fallible_report(StopReason::Stall));
            }
            let s: Vec<f64> = x.iter().zip(&self.x).map(|(a, b)| a - b).collect();
            let y: Vec<f64> = g.iter().zip(&self.g).map(|(a, b)| a - b).collect();
            let sy: f64 = s.iter().zip(&y).map(|(a, b)| a * b).sum();
            if sy.is_finite() && (1.0 / sy).is_finite()
                && s.iter().chain(&y).all(|v| v.is_finite())
                && sy > 1e-14 * inf_norm(&s) * inf_norm(&y).max(1e-30)
            {
                while self.pairs.len() >= self.memory {
                    self.pairs.pop_front();
                }
                self.pairs.push_back((s, y, 1.0 / sy));
            }
            self.x = x;
            self.f = f;
            self.g = g;
            self.iters += 1;
            self.history.push(f);
        }
        Ok(self.fallible_report(
            self.fallible_stop(rule, ceiling).unwrap_or(StopReason::IterationCap),
        ))
    }
}
