//! Actual L-BFGS solves over the live reverse-mode problem IR, without finite
//! differences or user-written callback adapters. This opt-in path leaves the
//! existing finite-difference Study and its replay semantics unchanged.

use crate::lbfgs::{LbfgsError, LbfgsReport, LbfgsState};
use crate::runner::Packing;
use crate::stop::{StopReason, StopRule};
use fs_exec::Cx;
use fs_opt::reverse::ReverseError;
use fs_opt::{Manifold, OptError, ReverseProblem, ReverseProblemError};

/// A structural refusal, cancellation, or original evaluator/optimizer error.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseStudyError {
    /// Unconstrained L-BFGS must not silently drop declared constraints.
    ConstraintsUnsupported,
    /// This Euclidean L-BFGS driver cannot transport manifold curvature pairs.
    NonEuclidean {
        /// Variable in declaration order.
        variable: usize,
    },
    /// The initial vector does not bind all declared variable points.
    PackedPointLength {
        /// Required point-storage length.
        expected: usize,
        /// Supplied point-storage length.
        actual: usize,
    },
    /// The original typed optimizer or reverse-evaluation refusal.
    Optimizer(LbfgsError<ReverseProblemError>),
    /// Cancellation was observed; the last accepted checkpoint remains usable.
    Cancelled,
}

impl From<LbfgsError<ReverseProblemError>> for ReverseStudyError {
    fn from(error: LbfgsError<ReverseProblemError>) -> Self {
        match error {
            LbfgsError::Evaluation(ReverseProblemError::Reverse(
                ReverseError::Evaluation(OptError::Cancelled),
            )) => Self::Cancelled,
            other => Self::Optimizer(other),
        }
    }
}

impl core::fmt::Display for ReverseStudyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ConstraintsUnsupported => {
                write!(f, "reverse L-BFGS does not solve constrained problems")
            }
            Self::NonEuclidean { variable } => {
                write!(f, "reverse L-BFGS variable {variable} is not Euclidean")
            }
            Self::PackedPointLength { expected, actual } => {
                write!(f, "reverse study needs {expected} point coordinates, received {actual}")
            }
            Self::Optimizer(error) => write!(f, "{error}"),
            Self::Cancelled => write!(f, "reverse study cancelled at its last accepted iterate"),
        }
    }
}

impl std::error::Error for ReverseStudyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Optimizer(error) => Some(error),
            _ => None,
        }
    }
}

fn poll(cx: Option<&Cx<'_>>) -> Result<(), ReverseStudyError> {
    if let Some(cx) = cx {
        cx.checkpoint().map_err(|_| ReverseStudyError::Cancelled)?;
    }
    Ok(())
}

fn value_gradient(
    oracle: &ReverseProblem<'_>,
    packing: &Packing,
    point: &[f64],
    cx: Option<&Cx<'_>>,
) -> Result<(f64, Vec<f64>), ReverseProblemError> {
    let bindings = packing.unpack(point);
    let tape = match cx {
        Some(cx) => oracle.evaluate_cancellable(&bindings, cx)?,
        None => oracle.evaluate(&bindings)?,
    };
    let gradient = match cx {
        Some(cx) => tape.objective_gradient_cancellable(cx)?,
        None => tape.objective_gradient()?,
    };
    // Rn point and parameter coordinates coincide. Packing supplies variable
    // order; the reverse program supplies one full gradient per variable.
    let mut packed = Vec::with_capacity(packing.dim);
    for block in gradient {
        packed.extend(block);
    }
    Ok((tape.objective_value(), packed))
}

fn unavailable_finite_trial(error: &ReverseProblemError) -> bool {
    matches!(error,
        ReverseProblemError::NonFiniteObjective { .. }
        | ReverseProblemError::Reverse(ReverseError::NonFiniteAdjoint { .. })
        | ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::EvalNonFinite { .. }))
    )
}

/// Resumable, budget-enforced L-BFGS over a compiled algebraic IR problem.
///
/// The borrowed immutable oracle fixes objective semantics for the lifetime of
/// the checkpoint. Resume does not take a new problem, preventing stale caches
/// from crossing problem meanings. Clone retains curvature, accepted history,
/// evaluation accounting and rejected-trial diagnostics. The compiled program
/// is shared rather than cloned or recompiled.
///
/// All variables must be Rn and constraints must be absent. Other manifolds
/// require Riemannian curvature transport, not flat coordinate updates. Unsupported
/// tags, kinks and physics/UQ execution are already refused by ReverseProblem.
#[derive(Debug, Clone)]
pub struct ReverseStudy<'oracle, 'problem> {
    oracle: &'oracle ReverseProblem<'problem>,
    packing: Packing,
    state: LbfgsState,
    evaluation_limit: usize,
    rejected_trials: usize,
    last_rejection: Option<ReverseProblemError>,
}

impl<'oracle, 'problem> ReverseStudy<'oracle, 'problem> {
    /// Bind a solve to a compiled problem and evaluate the initial point once.
    /// The problem's explicit positive/unlimited evaluation budget includes
    /// this initial primal+gradient evaluation. `None` opts out of cancellation.
    /// Invalid initial arithmetic is an error, never a line-search rejection.
    pub fn new(
        oracle: &'oracle ReverseProblem<'problem>,
        point: &[f64],
        memory: usize,
        cx: Option<&Cx<'_>>,
    ) -> Result<Self, ReverseStudyError> {
        poll(cx)?;
        if !oracle.problem().constraints().is_empty() {
            return Err(ReverseStudyError::ConstraintsUnsupported);
        }
        for (variable, declaration) in oracle.problem().vars().iter().enumerate() {
            if !matches!(declaration.manifold, Manifold::Rn { .. }) {
                return Err(ReverseStudyError::NonEuclidean { variable });
            }
        }
        let packing = Packing::new(oracle.problem());
        if point.len() != packing.dim {
            return Err(ReverseStudyError::PackedPointLength {
                expected: packing.dim, actual: point.len(),
            });
        }
        let evaluation_limit = oracle.problem().budget().limit.maximum().map_or(
            usize::MAX, |cap| usize::try_from(cap.get()).unwrap_or(usize::MAX),
        );
        let state = LbfgsState::try_new(point, memory, &mut |x| {
            value_gradient(oracle, &packing, x, cx)
        })?;
        poll(cx)?;
        Ok(Self {
            oracle, packing, state, evaluation_limit,
            rejected_trials: 0, last_rejection: None,
        })
    }

    /// Read-only accepted optimizer state, including gradient and work counts.
    /// No mutable access is exposed that could invalidate the retained IR value.
    #[must_use]
    pub fn optimizer(&self) -> &LbfgsState {
        &self.state
    }

    /// Count of non-finite trial valuations rejected by the line search.
    #[must_use]
    pub fn rejected_trials(&self) -> usize {
        self.rejected_trials
    }

    /// Most recent rejected trial's original numerical error, if any.
    #[must_use]
    pub fn last_rejection(&self) -> Option<&ReverseProblemError> {
        self.last_rejection.as_ref()
    }

    /// Run up to `additional_iters` accepted steps, using one primal sweep and
    /// one reverse sweep per objective callback. Both the problem's budget and
    /// all Budget leaves in `rule` are enforced before each trial, including zoom.
    ///
    /// A trial outside the finite primal/derivative domain is rejected using
    /// the existing +infinity barrier convention. Its zero placeholder gradient
    /// is never accepted; the original refusal and trial count remain inspectable.
    /// All other errors propagate. No interval/global-convergence claim is made.
    ///
    /// Cancellation is polled between iterations and inside reverse evaluations.
    /// A cancelled in-flight search is discarded, not refunded: accepted state
    /// and spent work remain intact. Resume restarts that search. Splits between
    /// completed iterations reproduce the uninterrupted trajectory and accounting;
    /// mid-search cancellation does not promise identical total work on replay.
    /// Existing packing/two-loop operations remain whole non-interruptible phases.
    pub fn run(
        &mut self,
        rule: &StopRule,
        additional_iters: usize,
        cx: Option<&Cx<'_>>,
    ) -> Result<LbfgsReport, ReverseStudyError> {
        poll(cx)?;
        let oracle = self.oracle;
        let packing = &self.packing;
        let rejected = &mut self.rejected_trials;
        let last_rejection = &mut self.last_rejection;
        let mut fg = |point: &[f64]| {
            match value_gradient(oracle, packing, point, cx) {
                Err(error) if unavailable_finite_trial(&error) => {
                    *rejected += 1;
                    *last_rejection = Some(error);
                    Ok((f64::INFINITY, vec![0.0; packing.dim]))
                }
                result => result,
            }
        };
        let mut report = self.state.try_run(&mut fg, rule, 0, self.evaluation_limit)?;
        for _ in 0..additional_iters {
            if report.reason != StopReason::IterationCap {
                break;
            }
            poll(cx)?;
            report = self.state.try_run(&mut fg, rule, 1, self.evaluation_limit)?;
        }
        poll(cx)?;
        Ok(report)
    }
}
