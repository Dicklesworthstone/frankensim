//! Constrained optimization of the live IR with shared reverse derivatives.

use crate::runner::Packing;
use crate::sqp::{SqpError, SqpReport, SqpRunReport, SqpSample, SqpState, SqpStop};
use fs_exec::Cx;
use fs_opt::reverse::ReverseError;
use fs_opt::{ConstraintKind, Manifold, OptError, ReverseProblem, ReverseProblemError};

/// Structural refusal or the original typed evaluation/optimizer failure.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseSqpError {
    /// Dense Euclidean SQP cannot use manifold point storage as flat parameters.
    NonEuclidean {
        /// Variable in declaration order.
        variable: usize,
    },
    /// The supplied point does not match the complete declared packing.
    PackedPointLength {
        /// Required scalar count.
        expected: usize,
        /// Supplied scalar count.
        actual: usize,
    },
    /// Decision plus constraint dimensions exceed the explicit dense-QP cap.
    DimensionCap,
    /// Original optimizer or callback error.
    Optimizer(SqpError<ReverseProblemError>),
    /// Cancellation preserves accepted state and charges any attempted sample.
    Cancelled,
}

impl From<SqpError<ReverseProblemError>> for ReverseSqpError {
    fn from(error: SqpError<ReverseProblemError>) -> Self {
        match error {
            SqpError::Cancelled | SqpError::Evaluation(ReverseProblemError::Reverse(
                ReverseError::Evaluation(OptError::Cancelled),
            )) => Self::Cancelled,
            other => Self::Optimizer(other),
        }
    }
}

impl core::fmt::Display for ReverseSqpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonEuclidean { variable } => write!(f, "reverse SQP variable {variable} is not Euclidean"),
            Self::PackedPointLength { expected, actual } => write!(f, "reverse SQP needs {expected} coordinates, received {actual}"),
            Self::DimensionCap => write!(f, "reverse SQP exceeds its explicit dense dimension cap"),
            Self::Optimizer(error) => write!(f, "{error}"),
            Self::Cancelled => write!(f, "reverse SQP cancelled with accepted state retained"),
        }
    }
}

impl std::error::Error for ReverseSqpError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Optimizer(error) => Some(error), _ => None }
    }
}

/// Solution plus multipliers in the original IR constraint declaration order.
#[derive(Debug, Clone)]
pub struct ReverseSqpReport {
    /// Segment stop attribution.
    pub stop: SqpStop,
    /// Accepted decision, weighted objective, grouped duals and KKT residuals.
    pub solution: SqpReport,
    /// Equality and inequality multipliers interleaved exactly as declared,
    /// suitable for `ProblemEvaluation::lagrangian_gradient`.
    pub constraint_multipliers: Vec<f64>,
}

fn capacity<T>(len: usize) -> Result<Vec<T>, ReverseProblemError> {
    let mut values = Vec::new();
    values.try_reserve_exact(len).map_err(|_| OptError::RuntimeAllocationRefused {
        path: "reverse-sqp/sample", node: None, variable: None,
        elements: len as u64, element_bytes: core::mem::size_of::<T>() as u64,
    })?;
    Ok(values)
}

fn evaluate(
    oracle: &ReverseProblem<'_>,
    packing: &Packing,
    point: &[f64],
    cx: Option<&Cx<'_>>,
) -> Result<SqpSample, ReverseProblemError> {
    let bindings = packing.unpack(point);
    let tape = match cx {
        Some(cx) => oracle.evaluate_cancellable(&bindings, cx)?,
        None => oracle.evaluate(&bindings)?,
    };
    let blocks = match cx {
        Some(cx) => tape.objective_gradient_cancellable(cx)?,
        None => tape.objective_gradient()?,
    };
    let mut gradient = capacity(packing.dim)?;
    for block in blocks { gradient.extend(block); }
    let constraints = oracle.problem().constraints();
    let ne = constraints.iter().filter(|c| c.kind == ConstraintKind::EqZero).count();
    let ni = constraints.len() - ne;
    // Constructor preflight bounds n + ne + ni and its square, hence all
    // products below, before the first primal or Jacobian allocation.
    let mut sample = SqpSample {
        f: tape.objective_value(), gradient,
        ce: capacity(ne)?, ci: capacity(ni)?,
        je: capacity(ne * packing.dim)?, ji: capacity(ni * packing.dim)?,
    };
    let mut seeds = capacity(constraints.len())?;
    seeds.resize(constraints.len(), 0.0);
    for (index, constraint) in constraints.iter().enumerate() {
        seeds[index] = 1.0;
        let blocks = match cx {
            Some(cx) => tape.constraint_pullback_cancellable(&seeds, cx)?,
            None => tape.constraint_pullback(&seeds)?,
        };
        seeds[index] = 0.0;
        let (values, rows) = match constraint.kind {
            ConstraintKind::EqZero => (&mut sample.ce, &mut sample.je),
            ConstraintKind::LeZero => (&mut sample.ci, &mut sample.ji),
        };
        values.push(tape.constraint_values()[index]);
        for block in blocks { rows.extend(block); }
    }
    Ok(sample)
}

fn numerical_trial_failure(error: &ReverseProblemError) -> bool {
    matches!(error,
        ReverseProblemError::NonFiniteObjective { .. }
        | ReverseProblemError::Reverse(ReverseError::NonFiniteAdjoint { .. })
        | ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::EvalNonFinite { .. }))
    )
}

/// A constrained solve bound to one immutable compiled problem. Clone retains
/// the SQP model, sample, history, work and diagnostics while sharing the oracle.
/// Callers supply no objective, constraint or derivative callbacks.
///
/// This is small-dense, smooth Euclidean SQP. Both constraint kinds are solved;
/// non-Euclidean variables refuse rather than having their geometry discarded.
#[derive(Debug, Clone)]
pub struct ReverseSqpStudy<'oracle, 'problem> {
    oracle: &'oracle ReverseProblem<'problem>,
    packing: Packing,
    state: SqpState,
    evaluation_limit: usize,
    domain_rejections: usize,
    last_rejection: Option<ReverseProblemError>,
}

impl<'oracle, 'problem> ReverseSqpStudy<'oracle, 'problem> {
    /// Validate dimensions before evaluating the initial point once. The
    /// problem's positive/unlimited evaluation budget includes that attempt.
    /// `max_kkt_dimension` bounds total decision plus declared constraints.
    /// The compiled oracle's own limits independently bound the reverse tape.
    pub fn new(
        oracle: &'oracle ReverseProblem<'problem>,
        point: &[f64],
        max_kkt_dimension: usize,
        cx: Option<&Cx<'_>>,
    ) -> Result<Self, ReverseSqpError> {
        if let Some(cx) = cx { cx.checkpoint().map_err(|_| ReverseSqpError::Cancelled)?; }
        let mut n = 0usize;
        for (variable, declaration) in oracle.problem().vars().iter().enumerate() {
            let Manifold::Rn { dim } = declaration.manifold else {
                return Err(ReverseSqpError::NonEuclidean { variable });
            };
            n = n.checked_add(dim as usize).ok_or(ReverseSqpError::DimensionCap)?;
        }
        let dim = n.checked_add(oracle.problem().constraints().len()).ok_or(ReverseSqpError::DimensionCap)?;
        if dim > max_kkt_dimension || dim.checked_mul(dim).is_none() {
            return Err(ReverseSqpError::DimensionCap);
        }
        if point.len() != n {
            return Err(ReverseSqpError::PackedPointLength { expected: n, actual: point.len() });
        }
        let packing = Packing::new(oracle.problem());
        let evaluation_limit = oracle.problem().budget().limit.maximum().map_or(
            usize::MAX, |cap| usize::try_from(cap.get()).unwrap_or(usize::MAX),
        );
        let state = SqpState::try_new(point, max_kkt_dimension, &mut |x| {
            evaluate(oracle, &packing, x, cx).map(Some)
        }, cx)?;
        Ok(Self { oracle, packing, state, evaluation_limit, domain_rejections: 0, last_rejection: None })
    }

    /// Read-only accepted point, cached derivative sample, history and work.
    #[must_use]
    pub fn optimizer(&self) -> &SqpState { &self.state }
    /// Numerical trial rejections, excluding ordinary merit rejections.
    #[must_use]
    pub fn domain_rejections(&self) -> usize { self.domain_rejections }
    /// Original numerical error from the most recently unavailable trial.
    #[must_use]
    pub fn last_rejection(&self) -> Option<&ReverseProblemError> { self.last_rejection.as_ref() }

    /// Continue within the immutable problem's cumulative evaluation limit.
    pub fn run(&mut self, tol: f64, additional_iters: usize, cx: Option<&Cx<'_>>)
        -> Result<ReverseSqpReport, ReverseSqpError> {
        self.run_with_budget(tol, additional_iters, self.evaluation_limit, cx)
    }

    /// Install an optional stricter cumulative ceiling without raising the
    /// problem's limit. Each attempt performs at most one primal sweep and
    /// `1 + constraint_count` reverse sweeps, with no coordinate perturbations.
    /// Reports and checkpoint reads perform no graph evaluations.
    ///
    /// Invalid initial arithmetic is an error. Only non-finite *trial* values
    /// or derivatives are backtracked; other typed failures propagate. A
    /// cancelled/incomplete search retains accepted state and spent work, but
    /// its probes may be repeated on resume. Dense kernels and Packing remain
    /// whole phases; reverse sweeps poll through the existing Cx path.
    pub fn run_with_budget(
        &mut self, tol: f64, additional_iters: usize, max_evals: usize, cx: Option<&Cx<'_>>,
    ) -> Result<ReverseSqpReport, ReverseSqpError> {
        let oracle = self.oracle;
        let packing = &self.packing;
        let rejected = &mut self.domain_rejections;
        let last = &mut self.last_rejection;
        let mut callback = |x: &[f64]| match evaluate(oracle, packing, x, cx) {
            Ok(sample) => Ok(Some(sample)),
            Err(error) if numerical_trial_failure(&error) => {
                *rejected += 1;
                *last = Some(error);
                Ok(None)
            }
            Err(error) => Err(error),
        };
        let SqpRunReport { stop, solution } = self.state.try_run(
            &mut callback, tol, additional_iters, max_evals.min(self.evaluation_limit), cx,
        )?;
        let mut constraint_multipliers = Vec::with_capacity(oracle.problem().constraints().len());
        let (mut e, mut i) = (0, 0);
        for constraint in oracle.problem().constraints() {
            constraint_multipliers.push(match constraint.kind {
                ConstraintKind::EqZero => { let v = solution.lambda[e]; e += 1; v }
                ConstraintKind::LeZero => { let v = solution.nu[i]; i += 1; v }
            });
        }
        Ok(ReverseSqpReport { stop, solution, constraint_multipliers })
    }
}
