//! Parameter/BC uncertainty propagation and sampling execution under explicit budgets.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`
//!
//! Runs bounded Philox Monte Carlo for declared probability measures, including
//! correlated Gaussian inputs. Other method/kind declarations refuse before
//! evaluating the model. Empirical statistics never mint Verified evidence.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;

/// Classification of uncertainty for a physical parameter or boundary condition.
#[derive(Debug, Clone, PartialEq)]
pub enum UncertaintyKind {
    /// Pure aleatory Gaussian distribution (known mean and standard deviation).
    AleatoryGaussian {
        /// Distribution mean.
        mean: f64,
        /// Distribution standard deviation.
        std_dev: f64,
    },
    /// Pure aleatory uniform distribution on [lo, hi].
    AleatoryUniform {
        /// Lower interval bound.
        lo: f64,
        /// Upper interval bound.
        hi: f64,
    },
    /// Epistemic interval enclosure (no probability measure claimed).
    EpistemicInterval {
        /// Lower interval bound.
        lo: f64,
        /// Upper interval bound.
        hi: f64,
    },
    /// Statistical confidence estimate with stated confidence level (e.g. 95%).
    StatisticalConfidence {
        /// Point estimate.
        estimate: f64,
        /// Half-width margin of error.
        half_width: f64,
        /// Confidence level in [0, 1].
        confidence: f64,
    },
    /// Unstated or uncharacterized uncertainty.
    Unstated,
}

impl UncertaintyKind {
    /// Whether this kind carries a well-defined probability measure.
    #[must_use]
    pub fn has_probability_measure(&self) -> bool {
        matches!(
            self,
            Self::AleatoryGaussian { .. } | Self::AleatoryUniform { .. }
        )
    }

    /// Machine-readable kind label.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        match self {
            Self::AleatoryGaussian { .. } => "aleatory-gaussian",
            Self::AleatoryUniform { .. } => "aleatory-uniform",
            Self::EpistemicInterval { .. } => "epistemic-interval",
            Self::StatisticalConfidence { .. } => "statistical-confidence",
            Self::Unstated => "unstated",
        }
    }
}

/// A declared parameter with typed uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterUncertainty {
    /// Parameter identifier name.
    pub name: String,
    /// Typed classification of uncertainty.
    pub kind: UncertaintyKind,
    /// Physical unit.
    pub unit: String,
}

impl ParameterUncertainty {
    /// Create a new parameter with aleatory Gaussian distribution.
    #[must_use]
    pub fn gaussian(
        name: impl Into<String>,
        mean: f64,
        std_dev: f64,
        unit: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            kind: UncertaintyKind::AleatoryGaussian { mean, std_dev },
            unit: unit.into(),
        }
    }

    /// Create a new parameter with aleatory uniform distribution.
    #[must_use]
    pub fn uniform(name: impl Into<String>, lo: f64, hi: f64, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            kind: UncertaintyKind::AleatoryUniform { lo, hi },
            unit: unit.into(),
        }
    }

    /// Create a new parameter with epistemic interval bounds.
    #[must_use]
    pub fn interval(name: impl Into<String>, lo: f64, hi: f64, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            kind: UncertaintyKind::EpistemicInterval { lo, hi },
            unit: unit.into(),
        }
    }
}

/// Dependence and correlation model between uncertain parameters.
#[derive(Debug, Clone, PartialEq)]
pub enum CorrelationModel {
    /// Mutually independent marginals.
    Independent,
    /// Correlations without a declared joint distribution. Sampling refuses:
    /// Gaussian marginals alone do not imply a jointly Gaussian law.
    Correlated {
        /// Normalized symmetric correlation matrix.
        matrix: Vec<Vec<f64>>,
    },
    /// Explicit multivariate Gaussian joint law with Gaussian marginals.
    JointGaussian {
        /// Correlation matrix in parameter declaration order; numerically PSD.
        matrix: Vec<Vec<f64>>,
    },
    /// Unknown correlation across parameters (forbids joint probability propagation).
    Unknown,
}

/// Propagation method for uncertainty analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropagationMethod {
    /// Standard Monte Carlo sampling via Philox PRNG.
    MonteCarlo,
    /// Quasi-Monte Carlo via low-discrepancy Sobol sequence (dimensions 1..=10).
    QuasiMonteCarlo,
    /// Polynomial Chaos Expansion with orthogonal polynomial regression.
    PolynomialChaos,
    /// Multilevel Monte Carlo across hierarchical discretization ladders.
    MultilevelMonteCarlo,
    /// Epistemic interval bounding (vertex / box evaluation).
    EpistemicBounding,
}

/// Status of the uncertainty propagation execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UqStatus {
    /// Completed all planned samples; no statistical stopping bound is claimed.
    Complete,
    /// Stopped honestly due to sample or time budget exhaustion.
    BudgetTruncated,
    /// Observed cancellation request before completion.
    Cancelled,
    /// Refused invalid input (e.g. non-PSD correlation, unknown dependence).
    Refused,
}

impl UqStatus {
    /// Machine-readable status label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::BudgetTruncated => "budget-truncated",
            Self::Cancelled => "cancelled",
            Self::Refused => "refused",
        }
    }
}

/// Specification plan for executing uncertainty propagation.
#[derive(Debug, Clone, PartialEq)]
pub struct UqPlan {
    /// Target physical Quantity of Interest.
    pub target_qoi: String,
    /// Optional upper compliance limit (QoI <= threshold).
    pub compliance_threshold: Option<f64>,
    /// Uncertain input parameters.
    pub parameters: Vec<ParameterUncertainty>,
    /// Parameter dependence / correlation model.
    pub correlation: CorrelationModel,
    /// Propagation algorithm to employ.
    pub method: PropagationMethod,
    /// Maximum sample evaluation budget.
    pub budget_max_samples: usize,
    /// Pseudo-random number generator seed.
    pub seed: u64,
}

impl UqPlan {
    /// Create a new UQ plan for a target QoI. Multivariate dependence must be
    /// declared explicitly; constructing a plan does not assume independence.
    #[must_use]
    pub fn new(
        target_qoi: impl Into<String>,
        method: PropagationMethod,
        max_samples: usize,
    ) -> Self {
        Self {
            target_qoi: target_qoi.into(),
            compliance_threshold: None,
            parameters: Vec::new(),
            correlation: CorrelationModel::Unknown,
            method,
            budget_max_samples: max_samples,
            seed: 0x0517,
        }
    }

    /// Add an uncertain parameter.
    #[must_use]
    pub fn with_parameter(mut self, param: ParameterUncertainty) -> Self {
        self.parameters.push(param);
        self
    }

    /// Set the compliance threshold (QoI <= threshold).
    #[must_use]
    pub fn with_compliance_threshold(mut self, threshold: f64) -> Self {
        self.compliance_threshold = Some(threshold);
        self
    }

    /// Set the correlation model.
    #[must_use]
    pub fn with_correlation(mut self, correlation: CorrelationModel) -> Self {
        self.correlation = correlation;
        self
    }
}

/// Structured outcome of uncertainty propagation.
#[derive(Debug, Clone, PartialEq)]
pub struct UqResult {
    /// Evaluated Quantity of Interest name.
    pub qoi_name: String,
    /// Algorithm used for propagation.
    pub method_used: PropagationMethod,
    /// Number of evaluator calls, including a terminal rejected evaluation.
    pub samples_evaluated: usize,
    /// Sample mean, absent when there are no observations or the run refused.
    pub mean: Option<f64>,
    /// Sample standard deviation, absent when fewer than two observations exist
    /// or a partial-run dispersion cannot be represented as a finite f64.
    pub std_dev: Option<f64>,
    /// Quantile percentiles [p05, p50, p95].
    pub percentiles: Option<[f64; 3]>,
    /// Empirical min/max; [0, 0] is a placeholder when statistics are absent.
    pub interval_bounds: [f64; 2],
    /// Empirical probability of meeting compliance ceiling.
    pub probability_of_compliance: Option<f64>,
    /// Descriptive standard error; unavailable (numeric placeholder 0) when
    /// `std_dev` is None. Not an optional-stopping confidence bound.
    pub sampling_error: f64,
    /// Assigned evidence color classification.
    pub evidence_color: Color,
    /// Final execution status.
    pub status: UqStatus,
    /// Detailed reason if propagation refused or failed.
    pub rejection_reason: Option<String>,
}

impl UqResult {
    /// Generate a deterministic BLAKE3 digest of the uncertainty evidence.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        // All result fields, including authority and refusal, participate. Debug
        // formatting escapes strings and preserves round-trippable f64 values;
        // the old six-decimal range omitted statistically meaningful changes.
        hash_domain(
            "org.frankensim.uq.result.v2",
            format!("{self:?}").as_bytes(),
        )
    }
}

/// Evaluator and runner for product uncertainty propagation plans.
pub struct UqPropagator;

impl UqPropagator {
    /// Execute the entire admitted plan. For fallible solvers, cancellation,
    /// bounded work chunks, and resumption use [`crate::UqExecution`].
    pub fn run<F: Fn(&[f64]) -> f64>(plan: &UqPlan, evaluator: F) -> UqResult {
        match crate::UqExecution::new(plan) {
            Ok(mut execution) => execution.advance(
                plan.budget_max_samples,
                || false,
                |parameters| Ok::<f64, core::convert::Infallible>(evaluator(parameters)),
            ),
            Err(reason) => refused(plan, reason, 0),
        }
    }
}

pub(crate) fn refused(plan: &UqPlan, reason: impl Into<String>, evaluated: usize) -> UqResult {
    UqResult {
        qoi_name: plan.target_qoi.clone(),
        method_used: plan.method,
        samples_evaluated: evaluated,
        mean: None,
        std_dev: None,
        percentiles: None,
        interval_bounds: [0.0, 0.0],
        probability_of_compliance: None,
        sampling_error: 0.0,
        evidence_color: Color::Estimated {
            estimator: "refused; no uncertainty result".into(),
            dispersion: 0.0,
        },
        status: UqStatus::Refused,
        rejection_reason: Some(reason.into()),
    }
}

/// Returns a numerical PSD factor for an explicitly joint Gaussian model.
/// Correlation alone does not specify a non-Gaussian copula, so those models
/// refuse. This is floating-point admission, not an interval PSD certificate.
pub(crate) fn admit_plan(plan: &UqPlan) -> Result<Option<Vec<Vec<f64>>>, &'static str> {
    let dim = plan.parameters.len();
    if dim == 0 || dim > 256 || !(2..=1_000_000).contains(&plan.budget_max_samples) {
        return Err("supported envelope: 1..=256 parameters and 2..=1000000 samples");
    }
    if plan.method != PropagationMethod::MonteCarlo {
        return Err(
            "this driver implements MonteCarlo only; use an admitted method-specific producer",
        );
    }
    if plan.target_qoi.is_empty() || plan.compliance_threshold.is_some_and(|x| !x.is_finite()) {
        return Err("QoI name and finite compliance threshold are required");
    }
    let mut names = std::collections::BTreeSet::new();
    for p in &plan.parameters {
        if p.name.is_empty() || p.unit.is_empty() || !names.insert(&p.name) {
            return Err("parameters require distinct nonempty names and explicit units");
        }
        match p.kind {
            UncertaintyKind::AleatoryGaussian { mean, std_dev }
                if mean.is_finite() && std_dev.is_finite() && std_dev >= 0.0 => {}
            UncertaintyKind::AleatoryUniform { lo, hi }
                if lo.is_finite() && hi.is_finite() && lo <= hi => {}
            _ => {
                return Err(
                    "declare a finite Gaussian or uniform probability measure; bands and unstated inputs are not distributions",
                );
            }
        }
    }
    let random_dimensions = plan
        .parameters
        .iter()
        .filter(|p| match p.kind {
            UncertaintyKind::AleatoryGaussian { std_dev, .. } => std_dev > 0.0,
            UncertaintyKind::AleatoryUniform { lo, hi } => lo < hi,
            _ => false,
        })
        .count();
    let matrix = match &plan.correlation {
        CorrelationModel::Unknown if random_dimensions > 1 => {
            return Err(
                "joint probability requires explicit dependence; independence is not a default",
            );
        }
        CorrelationModel::Unknown | CorrelationModel::Independent => return Ok(None),
        CorrelationModel::Correlated { .. } => {
            return Err(
                "correlations alone do not specify a joint measure; declare JointGaussian or an implemented copula",
            );
        }
        CorrelationModel::JointGaussian { matrix } => matrix,
    };
    if plan
        .parameters
        .iter()
        .any(|p| !matches!(p.kind, UncertaintyKind::AleatoryGaussian { .. }))
    {
        return Err(
            "correlated sampling requires joint Gaussian marginals; a correlation matrix does not specify a general copula",
        );
    }
    if matrix.len() != dim || matrix.iter().any(|r| r.len() != dim) {
        return Err("correlation matrix dimensions must match parameter declaration order");
    }
    for i in 0..dim {
        for j in 0..dim {
            let x = matrix[i][j];
            if !x.is_finite() || x.abs() > 1.0 || x != matrix[j][i] || (i == j && x != 1.0) {
                return Err("correlation must be finite, symmetric, in [-1,1], with unit diagonal");
            }
        }
    }
    let tolerance = 64.0 * f64::EPSILON * dim as f64;
    let mut l = vec![vec![0.0; dim]; dim];
    for i in 0..dim {
        for j in 0..=i {
            let residual = matrix[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if residual < -tolerance {
                    return Err("correlation matrix is not positive semidefinite");
                }
                l[i][j] = residual.max(0.0).sqrt();
            } else if l[j][j] > 0.0 {
                l[i][j] = residual / l[j][j];
            } else if residual.abs() > tolerance {
                return Err("correlation matrix has an inconsistent singular pivot");
            }
        }
    }
    Ok(Some(l))
}
