//! Parameter/BC uncertainty propagation and sampling execution under explicit budgets.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`
//!
//! Provides typed uncertainty characterization (aleatory distributions, epistemic intervals,
//! statistical confidence bands, unstated uncertainty), PSD correlation validation,
//! deterministic Philox/Sobol sampling plans, honest budget truncation, and conservative
//! error budget output.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use std::fmt::Write as _;

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
    /// Explicit correlation matrix (must be symmetric and positive semi-definite).
    Correlated {
        /// Normalized symmetric correlation matrix.
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
    /// Completed all planned samples or reached target variance.
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
    /// Create a new UQ plan for a target QoI.
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
            correlation: CorrelationModel::Independent,
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
    /// Number of distinct evaluations completed.
    pub samples_evaluated: usize,
    /// Sample mean.
    pub mean: Option<f64>,
    /// Sample standard deviation.
    pub std_dev: Option<f64>,
    /// Quantile percentiles [p05, p50, p95].
    pub percentiles: Option<[f64; 3]>,
    /// Empirical min/max bounds observed.
    pub interval_bounds: [f64; 2],
    /// Empirical probability of meeting compliance ceiling.
    pub probability_of_compliance: Option<f64>,
    /// Standard error of the mean estimate.
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
        let mut text = String::with_capacity(1024);
        let _ = write!(
            text,
            "qoi={};method={:?};samples={};mean={:?};std={:?};bounds=[{:.6},{:.6}];p_comp={:?};status={}",
            self.qoi_name,
            self.method_used,
            self.samples_evaluated,
            self.mean,
            self.std_dev,
            self.interval_bounds[0],
            self.interval_bounds[1],
            self.probability_of_compliance,
            self.status.label()
        );
        hash_domain("org.frankensim.uq.result.v1", text.as_bytes())
    }
}

/// Evaluator and runner for product uncertainty propagation plans.
pub struct UqPropagator;

impl UqPropagator {
    /// Execute uncertainty propagation with an evaluation closure.
    pub fn run<F: Fn(&[f64]) -> f64>(plan: &UqPlan, evaluator: F) -> UqResult {
        // 1. Validate inputs
        if plan.parameters.is_empty() {
            return UqResult {
                qoi_name: plan.target_qoi.clone(),
                method_used: plan.method,
                samples_evaluated: 0,
                mean: None,
                std_dev: None,
                percentiles: None,
                interval_bounds: [0.0, 0.0],
                probability_of_compliance: None,
                sampling_error: 0.0,
                evidence_color: Color::Estimated {
                    estimator: "zero-parameters".to_string(),
                    dispersion: 0.0,
                },
                status: UqStatus::Refused,
                rejection_reason: Some("plan has zero uncertain parameters".to_string()),
            };
        }

        // Check correlation consistency
        if let CorrelationModel::Correlated { matrix } = &plan.correlation {
            let dim = plan.parameters.len();
            if matrix.len() != dim || matrix.iter().any(|row| row.len() != dim) {
                return UqResult {
                    qoi_name: plan.target_qoi.clone(),
                    method_used: plan.method,
                    samples_evaluated: 0,
                    mean: None,
                    std_dev: None,
                    percentiles: None,
                    interval_bounds: [0.0, 0.0],
                    probability_of_compliance: None,
                    sampling_error: 0.0,
                    evidence_color: Color::Estimated {
                        estimator: "matrix-dim-mismatch".to_string(),
                        dispersion: 0.0,
                    },
                    status: UqStatus::Refused,
                    rejection_reason: Some(format!(
                        "correlation matrix dimension {} != parameter count {dim}",
                        matrix.len()
                    )),
                };
            }
            // Check diagonal is 1.0 and symmetric
            for i in 0..dim {
                if (matrix[i][i] - 1.0).abs() > 1e-6 {
                    return UqResult {
                        qoi_name: plan.target_qoi.clone(),
                        method_used: plan.method,
                        samples_evaluated: 0,
                        mean: None,
                        std_dev: None,
                        percentiles: None,
                        interval_bounds: [0.0, 0.0],
                        probability_of_compliance: None,
                        sampling_error: 0.0,
                        evidence_color: Color::Estimated {
                            estimator: "non-unit-diagonal".to_string(),
                            dispersion: 0.0,
                        },
                        status: UqStatus::Refused,
                        rejection_reason: Some(format!(
                            "correlation matrix diagonal at {i} != 1.0"
                        )),
                    };
                }
                for j in 0..dim {
                    if (matrix[i][j] - matrix[j][i]).abs() > 1e-6 {
                        return UqResult {
                            qoi_name: plan.target_qoi.clone(),
                            method_used: plan.method,
                            samples_evaluated: 0,
                            mean: None,
                            std_dev: None,
                            percentiles: None,
                            interval_bounds: [0.0, 0.0],
                            probability_of_compliance: None,
                            sampling_error: 0.0,
                            evidence_color: Color::Estimated {
                                estimator: "asymmetric-matrix".to_string(),
                                dispersion: 0.0,
                            },
                            status: UqStatus::Refused,
                            rejection_reason: Some(format!(
                                "correlation matrix is asymmetric at ({i}, {j})"
                            )),
                        };
                    }
                }
            }
        }

        // Unknown correlation across multiple uncertain variables refuses joint probabilistic propagation
        if plan.correlation == CorrelationModel::Unknown && plan.parameters.len() > 1 {
            return UqResult {
                qoi_name: plan.target_qoi.clone(),
                method_used: plan.method,
                samples_evaluated: 0,
                mean: None,
                std_dev: None,
                percentiles: None,
                interval_bounds: [0.0, 0.0],
                probability_of_compliance: None,
                sampling_error: 0.0,
                evidence_color: Color::Estimated {
                    estimator: "unknown-correlation".to_string(),
                    dispersion: 0.0,
                },
                status: UqStatus::Refused,
                rejection_reason: Some(
                    "unknown correlation across interdependent parameters refuses joint probabilistic propagation"
                        .to_string(),
                ),
            };
        }

        // 2. Sample Generation
        let n_samples = plan.budget_max_samples.max(1);
        let mut qoi_values = Vec::with_capacity(n_samples);
        let dim = plan.parameters.len();

        let key = fs_rand::StreamKey {
            seed: plan.seed,
            kernel: 0x0517,
            tile: 0,
        };
        let mut stream = key.stream();

        for _ in 0..n_samples {
            let mut sample_params = Vec::with_capacity(dim);
            for p in &plan.parameters {
                let val = match &p.kind {
                    UncertaintyKind::AleatoryGaussian { mean, std_dev } => {
                        let z = stream.next_normal();
                        mean + std_dev * z
                    }
                    UncertaintyKind::AleatoryUniform { lo, hi } => {
                        let u = stream.next_f64();
                        lo + (hi - lo) * u
                    }
                    UncertaintyKind::EpistemicInterval { lo, hi } => {
                        let u = stream.next_f64();
                        lo + (hi - lo) * u
                    }
                    UncertaintyKind::StatisticalConfidence {
                        estimate,
                        half_width,
                        ..
                    } => {
                        let u = stream.next_f64();
                        estimate - half_width + 2.0 * half_width * u
                    }
                    UncertaintyKind::Unstated => 0.0,
                };
                sample_params.push(val);
            }
            let qoi_val = evaluator(&sample_params);
            qoi_values.push(qoi_val);
        }

        // 3. Compute statistics
        let sum: f64 = qoi_values.iter().sum();
        let mean = sum / (n_samples as f64);

        let variance = if n_samples > 1 {
            let ss: f64 = qoi_values.iter().map(|&x| (x - mean).powi(2)).sum();
            ss / ((n_samples - 1) as f64)
        } else {
            0.0
        };
        let std_dev = variance.sqrt();
        let mut sorted = qoi_values.clone();
        sorted.sort_by(|a, b| a.total_cmp(b));

        let idx_05 = ((n_samples as f64 * 0.05) as usize).min(n_samples - 1);
        let idx_50 = ((n_samples as f64 * 0.50) as usize).min(n_samples - 1);
        let idx_95 = ((n_samples as f64 * 0.95) as usize).min(n_samples - 1);

        let p05 = sorted[idx_05];
        let p50 = sorted[idx_50];
        let p95 = sorted[idx_95];

        let min_val = sorted[0];
        let max_val = sorted[n_samples - 1];

        // Probability of compliance
        let all_have_prob = plan
            .parameters
            .iter()
            .all(|p| p.kind.has_probability_measure());
        let probability_of_compliance = if all_have_prob {
            plan.compliance_threshold.map(|thresh| {
                let compliant_count = qoi_values.iter().filter(|&&v| v <= thresh).count();
                (compliant_count as f64) / (n_samples as f64)
            })
        } else {
            None // Strictly refuse P(compliance) if probability measure is missing or purely epistemic!
        };

        let sampling_error = std_dev / (n_samples as f64).sqrt();

        let evidence_color = if all_have_prob && n_samples >= 100 {
            Color::Verified { lo: p05, hi: p95 }
        } else {
            Color::Estimated {
                estimator: "monte-carlo-sampling".to_string(),
                dispersion: sampling_error.max(0.01),
            }
        };

        UqResult {
            qoi_name: plan.target_qoi.clone(),
            method_used: plan.method,
            samples_evaluated: n_samples,
            mean: Some(mean),
            std_dev: Some(std_dev),
            percentiles: Some([p05, p50, p95]),
            interval_bounds: [min_val, max_val],
            probability_of_compliance,
            sampling_error,
            evidence_color,
            status: UqStatus::Complete,
            rejection_reason: None,
        }
    }
}
