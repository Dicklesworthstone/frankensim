//! fs-toleralloc — adjoint-driven tolerance allocation (plan addendum,
//! Proposal 11's commercial kicker). Layer: L4.
//!
//! GD&T (geometric dimensioning and tolerancing) today is assigned by
//! convention and fear. This replaces fear with a CERTIFIED SENSITIVITY: spend
//! tight manufacturing tolerances ONLY where `∂QoI/∂geometry` is large, and
//! provably LOOSEN everywhere else — delivered as savings a CFO understands.
//!
//! The allocation minimizes manufacturing cost subject to
//! `P(performance ∈ spec) ≥ target`, propagated FIRST-ORDER: the QoI variance
//! from independent feature tolerances is `Σ sᵢ² σᵢ²` with `σᵢ = tᵢ / k`. The
//! cost-optimal solution (Lagrange) allocates `tᵢ ∝ (cᵢ / sᵢ²)^{1/3}` — LOOSE
//! where sensitivity is small, TIGHT where it is large — normalized so the
//! variance budget is exactly met.
//!
//! First-order propagation is a LINEARIZATION, so [`robustness_check`] compares
//! it against the QoI evaluated at sampled tolerance-band EXTREMES and flags
//! where the linearization fails. Every loosened tolerance in the
//! [`gdt_report`] carries the certified sensitivity (with its color) that
//! justifies it. Deterministic; depends only on `fs-evidence` and the
//! `fs-math` deterministic scalar kernels.
//!
//! The additive correlated-stack lane admits a bounded lower-triangular
//! correlation factor `L` with binary64-near-unit rows, so `C = L Lᵀ` is
//! positive semidefinite by construction. [`propagate_correlated_stack`]
//! evaluates the signed first-order variance `aᵀ C a`, where
//! `aᵢ = sensitivityᵢ · σᵢ`, and retains the exact external model identity and
//! caller-supplied positional terms in its receipt.
//!
//! DETERMINISM DOCTRINE (bead frankensim-lyms): every transcendental in
//! this crate routes through `fs_math::det` so the "fully deterministic"
//! contract holds cross-ISA by construction — platform libm `ln`/`exp`
//! differ by ≥1 ULP across ISAs and libm versions. `sqrt` stays primitive
//! (IEEE-754 requires correct rounding for it).

use std::{collections::BTreeMap, num::NonZeroU64};

use fs_math::det;

pub use dependency_allocation::{
    DependencyGroup, DependencyGroupId, DependencyGroupReceipt, GROUPED_ALLOCATION_SCHEMA_V1,
    GroupedAllocationError, GroupedAllocationItem, GroupedAllocationReceipt,
    GroupedAllocationResource, GroupedDependenceIdentity, GroupedDependenceModel,
    GroupedDerivedQuantity, GroupedFeature, MAX_DEPENDENCY_GROUPS_V1, MAX_DEPENDENCY_NAME_BYTES_V1,
    MAX_DEPENDENCY_NAMESPACE_BYTES_V1, MAX_GROUPED_FEATURES_V1, allocate_grouped,
};
pub use fs_evidence::ColorRank;
pub use gear_backlash::{
    GEAR_BACKLASH_CONSUMER_SCHEMA_V1, GearBacklashConsumerDraftV1, GearBacklashConsumptionErrorV1,
    GearBacklashConsumptionReceiptV1, GearBacklashLengthUnitV1, GearBacklashModeShareV1,
    GearBacklashNumericIssueV1, GearBacklashNumericQuantityV1, GearBacklashProbabilityErrorV1,
    GearBacklashProbabilityV1, GearBacklashQuantileV1, GearBacklashSupportPointV1,
    GearBacklashValueV1, GearBacklashVarianceV1, MAX_GEAR_BACKLASH_QUANTILES_V1,
    consume_structured_backlash,
};
pub use structured_propagation::{
    ClampDisposition, InteriorKnotOwner, MAX_EXACT_STRUCTURED_WEIGHT_V1, MAX_STRUCTURED_DEPTH_V1,
    MAX_STRUCTURED_KEY_BYTES_V1, MAX_STRUCTURED_LAWS_V1, MAX_STRUCTURED_LEAVES_V1,
    MAX_STRUCTURED_NAMESPACE_BYTES_V1, MAX_STRUCTURED_NODES_V1, MAX_STRUCTURED_PIECES_PER_LAW_V1,
    PiecewiseQuadraticLaw, QuadraticResponsePiece, STRUCTURED_PROPAGATION_SCHEMA_V1,
    StructuredCoefficient, StructuredEvaluationStage, StructuredKeyRole, StructuredLawId,
    StructuredLawIssue, StructuredLeafReceipt, StructuredModeMomentReceipt,
    StructuredModelIdentity, StructuredMomentQuantity, StructuredMomentScope, StructuredNodeId,
    StructuredNodeMomentReceipt, StructuredNodeSpec, StructuredNumericIssue,
    StructuredPopulationModel, StructuredPropagationError, StructuredPropagationReceipt,
    StructuredResource, StructuredScalarLocation, StructuredTopologyIssue,
    propagate_structured_population,
};

mod dependency_allocation;
mod gear_backlash;
mod structured_propagation;

/// Maximum axis count admitted by the version-one correlated-stack lane.
pub const MAX_CORRELATED_STACK_TERMS_V1: usize = 128;

/// Maximum byte length of a correlation-model namespace.
pub const MAX_CORRELATION_MODEL_NAMESPACE_BYTES_V1: usize = 256;

/// Maximum UTF-8 byte length of one correlated-stack term name.
pub const MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1: usize = 256;

/// A geometric feature whose tolerance is being allocated.
#[derive(Debug, Clone, PartialEq)]
pub struct Feature {
    /// A stable name.
    pub name: String,
    /// `|∂QoI/∂geometry|` at this feature (the certified sensitivity, > 0).
    pub sensitivity: f64,
    /// The color of that sensitivity (verified for an adjoint-derived one).
    pub sensitivity_color: ColorRank,
    /// The cost coefficient `cᵢ` (cost `≈ cᵢ / tolerance`; tighter is costlier).
    pub cost_coeff: f64,
    /// The baseline (convention-assigned) tolerance, for tighten/loosen labels.
    pub baseline_tolerance: f64,
}

/// Why one supplied correlation-factor entry was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrelationFactorIssue {
    /// NaN and infinities are not correlation coefficients.
    NonFinite,
    /// Negative zero is not the canonical encoding of an exact zero.
    NonCanonicalNegativeZero,
    /// A row-major entry above the lower triangle must be canonical `+0.0`.
    AboveDiagonalNonZero,
    /// Canonical lower-triangular factors use nonnegative diagonal entries.
    NegativeDiagonal,
    /// A row admitted near unit norm cannot contain magnitude above one.
    MagnitudeAboveOne,
}

/// Refusal from admitting one externally identified correlation factor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorrelationAdmissionError {
    /// The namespace is not a bounded canonical slash-separated key.
    InvalidNamespace {
        /// Bounded UTF-8 prefix of the rejected namespace.
        namespace: String,
        /// Stable grammar explanation.
        reason: &'static str,
    },
    /// An all-zero semantic digest cannot identify an external model.
    ZeroDigest,
    /// The declared axis count is empty or exceeds the versioned cap.
    InvalidDimension {
        /// Supplied axis count.
        dimension: usize,
        /// Versioned maximum.
        max: usize,
    },
    /// The row-major factor length is not exactly `dimension²`.
    FactorLength {
        /// Declared matrix dimension.
        dimension: usize,
        /// Required scalar count.
        expected: usize,
        /// Supplied scalar count.
        actual: usize,
    },
    /// One factor entry violates the canonical lower-triangular grammar.
    InvalidFactorEntry {
        /// Zero-based row.
        row: usize,
        /// Zero-based column.
        column: usize,
        /// Stable refusal class.
        issue: CorrelationFactorIssue,
    },
    /// One binary64-computed factor-row norm is not near one within the
    /// admitted deterministic roundoff envelope.
    NonUnitRow {
        /// Zero-based factor row.
        row: usize,
        /// Exact IEEE-754 bits of the binary64-computed squared row norm.
        norm_squared_bits: u64,
        /// Exact IEEE-754 bits of the admitted absolute defect.
        tolerance_bits: u64,
    },
}

/// A bounded positive-semidefinite correlation model admitted from its factor.
///
/// Admission proves the finite lower-triangular factor grammar and PSD
/// construction. Its binary64 row-norm check is not an exact-real enclosure of
/// the implied diagonal. The external owner retains population, process,
/// calibration, and model-form authority. The exact factor is
/// representation-semantic: an equivalent singular correlation matrix need
/// not have a unique factor, and this crate neither derives nor authenticates
/// the caller-supplied semantic digest.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedCorrelationModel {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    semantic_digest: [u8; 32],
    dimension: usize,
    lower_factor: Box<[f64]>,
    max_row_norm_defect: f64,
}

impl AdmittedCorrelationModel {
    /// Admit an externally identified lower-triangular factor `L` for the
    /// correlation matrix `C = L Lᵀ`.
    ///
    /// `lower_factor` is row-major and must contain canonical `+0.0` above the
    /// diagonal. Every binary64-computed row squared norm must be within a
    /// deterministic `64 · dimension · ε` admission envelope around one. This
    /// is a measured structural check, not an exact-real diagonal enclosure.
    ///
    /// # Errors
    ///
    /// Refuses an unstable identity, empty/oversized dimension, malformed
    /// factor layout, non-finite/noncanonical entries, or a row outside the
    /// measured near-unit envelope.
    pub fn try_new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        semantic_digest: [u8; 32],
        dimension: usize,
        lower_factor: Vec<f64>,
    ) -> Result<Self, CorrelationAdmissionError> {
        let namespace = namespace.into();
        validate_correlation_namespace(&namespace)?;
        if semantic_digest == [0; 32] {
            return Err(CorrelationAdmissionError::ZeroDigest);
        }
        if dimension == 0 || dimension > MAX_CORRELATED_STACK_TERMS_V1 {
            return Err(CorrelationAdmissionError::InvalidDimension {
                dimension,
                max: MAX_CORRELATED_STACK_TERMS_V1,
            });
        }
        let expected = dimension * dimension;
        if lower_factor.len() != expected {
            return Err(CorrelationAdmissionError::FactorLength {
                dimension,
                expected,
                actual: lower_factor.len(),
            });
        }

        let row_tolerance = 64.0 * f64::EPSILON * dimension as f64;
        let mut max_row_norm_defect = 0.0_f64;
        for row in 0..dimension {
            let mut norm_squared = 0.0_f64;
            for column in 0..dimension {
                let value = lower_factor[row * dimension + column];
                let issue = if !value.is_finite() {
                    Some(CorrelationFactorIssue::NonFinite)
                } else if value == 0.0 && value.is_sign_negative() {
                    Some(CorrelationFactorIssue::NonCanonicalNegativeZero)
                } else if column > row && value != 0.0 {
                    Some(CorrelationFactorIssue::AboveDiagonalNonZero)
                } else if column == row && value < 0.0 {
                    Some(CorrelationFactorIssue::NegativeDiagonal)
                } else if value.abs() > 1.0 {
                    Some(CorrelationFactorIssue::MagnitudeAboveOne)
                } else {
                    None
                };
                if let Some(issue) = issue {
                    return Err(CorrelationAdmissionError::InvalidFactorEntry {
                        row,
                        column,
                        issue,
                    });
                }
                if column <= row {
                    norm_squared += value * value;
                }
            }
            let defect = (norm_squared - 1.0).abs();
            if defect > row_tolerance {
                return Err(CorrelationAdmissionError::NonUnitRow {
                    row,
                    norm_squared_bits: norm_squared.to_bits(),
                    tolerance_bits: row_tolerance.to_bits(),
                });
            }
            max_row_norm_defect = max_row_norm_defect.max(defect);
        }

        Ok(Self {
            namespace: namespace.into_boxed_str(),
            schema_version,
            semantic_digest,
            dimension,
            lower_factor: lower_factor.into_boxed_slice(),
            max_row_norm_defect,
        })
    }

    /// External model namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Explicit external schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact semantic digest supplied by the external owner.
    #[must_use]
    pub const fn semantic_digest(&self) -> [u8; 32] {
        self.semantic_digest
    }

    /// Number of positional factor axes.
    ///
    /// External axis identifiers/order are not carried or authenticated by
    /// this seed type; callers bind terms positionally in the receipt.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Exact admitted row-major lower-triangular factor.
    #[must_use]
    pub fn lower_factor(&self) -> &[f64] {
        &self.lower_factor
    }

    /// Largest binary64-computed absolute defect in a factor row's squared norm.
    #[must_use]
    pub const fn max_row_norm_defect(&self) -> f64 {
        self.max_row_norm_defect
    }
}

/// One signed first-order term in a correlated manufacturing stack.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelatedStackTerm {
    /// Bounded caller label in positional factor order.
    pub name: String,
    /// Signed derivative `∂QoI/∂axis`; zero is permitted and remains explicit.
    pub signed_sensitivity: f64,
    /// Evidence color carried by the supplied sensitivity.
    pub sensitivity_color: ColorRank,
    /// Strictly positive standard deviation for this manufacturing axis.
    pub standard_deviation: f64,
}

/// A correlated-stack derived quantity that could not be represented safely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrelatedDerivedQuantity {
    /// One signed `sensitivity · standard_deviation` term.
    ScaledSensitivity,
    /// One nonzero scaled term normalized by the stack scale.
    NormalizedSensitivity,
    /// One nonzero factor-times-normalized-term product.
    CorrelationProjectionProduct,
    /// One factor-column projection whose numerical zero was ambiguous.
    CorrelationProjection,
    /// The standard deviation under independent axes.
    IndependentStandardDeviation,
    /// The variance under independent axes.
    IndependentVariance,
    /// The standard deviation under the admitted correlation factor.
    CorrelatedStandardDeviation,
    /// The variance under the admitted correlation factor.
    CorrelatedVariance,
}

/// Refusal from evaluating one admitted correlation model and term stack.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorrelatedStackError {
    /// A stack with no axes has no manufacturing semantics.
    NoTerms,
    /// The supplied stack exceeds the versioned work/retention cap.
    TooManyTerms {
        /// Supplied term count.
        actual: usize,
        /// Versioned maximum.
        max: usize,
    },
    /// Model and term axis counts differ.
    DimensionMismatch {
        /// Model axis count.
        model: usize,
        /// Supplied term count.
        terms: usize,
    },
    /// One term name is empty or unstable.
    InvalidTermName {
        /// Position in model order.
        index: usize,
        /// Bounded UTF-8 prefix of the rejected spelling.
        name: String,
        /// Stable explanation.
        reason: &'static str,
    },
    /// Two term names collide under deterministic lowercase comparison.
    AmbiguousTermName {
        /// First position.
        first_index: usize,
        /// Colliding position.
        duplicate_index: usize,
        /// Canonical comparison key.
        canonical_name: String,
    },
    /// One term scalar is outside its declared domain.
    InvalidTermField {
        /// Position in model order.
        index: usize,
        /// Term name.
        name: String,
        /// Rejected field.
        field: &'static str,
        /// Domain violation.
        issue: ScalarIssue,
    },
    /// Finite admitted inputs produced an unrepresentable result.
    InvalidDerived {
        /// Failed quantity.
        quantity: CorrelatedDerivedQuantity,
        /// Term position when the failure belongs to one axis.
        term_index: Option<usize>,
        /// Numeric failure class.
        issue: ScalarIssue,
    },
}

/// Non-forgeable result of one correlated first-order stack evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelatedStackReceipt {
    model: AdmittedCorrelationModel,
    terms: Box<[CorrelatedStackTerm]>,
    independent_standard_deviation: f64,
    independent_variance: f64,
    correlated_standard_deviation: f64,
    correlated_variance: f64,
    correlation_variance_delta: f64,
}

impl CorrelatedStackReceipt {
    /// Exact admitted correlation model used by the evaluation.
    #[must_use]
    pub const fn model(&self) -> &AdmittedCorrelationModel {
        &self.model
    }

    /// Exact caller-supplied terms in positional factor order.
    #[must_use]
    pub fn terms(&self) -> &[CorrelatedStackTerm] {
        &self.terms
    }

    /// First-order standard deviation under an independence assumption.
    #[must_use]
    pub const fn independent_standard_deviation(&self) -> f64 {
        self.independent_standard_deviation
    }

    /// First-order variance under an independence assumption.
    #[must_use]
    pub const fn independent_variance(&self) -> f64 {
        self.independent_variance
    }

    /// First-order standard deviation under the admitted factor.
    #[must_use]
    pub const fn correlated_standard_deviation(&self) -> f64 {
        self.correlated_standard_deviation
    }

    /// First-order variance under the admitted factor.
    #[must_use]
    pub const fn correlated_variance(&self) -> f64 {
        self.correlated_variance
    }

    /// Signed binary64 `correlated_variance - independent_variance`;
    /// correlation may increase or decrease the propagated variance. A zero
    /// delta does not certify independence or absence of an exact-real effect.
    #[must_use]
    pub const fn correlation_variance_delta(&self) -> f64 {
        self.correlation_variance_delta
    }
}

/// Why a scalar was rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarIssue {
    /// The value was NaN or infinite.
    NonFinite,
    /// The value was zero or negative where strict positivity is required.
    NonPositive,
    /// The value was negative where zero is permitted.
    Negative,
    /// The value was outside the open unit interval `(0, 1)`.
    OutsideOpenUnitInterval,
    /// A semantic exact zero used the noncanonical negative-zero encoding.
    NonCanonicalNegativeZero,
    /// A nonzero mathematical result rounded to zero in binary64.
    Underflow,
    /// Nonzero inputs produced a numerical zero whose exactness is unknown.
    AmbiguousZero,
}

/// A derived quantity that could not be represented safely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedQuantity {
    /// The logarithmic normalization for the allocation.
    AllocationNormalization,
    /// One feature's allocated tolerance.
    Tolerance,
    /// One feature's manufacturing-cost contribution.
    CostContribution,
    /// One feature's QoI-variance contribution.
    VarianceContribution,
    /// The accumulated manufacturing cost.
    TotalCost,
    /// The accumulated QoI variance.
    AchievedVariance,
    /// The linearized standard deviation.
    LinearizedStandardDeviation,
    /// A sampled absolute deviation from the nominal QoI.
    SampledDeviation,
    /// The admissible sampled-extreme bound.
    RobustnessBound,
    /// The normal quantile used to derive a variance budget.
    NormalQuantile,
    /// The variance budget derived from a probability target.
    VarianceBudget,
}

/// A structured allocation or robustness failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ToleranceError {
    /// No features.
    NoFeatures,
    /// A feature name is empty or not a stable canonical spelling.
    InvalidFeatureName {
        /// Position in the input feature slice.
        index: usize,
        /// The rejected spelling.
        name: String,
        /// Stable explanation of the naming violation.
        reason: &'static str,
    },
    /// Two names collapse to the same canonical comparison key.
    AmbiguousFeatureName {
        /// Position of the first spelling.
        first_index: usize,
        /// Position of the colliding spelling.
        duplicate_index: usize,
        /// Locale-independent lowercase comparison key.
        canonical_name: String,
    },
    /// A feature scalar is outside its declared domain.
    InvalidFeatureField {
        /// Position in the input feature slice.
        index: usize,
        /// The offending feature name.
        feature: String,
        /// Field that was rejected.
        field: &'static str,
        /// Domain violation.
        issue: ScalarIssue,
    },
    /// A caller-supplied allocation item is unsafe to publish in a report.
    InvalidAllocationItem {
        /// Position in `Allocation::items`.
        index: usize,
        /// Item name.
        name: String,
        /// Field that was rejected.
        field: &'static str,
        /// Domain violation.
        issue: ScalarIssue,
    },
    /// A scalar API argument is outside its declared domain.
    InvalidArgument {
        /// Argument that was rejected.
        argument: &'static str,
        /// Domain violation.
        issue: ScalarIssue,
    },
    /// Finite admitted inputs produced an unrepresentable result.
    InvalidDerived {
        /// Quantity that failed.
        quantity: DerivedQuantity,
        /// Feature position, when the quantity belongs to one feature.
        feature_index: Option<usize>,
        /// Domain violation.
        issue: ScalarIssue,
    },
    /// A robustness claim requires at least one sampled band extreme.
    NoExtremeSamples,
    /// One sampled extreme is non-finite.
    InvalidExtremeQoi {
        /// Position in `extreme_qois`.
        index: usize,
        /// Domain violation.
        issue: ScalarIssue,
    },
    /// Arithmetic on one finite sampled extreme was unrepresentable.
    InvalidExtremeDerived {
        /// Position in `extreme_qois`.
        index: usize,
        /// Quantity that failed.
        quantity: DerivedQuantity,
        /// Domain violation.
        issue: ScalarIssue,
    },
}

/// What the allocator did to a feature's tolerance relative to its baseline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// Tolerance reduced (high sensitivity).
    Tighten,
    /// Tolerance widened (low sensitivity) — the savings.
    Loosen,
    /// Unchanged within rounding.
    Unchanged,
}

/// One feature's allocated tolerance.
#[derive(Debug, Clone, PartialEq)]
pub struct TolItem {
    /// The feature name.
    pub name: String,
    /// The allocated manufacturing tolerance.
    pub tolerance: f64,
    /// The certified sensitivity that justified it.
    pub sensitivity: f64,
    /// The sensitivity's color.
    pub sensitivity_color: ColorRank,
    /// Tighten / loosen / unchanged vs the baseline.
    pub action: Action,
}

/// The result of a tolerance allocation.
#[derive(Debug, Clone, PartialEq)]
pub struct Allocation {
    /// Per-feature allocation.
    pub items: Vec<TolItem>,
    /// Total manufacturing cost `Σ cᵢ / tᵢ` (lower is cheaper).
    pub total_cost: f64,
    /// The achieved QoI variance (== the budget, by construction).
    pub achieved_variance: f64,
}

/// Allocate cost-optimal tolerances that meet a QoI variance budget. `k` is the
/// tolerance-to-σ factor (`σ = t / k`, e.g. `k = 3` for a 3σ band).
///
/// # Errors
/// [`ToleranceError`] on empty input, ambiguous names, invalid feature fields,
/// invalid budget / `k`, or unrepresentable derived outputs.
pub fn allocate(
    features: &[Feature],
    variance_budget: f64,
    k: f64,
) -> Result<Allocation, ToleranceError> {
    if features.is_empty() {
        return Err(ToleranceError::NoFeatures);
    }
    validate_positive_argument("variance_budget", variance_budget)?;
    validate_positive_argument("k", k)?;

    let mut canonical_names = BTreeMap::new();
    for (index, feature) in features.iter().enumerate() {
        let canonical_name = canonical_feature_name(index, &feature.name)?;
        if let Some(&first_index) = canonical_names.get(&canonical_name) {
            return Err(ToleranceError::AmbiguousFeatureName {
                first_index,
                duplicate_index: index,
                canonical_name,
            });
        }
        canonical_names.insert(canonical_name, index);
        validate_positive_feature(index, feature, "sensitivity", feature.sensitivity)?;
        validate_positive_feature(index, feature, "cost_coeff", feature.cost_coeff)?;
        validate_positive_feature(
            index,
            feature,
            "baseline_tolerance",
            feature.baseline_tolerance,
        )?;
    }

    // Work in log space so finite, positive values do not overflow merely from
    // squaring a sensitivity or k. The public tolerance is still refused if
    // its mathematically required value is not representable as a positive
    // finite f64.
    let log_k = det::ln(k);
    let log_shapes: Vec<f64> = features
        .iter()
        .map(|feature| (det::ln(feature.cost_coeff) - 2.0 * det::ln(feature.sensitivity)) / 3.0)
        .collect();
    let log_variance_terms: Vec<f64> = features
        .iter()
        .zip(&log_shapes)
        .map(|(feature, &log_shape)| 2.0 * (det::ln(feature.sensitivity) - log_k + log_shape))
        .collect();
    if log_shapes.iter().any(|value| !value.is_finite())
        || log_variance_terms.iter().any(|value| !value.is_finite())
    {
        return Err(invalid_derived(
            DerivedQuantity::AllocationNormalization,
            None,
            ScalarIssue::NonFinite,
        ));
    }
    let max_log_variance = log_variance_terms
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let scaled_variance_sum: f64 = log_variance_terms
        .iter()
        .map(|term| det::exp(term - max_log_variance))
        .sum();
    let log_shape_variance = max_log_variance + det::ln(scaled_variance_sum);
    let log_scale = 0.5 * (det::ln(variance_budget) - log_shape_variance);
    if !log_scale.is_finite() {
        return Err(invalid_derived(
            DerivedQuantity::AllocationNormalization,
            None,
            ScalarIssue::NonFinite,
        ));
    }

    let mut items = Vec::with_capacity(features.len());
    let mut total_cost = 0.0;
    let mut achieved_variance = 0.0;
    for (index, (feature, &log_shape)) in features.iter().zip(&log_shapes).enumerate() {
        let log_tolerance = log_shape + log_scale;
        let tolerance = det::exp(log_tolerance);
        validate_positive_derived(DerivedQuantity::Tolerance, Some(index), tolerance)?;

        let cost_contribution = det::exp(det::ln(feature.cost_coeff) - log_tolerance);
        validate_positive_derived(
            DerivedQuantity::CostContribution,
            Some(index),
            cost_contribution,
        )?;
        total_cost += cost_contribution;
        validate_positive_derived(DerivedQuantity::TotalCost, None, total_cost)?;

        let log_variance = 2.0 * (det::ln(feature.sensitivity) - log_k + log_tolerance);
        let variance_contribution = det::exp(log_variance);
        validate_nonnegative_derived(
            DerivedQuantity::VarianceContribution,
            Some(index),
            variance_contribution,
        )?;
        achieved_variance += variance_contribution;
        validate_nonnegative_derived(DerivedQuantity::AchievedVariance, None, achieved_variance)?;

        let action = action_for(tolerance, feature.baseline_tolerance);
        items.push(TolItem {
            name: feature.name.clone(),
            tolerance,
            sensitivity: feature.sensitivity,
            sensitivity_color: feature.sensitivity_color,
            action,
        });
    }
    validate_positive_derived(DerivedQuantity::AchievedVariance, None, achieved_variance)?;
    Ok(Allocation {
        items,
        total_cost,
        achieved_variance,
    })
}

fn action_for(tolerance: f64, baseline: f64) -> Action {
    let log_ratio = det::ln(tolerance) - det::ln(baseline);
    if log_ratio > det::ln(1.01) {
        Action::Loosen
    } else if log_ratio < det::ln(0.99) {
        Action::Tighten
    } else {
        Action::Unchanged
    }
}

fn canonical_feature_name(index: usize, name: &str) -> Result<String, ToleranceError> {
    if name.is_empty() {
        return Err(ToleranceError::InvalidFeatureName {
            index,
            name: name.to_string(),
            reason: "name must not be empty",
        });
    }
    if name.trim() != name {
        return Err(ToleranceError::InvalidFeatureName {
            index,
            name: name.to_string(),
            reason: "name must not have leading or trailing whitespace",
        });
    }
    if name.chars().any(char::is_control) {
        return Err(ToleranceError::InvalidFeatureName {
            index,
            name: name.to_string(),
            reason: "name must not contain control characters",
        });
    }
    Ok(name.to_lowercase())
}

fn validate_positive_feature(
    index: usize,
    feature: &Feature,
    field: &'static str,
    value: f64,
) -> Result<(), ToleranceError> {
    let issue = if !value.is_finite() {
        Some(ScalarIssue::NonFinite)
    } else if value <= 0.0 {
        Some(ScalarIssue::NonPositive)
    } else {
        None
    };
    if let Some(issue) = issue {
        return Err(ToleranceError::InvalidFeatureField {
            index,
            feature: feature.name.clone(),
            field,
            issue,
        });
    }
    Ok(())
}

fn validate_positive_argument(argument: &'static str, value: f64) -> Result<(), ToleranceError> {
    let issue = if !value.is_finite() {
        Some(ScalarIssue::NonFinite)
    } else if value <= 0.0 {
        Some(ScalarIssue::NonPositive)
    } else {
        None
    };
    if let Some(issue) = issue {
        return Err(ToleranceError::InvalidArgument { argument, issue });
    }
    Ok(())
}

fn invalid_derived(
    quantity: DerivedQuantity,
    feature_index: Option<usize>,
    issue: ScalarIssue,
) -> ToleranceError {
    ToleranceError::InvalidDerived {
        quantity,
        feature_index,
        issue,
    }
}

fn validate_positive_derived(
    quantity: DerivedQuantity,
    feature_index: Option<usize>,
    value: f64,
) -> Result<(), ToleranceError> {
    if !value.is_finite() {
        Err(invalid_derived(
            quantity,
            feature_index,
            ScalarIssue::NonFinite,
        ))
    } else if value <= 0.0 {
        Err(invalid_derived(
            quantity,
            feature_index,
            ScalarIssue::NonPositive,
        ))
    } else {
        Ok(())
    }
}

fn validate_nonnegative_derived(
    quantity: DerivedQuantity,
    feature_index: Option<usize>,
    value: f64,
) -> Result<(), ToleranceError> {
    if !value.is_finite() {
        Err(invalid_derived(
            quantity,
            feature_index,
            ScalarIssue::NonFinite,
        ))
    } else if value < 0.0 {
        Err(invalid_derived(
            quantity,
            feature_index,
            ScalarIssue::Negative,
        ))
    } else {
        Ok(())
    }
}

fn validate_correlation_namespace(namespace: &str) -> Result<(), CorrelationAdmissionError> {
    let reason = if namespace.is_empty() {
        Some("namespace must not be empty")
    } else if namespace.len() > MAX_CORRELATION_MODEL_NAMESPACE_BYTES_V1 {
        Some("namespace exceeds the versioned byte cap")
    } else if namespace
        .as_bytes()
        .split(|byte| *byte == b'/')
        .any(|segment| segment.is_empty())
    {
        Some("namespace segments must not be empty")
    } else if namespace
        .as_bytes()
        .split(|byte| *byte == b'/')
        .any(|segment| {
            !segment[0].is_ascii_lowercase()
                || segment.last() == Some(&b'-')
                || segment.windows(2).any(|pair| pair == b"--")
                || segment.iter().any(|byte| {
                    !byte.is_ascii_lowercase() && !byte.is_ascii_digit() && *byte != b'-'
                })
        })
    {
        Some(
            "namespace segments must start lowercase and use lowercase ASCII letters, digits, or single interior hyphens",
        )
    } else {
        None
    };
    if let Some(reason) = reason {
        Err(CorrelationAdmissionError::InvalidNamespace {
            namespace: bounded_utf8_prefix(namespace, MAX_CORRELATION_MODEL_NAMESPACE_BYTES_V1),
            reason,
        })
    } else {
        Ok(())
    }
}

fn bounded_utf8_prefix(value: &str, max_bytes: usize) -> String {
    let mut end = value.len().min(max_bytes);
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    value[..end].to_string()
}

fn canonical_stack_term_name(index: usize, name: &str) -> Result<String, CorrelatedStackError> {
    let reason = if name.is_empty() {
        Some("name must not be empty")
    } else if name.len() > MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1 {
        Some("name exceeds the versioned byte cap")
    } else if name.trim() != name {
        Some("name must not have leading or trailing whitespace")
    } else if name.chars().any(char::is_control) {
        Some("name must not contain control characters")
    } else {
        None
    };
    if let Some(reason) = reason {
        Err(CorrelatedStackError::InvalidTermName {
            index,
            name: bounded_utf8_prefix(name, MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1),
            reason,
        })
    } else {
        let canonical_name = name.to_lowercase();
        if canonical_name.len() > MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1 {
            Err(CorrelatedStackError::InvalidTermName {
                index,
                name: name.to_string(),
                reason: "lowercase comparison key exceeds the versioned byte cap",
            })
        } else {
            Ok(canonical_name)
        }
    }
}

fn correlated_invalid_derived(
    quantity: CorrelatedDerivedQuantity,
    term_index: Option<usize>,
    issue: ScalarIssue,
) -> CorrelatedStackError {
    CorrelatedStackError::InvalidDerived {
        quantity,
        term_index,
        issue,
    }
}

fn normalized_l2(values: &[f64]) -> f64 {
    let mut scale = 0.0_f64;
    let mut sum_squares = 1.0_f64;
    for value in values.iter().copied().map(f64::abs) {
        if value == 0.0 {
            continue;
        }
        if scale < value {
            let ratio = scale / value;
            sum_squares = 1.0 + sum_squares * ratio * ratio;
            scale = value;
        } else {
            let ratio = value / scale;
            sum_squares += ratio * ratio;
        }
    }
    if scale == 0.0 {
        0.0
    } else {
        scale * sum_squares.sqrt()
    }
}

fn rescale_stack_deviation(
    scale: f64,
    normalized_standard_deviation: f64,
    standard_deviation_quantity: CorrelatedDerivedQuantity,
    variance_quantity: CorrelatedDerivedQuantity,
) -> Result<(f64, f64), CorrelatedStackError> {
    if normalized_standard_deviation == 0.0 {
        return Err(correlated_invalid_derived(
            standard_deviation_quantity,
            None,
            ScalarIssue::Underflow,
        ));
    }
    let standard_deviation = scale * normalized_standard_deviation;
    if !standard_deviation.is_finite() {
        return Err(correlated_invalid_derived(
            standard_deviation_quantity,
            None,
            ScalarIssue::NonFinite,
        ));
    }
    if standard_deviation == 0.0 {
        return Err(correlated_invalid_derived(
            standard_deviation_quantity,
            None,
            ScalarIssue::Underflow,
        ));
    }
    let variance = standard_deviation * standard_deviation;
    if !variance.is_finite() {
        return Err(correlated_invalid_derived(
            variance_quantity,
            None,
            ScalarIssue::NonFinite,
        ));
    }
    if variance == 0.0 {
        return Err(correlated_invalid_derived(
            variance_quantity,
            None,
            ScalarIssue::Underflow,
        ));
    }
    Ok((standard_deviation, variance))
}

/// Propagate signed first-order tolerance terms through an admitted
/// positive-semidefinite correlation factor.
///
/// For `aᵢ = signed_sensitivityᵢ · standard_deviationᵢ` and `C = L Lᵀ`, this
/// evaluates `variance = aᵀ C a = ||Lᵀ a||²`. The receipt also retains the
/// counterfactual independent-axis variance `Σ aᵢ²`, so an invalid independence
/// assumption is explicit rather than silently substituted.
///
/// # Errors
///
/// Refuses an empty/oversized or dimension-mismatched stack, unstable or
/// ambiguous term names, invalid scalars, and non-representable products,
/// standard deviations, or variances.
pub fn propagate_correlated_stack(
    model: &AdmittedCorrelationModel,
    terms: &[CorrelatedStackTerm],
) -> Result<CorrelatedStackReceipt, CorrelatedStackError> {
    if terms.is_empty() {
        return Err(CorrelatedStackError::NoTerms);
    }
    if terms.len() > MAX_CORRELATED_STACK_TERMS_V1 {
        return Err(CorrelatedStackError::TooManyTerms {
            actual: terms.len(),
            max: MAX_CORRELATED_STACK_TERMS_V1,
        });
    }
    if model.dimension != terms.len() {
        return Err(CorrelatedStackError::DimensionMismatch {
            model: model.dimension,
            terms: terms.len(),
        });
    }

    let mut canonical_names = BTreeMap::new();
    let mut scaled_terms = Vec::with_capacity(terms.len());
    let mut scale = 0.0_f64;
    for (index, term) in terms.iter().enumerate() {
        let canonical_name = canonical_stack_term_name(index, &term.name)?;
        if let Some(&first_index) = canonical_names.get(&canonical_name) {
            return Err(CorrelatedStackError::AmbiguousTermName {
                first_index,
                duplicate_index: index,
                canonical_name,
            });
        }
        canonical_names.insert(canonical_name, index);
        if !term.signed_sensitivity.is_finite() {
            return Err(CorrelatedStackError::InvalidTermField {
                index,
                name: term.name.clone(),
                field: "signed_sensitivity",
                issue: ScalarIssue::NonFinite,
            });
        }
        if term.signed_sensitivity == 0.0 && term.signed_sensitivity.is_sign_negative() {
            return Err(CorrelatedStackError::InvalidTermField {
                index,
                name: term.name.clone(),
                field: "signed_sensitivity",
                issue: ScalarIssue::NonCanonicalNegativeZero,
            });
        }
        let standard_deviation_issue = if !term.standard_deviation.is_finite() {
            Some(ScalarIssue::NonFinite)
        } else if term.standard_deviation <= 0.0 {
            Some(ScalarIssue::NonPositive)
        } else {
            None
        };
        if let Some(issue) = standard_deviation_issue {
            return Err(CorrelatedStackError::InvalidTermField {
                index,
                name: term.name.clone(),
                field: "standard_deviation",
                issue,
            });
        }
        let scaled = term.signed_sensitivity * term.standard_deviation;
        if !scaled.is_finite() {
            return Err(correlated_invalid_derived(
                CorrelatedDerivedQuantity::ScaledSensitivity,
                Some(index),
                ScalarIssue::NonFinite,
            ));
        }
        if scaled == 0.0 && term.signed_sensitivity != 0.0 {
            return Err(correlated_invalid_derived(
                CorrelatedDerivedQuantity::ScaledSensitivity,
                Some(index),
                ScalarIssue::Underflow,
            ));
        }
        scale = scale.max(scaled.abs());
        scaled_terms.push(scaled);
    }

    if scale == 0.0 {
        return Ok(CorrelatedStackReceipt {
            model: model.clone(),
            terms: terms.to_vec().into_boxed_slice(),
            independent_standard_deviation: 0.0,
            independent_variance: 0.0,
            correlated_standard_deviation: 0.0,
            correlated_variance: 0.0,
            correlation_variance_delta: 0.0,
        });
    }

    let mut normalized_terms = Vec::with_capacity(scaled_terms.len());
    for (index, scaled) in scaled_terms.iter().copied().enumerate() {
        let normalized = scaled / scale;
        if scaled != 0.0 && normalized == 0.0 {
            return Err(correlated_invalid_derived(
                CorrelatedDerivedQuantity::NormalizedSensitivity,
                Some(index),
                ScalarIssue::Underflow,
            ));
        }
        normalized_terms.push(normalized);
    }
    let independent_normalized_standard_deviation = normalized_l2(&normalized_terms);
    let (independent_standard_deviation, independent_variance) = rescale_stack_deviation(
        scale,
        independent_normalized_standard_deviation,
        CorrelatedDerivedQuantity::IndependentStandardDeviation,
        CorrelatedDerivedQuantity::IndependentVariance,
    )?;

    let mut normalized_projections = Vec::with_capacity(model.dimension);
    for column in 0..model.dimension {
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        let mut had_nonzero_product = false;
        for (row, normalized) in normalized_terms.iter().copied().enumerate().skip(column) {
            let factor = model.lower_factor[row * model.dimension + column];
            let value = factor * normalized;
            if factor != 0.0 && normalized != 0.0 {
                if value == 0.0 {
                    return Err(correlated_invalid_derived(
                        CorrelatedDerivedQuantity::CorrelationProjectionProduct,
                        Some(row),
                        ScalarIssue::Underflow,
                    ));
                }
                had_nonzero_product = true;
            }
            let next = sum + value;
            correction += if sum.abs() >= value.abs() {
                (sum - next) + value
            } else {
                (value - next) + sum
            };
            sum = next;
        }
        let projection = sum + correction;
        if !projection.is_finite() {
            return Err(correlated_invalid_derived(
                CorrelatedDerivedQuantity::CorrelationProjection,
                None,
                ScalarIssue::NonFinite,
            ));
        }
        if projection == 0.0 && had_nonzero_product {
            return Err(correlated_invalid_derived(
                CorrelatedDerivedQuantity::CorrelationProjection,
                None,
                ScalarIssue::AmbiguousZero,
            ));
        }
        normalized_projections.push(if projection == 0.0 { 0.0 } else { projection });
    }
    let correlated_normalized_standard_deviation = normalized_l2(&normalized_projections);
    let (correlated_standard_deviation, correlated_variance) = rescale_stack_deviation(
        scale,
        correlated_normalized_standard_deviation,
        CorrelatedDerivedQuantity::CorrelatedStandardDeviation,
        CorrelatedDerivedQuantity::CorrelatedVariance,
    )?;

    Ok(CorrelatedStackReceipt {
        model: model.clone(),
        terms: terms.to_vec().into_boxed_slice(),
        independent_standard_deviation,
        independent_variance,
        correlated_standard_deviation,
        correlated_variance,
        correlation_variance_delta: correlated_variance - independent_variance,
    })
}

/// The linearization's verdict against sampled tolerance-band extremes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RobustnessVerdict {
    /// The first-order predicted QoI standard deviation (`√budget`).
    pub linearized_std: f64,
    /// The largest `|QoI − nominal|` observed at a sampled band extreme.
    pub sampled_max_deviation: f64,
    /// Did the extremes stay within `k · linearized_std · (1 + margin)` (the
    /// linearization held)?
    pub confirmed: bool,
}

/// Check the linearized allocation against the QoI evaluated at sampled
/// tolerance-band EXTREMES. `extreme_qois` are QoI values at `±t` corners;
/// `nominal_qoi` is the on-design value. `margin` is non-negative slack.
///
/// # Errors
///
/// Refuses an empty sample set, non-finite inputs or allocation variance,
/// negative `margin`, non-positive `k`, and unrepresentable derived values.
pub fn robustness_check(
    allocation: &Allocation,
    extreme_qois: &[f64],
    nominal_qoi: f64,
    k: f64,
    margin: f64,
) -> Result<RobustnessVerdict, ToleranceError> {
    if extreme_qois.is_empty() {
        return Err(ToleranceError::NoExtremeSamples);
    }
    validate_allocation(allocation)?;
    validate_positive_argument("k", k)?;
    if !nominal_qoi.is_finite() {
        return Err(ToleranceError::InvalidArgument {
            argument: "nominal_qoi",
            issue: ScalarIssue::NonFinite,
        });
    }
    if !margin.is_finite() {
        return Err(ToleranceError::InvalidArgument {
            argument: "margin",
            issue: ScalarIssue::NonFinite,
        });
    }
    if margin < 0.0 {
        return Err(ToleranceError::InvalidArgument {
            argument: "margin",
            issue: ScalarIssue::Negative,
        });
    }
    let linearized_std = allocation.achieved_variance.sqrt();
    validate_nonnegative_derived(
        DerivedQuantity::LinearizedStandardDeviation,
        None,
        linearized_std,
    )?;
    let mut sampled_max_deviation = 0.0_f64;
    for (index, &qoi) in extreme_qois.iter().enumerate() {
        if !qoi.is_finite() {
            return Err(ToleranceError::InvalidExtremeQoi {
                index,
                issue: ScalarIssue::NonFinite,
            });
        }
        let deviation = (qoi - nominal_qoi).abs();
        if !deviation.is_finite() {
            return Err(ToleranceError::InvalidExtremeDerived {
                index,
                quantity: DerivedQuantity::SampledDeviation,
                issue: ScalarIssue::NonFinite,
            });
        }
        if deviation > sampled_max_deviation {
            sampled_max_deviation = deviation;
        }
    }
    // an extreme lives at ~k·σ; the linearization holds if the observed extreme
    // does not exceed that by more than the margin.
    let bound = k * linearized_std * (1.0 + margin);
    validate_nonnegative_derived(DerivedQuantity::RobustnessBound, None, bound)?;
    Ok(RobustnessVerdict {
        linearized_std,
        sampled_max_deviation,
        confirmed: sampled_max_deviation <= bound,
    })
}

/// One GD&T suggestion carrying its justification.
#[derive(Debug, Clone, PartialEq)]
pub struct Suggestion {
    /// The feature.
    pub name: String,
    /// The suggested tolerance.
    pub tolerance: f64,
    /// Tighten / loosen / unchanged.
    pub action: Action,
    /// The certified sensitivity that justifies it.
    pub certified_sensitivity: f64,
    /// The color of that sensitivity.
    pub color: ColorRank,
}

/// Build the GD&T suggestion report: every entry (and in particular every
/// LOOSENED tolerance) carries the certified sensitivity that justifies it.
///
/// # Errors
///
/// Refuses forged or deserialized allocations containing unstable/ambiguous
/// names or non-positive/non-finite tolerance and sensitivity fields.
pub fn gdt_report(allocation: &Allocation) -> Result<Vec<Suggestion>, ToleranceError> {
    validate_allocation(allocation)?;
    Ok(allocation
        .items
        .iter()
        .map(|item| Suggestion {
            name: item.name.clone(),
            tolerance: item.tolerance,
            action: item.action,
            certified_sensitivity: item.sensitivity,
            color: item.sensitivity_color,
        })
        .collect())
}

fn validate_allocation(allocation: &Allocation) -> Result<(), ToleranceError> {
    if allocation.items.is_empty() {
        return Err(ToleranceError::NoFeatures);
    }
    validate_positive_argument("allocation.total_cost", allocation.total_cost)?;
    validate_positive_argument("allocation.achieved_variance", allocation.achieved_variance)?;
    let mut canonical_names = BTreeMap::new();
    for (index, item) in allocation.items.iter().enumerate() {
        let canonical_name = canonical_feature_name(index, &item.name)?;
        if let Some(&first_index) = canonical_names.get(&canonical_name) {
            return Err(ToleranceError::AmbiguousFeatureName {
                first_index,
                duplicate_index: index,
                canonical_name,
            });
        }
        canonical_names.insert(canonical_name, index);
        validate_allocation_item(index, &item.name, "tolerance", item.tolerance)?;
        validate_allocation_item(index, &item.name, "sensitivity", item.sensitivity)?;
    }
    Ok(())
}

fn validate_allocation_item(
    index: usize,
    name: &str,
    field: &'static str,
    value: f64,
) -> Result<(), ToleranceError> {
    let issue = if !value.is_finite() {
        Some(ScalarIssue::NonFinite)
    } else if value <= 0.0 {
        Some(ScalarIssue::NonPositive)
    } else {
        None
    };
    if let Some(issue) = issue {
        return Err(ToleranceError::InvalidAllocationItem {
            index,
            name: name.to_string(),
            field,
            issue,
        });
    }
    Ok(())
}

/// The QoI variance budget for a two-sided `P(|QoI − nominal| ≤ spec_margin) ≥
/// target`: `budget = (spec_margin / z)²` with `z = Φ⁻¹((1 + target) / 2)`.
///
/// # Errors
/// [`ToleranceError`] if `target ∉ (0, 1)`, `spec_margin ≤ 0`, any argument is
/// non-finite, or the derived quantile/budget is not positive and finite.
pub fn variance_budget(spec_margin: f64, target: f64) -> Result<f64, ToleranceError> {
    validate_positive_argument("spec_margin", spec_margin)?;
    if !target.is_finite() {
        return Err(ToleranceError::InvalidArgument {
            argument: "target",
            issue: ScalarIssue::NonFinite,
        });
    }
    if !(target > 0.0 && target < 1.0) {
        return Err(ToleranceError::InvalidArgument {
            argument: "target",
            issue: ScalarIssue::OutsideOpenUnitInterval,
        });
    }
    let z = two_sided_normal_quantile(target);
    validate_positive_derived(DerivedQuantity::NormalQuantile, None, z)?;
    let sigma = spec_margin / z;
    validate_positive_derived(DerivedQuantity::VarianceBudget, None, sigma)?;
    let budget = sigma * sigma;
    validate_positive_derived(DerivedQuantity::VarianceBudget, None, budget)?;
    Ok(budget)
}

/// Positive normal quantile `Φ⁻¹((1 + target) / 2)` for a two-sided central
/// probability, using Acklam's rational approximation. It evaluates directly
/// from `target / 2` in the central region and `(1 - target) / 2` in the upper
/// tail so representable targets adjacent to zero or one never round the CDF
/// probability to exactly `0.5` or `1.0` first.
#[allow(clippy::unreadable_literal, clippy::excessive_precision)]
fn two_sided_normal_quantile(target: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const CENTRAL_TARGET_LIMIT: f64 = 0.9515;
    if target <= CENTRAL_TARGET_LIMIT {
        let q = target * 0.5;
        let r = q * q;
        let numerator = ((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5];
        let denominator = ((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0;
        // Reassociate the final multiply so the smallest positive target is
        // not lost by computing target/2 before multiplying by a ~2.5 factor.
        target * (0.5 * numerator / denominator)
    } else {
        let upper_tail = (1.0 - target) * 0.5;
        let q = (-2.0 * det::ln(upper_tail)).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}
