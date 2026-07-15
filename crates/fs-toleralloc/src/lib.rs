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
//! DETERMINISM DOCTRINE (bead frankensim-lyms): every transcendental in
//! this crate routes through `fs_math::det` so the "fully deterministic"
//! contract holds cross-ISA by construction — platform libm `ln`/`exp`
//! differ by ≥1 ULP across ISAs and libm versions. `sqrt` stays primitive
//! (IEEE-754 requires correct rounding for it).

use std::collections::BTreeMap;

use fs_math::det;

pub use fs_evidence::ColorRank;

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
