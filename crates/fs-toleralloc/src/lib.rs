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
//! justifies it. Deterministic; depends only on `fs-evidence`.

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

/// A structured allocation failure.
#[derive(Debug, Clone, PartialEq)]
pub enum ToleranceError {
    /// No features.
    NoFeatures,
    /// A non-positive sensitivity, cost, or baseline.
    NonPositive {
        /// The offending feature.
        feature: String,
        /// What was non-positive.
        what: &'static str,
    },
    /// The variance budget or `k` is not positive.
    BadBudget,
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
/// [`ToleranceError`] on empty input, a non-positive feature field, or a
/// non-positive budget / `k`.
pub fn allocate(
    features: &[Feature],
    variance_budget: f64,
    k: f64,
) -> Result<Allocation, ToleranceError> {
    if features.is_empty() {
        return Err(ToleranceError::NoFeatures);
    }
    if !(variance_budget > 0.0 && k > 0.0) {
        return Err(ToleranceError::BadBudget);
    }
    for f in features {
        let bad = if f.sensitivity <= 0.0 {
            Some("sensitivity")
        } else if f.cost_coeff <= 0.0 {
            Some("cost_coeff")
        } else if f.baseline_tolerance <= 0.0 {
            Some("baseline_tolerance")
        } else {
            None
        };
        if let Some(what) = bad {
            return Err(ToleranceError::NonPositive {
                feature: f.name.clone(),
                what,
            });
        }
    }

    // shape: tᵢ ∝ (cᵢ / sᵢ²)^{1/3} — loose where sensitivity is small.
    let raw: Vec<f64> = features
        .iter()
        .map(|f| (f.cost_coeff / (f.sensitivity * f.sensitivity)).cbrt())
        .collect();
    // variance with the raw shape, then rescale so it exactly meets the budget.
    let var_raw: f64 = features
        .iter()
        .zip(&raw)
        .map(|(f, &t)| (f.sensitivity * f.sensitivity / (k * k)) * t * t)
        .sum();
    let scale = (variance_budget / var_raw).sqrt();

    let mut items = Vec::with_capacity(features.len());
    let mut total_cost = 0.0;
    let mut achieved_variance = 0.0;
    for (f, &t_raw) in features.iter().zip(&raw) {
        let tolerance = t_raw * scale;
        total_cost += f.cost_coeff / tolerance;
        achieved_variance += (f.sensitivity * f.sensitivity / (k * k)) * tolerance * tolerance;
        let action = action_for(tolerance, f.baseline_tolerance);
        items.push(TolItem {
            name: f.name.clone(),
            tolerance,
            sensitivity: f.sensitivity,
            sensitivity_color: f.sensitivity_color,
            action,
        });
    }
    Ok(Allocation {
        items,
        total_cost,
        achieved_variance,
    })
}

fn action_for(tolerance: f64, baseline: f64) -> Action {
    let ratio = tolerance / baseline;
    if ratio > 1.01 {
        Action::Loosen
    } else if ratio < 0.99 {
        Action::Tighten
    } else {
        Action::Unchanged
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
/// `nominal_qoi` is the on-design value. `margin` is the allowed slack.
#[must_use]
pub fn robustness_check(
    allocation: &Allocation,
    extreme_qois: &[f64],
    nominal_qoi: f64,
    k: f64,
    margin: f64,
) -> RobustnessVerdict {
    let linearized_std = allocation.achieved_variance.max(0.0).sqrt();
    let sampled_max_deviation = extreme_qois
        .iter()
        .map(|q| (q - nominal_qoi).abs())
        .fold(0.0_f64, f64::max);
    // an extreme lives at ~k·σ; the linearization holds if the observed extreme
    // does not exceed that by more than the margin.
    let bound = k * linearized_std * (1.0 + margin);
    RobustnessVerdict {
        linearized_std,
        sampled_max_deviation,
        confirmed: sampled_max_deviation <= bound,
    }
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
#[must_use]
pub fn gdt_report(allocation: &Allocation) -> Vec<Suggestion> {
    allocation
        .items
        .iter()
        .map(|i| Suggestion {
            name: i.name.clone(),
            tolerance: i.tolerance,
            action: i.action,
            certified_sensitivity: i.sensitivity,
            color: i.sensitivity_color,
        })
        .collect()
}

/// The QoI variance budget for a two-sided `P(|QoI − nominal| ≤ spec_margin) ≥
/// target`: `budget = (spec_margin / z)²` with `z = Φ⁻¹((1 + target) / 2)`.
///
/// # Errors
/// [`ToleranceError::BadBudget`] if `target ∉ (0, 1)` or `spec_margin ≤ 0`.
pub fn variance_budget(spec_margin: f64, target: f64) -> Result<f64, ToleranceError> {
    if !(target > 0.0 && target < 1.0) || spec_margin <= 0.0 {
        return Err(ToleranceError::BadBudget);
    }
    let z = inverse_normal_cdf(f64::midpoint(1.0, target));
    let sigma = spec_margin / z;
    Ok(sigma * sigma)
}

/// Standard-normal inverse CDF `Φ⁻¹(p)` (Acklam's rational approximation). The
/// coefficient tables are published numerical constants.
#[allow(clippy::unreadable_literal, clippy::excessive_precision)]
fn inverse_normal_cdf(p: f64) -> f64 {
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
    let plow = 0.02425;
    let phigh = 1.0 - plow;
    if p < plow {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= phigh {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}
