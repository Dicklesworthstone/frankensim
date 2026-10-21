//! Cost-aware fixed production allocation from pilot correction dispersion.
//!
//! Uses Giles's N_l proportional to sqrt(V_l/C_l), with an explicit minimum
//! and caller-supplied dispersion floor. Costs price PHYSICAL level solves;
//! a correction above level zero pays for both fine and coarse solves.
//! The predicted sampling error is a pilot-based estimate, not a bound on
//! future observations, discretization bias, or model error.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`.

use crate::product_sensitivity::Sum;

/// Dimensional allocation policy. There is deliberately no unit-free default.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MlmcAllocationPolicy {
    /// Target descriptive standard error in the physical QoI's units.
    pub target_standard_error: f64,
    /// At least two NEW production corrections on every declared level.
    pub minimum_samples_per_level: usize,
    /// Positive QoI-unit safeguard: a constant pilot does not prove zero variance.
    /// This floor affects allocation only, never reported observed statistics.
    pub correction_std_dev_floor: f64,
    /// Remaining physical-call cap, including both halves of higher-level pairs.
    pub maximum_model_calls: usize,
}

/// A complete fixed allocation. Increasing minimum counts can over-allocate
/// relative to the continuous optimum; exact integer cost optimality is not claimed.
#[derive(Debug, Clone, PartialEq)]
pub struct MlmcAllocation {
    pub samples_per_level: Vec<usize>,
    pub model_calls: usize,
    /// Pilot standard deviations after applying the declared allocation floor.
    pub allocation_standard_deviations: Vec<f64>,
    pub predicted_standard_error: f64,
    /// Sum of allocated corrections times their declared fine/coarse cost.
    /// None means the dimensional total overflows, not that work is free.
    pub declared_cost: Option<f64>,
}

/// Planning refuses before any production solve or unbounded allocation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MlmcAllocationError {
    Invalid(&'static str),
    BudgetExceeded { maximum_model_calls: usize },
    Numerical(&'static str),
}
impl core::fmt::Display for MlmcAllocationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) | Self::Numerical(reason) => f.write_str(reason),
            Self::BudgetExceeded { maximum_model_calls } => write!(
                f, "MLMC allocation cannot fit the {maximum_model_calls} remaining physical calls"
            ),
        }
    }
}
impl std::error::Error for MlmcAllocationError {}

/// Plan a fixed finest-level estimator from per-correction pilot standard
/// deviations and caller-declared costs of ONE physical solve at each level.
/// The same declared ladder is retained; no level is dropped or extrapolated.
///
/// Counts are rounded up from the continuous Giles allocation, then bounded
/// below by the explicit minimum. A bounded rounding repair checks the actual
/// represented predicted standard error. Costs are normalized before summing
/// fine/coarse prices, and error norms are scaled before squaring. The target
/// itself is never squared, so tiny or huge dimensional tolerances are usable.
/// Invalid data, unsupported numerical range, and an insufficient call cap
/// return an error rather than a smaller allocation masquerading as success.
pub fn allocate_mlmc(
    correction_std_devs: &[f64],
    physical_model_costs: &[f64],
    policy: MlmcAllocationPolicy,
) -> Result<MlmcAllocation, MlmcAllocationError> {
    validate_policy(physical_model_costs, policy)?;
    if correction_std_devs.len() != physical_model_costs.len()
        || correction_std_devs.iter().any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(MlmcAllocationError::Invalid(
            "MLMC needs one finite nonnegative correction deviation per physical level",
        ));
    }
    let budget_error = || MlmcAllocationError::BudgetExceeded {
        maximum_model_calls: policy.maximum_model_calls,
    };
    let deviations: Vec<_> = correction_std_devs.iter()
        .map(|s| s.max(policy.correction_std_dev_floor)).collect();
    let cost_scale = physical_model_costs.iter().copied().fold(0.0_f64, f64::max);
    let costs: Vec<_> = physical_model_costs.iter().enumerate().map(|(level, cost)| {
        cost / cost_scale + if level == 0 { 0.0 } else {
            physical_model_costs[level - 1] / cost_scale
        }
    }).collect();
    if costs.iter().any(|c| *c == 0.0) {
        return Err(MlmcAllocationError::Numerical(
            "MLMC relative correction costs underflow; respecify the cost model",
        ));
    }
    let roots: Vec<_> = costs.iter().map(|&c| fs_math::det::sqrt(c)).collect();
    let ratios: Vec<_> = deviations.iter().map(|s| s / policy.target_standard_error).collect();
    // Even spending the entire physical budget on one level could not meet
    // this necessary condition. Check before multiplying allocation weights.
    if ratios.iter().any(|r| !r.is_finite() || *r > fs_math::det::sqrt(policy.maximum_model_calls as f64)) {
        return Err(budget_error());
    }
    let mut sum = Sum::default();
    for (&ratio, &root) in ratios.iter().zip(&roots) {
        sum.add(ratio * root);
    }
    let mut counts = Vec::with_capacity(deviations.len());
    for (&ratio, &root) in ratios.iter().zip(&roots) {
        let raw = (ratio / root) * sum.value;
        if !raw.is_finite() || raw > policy.maximum_model_calls as f64 {
            return Err(budget_error());
        }
        counts.push((raw.ceil() as usize).max(policy.minimum_samples_per_level));
    }
    let mut calls = physical_calls(&counts).ok_or_else(budget_error)?;
    if calls > policy.maximum_model_calls { return Err(budget_error()); }
    let scale = deviations.iter().copied().fold(0.0_f64, f64::max);
    // Only repair f64 allocation/ceiling roundoff, not an adaptive data-dependent
    // production stop. At most one extra correction per level plus one is allowed.
    let mut predicted = predicted_error(&deviations, &counts, scale)?;
    for _ in 0..=counts.len() {
        if predicted <= policy.target_standard_error { break; }
        let best = (0..counts.len()).max_by(|&a, &b| {
            let gain = |i: usize| {
                let n = counts[i] as f64;
                (deviations[i] / scale).powi(2) / n / (n + 1.0)
            };
            gain(a).total_cmp(&gain(b)).then_with(|| b.cmp(&a))
        }).ok_or(MlmcAllocationError::Numerical("MLMC allocation has no levels"))?;
        calls += if best == 0 { 1 } else { 2 };
        if calls > policy.maximum_model_calls { return Err(budget_error()); }
        counts[best] += 1;
        predicted = predicted_error(&deviations, &counts, scale)?;
    }
    if predicted > policy.target_standard_error {
        return Err(MlmcAllocationError::Numerical("MLMC allocation rounding did not meet its predicted target"));
    }
    let mut cost = Sum::default();
    for (&n, &c) in counts.iter().zip(&costs) { cost.add(n as f64 * c); }
    let cost = cost.value * cost_scale;
    Ok(MlmcAllocation {
        samples_per_level: counts, model_calls: calls,
        allocation_standard_deviations: deviations,
        predicted_standard_error: predicted,
        declared_cost: cost.is_finite().then_some(cost),
    })
}

pub(crate) fn physical_calls(counts: &[usize]) -> Option<usize> {
    counts.iter().enumerate().try_fold(0usize, |sum, (level, &n)| {
        n.checked_mul(if level == 0 { 1 } else { 2 }).and_then(|v| sum.checked_add(v))
    })
}

pub(crate) fn validate_policy(
    costs: &[f64], policy: MlmcAllocationPolicy,
) -> Result<(), MlmcAllocationError> {
    if costs.is_empty() || costs.len() > 64 || costs.iter().any(|c| !c.is_finite() || *c <= 0.0)
        || !policy.target_standard_error.is_finite() || policy.target_standard_error <= 0.0
        || !policy.correction_std_dev_floor.is_finite() || policy.correction_std_dev_floor <= 0.0
        || policy.minimum_samples_per_level < 2 || policy.maximum_model_calls > 1_000_000
    {
        return Err(MlmcAllocationError::Invalid(
            "MLMC needs 1..=64 positive physical costs, positive finite error/floor, at least two corrections, and at most one million calls",
        ));
    }
    let minimum = policy.minimum_samples_per_level.checked_mul(2 * costs.len() - 1);
    if minimum.is_none_or(|n| n > policy.maximum_model_calls) {
        return Err(MlmcAllocationError::BudgetExceeded { maximum_model_calls: policy.maximum_model_calls });
    }
    Ok(())
}

fn predicted_error(deviations: &[f64], counts: &[usize], scale: f64) -> Result<f64, MlmcAllocationError> {
    let mut variance = Sum::default();
    for (&deviation, &count) in deviations.iter().zip(counts) {
        variance.add((deviation / scale).powi(2) / count as f64);
    }
    let error = fs_math::det::sqrt(variance.value) * scale;
    if !error.is_finite() || error == 0.0 {
        Err(MlmcAllocationError::Numerical("MLMC predicted standard error is not representable"))
    } else { Ok(error) }
}
