//! Decision-plan dimensional adoption (bead sj31i.7.4): typed decision
//! plans over OED campaigns with the shared inference dimensional core.
//!
//! A [`DecisionPlan`] binds the campaign's objective schema, the canonical
//! alternative set, dimensionless allocation coefficients, typed cost and
//! utility measures, and the dimensionless information gain into one
//! checked record with a domain-separated receipt identity. The already
//! landed scalar campaign slice stays untouched as the implementation
//! contribution; [`verify_decision_plan`] is the independent proof lane
//! that recomputes the plan's quantitative claims from the report.
//!
//! Distinctions enforced mechanically:
//!
//! - information gain is dimensionless (prior/posterior variance share
//!   dimensions, so their reduction ratio carries none);
//! - utility terms carry the objective schema (physical or explicitly
//!   dimensionless), never silently interchangeable with cost;
//! - costs are [`DecisionMeasure::Cost`] and utilities
//!   [`DecisionMeasure::Utility`]; cross-kind algebra refuses by type;
//! - allocation coefficients are finite, non-negative, dimensionless, and
//!   sum to one in deterministic canonical order.

use fs_qty::inference::{DecisionMeasure, InferenceError};

use crate::{OedError, OedReport};

/// Wire prefix of the decision-plan receipt identity.
pub const DECISION_PLAN_RECEIPT_PREFIX: &str = "oed-decision-plan:v1:";

const PLAN_RECEIPT_DOMAIN: &str = "org.frankensim.fs-oed-e2e.decision-plan.v1";

/// Structured decision-plan refusals; every mismatch names both sides.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanError {
    /// A plan needs at least one alternative.
    EmptyAlternatives,
    /// Too many alternatives for the checked envelope.
    AlternativeLimit {
        /// Requested count.
        count: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// Alternatives must be unique after canonical ordering.
    DuplicateAlternative {
        /// The duplicated name.
        name: String,
    },
    /// An allocation coefficient must be finite, non-negative, and
    /// dimensionless.
    InvalidCoefficient {
        /// Offending alternative.
        alternative: String,
        /// Why the coefficient was refused.
        reason: &'static str,
    },
    /// Allocation coefficients must sum to one within exact tolerance.
    AllocationNotNormalized {
        /// Computed sum in deterministic canonical order.
        sum: f64,
    },
    /// Information gain must be finite and dimensionless: the prior and
    /// posterior variance dimensions disagreed.
    InformationGainNotDimensionless {
        /// Prior variance dimensions rendered.
        prior: String,
        /// Posterior variance dimensions rendered.
        posterior: String,
    },
    /// A negative or non-finite information gain is not a decision input.
    InvalidInformationGain,
    /// The report and the plan disagree on the objective schema.
    ObjectiveSchemaMismatch,
    /// The campaign report could not be read for plan extraction.
    Report(crate::OedError),
}

impl core::fmt::Display for PlanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyAlternatives => write!(f, "a decision plan needs at least one alternative"),
            Self::AlternativeLimit { count, max } => {
                write!(
                    f,
                    "alternative count {count} exceeds the admitted envelope {max}"
                )
            }
            Self::DuplicateAlternative { name } => {
                write!(
                    f,
                    "alternative {name:?} is duplicated after canonical ordering"
                )
            }
            Self::InvalidCoefficient {
                alternative,
                reason,
            } => write!(
                f,
                "allocation coefficient for {alternative:?} is invalid: {reason}"
            ),
            Self::AllocationNotNormalized { sum } => {
                write!(f, "allocation coefficients sum to {sum}, not one")
            }
            Self::InformationGainNotDimensionless { prior, posterior } => write!(
                f,
                "information gain is not dimensionless: prior variance {prior} vs posterior variance {posterior}"
            ),
            Self::InvalidInformationGain => {
                write!(f, "information gain must be finite and non-negative")
            }
            Self::ObjectiveSchemaMismatch => {
                write!(f, "report and plan disagree on the objective schema")
            }
            Self::Report(error) => write!(f, "report extraction failed: {error}"),
        }
    }
}

impl std::error::Error for PlanError {}

impl From<OedError> for PlanError {
    fn from(error: OedError) -> Self {
        Self::Report(error)
    }
}

impl From<InferenceError> for PlanError {
    fn from(error: InferenceError) -> Self {
        match error {
            InferenceError::NonFiniteDecisionValue => Self::InvalidInformationGain,
            other => Self::InvalidCoefficient {
                alternative: "decision-measure".to_string(),
                reason: match other {
                    InferenceError::DecisionMeasureMismatch { .. } => "cost/utility measures mixed",
                    _ => "decision measure refused",
                },
            },
        }
    }
}

/// One typed allocation coefficient: a finite, non-negative, dimensionless
/// weight bound to one named alternative.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AllocationCoefficient {
    alternative_index: usize,
    weight: f64,
}

impl AllocationCoefficient {
    /// The alternative's canonical index.
    #[must_use]
    pub const fn alternative_index(self) -> usize {
        self.alternative_index
    }

    /// The dimensionless weight.
    #[must_use]
    pub const fn weight(self) -> f64 {
        self.weight
    }
}

/// A checked decision plan: canonical alternatives, dimensionless
/// normalized allocation coefficients, typed cost and utility measures,
/// the dimensionless information gain, and a domain-separated receipt
/// identity binding the objective schema and numeric policy.
#[derive(Debug, Clone, PartialEq)]
pub struct DecisionPlan {
    alternatives: Vec<String>,
    coefficients: Vec<AllocationCoefficient>,
    cost: DecisionMeasure,
    utility: DecisionMeasure,
    information_gain: f64,
    receipt_identity: String,
}

impl DecisionPlan {
    /// Canonical alternatives (sorted, deduplicated).
    #[must_use]
    pub fn alternatives(&self) -> &[String] {
        &self.alternatives
    }

    /// Allocation coefficients in canonical alternative order.
    #[must_use]
    pub fn coefficients(&self) -> &[AllocationCoefficient] {
        &self.coefficients
    }

    /// The typed cost measure.
    #[must_use]
    pub const fn cost(&self) -> DecisionMeasure {
        self.cost
    }

    /// The typed utility measure.
    #[must_use]
    pub const fn utility(&self) -> DecisionMeasure {
        self.utility
    }

    /// The dimensionless information gain (variance-reduction ratio).
    #[must_use]
    pub const fn information_gain(&self) -> f64 {
        self.information_gain
    }

    /// The plan receipt identity `oed-decision-plan:v1:<64 lowercase hex>`.
    #[must_use]
    pub fn receipt_identity(&self) -> &str {
        &self.receipt_identity
    }

    /// Admit a plan from typed parts. Alternatives are canonically sorted;
    /// coefficients are checked dimensionless and normalized in that
    /// canonical order; cost and utility are typed measures that cannot
    /// mix by construction.
    ///
    /// # Errors
    /// Returns [`PlanError`] for empty/oversized/duplicated alternatives,
    /// non-dimensionless or non-normalized coefficients, non-finite
    /// measures, or a negative information gain.
    pub fn try_new(
        alternatives: Vec<String>,
        weights: Vec<f64>,
        cost: DecisionMeasure,
        utility: DecisionMeasure,
        information_gain: f64,
        objective_spec: crate::ObjectiveSpec,
    ) -> Result<Self, PlanError> {
        if alternatives.is_empty() {
            return Err(PlanError::EmptyAlternatives);
        }
        if alternatives.len() > crate::MAX_CAMPAIGN_CANDIDATES {
            return Err(PlanError::AlternativeLimit {
                count: alternatives.len(),
                max: crate::MAX_CAMPAIGN_CANDIDATES,
            });
        }
        if alternatives.len() != weights.len() {
            return Err(PlanError::InvalidCoefficient {
                alternative: "plan".to_string(),
                reason: "weights must cover every alternative",
            });
        }
        if !matches!(cost, DecisionMeasure::Cost(_)) {
            return Err(PlanError::InvalidCoefficient {
                alternative: "cost".to_string(),
                reason: "cost slot requires a DecisionMeasure::Cost",
            });
        }
        if !matches!(utility, DecisionMeasure::Utility(_)) {
            return Err(PlanError::InvalidCoefficient {
                alternative: "utility".to_string(),
                reason: "utility slot requires a DecisionMeasure::Utility",
            });
        }
        if !information_gain.is_finite() || information_gain < 0.0 {
            return Err(PlanError::InvalidInformationGain);
        }
        let mut ordered: Vec<(String, f64)> = alternatives.into_iter().zip(weights).collect();
        ordered.sort_by(|left, right| left.0.cmp(&right.0));
        for window in ordered.windows(2) {
            if window[0].0 == window[1].0 {
                return Err(PlanError::DuplicateAlternative {
                    name: window[0].0.clone(),
                });
            }
        }
        let mut sum = 0.0_f64;
        let mut coefficients = Vec::with_capacity(ordered.len());
        let mut alternatives = Vec::with_capacity(ordered.len());
        for (index, (alternative, weight)) in ordered.into_iter().enumerate() {
            if !weight.is_finite() {
                return Err(PlanError::InvalidCoefficient {
                    alternative,
                    reason: "not finite",
                });
            }
            if weight < 0.0 {
                return Err(PlanError::InvalidCoefficient {
                    alternative,
                    reason: "negative",
                });
            }
            sum += weight;
            coefficients.push(AllocationCoefficient {
                alternative_index: index,
                weight,
            });
            alternatives.push(alternative);
        }
        if (sum - 1.0).abs() > 1e-9 {
            return Err(PlanError::AllocationNotNormalized { sum });
        }
        let receipt_identity = plan_identity(
            &alternatives,
            &coefficients,
            cost,
            utility,
            information_gain,
            objective_spec,
        );
        Ok(Self {
            alternatives,
            coefficients,
            cost,
            utility,
            information_gain,
            receipt_identity,
        })
    }

    /// Extract a checked decision plan from a finished campaign report:
    /// allocation from the report's recorded allocation, cost from the
    /// acquisition-cost payload, utility from the final EVPI, and
    /// information gain from the verified dimensionless variance-reduction
    /// ratio.
    ///
    /// # Errors
    /// Returns [`PlanError`] when the report's variance dimensions disagree
    /// (the gain would not be dimensionless), the gain is invalid, or the
    /// allocation is not a normalized dimensionless weighting.
    pub fn from_report(report: &OedReport) -> Result<Self, PlanError> {
        let prior = report.prior_total_variance();
        let posterior = report.posterior_total_variance();
        if prior.dims != posterior.dims {
            return Err(PlanError::InformationGainNotDimensionless {
                prior: prior.dims.unit_string(),
                posterior: posterior.dims.unit_string(),
            });
        }
        let information_gain = report.variance_reduction();
        let cost_value = report
            .allocation()
            .iter()
            .map(|(_, weight)| weight)
            .sum::<f64>();
        let cost = DecisionMeasure::cost(cost_value)?;
        let utility = DecisionMeasure::utility(report.final_evpi().value())?;
        // The report's allocation is budget amounts per alternative; the
        // plan's coefficients are the normalized dimensionless weights.
        let mut alternatives = Vec::new();
        let mut weights = Vec::new();
        if cost_value > 0.0 {
            for (alternative, amount) in report.allocation() {
                alternatives.push(alternative.clone());
                weights.push(*amount / cost_value);
            }
        } else {
            for (alternative, _) in report.allocation() {
                alternatives.push(alternative.clone());
                weights.push(0.0);
            }
        }
        // An all-zero allocation (a no-sensor campaign) is uniform over the
        // admitted alternatives, never a divide-by-zero.
        if cost_value == 0.0 && !alternatives.is_empty() {
            let uniform = 1.0 / alternatives.len() as f64;
            weights.fill(uniform);
        }
        Self::try_new(
            alternatives,
            weights,
            cost,
            utility,
            information_gain,
            report.objective_spec(),
        )
    }
}

fn plan_identity(
    alternatives: &[String],
    coefficients: &[AllocationCoefficient],
    cost: DecisionMeasure,
    utility: DecisionMeasure,
    information_gain: f64,
    objective_spec: crate::ObjectiveSpec,
) -> String {
    let mut hasher = fs_blake3::DomainHasher::new(PLAN_RECEIPT_DOMAIN);
    hasher.update(&crate::OED_REPORT_IDENTITY_VERSION.to_le_bytes());
    hasher.update(&objective_spec.quantity_spec().canonical_bytes());
    hasher.update(&(alternatives.len() as u64).to_le_bytes());
    for (alternative, coefficient) in alternatives.iter().zip(coefficients) {
        hasher.update(&(alternative.len() as u64).to_le_bytes());
        hasher.update(alternative.as_bytes());
        hasher.update(&coefficient.weight.to_le_bytes());
    }
    hasher.update(cost.kind_name().as_bytes());
    hasher.update(&cost.value().to_le_bytes());
    hasher.update(utility.kind_name().as_bytes());
    hasher.update(&utility.value().to_le_bytes());
    hasher.update(&information_gain.to_le_bytes());
    format!("{DECISION_PLAN_RECEIPT_PREFIX}{}", hasher.finalize())
}

/// The independent plan-verification verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanVerdict {
    /// The independent recomputation agrees with the report's recorded
    /// plan quantities.
    Verified,
    /// The recomputation decisively disagrees with the report.
    Discrepancy,
}

impl PlanVerdict {
    /// Stable wire name for receipts and structured logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Discrepancy => "discrepancy",
        }
    }
}

/// Independent check of a decision plan against its campaign report.
#[derive(Debug, Clone, PartialEq)]
pub struct PlanCheck {
    verdict: PlanVerdict,
    recomputed_information_gain: f64,
    reported_information_gain: f64,
    allocation_matches: bool,
}

impl PlanCheck {
    /// The verdict.
    #[must_use]
    pub const fn verdict(&self) -> PlanVerdict {
        self.verdict
    }

    /// Independently recomputed dimensionless information gain.
    #[must_use]
    pub const fn recomputed_information_gain(&self) -> f64 {
        self.recomputed_information_gain
    }

    /// The report's recorded information gain.
    #[must_use]
    pub const fn reported_information_gain(&self) -> f64 {
        self.reported_information_gain
    }

    /// True when the plan's allocation matches the report's record in
    /// canonical order.
    #[must_use]
    pub const fn allocation_matches(&self) -> bool {
        self.allocation_matches
    }
}

/// Independently verify a decision plan against its campaign report:
/// recompute the information gain from the recorded prior/posterior
/// variance totals, and compare the allocation in canonical order. This is
/// the independent-proof lane for the landed scalar campaign slice; a
/// `Discrepancy` is data about the plan, not an error.
///
/// # Errors
/// Returns [`PlanError`] when the report cannot be read.
pub fn verify_decision_plan(
    report: &OedReport,
    plan: &DecisionPlan,
) -> Result<PlanCheck, PlanError> {
    let prior = report.prior_total_variance();
    let posterior = report.posterior_total_variance();
    let recomputed = if prior.value > 0.0 {
        ((prior.value - posterior.value) / prior.value).max(0.0)
    } else {
        0.0
    };
    let reported = plan.information_gain();
    let gain_agrees =
        (recomputed - reported).abs() <= 64.0_f64 * f64::EPSILON * recomputed.abs().max(1.0);
    let mut canonical: Vec<(String, f64)> = report
        .allocation()
        .iter()
        .map(|(name, weight)| (name.clone(), *weight))
        .collect();
    canonical.sort_by(|left, right| left.0.cmp(&right.0));
    let total: f64 = canonical.iter().map(|(_, weight)| *weight).sum();
    let normalized: Vec<f64> = if total > 0.0 {
        canonical
            .iter()
            .map(|(_, weight)| *weight / total)
            .collect()
    } else if canonical.is_empty() {
        Vec::new()
    } else {
        vec![1.0 / canonical.len() as f64; canonical.len()]
    };
    let allocation_matches = canonical.len() == plan.alternatives().len()
        && canonical
            .iter()
            .zip(plan.alternatives())
            .all(|((name, _), alt)| name == alt)
        && normalized
            .iter()
            .zip(plan.coefficients())
            .all(|(weight, coefficient)| (weight - coefficient.weight()).abs() <= 1e-12);
    let verdict = if gain_agrees && allocation_matches {
        PlanVerdict::Verified
    } else {
        PlanVerdict::Discrepancy
    };
    Ok(PlanCheck {
        verdict,
        recomputed_information_gain: recomputed,
        reported_information_gain: reported,
        allocation_matches,
    })
}
