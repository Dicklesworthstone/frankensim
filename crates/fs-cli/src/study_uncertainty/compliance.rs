//! Predeclared Bernoulli decisions for the native numerical cooling event.
//! The statistical owner replays raw observations; this module only applies
//! the retained policy and projects its estimated sampling evidence.

use fs_blake3::ContentHash;
use fs_package::Claim;
use fs_project::uncertainty::CompliancePolicy;
use fs_uq::{AnytimeEstimate, UqExecution, UqStatus};

use super::{Model, Result, fail, input_law, optional, quoted};

const SCOPE: &str = "The predeclared Beta(1/2,1/2) Bernoulli likelihood mixture gives a time-uniform mathematical confidence sequence under a fixed conditional success probability. Alpha, threshold, input laws and stopping policy remain fixed across resume. Floating-point inversion is not an outward-rounded certificate. No fixed-count temperature distribution is claimed from a policy-stopped sample. Child engineering uncertainty budgets and verdicts remain unchanged; numerical, geometric, material and physical-model errors are not bounded by this sampling interval. Refused observations are never replaced or skipped, and a refused execution has no probability-confidence claim.";

pub(super) fn scope(model: &Model) -> String {
    format!("Estimated probability of the declared native numerical temperature event under {}. {SCOPE}", input_law(model))
}

pub(super) struct Assessment {
    policy: CompliancePolicy,
    count: usize,
    estimate: Option<AnytimeEstimate>,
    decision: &'static str,
    scope: String,
}

pub(super) fn assess(model: &Model, execution: &UqExecution) -> Result<Option<Assessment>> {
    let Some(policy) = model.bound.study().compliance().copied() else {
        return Ok(None);
    };
    let count = execution.observations().len();
    // A terminal failed execution must not present its accepted prefix as an
    // uncensored Bernoulli sample, even if that prefix would cross a boundary.
    let estimate = if execution.report().status == UqStatus::Refused {
        None
    } else {
        execution
            .assess_bernoulli_compliance(policy.alpha, 0.0)
            .map_err(|error| fail("cli-uncertainty-compliance", error.to_string()))?
    };
    let decision = match estimate {
        Some(interval)
            if count >= policy.min_samples && interval.lo >= policy.required_probability =>
        {
            "meets-probability-target"
        }
        Some(interval)
            if count >= policy.min_samples && interval.hi < policy.required_probability =>
        {
            "below-probability-target"
        }
        _ => "indeterminate",
    };
    Ok(Some(Assessment {
        policy,
        count,
        estimate,
        decision,
        scope: scope(model),
    }))
}

impl Assessment {
    pub(super) fn resolved(&self) -> bool {
        self.decision != "indeterminate"
    }

    pub(super) fn json(&self) -> String {
        let interval = self.estimate.map_or_else(
            || "null".into(),
            |estimate| format!("[{},{}]", estimate.lo, estimate.hi),
        );
        format!(
            "{{\"method\":\"bernoulli-beta-half-mixture\",\"required_probability\":{},\"alpha\":{},\"min_decision_samples\":{},\"samples_evaluated\":{},\"empirical_probability_of_compliance\":{},\"probability_confidence_sequence\":{interval},\"decision\":{}}}",
            self.policy.required_probability,
            self.policy.alpha,
            self.policy.min_samples,
            self.count,
            optional(self.estimate.map(|estimate| estimate.mean)),
            quoted(self.decision)
        )
    }

    pub(super) fn html(&self) -> String {
        let interval = self.estimate.map_or_else(
            || "unavailable".into(),
            |estimate| {
                format!(
                    "[{}, {}]; empirical frequency {}",
                    estimate.lo, estimate.hi, estimate.mean
                )
            },
        );
        format!(
            "<h2>Probability decision</h2><p>Decision: {}. Required probability: {}; alpha: {}; minimum decision samples: {}.</p><p>Bernoulli Beta(1/2,1/2) confidence sequence: {interval}.</p>",
            self.decision,
            self.policy.required_probability,
            self.policy.alpha,
            self.policy.min_samples
        )
    }

    pub(super) fn claim(&self, report: ContentHash) -> Option<Claim> {
        self.estimate.map(|estimate| {
            Claim::estimated(
                "cooling.uncertainty.compliance-probability",
                format!(
                    "Empirical frequency {} across {} completed native observations; probability confidence sequence [{}, {}] at predeclared alpha {}. Required probability {}; decision {} after minimum {}. Retained result {}. {}",
                    estimate.mean, self.count, estimate.lo, estimate.hi, self.policy.alpha,
                    self.policy.required_probability, self.decision, self.policy.min_samples, report.to_hex(), self.scope
                ),
                "predeclared-bernoulli-beta-half-mixture-confidence-sequence",
                (estimate.mean - estimate.lo).max(estimate.hi - estimate.mean),
            )
        })
    }
}
