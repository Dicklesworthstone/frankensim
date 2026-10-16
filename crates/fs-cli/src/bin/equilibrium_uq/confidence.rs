//! A stopping policy over the existing MC execution and compliance CS owners.
//! No new sampler or confidence-sequence formula lives here. QMC is refused.
use fs_uq::{AnytimeEstimate, PropagationMethod, UqExecution, UqPlan, UqResult, UqStatus};

/// Fixed before sampling; never take the tighter bound after looking at data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Method { GaussianMixture, BernoulliMixture }
impl Method {
    pub fn parse(value: &str) -> Result<Self, &'static str> {
        match value {
            "gaussian-mixture" => Ok(Self::GaussianMixture),
            "bernoulli-mixture" => Ok(Self::BernoulliMixture),
            _ => Err("--confidence-method must be gaussian-mixture or bernoulli-mixture"),
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::GaussianMixture => "gaussian-mixture-indicator-confidence-sequence",
            Self::BernoulliMixture => "beta-half-bernoulli-mixture-confidence-sequence",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Policy {
    pub probability: f64,
    pub alpha: f64,
    pub method: Method,
}
impl Policy {
    pub fn new(probability: f64, alpha: f64) -> Result<Self, &'static str> {
        if !probability.is_finite() || probability <= 0.0 || probability >= 1.0 {
            return Err("--require-probability must be strictly inside (0,1)");
        }
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 || !alpha.recip().is_finite() {
            return Err("--confidence-alpha must be inside (0,1) with a finite reciprocal");
        }
        Ok(Self { probability, alpha, method: Method::GaussianMixture })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Decision { Satisfied, Violated, Inconclusive }
impl Decision {
    pub fn label(self) -> &'static str {
        match self {
            Self::Satisfied => "satisfied", Self::Violated => "violated", Self::Inconclusive => "inconclusive",
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct Assessment {
    pub policy: Policy,
    pub lower: f64,
    pub upper: f64,
    pub samples: u64,
    pub decision: Decision,
}
impl Assessment {
    fn from_interval(policy: Policy, interval: AnytimeEstimate) -> Self {
        // Equality is included in the requirement p >= required_probability.
        // Failure therefore requires a STRICT upper-bound separation.
        let decision = if interval.lo >= policy.probability { Decision::Satisfied }
            else if interval.hi < policy.probability { Decision::Violated }
            else { Decision::Inconclusive };
        Self { policy, lower: interval.lo, upper: interval.hi, samples: interval.n, decision }
    }
}

pub(super) struct McRun {
    pub result: UqResult,
    pub assessment: Option<Assessment>,
}

/// Assess only at completed prefixes 2,4,8,... and the original sample cap.
/// Replaying retained indicators then costs O(n) in total rather than O(n^2).
/// Original errors/cancellation return with NO assessment of a censored prefix.
/// Evaluators, laws, threshold and policy must retain the same meaning throughout.
pub(super) fn run<F, E, C>(plan: &UqPlan, policy: Option<Policy>, mut cancelled: C, mut evaluate: F)
    -> Result<McRun, Box<dyn std::error::Error>>
where F: FnMut(&[f64]) -> Result<f64, E>, E: core::fmt::Display, C: FnMut() -> bool,
{
    if plan.method != PropagationMethod::MonteCarlo {
        return Err("MC confidence decisions cannot be applied to dependent QMC points".into());
    }
    if let Some(policy) = policy {
        Policy::new(policy.probability, policy.alpha)?;
        if plan.compliance_threshold.is_none() { return Err("compliance decision needs a fixed threshold".into()); }
    }
    let mut execution = UqExecution::new(plan)?;
    let mut checkpoint = if policy.is_some() { 2 } else { plan.budget_max_samples };
    loop {
        let result = execution.advance(checkpoint - execution.observations().len(), &mut cancelled, &mut evaluate);
        if matches!(result.status, UqStatus::Cancelled | UqStatus::Refused) {
            return Ok(McRun { result, assessment: None });
        }
        let assessment = match policy {
            Some(policy) => {
                let interval = match policy.method {
                    Method::GaussianMixture => execution.assess_compliance(policy.alpha, 0.0)?,
                    Method::BernoulliMixture => execution.assess_bernoulli_compliance(policy.alpha, 0.0)?,
                }.ok_or("missing completed MC confidence prefix")?;
                Some(Assessment::from_interval(policy, interval))
            }
            None => None,
        };
        if result.status == UqStatus::Complete
            || assessment.as_ref().is_some_and(|a| a.decision != Decision::Inconclusive) {
            return Ok(McRun { result, assessment });
        }
        checkpoint = checkpoint.saturating_mul(2).min(plan.budget_max_samples);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_uq::{CorrelationModel, ParameterUncertainty};
    fn plan(n: usize) -> UqPlan {
        UqPlan::new("test-displacement", PropagationMethod::MonteCarlo, n)
            .with_correlation(CorrelationModel::Independent)
            .with_parameter(ParameterUncertainty::uniform("force", 0.0, 1.0, "N"))
            .with_compliance_threshold(0.5)
    }
    #[test]
    fn constant_events_stop_on_the_existing_cs_not_on_zero_empirical_error() {
        for (value, expected) in [(0.0, Decision::Satisfied), (1.0, Decision::Violated)] {
            let plan = plan(4096);
            let mut calls = 0;
            let report = run(&plan, Some(Policy::new(0.5, 0.05).unwrap()), || false, |_| {
                calls += 1; Ok::<_, &str>(value)
            }).unwrap();
            let assessment = report.assessment.unwrap();
            assert_eq!(assessment.decision, expected);
            assert!(calls > 2 && calls < plan.budget_max_samples);
            assert_eq!(report.result.status, UqStatus::BudgetTruncated);
            assert_eq!(report.result.samples_evaluated, calls);
            assert_eq!(assessment.samples as usize, calls);
            // Independent call into the SAME retained-observation owner.
            let mut direct = UqExecution::new(&plan).unwrap();
            direct.advance(calls, || false, |_| Ok::<_, &str>(value));
            let actual = direct.assess_compliance(0.05, 0.0).unwrap().unwrap();
            assert_eq!(assessment.lower, actual.lo); assert_eq!(assessment.upper, actual.hi);
            if expected == Decision::Satisfied { assert!(assessment.lower < 1.0); }
            else { assert!(assessment.upper > 0.0); }
        }
    }
    #[test]
    fn an_undecided_non_power_of_two_budget_retains_every_sample_and_original_statistics() {
        let plan = plan(19);
        // At this alpha and small cap the CS cannot separate p>=1/2 for ANY stream.
        let report = run(&plan, Some(Policy::new(0.5, 1e-100).unwrap()), || false,
            |x| Ok::<_, &str>(x[0])).unwrap();
        assert_eq!(report.assessment.unwrap().decision, Decision::Inconclusive);
        let mut direct = UqExecution::new(&plan).unwrap();
        let expected = direct.advance(19, || false, |x| Ok::<_, &str>(x[0]));
        assert_eq!(report.result, expected);
        let no_policy = run(&plan, None, || false, |x| Ok::<_, &str>(x[0])).unwrap();
        assert_eq!(no_policy.result, expected); assert!(no_policy.assessment.is_none());
    }
    #[test]
    fn failures_and_cancellation_do_not_publish_a_previous_prefix_decision() {
        for nonfinite in [false, true] {
            let mut calls = 0;
            let result = run(&plan(64), Some(Policy::new(0.5, 0.05).unwrap()), || false, |_| {
                calls += 1;
                if calls == 3 { if nonfinite { Ok(f64::NAN) } else { Err("physical refusal") } }
                else { Ok(0.0) }
            }).unwrap();
            assert_eq!(result.result.status, UqStatus::Refused);
            assert_eq!(calls, 3); assert!(result.assessment.is_none());
        }
        let report = run(&plan(64), Some(Policy::new(0.5, 0.05).unwrap()), || true,
            |_| -> Result<f64, &str> { panic!("pre-cancelled model must not execute") }).unwrap();
        assert_eq!(report.result.status, UqStatus::Cancelled); assert!(report.assessment.is_none());
    }
    #[test]
    fn boundary_equality_cannot_be_mislabeled_as_violated() {
        let interval = |lo, hi| AnytimeEstimate { mean: 0.5, lo, hi, n: 4, converged: false };
        let policy = Policy::new(0.5, 0.05).unwrap();
        assert_eq!(Assessment::from_interval(policy, interval(0.0, 0.5)).decision, Decision::Inconclusive);
        assert_eq!(Assessment::from_interval(policy, interval(0.5, 1.0)).decision, Decision::Satisfied);
    }
    #[test]
    fn malformed_policy_or_non_mc_sampling_refuses_before_any_model_call() {
        for (p, alpha) in [(0.0, 0.05), (1.0, 0.05), (f64::NAN, 0.05),
            (0.5, 0.0), (0.5, 1.0), (0.5, f64::INFINITY), (0.5, f64::from_bits(1))] {
            assert!(Policy::new(p, alpha).is_err());
        }
        let mut qmc = plan(64); qmc.method = PropagationMethod::QuasiMonteCarlo;
        assert!(run(&qmc, Some(Policy::new(0.5, 0.05).unwrap()), || false,
            |_| -> Result<f64, &str> { panic!("QMC must not enter the MC confidence owner") }).is_err());
    }
}


#[cfg(test)]
#[path = "confidence/bernoulli_tests.rs"]
mod bernoulli_tests;
