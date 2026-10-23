//! Sequential target decisions using the EXISTING weighted-failure CS owner.
//! A proposal frequency or descriptive standard error never licenses stopping.
use super::*;
use super::super::compliance::Decision;
use fs_uq::AnytimeEstimate;
use std::path::Path;

pub(super) struct Assessment {
    pub(super) failure: Option<AnytimeEstimate>,
    pub(super) decision: Decision,
}
impl Assessment {
    pub(super) fn reached(&self) -> bool { self.decision != Decision::Indeterminate }
}

/// None denotes a study that has not entered independently sampled production.
/// The minimum and confidence sequence count PRODUCTION calls, never pilots.
pub(super) fn assess(policy: Policy, production: Option<&ImportanceExecution>) -> Result<Assessment> {
    let failure = if let Some(execution) = production {
        policy.validate_plan(execution.plan())?;
        execution.assess_failure_probability(policy.alpha, 0.0).map_err(|error| Failure {
            code: "cooling-network-uq-confidence", message: error.to_string(),
        })?
    } else { None };
    let maximum_failure = 1.0 - policy.required_probability;
    let decision = match failure.as_ref() {
        Some(interval) if interval.n >= policy.min_samples as u64 => {
            if interval.hi <= maximum_failure { Decision::MeetsTarget }
            else if interval.lo > maximum_failure { Decision::BelowTarget }
            else { Decision::Indeterminate }
        }
        _ => Decision::Indeterminate,
    };
    Ok(Assessment { failure, decision })
}

/// Separate from fixed-count bindings AND ordinary-MC compliance. Changing any
/// policy choice on resume must refuse before a physical call or output write.
pub(super) fn binding(policy: Policy, parameters: &str) -> String {
    format!("importance-weighted-failure-decision-v1\n{}", policy.checkpoint_binding(parameters))
}

pub(super) fn termination(assessment: &Assessment, status: UqStatus) -> &'static str {
    if assessment.reached() { "probability-target" }
    else if status == UqStatus::Complete { "sample-budget" }
    else if status == UqStatus::Cancelled { "wall-time-budget" }
    else { "sample-chunk" }
}

pub(super) fn exit_code(assessment: &Assessment) -> u8 {
    if assessment.reached() { exit::SUCCESS } else { exit::BUDGET }
}

/// Only sampling-confidence fields. The caller adds physical work and proposal
/// identity, rather than reporting a biased proposal temperature distribution.
pub(super) fn fields(
    config: &Config, policy: Policy, assessment: &Assessment,
    termination: &str, checkpoint: Option<&Path>,
) -> Result<String> {
    let (failure_interval, compliance_interval, samples) = match &assessment.failure {
        Some(interval) => (
            format!("[{},{}]", number_json(interval.lo)?, number_json(interval.hi)?),
            format!("[{},{}]", number_json(1.0 - interval.hi)?, number_json(1.0 - interval.lo)?),
            interval.n,
        ),
        None => ("null".into(), "null".into(), 0),
    };
    // A resolved result has no invocation-specific path, enabling exact terminal
    // replay. Unresolved exits expose where the paid state was retained, if any.
    let recovery = if assessment.reached() { String::new() } else {
        format!(",\"checkpoint\":{}", checkpoint.map_or_else(|| "null".into(),
            |path| quote(&path.to_string_lossy())))
    };
    let no_claim = format!(
        "{} Sampling confidence for the declared numerical model under the weighted bounded-mean confidence sequence's assumptions. The proposal, failure threshold, target, alpha and minimum must be predeclared independently of production outcomes. No multiplicity control across studies, models, seeds or policies. Bounds are not outward-rounded and exclude generator-grid, likelihood-rounding, solver, model-form and physical-population error. The defensive range bound may be much wider than descriptive standard error, especially for rare events. Zero failures do not prove zero risk. A resolved below-target result is a successful calculation, not approval. Stopped weighted means are not claimed unbiased; no temperature distribution or physical safety certificate is published.",
        config.qoi.no_claim(),
    );
    Ok(format!(
        "\"authority\":\"estimated-model-sampling-confidence\",\"status\":{},\"termination\":{},\"decision\":{},\"scope\":\"probability-of-declared-numerical-model\",\"failure_event\":\"qoi > temperature_limit_k\",\"required_probability\":{},\"maximum_failure_probability\":{},\"alpha\":{},\"min_decision_samples\":{},\"confidence_samples\":{},\"estimated_failure_probability\":{},\"failure_probability_confidence_sequence\":{},\"compliance_probability_confidence_sequence\":{},\"method\":\"defensive-weighted-gaussian-mixture-cs-sigma-half-rho-one\"{},\"no_claim\":{}",
        quote(if assessment.reached() { "decision-reached" } else { "inconclusive" }),
        quote(termination), quote(assessment.decision.label()),
        number_json(policy.required_probability)?, number_json(1.0 - policy.required_probability)?,
        number_json(policy.alpha)?, policy.min_samples, samples,
        optional_number(assessment.failure.as_ref().map(|x| x.mean))?,
        failure_interval, compliance_interval, recovery, quote(&no_claim),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn plan(samples: usize) -> UqPlan {
        UqPlan::new("temperature", PropagationMethod::MonteCarlo, samples)
            .with_parameter(ParameterUncertainty::gaussian("ambient", 300.0, 1.0, "K"))
            .with_correlation(CorrelationModel::Independent).with_compliance_threshold(310.0)
    }
    fn run(shift: f64) -> ImportanceExecution {
        ImportanceExecution::new(&plan(128), ImportanceConfig {
            germ_shift: vec![shift], defensive_fraction: 0.2,
        }).unwrap()
    }

    #[test]
    fn both_targets_need_confidence_and_the_production_minimum() {
        let policy = Policy::new(0.5, 0.05, 32).unwrap();
        assert!(!assess(policy, None).unwrap().reached());
        for (value, expected) in [(310.0, Decision::MeetsTarget), (320.0, Decision::BelowTarget)] {
            let mut execution = run(0.0);
            execution.advance(31, || false, |_| Ok::<_, &str>(value));
            assert!(!assess(policy, Some(&execution)).unwrap().reached());
            execution.advance(1, || false, |_| Ok::<_, &str>(value));
            let assessment = assess(policy, Some(&execution)).unwrap();
            assert_eq!(assessment.decision, expected);
            let original = execution.assess_failure_probability(0.05, 0.0).unwrap().unwrap();
            let actual = assessment.failure.unwrap();
            assert_eq!((actual.mean, actual.lo, actual.hi, actual.n),
                (original.mean, original.lo, original.hi, original.n));
        }
    }

    #[test]
    fn no_failures_and_zero_empirical_error_do_not_approve_a_tilted_run() {
        let mut execution = run(2.0);
        execution.advance(32, || false, |_| Ok::<_, &str>(300.0));
        assert_eq!(execution.report().estimate.unwrap().standard_error, Some(0.0));
        let assessment = assess(Policy::new(0.9, 0.05, 2).unwrap(), Some(&execution)).unwrap();
        assert_eq!(assessment.decision, Decision::Indeterminate);
        assert!(assessment.failure.unwrap().hi > 0.1);
        execution.advance(1, || false, |_| Err::<f64, _>("physical refusal"));
        assert!(assess(Policy::new(0.9, 0.05, 2).unwrap(), Some(&execution)).is_err());
    }

    #[test]
    fn stopping_policy_and_sampling_semantics_are_bound_before_recovery() {
        let policy = Policy::new(0.5, 0.05, 32).unwrap();
        let original = binding(policy, "[]");
        assert_ne!(original, policy.checkpoint_binding("[]"));
        for changed in [Policy::new(0.6, 0.05, 32).unwrap(),
            Policy::new(0.5, 0.01, 32).unwrap(), Policy::new(0.5, 0.05, 33).unwrap()] {
            assert_ne!(original, binding(changed, "[]"));
        }
        let flags: Vec<OsString> = ["--importance", "0.2,2", "--compliance-probability", "0.5",
            "--confidence-alpha", "0.05", "--min-decision-samples", "32",
            "--checkpoint", "new.imp", "--max-new-samples", "7"]
            .iter().map(OsString::from).collect();
        assert!(Options::parse(&flags).is_ok());
        assert!(Policy::new(0.5, 0.05, 32).unwrap().validate_plan(&plan(8)).is_err());
    }
}
