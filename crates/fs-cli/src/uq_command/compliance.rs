//! Sequential decisions about P(numerical cooling QoI <= declared ceiling).
//! The statistical producer is fs-uq::UqExecution::assess_compliance; this
//! adapter neither invents a second bound nor uses Monte Carlo standard error.

use super::{Config, Failure, Result, UqExecution, UqPlan, bad, number_json, optional_number, quote};
use fs_uq::AnytimeEstimate;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Policy {
    pub(super) required_probability: f64,
    pub(super) alpha: f64,
    pub(super) min_samples: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Decision {
    MeetsTarget,
    BelowTarget,
    Indeterminate,
}

impl Decision {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::MeetsTarget => "meets-probability-target",
            Self::BelowTarget => "below-probability-target",
            Self::Indeterminate => "indeterminate",
        }
    }
}

pub(super) struct Assessment {
    pub(super) estimate: Option<AnytimeEstimate>,
    pub(super) decision: Decision,
}

impl Assessment {
    pub(super) fn reached(&self) -> bool {
        self.decision != Decision::Indeterminate
    }
}

impl Policy {
    pub(super) fn new(required_probability: f64, alpha: f64, min_samples: usize) -> Result<Self> {
        if !required_probability.is_finite() || required_probability <= 0.0 || required_probability >= 1.0 {
            return Err(bad("--compliance-probability must be strictly inside (0,1)"));
        }
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 || !(1.0 / alpha).is_finite() {
            return Err(bad("--confidence-alpha must be strictly inside (0,1) with a finite reciprocal"));
        }
        if !(2..=super::MAX_PRODUCT_SAMPLES).contains(&min_samples) {
            return Err(bad("--min-decision-samples must be between two and the product sample cap"));
        }
        Ok(Self { required_probability, alpha, min_samples })
    }

    pub(super) fn validate_plan(&self, plan: &UqPlan) -> Result<()> {
        if plan.compliance_threshold.is_none() {
            return Err(bad("sequential compliance requires temperature_limit_k in the UQ request"));
        }
        if self.min_samples > plan.budget_max_samples {
            return Err(bad("--min-decision-samples exceeds the original UQ sample budget"));
        }
        Ok(())
    }

    pub(super) fn assess(&self, execution: &UqExecution) -> Result<Assessment> {
        self.validate_plan(execution.plan())?;
        // Width is not a stop rule here. Only separation from the declared
        // probability target licenses a decision; the empirical mean alone
        // and clipping at zero/one cannot manufacture one.
        let estimate = execution.assess_compliance(self.alpha, 0.0).map_err(|error| Failure {
            code: "cooling-network-uq-confidence",
            message: error.to_string(),
        })?;
        let decision = match estimate.as_ref() {
            Some(interval) if interval.n >= self.min_samples as u64 => {
                if interval.lo >= self.required_probability {
                    Decision::MeetsTarget
                } else if interval.hi < self.required_probability {
                    Decision::BelowTarget
                } else {
                    Decision::Indeterminate
                }
            }
            _ => Decision::Indeterminate,
        };
        Ok(Assessment { estimate, decision })
    }

    /// Versioned, length-unambiguous identity input, not just a display label.
    /// Base model, executable, threshold, sampler and distributions are also
    /// bound by the existing file adapter and fs-uq checkpoint envelope.
    pub(super) fn checkpoint_binding(&self, parameters: &str) -> String {
        format!(
            "{{\"parameters\":[{parameters}],\"compliance_policy\":{{\"version\":1,\"method\":\"gaussian-mixture-cs-sigma-half-rho-one\",\"probability_bits\":\"{:016x}\",\"alpha_bits\":\"{:016x}\",\"min_samples\":{}}}}}",
            self.required_probability.to_bits(), self.alpha.to_bits(), self.min_samples,
        )
    }

    pub(super) fn render(
        &self,
        config: &Config,
        qoi_kind: &str,
        execution: &UqExecution,
        assessment: &Assessment,
        termination: &str,
        checkpoint: Option<&Path>,
    ) -> Result<String> {
        let probability_interval = match assessment.estimate.as_ref() {
            Some(interval) => format!("[{},{}]", number_json(interval.lo)?, number_json(interval.hi)?),
            None => "null".to_string(),
        };
        // A completed decision has no invocation-specific fields, so resuming
        // a checkpoint at the stopping ordinal reproduces the same result bytes.
        let recovery = if assessment.reached() {
            String::new()
        } else {
            format!(",\"checkpoint\":{}", checkpoint.map_or_else(
                || "null".to_string(), |path| quote(&path.to_string_lossy()),
            ))
        };
        Ok(format!(
            "{{\"schema\":\"frankensim.cooling-network-uq.compliance.v1\",\"authority\":\"estimated-model-sampling-confidence\",\"status\":{},\"termination\":{},\"decision\":{},\"scope\":\"probability-of-declared-numerical-model\",\"qoi\":{{\"kind\":{},\"unit\":\"K\"}},\"seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"temperature_limit_k\":{},\"required_probability\":{},\"alpha\":{},\"min_decision_samples\":{},\"empirical_probability_of_compliance\":{},\"probability_confidence_sequence\":{},\"method\":\"gaussian-mixture-cs-sigma-half-rho-one\",\"correlation\":{},\"parameters\":[{}]{},\"no_claim\":\"sampling confidence for one fixed model, probability law, threshold and predeclared alpha under the confidence sequence's fixed-conditional-mean assumptions; no multiplicity control across models, seeds or policies; numerical bounds are not outward-rounded; no physical/model-form uncertainty bound, mesh certificate, experimental validation or physical safety certification; stopped temperature means and quantiles are not confidence bounds\"}}\n",
            quote(if assessment.reached() { "decision-reached" } else { "inconclusive" }),
            quote(termination), quote(assessment.decision.label()), quote(qoi_kind), quote(&config.seed.to_string()),
            config.samples, execution.observations().len(), optional_number(config.threshold_k)?,
            number_json(self.required_probability)?, number_json(self.alpha)?, self.min_samples,
            optional_number(assessment.estimate.as_ref().map(|interval| interval.mean))?,
            probability_interval, quote(config.correlation_label), config.render_parameters()?, recovery,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_blake3::ContentHash;
    use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod};

    fn plan(samples: usize) -> UqPlan {
        UqPlan::new("temperature", PropagationMethod::MonteCarlo, samples)
            .with_parameter(ParameterUncertainty::uniform("ambient", 300.0, 300.0, "K"))
            .with_correlation(CorrelationModel::Independent)
            .with_compliance_threshold(310.0)
    }

    #[test]
    fn both_decisions_require_confidence_and_the_declared_minimum() {
        let policy = Policy::new(0.5, 0.05, 32).unwrap();
        for (value, expected) in [(300.0, Decision::MeetsTarget), (320.0, Decision::BelowTarget)] {
            let mut execution = UqExecution::new(&plan(128)).unwrap();
            assert!(!policy.assess(&execution).unwrap().reached());
            execution.advance(31, || false, |_| Ok::<_, &str>(value));
            assert!(!policy.assess(&execution).unwrap().reached());
            execution.advance(1, || false, |_| Ok::<_, &str>(value));
            assert_eq!(policy.assess(&execution).unwrap().decision, expected);
        }
    }

    #[test]
    fn probability_at_the_target_is_not_decided_from_its_empirical_mean() {
        let policy = Policy::new(0.5, 0.05, 2).unwrap();
        let mut execution = UqExecution::new(&plan(128)).unwrap();
        for i in 0..128 {
            execution.advance(1, || false, |_| Ok::<_, &str>(if i % 2 == 0 { 300.0 } else { 320.0 }));
            assert_eq!(policy.assess(&execution).unwrap().decision, Decision::Indeterminate);
        }
    }

    #[test]
    fn restored_prefix_reaches_the_same_first_decision_and_observation_bits() {
        let policy = Policy::new(0.5, 0.05, 32).unwrap();
        let plan = plan(128);
        let model = ContentHash([17; 32]);
        let mut full = UqExecution::new(&plan).unwrap();
        while !policy.assess(&full).unwrap().reached() {
            full.advance(1, || false, |values| Ok::<_, &str>(values[0]));
        }
        let mut prefix = UqExecution::new(&plan).unwrap();
        prefix.advance(7, || false, |values| Ok::<_, &str>(values[0]));
        let mut restored = UqExecution::restore(&plan, model, &prefix.checkpoint(model).unwrap()).unwrap();
        let mut calls = 0;
        while !policy.assess(&restored).unwrap().reached() {
            restored.advance(1, || false, |values| { calls += 1; Ok::<_, &str>(values[0]) });
        }
        assert_eq!(calls, full.observations().len() - 7);
        assert_eq!(full.checkpoint(model).unwrap(), restored.checkpoint(model).unwrap());
    }

    #[test]
    fn failed_samples_cannot_leave_a_successful_compliance_prefix() {
        let policy = Policy::new(0.5, 0.05, 2).unwrap();
        let mut execution = UqExecution::new(&plan(128)).unwrap();
        execution.advance(40, || false, |_| Ok::<_, &str>(300.0));
        execution.advance(1, || false, |_| Err::<f64, _>("physical solver refused"));
        assert!(policy.assess(&execution).is_err());
    }

    #[test]
    fn admission_and_identity_bind_every_stopping_choice() {
        for invalid in [0.0, 1.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(Policy::new(invalid, 0.05, 2).is_err());
            assert!(Policy::new(0.5, invalid, 2).is_err());
        }
        assert!(Policy::new(0.5, f64::from_bits(1), 2).is_err());
        assert!(Policy::new(0.5, 0.05, 1).is_err());
        let policy = Policy::new(0.5, 0.05, 32).unwrap();
        assert!(policy.validate_plan(&plan(8)).is_err());
        let mut missing = plan(128); missing.compliance_threshold = None;
        assert!(policy.validate_plan(&missing).is_err());
        let identity = policy.checkpoint_binding("{}");
        for changed in [
            Policy::new(0.6, 0.05, 32).unwrap(),
            Policy::new(0.5, 0.01, 32).unwrap(),
            Policy::new(0.5, 0.05, 33).unwrap(),
        ] {
            assert_ne!(identity, changed.checkpoint_binding("{}"));
        }
    }
}
