//! Learned rare-event proposals evaluated by the original cooling product.
//! Pilots select a proposal; only independent, fixed production estimates risk.
mod recovery;

use super::*;
use fs_uq::{
    ImportanceConfig, ImportanceFitConfig, PilotImportanceConfig,
    PilotImportanceExecution, PilotImportanceReport,
};

const USAGE: &str = "--importance-pilot requires ROUNDS,SAMPLES_PER_PILOT,DEFENSIVE_FRACTION,ELITE_FRACTION,RELAXATION,MAX_SHIFT_NORM";

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Settings {
    rounds: usize,
    samples_per_pilot: usize,
    defensive_fraction: f64,
    elite_fraction: f64,
    relaxation: f64,
    maximum_shift_norm: f64,
}

impl Settings {
    pub(super) fn parse(value: &OsString) -> Result<Self> {
        let text = value.to_str().filter(|text| text.len() <= 256)
            .ok_or_else(|| bad(USAGE))?;
        let parts: Vec<_> = text.split(',').collect();
        if parts.len() != 6 || parts.iter().any(|part| part.is_empty()) {
            return Err(bad(USAGE));
        }
        let rounds = count_option(&OsString::from(parts[0]), "importance pilot rounds")?;
        let samples_per_pilot = count_option(&OsString::from(parts[1]), "importance pilot samples")?;
        let defensive_fraction = float_option(&OsString::from(parts[2]), "importance defensive fraction")?;
        let elite_fraction = float_option(&OsString::from(parts[3]), "importance elite fraction")?;
        let relaxation = float_option(&OsString::from(parts[4]), "importance relaxation")?;
        let maximum_shift_norm = float_option(&OsString::from(parts[5]), "importance maximum shift norm")?;
        if !(1..=16).contains(&rounds) || samples_per_pilot < 2 {
            return Err(bad("importance learning requires 1..=16 pilots with at least two samples each"));
        }
        if defensive_fraction <= 0.0 || defensive_fraction > 1.0
            || !defensive_fraction.recip().is_finite()
            || elite_fraction <= 0.0 || elite_fraction > 1.0
            || relaxation <= 0.0 || relaxation > 1.0
            || maximum_shift_norm <= 0.0 || maximum_shift_norm > 1e6
        {
            return Err(bad("importance fractions and relaxation must lie in (0,1], defensive reciprocal must be finite, and maximum shift norm must lie in (0,1e6]"));
        }
        if rounds.checked_mul(samples_per_pilot).and_then(|n| n.checked_add(2))
            .is_none_or(|n| n > MAX_PRODUCT_SAMPLES)
        {
            return Err(bad("importance pilots must leave at least two production calls within the 10000-call ceiling"));
        }
        Ok(Self {
            rounds, samples_per_pilot, defensive_fraction, elite_fraction,
            relaxation, maximum_shift_norm,
        })
    }

    fn library_config(self, dimension: usize) -> PilotImportanceConfig {
        PilotImportanceConfig {
            pilot_rounds: self.rounds,
            samples_per_pilot: self.samples_per_pilot,
            initial_proposal: ImportanceConfig {
                germ_shift: vec![0.0; dimension],
                defensive_fraction: self.defensive_fraction,
            },
            fit: ImportanceFitConfig {
                elite_fraction: self.elite_fraction,
                relaxation: self.relaxation,
                maximum_shift_norm: self.maximum_shift_norm,
            },
        }
    }

    fn validate_decision(self, plan: &UqPlan, policy: Policy) -> Result<usize> {
        policy.validate_plan(plan)?;
        let pilots = self.rounds.checked_mul(self.samples_per_pilot)
            .ok_or_else(|| bad("importance pilot work overflows"))?;
        plan.budget_max_samples.checked_sub(pilots)
            .filter(|&n| n >= policy.min_samples)
            .ok_or_else(|| bad("--min-decision-samples must fit independent production AFTER every declared pilot"))
    }
}

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options, settings: Settings,
) -> Result<ExecutionOutput> {
    if options.checkpoint.is_some() || options.resume.is_some() || options.compliance.is_some() {
        return recovery::execute(worker, base_text, base, config, options, settings);
    }
    // Setup, every pilot, production, fitting and publication share ONE allowance.
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let plan = config.plan();
    let mut execution = PilotImportanceExecution::new(
        &plan, settings.library_config(plan.parameters.len()),
    ).map_err(bad)?;
    let report = execution.advance_interruptible(
        config.samples,
        || Instant::now() >= deadline,
        |values| super::evaluate(worker, base, config, values, deadline),
    );
    if report.status == UqStatus::Refused {
        return Err(model_failure(report.rejection_reason.clone()
            .unwrap_or_else(|| "importance pilot or production model refused".into())));
    }
    if report.status != UqStatus::Complete || Instant::now() >= deadline {
        return Err(budget(format!(
            "importance study stopped after {} of {} physical evaluations and {} completed pilots; no partial or pilot-pooled probability published; use --checkpoint to retain completed physics and learned proposals",
            report.model_calls_evaluated, report.model_calls_planned, report.completed_pilots.len(),
        )));
    }
    let stdout = render(base, config, settings, &report)?;
    if Instant::now() >= deadline {
        return Err(budget("importance result preparation exhausted the shared wall allowance; no estimate published"));
    }
    Ok(ExecutionOutput { stdout, exit_code: exit::SUCCESS })
}

fn pilots_json(report: &PilotImportanceReport) -> Result<String> {
    report.completed_pilots.iter().enumerate().map(|(round, pilot)| {
        let fit = &pilot.fit;
        Ok(format!(
            "{{\"round\":{},\"seed\":{},\"samples\":{},\"sampled_proposal\":{},\"selected\":{},\"observed_failures\":{},\"cutoff_k\":{},\"targets_failure_event\":{},\"effective_samples\":{},\"norm_before_limit\":{},\"norm_limited\":{},\"next_proposal\":{}}}",
            round, quote(&pilot.seed.to_string()), fit.samples, proposal_json(&pilot.sampled_proposal)?,
            fit.selected, fit.observed_failures, number_json(fit.cutoff)?, fit.targets_failure_event,
            number_json(fit.effective_samples)?, number_json(fit.norm_before_limit)?,
            fit.norm_limited, proposal_json(&fit.proposal)?,
        ))
    }).collect::<Result<Vec<_>>>().map(|rows| rows.join(","))
}

fn render(base: &J, config: &Config, settings: Settings, report: &PilotImportanceReport) -> Result<String> {
    if report.status != UqStatus::Complete {
        return Err(model_failure("only a complete fixed importance study may publish a probability"));
    }
    let production = report.production.as_ref()
        .filter(|production| production.status == UqStatus::Complete)
        .ok_or_else(|| model_failure("missing completed independent importance production"))?;
    let estimate = report.estimate.as_ref().ok_or_else(|| model_failure(
        production.unavailable_reason.unwrap_or("completed importance statistic is unavailable"),
    ))?;
    let pilots = pilots_json(report)?;
    let no_claim = format!(
        "{} The strict event is QoI > temperature_limit_k. Only independently keyed final production contributes target/proposal-weighted failure indicators. Pilot failure counts and proposal failure frequency are not target probabilities. The estimate is neither self-normalized nor clipped to [0,1]. Standard error and effective sample count are descriptive, not confidence bounds or stopping rules. Zero observed failures and zero empirical error do not prove zero risk. Complete means the fixed call design finished, not achieved accuracy, variance reduction or physical approval. A single shift may miss disconnected failure modes. Finite-generator, floating-point, solver and physical-model errors are not bounded.",
        config.qoi.no_claim(),
    );
    Ok(format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.importance-pilot.v1\",\"authority\":\"estimated-model-importance-sampling\",\"method\":\"pilot-learned-defensive-gaussian\",\"status\":\"complete\",\"qoi\":{},\"seed\":{},\"production_seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"pilot_rounds\":{},\"samples_per_pilot\":{},\"pilot_samples_evaluated\":{},\"production_samples\":{},\"temperature_limit_k\":{},\"failure_event\":\"qoi>temperature_limit_k\",\"estimated_failure_probability\":{},\"failure_probability_standard_error\":{},\"failure_effective_samples\":{},\"mean_likelihood_weight\":{},\"production_proposal_failures\":{},\"production_proposal\":{},\"fit_policy\":{{\"elite_fraction\":{},\"relaxation\":{},\"maximum_shift_norm\":{}}},\"pilots\":[{}],\"correlation\":{},\"parameters\":[{}],\"no_claim\":{}}}\n",
        config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
        quote(&report.current_seed.to_string()), report.model_calls_planned, report.model_calls_evaluated,
        settings.rounds, settings.samples_per_pilot,
        report.model_calls_evaluated - production.evaluations_attempted, production.evaluations_attempted,
        number_json(production.threshold)?, number_json(estimate.failure_probability)?,
        optional_number(estimate.standard_error)?, optional_number(estimate.failure_effective_samples)?,
        number_json(estimate.mean_likelihood_weight)?, production.proposal_failures,
        proposal_json(&report.current_proposal)?, number_json(settings.elite_fraction)?,
        number_json(settings.relaxation)?, number_json(settings.maximum_shift_norm)?, pilots,
        quote(config.correlation_label), config.render_parameters()?, quote(&no_claim),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }
    const VALID: &str = "2,8,0.2,0.25,1,8";

    #[test]
    fn all_pilot_choices_are_explicit_and_numerically_admitted() {
        let settings = Settings::parse(&OsString::from(VALID)).unwrap();
        assert_eq!(settings.rounds, 2);
        assert_eq!(settings.samples_per_pilot, 8);
        let config = settings.library_config(3);
        assert_eq!(config.initial_proposal.germ_shift, vec![0.0; 3]);
        assert_eq!(config.initial_proposal.defensive_fraction, 0.2);
        assert_eq!(config.fit.elite_fraction, 0.25);
        for value in [
            "", "2,8,0.2", "2,8,0.2,0.25,1,8,9", "2,,0.2,0.25,1,8",
            "0,8,0.2,0.25,1,8", "17,8,0.2,0.25,1,8", "2,1,0.2,0.25,1,8",
            "2.0,8,0.2,0.25,1,8", "2,-8,0.2,0.25,1,8", "2,5000,0.2,0.25,1,8",
            "2,8,0,0.25,1,8", "2,8,1.1,0.25,1,8", "2,8,1e-320,0.25,1,8",
            "2,8,NaN,0.25,1,8", "2,8,0.2,0,1,8", "2,8,0.2,1.1,1,8",
            "2,8,0.2,0.25,0,8", "2,8,0.2,0.25,1.1,8", "2,8,0.2,0.25,1,0",
            "2,8,0.2,0.25,1,1000001", "2,8,0.2,0.25,1,inf",
        ] {
            assert!(Settings::parse(&OsString::from(value)).is_err(), "accepted {value}");
        }
    }

    #[test]
    fn pilot_mode_cannot_silently_override_another_method() {
        assert!(Options::parse(&args(&["--importance-pilot", VALID])).unwrap().importance.is_some());
        for other in [
            vec!["--importance-pilot", VALID], vec!["--importance", "0.2,0"], vec!["--pce", "0,2"],
            vec!["--pce-mean", "0,2,4"], vec!["--qmc-replicates", "2"],
            vec!["--sensitivity", "sobol"], vec!["--sensitivity", "shapley"],
            vec!["--mean-control", "adjoint"], vec!["--max-new-samples", "0"],
        ] {
            let own = args(&["--importance-pilot", VALID]);
            let extra = args(&other);
            for values in [
                own.iter().chain(&extra).cloned().collect::<Vec<_>>(),
                extra.iter().chain(&own).cloned().collect::<Vec<_>>(),
            ] {
                assert!(Options::parse(&values).is_err(), "accepted {values:?}");
            }
        }
    }

    #[test]
    fn the_original_budget_must_cover_every_pilot_and_independent_production() {
        let settings = Settings::parse(&OsString::from(VALID)).unwrap();
        for count in [2, 15, 16, 17, 18] {
            let plan = fs_uq::UqPlan::new("temperature", PropagationMethod::MonteCarlo, count)
                .with_parameter(ParameterUncertainty::gaussian("inlet", 300.0, 0.2, "K"))
                .with_correlation(CorrelationModel::Independent)
                .with_compliance_threshold(310.0);
            let result = PilotImportanceExecution::new(&plan, settings.library_config(1));
            assert_eq!(result.is_ok(), count == 18);
        }
    }

    #[test]
    fn decision_minimum_cannot_be_satisfied_by_learning_observations() {
        let settings = Settings::parse(&OsString::from(VALID)).unwrap();
        let policy = Policy::new(0.5, 0.05, 16).unwrap();
        for count in [16, 17, 18, 31, 32] {
            let plan = UqPlan::new("temperature", PropagationMethod::MonteCarlo, count)
                .with_parameter(ParameterUncertainty::gaussian("inlet", 300.0, 0.2, "K"))
                .with_correlation(CorrelationModel::Independent)
                .with_compliance_threshold(310.0);
            assert_eq!(settings.validate_decision(&plan, policy).is_ok(), count == 32);
        }
        let values = args(&["--importance-pilot", VALID, "--compliance-probability", "0.5",
            "--confidence-alpha", "0.05", "--min-decision-samples", "16"]);
        assert!(Options::parse(&values).is_ok());
    }

    #[test]
    fn learned_proposal_confidence_reuses_only_the_original_production_owner() {
        let config = PilotImportanceConfig {
            pilot_rounds: 2, samples_per_pilot: 4,
            initial_proposal: ImportanceConfig { germ_shift: vec![0.0], defensive_fraction: 0.5 },
            fit: ImportanceFitConfig {
                elite_fraction: 0.5, relaxation: 1.0, maximum_shift_norm: 6.0,
            },
        };
        let plan = UqPlan::new("temperature", PropagationMethod::MonteCarlo, 32)
            .with_parameter(ParameterUncertainty::gaussian("inlet", 300.0, 1.0, "K"))
            .with_correlation(CorrelationModel::Independent).with_compliance_threshold(310.0);
        let mut learned = PilotImportanceExecution::new(&plan, config).unwrap();
        let policy = Policy::new(0.5, 0.05, 8).unwrap();
        learned.advance(7, || false, |x| Ok::<_, &str>(x[0]));
        assert!(decision::assess(policy, learned.production()).unwrap().failure.is_none());
        learned.advance(5, || false, |x| Ok::<_, &str>(x[0]));
        assert_eq!(learned.evaluations_attempted(), 12);
        let assessment = decision::assess(policy, learned.production()).unwrap();
        assert!(!assessment.reached(), "eight pilots are not eight production observations");
        let actual = assessment.failure.unwrap();
        let expected = learned.production().unwrap().assess_failure_probability(0.05, 0.0).unwrap().unwrap();
        assert_eq!(actual.n, 4);
        assert_eq!((actual.mean, actual.lo, actual.hi), (expected.mean, expected.lo, expected.hi));
    }
}
