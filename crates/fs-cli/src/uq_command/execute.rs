use super::*;
use child::{EvaluationError, evaluate_sample};
use compliance::Policy;
use model::Config;
use std::path::PathBuf;
use std::time::{Duration, Instant};

#[path = "checkpoint.rs"]
mod checkpoint;
#[path = "compliance.rs"]
mod compliance;

#[derive(Debug, Default)]
pub(super) struct Options {
    checkpoint: Option<PathBuf>,
    resume: Option<PathBuf>,
    max_new_samples: Option<usize>,
    compliance: Option<Policy>,
}

impl Options {
    pub(super) fn parse(args: &[OsString]) -> Result<Self> {
        let mut options = Self::default();
        let mut probability = None;
        let mut alpha = None;
        let mut minimum = None;
        let mut index = 0;
        while index < args.len() {
            let flag = &args[index];
            let value = args.get(index + 1).ok_or_else(|| bad("each UQ execution option requires a value"))?;
            if flag == "--checkpoint" && options.checkpoint.is_none() {
                options.checkpoint = Some(PathBuf::from(value.as_os_str()));
            } else if flag == "--resume" && options.resume.is_none() {
                options.resume = Some(PathBuf::from(value.as_os_str()));
            } else if flag == "--max-new-samples" && options.max_new_samples.is_none() {
                options.max_new_samples = Some(count_option(value, "--max-new-samples")?);
            } else if flag == "--compliance-probability" && probability.is_none() {
                probability = Some(float_option(value, "--compliance-probability")?);
            } else if flag == "--confidence-alpha" && alpha.is_none() {
                alpha = Some(float_option(value, "--confidence-alpha")?);
            } else if flag == "--min-decision-samples" && minimum.is_none() {
                minimum = Some(count_option(value, "--min-decision-samples")?);
            } else {
                return Err(bad(format!("unknown or duplicate UQ option {}", flag.to_string_lossy())));
            }
            index += 2;
        }
        options.compliance = match (probability, alpha, minimum) {
            (None, None, None) => None,
            (Some(probability), Some(alpha), Some(minimum)) => Some(Policy::new(probability, alpha, minimum)?),
            _ => return Err(bad("sequential compliance requires --compliance-probability, --confidence-alpha and --min-decision-samples together")),
        };
        if options.max_new_samples.is_some() && options.checkpoint.is_none() {
            return Err(bad("--max-new-samples requires --checkpoint so the unfinished prefix is retained"));
        }
        Ok(options)
    }

    fn checkpoint_binding(&self, config: &Config) -> Result<String> {
        let parameters = config.render_parameters()?;
        // Preserve the existing fixed-count identity exactly. Adding, removing
        // or changing a sequential policy yields a different identity and is
        // refused BEFORE reserving an output or evaluating a cooling sample.
        Ok(self.compliance.map_or_else(
            || parameters.clone(), |policy| policy.checkpoint_binding(&parameters),
        ))
    }
}

fn count_option(value: &OsString, flag: &str) -> Result<usize> {
    let value = value.to_str().and_then(|text| {
        (!text.is_empty() && text.bytes().all(|b| b.is_ascii_digit()))
            .then(|| text.parse::<usize>().ok()).flatten()
    }).ok_or_else(|| bad(format!("{flag} must be a nonnegative decimal integer")))?;
    if value > MAX_PRODUCT_SAMPLES {
        return Err(bad(format!("{flag} exceeds {MAX_PRODUCT_SAMPLES}")));
    }
    Ok(value)
}

fn float_option(value: &OsString, flag: &str) -> Result<f64> {
    value.to_str().and_then(|text| text.parse::<f64>().ok())
        .filter(|value| value.is_finite())
        .ok_or_else(|| bad(format!("{flag} must be a finite number")))
}

pub(super) struct ExecutionOutput {
    pub(super) stdout: String,
    pub(super) exit_code: u8,
}

// Preserve the original two-argument caller and its no-partial-output behavior.
pub(super) fn execute(base_text: &str, uq_text: &str) -> Result<String> {
    execute_with_options(base_text, uq_text, &Options::default()).map(|output| output.stdout)
}

pub(super) fn execute_with_options(base_text: &str, uq_text: &str, options: &Options) -> Result<ExecutionOutput> {
    let base = J::parse(base_text).map_err(|error| bad(format!("invalid base JSON: {error}")))?;
    let config = Config::parse(uq_text, &base)?;
    let plan = config.plan();
    if let Some(policy) = options.compliance {
        policy.validate_plan(&plan)?;
    }
    let mut execution = UqExecution::new(&plan).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(base_text, &options.checkpoint_binding(&config)?)?)
    } else { None };
    if let Some(path) = &options.resume {
        execution = checkpoint::restore(path, &plan, identity.expect("resume requested identity"))?;
    }
    // Restore/admit BEFORE reserving output. Existing paths, including the
    // resume source itself, are never overwritten by this invocation.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    if let Some(output) = &output {
        output.save(&execution, identity.expect("output requested identity"))?;
    }
    let initial_count = execution.observations().len();
    let allowance = options.max_new_samples.unwrap_or(config.samples)
        .min(config.samples - initial_count);
    let mut assessment = options.compliance
        .map(|policy| policy.assess(&execution)).transpose()?;
    // This is a fresh evaluation-time allowance per invocation. It does not
    // reset the immutable lifetime sample budget retained in the checkpoint.
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    for _ in 0..allowance {
        // A checkpoint already at its stopping ordinal must not launch another
        // child or alter the accepted prefix when resumed.
        if assessment.as_ref().is_some_and(|value| value.reached()) { break; }
        let report = execution.advance_interruptible(1, || Instant::now() >= deadline, |values| {
            let request = config.sample_request(&base, values)?;
            match evaluate_sample(&request, deadline) {
                Ok(value) => Ok(Some(value)),
                Err(EvaluationError::Budget) => Ok(None),
                Err(EvaluationError::Child(message)) => Err(model_failure(message)),
            }
        });
        if report.status == UqStatus::Refused {
            let refusal = Failure {
                code: "cooling-network-uq-refused",
                message: report.rejection_reason.clone().unwrap_or_else(|| "model evaluation refused".into()),
            };
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}"));
                    error
                })?;
            }
            return Err(refusal);
        }
        if let Some(output) = &output {
            output.save(&execution, identity.expect("output requested identity"))?;
        }
        assessment = options.compliance
            .map(|policy| policy.assess(&execution)).transpose()?;
        if matches!(report.status, UqStatus::Complete | UqStatus::Cancelled) { break; }
    }
    let report = execution.report();
    if let (Some(policy), Some(assessment)) = (options.compliance, assessment) {
        let termination = if assessment.reached() {
            "probability-target"
        } else if report.status == UqStatus::Complete {
            "sample-budget"
        } else if report.status == UqStatus::Cancelled {
            "wall-time-budget"
        } else {
            "sample-chunk"
        };
        let stdout = policy.render(
            &config, objective_kind(&base)?, &execution, &assessment,
            termination, options.checkpoint.as_deref(),
        )?;
        // A below-target decision is a successful statistical computation,
        // NOT an engineering approval. An unresolved target never exits zero,
        // even after every planned sample was evaluated.
        let exit_code = if assessment.reached() { exit::SUCCESS } else { exit::BUDGET };
        return Ok(ExecutionOutput { stdout, exit_code });
    }
    if report.status == UqStatus::Complete {
        return Ok(ExecutionOutput { stdout: render_result(&config, &base, &report)?, exit_code: exit::SUCCESS });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!("UQ wall-time budget exhausted after {} model evaluations; no partial distribution published; use --checkpoint to retain completed samples", report.samples_evaluated)));
    };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"next_sample_ordinal\":{},\"checkpoint\":{},\"no_claim\":\"retained prefix only; no completed distribution or compliance decision; resume with the same base request, UQ plan and executable on the same deterministic runtime profile\"}}\n",
        quote(termination), config.samples, report.samples_evaluated,
        report.samples_evaluated - initial_count, report.samples_evaluated,
        quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}

fn objective_kind(base: &J) -> Result<&'static str> {
    let objective = field(base, "objective")?;
    Ok(if objective.get("mean_wall_region").is_some() {
        "mean_wall_region"
    } else if objective.get("max_wall_region").is_some() {
        "max_wall_region"
    } else if objective.get("max_solid_temperature").is_some() {
        "max_solid_temperature"
    } else { "max_vertices" })
}

fn render_result(config: &Config, base: &J, result: &fs_uq::UqResult) -> Result<String> {
    let percentiles = match result.percentiles {
        Some(values) => format!("[{},{},{}]", number_json(values[0])?, number_json(values[1])?, number_json(values[2])?),
        None => "null".into(),
    };
    let bounds = if result.mean.is_some() {
        format!("[{},{}]", number_json(result.interval_bounds[0])?, number_json(result.interval_bounds[1])?)
    } else { "null".into() };
    let objective_kind = objective_kind(base)?;
    Ok(format!(
        "{{\"schema\":{},\"authority\":\"estimated-empirical-monte-carlo\",\"status\":{},\"qoi\":{{\"kind\":{},\"unit\":\"K\"}},\"seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"mean_k\":{},\"std_dev_k\":{},\"percentiles_p05_p50_p95_k\":{},\"empirical_bounds_k\":{},\"temperature_limit_k\":{},\"empirical_probability_of_compliance\":{},\"sampling_standard_error_k\":{},\"correlation\":{},\"parameters\":[{}],\"no_claim\":\"fixed-count empirical propagation through actual cooling-network child solves; no confidence sequence, optional-stopping guarantee, physical/model-form uncertainty bound, mesh-convergence certificate, experimental validation, or native .fsim/ledger package claim\"}}\n",
        quote(RESULT_SCHEMA), quote(result.status.label()), quote(objective_kind), quote(&config.seed.to_string()),
        config.samples, result.samples_evaluated, optional_number(result.mean)?, optional_number(result.std_dev)?,
        percentiles, bounds, optional_number(config.threshold_k)?, optional_number(result.probability_of_compliance)?,
        number_json(result.sampling_error)?, quote(config.correlation_label), config.render_parameters()?
    ))
}

#[cfg(test)]
mod options_tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }

    #[test]
    fn options_require_recoverable_chunks_and_reject_ambiguous_flags() {
        assert!(Options::parse(&args(&["--checkpoint", "a.bin", "--max-new-samples", "0"])).is_ok());
        assert!(Options::parse(&args(&["--resume", "a.bin", "--checkpoint", "b.bin", "--max-new-samples", "3"])).is_ok());
        for values in [
            vec!["--max-new-samples", "3"],
            vec!["--checkpoint"],
            vec!["--unknown", "a.bin"],
            vec!["--checkpoint", "a.bin", "--checkpoint", "b.bin"],
            vec!["--checkpoint", "a.bin", "--max-new-samples", "-1"],
            vec!["--checkpoint", "a.bin", "--max-new-samples", "10001"],
        ] {
            assert!(Options::parse(&args(&values)).is_err(), "accepted {values:?}");
        }
    }

    #[test]
    fn sequential_options_are_an_explicit_all_or_nothing_policy() {
        let complete = ["--compliance-probability", "0.95", "--confidence-alpha", "0.05", "--min-decision-samples", "32"];
        let options = Options::parse(&args(&complete)).unwrap();
        assert_eq!(options.compliance, Some(Policy::new(0.95, 0.05, 32).unwrap()));
        for values in [
            vec!["--compliance-probability", "0.95"],
            vec!["--confidence-alpha", "0.05", "--min-decision-samples", "32"],
            vec!["--compliance-probability", "NaN", "--confidence-alpha", "0.05", "--min-decision-samples", "32"],
            vec!["--compliance-probability", "0.95", "--confidence-alpha", "1", "--min-decision-samples", "32"],
            vec!["--compliance-probability", "0.95", "--confidence-alpha", "0.05", "--min-decision-samples", "1"],
        ] {
            assert!(Options::parse(&args(&values)).is_err(), "accepted {values:?}");
        }
        let mut duplicate = args(&complete);
        duplicate.extend(args(&["--confidence-alpha", "0.01"]));
        assert!(Options::parse(&duplicate).is_err());
    }
}
