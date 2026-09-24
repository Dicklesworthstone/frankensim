//! Direct physical-input variance sensitivity through the existing cooling child.
//! The sampler owns A/B/hybrid pairing; every observation is a complete solve.
use super::*;
use fs_uq::{SobolExecution, SobolReport};

pub(super) fn execute(
    base_text: &str, base: &J, config: &Config, options: &Options,
) -> Result<ExecutionOutput> {
    let plan = config.plan();
    let mut execution = SobolExecution::new(&plan).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(base_text, &config.render_parameters()?)?)
    } else { None };
    if let Some(path) = &options.resume {
        execution = checkpoint::restore_sobol(path, &plan, identity.expect("resume identity"))?;
    }
    // Re-admit the exact model and layout before reserving a new output.
    // Preserve the source checkpoint and publish each paid evaluation atomically.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    if let Some(output) = &output {
        output.save_sobol(&execution, identity.expect("output identity"))?;
    }
    let initial_count = execution.observations().len();
    let allowance = options.max_new_samples.unwrap_or(config.samples)
        .min(config.samples - initial_count);
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    for _ in 0..allowance {
        let report = execution.advance_interruptible(1, || Instant::now() >= deadline, |values| {
            let request = config.sample_request(base, values)?;
            let evaluated = match config.qoi {
                Qoi::Steady => evaluate_sample(&request, deadline),
                qoi => evaluate_sample_for(&request, deadline, qoi),
            };
            match evaluated {
                Ok(value) => Ok(Some(value)),
                Err(EvaluationError::Budget) => Ok(None),
                Err(EvaluationError::Child(message)) => Err(model_failure(message)),
            }
        });
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason.clone()
                .unwrap_or_else(|| "Sobol physical evaluation refused".into()));
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}"));
                    error
                })?;
            }
            return Err(refusal);
        }
        if let Some(output) = &output {
            output.save_sobol(&execution, identity.expect("output identity"))?;
        }
        if matches!(report.status, UqStatus::Complete | UqStatus::Cancelled) { break; }
    }
    let report = execution.report();
    if report.status == UqStatus::Complete {
        return Ok(ExecutionOutput { stdout: render(base, config, &report)?, exit_code: exit::SUCCESS });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!(
            "Sobol sensitivity stopped after {} accepted evaluations and {} complete rows; no shortened-design sensitivity published; use --checkpoint to retain completed evaluations",
            report.evaluations_accepted, report.completed_rows,
        )));
    };
    let width = plan.parameters.len() + 2;
    let slot = report.evaluations_accepted % width;
    let next_kind = match slot { 0 => "base-a", 1 => "base-b", _ => "hybrid" };
    let parameter = if slot < 2 { "null".into() }
        else { quote(&plan.parameters[slot - 2].name) };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.sensitivity.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"base_samples\":{},\"evaluations_per_row\":{},\"completed_rows\":{},\"next_sample_ordinal\":{},\"next_row_ordinal\":{},\"next_row_slot\":{},\"next_evaluation_kind\":{},\"next_hybrid_parameter\":{},\"checkpoint\":{},\"no_claim\":\"retained direct-model Sobol evaluation prefix, including incomplete A/B/hybrid rows; no sensitivity indices, shortened-design statistics or optional-stopping inference; resume the identical base, plan, executable and fixed pairing\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?),
        report.evaluations_planned, report.evaluations_attempted,
        report.evaluations_attempted - initial_count, report.base_samples, width,
        report.completed_rows, report.evaluations_accepted, report.completed_rows,
        slot, quote(next_kind), parameter, quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}

fn render(base: &J, config: &Config, report: &SobolReport) -> Result<String> {
    if report.status != UqStatus::Complete {
        return Err(model_failure("only complete fixed-design Sobol reports may be published"));
    }
    let estimate = report.estimate.as_ref().ok_or_else(|| model_failure(format!(
        "Sobol indices unavailable after {} evaluations: {}",
        report.evaluations_accepted, report.unavailable_reason.unwrap_or("missing complete estimate"),
    )))?;
    let effects = estimate.effects.iter().map(|effect| Ok(format!(
        "{{\"parameter\":{},\"unit\":{},\"first_order\":{},\"total_order\":{}}}",
        quote(&effect.parameter), quote(&effect.unit),
        number_json(effect.first_order)?, number_json(effect.total_order)?,
    ))).collect::<Result<Vec<_>>>()?.join(",");
    let no_claim = format!(
        "{} Direct-model variance sensitivity under declared independent inputs, not surrogate coefficients or local derivatives. Main and total effects are noisy fixed-design estimates, not confidence intervals or causal effects. Values are not clipped to [0,1]; total effects can sum above one because interactions overlap. No convergence, continuous-time peak bound, model validation, compliance decision or optional-stopping guarantee is inferred. A temperature limit does not change this observable into a compliance indicator.",
        config.qoi.no_claim(),
    );
    Ok(format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.sensitivity.v1\",\"authority\":\"estimated-direct-model-sobol\",\"method\":\"jansen-pick-freeze-philox\",\"status\":\"complete\",\"qoi\":{},\"seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"base_samples\":{},\"evaluations_per_row\":{},\"base_output_std_dev_k\":{},\"correlation\":{},\"parameters\":[{}],\"effects\":[{}],\"no_claim\":{}}}\n",
        config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
        report.evaluations_planned, report.evaluations_attempted, report.base_samples,
        estimate.effects.len() + 2, optional_number(estimate.base_output_std_dev)?,
        quote(config.correlation_label), config.render_parameters()?, effects, quote(&no_claim),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }
    #[test]
    fn sensitivity_is_explicit_and_rejects_incompatible_execution_semantics() {
        assert!(Options::parse(&args(&["--sensitivity", "sobol"])).unwrap().sobol_sensitivity);
        assert!(!Options::parse(&[]).unwrap().sobol_sensitivity);
        assert!(Options::parse(&args(&["--sensitivity", "sobol", "--checkpoint", "a.bin", "--max-new-samples", "3"])).is_ok());
        assert!(Options::parse(&args(&["--sensitivity", "sobol", "--resume", "a.bin", "--checkpoint", "b.bin"])).is_ok());
        for tail in [
            vec!["--sensitivity", "sobol"], vec!["--qmc-replicates", "2"],
            vec!["--max-new-samples", "3"],
            vec!["--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "16"],
        ] {
            let mut values = args(&["--sensitivity", "sobol"]); values.extend(args(&tail));
            assert!(Options::parse(&values).is_err(), "accepted {tail:?}");
        }
        assert!(Options::parse(&args(&["--sensitivity", "unknown"])).is_err());
    }
}
