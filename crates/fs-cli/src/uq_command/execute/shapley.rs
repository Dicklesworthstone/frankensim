//! Physical-input attribution through the original cooling producer.
//! Conditional Gaussian dependence or explicitly independent mixed marginals.
use super::*;
use fs_uq::{ShapleyExecution, ShapleyReport};

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options,
) -> Result<ExecutionOutput> {
    let plan = config.plan();
    let mut execution = ShapleyExecution::new(&plan).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(worker, base_text, &config.render_parameters()?)?)
    } else { None };
    if let Some(path) = &options.resume {
        execution = checkpoint::restore_shapley(path, &plan, identity.expect("resume identity"))?;
    }
    // Reject changed law/model/format before reserving a destination. Every
    // accepted conditional half-pair is persisted, not only complete rows.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    if let Some(output) = &output {
        output.save_shapley(&execution, identity.expect("output identity"))?;
    }
    let initial = execution.observations().len();
    let allowance = options.max_new_samples.unwrap_or(config.samples).min(config.samples - initial);
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    for _ in 0..allowance {
        let report = execution.advance_interruptible(1, || Instant::now() >= deadline, |values| {
            let request = config.sample_request(base, values)?;
            let evaluated = match config.qoi {
                Qoi::Steady => evaluate_sample(worker, &request, deadline),
                qoi => evaluate_sample_for(worker, &request, deadline, qoi),
            };
            match evaluated {
                Ok(value) => Ok(Some(value)),
                Err(EvaluationError::Budget) => Ok(None),
                Err(EvaluationError::Child(message)) => Err(model_failure(message)),
            }
        });
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason.clone()
                .unwrap_or_else(|| "Shapley physical evaluation refused".into()));
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}")); error
                })?;
            }
            return Err(refusal);
        }
        if let Some(output) = &output {
            output.save_shapley(&execution, identity.expect("output identity"))?;
        }
        if matches!(report.status, UqStatus::Complete | UqStatus::Cancelled) { break; }
    }
    let report = execution.report();
    if report.status == UqStatus::Complete {
        return Ok(ExecutionOutput { stdout: render(base, config, &report)?, exit_code: exit::SUCCESS });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!(
            "Shapley sensitivity stopped after {} accepted physical evaluations; no shortened-design attribution published; use --checkpoint to retain paid conditional pairs",
            report.evaluations_accepted,
        )));
    };
    let width = 2 * ((1usize << plan.parameters.len()) - 1);
    let slot = report.evaluations_accepted % width;
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.shapley.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"base_samples\":{},\"evaluations_per_row\":{},\"completed_rows\":{},\"next_sample_ordinal\":{},\"next_row_ordinal\":{},\"next_subset_mask\":{},\"next_pair_half\":{},\"checkpoint\":{},\"no_claim\":\"retained conditional-pair prefix only; no sensitivity estimate or shortened-design inference; resume the same model, law, executable and fixed allocation\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?),
        report.evaluations_planned, report.evaluations_attempted,
        report.evaluations_attempted - initial, report.base_samples, width,
        report.completed_rows, report.evaluations_accepted, report.completed_rows,
        slot / 2, slot % 2, quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}

pub(super) fn render(base: &J, config: &Config, report: &ShapleyReport) -> Result<String> {
    if report.status != UqStatus::Complete {
        return Err(model_failure("only complete fixed-design Shapley reports may be published"));
    }
    let estimate = report.estimate.as_ref().ok_or_else(|| model_failure(format!(
        "physical Shapley effects unavailable after {} evaluations: {}",
        report.evaluations_accepted, report.unavailable_reason.unwrap_or("missing complete estimate"),
    )))?;
    let effects = estimate.effects.iter().map(|effect| Ok(format!(
        "{{\"parameter\":{},\"unit\":{},\"index\":{}}}",
        quote(&effect.parameter), quote(&effect.unit), number_json(effect.index)?,
    ))).collect::<Result<Vec<_>>>()?.join(",");
    let fractions = estimate.conditional_variance_fractions.iter().map(|&v| number_json(v))
        .collect::<Result<Vec<_>>>()?.join(",");
    let width = 2 * (estimate.conditional_variance_fractions.len() - 1);
    let mixed = config.plan().parameters.iter().any(|parameter|
        matches!(parameter.kind, fs_uq::UncertaintyKind::AleatoryUniform { .. }));
    let (method, law) = if mixed {
        ("conditional-product-philox", "explicitly independent Gaussian/uniform product law")
    } else {
        ("conditional-gaussian-philox", "declared Gaussian joint law")
    };
    let no_claim = format!(
        "{} Physical-input Shapley variance attribution under the {law}, not independent-germ Sobol effects, local derivatives or causal effects. Dependence and interactions both contribute. Noisy effects and conditional variance fractions remain unclipped; summing to one does not prove accuracy. Fixed-design sampling supplies no confidence bound, automatic convergence or model validation. A temperature ceiling does not turn this response into a compliance indicator.",
        config.qoi.no_claim(),
    );
    Ok(format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.shapley.v1\",\"authority\":\"estimated-direct-model-shapley\",\"method\":{},\"status\":\"complete\",\"qoi\":{},\"seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"base_samples\":{},\"evaluations_per_row\":{},\"base_output_std_dev_k\":{},\"correlation\":{},\"parameters\":[{}],\"effects\":[{}],\"conditional_variance_fractions\":[{}],\"no_claim\":{}}}\n",
        quote(method), config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
        report.evaluations_planned, report.evaluations_attempted, report.base_samples, width,
        optional_number(estimate.base_output_std_dev)?, quote(config.correlation_label),
        config.render_parameters()?, effects, fractions, quote(&no_claim),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }

    #[test]
    fn shapley_is_explicit_recoverable_and_preserves_the_sobol_selector() {
        let options = Options::parse(&args(&["--sensitivity", "shapley"])).unwrap();
        assert!(options.shapley_sensitivity && !options.sobol_sensitivity);
        let options = Options::parse(&args(&["--sensitivity", "sobol"])).unwrap();
        assert!(options.sobol_sensitivity && !options.shapley_sensitivity);
        assert!(!Options::parse(&[]).unwrap().shapley_sensitivity);
        assert!(Options::parse(&args(&["--sensitivity", "shapley", "--resume", "a.bin",
            "--checkpoint", "b.bin", "--max-new-samples", "1"])).is_ok());
    }

    #[test]
    fn incompatible_modes_refuse_in_either_argument_order() {
        for tail in [vec!["--sensitivity", "sobol"], vec!["--sensitivity", "shapley"],
            vec!["--qmc-replicates", "2"], vec!["--mean-control", "adjoint"],
            vec!["--pce", "1,6"], vec!["--pce-mean", "0,2,4"],
            vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "2"]]
        {
            for reversed in [false, true] {
                let selector = args(&["--sensitivity", "shapley"]);
                let extra = args(&tail);
                let combined = if reversed { [extra, selector].concat() } else { [selector, extra].concat() };
                assert!(Options::parse(&combined).is_err(), "accepted {combined:?}");
            }
        }
    }
}
