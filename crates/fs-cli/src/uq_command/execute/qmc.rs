//! Fixed-layout replicated QMC through the existing cooling parser/solver.
//! Recovery retains incomplete nets; estimates require the entire fixed layout.
//! This adapter never treats within-net points as independent MC observations.
use super::*;
use fs_uq::{QmcConfig, QmcExecution, QmcReport};

pub(super) fn execute(
    base_text: &str, base: &J, config: &Config, options: &Options, replicates: usize,
) -> Result<ExecutionOutput> {
    // No flooring an incomplete net or increasing the caller's sample budget.
    if replicates == 0 || config.samples % replicates != 0 {
        return Err(bad("samples must divide exactly into the declared QMC replicates"));
    }
    let layout = QmcConfig {
        replicates, samples_per_replicate: config.samples / replicates,
    };
    let mut plan = config.plan();
    plan.method = PropagationMethod::QuasiMonteCarlo;
    let mut execution = QmcExecution::new(&plan, layout).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(base_text, &config.render_parameters()?)?)
    } else { None };
    if let Some(path) = &options.resume {
        execution = checkpoint::restore_qmc(path, &plan, layout, identity.expect("resume identity"))?;
    }
    // Reject changed models/layouts and invalid checkpoints before reserving a
    // destination. Reuse the same fresh-path, bounded, atomic file transport.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    if let Some(output) = &output {
        output.save_qmc(&execution, identity.expect("output identity"))?;
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
            let refusal = Failure {
                code: "cooling-network-uq-refused",
                message: report.rejection_reason.clone().unwrap_or_else(|| "QMC model evaluation refused".into()),
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
            output.save_qmc(&execution, identity.expect("output identity"))?;
        }
        if matches!(report.status, UqStatus::Complete | UqStatus::Cancelled) { break; }
    }
    let report = execution.report();
    if report.status == UqStatus::Complete {
        return Ok(ExecutionOutput { stdout: render(config, base, &report)?, exit_code: exit::SUCCESS });
    }
    // A deadline can depend on the sampled physics. Even complete replicates
    // from a shortened run do not license publishing a fixed-layout estimate.
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!(
            "QMC wall-time budget exhausted after {} accepted samples and {} complete replicates; no partial distribution published; use --checkpoint to retain completed samples",
            report.samples_accepted, report.completed_replicates,
        )));
    };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.qmc.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"replicates\":{},\"samples_per_replicate\":{},\"completed_replicates\":{},\"next_sample_ordinal\":{},\"next_replicate_ordinal\":{},\"next_point_ordinal\":{},\"checkpoint\":{},\"no_claim\":\"retained QMC prefix, including incomplete nets; no completed distribution, probability estimate or optional-stopping inference; resume the identical base, plan, executable and fixed net layout\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?),
        config.samples, report.samples_evaluated, report.samples_evaluated - initial_count,
        layout.replicates, layout.samples_per_replicate, report.completed_replicates,
        report.samples_accepted, report.samples_accepted / layout.samples_per_replicate,
        report.samples_accepted % layout.samples_per_replicate, quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}

fn render(config: &Config, base: &J, report: &QmcReport) -> Result<String> {
    if report.status != UqStatus::Complete {
        return Err(model_failure("only complete fixed-layout QMC reports may be published"));
    }
    let estimate = report.estimate.as_ref().ok_or_else(|| model_failure("missing QMC estimate"))?;
    let replicate_means = report.replicate_means.iter().map(|&x| number_json(x))
        .collect::<Result<Vec<_>>>()?.join(",");
    let probability = report.compliance.as_ref().map(|x| x.mean);
    let probability_error = report.compliance.as_ref().and_then(|x| x.standard_error);
    let no_claim = format!(
        "{} Replicate standard errors describe variation between independently keyed Sobol scrambles, not independent net points. They are not confidence intervals, sequential stopping guarantees or physical compliance approval. The 32-bit midpoint grid and approximate inverse-normal transform have unbounded discretization bias here; zero replicate error does not prove exactness.",
        config.qoi.no_claim(),
    );
    Ok(format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.qmc.v1\",\"authority\":\"estimated-randomized-quadrature\",\"method\":\"owen-sobol-replicated\",\"status\":\"complete\",\"qoi\":{},\"seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"replicates\":{},\"samples_per_replicate\":{},\"completed_replicates\":{},\"mean_k\":{},\"between_replicate_standard_error_k\":{},\"replicate_means_k\":[{}],\"temperature_limit_k\":{},\"estimated_probability_of_compliance\":{},\"between_replicate_probability_standard_error\":{},\"correlation\":{},\"parameters\":[{}],\"no_claim\":{}}}\n",
        config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
        report.samples_planned, report.samples_evaluated, report.config.replicates,
        report.config.samples_per_replicate, report.completed_replicates,
        number_json(estimate.mean)?, optional_number(estimate.standard_error)?, replicate_means,
        optional_number(config.threshold_k)?, optional_number(probability)?, optional_number(probability_error)?,
        quote(config.correlation_label), config.render_parameters()?, quote(&no_claim),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }

    #[test]
    fn qmc_recovery_is_explicit_and_never_reuses_iid_stopping() {
        assert_eq!(Options::parse(&args(&["--qmc-replicates", "3"])).unwrap().qmc_replicates, Some(3));
        assert!(Options::parse(&[]).unwrap().qmc_replicates.is_none());
        assert!(Options::parse(&args(&["--qmc-replicates", "2", "--checkpoint", "a.bin", "--max-new-samples", "3"])).is_ok());
        assert!(Options::parse(&args(&["--qmc-replicates", "2", "--resume", "a.bin", "--checkpoint", "b.bin"])).is_ok());
        for tail in [
            vec!["--max-new-samples", "0"],
            vec!["--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "16"],
            vec!["--qmc-replicates", "4"],
        ] {
            let mut values = args(&["--qmc-replicates", "3"]); values.extend(args(&tail));
            assert!(Options::parse(&values).is_err(), "accepted {tail:?}");
        }
        for count in ["0", "1", "257", "-1", "NaN"] {
            assert!(Options::parse(&args(&["--qmc-replicates", count])).is_err());
        }
    }
}
