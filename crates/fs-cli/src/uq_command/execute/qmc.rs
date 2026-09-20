//! Fixed-layout replicated QMC through the existing cooling parser/solver.
//! This adapter never calls MC confidence-sequence or checkpoint code.
use super::*;
use fs_uq::{QmcConfig, QmcExecution, QmcReport};

pub(super) fn execute(base: &J, config: &Config, replicates: usize) -> Result<ExecutionOutput> {
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
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let report = execution.advance_interruptible(config.samples, || Instant::now() >= deadline, |values| {
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
        return Err(Failure {
            code: "cooling-network-uq-refused",
            message: report.rejection_reason.clone().unwrap_or_else(|| "QMC model evaluation refused".into()),
        });
    }
    if report.status != UqStatus::Complete {
        // A deadline is potentially value-dependent: do not publish a shortened
        // estimate as a successful fixed-budget result or infer convergence.
        return Err(budget(format!(
            "QMC wall-time budget exhausted after {} accepted samples and {} complete replicates; no partial distribution published; durable QMC resume is not implemented",
            report.samples_accepted, report.completed_replicates,
        )));
    }
    Ok(ExecutionOutput { stdout: render(config, base, &report)?, exit_code: exit::SUCCESS })
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
    fn qmc_is_explicit_and_never_reuses_iid_stopping_or_mc_recovery() {
        assert_eq!(Options::parse(&args(&["--qmc-replicates", "3"])).unwrap().qmc_replicates, Some(3));
        assert!(Options::parse(&[]).unwrap().qmc_replicates.is_none());
        for tail in [
            vec!["--checkpoint", "unused.bin"],
            vec!["--resume", "unused.bin"],
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
