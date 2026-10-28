//! Polynomial uncertainty through the original cooling parser and physics.
//! A complete fixed fit exposes surrogate moments, not a physical mean bound.
use super::*;
use fs_uq::{PceBasis, PceConfig, PceExecution, PceReport, ProductPceModel};
mod mean;

#[derive(Debug, Clone, Copy)]
pub(super) struct Settings {
    degree: usize,
    training: usize,
    fit_samples: Option<usize>,
}
impl Settings {
    pub(super) fn parse(value: &OsString) -> Result<Self> {
        let text = value.to_str().ok_or_else(|| bad("--pce requires DEGREE,TRAINING"))?;
        let (degree, training) = text.split_once(',')
            .ok_or_else(|| bad("--pce requires DEGREE,TRAINING"))?;
        let degree = count_option(&OsString::from(degree), "PCE degree")?;
        let training = count_option(&OsString::from(training), "PCE training count")?;
        if degree > 8 || training < 2 {
            return Err(bad("PCE degree must be in 0..=8 and training count at least two"));
        }
        Ok(Self { degree, training, fit_samples: None })
    }
    pub(super) fn parse_mean(value: &OsString) -> Result<Self> {
        let text = value.to_str().ok_or_else(|| bad("--pce-mean requires DEGREE,TRAINING,FIT"))?;
        let (regression, fit) = text.rsplit_once(',')
            .ok_or_else(|| bad("--pce-mean requires DEGREE,TRAINING,FIT"))?;
        let mut settings = Self::parse(&OsString::from(regression))?;
        let fit = count_option(&OsString::from(fit), "PCE fit count")?;
        if settings.training.checked_add(2).is_none_or(|n| n > fit) {
            return Err(bad("PCE fit count must include training and at least two held-out calls"));
        }
        settings.fit_samples = Some(fit);
        Ok(settings)
    }
    fn regression(self) -> PceConfig { PceConfig::new(self.degree, self.training) }
}

pub(super) fn validate_options(options: &Options) -> Result<()> {
    if options.polynomial.is_some() && (options.qmc_replicates.is_some()
        || options.sobol_sensitivity || options.adjoint_mean_control
        || options.compliance.is_some() || options.design.is_some()
        || options.checkpoint.is_some() || options.resume.is_some()
        || options.max_new_samples.is_some())
    {
        return Err(bad("PCE requires a fixed design; QMC, sensitivity, adjoint control, sequential decisions, candidate selection and checkpoint/chunk options are unavailable"));
    }
    Ok(())
}

pub(super) fn execute(worker: &Path, base: &J, config: &Config, settings: Settings) -> Result<ExecutionOutput> {
    if let Some(fit_samples) = settings.fit_samples {
        return mean::execute(worker, base, config, settings, fit_samples);
    }
    // Include design setup and fitting in the same soft wall allowance as all
    // child solves. Bounded synchronous library reductions are not preemptible.
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let mut plan = config.plan();
    plan.method = PropagationMethod::PolynomialChaos;
    let mut execution = PceExecution::new(&plan, settings.regression()).map_err(bad)?;
    let report = execution.advance_interruptible(config.samples,
        || Instant::now() >= deadline,
        |values| evaluate(worker, config, base, values, deadline));
    complete(report.status, report.rejection_reason.as_deref(), report.samples_evaluated, deadline)?;
    let model = execution.model().ok_or_else(|| model_failure("completed PCE has no model"))?;
    let stdout = render(config, base, &report, model)?;
    // Do not turn a data-dependent timeout into a completed selected-prefix fit.
    complete(report.status, None, report.samples_evaluated, deadline)?;
    Ok(ExecutionOutput { stdout, exit_code: exit::SUCCESS })
}

fn evaluate(worker: &Path, config: &Config, base: &J, values: &[f64], deadline: Instant) -> Result<Option<f64>> {
    let request = config.sample_request(base, values)?;
    let result = match config.qoi {
        Qoi::Steady => evaluate_sample(worker, &request, deadline),
        qoi => evaluate_sample_for(worker, &request, deadline, qoi),
    };
    match result {
        Ok(value) => Ok(Some(value)),
        Err(EvaluationError::Budget) => Ok(None),
        Err(EvaluationError::Child(message)) => Err(model_failure(message)),
    }
}

fn complete(status: UqStatus, reason: Option<&str>, calls: usize, deadline: Instant) -> Result<()> {
    if status == UqStatus::Refused {
        return Err(model_failure(reason.unwrap_or("polynomial physical study refused")));
    }
    if status != UqStatus::Complete || Instant::now() >= deadline {
        return Err(budget(format!("PCE wall allowance exhausted after {calls} completed model calls; no partial estimate published; this mode does not support durable recovery")));
    }
    Ok(())
}

fn render(config: &Config, base: &J, report: &PceReport, model: &ProductPceModel) -> Result<String> {
    let summary = report.summary.as_ref().ok_or_else(|| model_failure("missing completed PCE summary"))?;
    let basis = model.basis().iter().map(|family| quote(match family {
        PceBasis::StandardNormal => "standard-normal-hermite",
        PceBasis::Uniform01 => "uniform-01-shifted-legendre",
    })).collect::<Vec<_>>().join(",");
    let indices = model.indices().iter().map(|index| {
        format!("[{}]", index.iter().map(usize::to_string).collect::<Vec<_>>().join(","))
    }).collect::<Vec<_>>().join(",");
    let coefficients = model.coefficients().iter().map(|&c| number_json(c))
        .collect::<Result<Vec<_>>>()?.join(",");
    let no_claim = format!("{} Polynomial moments describe the fitted surrogate, not the physical-model expectation or a confidence bound. Held-out residuals and physical compliance frequency are diagnostics, not certified error or approval. Uniform coordinates are independent U(0,1) germs; Gaussian coordinates are independent N(0,1) germs before the admitted correlation map. No automatic degree selection or durable recovery.", config.qoi.no_claim());
    Ok(format!(concat!(
        "{{\"schema\":\"frankensim.cooling-network-uq.pce.v1\",",
        "\"authority\":\"estimated-polynomial-surrogate\",\"method\":\"hermite-legendre-regression\",",
        "\"status\":\"complete\",\"qoi\":{},\"seed\":{},",
        "\"samples_planned\":{},\"samples_evaluated\":{},\"training_samples\":{},\"validation_samples\":{},",
        "\"degree\":{},\"basis_terms\":{},\"r_condition_estimate\":{},",
        "\"surrogate_mean_k\":{},\"surrogate_variance_k2\":{},\"surrogate_std_dev_k\":{},",
        "\"training_rmse_k\":{},\"validation_rmse_k\":{},\"validation_mean_error_k\":{},",
        "\"temperature_limit_k\":{},\"held_out_probability_of_compliance\":{},",
        "\"held_out_probability_standard_error\":{},\"germ_basis\":[{}],\"indices\":[{}],",
        "\"coefficients_k\":[{}],\"correlation\":{},\"parameters\":[{}],\"no_claim\":{}}}\n"),
        config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
        report.samples_planned, report.samples_evaluated, report.training_samples_accepted,
        report.validation_samples_accepted, report.config.degree, report.basis_terms,
        number_json(report.r_condition_estimate)?, number_json(summary.surrogate_mean)?,
        optional_number(summary.surrogate_variance)?, optional_number(summary.surrogate_std_dev)?,
        number_json(summary.training_rmse)?, number_json(summary.validation_rmse)?,
        number_json(summary.validation_mean_error)?, optional_number(config.threshold_k)?,
        optional_number(summary.probability_of_compliance)?, optional_number(summary.compliance_standard_error)?,
        basis, indices, coefficients, quote(config.correlation_label), config.render_parameters()?, quote(&no_claim),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }

    #[test]
    fn pce_layout_is_explicit_bounded_and_not_a_stopping_policy() {
        let settings = Options::parse(&args(&["--pce", "2,24"])).unwrap().polynomial.unwrap();
        assert_eq!(settings.degree, 2);
        assert_eq!(settings.training, 24);
        for value in ["", "2", "2,", ",24", "2,24,8", "9,24", "-1,24", "2,1", "2,10001", "2, 24"] {
            assert!(Options::parse(&args(&["--pce", value])).is_err(), "{value}");
        }
        for tail in [
            vec!["--pce", "1,8"], vec!["--qmc-replicates", "2"],
            vec!["--sensitivity", "sobol"], vec!["--mean-control", "adjoint"],
            vec!["--checkpoint", "new.bin"], vec!["--resume", "old.bin"],
            vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "8"],
        ] {
            let mut values = args(&["--pce", "1,8"]);
            values.extend(args(&tail));
            assert!(Options::parse(&values).is_err(), "{tail:?}");
        }
    }

    #[test]
    fn physical_mean_mode_cannot_hide_fitting_work_or_combine_with_a_surrogate_run() {
        let settings = Options::parse(&args(&["--pce-mean", "1,6,8"])).unwrap().polynomial.unwrap();
        assert_eq!(settings.fit_samples, Some(8));
        for value in ["1,6", "1,6,7", "1,6,8,12", "1,6,-8", "1,6,10001"] {
            assert!(Options::parse(&args(&["--pce-mean", value])).is_err(), "{value}");
        }
        for tail in [vec!["--pce", "1,6"], vec!["--pce-mean", "1,6,8"],
            vec!["--mean-control", "adjoint"], vec!["--resume", "old.bin"]] {
            let mut values = args(&["--pce-mean", "1,6,8"]);
            values.extend(args(&tail));
            assert!(Options::parse(&values).is_err(), "{tail:?}");
        }
        assert!(Options::parse(&args(&["--pce", "1,6", "--pce-mean", "1,6,8"])).is_err());
    }
}
