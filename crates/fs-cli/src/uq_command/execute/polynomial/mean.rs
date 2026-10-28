//! Independent physical correction of a frozen cooling surrogate, one budget.
use super::*;
use fs_uq::{PceMeanConfig, PceMeanExecution, PceMeanReport};

pub(super) fn execute(worker: &Path, base: &J, config: &Config, settings: Settings, fit_samples: usize)
    -> Result<ExecutionOutput>
{
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let policy = PceMeanConfig {
        fit_evaluations: fit_samples,
        regression: settings.regression(),
        // This mode is explicitly a unit-weight frozen polynomial correction.
        // Neither the fit's validation data nor production tunes this weight.
        control_weight: 1.0,
    };
    let mut execution = PceMeanExecution::new(&config.plan(), policy).map_err(bad)?;
    let report = execution.advance_interruptible(config.samples,
        || Instant::now() >= deadline,
        |values| evaluate(worker, config, base, values, deadline));
    complete(report.status, report.rejection_reason.as_deref(), report.model_calls_evaluated, deadline)?;
    let stdout = render_mean(config, base, &execution, &report)?;
    complete(report.status, None, report.model_calls_evaluated, deadline)?;
    Ok(ExecutionOutput { stdout, exit_code: exit::SUCCESS })
}

fn render_mean(config: &Config, base: &J, execution: &PceMeanExecution, report: &PceMeanReport)
    -> Result<String>
{
    let estimate = report.estimate.as_ref().ok_or_else(|| model_failure("missing completed physical-mean estimate"))?;
    let raw = report.production.as_ref().ok_or_else(|| model_failure("missing independent production results"))?;
    if raw.status != UqStatus::Complete {
        return Err(model_failure("physical-mean publication needs complete production"));
    }
    let model = execution.fit().model().ok_or_else(|| model_failure("missing frozen cooling polynomial"))?;
    let fit = super::render(config, base, &report.fit, model)?;
    let percentiles = match raw.percentiles {
        Some(p) => format!("[{},{},{}]", number_json(p[0])?, number_json(p[1])?, number_json(p[2])?),
        None => "null".into(),
    };
    let no_claim = format!("{} The controlled mean uses only independent production cooling responses minus the centered frozen polynomial, at fixed weight one. Fit and held-out responses do not enter that mean. All physical ranges, quantiles and compliance frequencies use unadjusted production responses. Sampling standard errors and observed variance ratios are descriptive, not confidence bounds or stopping rules; a control can increase variance. No finite-generator, surrogate-rounding, solver or physical-model error bound, durable recovery or guaranteed speedup.", config.qoi.no_claim());
    Ok(format!(concat!(
        "{{\"schema\":\"frankensim.cooling-network-uq.pce-mean.v1\",",
        "\"authority\":\"estimated-controlled-physical-mean\",\"status\":\"complete\",",
        "\"method\":\"independent-pce-mean-control\",\"qoi\":{},\"control_weight\":1,",
        "\"model_calls_planned\":{},\"model_calls_evaluated\":{},\"fit_evaluations\":{},",
        "\"production_samples_planned\":{},\"production_samples_evaluated\":{},",
        "\"fit_seed\":{},\"production_seed\":{},\"mean_k\":{},\"raw_mean_k\":{},",
        "\"sampling_standard_error_k\":{},\"raw_sampling_standard_error_k\":{},",
        "\"observed_variance_ratio\":{},\"raw_production_std_dev_k\":{},",
        "\"raw_production_percentiles_p05_p50_p95_k\":{},\"raw_production_empirical_bounds_k\":[{},{}],",
        "\"temperature_limit_k\":{},\"raw_production_probability_of_compliance\":{},",
        "\"fit\":{},\"no_claim\":{}}}\n"),
        config.qoi.render(objective_kind(base)?), report.model_calls_planned,
        report.model_calls_evaluated, report.fit.samples_evaluated,
        execution.production().plan().budget_max_samples, estimate.samples,
        quote(&report.fit_seed.to_string()), quote(&report.production_seed.to_string()),
        number_json(estimate.mean)?, number_json(estimate.raw_mean)?,
        optional_number(estimate.standard_error)?, optional_number(estimate.raw_standard_error)?,
        optional_number(estimate.variance_ratio)?, optional_number(raw.std_dev)?, percentiles,
        number_json(raw.interval_bounds[0])?, number_json(raw.interval_bounds[1])?,
        optional_number(config.threshold_k)?, optional_number(raw.probability_of_compliance)?,
        fit.trim_end(), quote(&no_claim),
    ))
}
