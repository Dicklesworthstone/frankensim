use super::*;
use child::{EvaluationError, evaluate_sample};
use model::Config;
use std::time::{Duration, Instant};

pub(super) fn execute(base_text: &str, uq_text: &str) -> Result<String> {
    let base = J::parse(base_text).map_err(|error| bad(format!("invalid base JSON: {error}")))?;
    let config = Config::parse(uq_text, &base)?;
    let mut execution = UqExecution::new(&config.plan()).map_err(bad)?;
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let mut timed_out = false;
    let report = execution.advance(config.samples, || Instant::now() >= deadline, |values| {
        let request = config.sample_request(&base, values)?;
        match evaluate_sample(&request, deadline) {
            Ok(value) => Ok(value),
            Err(EvaluationError::Budget) => {
                timed_out = true;
                Err(budget("overall UQ wall-time budget expired during a cooling sample"))
            }
            Err(EvaluationError::Child(message)) => Err(model_failure(message)),
        }
    });
    if timed_out || report.status == UqStatus::Cancelled {
        return Err(budget(format!("UQ wall-time budget exhausted after {} model evaluations; no partial distribution published", report.samples_evaluated)));
    }
    if report.status != UqStatus::Complete {
        return Err(Failure { code: "cooling-network-uq-refused", message: report.rejection_reason.clone().unwrap_or_else(|| format!("UQ ended with {}", report.status.label())) });
    }
    render_result(&config, &base, &report)
}

fn render_result(config: &Config, base: &J, result: &fs_uq::UqResult) -> Result<String> {
    let percentiles = match result.percentiles {
        Some(values) => format!("[{},{},{}]", number_json(values[0])?, number_json(values[1])?, number_json(values[2])?),
        None => "null".into(),
    };
    let bounds = if result.mean.is_some() {
        format!("[{},{}]", number_json(result.interval_bounds[0])?, number_json(result.interval_bounds[1])?)
    } else { "null".into() };
    let objective = field(base, "objective")?;
    let objective_kind = if objective.get("mean_wall_region").is_some() {
        "mean_wall_region"
    } else if objective.get("max_wall_region").is_some() {
        "max_wall_region"
    } else if objective.get("max_solid_temperature").is_some() {
        "max_solid_temperature"
    } else { "max_vertices" };
    Ok(format!(
        "{{\"schema\":{},\"authority\":\"estimated-empirical-monte-carlo\",\"status\":{},\"qoi\":{{\"kind\":{},\"unit\":\"K\"}},\"seed\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"mean_k\":{},\"std_dev_k\":{},\"percentiles_p05_p50_p95_k\":{},\"empirical_bounds_k\":{},\"temperature_limit_k\":{},\"empirical_probability_of_compliance\":{},\"sampling_standard_error_k\":{},\"correlation\":{},\"parameters\":[{}],\"no_claim\":\"fixed-count empirical propagation through actual cooling-network child solves; no confidence sequence, optional-stopping guarantee, physical/model-form uncertainty bound, mesh-convergence certificate, experimental validation, or native .fsim/ledger package claim\"}}\n",
        quote(RESULT_SCHEMA), quote(result.status.label()), quote(objective_kind), quote(&config.seed.to_string()),
        config.samples, result.samples_evaluated, optional_number(result.mean)?, optional_number(result.std_dev)?,
        percentiles, bounds, optional_number(config.threshold_k)?, optional_number(result.probability_of_compliance)?,
        number_json(result.sampling_error)?, quote(config.correlation_label), config.render_parameters()?
    ))
}
