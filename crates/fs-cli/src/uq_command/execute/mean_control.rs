//! One real nominal coupled adjoint, fixed before the existing MC loop.
//! Raw samples and their compliance/quantile interpretation never change.
use super::*;
use fs_uq::{LinearControlVariate, UqControlError};

#[derive(Debug)]
enum Target {
    Inlet(usize),
    FanSpeed,
    SurfaceHtc(String),
    ContactResistance(String),
}

pub(super) struct MeanControl {
    frozen: LinearControlVariate,
    nominal_temperature: f64,
    adjoint_residual: f64,
}

pub(super) fn prepare(base: &J, config: &Config, execution: &UqExecution,
    deadline: Instant) -> Result<MeanControl>
{
    prepare_with(base, config, execution, |request| {
        child::evaluate_document(request, deadline).map_err(|error| match error {
            EvaluationError::Budget => budget("shared UQ wall allowance expired during the nominal forward/adjoint"),
            EvaluationError::Child(message) => model_failure(message),
        })
    })
}

fn prepare_with(base: &J, config: &Config, execution: &UqExecution,
    evaluate: impl FnOnce(&str) -> Result<J>) -> Result<MeanControl>
{
    if config.qoi != Qoi::Steady
        || ["transient", "radiation", "recirculation", "mesh_convergence"].iter()
            .any(|key| base.get(key).is_some()) {
        return Err(bad("adjoint mean control requires steady cooling without radiation, recirculation or mesh studies; unsupported derivative paths are not assigned zero"));
    }
    if config.samples >= MAX_PRODUCT_SAMPLES {
        return Err(bad("nominal forward/adjoint plus Monte Carlo samples must fit the total 10000-model-call cap"));
    }
    if execution.evaluations_attempted() != 0 || !execution.observations().is_empty() {
        return Err(bad("freeze the nominal control before sampling"));
    }
    if execution.plan() != &config.plan() {
        return Err(bad("nominal control and sample execution must share the same complete plan"));
    }
    // Reuse the adapter's canonical target declarations rather than maintain a
    // second uncertainty grammar or guess the physical meaning from plan names.
    let encoded = J::parse(&format!("[{}]", config.render_parameters()?))
        .map_err(|error| bad(error.to_string()))?;
    let targets = array(&encoded, "control parameters", 256)?.iter()
        .map(|row| Target::parse(field(row, "target")?)).collect::<Result<Vec<_>>>()?;
    let means = execution.parameter_means();
    if targets.len() != means.len() { return Err(bad("control parameter arity mismatch")); }
    let mut nominal_base = base.clone();
    enable_gradient(&mut nominal_base)?;
    // This sets every physical parameter to its actual declared expectation;
    // it does not assume that the base request is the distribution mean.
    let request = config.sample_request(&nominal_base, &means)?;
    let document = evaluate(&request)?;
    let nominal_temperature = config.qoi.extract(&document)?;
    let adjoint_residual = finite_field(&document, "adjoint_residual")?;
    if adjoint_residual < 0.0 { return Err(model_failure("negative nominal adjoint residual")); }
    let gradient = targets.iter().zip(&means).map(|(target, &mean)|
        target.derivative(&document, mean)).collect::<Result<Vec<_>>>()?;
    let frozen = execution.freeze_linear_control_variate(&gradient).map_err(control_error)?;
    Ok(MeanControl { frozen, nominal_temperature, adjoint_residual })
}

impl Target {
    fn parse(value: &J) -> Result<Self> {
        match value.str_field("kind") {
            Some("inlet-temperature") => Ok(Self::Inlet(integer_raw(field(value,"index")?,"inlet index")?)),
            Some("fan-speed-ratio") => Ok(Self::FanSpeed),
            Some("surface-htc") => Ok(Self::SurfaceHtc(string(field(value,"surface")?,"surface")?)),
            Some("contact-resistance") => Ok(Self::ContactResistance(string(field(value,"contact")?,"contact")?)),
            _ => Err(bad("adjoint mean control supports inlet-temperature, fan-speed-ratio, declared surface-htc and contact-resistance only; other UQ inputs retain the ordinary sampling path")),
        }
    }

    fn derivative(&self, document: &J, mean: f64) -> Result<f64> {
        let derivative = match self {
            Self::Inlet(index) => {
                let values = field(document,"dobjective_dinlet_k")?.as_array()
                    .ok_or_else(|| model_failure("nominal producer supplied no inlet adjoint"))?;
                number(values.get(*index).ok_or_else(||
                    model_failure("nominal inlet adjoint does not cover the declared inlet index"))?, "inlet derivative")?
            }
            Self::FanSpeed => {
                let row = field(document, "fan_speed_sensitivity")?;
                if row.str_field("status") != Some("available")
                    || row.str_field("method") != Some("fan-affinity-coupled-adjoint") {
                    return Err(model_failure("nominal producer has no admitted total fan-speed derivative"));
                }
                // The UQ variable is LINEAR speed, not log(speed). The total
                // adjoint includes the fan, capacity and convection chain rule.
                linear_derivative(finite_field(row,"dobjective_dlog_speed_ratio_k")?, mean)?
            }
            Self::SurfaceHtc(name) => {
                let row = named_row(field(document,"walls")?, "region", name)?;
                if finite_field(row,"htc_w_m2_k")? != mean {
                    return Err(model_failure("nominal surface coefficient differs from its declared parameter mean"));
                }
                linear_derivative(finite_field(row,"dobjective_dlog_htc")?, mean)?
            }
            Self::ContactResistance(name) => {
                let contacts = field(document,"contact_sensitivities")?;
                if contacts.str_field("method") != Some("coupled-adjoint-contact-bilinear-form") {
                    return Err(model_failure("nominal producer supplied no steady coupled contact adjoint"));
                }
                let row = named_row(field(contacts,"rows")?, "contact", name)?;
                if finite_field(row,"resistance_m2_k_w")? != mean {
                    return Err(model_failure("nominal contact resistance differs from its declared parameter mean"));
                }
                linear_derivative(finite_field(row,"dobjective_dlog_resistance_k")?, mean)?
            }
        };
        if derivative.is_finite() { Ok(derivative) }
        else { Err(model_failure("nominal linear-coordinate derivative is nonfinite")) }
    }
}

fn linear_derivative(log_derivative: f64, parameter: f64) -> Result<f64> {
    if !(parameter.is_finite() && parameter > 0.0) {
        return Err(model_failure("log-coordinate adjoint requires a positive finite parameter mean"));
    }
    let value = log_derivative / parameter;
    if value.is_finite() { Ok(value) }
    else { Err(model_failure("nominal linear-coordinate derivative overflowed")) }
}

fn finite_field(value: &J, name: &str) -> Result<f64> {
    value.get(name).and_then(J::as_f64).filter(|v| v.is_finite())
        .ok_or_else(|| model_failure(format!("nominal producer has no finite {name}")))
}

fn named_row<'a>(value: &'a J, key: &str, name: &str) -> Result<&'a J> {
    let rows = value.as_array().ok_or_else(|| model_failure("nominal derivative rows are absent"))?;
    let mut matching = rows.iter().filter(|row| row.str_field(key) == Some(name));
    let row = matching.next().ok_or_else(|| model_failure(format!("nominal derivative missing for {name}")))?;
    if matching.next().is_some() { return Err(model_failure(format!("ambiguous nominal derivative for {name}"))); }
    Ok(row)
}

fn enable_gradient(base: &mut J) -> Result<()> {
    let J::Object(rows) = base else { return Err(bad("nominal request must be an object")); };
    let Some((_, J::Object(objective))) = rows.iter_mut().find(|(key, _)| key == "objective")
        else { return Err(bad("nominal request has no objective")); };
    let gradient = objective.iter_mut().find(|(key, _)| key == "gradient")
        .ok_or_else(|| bad("nominal objective must explicitly declare gradient=false"))?;
    if gradient.1 != J::Bool(false) { return Err(bad("sampling base must keep gradient=false")); }
    gradient.1 = J::Bool(true);
    Ok(())
}

impl MeanControl {
    pub(super) fn attach(&self, raw_output: String, execution: &UqExecution,
        deadline: Instant) -> Result<String>
    {
        if execution.report().status != UqStatus::Complete {
            return Err(bad("a controlled mean requires the complete fixed sample count"));
        }
        let estimate = execution.assess_linear_control_variate_interruptible(&self.frozen,
            || Instant::now() >= deadline).map_err(control_error)?
            .ok_or_else(|| model_failure("complete controlled run has no observations"))?;
        let values = |values: &[f64]| -> Result<String> {
            Ok(format!("[{}]", values.iter().copied().map(number_json)
                .collect::<Result<Vec<_>>>()?.join(",")))
        };
        let controlled = format!(
            "{{\"method\":\"frozen-nominal-coupled-adjoint\",\"authority\":\"estimated-fixed-sample-mean\",\"nominal_forward_adjoint_evaluations\":1,\"sample_model_evaluations\":{},\"total_model_evaluations\":{},\"nominal_objective_k\":{},\"nominal_adjoint_residual\":{},\"parameter_means\":{},\"gradient_k_per_parameter_unit\":{},\"raw_mean_k\":{},\"mean_k\":{},\"adjusted_std_dev_k\":{},\"sampling_standard_error_k\":{},\"adjusted_to_raw_variance_ratio\":{},\"scope\":\"mean of the declared numerical model only; coefficients fixed before sampling, analytic marginal centering; raw quantiles, bounds and compliance indicators unchanged; standard errors are descriptive, not optional-stopping bounds; no guaranteed variance reduction or native .fsim integration\"}}",
            estimate.n, estimate.n + 1, number_json(self.nominal_temperature)?,
            number_json(self.adjoint_residual)?, values(self.frozen.parameter_means())?,
            values(self.frozen.gradient())?, number_json(estimate.raw_mean)?, number_json(estimate.mean)?,
            optional_number(estimate.std_dev)?, optional_number(estimate.standard_error)?,
            optional_number(estimate.variance_ratio)?,
        );
        let prefix = raw_output.strip_suffix("}\n")
            .ok_or_else(|| model_failure("unexpected complete UQ result framing"))?;
        if Instant::now() >= deadline { return Err(budget("shared UQ wall allowance expired before controlled mean publication")); }
        Ok(format!("{prefix},\"mean_control_variate\":{controlled}}}\n"))
    }
}

fn control_error(error: UqControlError) -> Failure {
    if error == UqControlError::Cancelled { budget(error.to_string()) }
    else { model_failure(error.to_string()) }
}

#[cfg(test)]
mod tests;
