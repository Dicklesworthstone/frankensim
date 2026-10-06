//! Project the shared Jansen estimator; do not reinterpret hybrids as iid data.
use fs_project::uncertainty::UniformParameter;
use fs_uq::{SobolReport, UqStatus};
use super::{Result, fail, optional, quoted};

pub(super) const SCOPE: &str = "Global Sobol main and total effects estimated from a predeclared independent-uniform pick-freeze design through actual native physical solves. Each row is A, B, then A with coordinate i replaced by B_i in parameter declaration order. Every base and hybrid solve consumes the original total evaluation budget. Hybrid observations are dependent: no Monte Carlo standard error, pass probability, quantiles, mean control or sequential compliance decision is inferred from this design. Only the complete fixed design produces indices; zero sampled base variance leaves them undefined, not zero. Finite-sample indices are not clipped to [0,1], and total effects need not sum to one. Effects describe this input law and numerical model, not local derivatives, causal effects, certified rankings, confidence intervals, continuum or physical validation.";

pub(super) fn json(report: &SobolReport, parameters: &[UniformParameter]) -> Result<String> {
    let width = parameters.len() + 2;
    let order = parameters.iter().map(|p| quoted(&p.name)).collect::<Vec<_>>().join(",");
    let estimate = if let Some(estimate) = &report.estimate {
        if report.status != UqStatus::Complete || estimate.effects.len() != parameters.len() {
            return Err(fail("cli-uncertainty-sensitivity", "incomplete or mismatched sensitivity estimate"));
        }
        let mut rows = Vec::with_capacity(parameters.len());
        for (effect, parameter) in estimate.effects.iter().zip(parameters) {
            if effect.parameter != parameter.name || effect.unit != parameter.target.unit()
                || !effect.first_order.is_finite() || !effect.total_order.is_finite() {
                return Err(fail("cli-uncertainty-sensitivity", "sensitivity effect differs from the declared physical input"));
            }
            rows.push(format!(
                "{{\"parameter\":{},\"entity\":{},\"parameter_unit\":{},\"index_unit\":\"1\",\"first_order\":{},\"total_order\":{}}}",
                quoted(&effect.parameter), quoted(&parameter.entity), quoted(&effect.unit),
                effect.first_order, effect.total_order));
        }
        format!("{{\"base_output_std_dev_k\":{},\"effects\":[{}]}}",
            optional(estimate.base_output_std_dev), rows.join(","))
    } else { "null".into() };
    Ok(format!(concat!("{{\"estimator\":\"jansen-pick-freeze\",\"sampler\":\"philox\",",
        "\"parameter_order\":[{}],\"row_order\":\"A,B,A_with_B_i_in_parameter_order\",",
        "\"base_rows_planned\":{},\"completed_rows\":{},\"evaluations_per_row\":{},",
        "\"evaluations_planned\":{},\"evaluations_accepted\":{},\"partial_row_evaluations\":{},",
        "\"estimate\":{},\"unavailable_reason\":{}}}"),
        order, report.base_samples, report.completed_rows, width,
        report.evaluations_planned, report.evaluations_accepted,
        report.evaluations_accepted % width, estimate,
        report.unavailable_reason.map_or_else(|| "null".into(), quoted)))
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
        .replace('"', "&quot;").replace('\'', "&#39;")
}

pub(super) fn html(report: &SobolReport) -> String {
    let Some(estimate) = &report.estimate else {
        return format!("<h2>Global parameter sensitivity</h2><p>{}/{} complete pick-freeze rows. Indices unavailable: {}.</p>",
            report.completed_rows, report.base_samples,
            escape(report.unavailable_reason.unwrap_or("the full fixed design has not completed successfully")));
    };
    let rows = estimate.effects.iter().map(|e| format!(
        "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
        escape(&e.parameter), escape(&e.unit), e.first_order, e.total_order))
        .collect::<String>();
    format!("<h2>Global parameter sensitivity</h2><p>{} complete pick-freeze rows, {} native evaluations. Estimated dimensionless variance shares; no confidence bounds or clipping.</p><table><tr><th>Input</th><th>Input unit</th><th>Main effect</th><th>Total effect</th></tr>{rows}</table>",
        report.completed_rows, report.evaluations_accepted)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_project::uncertainty::Target;
    use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, SobolExecution, UqPlan};
    use crate::json_read::JsonValue as J;

    fn parameters() -> Vec<UniformParameter> {
        vec![UniformParameter { name: "power<&>".into(), entity: "solid".into(),
            target: Target::Power, low: 1.0, high: 2.0 }]
    }
    fn plan() -> UqPlan {
        UqPlan::new("temperature-max", PropagationMethod::MonteCarlo, 6)
            .with_correlation(CorrelationModel::Independent)
            .with_parameter(ParameterUncertainty::uniform("power<&>", 1.0, 2.0, "W"))
    }
    #[test]
    fn partial_and_refused_designs_never_publish_indices_or_iid_statistics() {
        let mut run = SobolExecution::new(&plan()).unwrap();
        run.advance(2, || false, |x| Ok::<_, &str>(300.0 + x[0]));
        let report = run.report();
        let value = J::parse(&json(&report, &parameters()).unwrap()).unwrap();
        assert_eq!(value.f64_field("partial_row_evaluations"), Some(2.0));
        assert_eq!(value.get("estimate"), Some(&J::Null));
        run.advance(1, || false, |_| Err::<f64, _>("physical refusal"));
        assert_eq!(run.report().status, UqStatus::Refused);
        assert!(json(&run.report(), &parameters()).unwrap().contains("\"estimate\":null"));
        assert!(html(&run.report()).contains("unavailable"));
    }
    #[test]
    fn undefined_normalization_and_unclipped_effects_remain_explicit() {
        let mut run = SobolExecution::new(&plan()).unwrap();
        run.advance(6, || false, |_| Ok::<_, &str>(300.0));
        let report = run.report();
        assert_eq!(report.status, UqStatus::Complete);
        let value = J::parse(&json(&report, &parameters()).unwrap()).unwrap();
        assert_eq!(value.get("estimate"), Some(&J::Null));
        assert!(value.str_field("unavailable_reason").unwrap().contains("variance"));
        let mut run = SobolExecution::new(&plan()).unwrap();
        let mut report = run.advance(6, || false, |x| Ok::<_, &str>(300.0 + x[0]));
        // Projection must not sanitize legitimate finite-sample out-of-range indices.
        let e = &mut report.estimate.as_mut().unwrap().effects[0];
        e.first_order = -0.25; e.total_order = 1.5;
        let text = json(&report, &parameters()).unwrap();
        assert!(text.contains("\"first_order\":-0.25") && text.contains("\"total_order\":1.5"));
        assert!(html(&report).contains("power&lt;&amp;&gt;"));
        let mut wrong = parameters(); wrong[0].name = "another control".into();
        assert!(json(&report, &wrong).is_err());
    }
}
