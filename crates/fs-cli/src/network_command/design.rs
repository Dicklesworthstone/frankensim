//! A bounded, adjoint-guided target search over the SAME file-driven solid.
//! Only one effective h varies; hydraulics and geometry remain immutable.
//! A returned crossing is local: no global monotonicity/minimality or
//! infeasibility theorem is inferred for an arbitrary thermal network.

use super::*;
use fs_math::det;

#[derive(Debug)]
pub(super) struct DesignRequest {
    surface: String,
    limit: f64,
    minimum: f64,
    maximum: f64,
    temperature_tolerance: f64,
    log_tolerance: f64,
    evaluations: usize,
}

impl DesignRequest {
    pub(super) fn parse(value: &J, surfaces: &BTreeSet<String>) -> Result<Self> {
        object(value, &["surface", "mean_temperature_limit_k", "min_htc_w_m2_k", "max_htc_w_m2_k",
            "temperature_tolerance_k", "log_htc_tolerance", "max_evaluations"], "design")?;
        let result = Self {
            surface: string(get(value, "surface")?, "design.surface")?,
            limit: positive(get(value, "mean_temperature_limit_k")?, "mean_temperature_limit_k")?,
            minimum: positive(get(value, "min_htc_w_m2_k")?, "min_htc_w_m2_k")?,
            maximum: positive(get(value, "max_htc_w_m2_k")?, "max_htc_w_m2_k")?,
            temperature_tolerance: positive(get(value, "temperature_tolerance_k")?, "design.temperature_tolerance_k")?,
            log_tolerance: positive(get(value, "log_htc_tolerance")?, "log_htc_tolerance")?,
            evaluations: count(get(value, "max_evaluations")?, "max_evaluations", 4096)?,
        };
        if !surfaces.contains(&result.surface) || result.minimum >= result.maximum {
            return Err(bad("design requires an existing surface and strictly ordered positive h bounds"));
        }
        Ok(result)
    }
}

pub(super) struct Designed {
    pub passing: Evaluation,
    surface: String,
    limit: f64,
    passing_h: f64,
    failed_lower: Option<(f64, f64)>,
    log_width: f64,
    history: Vec<(f64, f64)>,
    solid_solves: usize,
}

fn evaluate(request: &Request, cx: &Cx<'_>, flow: &GraphSolution, design: &DesignRequest,
    h: f64, history: &mut Vec<(f64, f64)>, solid_solves: &mut usize) -> Result<Evaluation> {
    poll(cx)?;
    if history.len() >= design.evaluations {
        return Err(Failure { code: "cooling-network-design-budget",
            message: format!("{} coupled design evaluations exhausted; no design published", history.len()) });
    }
    let mut coefficients: BTreeMap<_, _> = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
    coefficients.insert(design.surface.clone(), h);
    let result = request.evaluate(cx, flow, &coefficients, true)?;
    if !result.objective.is_finite() { return Err(producer("nonfinite design objective")); }
    *solid_solves = solid_solves.checked_add(result.coupled.iterations).ok_or_else(|| bad("design work count overflow"))?;
    history.push((h, result.objective));
    poll(cx)?;
    Ok(result)
}

pub(super) fn solve(request: &Request, cx: &Cx<'_>, flow: &GraphSolution, design: &DesignRequest) -> Result<Designed> {
    let mut history = Vec::new();
    let mut solid_solves = 0;
    let lower = evaluate(request, cx, flow, design, design.minimum, &mut history, &mut solid_solves)?;
    if lower.objective <= design.limit {
        return Ok(Designed { passing: lower, surface: design.surface.clone(), limit: design.limit,
            passing_h: design.minimum, failed_lower: None, log_width: 0.0, history, solid_solves });
    }
    let mut passing = evaluate(request, cx, flow, design, design.maximum, &mut history, &mut solid_solves)?;
    if passing.objective > design.limit {
        return Err(Failure { code: "cooling-network-design-bracket", message: format!(
            "neither endpoint meets the {} K mean-wall limit: lower {} K, upper {} K; this is not a proof of infeasibility between the bounds",
            design.limit, lower.objective, passing.objective) });
    }
    let mut low_h = design.minimum;
    let mut low_value = lower.objective;
    let mut high_h = design.maximum;
    let mut low = det::ln(low_h);
    let mut high = det::ln(high_h);
    loop {
        poll(cx)?;
        let width = high - low;
        let slack = design.limit - passing.objective;
        if width <= design.log_tolerance && slack <= design.temperature_tolerance {
            return Ok(Designed { passing, surface: design.surface.clone(), limit: design.limit,
                passing_h: high_h, failed_lower: Some((low_h, low_value)), log_width: width, history, solid_solves });
        }
        let index = passing.coupled.solid.iter().position(|s| s.region == design.surface)
            .ok_or_else(|| bad("design surface missing from evaluated output"))?;
        let derivative = passing.gradient.as_ref().ok_or_else(|| bad("design adjoint was not computed"))?.log_htc[index];
        // Both endpoints are evaluated. A Newton proposal must lie strictly in
        // the central 80% of that log-space bracket; otherwise bisect. The
        // derivative guides work, not the meaning of a passing endpoint.
        let newton = high - (passing.objective - design.limit) / derivative;
        let candidate = if derivative != 0.0 && newton.is_finite()
            && newton > low + 0.1 * width && newton < high - 0.1 * width { newton }
            else { 0.5 * low + 0.5 * high };
        let h = det::exp(candidate);
        if !(candidate > low && candidate < high && h.is_finite() && h > low_h && h < high_h) {
            return Err(Failure { code: "cooling-network-design-resolution",
                message: "floating-point h resolution cannot meet both requested design tolerances".into() });
        }
        let evaluated = evaluate(request, cx, flow, design, h, &mut history, &mut solid_solves)?;
        if evaluated.objective <= design.limit {
            high = candidate;
            high_h = h;
            passing = evaluated;
        } else {
            low = candidate;
            low_h = h;
            low_value = evaluated.objective;
        }
    }
}

pub(super) fn attach(result: String, designed: &Designed) -> Result<String> {
    let prefix = result.strip_suffix("}\n").ok_or_else(|| bad("internal result framing mismatch"))?;
    let failed = match designed.failed_lower {
        None => "null".to_string(),
        Some((h, temperature)) => format!("{{\"htc_w_m2_k\":{},\"mean_temperature_k\":{}}}", num(h)?, num(temperature)?),
    };
    let trials = designed.history.iter().map(|&(h, temperature)| {
        Ok(format!("{{\"htc_w_m2_k\":{},\"mean_temperature_k\":{}}}", num(h)?, num(temperature)?))
    }).collect::<Result<Vec<_>>>()?.join(",");
    Ok(format!("{prefix},\"design\":{{\"surface\":{},\"limit_k\":{},\"selected_htc_w_m2_k\":{},\"status\":{},\"failed_lower\":{},\"log_bracket_width\":{},\"evaluations\":{},\"total_solid_solves\":{},\"history\":[{}],\"search_claim\":\"evaluated passing endpoint of a local target bracket, or feasible declared minimum; no global optimality or infeasibility proof\"}}}}\n",
        quote(&designed.surface), num(designed.limit)?, num(designed.passing_h)?,
        quote(if designed.failed_lower.is_none() { "minimum-feasible" } else { "target-bracketed" }),
        failed, num(designed.log_width)?, designed.history.len(), designed.solid_solves, trials))
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::{request, with_cx, close};

    fn design() -> DesignRequest {
        DesignRequest { surface: "last-face".into(), limit: 323.0, minimum: 10.0, maximum: 1000.0,
            temperature_tolerance: 1e-5, log_tolerance: 1e-4, evaluations: 80 }
    }

    #[test]
    fn file_driven_sizing_returns_a_passing_evaluated_solid() {
        let r = request();
        with_cx(|cx| {
            let flow = r.flow(cx).unwrap();
            let d = design();
            let result = solve(&r, cx, &flow, &d).unwrap();
            assert!(result.passing.objective <= d.limit);
            assert!(d.limit - result.passing.objective <= d.temperature_tolerance);
            assert!(result.log_width <= d.log_tolerance);
            assert!(result.failed_lower.unwrap().1 > d.limit);
            close(result.passing_h, 194.501441754, 0.03);
            assert!(result.history.len() <= d.evaluations);
            let rendered = attach(render(&r, &flow, &result.passing).unwrap(), &result).unwrap();
            let doc = J::parse(&rendered).unwrap();
            let summary = doc.get("design").unwrap();
            close(summary.f64_field("selected_htc_w_m2_k").unwrap(), result.passing_h, 1e-10);
            let wall = &doc.get("walls").unwrap().as_array().unwrap()[1];
            close(wall.f64_field("htc_w_m2_k").unwrap(), result.passing_h, 1e-10);
        });
    }

    #[test]
    fn design_minimum_bracket_and_budget_have_distinct_outcomes() {
        let r = request();
        with_cx(|cx| {
            let flow = r.flow(cx).unwrap();
            let mut d = design();
            d.limit = 400.0;
            d.evaluations = 1;
            let result = solve(&r, cx, &flow, &d).unwrap();
            assert!(result.failed_lower.is_none());
            assert_eq!(result.history.len(), 1);
            d = design(); d.limit = 290.0;
            assert_eq!(solve(&r, cx, &flow, &d).err().unwrap().code, "cooling-network-design-bracket");
            d = design(); d.evaluations = 1;
            assert_eq!(solve(&r, cx, &flow, &d).err().unwrap().code, "cooling-network-design-budget");
        });
    }

    #[test]
    fn optional_design_request_is_admitted_and_not_silently_ignored() {
        let text = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/size-mixed-slab.json"));
        let r = Request::parse(text).unwrap();
        assert_eq!(r.design.as_ref().unwrap().surface, "last-face");
        assert!(Request::parse(&text.replace("\"gradient\": true", "\"gradient\": false")).is_err());
        assert!(Request::parse(&text.replace("\"max_evaluations\": 80", "\"max_evaluations\": 0")).is_err());
    }
}
