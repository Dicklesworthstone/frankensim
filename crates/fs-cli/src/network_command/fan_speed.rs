//! Fan-speed target search over actual hydraulic, convection and FEM evaluations.
//! No fixed-h adjoint is mislabeled as dT/d(speed): bounded bisection re-solves
//! the complete model, including every correlation domain check, at each trial.

use super::*;

#[derive(Debug)]
pub(super) struct FanSpeedDesign {
    minimum: f64,
    maximum: f64,
    limit: f64,
    speed_tolerance: f64,
    temperature_tolerance: f64,
    evaluations: usize,
}

impl FanSpeedDesign {
    pub(super) fn parse(value: &J, fan: &fan_drive::FanDrive) -> Result<Self> {
        object(value, &["min_speed_ratio", "max_speed_ratio", "temperature_limit_k",
            "speed_ratio_tolerance", "temperature_tolerance_k", "max_evaluations"], "fan_speed_design")?;
        let design = Self {
            minimum: positive(get(value, "min_speed_ratio")?, "fan_speed_design.min_speed_ratio")?,
            maximum: positive(get(value, "max_speed_ratio")?, "fan_speed_design.max_speed_ratio")?,
            limit: positive(get(value, "temperature_limit_k")?, "fan_speed_design.temperature_limit_k")?,
            speed_tolerance: positive(get(value, "speed_ratio_tolerance")?, "speed_ratio_tolerance")?,
            temperature_tolerance: positive(get(value, "temperature_tolerance_k")?, "fan_speed_design.temperature_tolerance_k")?,
            evaluations: count(get(value, "max_evaluations")?, "fan_speed_design.max_evaluations", 4096)?,
        };
        if design.minimum >= design.maximum { return Err(bad("fan-speed bounds must be strictly ordered")); }
        // Domain admission precedes every expensive hydraulic or solid solve.
        fan.bank(design.minimum)?;
        fan.bank(design.maximum)?;
        Ok(design)
    }
}

#[derive(Debug, Clone)]
struct Trial {
    speed: f64,
    temperature: f64,
    flow: f64,
    active_vertex: Option<usize>,
}
impl Trial {
    fn render(&self) -> Result<String> {
        Ok(format!("{{\"speed_ratio\":{},\"temperature_k\":{},\"flow_m3_s\":{},\"active_vertex\":{}}}",
            num(self.speed)?, num(self.temperature)?, num(self.flow)?,
            self.active_vertex.map_or_else(|| "null".into(), |v| v.to_string())))
    }
}

struct Evaluator<'a> {
    request: &'a Request,
    design: &'a FanSpeedDesign,
    fan: &'a fan_drive::FanDrive,
    history: Vec<Trial>,
    solid_solves: usize,
}
impl Evaluator<'_> {
    fn trial(&mut self, cx: &Cx<'_>, speed: f64) -> Result<(GraphSolution, Evaluation)> {
        poll(cx)?;
        if self.history.len() >= self.design.evaluations {
            return Err(Failure { code: "cooling-network-design-budget",
                message: format!("{} complete fan-speed evaluations exhausted; no partial design published", self.history.len()) });
        }
        let flow = self.fan.solve(cx, &self.request.graph, self.request.limits, speed)?;
        let slots = self.request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let value = self.request.evaluate(cx, &flow, &slots, self.request.gradient)?;
        if !value.objective.is_finite() { return Err(producer("nonfinite fan-speed objective")); }
        self.solid_solves = self.solid_solves.checked_add(value.coupled.iterations)
            .ok_or_else(|| bad("fan-speed work count overflow"))?;
        self.history.push(Trial { speed, temperature: value.objective,
            flow: flow.node_outflows[self.fan.inlet].value(), active_vertex: value.objective_state.vertex });
        poll(cx)?;
        Ok((flow, value))
    }

    fn finish(&self, cx: &Cx<'_>, passing: (GraphSolution, Evaluation), speed: f64,
        failed: Option<&Trial>) -> Result<String> {
        let (flow, value) = passing;
        if value.objective > self.design.limit { return Err(producer("internal fan-speed feasibility mismatch")); }
        let output = self.fan.attach(render(self.request, &flow, &value)?, &flow, speed)?;
        let prefix = output.strip_suffix("}\n").ok_or_else(|| bad("internal result framing mismatch"))?;
        let trials = self.history.iter().map(Trial::render).collect::<Result<Vec<_>>>()?.join(",");
        let lower = failed.map(Trial::render).transpose()?.unwrap_or_else(|| "null".into());
        let width = failed.map_or(0.0, |lo| speed - lo.speed);
        let result = format!("{prefix},\"fan_speed_design\":{{\"selected_speed_ratio\":{},\"temperature_limit_k\":{},\"status\":{},\"failed_lower\":{lower},\"speed_bracket_width\":{},\"evaluations\":{},\"total_solid_solves\":{},\"history\":[{trials}],\"search_claim\":\"evaluated passing endpoint of a local speed bracket or feasible declared minimum; no global monotonicity, minimum-speed or hardware-compliance certificate\"}}}}\n",
            num(speed)?, num(self.design.limit)?, quote(if failed.is_some() { "target-bracketed" } else { "minimum-feasible" }),
            num(width)?, self.history.len(), self.solid_solves);
        poll(cx)?;
        Ok(result)
    }
}

pub(super) fn solve(request: &Request, cx: &Cx<'_>, design: &FanSpeedDesign) -> Result<String> {
    let fan = request.fan.as_ref().ok_or_else(|| bad("fan-speed design requires a fan drive"))?;
    let mut evaluator = Evaluator { request, design, fan, history: Vec::new(), solid_solves: 0 };
    let lower = evaluator.trial(cx, design.minimum)?;
    if lower.1.objective <= design.limit { return evaluator.finish(cx, lower, design.minimum, None); }
    let mut failed = evaluator.history.last().expect("lower evaluated").clone();
    let mut passing = evaluator.trial(cx, design.maximum)?;
    if passing.1.objective > design.limit {
        return Err(Failure { code: "cooling-network-design-bracket", message: format!(
            "neither fan-speed endpoint meets {} K: lower {} K, upper {} K; no interior infeasibility claim is made",
            design.limit, failed.temperature, passing.1.objective) });
    }
    let mut high = design.maximum;
    loop {
        poll(cx)?;
        let width = high - failed.speed;
        if width <= design.speed_tolerance && design.limit - passing.1.objective <= design.temperature_tolerance {
            return evaluator.finish(cx, passing, high, Some(&failed));
        }
        let middle = 0.5 * failed.speed + 0.5 * high;
        if !(middle > failed.speed && middle < high) {
            return Err(Failure { code: "cooling-network-design-resolution",
                message: "fan-speed floating-point resolution cannot meet both design tolerances".into() });
        }
        // Producer/domain failures propagate. They are not hot/cold verdicts.
        let value = evaluator.trial(cx, middle)?;
        if value.1.objective <= design.limit {
            high = middle;
            passing = value;
        } else {
            failed = evaluator.history.last().expect("middle evaluated").clone();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::{with_cx, close};
    const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/size-fan-hotspot.json"));

    #[test]
    fn fan_speed_search_returns_the_passing_recomputed_hotspot_field() {
        let r = Request::parse(FIXTURE).unwrap();
        let output = execute(&r, &CancelGate::new_clock_free()).unwrap();
        let result = J::parse(&output).unwrap();
        let design = result.get("fan_speed_design").unwrap();
        let selected = design.f64_field("selected_speed_ratio").unwrap();
        let temperature = result.path(&["objective", "value_k"]).unwrap().as_f64().unwrap();
        assert!(temperature <= 302.2 && 302.2 - temperature <= 1e-5);
        assert!(design.f64_field("speed_bracket_width").unwrap() <= 1e-4);
        assert!(design.path(&["failed_lower", "temperature_k"]).unwrap().as_f64().unwrap() > 302.2);
        assert!(result.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().all(|t| t.as_f64().unwrap() <= 302.2));
        close(result.path(&["fan", "speed_ratio"]).unwrap().as_f64().unwrap(), selected, 1e-14);
        close(result.path(&["fan", "flow_m3_s"]).unwrap().as_f64().unwrap(), 0.004 * selected, 1e-9);
        with_cx(|cx| {
            let flow = r.fan.as_ref().unwrap().solve(cx, &r.graph, r.limits, selected).unwrap();
            let slots = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            close(r.evaluate(cx, &flow, &slots, false).unwrap().objective, temperature, 1e-8);
        });
    }

    #[test]
    fn fan_speed_minimum_missing_bracket_and_exhaustion_are_distinct() {
        let mut r = Request::parse(FIXTURE).unwrap();
        r.fan_speed_design.as_mut().unwrap().limit = 400.0;
        r.fan_speed_design.as_mut().unwrap().evaluations = 1;
        let result = J::parse(&execute(&r, &CancelGate::new_clock_free()).unwrap()).unwrap();
        assert_eq!(result.path(&["fan_speed_design", "status"]).unwrap().as_str(), Some("minimum-feasible"));
        r.fan_speed_design.as_mut().unwrap().limit = 302.2;
        assert_eq!(execute(&r, &CancelGate::new_clock_free()).unwrap_err().code, "cooling-network-design-budget");
        r.fan_speed_design.as_mut().unwrap().evaluations = 64;
        r.fan_speed_design.as_mut().unwrap().limit = 290.0;
        assert_eq!(execute(&r, &CancelGate::new_clock_free()).unwrap_err().code, "cooling-network-design-bracket");
    }

    #[test]
    fn fan_speed_domain_conflicts_and_cancellation_refuse() {
        let r = Request::parse(FIXTURE).unwrap();
        let invalid = J::parse(r#"{"min_speed_ratio":0.1,"max_speed_ratio":2,"temperature_limit_k":302.2,"speed_ratio_tolerance":0.0001,"temperature_tolerance_k":0.00001,"max_evaluations":64}"#).unwrap();
        assert!(FanSpeedDesign::parse(&invalid, r.fan.as_ref().unwrap()).is_err());
        let gate = CancelGate::new_clock_free(); gate.request();
        assert!(execute(&r, &gate).is_err());
        let both = FIXTURE.replace("\"fan_speed_design\": {", "\"design\": {}, \"fan_speed_design\": {");
        assert!(Request::parse(&both).is_err());
    }
}
