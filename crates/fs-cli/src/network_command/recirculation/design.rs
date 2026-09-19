//! One return-fraction temperature target, evaluated through the complete model.
//! No monotonicity is assumed: the passing side can be at either end, and a
//! local bracket is not a global optimum or an infeasibility certificate.
use super::*;
use fan_speed::ThermalEvaluation;

#[derive(Debug)]
pub(super) struct Design {
    input: J,
    supply: usize,
    returned: usize,
    minimum: f64,
    maximum: f64,
    limit: f64,
    fraction_tolerance: f64,
    temperature_tolerance: f64,
    max_evaluations: usize,
}

impl Design {
    pub(super) fn parse(value: &J, root: &J, policy: &Policy) -> Result<Self> {
        object(value, &["supply_node", "return_node", "min_fraction", "max_fraction",
            "temperature_limit_k", "fraction_tolerance", "temperature_tolerance_k",
            "max_evaluations"], "recirculation.design")?;
        if ["transient", "mesh_convergence", "design", "fan_speed_design"].iter()
            .any(|key| root.get(key).is_some()) {
            return Err(bad("return-fraction design requires a steady request without another design or mesh study"));
        }
        let supply = integer_raw(get(value, "supply_node")?, "design.supply_node")?;
        let returned = integer_raw(get(value, "return_node")?, "design.return_node")?;
        if !policy.links.iter().any(|link| (link.supply_node, link.return_node) == (supply, returned)) {
            return Err(bad("return-fraction design must select one configured link"));
        }
        let minimum = number(get(value, "min_fraction")?, "min_fraction")?;
        let maximum = number(get(value, "max_fraction")?, "max_fraction")?;
        if !(0.0 <= minimum && minimum < maximum && maximum < 1.0) {
            return Err(bad("return design needs 0 <= min_fraction < max_fraction < 1"));
        }
        let mut total = 0.0;
        for link in policy.links.iter().filter(|link| link.supply_node == supply) {
            total += if link.return_node == returned { maximum } else { link.fraction };
        }
        if !total.is_finite() || total >= 1.0 {
            return Err(bad("the maximum designed fraction must preserve positive fresh makeup"));
        }
        let max_evaluations = count(get(value, "max_evaluations")?, "return design max_evaluations", 256)?;
        if max_evaluations < 2 { return Err(bad("return design needs at least two endpoint evaluations")); }
        let mut input = root.clone();
        members(member_mut(&mut input, "recirculation")?)?.retain(|(key, _)| key != "design");
        Ok(Self { input, supply, returned, minimum, maximum, max_evaluations,
            limit: positive(get(value, "temperature_limit_k")?, "return design temperature_limit_k")?,
            fraction_tolerance: positive(get(value, "fraction_tolerance")?, "fraction_tolerance")?,
            temperature_tolerance: positive(get(value, "temperature_tolerance_k")?, "return design temperature_tolerance_k")? })
    }

    fn candidate(&self, fraction: f64) -> Result<Request> {
        Request::parse(&encode(&self.resolved_input(fraction)?)?)
    }

    pub(super) fn solve(&self, cx: &Cx<'_>) -> Result<String> {
        let mut evaluator = Evaluator { design: self, history: Vec::new(), solid_solves: 0 };
        let lower = evaluator.trial(cx, self.minimum)?;
        let upper = evaluator.trial(cx, self.maximum)?;
        let lower_passes = lower.observation.temperature <= self.limit;
        let upper_passes = upper.observation.temperature <= self.limit;
        if lower_passes && upper_passes {
            return self.finish(cx, upper, None, &evaluator.history, evaluator.solid_solves);
        }
        if !lower_passes && !upper_passes {
            return Err(Failure { code: "cooling-network-design-bracket",
                message: "neither declared return-fraction endpoint meets the temperature limit; no interior infeasibility claim".into() });
        }
        let (mut passing, mut failed) = if lower_passes { (lower, upper) } else { (upper, lower) };
        loop {
            poll(cx)?;
            let low = passing.observation.fraction.min(failed.observation.fraction);
            let high = passing.observation.fraction.max(failed.observation.fraction);
            if high-low <= self.fraction_tolerance
                && self.limit-passing.observation.temperature <= self.temperature_tolerance {
                return self.finish(cx, passing, Some(failed.observation), &evaluator.history, evaluator.solid_solves);
            }
            let mut ends = [passing.observation, failed.observation];
            if (ends[1].temperature-self.limit).abs() < (ends[0].temperature-self.limit).abs() { ends.swap(0,1); }
            let next = ends.into_iter().find_map(|value| proposal(low, high, value, self.limit))
                .unwrap_or(0.5*low+0.5*high);
            if !(next > low && next < high) {
                return Err(Failure { code: "cooling-network-design-resolution",
                    message: "return-fraction resolution cannot meet both design tolerances".into() });
            }
            let candidate = evaluator.trial(cx, next)?;
            if candidate.observation.temperature <= self.limit { passing = candidate; }
            else { failed = candidate; }
        }
    }

    fn finish(&self, cx: &Cx<'_>, passing: Candidate, failed: Option<Observation>,
        history: &[Observation], solid_solves: usize) -> Result<String> {
        poll(cx)?;
        if passing.observation.temperature > self.limit { return Err(producer("return design attempted to publish a failing candidate")); }
        let mut output = passing.thermal.render(&passing.request, &passing.flow)?;
        if let Some(fan) = &passing.request.fan {
            output = fan.attach(output, &passing.flow, fan.speed_ratio)?;
        }
        let prefix = output.strip_suffix("}\n").ok_or_else(|| bad("return design result framing"))?;
        let trials = history.iter().map(Observation::render).collect::<Result<Vec<_>>>()?.join(",");
        let width = failed.map_or(0.0, |bad| (bad.fraction-passing.observation.fraction).abs());
        let failed_json = failed.map(|value| value.render()).transpose()?.unwrap_or_else(|| "null".into());
        let resolved = encode(&self.resolved_input(passing.observation.fraction)?)?;
        let result = format!("{prefix},\"recirculation_design\":{{\"supply_node\":{},\"return_node\":{},\"selected_fraction\":{},\"temperature_limit_k\":{},\"status\":{},\"failed_endpoint\":{failed_json},\"fraction_bracket_width\":{},\"evaluations\":{},\"total_solid_solves\":{},\"history\":[{trials}],\"resolved_request\":{resolved},\"search_claim\":\"actual passing endpoint of a local temperature-threshold bracket, or feasible declared maximum when both bounds pass; either slope sign is supported; derivatives only propose candidates; no global monotonicity, unique root, maximum-return or hardware-compliance certificate\"}}}}\n",
            self.supply, self.returned, num(passing.observation.fraction)?, num(self.limit)?,
            quote(if failed.is_some() { "target-bracketed" } else { "both-bounds-feasible" }),
            num(width)?, history.len(), solid_solves);
        poll(cx)?;
        Ok(result)
    }

    fn resolved_input(&self, fraction: f64) -> Result<J> {
        let mut input = self.input.clone();
        let J::Array(links) = member_mut(member_mut(&mut input, "recirculation")?, "links")?
            else { return Err(bad("return design lost its links")); };
        let mut changed = 0;
        for link in links {
            if integer_raw(get(link, "supply_node")?, "supply_node")? == self.supply
                && integer_raw(get(link, "return_node")?, "return_node")? == self.returned {
                *member_mut(link, "fraction")? = J::Number { value: fraction, raw: num(fraction)? };
                changed += 1;
            }
        }
        if changed != 1 { return Err(bad("return design target was not unique")); }
        Ok(input)
    }
}

struct Evaluator<'a> { design: &'a Design, history: Vec<Observation>, solid_solves: usize }
impl Evaluator<'_> {
    fn trial(&mut self, cx: &Cx<'_>, fraction: f64) -> Result<Candidate> {
        poll(cx)?;
        if self.history.len() >= self.design.max_evaluations {
            return Err(Failure { code: "cooling-network-design-budget",
                message: "return-fraction evaluation budget exhausted; no partial design published".into() });
        }
        let request = self.design.candidate(fraction)?;
        poll(cx)?;
        let flow = request.flow(cx)?;
        let h = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let thermal = match &request.radiation {
            Some(policy) => policy.evaluate(&request, cx, &flow, &h, request.gradient)?,
            None => {
                let value = request.evaluate(cx, &flow, &h, request.gradient)?;
                let solid_solves = value.coupled.iterations;
                ThermalEvaluation { value, radiation: None, solid_solves }
            }
        };
        let policy = request.recirculation.as_ref().ok_or_else(|| bad("return design lost its policy"))?;
        let index = policy.links.iter().position(|link|
            (link.supply_node, link.return_node) == (self.design.supply, self.design.returned))
            .ok_or_else(|| bad("return design lost its target"))?;
        let slope = policy.gradients(&request, &thermal.value)?.map(|values| values[index]);
        let observation = Observation { fraction, temperature: thermal.value.objective, slope };
        if !observation.temperature.is_finite() { return Err(producer("nonfinite return design objective")); }
        self.solid_solves = self.solid_solves.checked_add(thermal.solid_solves)
            .ok_or_else(|| bad("return design solid-work overflow"))?;
        self.history.push(observation);
        poll(cx)?;
        Ok(Candidate { request, flow, thermal, observation })
    }
}

struct Candidate { request: Request, flow: GraphSolution, thermal: ThermalEvaluation, observation: Observation }
#[derive(Clone, Copy)]
struct Observation { fraction: f64, temperature: f64, slope: Option<f64> }
impl Observation {
    fn render(&self) -> Result<String> {
        Ok(format!("{{\"fraction\":{},\"temperature_k\":{},\"dtemperature_dfraction_k\":{}}}",
            num(self.fraction)?, num(self.temperature)?, optional(self.slope)?))
    }
}
fn proposal(low: f64, high: f64, value: Observation, limit: f64) -> Option<f64> {
    let slope = value.slope.filter(|slope| slope.is_finite() && *slope != 0.0)?;
    let next = value.fraction - (value.temperature-limit)/slope;
    let guard = 0.1*(high-low);
    (next.is_finite() && next > low+guard && next < high-guard).then_some(next)
}
fn members(value: &mut J) -> Result<&mut Vec<(String,J)>> {
    if let J::Object(fields) = value { Ok(fields) } else { Err(bad("expected return design object")) }
}
fn member_mut<'a>(value: &'a mut J, key: &str) -> Result<&'a mut J> {
    members(value)?.iter_mut().find_map(|(name,value)| (name==key).then_some(value))
        .ok_or_else(|| bad(format!("missing return design field {key}")))
}
fn encode(value: &J) -> Result<String> {
    fn write(value: &J, out: &mut String) -> Result<()> {
        match value {
            J::Null => out.push_str("null"), J::Bool(v) => out.push_str(if *v {"true"} else {"false"}),
            J::Number {value,raw} => { if !value.is_finite() { return Err(bad("nonfinite return design input")); } out.push_str(raw); },
            J::Str(s) => out.push_str(&quote(s)),
            J::Array(values) => {
                out.push('['); for (i,v) in values.iter().enumerate() { if i>0 {out.push(',');} write(v,out)?; } out.push(']');
            }
            J::Object(fields) => {
                out.push('{'); for (i,(k,v)) in fields.iter().enumerate() { if i>0 {out.push(',');} out.push_str(&quote(k)); out.push(':'); write(v,out)?; } out.push('}');
            }
        }
        if out.len() as u64 > MAX_INPUT_BYTES { return Err(bad("resolved return request exceeds input limit")); }
        Ok(())
    }
    let mut out=String::new(); write(value,&mut out)?; Ok(out)
}

#[cfg(test)]
mod tests;
