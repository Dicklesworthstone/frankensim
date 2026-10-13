//! Finite operating-case envelopes for one shared cooling design.
//!
//! Cases vary declared source loads and external supply temperatures only.
//! Geometry, fluid properties and hydraulic design controls are shared. The
//! nominal case is always included. An unnormalized log-sum-exp upper envelope
//! closes the existing ONE-inequality optimizer without averaging violations.
//! Every case's actual selected temperature must also pass at replay. This is
//! a finite deterministic case set, not probabilistic or continuum robustness.

use super::*;
use fs_conduction::{ComponentPower, PowerAudit, PowerMap, PowerUncertainty};

#[derive(Debug, Clone)]
struct Loads {
    uniform: f64,
    nodal: Option<ScalarField>,
    audit: Option<PowerAudit>,
    map: Option<PowerMap>,
    inlets: Vec<TransportInlet>,
    objective: objective::Objective,
}

impl Loads {
    fn capture(request: &Request) -> Self {
        Self { uniform: request.source, nodal: request.solid_data.nodal_source.clone(),
            audit: request.solid_data.power.clone(), map: request.solid_data.component_map.clone(),
            inlets: request.inlets.clone(), objective: request.objective.clone() }
    }

    fn install(&self, request: &mut Request) {
        request.source = self.uniform;
        request.solid_data.nodal_source.clone_from(&self.nodal);
        request.solid_data.power.clone_from(&self.audit);
        request.solid_data.component_map.clone_from(&self.map);
        request.inlets.clone_from(&self.inlets);
        request.objective = self.objective.clone();
    }

    // Absolute component watts are a new declared load case. Unlisted rows
    // keep their original watts; changed rows do not inherit an uncertainty
    // band whose meaning the case never declared.
    fn with_component_powers(&self, value: &J, mesh: &ConductionMesh, tolerance: f64) -> Result<Self> {
        let map = self.map.as_ref().ok_or_else(|| bad("component_powers_w requires a component-power model"))?;
        let entries = value.as_object().ok_or_else(|| bad("component_powers_w must be a name-to-watts object"))?;
        if entries.is_empty() || entries.len() > map.components().len() {
            return Err(bad("component_powers_w must name a nonempty subset of the original components"));
        }
        let mut replacements = BTreeMap::new();
        for (name, value) in entries {
            let watts = number(value, "component case watts")?;
            if watts < 0.0 || !map.components().iter().any(|row| row.name() == name.as_str()) {
                return Err(bad("case component must exist and its watts must be nonnegative"));
            }
            if replacements.insert(name.as_str(), watts).is_some() {
                return Err(bad("case component is repeated"));
            }
        }
        let mut rows = Vec::with_capacity(map.components().len());
        let mut total = 0.0;
        for row in map.components() {
            let changed = replacements.get(row.name()).copied();
            let watts = changed.unwrap_or(row.watts());
            let uncertainty = if changed.is_some() { PowerUncertainty::Unstated } else { row.uncertainty() };
            rows.push(ComponentPower::new(row.name(), watts, uncertainty, row.vertices().to_vec()).map_err(producer)?);
            total = finite(total + watts)?;
        }
        let map = PowerMap::new(rows, total).map_err(producer)?;
        let (source, audit) = map.volumetric_source(mesh, tolerance).map_err(producer)?;
        let mut result = self.clone();
        result.uniform = 0.0;
        result.nodal = Some(source);
        result.audit = Some(audit);
        result.map = Some(map);
        Ok(result)
    }

    fn scaled(&self, factor: f64, mesh: &ConductionMesh, tolerance: f64) -> Result<Self> {
        let mut result = self.clone();
        if factor == 1.0 { return Ok(result); }
        if let Some(map) = &self.map {
            let mut rows = Vec::with_capacity(map.components().len());
            let mut total = 0.0;
            for row in map.components() {
                let watts = finite(row.watts() * factor)?;
                let uncertainty = match row.uncertainty() {
                    PowerUncertainty::Unstated => PowerUncertainty::Unstated,
                    PowerUncertainty::HalfWidth { half_width, confidence } =>
                        PowerUncertainty::HalfWidth { half_width: finite(half_width * factor)?, confidence },
                };
                rows.push(ComponentPower::new(row.name(), watts, uncertainty, row.vertices().to_vec()).map_err(producer)?);
                total = finite(total + watts)?;
            }
            let map = PowerMap::new(rows, total).map_err(producer)?;
            let (source, audit) = map.volumetric_source(mesh, tolerance).map_err(producer)?;
            result.uniform = 0.0;
            result.nodal = Some(source);
            result.audit = Some(audit);
            result.map = Some(map);
        } else {
            result.uniform = finite(self.uniform * factor)?;
            result.nodal = self.nodal.as_ref().map(|field| match field {
                ScalarField::Uniform(value) => finite(value * factor).map(ScalarField::Uniform),
                ScalarField::Nodal(values) => values.iter().map(|value| finite(value * factor))
                    .collect::<Result<Vec<_>>>().map(ScalarField::Nodal),
            }).transpose()?;
        }
        Ok(result)
    }
}

#[derive(Debug)]
struct Case { name: String, loads: Loads, limit: f64 }

#[derive(Debug)]
pub(super) struct Cases { smoothing_k: f64, cases: Vec<Case> }

impl Cases {
    pub(super) fn parse(value: &J, request: &Request, limit: f64, peak_smoothing: Option<f64>) -> Result<Self> {
        object(value, &["smoothing_k", "power_relative_tolerance", "cases"], "operating_envelope")?;
        let smoothing_k = positive(get(value, "smoothing_k")?, "operating_envelope.smoothing_k")?;
        let tolerance = number(get(value, "power_relative_tolerance")?, "power_relative_tolerance")?;
        if !(0.0..1.0).contains(&tolerance) { return Err(bad("power_relative_tolerance must be in [0,1)")); }
        let rows = array(get(value, "cases")?, "operating_envelope.cases", 7)?;
        if rows.is_empty() { return Err(bad("operating_envelope requires at least one additional case")); }
        let baseline = Loads::capture(request);
        let mut nominal = baseline.clone();
        prepare_objective(&mut nominal.objective, peak_smoothing)?;
        let mut cases = vec![Case { name: "nominal".into(), loads: nominal, limit }];
        let mut names = BTreeSet::from(["nominal".to_string()]);
        for row in rows {
            object(row, &["name", "source_multiplier", "component_powers_w", "inlet_temperatures",
                "objective", "temperature_limit_k", "maximum_smoothing_k"], "operating case")?;
            let name = string(get(row, "name")?, "case.name")?;
            if !names.insert(name.clone()) { return Err(bad("case names must be unique; nominal is reserved")); }
            let factor = row.get("source_multiplier").map(|v| number(v, "source_multiplier")).transpose()?.unwrap_or(1.0);
            if factor < 0.0 { return Err(bad("source_multiplier must be nonnegative")); }
            let mut loads = if let Some(powers) = row.get("component_powers_w") {
                if row.get("source_multiplier").is_some() {
                    return Err(bad("a case must choose source_multiplier or component_powers_w, not both"));
                }
                baseline.with_component_powers(powers, &request.mesh, tolerance)?
            } else { baseline.scaled(factor, &request.mesh, tolerance)? };
            if let Some(inlets) = row.get("inlet_temperatures") {
                let mut nodes = BTreeSet::new();
                for inlet in array(inlets, "inlet_temperatures", baseline.inlets.len())? {
                    object(inlet, &["node", "temperature_k"], "case inlet")?;
                    let node = integer(get(inlet, "node")?, "case inlet node", request.graph.node_count() - 1)?;
                    if !nodes.insert(node) { return Err(bad("case repeats an inlet node")); }
                    let source = loads.inlets.iter_mut().find(|entry| entry.node == node)
                        .ok_or_else(|| bad("case temperatures can change only existing external supplies"))?;
                    source.temperature = fs_qty::Temperature::new(positive(get(inlet, "temperature_k")?, "case temperature_k")?);
                }
            }
            if let Some(selection) = row.get("objective") {
                // This is a requirement selector, not a second gradient switch.
                object(selection, &["mean_wall_region", "max_wall_region", "max_solid_temperature", "max_vertices"], "case objective")?;
                loads.objective = objective::Objective::parse(selection, &request.surfaces, &request.mesh)?;
            }
            let local_smoothing = row.get("maximum_smoothing_k")
                .map(|value| positive(value, "case maximum_smoothing_k")).transpose()?;
            if loads.objective.is_mean() && local_smoothing.is_some() {
                return Err(bad("a case mean cannot declare maximum_smoothing_k"));
            }
            prepare_objective(&mut loads.objective, local_smoothing.or(peak_smoothing))?;
            let limit = row.get("temperature_limit_k").map(|value| positive(value, "case temperature_limit_k"))
                .transpose()?.unwrap_or(limit);
            cases.push(Case { name, loads, limit });
        }
        finite(smoothing_k * det::ln(cases.len() as f64))?;
        Ok(Self { smoothing_k, cases })
    }

    pub(super) fn len(&self) -> usize { self.cases.len() }
}

fn prepare_objective(objective: &mut objective::Objective, smoothing: Option<f64>) -> Result<()> {
    if !objective.is_mean() {
        objective.enable_upper_envelope(smoothing.ok_or_else(|| bad("each peak case requires explicit or inherited maximum_smoothing_k"))?)?;
    }
    Ok(())
}

/// Values and complete sample gradients for EVERY case, retained for replay.
/// Equality of only the combined gradient could conceal cancelling changes.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Signature {
    sample: ProjectedAlSample,
    constraint_temperature: f64,
    actual_temperature: f64,
    limit: f64,
}

pub(super) struct Candidate {
    pub sample: ProjectedAlSample,
    pub flow: GraphSolution,
    pub physical: Evaluation,
    pub signatures: Vec<Signature>,
    pub constraint_temperature: f64,
    pub actual_temperature: f64,
    pub passes: bool,
    pub solid_solves: usize,
    pub report: Option<String>,
}

/// Stable unnormalized upper envelope and derivative weights. No averaging by
/// case count, no subtraction of the bias, no normalization of a violation.
fn envelope(values: &[f64], tau: f64) -> Result<(f64, Vec<f64>)> {
    if values.is_empty() || !tau.is_finite() || tau <= 0.0
        || values.iter().any(|v| !v.is_finite())
    { return Err(bad("invalid operating envelope values or smoothing")); }
    let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut weights = Vec::with_capacity(values.len());
    let mut total = 0.0;
    for &value in values {
        // Negative overflow in the scaled exponent is harmless: exp(-inf)=0.
        // Form the difference first; the retained physical temperatures are finite.
        let weight = det::exp((value - maximum) / tau);
        total = finite(total + weight)?;
        weights.push(weight);
    }
    if total < 1.0 { return Err(producer("operating envelope lost its maximum term")); }
    for weight in &mut weights { *weight /= total; }
    Ok((finite(maximum + tau * det::ln(total))?, weights))
}

/// A complete case sweep. Sources always come from the retained nominal model,
/// never from the previous case/trial. On error the caller publishes nothing.
/// The nominal loads are restored on both success and refusal.
pub(super) fn evaluate_candidate(
    request: &mut Request, original: &LossGraph, design: &Specification,
    point: &[f64], cx: &Cx<'_>, render_cases: bool,
) -> Result<Candidate> {
    let result = evaluate_inner(request, original, design, point, cx, render_cases);
    if let Some(cases) = &design.operating_envelope { cases.cases[0].loads.install(request); }
    result
}

fn evaluate_inner(
    request: &mut Request, original: &LossGraph, design: &Specification,
    point: &[f64], cx: &Cx<'_>, render_cases: bool,
) -> Result<Candidate> {
    let count = design.case_count();
    let mut signatures = Vec::with_capacity(count);
    let mut nominal = None;
    let mut physical_max: f64 = f64::NEG_INFINITY;
    let mut all_pass = true;
    let mut solid_solves = 0_usize;
    let mut rendered = Vec::new();
    for index in 0..count {
        poll(cx)?;
        let (name, limit) = if let Some(cases) = &design.operating_envelope {
            let case = &cases.cases[index];
            case.loads.install(request);
            (case.name.as_str(), case.limit)
        } else { ("nominal", design.limit) };
        let (mut sample, flow, physical) = evaluate(request, original, design, point, cx)?;
        sample.constraint = finite((physical.objective - limit) / design.scale)?;
        let actual = physical.objective_state.physical_value();
        physical_max = physical_max.max(actual);
        all_pass &= physical.objective <= limit && actual <= limit;
        solid_solves = solid_solves.checked_add(physical.coupled.iterations)
            .ok_or_else(|| bad("operating-case work count overflow"))?;
        signatures.push(Signature { sample, constraint_temperature: physical.objective, actual_temperature: actual, limit });
        if render_cases && design.operating_envelope.is_some() {
            let mut output = render(request, &flow, &physical)?;
            if let Some(fan) = &request.fan { output = fan.attach(output, &flow, fan.speed_ratio)?; }
            let inlets = request.inlets.iter().map(|inlet| Ok(format!(
                "{{\"node\":{},\"temperature_k\":{}}}", inlet.node, num(inlet.temperature.value())?
            ))).collect::<Result<Vec<_>>>()?.join(",");
            rendered.push(format!("{{\"name\":{},\"temperature_limit_k\":{},\"inlets\":[{inlets}],\"result\":{}}}",
                quote(name), num(limit)?, output.trim_end()));
        }
        if index == 0 { nominal = Some((flow, physical)); }
    }
    let mut sample = signatures[0].sample.clone();
    let mut constraint_temperature = signatures[0].constraint_temperature;
    let report = if let Some(cases) = &design.operating_envelope {
        // Only heat sources and supply temperatures vary. Hydraulic air power
        // and design-change cost are therefore identical and charged ONCE.
        for signature in &signatures[1..] {
            if signature.sample.objective != sample.objective || signature.sample.gradient != sample.gradient {
                return Err(producer("operating cases changed the shared hydraulic design cost"));
            }
        }
        // Different requirements are compared by excess over THEIR OWN
        // limits. An absolute-temperature maximum would hide a cooler but
        // tighter-limited component, and averaging would hide any violation.
        let values = signatures.iter().map(|s| finite(s.constraint_temperature - s.limit))
            .collect::<Result<Vec<_>>>()?;
        let (excess, weights) = envelope(&values, cases.smoothing_k)?;
        constraint_temperature = finite(design.limit + excess)?;
        sample.constraint = finite(excess / design.scale)?;
        sample.constraint_gradient.fill(0.0);
        for (signature, weight) in signatures.iter().zip(&weights) {
            poll(cx)?;
            for (gradient, partial) in sample.constraint_gradient.iter_mut().zip(&signature.sample.constraint_gradient) {
                *gradient = finite(*gradient + finite(weight * partial)?)?;
            }
        }
        all_pass &= excess <= 0.0;
        if render_cases {
            Some(format!("{{\"reported_top_level_case\":\"nominal\",\"smoothing_k\":{},\"constraint_temperature_k\":{},\"constraint_excess_k\":{},\"maximum_actual_temperature_k\":{},\"analytic_case_bias_bound_k\":{},\"case_weights\":{},\"scaled_constraint_gradient\":{},\"cases\":[{}],\"scope\":\"finite declared operating cases, nominal always included; unnormalized upper envelope of own-limit excesses; constraint_temperature_k is nominal-limit plus excess, not a case temperature; not a probability, unsampled-load guarantee or continuum bound\"}}",
                num(cases.smoothing_k)?, num(constraint_temperature)?, num(excess)?, num(physical_max)?,
                num(cases.smoothing_k * det::ln(count as f64))?, numbers(&weights)?, numbers(&sample.constraint_gradient)?, rendered.join(",")))
        } else { None }
    } else { None };
    let (flow, physical) = nominal.ok_or_else(|| bad("operating sweep omitted nominal case"))?;
    poll(cx)?;
    Ok(Candidate { sample, flow, physical, signatures, constraint_temperature,
        actual_temperature: physical_max, passes: all_pass, solid_solves, report })
}

#[cfg(test)]
mod tests;
