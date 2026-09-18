//! Resolve independent absolute-watt controls without changing their footprints.
use super::*;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub(super) struct Axis {
    pub component: String,
    pub interval: usize,
    pub minimum: f64,
    pub maximum: f64,
}

pub(super) struct Plan {
    pub base: J,
    pub axes: Vec<Axis>,
    pub limit: f64,
    pub power_tolerance: f64,
    pub temperature_tolerance: f64,
    pub max_evaluations: usize,
    pub max_total_steps: usize,
    pub wall_seconds: f64,
    pub planned_steps: usize,
    pub cycles: usize,
    pub intervals: usize,
    pub adjoint: bool,
    pub qoi: model::Qoi,
}

impl Plan {
    pub fn parse(base: &J, spec: &J) -> Result<Self> {
        object(spec, &["schema", "units", "temperature_limit_k", "power_tolerance_w",
            "temperature_tolerance_k", "max_evaluations", "max_total_steps", "wall_seconds",
            "priority"], "component design")?;
        if spec.str_field("schema") != Some("frankensim.cooling-component-design.v1")
            || spec.str_field("units") != Some("SI")
            || base.str_field("schema") != Some("frankensim.cooling-network.v1")
            || base.str_field("units") != Some("SI") {
            return Err(bad("component design and cooling base require their v1 schemas and SI units"));
        }
        for key in ["design", "fan_speed_design", "mesh_convergence"] {
            if base.get(key).is_some() { return Err(bad("component allocation does not nest other design or mesh studies")); }
        }
        if base.path(&["objective", "gradient"]) != Some(&J::Bool(false)) {
            return Err(bad("component allocation requires a transient base with objective.gradient=false"));
        }
        let schedule = field(base, "transient")?;
        for key in ["adaptive", "time_convergence", "fan_speed_design", "power_design"] {
            if schedule.get(key).is_some() { return Err(bad("component allocation requires fixed timesteps without nested searches")); }
        }
        if schedule.get("repeat").is_some_and(|r| r.get("fan_controller").is_some()) {
            return Err(bad("component allocation requires a fixed fan schedule, not a thermostat"));
        }
        let qoi_spec = J::Object(vec![("kind".into(), J::Str("transient-sampled-peak".into()))]);
        let qoi = model::Qoi::parse(Some(&qoi_spec), base)?;
        let cycles = match qoi {
            model::Qoi::TransientPeak { repeated_cycles } => repeated_cycles.unwrap_or(1),
            _ => return Err(bad("component allocation requires a complete transient")),
        };
        let adjoint = match schedule.get("adjoint") {
            None => false,
            Some(policy) => {
                if policy.str_field("qoi") != Some("sampled-peak")
                    || policy.get("component_power") != Some(&J::Bool(true)) {
                    return Err(bad("allocation adjoints require qoi=sampled-peak and component_power=true"));
                }
                true
            }
        };
        let intervals = array(field(schedule, "intervals")?, "intervals", 4096)?;
        if intervals.is_empty() { return Err(bad("empty workload schedule")); }
        let max_dt = positive(field(schedule, "max_step_s")?, "max_step_s")?;
        let max_steps = positive_count(field(schedule, "max_steps")?, "max_steps", 10_000)?;
        let mut cycle_steps = 0usize;
        for interval in intervals {
            let duration = positive(field(interval, "duration_s")?, "duration_s")?;
            let required = (duration / max_dt).ceil().max(1.0);
            if !required.is_finite() || required > max_steps as f64 {
                return Err(budget("base trajectory exceeds its per-cycle step cap"));
            }
            let steps = interval.get("steps").map(|v| positive_count(v, "interval.steps", max_steps))
                .transpose()?.unwrap_or(required as usize);
            if steps < required as usize { return Err(bad("interval.steps violates max_step_s")); }
            cycle_steps = cycle_steps.checked_add(steps).ok_or_else(|| budget("step count overflow"))?;
        }
        if cycle_steps > max_steps { return Err(budget("base trajectory exceeds max_steps")); }
        let planned_steps = cycle_steps.checked_mul(cycles).ok_or_else(|| budget("cycle step count overflow"))?;
        if let Some(repeat) = schedule.get("repeat") {
            let cap = positive_count(field(repeat, "max_total_steps")?, "repeat.max_total_steps", 1_000_000)?;
            if planned_steps > cap { return Err(budget("base trajectory exceeds its repeated-step cap")); }
        }
        let mut component_watts = BTreeMap::new();
        let components = base.path(&["solid", "component_power", "components"])
            .ok_or_else(|| bad("allocation requires solid.component_power with fixed footprints"))?;
        for row in array(components, "components", 4096)? {
            let name = string(field(row, "name")?, "component name")?;
            let watts = nonnegative(field(row, "watts")?, "component watts")?;
            if component_watts.insert(name, watts).is_some() { return Err(bad("duplicate base component name")); }
        }
        if component_watts.is_empty() { return Err(bad("empty component map")); }
        let rows = array(field(spec, "priority")?, "priority", 64)?;
        if rows.is_empty() { return Err(bad("priority requires at least one component/interval control")); }
        let mut axes = Vec::new();
        let mut seen = BTreeSet::new();
        let mut touched = BTreeSet::new();
        for row in rows {
            object(row, &["component", "interval", "min_power_w", "max_power_w"], "priority control")?;
            let component = string(field(row, "component")?, "priority.component")?;
            if !component_watts.contains_key(&component) { return Err(bad("priority names an unknown component")); }
            let interval = integer(field(row, "interval")?, "priority.interval", intervals.len()-1)?;
            if !seen.insert((interval, component.clone())) { return Err(bad("duplicate component/interval control")); }
            let minimum = nonnegative(field(row, "min_power_w")?, "min_power_w")?;
            let maximum = positive(field(row, "max_power_w")?, "max_power_w")?;
            if minimum >= maximum { return Err(bad("component bounds require 0 <= minimum < maximum")); }
            touched.insert(interval);
            axes.push(Axis { component, interval, minimum, maximum });
        }
        if touched.len().checked_mul(component_watts.len()).is_none_or(|n| n > 65_536) {
            return Err(budget("expanded component schedules exceed 65536 workload rows"));
        }
        // Convert only controlled intervals. Preserve every other component's
        // actual watts, including explicit overrides and dormant components.
        let mut resolved = base.clone();
        let schedule_out = member_mut(&mut resolved, "transient")?;
        let output_intervals = array_mut(member_mut(schedule_out, "intervals")?)?;
        for ordinal in touched {
            let original = &intervals[ordinal];
            let powers = match (original.get("power_scale"), original.get("component_powers_w")) {
                (Some(scale), None) => {
                    let scale = nonnegative(scale, "power_scale")?;
                    component_watts.iter().map(|(name, watts)|
                        Ok((name.clone(), number_value(checked(watts*scale)?)?)))
                        .collect::<Result<Vec<_>>>()?
                }
                (None, Some(powers)) => {
                    let fields = powers.as_object().ok_or_else(|| bad("component_powers_w must be an object"))?;
                    let mut names = BTreeSet::new();
                    for (name, watts) in fields {
                        nonnegative(watts, "interval component watts")?;
                        if !component_watts.contains_key(name) || !names.insert(name.clone()) {
                            return Err(bad("absolute workload names must match the base component map"));
                        }
                    }
                    if names.len() != component_watts.len() { return Err(bad("absolute workload omits a base component")); }
                    fields.to_vec()
                }
                _ => return Err(bad("each interval requires exactly one power_scale or component_powers_w")),
            };
            remove(&mut output_intervals[ordinal], "power_scale")?;
            put(&mut output_intervals[ordinal], "component_powers_w", J::Object(powers))?;
        }
        let limit = positive(field(spec, "temperature_limit_k")?, "temperature_limit_k")?;
        put(schedule_out, "temperature_limit_k", number_value(limit)?)?;
        let wall_seconds = positive(field(spec, "wall_seconds")?, "wall_seconds")?;
        if wall_seconds > 86400.0 { return Err(bad("component design wall_seconds exceeds one day")); }
        Ok(Self {
            base: resolved, axes, limit,
            power_tolerance: positive(field(spec, "power_tolerance_w")?, "power_tolerance_w")?,
            temperature_tolerance: positive(field(spec, "temperature_tolerance_k")?, "temperature_tolerance_k")?,
            max_evaluations: positive_count(field(spec, "max_evaluations")?, "max_evaluations", 4096)?,
            max_total_steps: positive_count(field(spec, "max_total_steps")?, "max_total_steps", 100_000_000)?,
            wall_seconds, planned_steps, cycles, intervals: intervals.len(), adjoint, qoi,
        })
    }

    pub fn request(&self, powers: &[f64]) -> Result<J> {
        if powers.len() != self.axes.len() { return Err(bad("allocation vector length mismatch")); }
        let mut request = self.base.clone();
        let rows = array_mut(member_mut(member_mut(&mut request, "transient")?, "intervals")?)?;
        for (axis, &watts) in self.axes.iter().zip(powers) {
            if !watts.is_finite() || watts < axis.minimum || watts > axis.maximum {
                return Err(bad("allocation trial is outside its declared watt bounds"));
            }
            let map = member_mut(&mut rows[axis.interval], "component_powers_w")?;
            *member_mut(map, &axis.component)? = number_value(watts)?;
        }
        Ok(request)
    }
}

pub(super) fn positive_count(value: &J, name: &str, cap: usize) -> Result<usize> {
    let n = integer(value, name, cap)?;
    if n == 0 { Err(bad(format!("{name} must be positive"))) } else { Ok(n) }
}
fn nonnegative(value: &J, name: &str) -> Result<f64> {
    let value = number(value, name)?;
    if value < 0.0 { Err(bad(format!("{name} must be nonnegative"))) } else { Ok(value) }
}
pub(super) fn member_mut<'a>(value: &'a mut J, key: &str) -> Result<&'a mut J> {
    let J::Object(fields) = value else { return Err(bad("object required")); };
    fields.iter_mut().find(|(name,_)| name == key).map(|(_,value)| value)
        .ok_or_else(|| bad(format!("missing field {key}")))
}
fn array_mut(value: &mut J) -> Result<&mut Vec<J>> {
    match value { J::Array(rows) => Ok(rows), _ => Err(bad("array required")) }
}
pub(super) fn put(value: &mut J, key: &str, new: J) -> Result<()> {
    let J::Object(fields) = value else { return Err(bad("object required")); };
    if let Some((_, value)) = fields.iter_mut().find(|(name,_)| name == key) { *value = new; }
    else { fields.push((key.into(),new)); }
    Ok(())
}
fn remove(value: &mut J, key: &str) -> Result<()> {
    let J::Object(fields) = value else { return Err(bad("object required")); };
    fields.retain(|(name,_)| name != key); Ok(())
}
