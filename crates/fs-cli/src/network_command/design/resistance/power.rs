//! Nominal fan air-power cost. P = inlet volume flow * network pressure rise.
//! The existing fan/network adjoint closes operating-point feedback. This is
//! neither electrical/shaft power nor a time-integrated energy prediction.

use super::*;
use fs_airflow::graph::sensitivity::{FlowWeight, GraphSensitivityConfig, PressureWeight};
use fs_solver::CheckedCgConfig;

pub(super) fn parse(value: Option<&J>, request: &Request) -> Result<f64> {
    let Some(value) = value else { return Ok(0.0); };
    object(value, &["weight_per_w"], "air_power_cost")?;
    if request.fan.is_none() { return Err(bad("air_power_cost requires hydraulics.fan")); }
    positive(get(value, "weight_per_w")?, "air_power_cost.weight_per_w")
}

pub(super) fn value(request: &Request, flow: &GraphSolution) -> Result<f64> {
    let fan = request.fan.as_ref().ok_or_else(|| bad("air-power cost lost its fan drive"))?;
    let q = flow.node_outflows[fan.inlet].value();
    let pressure = finite(flow.pressures[fan.inlet].value() - flow.pressures[fan.outlet].value())?;
    let power = finite(q * pressure)?;
    if q <= 0.0 || pressure <= 0.0 || power <= 0.0 {
        return Err(producer("air-power cost needs representable positive fan flow and pressure"));
    }
    Ok(power)
}

pub(super) fn add_cost(
    request: &Request, flow: &GraphSolution, design: &Specification,
    cx: &Cx<'_>, sample: &mut ProjectedAlSample,
) -> Result<()> {
    if design.air_power_weight == 0.0 { return Ok(()); }
    poll(cx)?;
    let fan = request.fan.as_ref().ok_or_else(|| bad("air-power cost lost its fan drive"))?;
    let bank = fan.bank(fan.speed_ratio)?;
    let primal = GraphSolveConfig { max_sweeps: request.limits.graph, max_node_iterations: 80,
        absolute_flow_tolerance: VolumetricFlowRate::new(request.limits.flow), relative_flow_tolerance: 0.0 };
    let sensitivity = GraphSensitivityConfig {
        minimum_absolute_flow: VolumetricFlowRate::new(request.limits.flow),
        linear: CheckedCgConfig { tolerance: request.limits.relative,
            max_iterations: request.limits.linear, max_corrections: 2 },
    };
    let state = request.graph.solve_with_fan_linearized(
        fan.inlet, fan.outlet, &bank, primal, fan.config, sensitivity, cx,
    ).map_err(producer)?;
    if state.solution() != flow { return Err(producer("power and thermal hydraulic states differ")); }
    let power = value(request, flow)?;
    let q = flow.node_outflows[fan.inlet].value();
    let pressure = finite(flow.pressures[fan.inlet].value() - flow.pressures[fan.outlet].value())?;
    // Seed P/(1 W) in the dimensionless hydraulic API; numeric derivatives
    // below restore watts. Both product-rule terms are necessary.
    let mut wp = vec![PressureWeight::new(0.0); flow.pressures.len()];
    wp[fan.inlet] = PressureWeight::new(q);
    wp[fan.outlet] = PressureWeight::new(-q);
    let mut wq = Vec::with_capacity(flow.branches.len());
    for branch in &flow.branches {
        poll(cx)?;
        let sign = if branch.from == fan.inlet { 1.0 }
            else if branch.to == fan.inlet { -1.0 } else { 0.0 };
        wq.push(FlowWeight::new(sign * pressure));
    }
    let gradient = state.adjoint(&wp, &wq, cx).map_err(producer)?;
    sample.objective = finite(sample.objective + finite(design.air_power_weight * power)?)?;
    for (index, control) in design.controls.iter().enumerate() {
        poll(cx)?;
        sample.gradient[index] = finite(sample.gradient[index]
            + finite(design.air_power_weight * gradient.log_resistances[control.branch])?)?;
    }
    if design.fan_control.is_some() {
        let index = design.controls.len();
        sample.gradient[index] = finite(sample.gradient[index]
            + finite(design.air_power_weight * gradient.log_speed_ratio)?)?;
    }
    poll(cx)?;
    Ok(())
}

pub(super) fn render(request: &Request, flow: &GraphSolution, weight: f64) -> Result<String> {
    if weight == 0.0 { return Ok("null".into()); }
    let power = value(request, flow)?;
    Ok(format!("{{\"weight_per_w\":{},\"selected_air_power_w\":{},\"selected_cost_contribution\":{},\"electrical_power_w\":null,\"scope\":\"nominal volume flow times network pressure rise; no motor efficiency, fan heating, shaft-power or time-integrated energy model\"}}",
        num(weight)?, num(power)?, num(finite(weight * power)?)?))
}
