//! Compose accepted thermal, convection and hydraulic adjoints for the CLI.
//! Temperature objectives are normalized by 1 K before passing numeric weights
//! to the dimensionless hydraulic objective API; rendered gradients restore K.

use super::*;
use fs_airflow::graph::sensitivity::{
    FlowWeight, GraphFanLinearization, GraphFanSensitivityError, GraphLinearization,
    GraphSensitivityConfig, PressureWeight,
};
use fs_solver::CheckedCgConfig;

#[derive(Debug)]
pub(super) struct HydraulicGradient {
    branches: Vec<String>,
    flows: Vec<f64>,
    log_resistances: Vec<f64>,
    boundary_nodes: Vec<usize>,
    boundary_pressures: Option<Vec<f64>>,
    log_speed: Option<f64>,
    curve_pressures: Option<Vec<f64>>,
    residual: f64,
    iterations: usize,
}

impl HydraulicGradient {
    pub(super) fn render(&self) -> Result<String> {
        Ok(format!("{{\"status\":\"available\",\"method\":\"coupled-thermal-hydraulic-adjoint\",\"branch_order\":[{}],\"dobjective_dflow_k_s_m3\":{},\"dobjective_dlog_resistance_k\":{},\"boundary_nodes\":[{}],\"dobjective_dboundary_pressure_k_pa\":{},\"dobjective_dlog_fan_speed_ratio_k\":{},\"dobjective_dfan_curve_pressure_k_pa\":{},\"hydraulic_linear_relative_residual\":{},\"hydraulic_linear_iterations\":{},\"scope\":\"local selected temperature objective; moving hydraulic equilibrium, junction mixing, shared-solid feedback and admitted flow-derived convection; fixed geometry, material laws, loads and fluid properties; branch-flow weights are an intermediate chain-rule quantity, not independently feasible controls; no uncertainty, continuum-error or finite-design-change bound\"}}",
            self.branches.iter().map(|name| quote(name)).collect::<Vec<_>>().join(","),
            numbers(&self.flows)?, numbers(&self.log_resistances)?,
            self.boundary_nodes.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
            self.boundary_pressures.as_ref().map(|v| numbers(v)).transpose()?.unwrap_or_else(|| "null".into()),
            optional(self.log_speed)?,
            self.curve_pressures.as_ref().map(|v| numbers(v)).transpose()?.unwrap_or_else(|| "null".into()),
            num(self.residual)?, self.iterations))
    }
}

enum Drive<'g> {
    Pressure(GraphLinearization<'g>),
    Fan(GraphFanLinearization<'g>),
}
impl Drive<'_> {
    fn solution(&self) -> &GraphSolution {
        match self { Self::Pressure(state) => state.solution(), Self::Fan(state) => state.solution() }
    }
}

/// None means this extra control family is unsupported, not a zero gradient.
/// Numerical failure or cancellation after admission still refuses publication.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_pullback(
    request: &Request, cx: &Cx<'_>, binding: &CoupledLinearization<'_, '_>,
    objective: &CoupledObjective, names: &[&str], derivations: &[convection::Derived],
) -> Result<Option<CoolingGradient>> {
    poll(cx)?;
    // Fan-speed searches evaluate speeds other than Request's nominal setting;
    // their existing affinity path owns that trajectory. Radiation/return and
    // transient adjoints likewise retain their own complete feedback model.
    if request.recirculation.is_some() || request.fan_speed_design.is_some()
        || request.radiation.is_some() || request.transient.is_some()
    { return Ok(None); }
    let flow = binding.hydraulics();
    for branch in &flow.branches {
        poll(cx)?;
        if branch.flow.value().abs() <= request.limits.flow { return Ok(None); }
    }
    for boundary in &flow.boundaries {
        poll(cx)?;
        if flow.node_outflows[boundary.node].value().abs() <= request.limits.flow { return Ok(None); }
    }
    let mut slopes = BTreeMap::new();
    for derived in derivations {
        poll(cx)?;
        let Some(slope) = reynolds_elasticity(&derived.nu)? else { return Ok(None); };
        slopes.insert(derived.surface.as_str(), slope);
    }
    let primal = GraphSolveConfig { max_sweeps: request.limits.graph, max_node_iterations: 80,
        absolute_flow_tolerance: VolumetricFlowRate::new(request.limits.flow), relative_flow_tolerance: 0.0 };
    let config = GraphSensitivityConfig {
        minimum_absolute_flow: VolumetricFlowRate::new(request.limits.flow),
        linear: CheckedCgConfig { tolerance: request.limits.relative,
            max_iterations: request.limits.linear, max_corrections: 2 },
    };
    let drive = if let Some(fan) = &request.fan {
        let bank = fan.bank(fan.speed_ratio)?;
        match request.graph.solve_with_fan_linearized(fan.inlet, fan.outlet, &bank, primal, fan.config, config, cx) {
            Ok(state) => Drive::Fan(state),
            Err(GraphFanSensitivityError::NonDifferentiableFan(_)) => return Ok(None),
            Err(error) => return Err(producer(error)),
        }
    } else {
        Drive::Pressure(request.graph.solve_linearized(&request.boundaries, primal, config, cx).map_err(producer)?)
    };
    // Reusing a stale or differently driven hydraulic state would produce a
    // plausible but incorrect chain rule. Both solves use identical declared
    // budgets and deterministic arithmetic; require the exact original state.
    if drive.solution() != flow { return Err(producer("hydraulic and thermal derivative states differ")); }
    let interface = InterfaceSolveConfig {
        max_iterations: request.limits.derivative, absolute_tolerance: request.limits.relative,
        relative_tolerance: request.limits.relative, relaxation: request.limits.relaxation,
    };
    let (thermal, mut flows) = binding.pullback_branch_flows_iqn(
        cx, objective, interface, acceleration::POLICY,
    ).map_err(producer)?;
    if names.len() != thermal.log_htc.len() || flows.len() != flow.branches.len() {
        return Err(bad("hydraulic gradient control ordering mismatch"));
    }
    let mut capacity = 0.0;
    for (gradient, branch) in flows.iter().zip(&flow.branches) {
        poll(cx)?;
        capacity = checked(capacity + checked(gradient * branch.flow.value())?)?;
    }
    let mut owners = BTreeMap::new();
    for (edge, regions) in request.region_paths.iter().enumerate() {
        for region in regions { poll(cx)?; owners.insert(region.as_str(), edge); }
    }
    let mut convection = 0.0;
    for (name, gradient) in names.iter().zip(&thermal.log_htc) {
        poll(cx)?;
        let edge = *owners.get(name).ok_or_else(|| bad("missing thermal region hydraulic owner"))?;
        let contribution = checked(gradient * slopes.get(name).copied().unwrap_or(0.0))?;
        // ln(|q|) has derivative 1/q, including reversed declared edges.
        flows[edge] = checked(flows[edge] + checked(contribution / flow.branches[edge].flow.value())?)?;
        convection = checked(convection + contribution)?;
    }
    let mut output = HydraulicGradient {
        branches: flow.branches.iter().map(|b| b.loss.name.clone()).collect(),
        flows, log_resistances: Vec::new(), boundary_nodes: Vec::new(), boundary_pressures: None,
        log_speed: None, curve_pressures: None, residual: 0.0, iterations: 0,
    };
    let wp = vec![PressureWeight::new(0.0); flow.pressures.len()];
    let wq: Vec<_> = output.flows.iter().copied().map(FlowWeight::new).collect();
    match &drive {
        Drive::Pressure(state) => {
            let gradient = state.adjoint(&wp, &wq, cx).map_err(producer)?;
            output.log_resistances = gradient.log_resistances;
            output.boundary_nodes = flow.boundaries.iter().map(|b| b.node).collect();
            output.boundary_pressures = Some(gradient.boundary_pressures.iter().map(|p| p.value()).collect());
            output.residual = gradient.linear_report.rel_residual;
            output.iterations = gradient.linear_report.iters;
        }
        Drive::Fan(state) => {
            let gradient = state.adjoint(&wp, &wq, cx).map_err(producer)?;
            output.log_resistances = gradient.log_resistances;
            output.log_speed = Some(gradient.log_speed_ratio);
            output.curve_pressures = Some(gradient.curve_pressures.iter().map(|p| p.value()).collect());
            output.residual = gradient.objective_linear_report.rel_residual.max(gradient.inlet_linear_report.rel_residual);
            output.iterations = gradient.objective_linear_report.iters + gradient.inlet_linear_report.iters;
        }
    }
    let speed = output.log_speed.map(|total| SpeedGradient::Available { flow: capacity, convection, total });
    poll(cx)?;
    Ok(Some(CoolingGradient { thermal, speed, hydraulic: Some(output) }))
}

#[cfg(test)]
mod tests;
