//! Fan-bank boundary drive for the existing passive hydraulic graph.
//! The fan connects ambient to one supply node; it is not an internal edge,
//! recirculation loop, prescribed pressure, or electrical-power model.

use super::*;
use fs_airflow::{FanArrangement, FanBank, FanCurve, FanPoint};
use fs_airflow::graph::GraphFanConfig;

#[derive(Debug)]
pub(super) struct FanDrive {
    curve: FanCurve,
    count: usize,
    arrangement: FanArrangement,
    pub speed_ratio: f64,
    pub inlet: usize,
    outlet: usize,
    temperature_k: f64,
    config: GraphFanConfig,
}

impl FanDrive {
    pub(super) fn parse(value: &J, nodes: usize) -> Result<Self> {
        object(value, &["name", "source", "source_id", "inlet", "outlet", "temperature_k",
            "count", "arrangement", "speed_ratio", "min_speed_ratio", "max_speed_ratio",
            "min_flow_m3_s", "pressure_tolerance_rel", "points", "max_iterations",
            "flow_tolerance_m3_s", "pressure_tolerance_pa"], "hydraulics.fan")?;
        let name = string(get(value, "name")?, "fan.name")?;
        let source = string(get(value, "source")?, "fan.source")?;
        let source_id = string(get(value, "source_id")?, "fan.source_id")?;
        let inlet = integer(get(value, "inlet")?, "fan.inlet", nodes - 1)?;
        let outlet = integer(get(value, "outlet")?, "fan.outlet", nodes - 1)?;
        if inlet == outlet { return Err(bad("fan terminals must be distinct")); }
        let arrangement = match get(value, "arrangement")?.as_str() {
            Some("series") => FanArrangement::Series,
            Some("parallel") => FanArrangement::Parallel,
            _ => return Err(bad("fan arrangement must be series or parallel")),
        };
        let mut points = Vec::new();
        for point in array(get(value, "points")?, "fan.points", 4096)? {
            object(point, &["flow_m3_s", "static_pressure_pa"], "fan.point")?;
            points.push(FanPoint::new(
                VolumetricFlowRate::new(number(get(point, "flow_m3_s")?, "fan.flow_m3_s")?),
                Pressure::new(number(get(point, "static_pressure_pa")?, "fan.static_pressure_pa")?),
            ));
        }
        let curve = FanCurve::new(name, points, SourceProvenance::new(source, source_id),
            number(get(value, "pressure_tolerance_rel")?, "fan.pressure_tolerance_rel")?,
            ToleranceBasis::EngineeringAllowance,
            VolumetricFlowRate::new(number(get(value, "min_flow_m3_s")?, "fan.min_flow_m3_s")?),
            (positive(get(value, "min_speed_ratio")?, "fan.min_speed_ratio")?,
                positive(get(value, "max_speed_ratio")?, "fan.max_speed_ratio")?),
        ).map_err(producer)?;
        let drive = Self { curve, count: count(get(value, "count")?, "fan.count", 1024)?, arrangement,
            speed_ratio: positive(get(value, "speed_ratio")?, "fan.speed_ratio")?, inlet, outlet,
            temperature_k: positive(get(value, "temperature_k")?, "fan.temperature_k")?,
            config: GraphFanConfig {
                max_iterations: count(get(value, "max_iterations")?, "fan.max_iterations", 4096)?,
                absolute_flow_tolerance: VolumetricFlowRate::new(positive(get(value, "flow_tolerance_m3_s")?, "fan.flow_tolerance_m3_s")?),
                absolute_pressure_tolerance: Pressure::new(positive(get(value, "pressure_tolerance_pa")?, "fan.pressure_tolerance_pa")?),
            } };
        drive.bank(drive.speed_ratio)?;
        Ok(drive)
    }

    pub(super) fn bank(&self, speed: f64) -> Result<FanBank> {
        FanBank::new(self.curve.clone(), self.count, self.arrangement, speed).map_err(producer)
    }

    pub(super) fn solve(&self, cx: &Cx<'_>, graph: &LossGraph, limits: Limits, speed: f64) -> Result<GraphSolution> {
        let bank = self.bank(speed)?;
        let operating = graph.solve_with_fan(self.inlet, self.outlet, &bank,
            GraphSolveConfig { max_sweeps: limits.graph, max_node_iterations: 80,
                absolute_flow_tolerance: VolumetricFlowRate::new(limits.flow), relative_flow_tolerance: 0.0 },
            self.config, cx).map_err(producer)?;
        // A zero-flow fan operating point cannot provide the declared thermal supply.
        if operating.flow.value() <= 0.0 { return Err(producer("fan operating point has no thermal through-flow")); }
        Ok(operating.network)
    }

    pub(super) fn attach(&self, output: String, flow: &GraphSolution, speed: f64) -> Result<String> {
        let bank = self.bank(speed)?;
        let q = flow.node_outflows[self.inlet].value();
        let pressure = flow.pressures[self.inlet].value() - flow.pressures[self.outlet].value();
        let fan_pressure = bank.pressure_at(VolumetricFlowRate::new(q)).map_err(producer)?.value();
        let residual = fan_pressure - pressure;
        if !residual.is_finite() || residual.abs() > self.config.absolute_pressure_tolerance.value() {
            return Err(producer("published flow is not the declared fan operating point"));
        }
        let prefix = output.strip_suffix("}\n").ok_or_else(|| bad("internal result framing mismatch"))?;
        Ok(format!("{prefix},\"fan\":{{\"name\":{},\"source\":{},\"source_id\":{},\"count\":{},\"arrangement\":{},\"speed_ratio\":{},\"inlet\":{},\"outlet\":{},\"flow_m3_s\":{},\"static_pressure_pa\":{},\"pressure_residual_pa\":{},\"air_power_w\":{},\"electrical_power_w\":null,\"pressure_tolerance_rel\":{},\"authority\":\"nominal fan/graph intersection; caller-declared curve and allowance, not propagated uncertainty or hardware validation; no fan heating\"}}}}\n",
            quote(self.curve.name()), quote(&self.curve.source().citation), quote(&self.curve.source().identifier),
            self.count, quote(match self.arrangement { FanArrangement::Series => "series", FanArrangement::Parallel => "parallel" }),
            num(speed)?, self.inlet, self.outlet, num(q)?, num(pressure)?, num(residual)?, num(q * pressure)?,
            num(self.curve.pressure_tolerance_rel())?))
    }
}

/// Exactly one drive is required. A fan supplies fresh air at its discharge;
/// outlet temperature is computed, not reused as a recirculating inlet.
pub(super) fn parse(value: &J, nodes: usize) -> Result<(Vec<FixedPressure>, Vec<TransportInlet>, Option<FanDrive>)> {
    match (value.get("boundaries"), value.get("fan")) {
        (None, Some(fan)) => {
            let drive = FanDrive::parse(fan, nodes)?;
            let inlets = vec![TransportInlet { node: drive.inlet, temperature: Temperature::new(drive.temperature_k) }];
            Ok((Vec::new(), inlets, Some(drive)))
        }
        (Some(entries), None) => {
            let mut boundaries = Vec::new();
            let mut inlets = Vec::new();
            let mut fixed = BTreeSet::new();
            for entry in array(entries, "boundaries", nodes)? {
                object(entry, &["node", "pressure_pa", "temperature_k"], "boundary")?;
                let node = integer(get(entry, "node")?, "boundary.node", nodes - 1)?;
                if !fixed.insert(node) { return Err(bad("duplicate hydraulic boundary")); }
                boundaries.push(FixedPressure { node, pressure: Pressure::new(number(get(entry, "pressure_pa")?, "pressure_pa")?) });
                if let Some(temp) = entry.get("temperature_k") {
                    inlets.push(TransportInlet { node, temperature: Temperature::new(positive(temp, "temperature_k")?) });
                }
            }
            Ok((boundaries, inlets, None))
        }
        _ => Err(bad("hydraulics requires exactly one of boundaries or fan")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::{with_cx, close};
    const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-hotspot.json"));

    #[test]
    fn fan_curve_and_loss_graph_set_pressure_and_flow_with_speed_scaling() {
        let request = Request::parse(FIXTURE).unwrap();
        let fan = request.fan.as_ref().unwrap();
        with_cx(|cx| {
            for speed in [0.5, 1.0, 1.5] {
                let flow = fan.solve(cx, &request.graph, request.limits, speed).unwrap();
                close(flow.node_outflows[fan.inlet].value(), 0.004 * speed, 1e-9);
                close(flow.pressures[fan.inlet].value(), 20.0 * speed * speed, 1e-5);
                close(flow.branches[0].flow.value(), 0.003 * speed, 1e-9);
                let output = fan.attach("{}\n".into(), &flow, speed).unwrap();
                let json = J::parse(&output).unwrap();
                assert_eq!(json.get("fan").unwrap().get("electrical_power_w"), Some(&J::Null));
            }
        });
    }

    #[test]
    fn fan_driven_component_heat_reaches_the_exhaust() {
        let request = Request::parse(FIXTURE).unwrap();
        let output = execute(&request, &CancelGate::new_clock_free()).unwrap();
        let json = J::parse(&output).unwrap();
        close(json.f64_field("source_w").unwrap(), 1.0, 1e-8);
        let temperatures = json.get("node_temperatures_k").unwrap().as_array().unwrap();
        close(temperatures[2].as_f64().unwrap(), 300.0 + 1.0 / (1.2 * 0.004 * 1007.0), 1e-5);
        assert!(json.get("fan").is_some());
    }

    #[test]
    fn fan_input_ambiguity_domain_and_budget_refuse() {
        for (from, to) in [
            ("\"fan\": {", "\"boundaries\": [], \"fan\": {"),
            ("\"speed_ratio\": 1", "\"speed_ratio\": 3"),
            ("\"arrangement\": \"series\"", "\"arrangement\": \"guess\""),
            ("\"count\": 1", "\"count\": 0"),
        ] {
            assert!(FIXTURE.contains(from));
            assert!(Request::parse(&FIXTURE.replace(from, to)).is_err());
        }
        let mut request = Request::parse(FIXTURE).unwrap();
        let fan = request.fan.as_mut().unwrap();
        assert!(fan.bank(3.0).is_err());
        fan.config.max_iterations = 1;
        assert!(with_cx(|cx| request.flow(cx)).is_err());
        let gate = CancelGate::new_clock_free(); gate.request();
        assert!(execute(&request, &gate).is_err());
    }
}
