use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-hotspot.json"));

fn request(fan: bool) -> Request {
    // Positive final pressure allows two-sided perturbations of both samples.
    let text = FIXTURE.replace("\"static_pressure_pa\": 0", "\"static_pressure_pa\": 4");
    let mut request = Request::parse(&text).unwrap();
    request.gradient = true;
    request.objective = objective::Objective::MeanWall(request.surfaces[0].name.clone());
    if !fan {
        request.fan = None;
        request.boundaries = vec![FixedPressure { node: 2, pressure: Pressure::new(0.0) },
            FixedPressure { node: 0, pressure: Pressure::new(20.0) }];
    }
    request
}
fn run(request: &Request) -> J {
    J::parse(&execute(request, &CancelGate::new_clock_free()).unwrap()).unwrap()
}
fn value(request: &Request) -> f64 {
    run(request).get("objective").unwrap().f64_field("value_k").unwrap()
}
fn values(value: &J, key: &str) -> Vec<f64> {
    value.get(key).unwrap().as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect()
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 4e-4 * b.abs().max(1.0), "{a:.12e} != {b:.12e}");
}

#[test]
fn cli_resistance_gradients_match_complete_thermal_resolves_for_both_drives() {
    for fan in [false, true] {
        let base = request(fan);
        let output = run(&base);
        let sensitivity = output.get("hydraulic_sensitivity").unwrap();
        assert_eq!(sensitivity.str_field("status"), Some("available"));
        assert_eq!(sensitivity.get("branch_order").unwrap().as_array().unwrap().len(), 3);
        let gradient = values(sensitivity, "dobjective_dlog_resistance_k");
        let h = 1e-3;
        for (edge, &expected) in gradient.iter().enumerate() {
            let perturb = |step: f64| {
                let mut request = request(fan);
                request.gradient = false;
                let mut branches = request.graph.branches().to_vec();
                branches[edge].loss.resistance = LossResistance::new(
                    branches[edge].loss.resistance.value() * (1.0 + step));
                request.graph = LossGraph::new(request.graph.node_count(), branches).unwrap();
                value(&request)
            };
            close(expected, (perturb(h)-perturb(-h))/(2.0*h));
        }
    }
}

#[test]
fn cli_boundary_order_and_pressure_gauge_derivatives_are_preserved() {
    let output = run(&request(false));
    let sensitivity = output.get("hydraulic_sensitivity").unwrap();
    assert_eq!(values(sensitivity, "boundary_nodes"), vec![2.0, 0.0]);
    let gradient = values(sensitivity, "dobjective_dboundary_pressure_k_pa");
    close(gradient.iter().sum(), 0.0);
    for (boundary, &expected) in gradient.iter().enumerate() {
        let h = 1e-2;
        let perturb = |step: f64| {
            let mut request = request(false); request.gradient = false;
            request.boundaries[boundary].pressure = Pressure::new(request.boundaries[boundary].pressure.value()+step);
            value(&request)
        };
        close(expected, (perturb(h)-perturb(-h))/(2.0*h));
    }
    assert_eq!(sensitivity.get("dobjective_dlog_fan_speed_ratio_k"), Some(&J::Null));
}

#[test]
fn cli_fan_speed_and_curve_controls_include_moving_operating_point() {
    let output = run(&request(true));
    let sensitivity = output.get("hydraulic_sensitivity").unwrap();
    let expected = sensitivity.f64_field("dobjective_dlog_fan_speed_ratio_k").unwrap();
    let h = 1e-3;
    let perturb = |step: f64| {
        let mut request = request(true); request.gradient = false;
        request.fan.as_mut().unwrap().speed_ratio *= 1.0 + step;
        value(&request)
    };
    close(expected, (perturb(h)-perturb(-h))/(2.0*h));
    close(expected, output.get("fan_speed_sensitivity").unwrap().f64_field("dobjective_dlog_speed_ratio_k").unwrap());
    let curve = values(sensitivity, "dobjective_dfan_curve_pressure_k_pa");
    for (index, &gradient) in curve.iter().enumerate() {
        let run = |step: f64| {
            let (old, nominal) = if index == 0 { ("\"static_pressure_pa\": 40",40.0) }
                else { ("\"static_pressure_pa\": 0",4.0) };
            let text = FIXTURE.replace(old,&format!("\"static_pressure_pa\": {}",nominal+step));
            let text = if index == 0 { text.replace("\"static_pressure_pa\": 0","\"static_pressure_pa\": 4") } else { text };
            let mut request = Request::parse(&text).unwrap();
            request.objective = objective::Objective::MeanWall(request.surfaces[0].name.clone());
            value(&request)
        };
        close(gradient, (run(h)-run(-h))/(2.0*h));
    }
}

#[test]
fn unsupported_fan_boundary_keeps_conditional_gradients_without_fake_zeros() {
    let mut request = request(true);
    request.fan.as_mut().unwrap().speed_ratio = 2.0;
    let output = run(&request);
    assert_eq!(output.get("hydraulic_sensitivity").unwrap().str_field("status"),Some("unavailable"));
    assert_eq!(output.get("fan_speed_sensitivity").unwrap().str_field("status"),Some("available"));
    let gate = CancelGate::new_clock_free(); gate.request();
    assert!(execute(&request,&gate).is_err());
}
