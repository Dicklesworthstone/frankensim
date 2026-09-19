use super::*;
use super::super::tests::{request as base, close, with_cx};

fn configured() -> Request {
    Request::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/recirculated-slab.json"))).unwrap()
}

fn evaluate(request: &Request) -> Evaluation {
    with_cx(|cx| {
        let flow = request.flow(cx).unwrap();
        let htc = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        request.evaluate(cx, &flow, &htc, false).unwrap()
    })
}

#[test]
fn accepted_transient_endpoints_close_fresh_exhaust_energy() {
    let mut request = configured();
    request.gradient = false;
    request.source = 2000.0;
    let value = J::parse(r#"{"initial_temperature_k":300,
        "volumetric_heat_capacity_j_m3_k":1000,"max_step_s":0.05,"max_steps":8,
        "intervals":[{"duration_s":0.1,"power_scale":1},
                     {"duration_s":0.1,"power_scale":0}]}"#).unwrap();
    request.transient = Some(transient::Schedule::parse(&value,
        request.mesh.vertex_count(), request.mesh.element_count(), None).unwrap());
    let result = execute(&request, &CancelGate::new_clock_free()).unwrap();
    let root = J::parse(&result).unwrap();
    let trajectory = root.get("transient").unwrap();
    let history = trajectory.get("history").unwrap().as_array().unwrap();
    assert_eq!(history.len(), 5);
    let mut external = 0.0;
    let mut previous_mixed = None;
    let mut intake_changed = false;
    for sample in &history[1..] {
        let report = sample.get("recirculation").unwrap();
        let dt = sample.f64_field("dt_s").unwrap();
        external += dt * report.f64_field("external_heat_gain_w").unwrap();
        let mixed = report.get("mixed_supplies").unwrap().as_array().unwrap()[0]
            .f64_field("mixed_temperature_k").unwrap();
        if previous_mixed.is_some_and(|old: f64| (old - mixed).abs() > 1e-5) {
            intake_changed = true;
        }
        previous_mixed = Some(mixed);
        assert!(report.f64_field("max_mixing_residual_k").unwrap() <= 1e-9);
    }
    assert!(intake_changed, "return intake must follow evolving exhaust, not a frozen initial sample");
    close(trajectory.f64_field("fresh_exhaust_energy_gain_j").unwrap(), external, 1e-10);
    let residual = trajectory.f64_field("stored_energy_change_j").unwrap()
        - trajectory.f64_field("input_energy_j").unwrap() + external;
    close(trajectory.f64_field("fresh_exhaust_energy_residual_j").unwrap(), residual, 1e-10);
    assert!(residual.abs() <= request.limits.heat * 0.2);
}

#[test]
fn endpoint_history_refuses_a_lost_return_model() {
    let request = configured();
    let once = evaluate(&base());
    assert!(history_field(&request, &once.coupled.transport).is_err());
    let actual = evaluate(&request);
    let fragment = history_field(&request, &actual.coupled.transport).unwrap();
    J::parse(&format!("{{\"endpoint\":true{fragment}}}")).unwrap();
}
