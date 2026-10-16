use super::*;
use super::super::tests::{request as base, close, with_cx};

fn configured(fractions: [f64; 2]) -> Request {
    let mut request = base();
    let value = J::parse(&format!(r#"{{"model":"prescribed-adiabatic-return",
        "source":"declared regression fixture, not hardware data",
        "temperature_tolerance_k":1e-9,"links":[
        {{"supply_node":0,"return_node":3,"fraction":{}}},
        {{"supply_node":1,"return_node":3,"fraction":{}}}]}}"#,
        fractions[0], fractions[1])).unwrap();
    request.recirculation = Some(Policy::parse(&value, &request.inlets).unwrap());
    request
}

fn evaluate(request: &Request, gradient: bool) -> Evaluation {
    with_cx(|cx| {
        let flow = request.flow(cx).unwrap();
        let htc = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        request.evaluate(cx, &flow, &htc, gradient).unwrap()
    })
}

#[test]
fn file_request_reaches_the_real_return_air_producer() {
    let text = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/recirculated-slab.json"));
    let request = Request::parse(text).unwrap();
    let output = execute(&request, &CancelGate::new_clock_free()).unwrap();
    let result = J::parse(&output).unwrap();
    let report = result.get("recirculation").unwrap();
    assert_eq!(report.str_field("status"), Some("solved"));
    assert_eq!(report.get("supplies").unwrap().as_array().unwrap().len(), 2);
    assert_eq!(report.get("streams").unwrap().as_array().unwrap().len(), 2);
    assert!(report.f64_field("heat_imbalance_w").unwrap().abs() <= request.limits.heat);
}

#[test]
fn zero_returns_preserve_the_once_through_field_and_gradient() {
    let once = evaluate(&base(), true);
    let recycled = evaluate(&configured([0.0, 0.0]), true);
    assert_eq!(once.temperatures, recycled.temperatures);
    assert_eq!(once.gradient.unwrap().log_htc, recycled.gradient.unwrap().log_htc);
    assert!(recycled.coupled.transport.recirculation.is_none());
    let text = execute(&configured([0.0, 0.0]), &CancelGate::new_clock_free()).unwrap();
    assert_eq!(J::parse(&text).unwrap().get("recirculation").unwrap().str_field("status"),
        Some("once-through-zero-returns"));
}

#[test]
fn both_recycled_supplies_match_independent_outer_energy_solution() {
    let c = [1.2 * 0.003 * 1007.0, 1.2 * 0.001 * 1007.0];
    for fractions in [[0.4, 0.0], [0.6, 0.2], [0.0, 0.8], [0.9, 0.9]] {
        for source in [0.0, 2000.0] {
            let mut request = configured(fractions);
            request.source = source;
            let result = evaluate(&request, false);
            let fresh = [c[0] * (1.0 - fractions[0]), c[1] * (1.0 - fractions[1])];
            let power = source * 0.05 * 0.1 * 0.1;
            let exhaust = (fresh[0] * 330.0 + fresh[1] * 290.0 + power) / (fresh[0] + fresh[1]);
            close(result.coupled.transport.node_temperatures_k[3].unwrap(), exhaust, 2e-6);
            let report = result.coupled.transport.recirculation.as_ref().unwrap();
            close(report.external_heat_gain_w, power, 3.0 * request.limits.heat);
            for supply in &report.supplies {
                let i = supply.node;
                let temperature = if i == 0 { 330.0 } else { 290.0 };
                close(supply.mixed_temperature_k,
                    (1.0 - fractions[i]) * temperature + fractions[i] * exhaust, 2e-6);
                close(supply.fresh_capacity_w_per_k, fresh[i], 1e-7);
            }
        }
    }
}

#[test]
fn fresh_temperature_and_htc_adjoints_include_return_feedback() {
    let mut request = configured([0.6, 0.2]);
    let nominal = evaluate(&request, true);
    let gradient = nominal.gradient.unwrap();
    for node in 0..2 {
        let original = request.inlets[node].temperature.value();
        let delta = 1e-3;
        request.inlets[node].temperature = Temperature::new(original + delta);
        let plus = evaluate(&request, false).objective;
        request.inlets[node].temperature = Temperature::new(original - delta);
        let minus = evaluate(&request, false).objective;
        request.inlets[node].temperature = Temperature::new(original);
        close(gradient.inlets[node], (plus - minus) / (2.0 * delta), 3e-5);
    }
    with_cx(|cx| {
        let flow = request.flow(cx).unwrap();
        let mut h: BTreeMap<_, _> = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let delta = 1e-4_f64;
        h.insert("last-face".into(), 80.0 * delta.exp());
        let plus = request.evaluate(cx, &flow, &h, false).unwrap().objective;
        h.insert("last-face".into(), 80.0 * (-delta).exp());
        let minus = request.evaluate(cx, &flow, &h, false).unwrap().objective;
        close(gradient.log_htc[1], (plus - minus) / (2.0 * delta), 3e-5);
    });
}

#[test]
fn malformed_or_physically_unowned_returns_refuse() {
    let original = configured([0.4, 0.2]);
    let inlet = &original.inlets;
    for links in [
        r#"[{"supply_node":0,"return_node":3,"fraction":1}]"#,
        r#"[{"supply_node":0,"return_node":3,"fraction":-0.1}]"#,
        r#"[{"supply_node":0,"return_node":0,"fraction":0.4}]"#,
        r#"[{"supply_node":2,"return_node":3,"fraction":0.4}]"#,
        r#"[{"supply_node":0,"return_node":3,"fraction":0.2},{"supply_node":0,"return_node":3,"fraction":0.3}]"#,
        r#"[{"supply_node":0,"return_node":3,"fraction":0.6},{"supply_node":0,"return_node":2,"fraction":0.4}]"#,
        r#"[{"supply_node":0.0,"return_node":3,"fraction":0.4}]"#,
        r#"[{"supply_node":0,"return_node":3,"fraction":0.4,"ignored":true}]"#,
        "[]",
    ] {
        let value = J::parse(&format!(r#"{{"model":"prescribed-adiabatic-return","source":"test",
            "temperature_tolerance_k":1e-9,"links":{links}}}"#)).unwrap();
        assert!(Policy::parse(&value, inlet).is_err(), "accepted {links}");
    }
    for bad_return in [1, 2, 999] {
        let mut request = configured([0.4, 0.2]);
        request.recirculation.as_mut().unwrap().links[0].return_node = bad_return;
        assert!(execute(&request, &CancelGate::new_clock_free()).is_err());
    }
    let once = evaluate(&base(), false);
    assert!(original.recirculation.as_ref().unwrap().report(&once.coupled.transport).is_err());
}

#[test]
fn return_link_order_and_branch_orientation_do_not_change_the_solution() {
    let request = configured([0.6, 0.2]);
    let expected = evaluate(&request, false);
    let mut reversed = configured([0.6, 0.2]);
    let value = J::parse(r#"{"model":"prescribed-adiabatic-return","source":"test",
        "temperature_tolerance_k":1e-9,"links":[
        {"supply_node":1,"return_node":3,"fraction":0.2},
        {"supply_node":0,"return_node":3,"fraction":0.6}]}"#).unwrap();
    reversed.recirculation = Some(Policy::parse(&value, &reversed.inlets).unwrap());
    let mut branches = reversed.graph.branches().to_vec();
    let edge = &mut branches[0];
    std::mem::swap(&mut edge.from, &mut edge.to);
    reversed.region_paths[0].reverse();
    reversed.graph = LossGraph::new(4, branches).unwrap();
    let actual = evaluate(&reversed, false);
    close(actual.objective, expected.objective, 1e-7);
}

#[test]
fn cancelled_return_air_run_publishes_no_result() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(execute(&configured([0.4, 0.2]), &gate).is_err());
}
