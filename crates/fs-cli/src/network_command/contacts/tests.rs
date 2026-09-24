//! Direct integration checks for declared finite-resistance contacts.
use super::*;
use super::super::tests::{close, with_cx};

// Compacted so fixture edits are independent of the example's formatting.
static CONTACT: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| crate::network_command::json::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/size-contact-slab.json"))));
fn request() -> Request { Request::parse(&CONTACT).unwrap() }
fn coefficients(r: &Request) -> BTreeMap<String, f64> { r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect() }
fn oracle(first: f64, bypass: f64, h: f64) -> (f64, f64, f64) {
    let c1: f64 = 1.2 * 0.003 * 1007.0;
    let ct: f64 = 1.2 * 0.004 * 1007.0;
    let mean = 0.75*first + 0.25*bypass;
    let e1 = c1 * -(-0.5/c1).exp_m1();
    let e2 = ct * -(-h*0.01/ct).exp_m1();
    // Two 25 mm slabs and R''/A = 1 K/W, plus both exponential films
    // and downstream mixed-air feedback. This is not a welded slab.
    let q = (first-mean)/(0.5 + 1.0 + 1.0/e1 + 1.0/e2 - 1.0/ct);
    (first-q/e1, mean-q/ct+q/e2, q)
}

#[test]
fn coupled_contact_resolves_a_temperature_jump_and_the_continuous_slab_solution() {
    for (first, bypass) in [(330.0,290.0), (290.0,330.0), (310.0,310.0)] {
        let mut r = request();
        r.inlets[0].temperature = Temperature::new(first);
        r.inlets[1].temperature = Temperature::new(bypass);
        with_cx(|cx| {
            let flow = r.flow(cx).unwrap();
            let e = r.evaluate(cx, &flow, &coefficients(&r), false).unwrap();
            let (left, right, q) = oracle(first, bypass, 80.0);
            close(e.coupled.solid[0].mean_wall_temperature_k, left, 1e-5);
            close(e.coupled.solid[1].mean_wall_temperature_k, right, 1e-5);
            close(e.contact_fluxes[0].heat_rate_a_to_b_w, q, 1e-5);
            close(e.contact_fluxes[0].mean_jump_k, q, 1e-5);
            for (&a, &b) in [1,4,7,10].iter().zip(&[12,13,14,15]) {
                close(e.temperatures[a] - e.temperatures[b], q, 1e-5);
            }
            close(e.robin_total_w, 0.0, r.limits.heat);
            let output = J::parse(&render(&r, &flow, &e).unwrap()).unwrap();
            let row = &output.get("contacts").unwrap().as_array().unwrap()[0];
            close(row.f64_field("heat_a_to_b_w").unwrap(), q, 1e-5);
            assert_eq!(row.get("uncertainty"), Some(&J::Null));
        });
    }
}

#[test]
fn contact_keeps_heterogeneous_component_power_and_peak_gradients_connected() {
    let text = CONTACT.replace("\"conductivity_w_m_k\":10,", r#""materials":[
        {"name":"a","conductivity_w_m_k":20,"source":"declared"},
        {"name":"b","conductivity_w_m_k":2,"source":"declared"}],
        "element_materials":["a","a","a","a","a","a","b","b","b","b","b","b"],"#)
        .replace("\"source_w_m3\":0,", r#""component_power":{"total_w":1,"relative_tolerance":1e-12,
            "components":[{"name":"chip","watts":1,"vertices":[4]}]},"#);
    let mut r = Request::parse(&text).unwrap();
    for inlet in &mut r.inlets { inlet.temperature = Temperature::new(300.0); }
    r.objective = objective::Objective::parse(&J::parse(r#"{"max_solid_temperature":true}"#).unwrap(), &r.surfaces, &r.mesh).unwrap();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap(); let mut h = coefficients(&r);
        let e = r.evaluate(cx, &flow, &h, true).unwrap();
        close(e.source_total_w, 1.0, 1e-9);
        close(e.robin_total_w, 1.0, r.limits.heat);
        close(e.coupled.transport.external_heat_gain_w, 1.0, 2.0*r.limits.heat);
        assert_eq!(e.objective_state.vertex, Some(4));
        close(e.objective, 301.9936143601, 1e-5);
        let delta = 1e-4_f64;
        h.insert("last-face".into(), 80.0*delta.exp());
        let plus = r.evaluate(cx, &flow, &h, false).unwrap();
        h.insert("last-face".into(), 80.0*(-delta).exp());
        let minus = r.evaluate(cx, &flow, &h, false).unwrap();
        assert_eq!(plus.objective_state.vertex, e.objective_state.vertex);
        assert_eq!(minus.objective_state.vertex, e.objective_state.vertex);
        close(e.gradient.as_ref().unwrap().log_htc[1], (plus.objective-minus.objective)/(2.0*delta), 3e-5);
    });
}

#[test]
fn contact_orientation_reverses_reported_flux_not_the_temperature_field() {
    let r = request();
    let reversed = CONTACT.replace("\"side_a\":", "\"swapped\":")
        .replace("\"side_b\":", "\"side_a\":").replace("\"swapped\":", "\"side_b\":");
    let reversed = Request::parse(&reversed).unwrap();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let a = r.evaluate(cx, &flow, &coefficients(&r), false).unwrap();
        let b = reversed.evaluate(cx, &flow, &coefficients(&reversed), false).unwrap();
        for (&a, &b) in a.temperatures.iter().zip(&b.temperatures) { close(a,b,1e-6); }
        close(a.contact_fluxes[0].heat_rate_a_to_b_w, -b.contact_fluxes[0].heat_rate_a_to_b_w, 1e-6);
        close(a.contact_fluxes[0].mean_jump_k, -b.contact_fluxes[0].mean_jump_k, 1e-6);
    });
}

#[test]
fn missing_nonmatching_and_double_owned_contact_faces_refuse_before_solve() {
    let r = request();
    assert!(Contacts::parse(None, &r.mesh, &r.surfaces, true).is_err());
    for (from, to) in [
        ("\"resistance_m2_k_w\":0.01", "\"resistance_m2_k_w\":0"),
        ("\"side_a\":[1,4,10]", "\"side_a\":[0,3,9]"),
        ("\"side_a\":[1,4,10]", "\"side_a\":[0,1,4]"),
        ("\"side_b\":[12,14,15]", "\"side_b\":[12,13,15]"),
        ("\"adiabatic_remainder\":true", "\"adiabatic_remainder\":false"),
        ("\"contacts\":[", "\"unknown_contacts\":["),
    ] {
        assert!(CONTACT.contains(from));
        assert!(Request::parse(&CONTACT.replace(from,to)).is_err(), "accepted {from} -> {to}");
    }
    assert!(Contacts::parse(Some(&J::Array(Vec::new())), &r.mesh, &r.surfaces, true).is_err());
}

#[test]
fn target_sizing_retains_the_bond_and_returns_its_passing_field() {
    let r = request();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let baseline = r.evaluate(cx, &flow, &coefficients(&r), false).unwrap();
        assert!(baseline.objective > 325.0);
        let selected = design::solve(&r, cx, &flow, r.design.as_ref().unwrap()).unwrap();
        assert!(selected.passing.objective <= 325.0);
        assert!(325.0-selected.passing.objective <= 1e-5);
        let h = selected.passing.htc[1];
        let (left, _, q) = oracle(330.0, 290.0, h);
        close(selected.passing.objective, left, 1e-5);
        close(selected.passing.contact_fluxes[0].heat_rate_a_to_b_w, q, 1e-5);
        close(selected.passing.contact_fluxes[0].mean_jump_k, q, 1e-5);
        close(h, 135.1864, 0.03);
        J::parse(&design::attach(render(&r, &flow, &selected.passing).unwrap(), &selected).unwrap()).unwrap();
    });
}
