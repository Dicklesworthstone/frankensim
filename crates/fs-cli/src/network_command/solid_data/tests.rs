use super::*;
use super::super::tests::{close, request, with_cx};

fn assign(r: &mut Request, text: &str) {
    let (data, k, source) = SolidData::parse(&J::parse(text).unwrap(), &r.mesh).unwrap();
    r.solid_data = data; r.conductivity = k; r.source = source;
}
fn layered(r: &Request) -> String {
    let names = r.mesh.complex().tets.iter().map(|tet| {
        let x = tet.iter().map(|&v| r.mesh.positions()[v as usize][0]).sum::<f64>() / 4.0;
        quote(if x < 0.025 { "spreader" } else { "substrate" })
    }).collect::<Vec<_>>().join(",");
    format!(r#"{{"materials":[{{"name":"spreader","conductivity_w_m_k":20,"source":"fixture"}},{{"name":"substrate","conductivity_w_m_k":2,"source":"fixture"}}],"element_materials":[{names}],"source_w_m3":0}}"#)
}

#[test]
fn heterogeneous_solid_matches_two_layer_continuous_series_solution() {
    let mut r = request();
    let text = layered(&r); assign(&mut r, &text);
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let e = r.evaluate(cx, &flow, &h, false).unwrap();
        let c1 = 1.2 * 0.003 * 1007.0_f64;
        let ct = 1.2 * 0.004 * 1007.0_f64;
        let film1 = c1 * -(-0.5 / c1).exp_m1();
        let film2 = ct * -(-0.8 / ct).exp_m1();
        let heat = 10.0 / (0.025 / (20.0 * 0.01) + 0.025 / (2.0 * 0.01)
            + 1.0 / film1 + 1.0 / film2 - 1.0 / ct);
        let first = 330.0 - heat / film1;
        for (p, &t) in r.mesh.positions().iter().zip(&e.temperatures) {
            let resistance = p[0].min(0.025) / (20.0 * 0.01)
                + (p[0] - 0.025).max(0.0) / (2.0 * 0.01);
            close(t, first - heat * resistance, 1e-5);
        }
        close(e.coupled.solid[0].heat_rate_w, -heat, 1e-5);
        close(e.coupled.solid[1].heat_rate_w, heat, 1e-5);
        let doc = J::parse(&render(&r, &flow, &e).unwrap()).unwrap();
        assert_eq!(doc.path(&["solid_inputs", "constitutive"]).unwrap().str_field("mode"), Some("element-materials"));
    });
}

#[test]
fn overlapping_components_deliver_their_declared_watts_to_the_air() {
    let mut r = request();
    assign(&mut r, r#"{"conductivity_w_m_k":10,"component_power":{"total_w":1,"relative_tolerance":1e-12,"components":[{"name":"chip","watts":0.75,"vertices":[0]},{"name":"aux","watts":0.25,"vertices":[0]}]}}"#);
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let e = r.evaluate(cx, &flow, &h, false).unwrap();
        close(e.source_total_w, 1.0, 1e-10);
        close(e.robin_total_w, 1.0, r.limits.heat);
        close(e.coupled.transport.external_heat_gain_w, 1.0, 2.0 * r.limits.heat);
        let doc = J::parse(&render(&r, &flow, &e).unwrap()).unwrap();
        let audit = doc.path(&["solid_inputs", "heating"]).unwrap();
        close(audit.f64_field("delivered_total_w").unwrap(), 1.0, 1e-12);
        assert_eq!(audit.get("uncertainty_w"), Some(&J::Null));
        assert_eq!(audit.get("components").unwrap().as_array().unwrap().len(), 2);
    });
}

#[test]
fn full_footprint_component_reproduces_the_uniform_source_field() {
    let mut r = request(); r.source = 2000.0;
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let a = r.evaluate(cx, &flow, &h, false).unwrap();
        let vertices = (0..r.mesh.vertex_count()).map(|v| v.to_string()).collect::<Vec<_>>().join(",");
        assign(&mut r, &format!(r#"{{"conductivity_w_m_k":10,"component_power":{{"total_w":1,"relative_tolerance":1e-12,"components":[{{"name":"all","watts":1,"vertices":[{vertices}]}}]}}}}"#));
        let b = r.evaluate(cx, &flow, &h, false).unwrap();
        for (&a, &b) in a.temperatures.iter().zip(&b.temperatures) { close(a, b, 1e-6); }
    });
}

#[test]
fn heterogeneous_powered_adjoint_matches_perturbed_coupled_solve() {
    let mut r = request(); let text = layered(&r); assign(&mut r, &text);
    r.source = 1000.0;
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let mut h: BTreeMap<_, _> = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let a = r.evaluate(cx, &flow, &h, true).unwrap();
        let delta = 1e-4_f64;
        h.insert("last-face".into(), 80.0 * delta.exp());
        let plus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        h.insert("last-face".into(), 80.0 * (-delta).exp());
        let minus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        close(a.gradient.unwrap().log_htc[1], (plus - minus) / (2.0 * delta), 3e-5);
    });
}

#[test]
fn incomplete_ambiguous_and_unbalanced_solid_inputs_refuse() {
    let r = request();
    let valid = r#"{"conductivity_w_m_k":10,"component_power":{"total_w":1,"relative_tolerance":1e-12,"components":[{"name":"chip","watts":1,"vertices":[0]}]}}"#;
    for text in [
        "{}".to_string(),
        valid.replace("\"total_w\":1", "\"total_w\":2"),
        valid.replace("\"vertices\":[0]", "\"vertices\":[0,0]"),
        valid.replace("\"vertices\":[0]", "\"vertices\":[20000]"),
        valid.replace("\"watts\":1", "\"watts\":-1"),
        valid.replace("\"component_power\":", "\"source_w_m3\":0,\"component_power\":"),
        layered(&r).replace("\"source_w_m3\":0", "\"conductivity_w_m_k\":10,\"source_w_m3\":0"),
        layered(&r).replace("\"spreader\",", "\"missing\","),
        r#"{"materials":[],"element_materials":[],"source_w_m3":0}"#.to_string(),
    ] {
        assert!(SolidData::parse(&J::parse(&text).unwrap(), &r.mesh).is_err(), "accepted {text}");
    }
}
