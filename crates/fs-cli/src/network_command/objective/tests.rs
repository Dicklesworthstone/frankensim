use super::*;
use super::super::tests::{close, request, with_cx};

const HOTSPOT: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/size-heterogeneous-hotspot.json"));
fn maximum(r: &Request, selector: &str) -> Objective {
    Objective::parse(&J::parse(selector).unwrap(), &r.surfaces, &r.mesh).unwrap()
}

#[test]
fn exact_maximum_selects_lowest_id_at_ties_and_never_averages_them() {
    let r = request();
    let objective = maximum(&r, r#"{"max_vertices":[3,1,2],"gradient":true}"#);
    let mut t = vec![300.0; r.mesh.vertex_count()];
    t[1] = 350.0; t[3] = 350.0; t[2] = 349.999;
    with_cx(|cx| {
        let s = objective.evaluate(cx, &t, &[]).unwrap();
        assert_eq!(s.value, 350.0); assert_eq!(s.vertex, Some(1));
        assert_eq!(s.ties, 2); assert_eq!(s.separation_k, Some(0.0));
        t[3] += 1e-6;
        let s = objective.evaluate(cx, &t, &[]).unwrap();
        assert_eq!(s.vertex, Some(3)); assert_eq!(s.ties, 1);
        close(s.separation_k.unwrap(), 1e-6, 1e-12);
    });
}

#[test]
fn maximum_scope_is_explicit_and_invalid_selectors_refuse() {
    let r = request();
    with_cx(|cx| {
        let mut t = vec![300.0; r.mesh.vertex_count()]; t[4] = 400.0;
        let face = maximum(&r, r#"{"max_wall_region":"first-face"}"#).evaluate(cx, &t, &[]).unwrap();
        let all = maximum(&r, r#"{"max_solid_temperature":true}"#).evaluate(cx, &t, &[]).unwrap();
        assert_eq!(face.value, 300.0); assert_eq!(all.value, 400.0);
        let selected = maximum(&r, r#"{"max_vertices":[4]}"#);
        t[4] = f64::NAN;
        assert!(selected.evaluate(cx, &t, &[]).is_err());
    });
    for json in [r#"{}"#, r#"{"max_vertices":[]}"#, r#"{"max_vertices":[0,0]}"#,
        r#"{"max_solid_temperature":false}"#, r#"{"max_vertices":[99999]}"#,
        r#"{"max_wall_region":"missing"}"#,
        r#"{"mean_wall_region":"first-face","max_solid_temperature":true}"#] {
        assert!(Objective::parse(&J::parse(json).unwrap(), &r.surfaces, &r.mesh).is_err());
    }
}

#[test]
fn peak_coupled_adjoint_tracks_the_actual_hotspot_not_a_surface_mean() {
    let r = Request::parse(HOTSPOT).unwrap();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let mut h: BTreeMap<_, _> = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let nominal = r.evaluate(cx, &flow, &h, true).unwrap();
        assert_eq!(nominal.objective_state.vertex, Some(4));
        assert_eq!(nominal.objective_state.ties, 1);
        assert!(nominal.objective > 301.5);
        assert!(nominal.coupled.solid[0].mean_wall_temperature_k < 301.5);
        let delta = 1e-4_f64;
        h.insert("last-face".into(), 80.0 * delta.exp());
        let plus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        h.insert("last-face".into(), 80.0 * (-delta).exp());
        let minus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        close(nominal.gradient.as_ref().unwrap().log_htc[1], (plus - minus) / (2.0 * delta), 3e-5);
        let doc = J::parse(&render(&r, &flow, &nominal).unwrap()).unwrap();
        assert_eq!(doc.get("objective_mean_k"), Some(&J::Null));
        assert_eq!(doc.get("dmean_dinlet_k"), Some(&J::Null));
        assert_eq!(doc.get("objective").unwrap().f64_field("active_vertex"), Some(4.0));
        assert!(doc.get("dobjective_dinlet_k").unwrap().as_array().is_some());
    });
}

#[test]
fn hotspot_sizing_rechecks_every_vertex_even_when_the_active_vertex_changes() {
    let r = Request::parse(HOTSPOT).unwrap();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let mut h: BTreeMap<_, _> = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        h.insert("last-face".into(), 10.0);
        assert_eq!(r.evaluate(cx, &flow, &h, false).unwrap().objective_state.vertex, Some(5));
        let d = design::solve(&r, cx, &flow, r.design.as_ref().unwrap()).unwrap();
        assert!(d.passing.temperatures.iter().all(|&t| t <= 301.5));
        assert!(301.5 - d.passing.objective <= 1e-5);
        assert_eq!(d.passing.objective_state.vertex, Some(4));
        let doc = J::parse(&design::attach(render(&r, &flow, &d.passing).unwrap(), &d).unwrap()).unwrap();
        let lower = doc.path(&["design", "failed_lower"]).unwrap();
        assert!(lower.f64_field("temperature_k").unwrap() > 301.5);
        assert_eq!(lower.get("mean_temperature_k"), Some(&J::Null));
    });
}

#[test]
fn a_mean_limit_cannot_be_reinterpreted_as_a_peak_limit() {
    assert!(Request::parse(&HOTSPOT.replace("\"temperature_limit_k\"", "\"mean_temperature_limit_k\"")).is_err());
    assert!(Request::parse(&HOTSPOT.replace("\"temperature_limit_k\": 301.5", "\"temperature_limit_k\": 301.5,\"mean_temperature_limit_k\": 301.5")).is_err());
}
