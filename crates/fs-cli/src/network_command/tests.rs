use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/mixed-slab.json"));

fn request() -> Request { Request::parse(FIXTURE).expect("explicit tetrahedral request") }
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.16e} versus {expected:.16e}");
}
fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 41, kernel_id: 717, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn oracle(first: f64, bypass: f64, h_last: f64) -> (f64, f64, f64) {
    let c1: f64 = 1.2 * 0.003 * 1007.0;
    let ct: f64 = 1.2 * 0.004 * 1007.0;
    let inlet_mean = 0.75 * first + 0.25 * bypass;
    let e1 = c1 * -(-50.0 * 0.01 / c1).exp_m1();
    let e2 = ct * -(-h_last * 0.01 / ct).exp_m1();
    let heat = (first - inlet_mean) / (0.05 / (10.0 * 0.01) + 1.0 / e1 + 1.0 / e2 - 1.0 / ct);
    (first - heat / e1, inlet_mean - heat / ct + heat / e2, heat)
}

#[test]
fn file_driven_fem_matches_independent_slab_and_mixing_solution() {
    for (first, bypass) in [(330.0, 290.0), (290.0, 330.0), (310.0, 310.0)] {
        let mut r = request();
        r.inlets[0].temperature = Temperature::new(first);
        r.inlets[1].temperature = Temperature::new(bypass);
        r.gradient = false;
        let output = execute(&r, &CancelGate::new_clock_free()).expect("actual coupled solve");
        let json = J::parse(&output).expect("complete JSON result");
        let (t1, t2, heat) = oracle(first, bypass, 80.0);
        close(json.f64_field("objective_mean_k").unwrap(), t1, 1e-5);
        let field = json.get("solid_temperatures_k").unwrap().as_array().unwrap();
        for (position, value) in r.mesh.positions().iter().zip(field) {
            close(value.as_f64().unwrap(), t1 + (t2 - t1) * position[0] / 0.05, 1e-5);
        }
        let walls = json.get("walls").unwrap().as_array().unwrap();
        close(walls[0].f64_field("outward_heat_w").unwrap(), -heat, 1e-5);
        close(walls[1].f64_field("outward_heat_w").unwrap(), heat, 1e-5);
        let nodes = json.get("node_temperatures_k").unwrap().as_array().unwrap();
        close(nodes[3].as_f64().unwrap(), 0.75 * first + 0.25 * bypass, 1e-6);
        assert_eq!(json.get("dmean_dinlet_k"), Some(&J::Null));
    }
}

#[test]
fn file_driven_total_gradient_matches_perturbed_coupled_fem() {
    let r = request();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let mut h: BTreeMap<_, _> = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let nominal = r.evaluate(cx, &flow, &h, true).unwrap();
        let delta = 1e-4_f64;
        h.insert("last-face".into(), 80.0 * delta.exp());
        let plus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        h.insert("last-face".into(), 80.0 * (-delta).exp());
        let minus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        close(nominal.gradient.unwrap().log_htc[1], (plus - minus) / (2.0 * delta), 2e-5);
    });
}

#[test]
fn declared_source_reaches_the_external_air_energy_balance() {
    let mut r = request();
    r.source = 2000.0;
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let e = r.evaluate(cx, &flow, &h, false).unwrap();
        close(e.source_total_w, 1.0, 1e-8);
        close(e.robin_total_w, e.source_total_w, r.limits.heat);
        close(e.coupled.transport.external_heat_gain_w, e.source_total_w, 2.0 * r.limits.heat);
    });
}

#[test]
fn reversing_the_declared_branch_orientation_preserves_the_field() {
    let mut r = request();
    r.gradient = false;
    let a = execute(&r, &CancelGate::new_clock_free()).unwrap();
    let mut edges = r.graph.branches().to_vec();
    let edge = &mut edges[0];
    std::mem::swap(&mut edge.from, &mut edge.to);
    r.region_paths[0].reverse();
    r.graph = LossGraph::new(4, edges).unwrap();
    let b = execute(&r, &CancelGate::new_clock_free()).unwrap();
    let a = J::parse(&a).unwrap();
    let b = J::parse(&b).unwrap();
    close(a.f64_field("objective_mean_k").unwrap(), b.f64_field("objective_mean_k").unwrap(), 1e-7);
    let branch = &b.get("branches").unwrap().as_array().unwrap()[0];
    close(branch.f64_field("flow_m3_s").unwrap(), -0.003, 1e-10);
}

#[test]
fn schema_dimensions_ownership_and_budget_mistakes_refuse() {
    for (from, to) in [
        ("\"SI\"", "\"engineering\""),
        ("\"schema\":", "\"unexpected\":true,\"schema\":"),
        ("\"seed\": \"41\"", "\"seed\": 41"),
        ("\"graph_sweeps\": 4096", "\"graph_sweeps\": 0"),
        ("\"node_count\": 4", "\"node_count\": 4.0"),
        ("\"htc_w_m2_k\": 50", "\"htc_w_m2_k\": -50"),
        ("\"last-face\"", "\"first-face\""),
        ("\"adiabatic_remainder\": true", "\"adiabatic_remainder\": false"),
    ] {
        assert!(FIXTURE.contains(from), "fixture edit did not apply: {from}");
        assert!(Request::parse(&FIXTURE.replace(from, to)).is_err(), "accepted {from} -> {to}");
    }
    assert!(Request::parse("{}").is_err());
    assert!(J::parse("{\"units\":\"SI\",\"units\":\"SI\"}").is_err());
}

#[test]
fn cancellation_and_exhaustion_publish_no_partial_solution() {
    let mut r = request();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(execute(&r, &gate).is_err());
    r.limits.coupling = 1;
    assert!(execute(&r, &CancelGate::new_clock_free()).is_err());
    let out = diagnostic(exit::BUDGET, producer("exhausted"), true);
    assert!(out.stdout.is_empty());
    assert_eq!(out.exit_code, exit::BUDGET);
    J::parse(&out.stderr).unwrap();
}

#[test]
fn result_json_strings_escape_controls_and_unicode() {
    let text = "quote\"slash\\line\ncontrol\u{1f}é";
    assert_eq!(J::parse(&quote(text)).unwrap().as_str(), Some(text));
    assert!(num(f64::NAN).is_err());
    assert!(num(f64::INFINITY).is_err());
}
