use super::*;

const MODEL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-hotspot.json"));

fn model() -> Request {
    Request::parse(&MODEL.replace("\"max_solid_temperature\": true", "\"mean_wall_region\": \"first-face\"")
        .replace("\"gradient\": false", "\"gradient\": true")).unwrap()
}
fn policy(limit: f64, evaluations: usize) -> String {
    format!(r#"{{"schema":"frankensim.cooling-resistance-design.v1",
        "temperature_limit_k":{limit},"temperature_scale_k":0.1,
        "controls":[
            {{"branch":"first-stream","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":1}},
            {{"branch":"fresh-bypass","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":0.3}}],
        "max_evaluations":{evaluations},"max_iterations":80,"max_multiplier_updates":40,
        "max_backtracks":24,"kkt_tolerance":0.000001}}"#)
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 41, kernel_id: 718, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1.0), "{a:.14e} != {b:.14e}");
}

#[test]
fn design_sample_derivatives_match_actual_changed_resistance_solves() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let design = Specification::parse(&policy(302.0, 160), &request).unwrap();
        let point = [0.1, -0.2];
        let (nominal, _, _) = evaluate(&mut request, &original, &design, &point, cx).unwrap();
        for index in 0..2 {
            let h = 1e-3;
            let mut plus = point; plus[index] += h;
            let mut minus = point; minus[index] -= h;
            let (a, _, _) = evaluate(&mut request, &original, &design, &plus, cx).unwrap();
            let (b, _, _) = evaluate(&mut request, &original, &design, &minus, cx).unwrap();
            close(nominal.gradient[index], (a.objective-b.objective)/(2.0*h), 1e-8);
            close(nominal.constraint_gradient[index], (a.constraint-b.constraint)/(2.0*h), 3e-4);
        }
        // Restoring the same point must not compound prior multipliers.
        let (again, _, _) = evaluate(&mut request, &original, &design, &point, cx).unwrap();
        assert_eq!(nominal, again);
        assert_eq!(request.graph.branches()[2], original.branches()[2]);
    });
}

#[test]
fn initially_feasible_nominal_design_costs_zero_and_reserves_replay() {
    with_cx(|cx| {
        let request = model();
        let design = Specification::parse(&policy(1000.0, 2), &request).unwrap();
        let output = J::parse(&solve(request, design, cx).unwrap()).unwrap();
        let result = output.get("resistance_design").unwrap();
        assert_eq!(result.str_field("status"), Some("locally-converged"));
        assert_eq!(result.f64_field("selected_cost"), Some(0.0));
        assert_eq!(result.f64_field("complete_evaluations"), Some(2.0));
        assert_eq!(result.f64_field("replay_evaluations"), Some(1.0));
        for control in result.get("controls").unwrap().as_array().unwrap() {
            assert_eq!(control.f64_field("multiplier"), Some(1.0));
        }
        assert!(output.get("hydraulic_sensitivity").is_some());
        assert!(output.get("fan").is_some());
    });
}

#[test]
fn constrained_design_publishes_a_recomputed_passing_coupled_field() {
    with_cx(|cx| {
        let mut request = model();
        let original = request.graph.clone();
        let design = Specification::parse(&policy(302.0, 240), &request).unwrap();
        let (_, _, initial) = evaluate(&mut request, &original, &design, &[0.0; 2], cx).unwrap();
        let (_, _, improved) = evaluate(&mut request, &original, &design, &[-0.7, 0.7], cx).unwrap();
        assert!(improved.objective < initial.objective);
        let limit = 0.5 * initial.objective + 0.5 * improved.objective;
        let request = model();
        let design = Specification::parse(&policy(limit, 240), &request).unwrap();
        let output = J::parse(&solve(request, design, cx).unwrap()).unwrap();
        let result = output.get("resistance_design").unwrap();
        let value = output.path(&["objective", "value_k"]).unwrap().as_f64().unwrap();
        assert!(value <= limit);
        assert!(result.f64_field("selected_cost").unwrap() > 0.0);
        assert!(result.f64_field("complete_evaluations").unwrap() <= 240.0);
        let mut replay = model();
        let mut branches = replay.graph.branches().to_vec();
        for control in result.get("controls").unwrap().as_array().unwrap() {
            let multiplier = control.f64_field("multiplier").unwrap();
            assert!((0.25..=4.0).contains(&multiplier));
            let branch = branches.iter_mut().find(|b| Some(b.loss.name.as_str()) == control.str_field("branch")).unwrap();
            branch.loss.resistance = LossResistance::new(control.f64_field("selected_resistance_pa_s2_m6").unwrap());
        }
        replay.graph = LossGraph::new(replay.graph.node_count(), branches).unwrap();
        let flow = replay.flow(cx).unwrap();
        let coefficients = replay.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let solved = replay.evaluate(cx, &flow, &coefficients, false).unwrap();
        close(solved.objective, value, 1e-12);
        assert_eq!(output.get("solid_temperatures_k").unwrap().as_array().unwrap().len(), solved.temperatures.len());
        for (reported, actual) in output.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().zip(solved.temperatures) {
            close(reported.as_f64().unwrap(), actual, 1e-12);
        }
    });
}

#[test]
fn exhausted_work_without_a_passing_sample_is_not_infeasibility() {
    with_cx(|cx| {
        let request = model();
        let design = Specification::parse(&policy(1.0, 2), &request).unwrap();
        let error = solve(request, design, cx).unwrap_err();
        assert_eq!(error.code, "cooling-network-design-no-feasible-sample");
        assert!(error.message.contains("not a proof of infeasibility"));
    });
}

#[test]
fn unsupported_models_and_invalid_controls_refuse_before_optimizer_work() {
    let request = model();
    let text = policy(302.0, 20);
    for (from, to) in [
        ("\"max_evaluations\":20", "\"max_evaluations\":1"),
        ("\"first-stream\"", "\"not-a-branch\""),
        ("\"fresh-bypass\"", "\"first-stream\""),
        ("\"minimum_multiplier\":0.25", "\"minimum_multiplier\":1.5"),
        ("\"maximum_multiplier\":4", "\"maximum_multiplier\":0.5"),
        ("\"weight\":1", "\"weight\":0"),
        ("\"temperature_scale_k\":0.1", "\"temperature_scale_k\":0"),
        ("\"kkt_tolerance\":0.000001", "\"kkt_tolerance\":1"),
        ("\"max_backtracks\":24", "\"max_backtracks\":65"),
    ] {
        assert!(text.contains(from));
        assert!(Specification::parse(&text.replace(from, to), &request).is_err(), "{to}");
    }
    let mut disabled = model(); disabled.gradient = false;
    assert!(Specification::parse(&text, &disabled).is_err());
    let peak = Request::parse(&MODEL.replace("\"gradient\": false", "\"gradient\": true")).unwrap();
    assert!(Specification::parse(&text, &peak).is_err());
}

#[test]
fn cancellation_never_publishes_a_feasible_but_incomplete_run() {
    let request = model();
    let design = Specification::parse(&policy(1000.0, 2), &request).unwrap();
    let gate = CancelGate::new_clock_free(); gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 41, kernel_id: 718, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        assert_eq!(solve(request, design, &cx).unwrap_err().code, "cooling-network-cancelled");
    });
}
