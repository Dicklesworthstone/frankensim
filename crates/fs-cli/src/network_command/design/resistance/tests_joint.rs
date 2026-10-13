use super::*;

const MODEL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/resistance-model.json"));

fn model(peak: bool) -> Request {
    let text = if peak { MODEL.replace("\"mean_wall_region\": \"first-face\"", "\"max_solid_temperature\": true") }
        else { MODEL.to_string() };
    Request::parse(&text).unwrap()
}

fn policy(limit: f64, peak: bool, speed_only: bool) -> String {
    let controls = if speed_only { "[]" } else {
        r#"[{"branch":"first-stream","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":1},
            {"branch":"fresh-bypass","minimum_multiplier":0.25,"maximum_multiplier":4,"weight":0.3}]"#
    };
    let smoothing = if peak { ",\"maximum_smoothing_k\":0.04" } else { "" };
    format!(r#"{{"schema":"frankensim.cooling-resistance-design.v1",
        "temperature_limit_k":{limit},"temperature_scale_k":0.1,"controls":{controls},
        "fan_speed_control":{{"minimum_speed_ratio":0.6,"maximum_speed_ratio":1.8,"weight":0.2}},
        "max_evaluations":240,"max_iterations":100,"max_multiplier_updates":60,
        "max_backtracks":24,"kkt_tolerance":0.000001{smoothing}}}"#)
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
fn joint_temperature_and_cost_partials_match_full_changed_model_solves() {
    with_cx(|cx| for peak in [false, true] {
        let mut request = model(peak);
        let original = request.graph.clone();
        let design = Specification::parse(&policy(303.0, peak, false), &request).unwrap();
        assert_eq!(design.dimension(), 3);
        let point = [0.1, -0.2, 0.15];
        let (nominal, _, _) = evaluate(&mut request, &original, &design, &point, cx).unwrap();
        for i in 0..3 {
            let h = 1e-3;
            let mut plus = point; plus[i] += h;
            let mut minus = point; minus[i] -= h;
            let (a, _, _) = evaluate(&mut request, &original, &design, &plus, cx).unwrap();
            let (b, _, _) = evaluate(&mut request, &original, &design, &minus, cx).unwrap();
            close(nominal.gradient[i], (a.objective-b.objective)/(2.0*h), 1e-8);
            close(nominal.constraint_gradient[i], (a.constraint-b.constraint)/(2.0*h), 3e-4);
        }
        assert_eq!(evaluate(&mut request, &original, &design, &point, cx).unwrap().0, nominal);
    });
}

#[test]
fn speed_coordinates_are_relative_to_original_not_previous_trials() {
    let mut request = model(false);
    request.fan.as_mut().unwrap().speed_ratio = 0.9;
    let original = request.graph.clone();
    let design = Specification::parse(&policy(303.0, false, false), &request).unwrap();
    design.install(&mut request, &original, &[0.1, -0.2, 0.2]).unwrap();
    let first = request.graph.branches().to_vec();
    let speed = request.fan.as_ref().unwrap().speed_ratio;
    close(speed, 0.9 * det::exp(0.2), 1e-14);
    design.install(&mut request, &original, &[-0.4, 0.5, -0.1]).unwrap();
    design.install(&mut request, &original, &[0.1, -0.2, 0.2]).unwrap();
    assert_eq!(request.graph.branches(), first);
    assert_eq!(request.fan.as_ref().unwrap().speed_ratio.to_bits(), speed.to_bits());
    assert_eq!(request.graph.branches()[2], original.branches()[2]);
    assert!(design.install(&mut request, &original, &[0.0, 0.0, f64::NAN]).is_err());
    assert!(design.install(&mut request, &original, &[100.0, 0.0, 0.0]).is_err());
    assert!(design.install(&mut request, &original, &[0.0; 2]).is_err());
    assert_eq!(request.graph.branches(), first);
    assert_eq!(request.fan.as_ref().unwrap().speed_ratio.to_bits(), speed.to_bits());
}

#[test]
fn initially_feasible_joint_design_replays_every_coordinate() {
    with_cx(|cx| {
        let request = model(false);
        let text = policy(1000.0, false, false).replace("\"max_evaluations\":240", "\"max_evaluations\":2");
        let design = Specification::parse(&text, &request).unwrap();
        let output = J::parse(&solve(request, design, cx).unwrap()).unwrap();
        let report = output.get("resistance_design").unwrap();
        assert_eq!(report.str_field("status"), Some("locally-converged"));
        assert_eq!(report.f64_field("selected_cost"), Some(0.0));
        assert_eq!(report.f64_field("complete_evaluations"), Some(2.0));
        let fan = report.get("fan_speed_control").unwrap();
        assert_eq!(fan.f64_field("coordinate_index"), Some(2.0));
        assert_eq!(fan.f64_field("selected_speed_ratio"), Some(1.0));
        assert_eq!(output.path(&["fan", "speed_ratio"]).unwrap().as_f64(), Some(1.0));
        assert_eq!(report.path(&["optimizer_terminal", "point_log_multipliers"]).unwrap().as_array().unwrap().len(), 3);
    });
}

#[test]
fn speed_only_design_keeps_all_resistances_fixed() {
    with_cx(|cx| {
        let mut request = model(false);
        let original = request.graph.clone();
        let design = Specification::parse(&policy(1000.0, false, true), &request).unwrap();
        assert_eq!(design.dimension(), 1);
        let (sample, _, _) = evaluate(&mut request, &original, &design, &[0.2], cx).unwrap();
        assert_eq!(sample.gradient.len(), 1);
        assert_eq!(request.graph.branches(), original.branches());
        close(request.fan.as_ref().unwrap().speed_ratio, det::exp(0.2), 1e-14);
    });
}

