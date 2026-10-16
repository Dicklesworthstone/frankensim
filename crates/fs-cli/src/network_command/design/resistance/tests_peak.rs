use super::*;
use crate::network_command::tests::with_cx;

const MODEL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/resistance-hotspot-model.json"));
const POLICY: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/resistance-hotspot-design.json"));

fn model(selector: &str) -> Request {
    Request::parse(&MODEL.replace("\"max_solid_temperature\": true", selector)).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1.0), "{a:.14e} != {b:.14e}");
}

#[test]
fn hotspot_constraint_adjoint_matches_complete_perturbed_physics_for_all_selectors() {
    with_cx(|cx| {
        for selector in [r#""max_solid_temperature": true"#, r#""max_wall_region": "first-face""#, r#""max_vertices": [4,5]"#] {
            let mut request = model(selector);
            let original = request.graph.clone();
            let design = Specification::parse(POLICY, &request).unwrap();
            let point = [-0.15, 0.1];
            let (sample, _, physical) = evaluate(&mut request, &original, &design, &point, cx).unwrap();
            assert!(physical.objective >= physical.objective_state.physical_value());
            for index in 0..2 {
                let h = 1e-3;
                let mut plus = point; plus[index] += h;
                let mut minus = point; minus[index] -= h;
                let (a, _, _) = evaluate(&mut request, &original, &design, &plus, cx).unwrap();
                let (b, _, _) = evaluate(&mut request, &original, &design, &minus, cx).unwrap();
                close(sample.constraint_gradient[index], (a.constraint-b.constraint)/(2.0*h), 3e-4);
            }
            let (again, _, _) = evaluate(&mut request, &original, &design, &point, cx).unwrap();
            assert_eq!(sample, again);
        }
    });
}

#[test]
fn peak_design_replays_a_field_whose_every_selected_vertex_passes() {
    with_cx(|cx| {
        let selector = r#""max_solid_temperature": true"#;
        let mut request = model(selector);
        let original = request.graph.clone();
        let mut design = Specification::parse(POLICY, &request).unwrap();
        let (_, _, initial) = evaluate(&mut request, &original, &design, &[0.0; 2], cx).unwrap();
        let (_, _, cooler) = evaluate(&mut request, &original, &design, &[-0.7, 0.7], cx).unwrap();
        assert!(cooler.objective < initial.objective);
        design.limit = 0.5 * initial.objective + 0.5 * cooler.objective;
        let limit = design.limit;
        let output = J::parse(&solve(model(selector), design, cx).unwrap()).unwrap();
        let objective = output.get("objective").unwrap();
        let summary = output.get("resistance_design").unwrap();
        assert_eq!(objective.str_field("kind"), Some("smooth_upper_envelope"));
        assert_eq!(objective.str_field("maximum_selection"), Some("max_solid_temperature"));
        assert_eq!(summary.str_field("constraint_form"), Some("smooth-maximum-upper-envelope"));
        let maximum = objective.f64_field("actual_maximum_k").unwrap();
        let upper = objective.f64_field("value_k").unwrap();
        assert!(maximum <= upper && upper <= limit);
        assert!(summary.f64_field("selected_cost").unwrap() > 0.0);
        assert_eq!(summary.f64_field("replay_evaluations"), Some(1.0));
        assert_eq!(summary.f64_field("selected_temperature_k"), Some(maximum));
        for temperature in output.get("solid_temperatures_k").unwrap().as_array().unwrap() {
            assert!(temperature.as_f64().unwrap() <= limit);
        }
        // Rebuild the selected losses and evaluate the ORIGINAL unsmoothed model.
        let mut replay = model(selector);
        let mut branches = replay.graph.branches().to_vec();
        for control in summary.get("controls").unwrap().as_array().unwrap() {
            let branch = branches.iter_mut().find(|b| Some(b.loss.name.as_str()) == control.str_field("branch")).unwrap();
            branch.loss.resistance = LossResistance::new(control.f64_field("selected_resistance_pa_s2_m6").unwrap());
        }
        replay.graph = LossGraph::new(replay.graph.node_count(), branches).unwrap();
        let flow = replay.flow(cx).unwrap();
        let coefficients = replay.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let actual = replay.evaluate(cx, &flow, &coefficients, false).unwrap();
        close(actual.objective, maximum, 1e-12);
        for (reported, actual) in output.get("solid_temperatures_k").unwrap().as_array().unwrap().iter().zip(actual.temperatures) {
            close(reported.as_f64().unwrap(), actual, 1e-12);
        }
    });
}

#[test]
fn raw_peak_feasibility_does_not_silently_relax_the_declared_smooth_constraint() {
    with_cx(|cx| {
        let selector = r#""max_solid_temperature": true"#;
        let mut request = model(selector);
        let original = request.graph.clone();
        let mut design = Specification::parse(POLICY, &request).unwrap();
        let (_, _, initial) = evaluate(&mut request, &original, &design, &[0.0; 2], cx).unwrap();
        let maximum = initial.objective_state.physical_value();
        assert!(initial.objective > maximum);
        design.limit = 0.5 * maximum + 0.5 * initial.objective;
        design.evaluations = 2; design.options.max_evaluations = 1;
        assert!(maximum < design.limit);
        assert!(!design.passes(&initial));
        assert_eq!(solve(model(selector), design, cx).unwrap_err().code, "cooling-network-design-no-feasible-sample");
    });
}

#[test]
fn smoothing_requires_an_explicit_peak_selector_and_positive_policy() {
    let peak = model(r#""max_solid_temperature": true"#);
    let mean = model(r#""mean_wall_region": "first-face""#);
    assert!(Specification::parse(POLICY, &mean).is_err());
    let without = POLICY.replace("  \"maximum_smoothing_k\": 0.04,\n", "");
    assert_ne!(without, POLICY);
    assert!(Specification::parse(&without, &peak).is_err());
    assert!(Specification::parse(&without, &mean).is_ok());
    for tau in ["0", "-1", "1e999", "null"] {
        assert!(Specification::parse(&POLICY.replace("\"maximum_smoothing_k\": 0.04", &format!("\"maximum_smoothing_k\": {tau}")), &peak).is_err());
    }
}

#[test]
fn singleton_constraint_keeps_the_original_selected_vertex_value_and_gradient() {
    with_cx(|cx| {
        let mut request = model(r#""max_vertices": [4]"#);
        let original = request.graph.clone();
        let flow = request.flow(cx).unwrap();
        let coefficients = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let plain = request.evaluate(cx, &flow, &coefficients, true).unwrap();
        let design = Specification::parse(POLICY, &request).unwrap();
        let (_, _, smooth) = evaluate(&mut request, &original, &design, &[0.0; 2], cx).unwrap();
        assert_eq!(plain.objective, smooth.objective);
        for index in 0..request.graph.branches().len() {
            assert_eq!(plain.gradient.as_ref().unwrap().log_resistance(index),
                smooth.gradient.as_ref().unwrap().log_resistance(index));
        }
    });
}
