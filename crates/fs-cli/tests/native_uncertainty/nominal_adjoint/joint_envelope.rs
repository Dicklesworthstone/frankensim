//! Independent native ambient re-solves supply both complete coefficient ends.
use super::*;

fn value_after(detail: &str, label: &str) -> Vec<f64> {
    let marker = format!("{label} = ");
    detail.split(&marker).skip(1).map(|tail|
        tail.split(" K").next().unwrap().parse::<f64>().unwrap()).collect()
}

#[test]
fn native_full_design_matches_all_heated_and_cooled_coefficient_corners() {
    for cooled in [false, true] {
        let fixture = Fixture::new();
        let mut spec = fixture_project(true, 295.0, 305.0);
        if cooled {
            spec.power.as_mut().unwrap()[0].watts.value = 0.0;
            spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value = 280.0;
        }
        let original = spec.clone();
        let result = run(&fixture, &spec, 800);
        assert_eq!(spec, original, "propagation cannot mutate the declared model");
        let propagation = result.get("propagation").unwrap();
        let nominal = propagation.f64_field("nominal_base_k").unwrap();
        let boundary = propagation.get("boundary_conditions").unwrap();
        let model = propagation.get("model_form").unwrap();
        assert_eq!(boundary.str_field("state"), Some("measured"));
        assert_eq!(model.str_field("state"), Some("measured"));
        let detail = boundary.str_field("detail").unwrap();
        assert!(detail.contains("4/4 distinct joint corners"));
        let model_width = model.f64_field("half_width_k").unwrap();
        let boundary_width = boundary.get("vertices").unwrap().as_array().unwrap().iter()
            .map(|row| (row.f64_field("t_max_k").unwrap()-nominal).abs()).fold(0.0, f64::max);
        let allowance = fs_convection::correlation_catalog().into_iter().find(|card|
            card.id.name() == "convection.churchill-chu-vertical-plate").unwrap().model.discrepancy_rel;
        let mut joint_width = 0.0_f64;
        for (i, ambient) in [295.0,305.0].into_iter().enumerate() {
            let mut changed = spec.clone();
            perturb(&mut changed, "natural-convection-ambient", ambient-300.0);
            // Keep envelope, source card, reservoir and all physical parameters
            // unchanged. This independent native model produces both h(T) ends
            // at the chosen fluid ambient, not a frozen nominal Robin law.
            let actual = run(&fixture, &changed, 801+i);
            let term = actual.get("propagation").unwrap().get("model_form").unwrap();
            assert_eq!(term.str_field("state"), Some("measured"));
            let vertices = term.get("vertices").unwrap().as_array().unwrap();
            assert_eq!(vertices.len(),2);
            for (j, scale) in [1.0-allowance,1.0+allowance].into_iter().enumerate() {
                let expected = vertices[j].f64_field("t_max_k").unwrap();
                joint_width = joint_width.max((expected-nominal).abs());
                let label = format!("joint fluid {ambient} K, fan pressure x1, card coefficient x{scale}");
                let retained = value_after(detail,&label);
                assert!(!retained.is_empty(), "missing actual joint result: {label}");
                for value in retained {
                    assert!((value-expected).abs() < 1e-7,
                        "cooled={cooled}: retained {value} differs from independent native corner {expected}");
                }
            }
        }
        let excess = (joint_width-boundary_width-model_width).max(0.0);
        assert!((propagation.f64_field("interaction_excess_k").unwrap()-excess).abs() < 2e-7);
        assert!((boundary.f64_field("half_width_k").unwrap()-boundary_width-excess).abs() < 2e-7);
    }
}

#[test]
fn unavailable_model_allowance_is_not_reported_as_measured_zero_interaction() {
    let fixture = Fixture::new();
    let spec = fixture_project(false,295.0,305.0);
    let result = run(&fixture,&spec,820);
    let propagation = result.get("propagation").unwrap();
    assert_eq!(propagation.get("boundary_conditions").unwrap().str_field("state"),Some("measured"));
    assert_eq!(propagation.get("model_form").unwrap().str_field("state"),Some("no-data"));
    assert_eq!(propagation.get("interaction_excess_k"),Some(&JsonValue::Null));
}
