//! Natural uncertainty vertices must re-evaluate h(T), not perturb a frozen h.
use super::*;

fn fixture(fanless: bool) -> Fixture {
    let mut f = Fixture::new();
    f.project.power.as_mut().unwrap()[0].watts.value = 12.0;
    f.project.power.as_mut().unwrap()[0].duty = 0.37;
    f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    let cooling = f.project.cooling.as_mut().unwrap();
    if fanless {
        cooling.fans.clear();
        cooling.vents.clear();
        cooling.fan_system = None;
        cooling.airflow_leakage = None;
    }
    cooling.conduction.as_mut().unwrap().boundaries[0].condition =
        fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    f
}

#[test]
fn natural_card_budget_measures_full_fixed_point_vertices_not_frozen_robin_solves() {
    for fanless in [false, true] {
        let f = fixture(fanless);
        let result = solve(&f, &f.project, 400);
        let propagation = result.get("propagation").unwrap();
        let nominal = propagation.f64_field("nominal_base_k").unwrap();
        let model = propagation.get("model_form").unwrap();
        assert_eq!(model.str_field("state"), Some("measured"));
        assert_eq!(model.str_field("method"), Some("interval-vertex-resolve"));
        let vertices = model.get("vertices").unwrap().as_array().unwrap();
        assert_eq!(vertices.len(), 2);
        let hot = vertices[0].f64_field("t_max_k").unwrap();
        let cold = vertices[1].f64_field("t_max_k").unwrap();
        assert!(hot > nominal && cold < nominal, "the allowance must change the actual physical answer");
        let width = (hot - nominal).abs().max((cold - nominal).abs());
        assert_eq!(model.f64_field("half_width_k").unwrap().to_bits(), width.to_bits());
        assert!(model.str_field("detail").unwrap().contains("not a validated interval"));
        let natural = result.get("natural_convection").unwrap();
        let law = &natural.get("laws").unwrap().as_array().unwrap()[0];
        assert!(law.get("htc_multiplier").is_none(), "the published nominal law is not perturbed");
        let raw_h = law.f64_field("htc_w_m2_k").unwrap();
        let allowance = fs_convection::correlation_catalog().into_iter().find(|card|
            card.id.name() == "convection.churchill-chu-vertical-plate").unwrap().model.discrepancy_rel;
        // Negative physical control: an independently supplied constant Robin
        // coefficient at the nominal wall is a DIFFERENT model. If a vertex
        // drops feedback, these ordinary native re-solves would coincide.
        for (index, (scale, actual)) in [(1.0 - allowance, hot), (1.0 + allowance, cold)].into_iter().enumerate() {
            let mut frozen = f.project.clone();
            frozen.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
                fs_project::ThermalBoundaryCondition::Convection {
                    coefficient: fs_qty::QtyAny::new(raw_h * scale, fs_project::spec::dims::HEAT_TRANSFER_COEFFICIENT),
                    reference_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
                };
            let wrong = solve(&f, &frozen, 401 + index).get("temperature").unwrap().f64_field("max").unwrap();
            assert!((wrong - actual).abs() > 1e-6, "vertex silently used a frozen natural coefficient");
        }
    }
}

#[test]
fn natural_error_analysis_cannot_polish_or_certify_the_frozen_operator() {
    let f = fixture(true);
    let result = solve(&f, &f.project, 410);
    let control = result.get("solver_control").unwrap();
    assert_eq!(control.str_field("status"), Some("unsupported-model"));
    assert_eq!(control.get("goal_met"), Some(&JsonValue::Bool(false)));
    assert_eq!(control.f64_field("primal_iterations"), Some(0.0));
    assert!(control.str_field("reason").unwrap().contains("natural-convection"));
    assert!(result.get("solver_algebraic").is_none(), "no linear enclosure is applicable");
    let tolerance = result.get("propagation").unwrap().get("solver_algebraic").unwrap();
    assert_eq!(tolerance.str_field("method"), Some("tolerance-tightening-resolve"));
    let roundoff = result.get("roundoff").unwrap();
    assert_eq!(roundoff.str_field("state"), Some("no-data"));
    assert!(roundoff.str_field("reason").unwrap().contains("coefficient-feedback"));
    // The retained wall mean and heat balance belong to the unchanged field.
    let law = &result.get("natural_convection").unwrap().get("laws").unwrap().as_array().unwrap()[0];
    assert!((law.f64_field("heat_rate_w").unwrap()
        - result.get("energy").unwrap().f64_field("robin_out_w").unwrap()).abs() < 1e-9);
}
