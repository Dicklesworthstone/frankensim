//! Full native natural-convection adaptation with a nonlinear material card.
//! The reported two-mesh error remains Estimated, not a continuum bound.
use super::*;

#[test]
fn native_natural_adaptive_goal_uses_total_feedback_and_publishes_the_probed_field() {
    let mut fixture = Fixture::new();
    let project = &mut fixture.project;
    project.power.as_mut().unwrap()[0].watts.value = 5.0;
    project.power.as_mut().unwrap()[0].duty = 1.0;
    project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    project.solver.as_mut().unwrap().fidelity = "adaptive".into();
    project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    project.budgets.as_mut().unwrap().accuracy_rel = 0.2;
    let (receipt, field) = fixture.solve(&fixture.project, 700);
    let adaptive = receipt.get("adaptive").expect("native adaptive producer executed");
    let history = adaptive.get("history").unwrap().as_array().unwrap();
    assert!(!history.is_empty(), "at least one real enriched solve must be compared");
    let goal_tolerance = fixture.project.solver.as_ref().unwrap().tolerance_rel * 1e-2;
    for row in history {
        assert_eq!(row.get("natural_coefficient_feedback"), Some(&JsonValue::Bool(true)));
        assert_eq!(row.get("uses_nonlinear_jacobian"), Some(&JsonValue::Bool(true)));
        assert!(row.f64_field("primal_residual").unwrap() < goal_tolerance);
        assert!(row.f64_field("dual_residual").unwrap() < goal_tolerance);
        let signed = row.f64_field("signed_linear_change_k").unwrap();
        let nonlinear = row.f64_field("linearization_remainder_k").unwrap();
        let owner = row.f64_field("maximum_remainder_k").unwrap();
        let measured = row.f64_field("measured_change_k").unwrap();
        assert!((signed + nonlinear + owner - measured).abs() < 1e-9,
            "the total residual, nonlinear and maximum-owner terms must match the observed goal");
        assert!(row.f64_field("estimated_change_k").unwrap() + 1e-12 >= measured.abs());
        assert!(row.f64_field("probe_tets").unwrap() > row.f64_field("tets").unwrap());
    }
    let temperatures = field.get("temperature").unwrap().as_array().unwrap();
    let maximum = temperatures.iter().map(|v| v.as_f64().unwrap()).fold(f64::NEG_INFINITY, f64::max);
    let last = history.last().unwrap();
    assert!((maximum - last.f64_field("t_max_k").unwrap()).abs() < 1e-10,
        "publish the last independently probed mesh, not the provisional enrichment");
    assert!(maximum > 300.0 && maximum < 450.0);
    // Replaying the same canonical input must not replace the accepted field
    // with a different law, tolerance or provisional mesh.
    let (replayed, replayed_field) = fixture.solve(&fixture.project, 701);
    assert_eq!(receipt.get("adaptive"), replayed.get("adaptive"));
    assert_eq!(field.get("temperature"), replayed_field.get("temperature"));
}
