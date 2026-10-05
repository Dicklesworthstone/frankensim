//! Exercise the real native import, card, nonlinear solve and receipt path.
use super::*;

#[test]
fn native_natural_adjoint_preserves_the_primal_and_matches_physical_controls() {
    let mut fixture = Fixture::new();
    fixture.project.power.as_mut().unwrap()[0].watts.value = 12.0;
    fixture.project.power.as_mut().unwrap()[0].duty = 0.37;
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    let baseline = solve(&fixture, &fixture.project, 300);
    let mut requested = fixture.project.clone();
    request(&mut requested);
    let nominal = solve(&fixture, &requested, 301);
    assert_eq!(baseline.get("temperature"), nominal.get("temperature"));
    assert_eq!(baseline.get("energy"), nominal.get("energy"));
    assert_eq!(baseline.get("natural_convection"), nominal.get("natural_convection"));
    let ledger = fixture.dir.join("adjoint.db");
    let a = json_artifact(&ledger, baseline.str_field("solution_artifact").unwrap());
    let b = json_artifact(&ledger, nominal.str_field("solution_artifact").unwrap());
    assert_eq!(a.get("temperature"), b.get("temperature"), "differentiate the retained field, not a second primal");
    let gradient = nominal.get("nominal_adjoint").unwrap();
    assert_eq!(gradient.str_field("mode"), Some("natural-convection-full-wall-feedback"));
    assert_eq!(gradient.str_field("authority"), Some("Estimated"));
    assert!(gradient.get("unsupported").unwrap().as_array().unwrap().iter()
        .all(|row| row.str_field("target") != Some("natural-convection-ambient")));
    let maximum = |receipt: &JsonValue| receipt.get("temperature").unwrap().f64_field("max").unwrap();
    for (case, target, step) in [(0, "power", 0.02), (1, "natural-convection-ambient", 0.005)] {
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut project = fixture.project.clone();
            if case == 0 { project.power.as_mut().unwrap()[0].watts.value += sign*step; }
            else {
                let fs_project::ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. }
                    = &mut project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
                    else { panic!("natural fixture"); };
                ambient_temperature.value += sign*step;
            }
            values.push(maximum(&solve(&fixture, &project, 302+case*2+side)));
        }
        let expected = (values[1]-values[0])/(2.0*step);
        let actual = coefficient(&nominal, target);
        assert!(actual > 0.0);
        assert!((actual-expected).abs() < 1e-4*expected.abs().max(1e-3),
            "{target}: full natural adjoint {actual:e}, physical finite difference {expected:e}");
    }
}
