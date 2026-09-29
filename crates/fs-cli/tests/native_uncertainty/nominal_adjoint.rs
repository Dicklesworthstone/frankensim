//! Actual native import/material/solve/receipt comparisons, not a solver stub.
use super::*;

fn request(project: &mut fs_project::ProjectSpec) {
    project.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
        name: "temperature-max-adjoint".into(), kind: "report".into(),
    });
}

fn solve(fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize) -> JsonValue {
    let source = fixture.dir.join(format!("adjoint-{ordinal}.fsim"));
    std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
    let ledger = fixture.dir.join("adjoint.db");
    command(&["--json", "import", source.to_str().unwrap(),
        fixture.sources.join("plate.stl").to_str().unwrap(), ledger.to_str().unwrap(),
        "--unit", "m", "--max-hole-edges", "0"], fs_cli::exit::SUCCESS);
    let run = command(&["--json", "solve", source.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", fixture.sources.join("aa6061.fsmcdpk").to_str().unwrap()], fs_cli::exit::SUCCESS);
    let receipt = json_artifact(&ledger, run.str_field("run_receipt").unwrap());
    let conduction = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|stage| stage.str_field("stage") == Some("conduction")).unwrap();
    json_artifact(&ledger, conduction.str_field("receipt").unwrap())
}

fn coefficient(receipt: &JsonValue, target: &str) -> f64 {
    receipt.get("nominal_adjoint").unwrap().get("parameters").unwrap().as_array().unwrap()
        .iter().find(|row| row.str_field("target") == Some(target)).unwrap()
        .f64_field("derivative").unwrap()
}

#[test]
fn requested_native_adjoint_keeps_the_identical_primal_field() {
    let fixture = Fixture::new();
    let baseline = solve(&fixture, &fixture.project, 0);
    assert!(baseline.get("nominal_adjoint").is_none());
    let mut project = fixture.project.clone();
    request(&mut project);
    let with_adjoint = solve(&fixture, &project, 1);
    assert_eq!(baseline.get("temperature"), with_adjoint.get("temperature"));
    assert_eq!(baseline.get("energy"), with_adjoint.get("energy"));
    // Compare actual retained fields, not only their extrema. Run identity is
    // intentionally different because the output request is part of the model.
    let ledger = fixture.dir.join("adjoint.db");
    let a = json_artifact(&ledger, baseline.str_field("solution_artifact").unwrap());
    let b = json_artifact(&ledger, with_adjoint.str_field("solution_artifact").unwrap());
    assert_eq!(a.get("temperature"), b.get("temperature"));
    let gradient = with_adjoint.get("nominal_adjoint").unwrap();
    assert_eq!(gradient.str_field("authority"), Some("Estimated"));
    assert_eq!(gradient.str_field("functional"), Some("selected-nodal-temperature"));
}

#[test]
fn native_power_derivative_matches_fresh_physics_and_retains_duty_factor() {
    let fixture = Fixture::new();
    let mut base = fixture.project.clone();
    base.power.as_mut().unwrap()[0].duty = 0.37;
    let mut requested = base.clone();
    request(&mut requested);
    let nominal = solve(&fixture, &requested, 10);
    let derivative = coefficient(&nominal, "power");
    let h = 0.01;
    let mut low = base.clone();
    low.power.as_mut().unwrap()[0].watts.value -= h;
    let mut high = base;
    high.power.as_mut().unwrap()[0].watts.value += h;
    let a = solve(&fixture, &low, 11);
    let b = solve(&fixture, &high, 12);
    let difference = (b.get("temperature").unwrap().f64_field("max").unwrap()
        - a.get("temperature").unwrap().f64_field("max").unwrap()) / (2.0*h);
    assert!(derivative > 0.0);
    assert!((derivative-difference).abs() < 1e-6 * derivative.abs().max(1.0),
        "adjoint {derivative} vs native physical finite difference {difference}");
}

#[test]
fn native_convection_coefficients_match_fresh_boundary_resolves() {
    let mut fixture = Fixture::new();
    let fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. }
        = &mut fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
        else { panic!("native reference boundary"); };
    reference_temperature.value = 300.0;
    let mut requested = fixture.project.clone();
    request(&mut requested);
    let nominal = solve(&fixture, &requested, 20);
    for (case, target, h) in [(0, "convection-coefficient", 0.001), (1, "convection-temperature", 0.01)] {
        let derivative = coefficient(&nominal, target);
        let mut results = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut project = fixture.project.clone();
            let fs_project::ThermalBoundaryCondition::Convection { coefficient, reference_temperature }
                = &mut project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition
                else { panic!("native reference boundary"); };
            if case == 0 { coefficient.value += sign*h; }
            else { reference_temperature.value += sign*h; }
            results.push(solve(&fixture, &project, 21+case*2+side).get("temperature").unwrap()
                .f64_field("max").unwrap());
        }
        let difference = (results[1]-results[0])/(2.0*h);
        assert!((derivative-difference).abs() < 1e-5 * derivative.abs().max(1.0),
            "{target}: adjoint {derivative} vs native physical finite difference {difference}");
    }
}
