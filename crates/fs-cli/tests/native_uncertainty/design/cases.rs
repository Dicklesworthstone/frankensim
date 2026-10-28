//! Actual native case sweeps, measured feasible incumbents, and partial-case
//! recovery. The existing physical fixture/solver remains the oracle.
use super::*;

fn source_with_cases(source: &str, target: &str) -> String {
    let prefix = source.trim_end().strip_suffix(')').unwrap();
    format!(r#"{prefix}
      :operating-cases (
        (case :name "cool" :inputs (({target} :entity "air" :value 293.15K)))
        (case :name "hot" :inputs (({target} :entity "air" :value 303.15K)))))"#)
}

fn set_ambient(project: &mut fs_project::ProjectSpec, ambient: f64) {
    for boundary in &mut project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        match &mut boundary.condition {
            fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. } =>
                reference_temperature.value = ambient,
            fs_project::ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } =>
                inlet_temperature.value = ambient,
            _ => panic!("unexpected case fixture boundary"),
        }
    }
}

fn cases_fixture(fan: bool) -> Fixture {
    let mut fixture = if fan { fan_fixture(false) } else { fixture(false) };
    let mut oracle = fixture.project.clone();
    set_ambient(&mut oracle, 303.15);
    if fan {
        oracle.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = 1.0;
    } else {
        oracle.power.as_mut().unwrap()[0].watts.value = 16.0; // 8 delivered W.
    }
    let (_, _, target) = fixture.solve_project(900, &oracle);
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = target + 1.0;
    requirement.margin.value = 1.0;
    write(&fixture);
    let source = if fan { fan_source() } else { DESIGN.into() };
    let target = if fan { "air-inlet-temperature" } else { "convection-temperature" };
    std::fs::write(fixture.sources.join("study.fsim"), source_with_cases(&source, target)).unwrap();
    fixture
}

fn case_report(fixture: &Fixture, result: &JsonValue, fan: bool) -> JsonValue {
    let receipt = result.get("receipt").unwrap();
    assert_eq!(receipt.str_field("driver"), Some(if fan {
        "native-cooling-fan-case-design-v1"
    } else { "native-cooling-power-case-design-v1" }));
    json_artifact(&fixture.ledger, receipt.str_field("report_json").unwrap())
}

#[test]
fn finite_case_power_and_fan_designs_pass_every_authored_condition() {
    for fan in [false, true] {
        let fixture = cases_fixture(fan);
        let result = fixture.study(None, fs_cli::exit::SUCCESS);
        let report = case_report(&fixture, &result, fan);
        assert_eq!(report.str_field("termination"), Some("target-bracketed"));
        assert_eq!(report.str_field("robustness_scope"), Some("authored-finite-cases-only"));
        assert_eq!(report.str_field("selected_operating_case"), Some("hot"));
        assert_eq!(report.str_field("selected_project_scope"), Some("evaluated-governing-operating-case"));
        assert_eq!(report.get("pending_candidate"), Some(&JsonValue::Null));
        let selected = report.get("selected").unwrap();
        let key = if fan { "fan_speed_ratio" } else { "power_w" };
        let coordinate = selected.f64_field(key).unwrap();
        assert!((coordinate - if fan { 1.0 } else { 8.0 }).abs() <= if fan { 0.001 } else { 0.05 });
        let limit = report.f64_field("temperature_limit_k").unwrap();
        let history = report.get("history").unwrap().as_array().unwrap();
        let evaluations = report.get("case_evaluations").unwrap().as_array().unwrap();
        assert_eq!(evaluations.len(), history.len() * 2);
        assert_eq!(report.f64_field("evaluations_completed"), Some(evaluations.len() as f64));
        assert_eq!(report.f64_field("candidates_completed"), Some(history.len() as f64));
        // Independently author and solve both physical conditions at the chosen
        // coordinate, without calling the new case binder or aggregation.
        for (ordinal, ambient) in [293.15, 303.15].into_iter().enumerate() {
            let mut project = fixture.project.clone();
            set_ambient(&mut project, ambient);
            if fan {
                project.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = coordinate;
            } else {
                let row = &mut project.power.as_mut().unwrap()[0];
                row.watts.value = coordinate / row.duty;
            }
            let (_, run, value) = fixture.solve_project(910 + ordinal, &project);
            assert!(value <= limit, "authored case {ambient} must be feasible");
            assert!(evaluations.iter().any(|row|
                row.path(&["observation", "run"]).and_then(JsonValue::as_str) == Some(run.as_str())));
        }
        let bytes = artifact(&fixture.ledger,
            result.path(&["receipt", "selected_project"]).unwrap().as_str().unwrap());
        let exported = fs_project::parse_sexpr(std::str::from_utf8(&bytes).unwrap()).unwrap();
        let mut expected = fixture.project.clone();
        set_ambient(&mut expected, 303.15);
        if fan {
            expected.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = coordinate;
        } else {
            let row = &mut expected.power.as_mut().unwrap()[0]; row.watts.value = coordinate / row.duty;
        }
        assert_eq!(exported.spec, expected, "export must be the measured governing case");
        assert_eq!(selected.str_field("project_hash"), Some(exported.hash().to_hex().as_str()));
    }
}

#[test]
fn incomplete_case_groups_resume_without_repeating_cases_or_selecting_partial_candidates() {
    let full = cases_fixture(true);
    let completed = full.study(None, fs_cli::exit::SUCCESS);
    let expected = case_report(&full, &completed, true);
    let split = cases_fixture(true);
    let first = split.study(Some("1"), fs_cli::exit::BUDGET);
    let report = case_report(&split, &first, true);
    assert_eq!(report.f64_field("evaluations_completed"), Some(1.0));
    assert_eq!(report.f64_field("candidates_completed"), Some(0.0));
    assert_eq!(report.path(&["pending_candidate", "cases_completed"]).unwrap().as_f64(), Some(1.0));
    assert_eq!(report.get("selected"), Some(&JsonValue::Null));
    assert_eq!(first.path(&["receipt", "selected_project"]), Some(&JsonValue::Null));
    let second = split.resume(first.str_field("run").unwrap(), Some("2"), fs_cli::exit::BUDGET);
    let report = case_report(&split, &second, true);
    assert_eq!(report.f64_field("evaluations_completed"), Some(3.0));
    assert_eq!(report.f64_field("candidates_completed"), Some(1.0));
    assert_eq!(report.path(&["selected", "fan_speed_ratio"]).unwrap().as_f64(), Some(1.2));
    assert_eq!(report.path(&["pending_candidate", "cases_completed"]).unwrap().as_f64(), Some(1.0));
    std::fs::rename(&split.sources, split.dir.join("relocated-case-design-inputs")).unwrap();
    let final_result = split.resume(second.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(case_report(&split, &final_result, true), expected,
        "case ordering, child identities, governing cases and selection must replay exactly");
}

#[test]
fn a_passing_first_case_does_not_hide_an_infeasible_baseline() {
    let mut fixture = cases_fixture(false);
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 301.0;
    requirement.margin.value = 1.0;
    std::fs::write(fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap()).unwrap();
    let first = fixture.study(Some("1"), fs_cli::exit::BUDGET);
    let partial = case_report(&fixture, &first, false);
    let raw = partial.get("case_evaluations").unwrap().as_array().unwrap();
    assert!(raw[0].path(&["observation", "temperature_k"]).unwrap().as_f64().unwrap() <= 300.0);
    assert_eq!(partial.get("selected"), Some(&JsonValue::Null));
    let failed = fixture.resume(first.str_field("run").unwrap(), None, fs_cli::exit::REFUSED);
    let report = case_report(&fixture, &failed, false);
    assert_eq!(report.str_field("termination"), Some("minimum-power-infeasible"));
    assert_eq!(report.f64_field("evaluations_completed"), Some(2.0));
    assert_eq!(report.f64_field("candidates_completed"), Some(1.0));
    assert_eq!(report.get("selected"), Some(&JsonValue::Null));
}
