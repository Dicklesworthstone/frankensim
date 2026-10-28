//! Native enclosure derivatives and probability calibration share one physical owner.
use super::*;

fn fixture() -> Fixture {
    let mut fixture = Fixture::new();
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
    fixture.project = fs_project::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-enclosure.fsim"
    ))).unwrap().decoded.spec;
    fixture.project.envelope.as_mut().unwrap().ambient_lo.value = 280.0;
    fixture.project.envelope.as_mut().unwrap().ambient_hi.value = 320.0;
    let enclosure = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap();
    enclosure.temperature_tolerance.value = 1e-11;
    enclosure.heat_tolerance.value = 1e-10;
    std::fs::copy(reference.join("gray-surface.fsmcdpk"),
        fixture.sources.join("gray-surface.fsmcdpk")).unwrap();
    std::fs::write(fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap()).unwrap();
    fixture
}

fn execute(
    fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize, expected: u8,
) -> JsonValue {
    let source = fixture.dir.join(format!("enclosure-adjoint-{ordinal}.fsim"));
    std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
    let ledger = fixture.dir.join("enclosure-adjoint.db");
    command(&["--json", "import", source.to_str().unwrap(),
        fixture.sources.join("plate.stl").to_str().unwrap(), ledger.to_str().unwrap(),
        "--unit", "m", "--max-hole-edges", "0"], fs_cli::exit::SUCCESS);
    command(&["--json", "solve", source.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", fixture.sources.join("aa6061.fsmcdpk").to_str().unwrap(),
        "--materials", fixture.sources.join("gray-surface.fsmcdpk").to_str().unwrap()], expected)
}

fn native(fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize) -> JsonValue {
    let solved = execute(fixture, project, ordinal, fs_cli::exit::SUCCESS);
    let ledger = fixture.dir.join("enclosure-adjoint.db");
    let receipt = json_artifact(&ledger, solved.str_field("run_receipt").unwrap());
    let stage = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    json_artifact(&ledger, stage.str_field("receipt").unwrap())
}

fn field(fixture: &Fixture, receipt: &JsonValue) -> JsonValue {
    json_artifact(&fixture.dir.join("enclosure-adjoint.db"),
        receipt.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone()
}

#[test]
fn native_enclosure_maximum_adjoint_preserves_the_field_and_matches_physical_resolves() {
    let fixture = fixture();
    let mut base = fixture.project.clone();
    base.power.as_mut().unwrap()[0].duty = 0.37;
    let original = native(&fixture, &base, 0);
    let mut requested = base.clone();
    request(&mut requested);
    let nominal = native(&fixture, &requested, 1);
    assert_eq!(field(&fixture, &original), field(&fixture, &nominal));
    assert_eq!(original.get("energy"), nominal.get("energy"));
    let adjoint = nominal.get("nominal_adjoint").unwrap();
    assert_eq!(adjoint.str_field("mode"), Some("enclosure-full-radiosity-feedback"));
    assert_eq!(adjoint.str_field("authority"), Some("Estimated"));
    assert_eq!(adjoint.get("stability_iterations"), Some(&JsonValue::Null));
    let selected = adjoint.f64_field("selected_vertex").unwrap() as usize;
    let emissivities: Vec<_> = adjoint.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field("target") == Some("radiation-emissivity")).collect();
    assert_eq!(emissivities.len(), 2);
    assert!(emissivities.iter().all(|row| row.f64_field("derivative").unwrap().is_finite()));
    assert!(adjoint.get("parameters").unwrap().as_array().unwrap().iter().all(|row|
        row.str_field("target") != Some("radiation-reservoir-temperature")));
    for (case, target, step) in [(0, "power", 0.02),
        (1, "convection-coefficient", 0.001), (2, "convection-temperature", 0.01)]
    {
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut perturbed = base.clone();
            if case == 0 {
                perturbed.power.as_mut().unwrap()[0].watts.value += sign * step;
            } else {
                let fs_project::ThermalBoundaryCondition::Convection {
                    coefficient, reference_temperature,
                } = &mut perturbed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .boundaries[0].condition else { panic!("prescribed convection fixture"); };
                if case == 1 { coefficient.value += sign * step; }
                else { reference_temperature.value += sign * step; }
            }
            let receipt = native(&fixture, &perturbed, 2 + 2 * case + side);
            // Keep the selected nodal functional fixed, even if the maximum
            // were to relocate under one of the physical perturbations.
            values.push(field(&fixture, &receipt).as_array().unwrap()[selected].as_f64().unwrap());
        }
        let expected = (values[1] - values[0]) / (2.0 * step);
        let actual = coefficient(&nominal, target);
        assert!((actual - expected).abs() < 5e-5 * expected.abs().max(0.02),
            "{target}: enclosure adjoint {actual:e} versus native differences {expected:e}");
    }
}

#[test]
fn native_enclosure_volume_mean_uses_its_own_goal_and_complete_radiation_dual() {
    let fixture = fixture();
    let original = native(&fixture, &fixture.project, 10);
    let mut project = fixture.project.clone();
    project.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
        name: "temperature-volume-mean-adjoint".into(), kind: "report".into(), region: None,
    });
    let nominal = native(&fixture, &project, 11);
    assert_eq!(field(&fixture, &original), field(&fixture, &nominal));
    let report = nominal.get("nominal_adjoint").unwrap();
    assert_eq!(report.str_field("functional"), Some("region-volume-mean-temperature"));
    assert_eq!(report.str_field("mode"), Some("enclosure-full-radiosity-feedback"));
    assert_eq!(report.get("selected_vertex"), Some(&JsonValue::Null));
    let mean = report.f64_field("value_k").unwrap();
    let temperatures = field(&fixture, &nominal);
    let values: Vec<_> = temperatures.as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    assert!(mean >= values.iter().copied().fold(f64::INFINITY, f64::min));
    assert!(mean <= values.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    let step = 0.02;
    let mut values = Vec::new();
    for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let mut perturbed = project.clone();
        perturbed.power.as_mut().unwrap()[0].watts.value += sign * step;
        let receipt = native(&fixture, &perturbed, 12 + side);
        values.push(receipt.get("nominal_adjoint").unwrap().f64_field("value_k").unwrap());
    }
    let expected = (values[1] - values[0]) / (2.0 * step);
    let actual = coefficient(&nominal, "power");
    assert!((actual - expected).abs() < 5e-5 * expected.abs().max(0.02),
        "volume-mean enclosure derivative {actual:e} versus native differences {expected:e}");
}

fn enclosure_study(controlled: bool) -> String {
    study_text(false, false, controlled)
        .replace(":materials (\"aa6061.fsmcdpk\")",
            ":materials (\"aa6061.fsmcdpk\" \"gray-surface.fsmcdpk\")")
        .replace(":entity \"air\" :low 4W :high 6W", ":entity \"body\" :low 4W :high 6W")
        .replace(":entity \"air\" :low 294K :high 300K",
            ":entity \"cold-patch\" :low 289K :high 291K")
}

#[test]
fn native_enclosure_nominal_control_reuses_one_complete_dual_without_changing_probability_children() {
    let raw = fixture();
    std::fs::write(raw.sources.join("study.fsim"), enclosure_study(false)).unwrap();
    let raw_result = raw.study(None, fs_cli::exit::SUCCESS);
    let (_, raw_report) = raw.report(&raw_result);
    let controlled = fixture();
    std::fs::write(controlled.sources.join("study.fsim"), enclosure_study(true)).unwrap();
    let result = controlled.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = controlled.report(&result);
    assert_eq!(rows(&report), rows(&raw_report), "calibration must not change any physical probability child");
    assert_eq!(report.get("statistics"), raw_report.get("statistics"));
    let control = report.get("mean_control").unwrap();
    assert_eq!(control.str_field("method"), Some("nominal-adjoint"));
    assert_eq!(control.f64_field("probe_solves_planned"), Some(1.0));
    assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
    let gradient = control.get("gradient").unwrap().as_array().unwrap();
    assert_eq!(gradient.len(), 2);
    assert!(gradient.iter().all(|v| v.as_f64().unwrap().is_finite()));
    assert!(gradient[0].as_f64().unwrap() > 0.0);
    assert_eq!(control.get("estimate").unwrap().f64_field("samples_in_estimate"), Some(4.0));
}

#[test]
fn native_enclosure_nominal_adjoint_refuses_unconnected_natural_feedback() {
    let fixture = fixture();
    let mut project = fixture.project.clone();
    project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.1, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(290.0, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    request(&mut project);
    let refused = execute(&fixture, &project, 20, fs_cli::exit::REFUSED);
    assert!(format!("{refused:?}").contains("enclosure"),
        "a missing enclosure/natural tangent must refuse, never fall back to frozen radiation");
}
