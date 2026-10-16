//! Full CLI import -> material-card resolution -> coupled physical solve ->
//! adjoint receipt. Differences compare the SAME selected native vertex.
use super::*;

fn fixture(natural: bool) -> Fixture {
    let mut fixture = Fixture::new();
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
    fixture.project = fs_project::parse_sexpr_migrating(
        &std::fs::read_to_string(reference.join("cooling-radiation.fsim")).unwrap())
        .unwrap().decoded.spec;
    std::fs::copy(reference.join("gray-surface.fsmcdpk"), fixture.sources.join("gray-surface.fsmcdpk")).unwrap();
    let envelope = fixture.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 280.0;
    envelope.ambient_hi.value = 320.0;
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let radiation = setup.radiation.as_mut().unwrap();
    radiation.temperature_tolerance.value = 1e-11;
    radiation.heat_tolerance.value = 1e-10;
    radiation.surfaces[0].reservoir_temperature.value = 296.0;
    if natural {
        setup.boundaries[0].condition = fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(293.15, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    }
    fixture
}
fn native(fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize) -> JsonValue {
    let source = fixture.dir.join(format!("radiation-adjoint-{ordinal}.fsim"));
    std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
    let ledger = fixture.dir.join("radiation-adjoint.db");
    command(&["--json", "import", source.to_str().unwrap(),
        fixture.sources.join("plate.stl").to_str().unwrap(), ledger.to_str().unwrap(),
        "--unit", "m", "--max-hole-edges", "0"], fs_cli::exit::SUCCESS);
    let run = command(&["--json", "solve", source.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", fixture.sources.join("aa6061.fsmcdpk").to_str().unwrap(),
        "--materials", fixture.sources.join("gray-surface.fsmcdpk").to_str().unwrap()], fs_cli::exit::SUCCESS);
    let receipt = json_artifact(&ledger, run.str_field("run_receipt").unwrap());
    let stage = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|s| s.str_field("stage") == Some("conduction")).unwrap();
    json_artifact(&ledger, stage.str_field("receipt").unwrap())
}
fn field(fixture: &Fixture, receipt: &JsonValue) -> JsonValue {
    json_artifact(&fixture.dir.join("radiation-adjoint.db"),
        receipt.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone()
}

#[test]
fn native_radiation_and_natural_radiation_adjoint_preserve_fields_and_match_resolves() {
    for natural in [false, true] {
        let fixture = fixture(natural);
        let mut project = fixture.project.clone();
        project.power.as_mut().unwrap()[0].duty = 0.37;
        let baseline = native(&fixture, &project, 0);
        let mut requested = project.clone();
        request(&mut requested);
        let nominal = native(&fixture, &requested, 1);
        assert_eq!(field(&fixture, &baseline), field(&fixture, &nominal));
        assert_eq!(baseline.get("energy"), nominal.get("energy"));
        let adjoint = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(adjoint.str_field("mode"), Some("radiation-full-wall-feedback"));
        assert_eq!(adjoint.str_field("authority"), Some("Estimated"));
        let selected = adjoint.f64_field("selected_vertex").unwrap() as usize;
        let ambient_target = if natural { "natural-convection-ambient" } else { "convection-temperature" };
        for (case, target, step) in [(0, "power", 0.01),
            (1, "radiation-reservoir-temperature", 0.01), (2, ambient_target, 0.01)] {
            let derivative = coefficient(&nominal, target);
            let mut results = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut perturbed = project.clone();
                if case == 0 {
                    perturbed.power.as_mut().unwrap()[0].watts.value += sign*step;
                } else {
                    let setup = perturbed.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
                    if case == 1 {
                        setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += sign*step;
                    } else {
                        match &mut setup.boundaries[0].condition {
                            fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. } =>
                                reference_temperature.value += sign*step,
                            fs_project::ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } =>
                                ambient_temperature.value += sign*step,
                            _ => panic!("fixture boundary"),
                        }
                    }
                }
                let receipt = native(&fixture, &perturbed, 2+case*2+side);
                results.push(field(&fixture, &receipt).as_array().unwrap()[selected].as_f64().unwrap());
            }
            let difference = (results[1]-results[0])/(2.0*step);
            assert!((derivative-difference).abs() < 2e-5*derivative.abs().max(0.01),
                "natural={natural}, {target}: {derivative:e} versus {difference:e}");
        }
        assert!(coefficient(&nominal, "radiation-emissivity").is_finite());
    }
}
