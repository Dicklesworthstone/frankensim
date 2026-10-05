//! Fresh native imports, material-card queries, nonlinear solves and retained fields.
use super::*;

fn project(natural: bool) -> fs_project::ProjectSpec {
    let mut project = fs_project::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-radiation.fsim"
    ))).unwrap().decoded.spec;
    project.envelope.as_mut().unwrap().ambient_lo.value = 290.0;
    project.envelope.as_mut().unwrap().ambient_hi.value = 330.0;
    project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    if natural {
        setup.boundaries[0].condition = fs_project::ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: fs_qty::QtyAny::new(0.06, fs_project::spec::dims::LENGTH),
            ambient_temperature: fs_qty::QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    } else if let fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. }
        = &mut setup.boundaries[0].condition { reference_temperature.value = 300.0; }
    let radiation = setup.radiation.as_mut().unwrap();
    // Different convective/radiative references expose the missing weighted
    // reference derivative. This reservoir HEATS the surface; no cooling-only
    // sign convention is allowed in the adjoint.
    radiation.surfaces[0].reservoir_temperature.value = 310.0;
    radiation.temperature_tolerance.value = 1e-10;
    radiation.heat_tolerance.value = 1e-10;
    project
}

fn run(fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize) -> JsonValue {
    let source = fixture.dir.join(format!("radiation-adjoint-{ordinal}.fsim"));
    std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
    let ledger = fixture.dir.join("radiation-adjoint.db");
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
    command(&["--json", "import", source.to_str().unwrap(),
        reference.join("plate.stl").to_str().unwrap(), ledger.to_str().unwrap(),
        "--unit", "m", "--max-hole-edges", "0"], fs_cli::exit::SUCCESS);
    let result = command(&["--json", "solve", source.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", reference.join("aa6061.fsmcdpk").to_str().unwrap(),
        "--materials", reference.join("gray-surface.fsmcdpk").to_str().unwrap()], fs_cli::exit::SUCCESS);
    let receipt = json_artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|s| s.str_field("stage") == Some("conduction")).unwrap();
    json_artifact(&ledger, stage.str_field("receipt").unwrap())
}

fn field(fixture: &Fixture, receipt: &JsonValue) -> JsonValue {
    json_artifact(&fixture.dir.join("radiation-adjoint.db"),
        receipt.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone()
}

fn perturb(project: &mut fs_project::ProjectSpec, parameter: &str, delta: f64) {
    if parameter == "power" { project.power.as_mut().unwrap()[0].watts.value += delta; return; }
    let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    if parameter == "radiation-reservoir-temperature" {
        setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += delta;
        return;
    }
    match &mut setup.boundaries[0].condition {
        fs_project::ThermalBoundaryCondition::Convection { coefficient, reference_temperature } => {
            match parameter {
                "convection-coefficient" => coefficient.value += delta,
                "convection-temperature" => reference_temperature.value += delta,
                _ => panic!("unknown parameter {parameter}"),
            }
        }
        fs_project::ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } => {
            assert_eq!(parameter, "natural-convection-ambient");
            ambient_temperature.value += delta;
        }
        _ => panic!("fixture's physical law"),
    }
}

#[test]
fn native_radiation_and_passive_cooling_adjoints_match_real_parameter_resolves() {
    for natural in [false, true] {
        let fixture = Fixture::new();
        let base = project(natural);
        let original = run(&fixture, &base, 0);
        let mut requested = base.clone();
        request(&mut requested);
        let nominal = run(&fixture, &requested, 1);
        assert_eq!(field(&fixture, &original), field(&fixture, &nominal));
        assert_eq!(original.get("energy"), nominal.get("energy"));
        let adjoint = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(adjoint.str_field("mode"), Some("radiation-full-wall-feedback"));
        assert_eq!(adjoint.str_field("authority"), Some("Estimated"));
        assert!(adjoint.f64_field("true_relative_residual").unwrap() < 1e-10);
        assert!(coefficient(&nominal, "radiation-emissivity").is_finite());
        let selected = adjoint.f64_field("selected_vertex").unwrap() as usize;
        let ambient = if natural { "natural-convection-ambient" } else { "convection-temperature" };
        let mut controls = vec![("power", 0.02), (ambient, 0.01), ("radiation-reservoir-temperature", 0.01)];
        if !natural { controls.push(("convection-coefficient", 0.001)); }
        for (case, (parameter, step)) in controls.into_iter().enumerate() {
            let mut values = Vec::new();
            for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
                let mut perturbed = base.clone();
                perturb(&mut perturbed, parameter, sign*step);
                let solved = run(&fixture, &perturbed, 2+2*case+side);
                values.push(field(&fixture, &solved).as_array().unwrap()[selected].as_f64().unwrap());
            }
            // Use the SAME nodal functional, not a possibly relocated maximum.
            let expected = (values[1]-values[0])/(2.0*step);
            let actual = coefficient(&nominal, parameter);
            assert!((actual-expected).abs() < 5e-5*expected.abs().max(0.05),
                "natural={natural}, {parameter}: adjoint {actual:e}, native re-solves {expected:e}");
        }
    }
}

#[test]
fn mean_feedback_rebinding_checks_constitutive_residual_even_for_zero_goal() {
    use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel,
        ConductivityTable, ThermalBc, ThermalBoundaryBuilder, ScalarField, LinearConfig};
    use fs_conduction::adjoint::RobinResponse;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let (complex, positions) = fs_conduction::fixtures::box_grid([1,1,1], [0.1,0.1,0.1]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 8.0), (450.0, 8.0)]).unwrap());
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("wall", ThermalBc::robin(10.0,300.0).unwrap()).unwrap().finish().unwrap();
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 17, kernel_id: 821, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let linear = LinearConfig { tolerance: 1e-9, ..LinearConfig::default() };
        let n = mesh.vertex_count();
        let call = |points: &[[f64;4]]| RobinResponse::pullback_mean_robin_at(&cx, problem, None,
            linear, &vec![300.0;n], &["wall"], points, &vec![0.0;n], 2*n);
        let zero = call(&[[10.0,300.0,0.2,0.1]]).unwrap();
        assert_eq!(zero.iterations, 0);
        assert!(zero.nodal_load.iter().all(|&v| v == 0.0));
        assert!(matches!(call(&[[10.0,310.0,0.2,0.1]]),
            Err(fs_conduction::ConductionError::LinearSolveFailed { .. })),
            "convergence of the old frozen row cannot stand in for the full law");
        assert!(call(&[]).is_err());
        assert!(call(&[[10.0,300.0,0.2,f64::NAN]]).is_err());
    });
}
