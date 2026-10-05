//! Actual native import/material/solve/receipt comparisons, not a solver stub.
use super::*;

#[path = "nominal_adjoint/natural_feedback.rs"]
mod natural_feedback;
#[path = "nominal_adjoint/radiative_feedback.rs"]
mod radiative_feedback;
#[path = "nominal_adjoint/radiation_owner.rs"]
mod radiation_owner;

fn request(project: &mut fs_project::ProjectSpec) {
    project.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
        name: "temperature-max-adjoint".into(), kind: "report".into(),
        region: None,
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

fn study_text(copula: bool, qmc: bool, controlled: bool) -> String {
    let mut text = STUDY.to_string();
    if qmc { text = text.replace("monte-carlo", "quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)"); }
    if copula { text = text.replace("independent", "(gaussian-copula :latent-correlation ((1 0.5) (0.5 1)))"); }
    if controlled { text = text.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)"); }
    text
}

#[test]
fn nominal_adjoint_mean_control_uses_one_native_calibration_for_all_four_samplers() {
    for (copula, qmc) in [(false,false),(true,false),(false,true),(true,true)] {
        let raw = Fixture::new();
        std::fs::write(raw.sources.join("study.fsim"),study_text(copula,qmc,false)).unwrap();
        let raw_result = raw.study(None,fs_cli::exit::SUCCESS);
        let (_, raw_report) = raw.report(&raw_result);
        let fixture = Fixture::new();
        std::fs::write(fixture.sources.join("study.fsim"),study_text(copula,qmc,true)).unwrap();
        let result = fixture.study(None,fs_cli::exit::SUCCESS);
        let (_, report) = fixture.report(&result);
        assert_eq!(rows(&report),rows(&raw_report),"control must not change any raw physical solve");
        assert_eq!(report.get("statistics"),raw_report.get("statistics"));
        assert_eq!(report.get("qmc"),raw_report.get("qmc"));
        let control = report.get("mean_control").unwrap();
        assert_eq!(control.str_field("method"),Some("nominal-adjoint"));
        assert_eq!(control.f64_field("probe_solves_planned"),Some(1.0));
        assert_eq!(control.f64_field("probe_solves_completed"),Some(1.0));
        let estimate=control.get("estimate").unwrap();
        assert_eq!(estimate.f64_field("samples_in_estimate"),Some(4.0));
        // The reference model is affine in power and ambient, including its
        // spatial field. A correct nominal control removes its sample variation.
        assert!(estimate.f64_field("variance_ratio").unwrap() < 1e-4,"{estimate:?}");
    }
}

#[test]
fn nominal_adjoint_checkpoint_retains_calibration_and_resumes_without_original_files() {
    let full=Fixture::new();
    let text=study_text(true,true,true);
    std::fs::write(full.sources.join("study.fsim"),&text).unwrap();
    let complete=full.study(None,fs_cli::exit::SUCCESS);
    let (expected, _)=full.report(&complete);
    let fixture=Fixture::new();
    std::fs::write(fixture.sources.join("study.fsim"),text).unwrap();
    let prefix=fixture.study(Some("1"),fs_cli::exit::BUDGET);
    let (_, prefix_report)=fixture.report(&prefix);
    assert!(rows(&prefix_report).is_empty(),"the one allowed solve is calibration, not a random observation");
    let frozen=prefix_report.get("mean_control").unwrap().get("calibration").unwrap().clone();
    std::fs::rename(&fixture.sources,fixture.dir.join("relocated-sources")).unwrap();
    let one=fixture.resume(prefix.str_field("run").unwrap(),Some("1"),fs_cli::exit::BUDGET);
    let (_, one_report)=fixture.report(&one);
    assert_eq!(rows(&one_report).len(),1);
    assert_eq!(one_report.get("mean_control").unwrap().get("calibration"),Some(&frozen));
    let final_result=fixture.resume(one.str_field("run").unwrap(),None,fs_cli::exit::SUCCESS);
    let (actual, _)=fixture.report(&final_result);
    assert_eq!(actual,expected,"split quadrature must preserve the exact coefficient/sample history");
    assert_eq!(fixture.resume(final_result.str_field("run").unwrap(),None,fs_cli::exit::SUCCESS),final_result);
}

#[test]
fn nominal_adjoint_zero_width_inputs_need_no_calibration_and_no_fake_observations() {
    let fixture=Fixture::new();
    let text=study_text(false,false,true).replace(":max-solves 1",":max-solves 0")
        .replace(":high 6W",":high 4W").replace(":high 300K",":high 294K");
    std::fs::write(fixture.sources.join("study.fsim"),text).unwrap();
    let empty=fixture.study(Some("0"),fs_cli::exit::BUDGET);
    let (_, report)=fixture.report(&empty);
    assert!(rows(&report).is_empty());
    let control=report.get("mean_control").unwrap();
    assert_eq!(control.f64_field("probe_solves_completed"),Some(0.0));
    assert_eq!(control.get("gradient").unwrap().as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect::<Vec<_>>(),[0.0,0.0]);
    assert!(control.get("calibration").unwrap().get("probes").unwrap().as_array().unwrap().is_empty());
    let result=fixture.resume(empty.str_field("run").unwrap(),None,fs_cli::exit::SUCCESS);
    let (_, report)=fixture.report(&result);
    assert_eq!(rows(&report).len(),4);
}

#[test]
fn unsupported_nominal_fan_derivative_is_refused_before_ledger_creation() {
    let fixture=Fixture::new();
    let text=study_text(false,false,true).replace(
        "(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)",
        "(uniform :name \"fan\" :target fan-speed-ratio :entity \"fixture-bank\" :low 0.8 :high 1.2)");
    assert!(text.contains(":target fan-speed-ratio"));
    std::fs::write(fixture.sources.join("study.fsim"),text).unwrap();
    fixture.study(None,fs_cli::exit::REFUSED);
    assert!(!fixture.ledger.exists());
}

#[test]
fn native_surface_and_volume_power_adjoints_retain_distinct_loads_and_duties() {
    let fixture = Fixture::new();
    let mut base = fixture.project.clone();
    base.solver.as_mut().unwrap().tolerance_rel = 1e-9;
    base.assembly.as_mut().unwrap().push(fs_project::spec::EntityDecl::Surface {
        parent: "enclosure".into(), name: "chip".into(),
        display: "Heated chip footprint".into(), expect_id: None,
    });
    let assignments = base.assignments.as_mut().unwrap();
    assignments[0].allow_overlap = true;
    let mut patch = assignments[0].clone();
    patch.target = "chip".into();
    // One complete tetrahedron face, not a centroid selection or replacement
    // tet. The primal gives this Surface precedence over the surrounding air.
    patch.selector = fs_io::MeshSelector::HalfSpace {
        normal: [0.0,0.0,1.0], offset: 0.0,
        side: fs_io::HalfSpaceSide::AtMost, tolerance: 0.0,
    };
    assignments.push(patch);
    let power = base.power.as_mut().unwrap();
    power[0].duty = 0.61;
    let mut chip = power[0].clone();
    chip.region = "chip".into(); chip.watts.value = 7.0; chip.duty = 0.37;
    power.push(chip);
    let original = solve(&fixture, &base, 100);
    let mut requested = base.clone();
    request(&mut requested);
    let nominal = solve(&fixture, &requested, 101);
    let field = |receipt: &JsonValue| {
        json_artifact(&fixture.dir.join("adjoint.db"),
            receipt.str_field("solution_artifact").unwrap()).get("temperature").unwrap().clone()
    };
    assert_eq!(field(&original), field(&nominal));
    assert_eq!(original.get("energy"), nominal.get("energy"));
    let derivative = |receipt: &JsonValue, entity: &str| {
        receipt.get("nominal_adjoint").unwrap().get("parameters").unwrap().as_array().unwrap()
            .iter().find(|row| row.str_field("target") == Some("power")
                && row.str_field("entity") == Some(entity)).unwrap().f64_field("derivative").unwrap()
    };
    let selected = nominal.get("nominal_adjoint").unwrap().f64_field("selected_vertex").unwrap() as usize;
    for (case, entity) in ["air", "chip"].into_iter().enumerate() {
        let actual = derivative(&nominal, entity);
        let mut values = Vec::new();
        let step = 0.02;
        for (side, sign) in [-1.0,1.0].into_iter().enumerate() {
            let mut perturbed = base.clone();
            perturbed.power.as_mut().unwrap()[case].watts.value += sign*step;
            let receipt = solve(&fixture, &perturbed, 102+2*case+side);
            values.push(field(&receipt).as_array().unwrap()[selected].as_f64().unwrap());
        }
        let expected = (values[1]-values[0])/(2.0*step);
        assert!(actual > 0.0, "positive {entity} watts must inject heat, not remove it");
        assert!((actual-expected).abs() < 3e-6*actual.abs().max(0.01),
            "{entity}: adjoint {actual:e}, actual load re-solves {expected:e}");
    }
    assert!((derivative(&nominal, "chip")/0.37-derivative(&nominal, "air")/0.61).abs() > 1e-4,
        "the fixture must distinguish patch heat from volumetric smearing");
    requested.power.as_mut().unwrap()[1].duty = 0.0;
    let zero = solve(&fixture, &requested, 110);
    assert_eq!(derivative(&zero, "chip"), 0.0, "declared pre-duty watts have zero influence at zero duty");
}
