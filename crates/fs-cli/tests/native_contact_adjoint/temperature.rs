//! A real imported project must use each contact's current source temperature.
use super::*;

#[path = "temperature/transient.rs"]
mod transient;
#[path = "temperature/radiation.rs"]
mod radiation;
#[path = "temperature/enthalpy.rs"]
mod enthalpy;

fn temperature_pack() -> (Vec<u8>, String, String) {
    temperature_pack_from_knots(vec![(200.0, 0.2), (450.0, 0.45)])
}

fn temperature_pack_from_knots(knots: Vec<(f64, f64)>) -> (Vec<u8>, String, String) {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let claims = pack(
        "synthetic-temperature-contact",
        fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        dims,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_project::spec::dims::TEMPERATURE,
            knots,
            dims,
        },
        InterpolationPolicy::LinearInside,
    );
    let state = |chemistry: &str| MaterialStateId {
        chemistry: chemistry.into(),
        phase: "solid".into(),
        process: "fixture".into(),
        revision: 0,
    };
    let bytes = NormalizedInterfacePack::new(
        SurfaceSpec {
            material: state("cold-body"),
            texture_frame: "normal-plus-x".into(),
        },
        SurfaceSpec {
            material: state("hot-body"),
            texture_frame: "normal-minus-x".into(),
        },
        SystemContext {
            medium: "dry-contact".into(),
            third_body: None,
            environment: "fixture-air".into(),
            history: "unaged".into(),
        },
        claims,
    )
    .unwrap()
    .to_bytes();
    let set = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Interface,
        source: "synthetic-temperature-contact".into(),
        bytes: bytes.clone(),
        expect: None,
    }])
    .unwrap();
    let card = set.interfaces()[0].card().to_hex();
    let library = set.library();
    let selected = library
        .interface(&card)
        .unwrap()
        .claims()
        .claims_ordered()
        .next()
        .unwrap()
        .0
        .0
        .to_hex();
    (bytes, card, selected)
}

fn prepare(f: &Fixture, ordinal: usize, output: Option<&str>) -> (PathBuf, PathBuf, PathBuf) {
    prepare_with_pack(f, ordinal, output, temperature_pack())
}

fn prepare_with_pack(
    f: &Fixture, ordinal: usize, output: Option<&str>,
    (bytes, card, selected): (Vec<u8>, String, String),
) -> (PathBuf, PathBuf, PathBuf) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut project = f.project.clone();
    let binding = &mut project.interface_cards.as_mut().unwrap()[0];
    binding.card = card;
    binding.claim = Some(selected);
    if let Some(name) = output {
        project
            .outputs
            .get_or_insert_with(Vec::new)
            .push(fs_project::spec::OutputRequest {
                name: name.into(),
                kind: "report".into(),
                region: None,
            });
    }
    let source = f.dir.join(format!("temperature-{ordinal}.fsim"));
    let contact = f.dir.join("temperature.fsintpk");
    let ledger = f.dir.join("temperature.db");
    std::fs::write(&source, fs_project::print_sexpr(&project).unwrap()).unwrap();
    std::fs::write(&contact, bytes).unwrap();
    command(&[
        "--json",
        "import",
        source.to_str().unwrap(),
        root.join("examples/contact-pair/cold-body.stl")
            .to_str()
            .unwrap(),
        root.join("examples/contact-pair/hot-body.stl")
            .to_str()
            .unwrap(),
        ledger.to_str().unwrap(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]);
    (source, contact, ledger)
}

fn run(f: &Fixture, source: &Path, contact: &Path, ledger: &Path) -> fs_cli::CommandOutput {
    let material = f.dir.join("solid.fsmcdpk");
    let args = [
        "--json",
        "solve",
        source.to_str().unwrap(),
        ledger.to_str().unwrap(),
        "--materials",
        material.to_str().unwrap(),
        "--interfaces",
        contact.to_str().unwrap(),
    ];
    fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
}

#[test]
fn native_temperature_contact_queries_the_accepted_faces_and_preserves_heat() {
    for nonlinear in [false, true] {
        for boundary in 0..3 {
            let airflow = boundary == 1;
            let mut f = Fixture::new(nonlinear, airflow);
            f.project
                .cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .radiation = None;
            if boundary == 2 {
                f.project
                    .cooling
                    .as_mut()
                    .unwrap()
                    .conduction
                    .as_mut()
                    .unwrap()
                    .boundaries[0]
                    .condition = fs_project::ThermalBoundaryCondition::NaturalConvection {
                    characteristic_length: fs_qty::QtyAny::new(1.0, fs_project::spec::dims::LENGTH),
                    ambient_temperature: fs_qty::QtyAny::new(
                        300.0,
                        fs_project::spec::dims::TEMPERATURE,
                    ),
                    correlation: "convection.churchill-chu-vertical-plate".into(),
                };
            }
            let (source, contact, ledger) = prepare(&f, 0, None);
            let result = run(&f, &source, &contact, &ledger);
            assert_eq!(
                result.exit_code,
                fs_cli::exit::SUCCESS,
                "{}\n{}",
                result.stdout,
                result.stderr
            );
            let result = J::parse(&result.stdout).unwrap();
            let run_receipt = artifact(&ledger, result.str_field("run_receipt").unwrap());
            let stage = run_receipt
                .get("stages")
                .unwrap()
                .as_array()
                .unwrap()
                .iter()
                .find(|row| row.str_field("stage") == Some("conduction"))
                .unwrap();
            let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
            let physical = receipt.get("temperature_contact").unwrap();
            assert_eq!(physical.str_field("algorithm"), Some("coupled-newton"));
            assert_eq!(physical.f64_field("thermal_responses"), Some(1.0));
            assert!(
                physical.f64_field("total_thermal_responses").unwrap()
                    >= physical.f64_field("thermal_responses").unwrap()
            );
            let history = physical.get("newton_history").unwrap().as_array().unwrap();
            assert_eq!(
                physical.f64_field("nonlinear_iterations"),
                Some(history.len() as f64)
            );
            assert_eq!(
                physical.f64_field("jacobian_assemblies"),
                Some(history.len() as f64)
            );
            assert!(
                physical.f64_field("residual_evaluations").unwrap()
                    >= history.len() as f64 + 1.0
            );
            assert_eq!(
                physical.f64_field("krylov_iterations").unwrap(),
                history.iter().map(|row| row.f64_field("linear_iterations").unwrap()).sum::<f64>()
            );
            for row in history {
                let step = row.f64_field("accepted_step_length").unwrap();
                assert!(step > 0.0 && step <= 1.0);
                assert!(row.f64_field("residual_after_w").unwrap()
                    < row.f64_field("residual_before_w").unwrap());
                assert!(row.f64_field("linear_relative_residual").unwrap()
                    <= row.f64_field("forcing").unwrap());
            }
            assert!(
                physical.f64_field("residual_norm_w").unwrap()
                    <= physical.f64_field("residual_tolerance_w").unwrap()
            );
            assert!(
                physical.f64_field("energy_closure_w").unwrap().abs()
                    <= physical.f64_field("energy_tolerance_w").unwrap()
            );
            assert_eq!(
                receipt.get("solver").unwrap().f64_field("final_residual"),
                physical.f64_field("residual_norm_w")
            );
            assert_eq!(
                receipt.get("energy").unwrap().f64_field("closure_w"),
                physical.f64_field("energy_closure_w")
            );
            let states = physical.get("face_states").unwrap().as_array().unwrap();
            assert!(!states.is_empty());
            for row in states {
                let query = row.f64_field("query_temperature_k").unwrap();
                let a = row.f64_field("side_a_temperature_k").unwrap();
                let b = row.f64_field("side_b_temperature_k").unwrap();
                assert_eq!(query, 0.5 * a + 0.5 * b);
                assert!(
                    (row.f64_field("resistance_m2_k_w").unwrap() - 0.001 * query).abs() < 1e-14
                );
                assert_eq!(
                    row.str_field("selected_claim"),
                    Some(temperature_pack().2.as_str())
                );
                assert_eq!(row.str_field("usage_receipt").unwrap().len(), 64);
            }
            let interface_heat = receipt
                .get("interfaces")
                .unwrap()
                .get("fluxes")
                .unwrap()
                .as_array()
                .unwrap()
                .iter()
                .map(|row| row.f64_field("heat_rate_a_to_b_w").unwrap())
                .sum::<f64>();
            assert!((interface_heat.abs() - 12.0 * 0.37).abs() < 1e-7);
            assert_eq!(
                receipt.get("solver_control"),
                Some(&J::Null),
                "a frozen-contact inverse is not a nonlinear certificate"
            );
            let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
            let maximum = field
                .get("temperature")
                .unwrap()
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .fold(f64::NEG_INFINITY, f64::max);
            let low = states
                .iter()
                .map(|row| row.f64_field("resistance_m2_k_w").unwrap())
                .fold(f64::INFINITY, f64::min);
            let high = states
                .iter()
                .map(|row| row.f64_field("resistance_m2_k_w").unwrap())
                .fold(f64::NEG_INFINITY, f64::max);
            let max_fixed = |resistance, ordinal| {
                f.solve(resistance, None, ordinal)
                    .1
                    .get("temperature")
                    .unwrap()
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_f64().unwrap())
                    .fold(f64::NEG_INFINITY, f64::max)
            };
            assert!(maximum >= max_fixed(low * 0.99, 1) - 1e-7);
            assert!(maximum <= max_fixed(high * 1.01, 2) + 1e-7);
            if airflow {
                assert!(
                    receipt.get("conjugate").unwrap().get("outcome").is_some()
                        || receipt
                            .get("conjugate")
                            .unwrap()
                            .get("iterations")
                            .is_some()
                );
            }
        }
    }
}

#[test]
fn temperature_contact_is_rebound_on_each_native_uniform_rung() {
    let mut f = Fixture::new(false, false);
    f.project.solver.as_mut().unwrap().fidelity = "ladder".into();
    f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    let (source, contact, ledger) = prepare(&f, 30, None);
    let result = run(&f, &source, &contact, &ledger);
    assert_eq!(
        result.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        result.stdout,
        result.stderr
    );
    let result = J::parse(&result.stdout).unwrap();
    let run_receipt = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = run_receipt
        .get("stages")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("stage") == Some("conduction"))
        .unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let rungs = receipt
        .get("ladder")
        .unwrap()
        .get("rungs")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(rungs.len(), 3);
    let states = receipt
        .get("temperature_contact")
        .unwrap()
        .get("face_states")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(
        states.len() as f64,
        receipt
            .get("interfaces")
            .unwrap()
            .f64_field("pair_count")
            .unwrap()
    );
    assert!(
        states.len() > 2,
        "the published contact must use refined faces, not the imported coarse pair"
    );
}

#[test]
fn temperature_contact_never_substitutes_frozen_derivatives_or_adaptive_goals() {
    let f = Fixture::new(false, false);
    let (source, contact, ledger) = prepare(&f, 20, Some(OUTPUT));
    let result = run(&f, &source, &contact, &ledger);
    assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
    assert!(
        format!("{}{}", result.stdout, result.stderr)
            .contains("cli-solve-conduction-contact-composition")
    );
    let mut f = Fixture::new(false, true);
    f.project.solver.as_mut().unwrap().fidelity = "adaptive".into();
    let (source, contact, ledger) = prepare(&f, 21, None);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let result = fs_cli::run(vec![
        "--json".into(),
        "solve".into(),
        source.display().to_string(),
        ledger.display().to_string(),
        "--materials".into(),
        f.dir.join("solid.fsmcdpk").display().to_string(),
        "--materials".into(),
        root.join("data/reference-project/gray-surface.fsmcdpk")
            .display()
            .to_string(),
        "--interfaces".into(),
        contact.display().to_string(),
    ]);
    assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
    assert!(
        format!("{}{}", result.stdout, result.stderr)
            .contains("cli-solve-conduction-contact-composition")
    );
}

#[path = "temperature/transient_storage.rs"]
mod transient_storage;
