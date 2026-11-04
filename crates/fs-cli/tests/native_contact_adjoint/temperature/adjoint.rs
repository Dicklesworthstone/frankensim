//! G3: imported source-backed contacts, complete adjoints and physical resolves.
use super::*;
use fs_project::{ThermalBoundaryCondition as B, spec::dims};
use fs_qty::QtyAny;

#[path = "adjoint/natural.rs"]
mod natural;

fn scaled_pack(scale: f64, surface: bool) -> (Vec<u8>, String, String) {
    if !surface {
        return temperature_pack_from_knots(vec![(200.0, 0.2 * scale), (450.0, 0.45 * scale)]);
    }
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: vec![([200.0, 0.8e6], 0.2 * scale), ([450.0, 0.8e6], 0.45 * scale),
            ([200.0, 1.2e6], 0.2 * scale), ([450.0, 1.2e6], 0.45 * scale)],
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

/// Author independently scaled source cards, then run the ordinary native
/// import and solve. No adjoint data or previously solved field is injected.
fn solve(f: &Fixture, scales: &[f64], surface: bool, output: Option<&str>, ordinal: usize) -> (J, J) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut project = f.project.clone();
    let bindings = project.interface_cards.as_mut().unwrap();
    assert_eq!(bindings.len(), scales.len());
    let mut cards = Vec::new();
    for (i, (binding, &scale)) in bindings.iter_mut().zip(scales).enumerate() {
        let (bytes, card, claim) = scaled_pack(scale, surface);
        binding.card = card;
        binding.claim = Some(claim);
        let path = f.dir.join(format!("adjoint-contact-{ordinal}-{i}.fsintpk"));
        std::fs::write(&path, bytes).unwrap();
        cards.push(path);
    }
    if let Some(name) = output {
        project.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let source = f.dir.join(format!("adjoint-project-{ordinal}.fsim"));
    let ledger = f.dir.join("source-contact-adjoint.db");
    std::fs::write(&source, fs_project::print_sexpr(&project).unwrap()).unwrap();
    command(&["--json", "import", source.to_str().unwrap(),
        root.join("examples/contact-pair/cold-body.stl").to_str().unwrap(),
        root.join("examples/contact-pair/hot-body.stl").to_str().unwrap(),
        ledger.to_str().unwrap(), "--unit", "m", "--max-hole-edges", "0"]);
    let mut args = vec!["--json".to_string(), "solve".into(), source.display().to_string(),
        ledger.display().to_string(), "--materials".into(), f.dir.join("solid.fsmcdpk").display().to_string()];
    for card in cards { args.extend(["--interfaces".into(), card.display().to_string()]); }
    let result = fs_cli::run(args);
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", result.stdout, result.stderr);
    let result = J::parse(&result.stdout).unwrap();
    let run = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = run.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

fn value(field: &J, vertex: usize) -> f64 {
    field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap()
}

fn derivative(report: &J, target: &str, entity: &str) -> f64 {
    let rows: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|r| r.str_field("target") == Some(target) && r.str_field("entity") == Some(entity)).collect();
    assert_eq!(rows.len(), 1, "each physical control has one declared owner");
    rows[0].f64_field("derivative").unwrap()
}

fn agrees(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 5e-5 * expected.abs().max(1e-3),
        "complete source-backed adjoint {actual:e} versus independent physical resolve {expected:e}");
}

#[test]
fn live_contact_and_power_controls_match_physical_resolves_without_changing_the_primal() {
    for (nonlinear, surface) in [(false, false), (true, false), (true, true)] {
        let mut f = Fixture::new(nonlinear, false);
        let (plain, plain_field) = solve(&f, &[1.0], surface, None, 0);
        let (nominal, field) = solve(&f, &[1.0], surface, Some(OUTPUT), 1);
        assert_eq!(plain_field.get("temperature"), field.get("temperature"));
        for key in ["energy", "interfaces", "temperature_contact", "solver"] {
            assert_eq!(plain.get(key), nominal.get(key), "{key}");
        }
        assert_eq!(nominal.get("temperature_contact").unwrap().str_field("algorithm"),
            Some(if surface { "bounded-contact-picard" } else { "coupled-newton" }));
        let report = nominal.get("nominal_adjoint").unwrap();
        assert_eq!(report.str_field("mode"), Some("nonlinear-contact-full-material-feedback"));
        assert_eq!(report.str_field("authority"), Some("Estimated"));
        assert!(report.f64_field("true_relative_residual").unwrap() <= 1e-12);
        assert_eq!(report.get("stability_iterations"), Some(&J::Null));
        let selected = report.f64_field("selected_vertex").unwrap() as usize;
        let contact = derivative(report, TARGET, "cold-hot-joint");
        let power = derivative(report, "power", "hot");
        assert!(contact > 0.0 && power > 0.0);
        let delta = 0.002;
        let low = solve(&f, &[1.0 - delta], surface, None, 2).1;
        let high = solve(&f, &[1.0 + delta], surface, None, 3).1;
        agrees(contact, (value(&high, selected) - value(&low, selected)) / (2.0 * delta));
        f.project.power.as_mut().unwrap()[0].watts.value = 12.0 - 0.01;
        let low = solve(&f, &[1.0], surface, None, 4).1;
        f.project.power.as_mut().unwrap()[0].watts.value = 12.0 + 0.01;
        let high = solve(&f, &[1.0], surface, None, 5).1;
        agrees(power, (value(&high, selected) - value(&low, selected)) / 0.02);
        f.project.power.as_mut().unwrap()[0].watts.value = 12.0;
        let (legacy, legacy_field) = solve(&f, &[1.0], surface, Some("temperature-max-adjoint"), 6);
        assert_eq!(legacy_field.get("temperature"), field.get("temperature"));
        let legacy = legacy.get("nominal_adjoint").unwrap();
        let retained: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
            .filter(|row| row.str_field("target") != Some(TARGET)).collect();
        assert_eq!(retained, legacy.get("parameters").unwrap().as_array().unwrap().iter().collect::<Vec<_>>());
        for key in ["selected_vertex", "value_k", "mode", "dual_iterations", "true_relative_residual"] {
            assert_eq!(report.get(key), legacy.get(key), "one unchanged complete dual: {key}");
        }
    }
}

#[test]
fn live_contact_robin_controls_keep_the_full_thermal_feedback() {
    let mut f = Fixture::new(true, false);
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition = B::Convection {
        coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
    };
    let (nominal, _) = solve(&f, &[1.0], false, Some(OUTPUT), 20);
    let report = nominal.get("nominal_adjoint").unwrap();
    let selected = report.f64_field("selected_vertex").unwrap() as usize;
    for (case, (target, step)) in [("convection-coefficient", 0.01), ("convection-temperature", 0.002)]
        .into_iter().enumerate()
    {
        let actual = derivative(report, target, "cold");
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let B::Convection { coefficient, reference_temperature } = &mut f.project.cooling.as_mut()
                .unwrap().conduction.as_mut().unwrap().boundaries[0].condition else { panic!("Robin law"); };
            coefficient.value = 8.0 + if case == 0 { sign * step } else { 0.0 };
            reference_temperature.value = 293.15 + if case == 1 { sign * step } else { 0.0 };
            values.push(value(&solve(&f, &[1.0], false, None, 21 + case * 2 + side).1, selected));
        }
        agrees(actual, (values[1] - values[0]) / (2.0 * step));
    }
}

#[test]
fn named_source_controls_remain_independent_when_two_interfaces_share_one_card() {
    let mut f = Fixture::new(true, false);
    let second = "other-joint";
    f.project.assembly.as_mut().unwrap().push(fs_project::spec::EntityDecl::Interface {
        parent: "assembly".into(), name: second.into(), display: "Other half of the joint".into(),
        from: "cold".into(), to: "hot".into(), expect_id: None,
    });
    let assignments = f.project.assignments.as_mut().unwrap();
    let first = assignments.iter_mut().find(|row| row.target == "cold-hot-joint").unwrap();
    first.selector = fs_io::MeshSelector::NearestDatum {
        point: [1.0, 1.0 / 3.0, 2.0 / 3.0], max_distance: 1e-10, tolerance: 0.0,
    };
    let mut other = first.clone();
    other.target = second.into();
    other.selector = fs_io::MeshSelector::NearestDatum {
        point: [1.0, 2.0 / 3.0, 1.0 / 3.0], max_distance: 1e-10, tolerance: 0.0,
    };
    assignments.push(other);
    let bindings = f.project.interface_cards.as_mut().unwrap();
    let mut other = bindings[0].clone();
    other.interface = second.into();
    bindings.push(other);
    let (nominal, _) = solve(&f, &[1.0, 1.0], false, Some(OUTPUT), 30);
    let report = nominal.get("nominal_adjoint").unwrap();
    let selected = report.f64_field("selected_vertex").unwrap() as usize;
    let rows: Vec<_> = report.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|row| row.str_field("target") == Some(TARGET)).collect();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].str_field("interface_card"), rows[1].str_field("interface_card"));
    assert_ne!(rows[0].str_field("entity"), rows[1].str_field("entity"));
    for row in &rows { assert_eq!(row.f64_field("face_pairs"), Some(1.0)); }
    for (control, name) in ["cold-hot-joint", second].into_iter().enumerate() {
        let delta = 0.002;
        let mut low_scale = [1.0, 1.0]; low_scale[control] -= delta;
        let mut high_scale = [1.0, 1.0]; high_scale[control] += delta;
        let low = solve(&f, &low_scale, false, None, 31 + 2 * control).1;
        let high = solve(&f, &high_scale, false, None, 32 + 2 * control).1;
        let actual = derivative(report, TARGET, name);
        assert!(actual > 0.0, "each real half-interface conducts heat from the hot body");
        agrees(actual, (value(&high, selected) - value(&low, selected)) / (2.0 * delta));
    }
}

#[test]
fn source_corners_do_not_inherit_the_primal_generalized_tangent_as_an_adjoint() {
    let mut f = Fixture::new(false, false);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition = B::FixedTemperature {
        temperature: QtyAny::new(300.0, dims::TEMPERATURE),
    };
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 300.0; envelope.ambient_hi.value = 300.0;
    let kink = || temperature_pack_from_knots(vec![(200.0, 0.2), (300.0, 0.25), (450.0, 0.55)]);
    let (source, card, ledger) = prepare_with_pack(&f, 40, None, kink());
    let primal = run(&f, &source, &card, &ledger);
    assert_eq!(primal.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", primal.stdout, primal.stderr);
    let (source, card, ledger) = prepare_with_pack(&f, 41, Some(OUTPUT), kink());
    let refused = run(&f, &source, &card, &ledger);
    assert_ne!(refused.exit_code, fs_cli::exit::SUCCESS);
    let diagnostic = format!("{}{}", refused.stdout, refused.stderr);
    assert!(diagnostic.contains("derivative") || diagnostic.contains("slope"), "{diagnostic}");
}

#[test]
fn source_contact_boundary_lift_reports_stay_explicitly_unsupported() {
    let f = Fixture::new(false, false);
    for (ordinal, output) in ["temperature-max-boundary-adjoint", "temperature-max-contact-boundary-adjoint",
        "temperature-volume-mean-adjoint"].into_iter().enumerate()
    {
        let (source, card, ledger) = prepare(&f, 50 + ordinal, Some(output));
        let refused = run(&f, &source, &card, &ledger);
        assert_ne!(refused.exit_code, fs_cli::exit::SUCCESS);
        assert!(format!("{}{}", refused.stdout, refused.stderr).contains("cli-solve-conduction-contact-composition"));
    }
}
