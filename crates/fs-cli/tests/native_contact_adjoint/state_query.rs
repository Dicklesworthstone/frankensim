//! One immutable pressure-response card, compared with independently fixed R''.
use super::*;

fn pressure_pack() -> (Vec<u8>, String) {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    // Interior pressure knot 300 deliberately overlaps the numerical thermal
    // range. It must not be mistaken for a temperature-dependent resistance.
    let claims = pack("synthetic-pressure-response", fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        dims, PropertyValue::Curve { abscissa: "normal_pressure".into(),
            abscissa_dims: fs_project::spec::dims::PRESSURE,
            knots: vec![(200.0, 0.20), (300.0, 0.15), (400.0, 0.10)], dims },
        InterpolationPolicy::LinearInside);
    let state = |chemistry: &str| MaterialStateId { chemistry: chemistry.into(),
        phase: "solid".into(), process: "fixture".into(), revision: 0 };
    let bytes = NormalizedInterfacePack::new(
        SurfaceSpec { material: state("cold-body"), texture_frame: "normal-plus-x".into() },
        SurfaceSpec { material: state("hot-body"), texture_frame: "normal-minus-x".into() },
        SystemContext { medium: "dry-contact".into(), third_body: None,
            environment: "fixture-air".into(), history: "unaged".into() }, claims).unwrap().to_bytes();
    let set = CardPackSet::admit(vec![RawCardPack { kind: CardPackKind::Interface,
        source: "synthetic-pressure-response".into(), bytes: bytes.clone(), expect: None }]).unwrap();
    (bytes, set.interfaces()[0].card().to_hex())
}

fn state_project(f: &Fixture, pressure: f64, width: f64) -> (fs_project::ProjectSpec, CardPackSet) {
    let mut project = f.project.clone();
    let (bytes, card) = pressure_pack();
    let binding = &mut project.interface_cards.as_mut().unwrap()[0];
    binding.card = card;
    binding.claim = None;
    binding.state = fs_project::InterfaceState::DryContact {
        pressure: fs_qty::QtyAny::new(pressure, fs_project::spec::dims::PRESSURE),
        pressure_half_width: fs_qty::QtyAny::new(width, fs_project::spec::dims::PRESSURE),
        finish: "fixture-machined".into(),
    };
    let cards = CardPackSet::admit(vec![
        RawCardPack { kind: CardPackKind::Interface, source: "pressure-response".into(), bytes, expect: None },
        RawCardPack { kind: CardPackKind::Material, source: "solid".into(),
            bytes: std::fs::read(f.dir.join("solid.fsmcdpk")).unwrap(), expect: None },
    ]).unwrap();
    (project, cards)
}

fn solve_pressure(f: &Fixture, pressure: f64, width: f64, ordinal: usize) -> (J, J) {
    let (project, _) = state_project(f, pressure, width);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let source = f.dir.join(format!("pressure-{ordinal}.fsim"));
    let contact = f.dir.join("pressure.fsintpk");
    std::fs::write(&source, fs_project::print_sexpr(&project).unwrap()).unwrap();
    std::fs::write(&contact, pressure_pack().0).unwrap();
    let ledger = f.dir.join("pressure.db");
    command(&["--json", "import", source.to_str().unwrap(),
        root.join("examples/contact-pair/cold-body.stl").to_str().unwrap(),
        root.join("examples/contact-pair/hot-body.stl").to_str().unwrap(),
        ledger.to_str().unwrap(), "--unit", "m", "--max-hole-edges", "0"]);
    let result = command(&["--json", "solve", source.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", f.dir.join("solid.fsmcdpk").to_str().unwrap(),
        "--materials", root.join("data/reference-project/gray-surface.fsmcdpk").to_str().unwrap(),
        "--interfaces", contact.to_str().unwrap()]);
    let receipt = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|s| s.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
}

#[test]
fn manufactured_pressure_reaches_actual_contact_solves_without_freezing_other_physics() {
    for nonlinear in [false, true] { for coupled in [false, true] {
        let f = Fixture::new(nonlinear, coupled);
        let mut maxima = Vec::new();
        for (i, pressure) in [250.0, 350.0].into_iter().enumerate() {
            let resistance = 0.30 - 0.0005 * pressure;
            let (project, cards) = state_project(&f, pressure, 25.0);
            let bound = fs_project::resolve_bindings(&project, &cards.library(),
                &fs_project::BindingRequirements::thermal_steady_v1());
            assert!(bound.admissible(), "{:?}", bound.violations);
            let interface = bound.bindings.iter().find(|b|
                matches!(&b.target, fs_project::BindingTarget::Interface(name) if name == "cold-hot-joint")).unwrap();
            let property = &interface.properties[0];
            assert!((property.value_lo - resistance).abs() < 1e-14);
            assert_eq!(property.value_lo, property.value_hi, "nominal state is constant in temperature");
            let (_, field) = solve_pressure(&f, pressure, 25.0, i);
            let (_, oracle, _) = f.solve(resistance, None, i + 40);
            let actual = field.get("temperature").unwrap().as_array().unwrap();
            let expected = oracle.get("temperature").unwrap().as_array().unwrap();
            assert_eq!(actual.len(), expected.len());
            for (a, b) in actual.iter().zip(expected) {
                assert!((a.as_f64().unwrap() - b.as_f64().unwrap()).abs() < 1e-7,
                    "state-dependent and independently fixed R'' solves disagree");
            }
            maxima.push(actual.iter().map(|v| v.as_f64().unwrap()).fold(f64::NEG_INFINITY, f64::max));
        }
        assert!(maxima[0] > maxima[1] + 1e-4, "changing only clamping pressure must change this physical answer");
    } }
}

#[test]
fn contact_band_support_is_not_a_nominal_temperature_perturbation() {
    let f = Fixture::new(false, false);
    let (_, narrow) = solve_pressure(&f, 300.0, 5.0, 50);
    let (_, wide) = solve_pressure(&f, 300.0, 75.0, 51);
    assert_eq!(narrow.get("temperature"), wide.get("temperature"),
        "admitted uncertainty width must not replace the nominal contact state");
    let (project, cards) = state_project(&f, 300.0, 110.0);
    let refused = fs_project::resolve_bindings(&project, &cards.library(),
        &fs_project::BindingRequirements::thermal_steady_v1());
    assert!(!refused.admissible(), "nominal pressure is valid but the manufactured band exceeds the card");
    assert!(!refused.bindings.iter().any(|b| matches!(b.target, fs_project::BindingTarget::Interface(_))));
}
