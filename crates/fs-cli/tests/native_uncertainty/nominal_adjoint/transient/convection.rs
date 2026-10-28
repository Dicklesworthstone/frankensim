//! G1/G3/G5: one true enthalpy adjoint calibrates both original Robin inputs.
//! Physical finite differences use independent native projects and full fields.

use super::*;

fn conductivity_card() -> fs_matdb::NormalizedMaterialCardPack {
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
        NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, UncertaintyModel,
    };
    let source = fs_blake3::hash_bytes(b"synthetic native UQ nonlinear conductivity fixture");
    let provenance = || Provenance {
        source: "synthetic k(T) for an independent native UQ regression".into(),
        license: "CC-BY-4.0; redistribution permitted with attribution".into(),
        artifact: Some(source),
    };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic conductivity coupon".into(),
        method: "authored regression curve; not measured data".into(),
        artifact: source,
        caveats: "test oracle only; no material validation".into(),
        provenance: provenance(),
    }).unwrap();
    let conductivity = fs_qty::Dims([1, 1, -3, -1, 0, 0]);
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("thermal-conductivity", conductivity),
        value: PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: dims::TEMPERATURE,
            knots: vec![(233.15, 1.0), (398.15, 4.3)], dims: conductivity,
        },
        validity: fs_evidence::ValidityDomain::unconstrained().with("T", 233.15, 398.15),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![observation],
        provenance: provenance(),
    }).unwrap();
    NormalizedMaterialCardPack::new(
        MaterialStateId { chemistry: "synthetic-uq".into(), phase: "solid".into(),
            process: "declared-curve".into(), revision: 0 },
        NormalizedPack::new("native-uq-enthalpy-convection",
            "frankensim-material-card-pack-compiler-v1", source,
            "CC-BY-4.0: redistribution permitted with attribution",
            claims, Vec::new(), Vec::new()).unwrap(),
    ).unwrap()
}

fn fixture(controlled: bool) -> Fixture {
    let mut fixture = enthalpy_fixture();
    let mut time = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.take().unwrap();
    fixture.project = fs_project::parse_sexpr_migrating(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-enclosure.fsim"
    ))).unwrap().decoded.spec;
    let card = conductivity_card();
    std::fs::write(fixture.sources.join("nonlinear.fsmcdpk"), card.to_bytes()).unwrap();
    let binding = &mut fixture.project.materials.as_mut().unwrap()[0];
    binding.card = card.card().content_hash().to_hex();
    binding.state = card.card().id().to_string();
    let material = &mut time.enthalpy.as_mut().unwrap().materials[0];
    material.region = "body".into();
    material.material_card_identity = fs_blake3::ContentHash::from_hex(&binding.card).unwrap();
    time.energy_tolerance.value = 1e-7;
    let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient = Some(time);
    let enclosure = setup.enclosure_radiation.as_mut().unwrap();
    enclosure.temperature_tolerance.value = 1e-10;
    enclosure.heat_tolerance.value = 1e-9;
    fixture.project.envelope.as_mut().unwrap().ambient_lo.value = 280.0;
    fixture.project.envelope.as_mut().unwrap().ambient_hi.value = 320.0;
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-11;
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
    std::fs::copy(reference.join("gray-surface.fsmcdpk"),
        fixture.sources.join("gray-surface.fsmcdpk")).unwrap();
    write_project(&fixture);
    let text = study_text(false, false, controlled)
        .replace(":wall-time 120s", ":wall-time 240s")
        .replace(":materials (\"aa6061.fsmcdpk\")",
            ":materials (\"nonlinear.fsmcdpk\" \"gray-surface.fsmcdpk\")")
        .replace("(uniform :name \"power\" :target power :entity \"air\" :low 4W :high 6W)",
            "(uniform :name \"htc\" :target convection-coefficient :entity \"cold-patch\" :low 12.0kg·s^-3·K^-1 :high 16.0kg·s^-3·K^-1)")
        .replace("(uniform :name \"ambient\" :target convection-temperature :entity \"air\" :low 294K :high 300K)",
            "(uniform :name \"reference\" :target convection-temperature :entity \"cold-patch\" :low 294K :high 300K)");
    std::fs::write(fixture.sources.join("study.fsim"), text).unwrap();
    fixture
}

fn set_boundary(project: &mut fs_project::ProjectSpec, coefficient_value: f64, reference: f64) {
    let boundary = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .boundaries.iter_mut().find(|row| row.target == "cold-patch").unwrap();
    let fs_project::ThermalBoundaryCondition::Convection { coefficient, reference_temperature }
        = &mut boundary.condition else { panic!("original prescribed convection"); };
    coefficient.value = coefficient_value;
    reference_temperature.value = reference;
}

fn physical_field(fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize) -> JsonValue {
    let source = fixture.dir.join(format!("enthalpy-convection-probe-{ordinal}.fsim"));
    std::fs::write(&source, fs_project::print_sexpr(project).unwrap()).unwrap();
    let ledger = fixture.dir.join("enthalpy-convection-probes.db");
    command(&["--json", "import", source.to_str().unwrap(),
        fixture.sources.join("plate.stl").to_str().unwrap(), ledger.to_str().unwrap(),
        "--unit", "m", "--max-hole-edges", "0"], fs_cli::exit::SUCCESS);
    let solved = command(&["--json", "solve", source.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", fixture.sources.join("nonlinear.fsmcdpk").to_str().unwrap(),
        "--materials", fixture.sources.join("gray-surface.fsmcdpk").to_str().unwrap()],
        fs_cli::exit::SUCCESS);
    let receipt = json_artifact(&ledger, solved.str_field("run_receipt").unwrap());
    let conduction = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let conduction = json_artifact(&ledger, conduction.str_field("receipt").unwrap());
    json_artifact(&ledger, conduction.str_field("solution_artifact").unwrap())
        .get("temperature").unwrap().clone()
}

fn calibration_receipt(fixture: &Fixture, report: &JsonValue) -> JsonValue {
    let probes = report.path(&["mean_control", "calibration", "probes"]).unwrap()
        .as_array().unwrap();
    assert_eq!(probes.len(), 1);
    let qoi = json_artifact(&fixture.ledger, probes[0].str_field("qoi_receipt").unwrap());
    json_artifact(&fixture.ledger,
        qoi.path(&["lineage", "conduction_receipt"]).unwrap().as_str().unwrap())
}

#[test]
fn g1_g3_enthalpy_convection_calibration_matches_nonlinear_enclosure_resolves_and_preserves_samples() {
    let raw = fixture(false);
    let result = raw.study(None, fs_cli::exit::SUCCESS);
    let (_, baseline) = raw.report(&result);
    let controlled = fixture(true);
    let result = controlled.study(None, fs_cli::exit::SUCCESS);
    let (_, report) = controlled.report(&result);
    assert_eq!(rows(&report), rows(&baseline),
        "calibration cannot change any raw physical probability project or trajectory");
    assert_eq!(report.get("statistics"), baseline.get("statistics"));
    let control = report.get("mean_control").unwrap();
    assert_eq!(control.str_field("status"), Some("frozen"));
    assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
    assert_eq!(control.get("estimate").unwrap().f64_field("samples_in_estimate"), Some(4.0));
    let coefficients = gradient(&report);
    assert_eq!(coefficients.len(), 2);
    assert!(coefficients[0].abs() > 1e-6 && coefficients[1] > 0.0);

    let receipt = calibration_receipt(&controlled, &report);
    let history = receipt.get("transient").unwrap();
    let adjoint = history.get("adjoint").unwrap();
    assert_eq!(adjoint.str_field("output"), Some("temperature-final-adjoint"));
    assert_eq!(adjoint.f64_field("state_index"), history.f64_field("fine_steps"));
    let controls = adjoint.get("prescribed_convection").unwrap();
    assert_eq!(controls.f64_field("control_steps"), history.f64_field("fine_steps"));
    let boundary = controls.get("boundaries").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("target") == Some("cold-patch")).unwrap();
    // Means differ from the authored h=10, Tref=290; identities must bind the
    // actual calibration project rather than the original unsampled project.
    assert_eq!(boundary.f64_field("coefficient_w_m2_k"), Some(14.0));
    assert_eq!(boundary.f64_field("reference_temperature_k"), Some(297.0));
    assert_eq!(coefficients[0].to_bits(), boundary
        .f64_field("dtemperature_dcoefficient_k_per_w_m2_k").unwrap().to_bits());
    assert_eq!(coefficients[1].to_bits(), boundary
        .f64_field("dtemperature_dreference_temperature_k_per_k").unwrap().to_bits());
    assert!(history.get("nonlinear").unwrap().f64_field("total_updates").unwrap() > 0.0);
    let surfaces = receipt.get("radiation").unwrap().get("surfaces").unwrap().as_array().unwrap();
    assert!(surfaces[0].f64_field("nonlinear_heat_w").unwrap().abs() > 1e-5);

    let vertex = adjoint.f64_field("active_vertex").unwrap() as usize;
    let mut nominal = controlled.project.clone();
    set_boundary(&mut nominal, 14.0, 297.0);
    let nominal_field = physical_field(&controlled, &nominal, 0);
    let retained_field = json_artifact(&controlled.ledger,
        receipt.str_field("solution_artifact").unwrap());
    assert_eq!(Some(&nominal_field), retained_field.get("temperature"),
        "requesting the calibration transpose leaves the accepted primal untouched");
    for (coordinate, step) in [(0, 0.01), (1, 0.02)] {
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut perturbed = nominal.clone();
            set_boundary(&mut perturbed,
                14.0 + if coordinate == 0 { sign*step } else { 0.0 },
                297.0 + if coordinate == 1 { sign*step } else { 0.0 });
            let field = physical_field(&controlled, &perturbed, 1 + 2*coordinate + side);
            values.push(field.as_array().unwrap()[vertex].as_f64().unwrap());
        }
        let difference = (values[1] - values[0]) / (2.0*step);
        assert!((coefficients[coordinate] - difference).abs() <= 3e-5*difference.abs().max(0.02),
            "coordinate {coordinate}: frozen coefficient {} versus physical difference {difference}",
            coefficients[coordinate]);
    }
}

#[test]
fn g5_enthalpy_convection_calibration_resumes_without_source_files_or_new_probes() {
    let full = fixture(true);
    let result = full.study(None, fs_cli::exit::SUCCESS);
    let (expected, _) = full.report(&result);
    let split = fixture(true);
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    let (_, report) = split.report(&prefix);
    assert!(rows(&report).is_empty(), "the adjoint probe is not a probability observation");
    let frozen = report.path(&["mean_control", "calibration"]).unwrap().clone();
    std::fs::rename(&split.sources, split.dir.join("retained-convection-inputs")).unwrap();
    let finished = split.resume(prefix.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let (actual, restored) = split.report(&finished);
    assert_eq!(restored.path(&["mean_control", "calibration"]), Some(&frozen));
    assert_eq!(actual, expected, "recovered coefficients and raw trajectories replay exactly");
}

