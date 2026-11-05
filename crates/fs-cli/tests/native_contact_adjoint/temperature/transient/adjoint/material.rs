//! Whole-law material controls use the retained live-contact capacity history.
use super::*;

const FINAL_MATERIAL: &str = "temperature-final-material-adjoint";
const PEAK_MATERIAL: &str = "temperature-sampled-peak-material-adjoint";
const DK: &str = "dtemperature_dconductivity_multiplier_k";
const DR: &str = "dtemperature_dresistance_multiplier_k";

fn rows<'a>(value: &'a J, key: &str) -> &'a [J] {
    value.get(key).unwrap().as_array().unwrap()
}

// Change every ordinate, retaining the original temperature/pressure axes and
// interpolation. A perturbation must scale k'(T) and R'(T) along with k and R.
fn scale_value(value: &mut PropertyValue, scale: f64) {
    match value {
        PropertyValue::Scalar { value, .. } => *value *= scale,
        PropertyValue::Curve { knots, .. } => {
            for (_, value) in knots {
                *value *= scale;
            }
        }
        PropertyValue::Surface2D { samples, .. } => {
            for (_, value) in samples {
                *value *= scale;
            }
        }
        _ => panic!("scalar, curve or surface fixture"),
    }
}

fn scaled_material(nonlinear: bool, scale: f64) -> (Vec<u8>, String, String) {
    let (bytes, _, _) = material_pack(nonlinear);
    let original = NormalizedMaterialCardPack::from_bytes(&bytes).unwrap();
    let claim = original.card().claims().claims_ordered().next().unwrap().1;
    let mut value = claim.value.clone();
    scale_value(&mut value, scale);
    let claims = pack(
        original.pack_id(),
        "thermal-conductivity",
        fs_conduction::CONDUCTIVITY_DIMS,
        value,
        claim.interpolation,
    );
    let scaled = NormalizedMaterialCardPack::new(original.card().id().clone(), claims).unwrap();
    let result = (
        scaled.to_bytes(),
        scaled.card().content_hash().to_hex(),
        scaled.card().id().to_string(),
    );
    if scale == 1.0 {
        assert_eq!(bytes, result.0);
    }
    result
}

fn scaled_contact(surface: bool, scale: f64) -> (Vec<u8>, String, String) {
    let original = contact_pack(surface);
    let decoded = NormalizedInterfacePack::from_bytes(&original.0).unwrap();
    let mut value = decoded
        .card()
        .claims()
        .claims_ordered()
        .next()
        .unwrap()
        .1
        .value
        .clone();
    scale_value(&mut value, scale);
    let result = temperature_pack_from_value(value);
    if scale == 1.0 {
        assert_eq!(original, result);
    } else {
        assert_ne!(original.1, result.1);
    }
    result
}

/// Fresh native import/solve with independent cold k, hot k and joint R cards.
/// No derivative output, retained temperature or nominal response is supplied.
fn physical_probe(
    f: &Fixture,
    ordinal: usize,
    nonlinear: bool,
    surface: bool,
    scales: [f64; 3],
) -> J {
    let mut probe = Fixture {
        dir: f.dir.clone(),
        project: f.project.clone(),
    };
    let mut materials = Vec::new();
    for binding in probe.project.materials.as_mut().unwrap() {
        let index = match binding.region.as_str() {
            "cold" => 0,
            "hot" => 1,
            _ => panic!("contact-pair region"),
        };
        let (bytes, card, state) = scaled_material(nonlinear, scales[index]);
        if scales[index] == 1.0 {
            assert_eq!(card, binding.card);
        } else {
            assert_ne!(card, binding.card);
        }
        binding.card = card;
        binding.state = state;
        let path = f
            .dir
            .join(format!("material-probe-{ordinal}-{index}.fsmcdpk"));
        std::fs::write(&path, bytes).unwrap();
        materials.push(path);
    }
    let (source, contact, ledger) =
        prepare_with_pack(&probe, ordinal, None, scaled_contact(surface, scales[2]));
    let mut args = vec![
        "--json".to_string(),
        "solve".into(),
        source.display().to_string(),
        ledger.display().to_string(),
        "--interfaces".into(),
        contact.display().to_string(),
    ];
    for path in materials {
        args.extend(["--materials".into(), path.display().to_string()]);
    }
    let result = fs_cli::run(args);
    assert_eq!(
        result.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        result.stdout,
        result.stderr
    );
    let result = J::parse(&result.stdout).unwrap();
    let run = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = rows(&run, "stages")
        .iter()
        .find(|row| row.str_field("stage") == Some("conduction"))
        .unwrap();
    artifact(&ledger, stage.str_field("receipt").unwrap())
}

fn controls(f: &Fixture, receipt: &J, surface: bool) -> [f64; 3] {
    let adjoint = report(receipt);
    let material = adjoint.get("material_controls").unwrap();
    assert_eq!(n(material, "control_steps"), n(adjoint, "state_index"));
    assert_eq!(material.str_field("units"), Some("K"));
    let conductivity = rows(material, "conductivity_regions");
    assert_eq!(conductivity.len(), 2);
    // Sharing one original card must retain independent declared coordinates.
    assert_eq!(
        conductivity[0].str_field("material_card"),
        conductivity[1].str_field("material_card")
    );
    let mut result = [0.0; 3];
    for (index, name) in ["cold", "hot"].into_iter().enumerate() {
        let row = conductivity
            .iter()
            .find(|row| row.str_field("region") == Some(name))
            .unwrap();
        let binding = f
            .project
            .materials
            .as_ref()
            .unwrap()
            .iter()
            .find(|binding| binding.region == name)
            .unwrap();
        assert_eq!(row.str_field("material_card"), Some(binding.card.as_str()));
        assert_eq!(
            row.str_field("material_state"),
            Some(binding.state.as_str())
        );
        assert!(n(row, "elements") > 0.0);
        assert_eq!(n(row, "reference_value"), 1.0);
        assert_eq!(n(row, "log_reference_value"), 0.0);
        result[index] = n(row, DK);
    }
    let contacts = rows(material, "contact_interfaces");
    assert_eq!(contacts.len(), 1);
    let contact = &contacts[0];
    let (_, card, claim) = contact_pack(surface);
    assert_eq!(contact.str_field("interface"), Some("cold-hot-joint"));
    assert_eq!(contact.str_field("interface_card"), Some(card.as_str()));
    assert_eq!(contact.get("mapped"), Some(&J::Bool(true)));
    assert!(n(contact, "face_pairs") > 0.0);
    assert_eq!(n(contact, "reference_value"), 1.0);
    assert_eq!(n(contact, "log_reference_value"), 0.0);
    result[2] = n(contact, DR);
    for grid in ["coarse", "fine"] {
        for step in rows(receipt.get("transient").unwrap(), grid) {
            for state in rows(step.get("temperature_contact").unwrap(), "face_states") {
                assert_eq!(state.str_field("interface_card"), Some(card.as_str()));
                assert_eq!(state.str_field("selected_claim"), Some(claim.as_str()));
            }
        }
    }
    result
}

fn material_fixture(nonlinear: bool) -> Fixture {
    let mut f = fixture(nonlinear, 0.5);
    f.project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = ThermalBoundaryCondition::Convection {
        coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
    };
    f.project.power.as_mut().unwrap()[0].watts.value = 24.0;
    f.project.power.as_mut().unwrap()[0].duty = 1.0;
    f
}

fn compare_histories(f: &Fixture, nonlinear: bool, surface: bool, peak: bool) {
    let output = if peak { PEAK_MATERIAL } else { FINAL_MATERIAL };
    let (plain, plain_field) = solve(f, 400, None, surface);
    let (receipt, field) = solve(f, 401, Some(output), surface);
    unchanged(&plain, &plain_field, &receipt, &field);
    let adjoint = report(&receipt);
    assert_eq!(adjoint.str_field("output"), Some(output));
    assert_eq!(n(adjoint, "value_k"), qoi(&plain, peak));
    let count = rows(plain.get("transient").unwrap(), "fine").len();
    if peak {
        assert!(n(adjoint, "state_index") > 0.0 && n(adjoint, "state_index") < count as f64);
        assert!(qoi(&plain, true) > qoi(&plain, false) + 1e-3);
    } else {
        assert_eq!(n(adjoint, "state_index"), count as f64);
    }
    let gradients = controls(f, &receipt, surface);
    assert!(
        gradients[..2].iter().all(|value| value.abs() > 1e-6),
        "each regional material must influence the selected history: {gradients:?}"
    );
    assert!(
        gradients[2] > 1e-6,
        "joint resistance retains heat in the hot body"
    );
    // Nonunit fine steps make a missing or duplicate dt in the material/contact
    // contractions observable. Each probe restarts the entire physical history.
    assert!(
        rows(plain.get("transient").unwrap(), "fine")
            .iter()
            .all(|row| n(row, "dt_s") < 0.5)
    );
    let epsilon = 0.002;
    for (coordinate, actual) in gradients.into_iter().enumerate() {
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut scales = [1.0; 3];
            scales[coordinate] += sign * epsilon;
            values.push(qoi(
                &physical_probe(f, 402 + 2 * coordinate + side, nonlinear, surface, scales),
                peak,
            ));
        }
        agrees(actual, (values[1] - values[0]) / (2.0 * epsilon));
    }
    let (legacy, legacy_field) = solve(f, 408, Some(if peak { PEAK } else { FINAL }), surface);
    assert_eq!(legacy_field.get("temperature"), field.get("temperature"));
    let legacy = report(&legacy);
    assert!(legacy.get("material_controls").is_none());
    for key in [
        "value_k",
        "state_index",
        "time_s",
        "active_vertex",
        "active_ties",
        "dtemperature_duniform_initial_k_per_k",
        "dtemperature_dcapacity_multiplier_k",
        "regional_power",
        "intervals",
        "robin",
        "prescribed_temperature",
        "adjoint_krylov_iterations",
        "maximum_dual_relative_residual",
        "temperature_contact",
    ] {
        assert_eq!(
            legacy.get(key),
            adjoint.get(key),
            "material controls changed the shared dual: {key}"
        );
    }
}

#[test]
fn live_contact_final_material_adjoint_matches_regional_and_whole_joint_law_resolves() {
    for nonlinear in [false, true] {
        let mut f = material_fixture(nonlinear);
        time(&mut f).power_schedules.clear();
        compare_histories(&f, nonlinear, false, false);
    }
}

#[test]
fn live_contact_sampled_peak_material_adjoint_matches_the_selected_history_prefix() {
    let mut f = material_fixture(true);
    let policy = time(&mut f);
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(300.0, dims::TEMPERATURE));
    policy.power_schedules[0].steps[0].watts.value = 80.0;
    policy.power_schedules[0].steps[1].watts.value = 0.0;
    compare_histories(&f, true, true, true);
}

#[test]
fn live_contact_initial_peak_has_zero_future_material_and_joint_controls() {
    let mut f = material_fixture(true);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    let policy = time(&mut f);
    policy.power_schedules.clear();
    policy.initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(320.0, dims::TEMPERATURE));
    let (plain, plain_field) = solve(&f, 420, None, false);
    let (receipt, field) = solve(&f, 421, Some(PEAK_MATERIAL), false);
    unchanged(&plain, &plain_field, &receipt, &field);
    let adjoint = report(&receipt);
    assert_eq!(n(adjoint, "state_index"), 0.0);
    assert_eq!(n(adjoint, "value_k"), 320.0);
    assert!(qoi(&plain, false) < 320.0);
    assert_eq!(n(adjoint, "adjoint_krylov_iterations"), 0.0);
    assert_eq!(n(adjoint, "dtemperature_duniform_initial_k_per_k"), 1.0);
    assert_eq!(controls(&f, &receipt, false), [0.0; 3]);
}
