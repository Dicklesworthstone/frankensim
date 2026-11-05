//! Full temperature x pressure support reaches imported nominal contact solves.
use super::*;
use fs_matdb::{ClaimSelection, QueryPoint, SelectionPolicy};
use fs_project::{InterfaceState, spec::dims};

fn resistance(temperature: f64, pressure: f64) -> f64 {
    0.001 * temperature - 1e-7 * (pressure - 1.0e6)
}

fn surface_pack(hole: bool) -> (Vec<u8>, String, String) {
    // The optional hole lies inside the original ±100 kPa band, away from
    // its nominal pressure and all four declared temperature/pressure corners.
    let temperatures = [200.0, 280.0, 320.0, 450.0];
    let pressures = [0.5e6, 1.03e6, 1.07e6, 1.5e6];
    let samples = pressures.into_iter().flat_map(|pressure|
        temperatures.into_iter().map(move |temperature|
            ([temperature, pressure], resistance(temperature, pressure)))).collect();
    let mut triangles = Vec::new();
    for pressure in 0..3_u32 {
        for temperature in 0..3_u32 {
            if hole && pressure == 1 && temperature == 1 { continue; }
            let a = 4 * pressure + temperature;
            triangles.extend([[a, a + 1, a + 5], [a, a + 5, a + 4]]);
        }
    }
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples, triangles, dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

fn pressure_band(f: &Fixture) -> (f64, f64) {
    let InterfaceState::DryContact { pressure, pressure_half_width, .. } =
        &f.project.interface_cards.as_ref().unwrap()[0].state else { panic!("dry contact") };
    (pressure.value, pressure_half_width.value)
}

fn with_pressure_band(f: &Fixture, nominal: f64, width: f64) -> Fixture {
    let mut changed = Fixture { dir: f.dir.clone(), project: f.project.clone() };
    let InterfaceState::DryContact { pressure, pressure_half_width, .. } =
        &mut changed.project.interface_cards.as_mut().unwrap()[0].state else { panic!("dry contact") };
    pressure.value = nominal;
    pressure_half_width.value = width;
    changed
}

fn full_temperature_range(f: &Fixture) -> (f64, f64) {
    for binding in f.project.materials.as_ref().unwrap() {
        assert_eq!((binding.temp_lo.value, binding.temp_hi.value), (233.15, 398.15));
    }
    (233.15, 398.15)
}

fn solve(f: &Fixture, ordinal: usize, source_pack: (Vec<u8>, String, String)) -> (J, J) {
    full_temperature_range(f);
    let (card, claim) = (source_pack.1.clone(), source_pack.2.clone());
    // Deliberately use the generic writer: the authored manufacturing band
    // must survive source admission and cannot become a nominal line.
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, source_pack);
    let output = run(f, &source, &contact, &ledger);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", output.stdout, output.stderr);
    let result = J::parse(&output.stdout).unwrap();
    let completed = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    let physical = receipt.get("temperature_contact").unwrap();
    assert!(physical.f64_field("residual_norm_w").unwrap()
        <= physical.f64_field("residual_tolerance_w").unwrap());
    assert!(physical.f64_field("energy_closure_w").unwrap().abs()
        <= physical.f64_field("energy_tolerance_w").unwrap());
    let faces = physical.get("face_states").unwrap().as_array().unwrap();
    assert!(!faces.is_empty());
    for face in faces {
        assert_eq!(face.str_field("selected_claim"), Some(claim.as_str()));
        assert_eq!(face.str_field("usage_receipt").unwrap().len(), 64);
        let query = face.f64_field("query_temperature_k").unwrap();
        let actual = face.f64_field("resistance_m2_k_w").unwrap();
        assert!((actual - resistance(query, pressure_band(f).0)).abs() < 1e-13);
    }
    let fluxes = receipt.get("interfaces").unwrap().get("fluxes").unwrap()
        .as_array().unwrap();
    assert!(!fluxes.is_empty());
    for flux in fluxes {
        assert_eq!(flux.str_field("card"), Some(card.as_str()));
    }
    let heat: f64 = fluxes.iter()
        .map(|row| row.f64_field("heat_rate_a_to_b_w").unwrap()).sum();
    assert!((heat.abs() - 12.0 * 0.37).abs() < 1e-7,
        "the hot body's original power still leaves through this contact");
    (receipt, field)
}

fn temperatures(field: &J) -> Vec<f64> {
    field.get("temperature").unwrap().as_array().unwrap().iter()
        .map(|value| value.as_f64().unwrap()).collect()
}

#[test]
fn full_surface_manufacturing_band_preserves_nominal_physics_and_pressure_response() {
    let f = Fixture::new(true, false);
    assert_eq!(pressure_band(&f), (1.0e6, 0.1e6));
    let (surface, actual) = solve(&f, 610, surface_pack(false));
    let (_, curve) = solve(&f, 611, temperature_pack_from_knots(
        vec![(200.0, 0.2), (450.0, 0.45)]));
    let actual_t = temperatures(&actual);
    let curve_t = temperatures(&curve);
    assert_eq!(actual_t.len(), curve_t.len());
    for (a, b) in actual_t.iter().zip(curve_t) {
        assert!((a - b).abs() < 2e-7, "surface at nominal pressure differs from independent R(T)");
    }
    assert_eq!(surface.get("temperature_contact").unwrap().str_field("algorithm"),
        Some("bounded-contact-picard"));

    let wider = with_pressure_band(&f, 1.0e6, 0.3e6);
    let (wide_receipt, wide_field) = solve(&wider, 612, surface_pack(false));
    assert_eq!(actual.get("temperature"), wide_field.get("temperature"),
        "admitted manufacturing uncertainty must not replace the nominal state");
    assert_eq!(surface.get("energy"), wide_receipt.get("energy"));

    let tighter = with_pressure_band(&f, 1.2e6, 0.1e6);
    let (_, compressed) = solve(&tighter, 613, surface_pack(false));
    let maximum = |field: &J| temperatures(field).into_iter().fold(f64::NEG_INFINITY, f64::max);
    assert!(maximum(&actual) > maximum(&compressed) + 1e-3,
        "higher nominal pressure lowers R(T,p) and the actual solved hot temperature");
}

fn assert_refused(f: &Fixture, ordinal: usize, source_pack: (Vec<u8>, String, String)) {
    full_temperature_range(f);
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, source_pack);
    let output = run(f, &source, &contact, &ledger);
    assert_eq!(output.exit_code, fs_cli::exit::REFUSED, "{}\n{}", output.stdout, output.stderr);
    assert!(output.stderr.contains("project-binding-"), "{}", output.stderr);
    let result = J::parse(&output.stdout).unwrap();
    assert_eq!(result.str_field("status"), Some("refused"));
    assert_eq!(result.str_field("stage"), Some("material-resolve"));
    assert!(result.get("run_receipt").is_none());
    assert!(result.get("solution_artifact").is_none(),
        "unsupported source coverage cannot publish a physical solution");
}

#[test]
fn surface_band_refuses_extrapolation_and_a_hole_hidden_from_nominal_and_corner_queries() {
    let f = Fixture::new(false, false);
    assert_eq!(pressure_band(&f), (1.0e6, 0.1e6));
    let outside = with_pressure_band(&f, 1.0e6, 0.6e6);
    assert_refused(&outside, 620, surface_pack(false));

    let holed = surface_pack(true);
    let original = NormalizedInterfacePack::from_bytes(&holed.0).unwrap();
    let claims = original.card().claims();
    let (selected, claim) = claims.claims_ordered()
        .find(|(id, _)| id.0.to_hex() == holed.2).unwrap();
    let point = |temperature, pressure| QueryPoint::new().with("T", temperature).unwrap()
        .with("normal_pressure", pressure).unwrap();
    let (low, high) = full_temperature_range(&f);
    claims.query_line_typed(&claim.key, "T", &point(low, 1.0e6), &point(high, 1.0e6),
        ClaimSelection::Pinned(selected)).expect("the entire nominal T line is supported");
    for temperature in [low, high] {
        for pressure in [0.9e6, 1.1e6] {
            let corner = claims.query(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                &point(temperature, pressure), SelectionPolicy::SingleClaimOnly).unwrap();
            assert_eq!(corner.receipt.selected, selected);
            claims.verify_receipt(&corner.receipt).unwrap();
        }
    }
    // The manufactured band includes the unmeasured interior cell despite its
    // valid nominal line and four corners. An explicit claim pin cannot fill it.
    assert_refused(&f, 621, holed);
}
