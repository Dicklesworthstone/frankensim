//! Original material-card emissivity is evaluated at the evolving enthalpy state.
//! The imported regular tetrahedron admits an independent scalar energy oracle.

use super::*;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};

const SOURCE_KNOTS: [(f64, f64); 3] = [(200.0, 0.15), (300.0, 0.40), (450.0, 0.85)];
const MATERIAL_MODEL: &str = "pointwise-p1-material-gray-surface-to-black-reservoir";

struct MaterialFixture {
    cards: CardPackSet,
    card: String,
    state: String,
    claim: String,
}

fn material(knots: &[(f64, f64)], specimen: &str) -> MaterialFixture {
    let raw = format!("synthetic original material; specimen={specimen}; epsilon(T)={knots:?}");
    let source = hash_bytes(raw.as_bytes());
    let provenance = || Provenance {
        source: "synthetic combined conductivity and emissivity; not experimental data".into(),
        license: "CC-BY-4.0; redistribution permitted".into(),
        artifact: Some(source),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: specimen.into(),
            method: "authored scalar conductivity and original emissivity ordinates".into(),
            artifact: source,
            caveats: "numerical verification only; no measured material qualification".into(),
            provenance: provenance(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal-conductivity", CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 6.0,
                dims: CONDUCTIVITY_DIMS,
            },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, 450.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let claim = claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::SURFACE_EMISSIVITY_PROPERTY,
                fs_qty::Dims::NONE,
            ),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: dims::TEMPERATURE,
                knots: knots.to_vec(),
                dims: fs_qty::Dims::NONE,
            },
            validity: fs_evidence::ValidityDomain::unconstrained().with(
                "T",
                knots.first().unwrap().0,
                knots.last().unwrap().0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap()
        .0
        .to_hex();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "native-material-emissivity".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        NormalizedPack::new(
            "native-material-emissivity",
            "synthetic-v1",
            source,
            "CC-BY-4.0; redistribution permitted",
            claims,
            Vec::new(),
            Vec::new(),
        )
        .unwrap(),
    )
    .unwrap();
    let cards = CardPackSet::admit(vec![raw_pack(
        CardPackKind::Material,
        "synthetic-original-material.fsmcdpk",
        pack.to_bytes(),
    )])
    .unwrap();
    let card = cards.materials()[0].card().to_hex();
    let state = cards.materials()[0].identity().to_string();
    MaterialFixture {
        cards,
        card,
        state,
        claim,
    }
}

fn fixture(material: &MaterialFixture, reservoir: f64) -> ProjectSpec {
    let mut spec = project(&material.card, reservoir);
    rebind_to(&mut spec, &material.cards);
    let binding = &mut spec.materials.as_mut().unwrap()[0];
    binding.claim = None;
    binding.source = "synthetic original material shared by storage and radiation".into();
    // Retain the full original material and caloric ranges, including source
    // refusal tests. The actual source must support the states being solved.
    close(binding.temp_lo.value, 233.15, 0.0);
    close(binding.temp_hi.value, 398.15, 0.0);
    let radiation = spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .radiation
        .as_mut()
        .unwrap();
    radiation.spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
    radiation.emissivity_response = Some(fs_project::RadiationEmissivityResponse::LocalTemperature);
    radiation.max_iterations = 1;
    radiation.surfaces[0].claim = Some(material.claim.clone());
    let policy = time(&mut spec);
    policy.initial_temperature = fs_project::TransientInitialTemperature::FromRegionalEnthalpy;
    policy.enthalpy.as_mut().unwrap().materials[0].material_card_identity =
        fs_blake3::ContentHash::from_hex(&material.card).unwrap();
    policy.power_schedules = vec![TransientRegionPower {
        region: "air".into(),
        source: "synthetic constant delivered power; independent emissivity oracle".into(),
        steps: vec![TransientPowerStep {
            until: policy.horizon,
            watts: QtyAny::new(5.0, dims::POWER),
        }],
    }];
    spec
}

fn epsilon(temperature: f64, knots: &[(f64, f64)]) -> f64 {
    let pair = knots
        .windows(2)
        .find(|pair| pair[0].0 <= temperature && temperature <= pair[1].0)
        .expect("independent oracle does not extrapolate the original source");
    let fraction = (temperature - pair[0].0) / (pair[1].0 - pair[0].0);
    pair[0].1 + fraction * (pair[1].1 - pair[0].1)
}

fn emitted_w(temperature: f64, reservoir: f64, emissivity: f64) -> f64 {
    emissivity * SIGMA * area_m2() * (temperature.powi(4) - reservoir.powi(4))
}

fn step(old: f64, dt: f64, reservoir: f64, knots: &[(f64, f64)], local: bool) -> f64 {
    let residual = |temperature: f64| {
        let emissivity = epsilon(if local { temperature } else { 300.0 }, knots);
        mass_kg(DENSITY) * CP * (temperature - old) / dt
            + convective_w(temperature)
            + emitted_w(temperature, reservoir, emissivity)
            - 5.0
    };
    let (mut low, mut high) = (233.15, 398.15);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (low + high);
        if residual(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    0.5 * (low + high)
}

fn check_source_states(report: &J, material: &MaterialFixture, curve: &str, temperature: f64) {
    assert_eq!(report.str_field("model"), Some(MATERIAL_MODEL));
    assert_eq!(report.str_field("spatial_model"), Some("pointwise"));
    assert_eq!(
        report.str_field("emissivity_response"),
        Some("local-temperature")
    );
    let surfaces = report.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 1);
    let surface = &surfaces[0];
    assert_eq!(surface.str_field("card"), Some(material.card.as_str()));
    assert_eq!(
        surface.str_field("material_state"),
        Some(material.state.as_str())
    );
    assert_eq!(
        surface.str_field("selected_claim"),
        Some(material.claim.as_str())
    );
    close(
        n(surface, "reference_emissivity"),
        epsilon(300.0, &SOURCE_KNOTS),
        1e-15,
    );
    close(n(surface, "query_temperature_k"), 300.0, 0.0);
    let states = surface.get("material_states").unwrap().as_array().unwrap();
    assert_eq!(
        states.len(),
        4,
        "one source query per original surface vertex"
    );
    let mut vertices = Vec::new();
    for state in states {
        vertices.push(n(state, "vertex") as usize);
        close(n(state, "temperature_k"), temperature, 2e-8);
        close(n(state, "liquid_mass_fraction"), 0.0, 0.0);
        close(
            n(state, "emissivity"),
            epsilon(temperature, &SOURCE_KNOTS),
            1e-10,
        );
        close(
            n(state, "specific_enthalpy_j_kg"),
            CP * (temperature - INITIAL_K),
            3e-7,
        );
        assert_eq!(state.str_field("phase_curve_identity"), Some(curve));
        assert!(
            fs_blake3::ContentHash::from_hex(state.str_field("material_state_identity").unwrap())
                .is_some()
        );
        assert_eq!(
            state.str_field("selected_claim"),
            Some(material.claim.as_str())
        );
        assert!(
            fs_blake3::ContentHash::from_hex(state.str_field("usage_receipt").unwrap()).is_some()
        );
    }
    vertices.sort_unstable();
    assert_eq!(vertices, vec![0, 1, 2, 3]);
}

#[test]
fn g1_local_emissivity_uses_original_material_at_each_accepted_enthalpy_endpoint() {
    let material = material(&SOURCE_KNOTS, "original-specimen");
    for reservoir in [270.0, 400.0] {
        let spec = fixture(&material, reservoir);
        let canonical = print_sexpr(&spec).unwrap();
        assert!(canonical.contains("local-temperature"));
        let solved = run(&spec, &material.cards);
        assert_eq!(solved.temperature.len(), 4);
        let history = solved.transient();
        let source = &history.get("materials").unwrap().as_array().unwrap()[0];
        let curve = source.str_field("phase_curve_identity").unwrap();
        assert_eq!(
            source.str_field("material_card_identity"),
            Some(material.card.as_str())
        );
        for grid in ["coarse", "fine"] {
            let mut old = INITIAL_K;
            for row in solved.rows(grid) {
                let dt = n(row, "dt_s");
                let expected = step(old, dt, reservoir, &SOURCE_KNOTS, true);
                let heat = emitted_w(expected, reservoir, epsilon(expected, &SOURCE_KNOTS));
                let storage = mass_kg(DENSITY) * CP * (expected - old);
                close(n(row, "final_region_max_k"), expected, 2e-8);
                close(n(row, "source_w"), 5.0, 1e-12);
                close(n(row, "stored_energy_change_j"), storage, 3e-6);
                close(
                    n(row, "net_input_w"),
                    5.0 - convective_w(expected) - heat,
                    3e-6,
                );
                let radiation = row.get("radiation").unwrap();
                close(n(radiation, "solid_solves"), 1.0, 0.0);
                close(n(radiation, "nonlinear_radiative_out_w"), heat, 3e-6);
                close(
                    n(radiation, "convective_out_w"),
                    convective_w(expected),
                    3e-6,
                );
                close(
                    storage - dt * (5.0 - convective_w(expected) - heat),
                    0.0,
                    3e-6,
                );
                assert!(
                    n(radiation, "physical_residual_j")
                        <= n(radiation, "physical_residual_threshold_j")
                );
                assert!(
                    n(radiation, "physical_energy_residual_j").abs()
                        <= n(history, "energy_tolerance_j")
                );
                check_source_states(
                    radiation.get("exchange").unwrap(),
                    &material,
                    curve,
                    expected,
                );
                old = expected;
            }
            if grid == "fine" {
                for &temperature in &solved.temperature {
                    close(temperature, old, 2e-8);
                }
                let report = solved.receipt.get("radiation").unwrap();
                check_source_states(report, &material, curve, old);
                close(
                    n(report, "nonlinear_radiative_out_w"),
                    emitted_w(old, reservoir, epsilon(old, &SOURCE_KNOTS)),
                    3e-6,
                );
            }
        }
        assert!(solved.html.contains("local-temperature"));
        assert!(solved.html.contains("Estimated"));
    }
}

#[test]
fn g1_omitted_emissivity_response_keeps_fixed_query_physics_on_the_same_original_card() {
    let material = material(&SOURCE_KNOTS, "original-specimen");
    let mut spec = fixture(&material, 400.0);
    let local = run(&spec, &material.cards);
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .radiation
        .as_mut()
        .unwrap()
        .emissivity_response = None;
    assert!(!print_sexpr(&spec).unwrap().contains(":emissivity-response"));
    let fixed = run(&spec, &material.cards);
    let report = fixed.receipt.get("radiation").unwrap();
    assert_eq!(
        report.str_field("model"),
        Some("pointwise-p1-gray-surface-to-black-reservoir")
    );
    let surface = &report.get("surfaces").unwrap().as_array().unwrap()[0];
    close(n(surface, "emissivity"), epsilon(300.0, &SOURCE_KNOTS), 0.0);
    let mut old = INITIAL_K;
    for row in fixed.rows("fine") {
        old = step(old, n(row, "dt_s"), 400.0, &SOURCE_KNOTS, false);
        close(n(row, "final_region_max_k"), old, 2e-8);
    }
    for &temperature in &fixed.temperature {
        close(temperature, old, 2e-8);
    }
    assert!(
        (local.temperature[0] - fixed.temperature[0]).abs() > 0.05,
        "the actual epsilon(T) trajectory must differ from the frozen 300 K source query"
    );
}

fn refused_before_solution(spec: &ProjectSpec, cards: &CardPackSet) -> SolveRefusal {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, transient_radiation::regular_tetra_stl());
    let decoded = decode(spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let refusal = run_solve(
        &ledger,
        &CancelGate::new_clock_free(),
        &mut benign_clock(),
        &decoded,
        cards,
        &mut Vec::new(),
    )
    .unwrap_err();
    assert_eq!(refusal.stage, Some("conduction"), "{refusal:?}");
    assert_eq!(
        refusal.code, "cli-solve-conduction-transient",
        "{refusal:?}"
    );
    assert_eq!(
        stage_receipt_hashes(&ledger, refusal.run.as_ref().unwrap()).len(),
        4,
        "source refusal must publish no conduction field, QoI or report"
    );
    refusal
}

#[test]
fn g0_local_emissivity_refuses_a_different_original_material_card() {
    let original = material(&SOURCE_KNOTS, "original-specimen");
    let other = material(&SOURCE_KNOTS, "different-specimen-same-emissivity-law");
    assert_ne!(original.card, other.card);
    let cards = CardPackSet::admit(
        original
            .cards
            .iter()
            .chain(other.cards.iter())
            .map(|entry| {
                raw_pack(
                    entry.kind(),
                    "separate-original-material-owner",
                    entry.bytes().to_vec(),
                )
            })
            .collect(),
    )
    .unwrap();
    let mut spec = fixture(&original, 400.0);
    let surface = &mut spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .radiation
        .as_mut()
        .unwrap()
        .surfaces[0];
    surface.card = other.card;
    surface.claim = Some(other.claim);
    let refusal = refused_before_solution(&spec, &cards);
    assert!(
        refusal.what.contains("same original material card"),
        "{refusal:?}"
    );
}

#[test]
fn g0_local_emissivity_refuses_runtime_source_support_without_shrinking_storage_domain() {
    let narrow = material(&[(200.0, 0.15), (320.0, 0.80)], "narrow-emissivity-source");
    let spec = fixture(&narrow, 400.0);
    // The authored 300 K source sample passes, whereas the actual initialized
    // 330 K enthalpy state has no source support. Conductivity retains its full
    // 200..450 K support and neither project nor caloric ranges are narrowed.
    let refusal = refused_before_solution(&spec, &narrow.cards);
    assert!(refusal.what.contains("emissivity"), "{refusal:?}");
    assert!(
        refusal.what.contains("does not cover") || refusal.what.contains("extrapolat"),
        "{refusal:?}"
    );
}

#[path = "material_emissivity/adjoint.rs"]
mod adjoint;
