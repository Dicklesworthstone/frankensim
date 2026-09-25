//! Real .fsim import -> solve -> QoI -> retained report radiation coverage.

use super::*;
use fs_project::{ConductionRadiation, RadiatingSurface};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const SIGMA: f64 = 5.670_374_419e-8;

fn emissivity_cards(base: &CardPackSet, value: f64) -> (CardPackSet, String) {
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
        ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
        UncertaintyModel,
    };
    let raw_fixture = format!("synthetic gray-surface fixture; T_K=300; emissivity={value}");
    let raw_identity = hash_bytes(raw_fixture.as_bytes());
    let provenance = || Provenance {
        source: "G1 manufactured gray surface; synthetic fixture data".to_string(),
        license: "CC-BY-4.0; redistribution permitted with attribution".to_string(),
        artifact: Some(raw_identity),
    };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured gray surface".to_string(),
        method: "synthetic constant-emissivity fixture declaration".to_string(),
        artifact: raw_identity,
        caveats: "illustrative fixture value, not an experimental observation or validated material property".to_string(),
        provenance: provenance(),
    }).expect("licensed synthetic observation fixture inserts");
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                fs_conduction::SURFACE_EMISSIVITY_PROPERTY,
                fs_qty::Dims::NONE,
            ),
            value: PropertyValue::Scalar {
                value,
                dims: fs_qty::Dims::NONE,
            },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, 450.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "radiation-fixture".to_string(),
            phase: "solid".to_string(),
            process: "gray-finish".to_string(),
            revision: 0,
        },
        NormalizedPack::new(
            "radiation-fixture",
            "manufactured-fixture-v1",
            hash_bytes(b"manufactured gray surface"),
            "CC-BY-4.0; redistribution permitted",
            claims,
            Vec::new(),
            Vec::new(),
        )
        .unwrap(),
    )
    .unwrap();
    let identity = pack.card().content_hash().to_hex();
    let mut raw: Vec<RawCardPack> = base
        .iter()
        .map(|entry| {
            raw_pack(
                entry.kind(),
                "fixture-retained-card",
                entry.bytes().to_vec(),
            )
        })
        .collect();
    raw.push(raw_pack(
        CardPackKind::Material,
        "fixture-emissivity",
        pack.to_bytes(),
    ));
    (CardPackSet::admit(raw).unwrap(), identity)
}

fn declare(spec: &mut ProjectSpec, card: &str, target: &str, reservoir: f64) {
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .radiation = Some(ConductionRadiation {
        surfaces: vec![RadiatingSurface {
            name: "gray-exterior".to_string(),
            target: target.to_string(),
            card: card.to_string(),
            claim: None,
            query_temperature: QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            reservoir_temperature: QtyAny::new(reservoir, fs_project::spec::dims::TEMPERATURE),
        }],
        max_iterations: 128,
        temperature_tolerance: QtyAny::new(1e-8, fs_project::spec::dims::TEMPERATURE),
        heat_tolerance: QtyAny::new(1e-7, fs_project::spec::dims::POWER),
        relaxation: 0.5,
    });
    // The solid residual is scaled by the absolute-temperature load vector,
    // not the 5 W net source. A 1e-9 relative gate can admit O(1e-5 W)
    // source imbalance on these Robin/contact fixtures; 1e-12 puts that
    // algebraic error below the independent 2e-6 W conservation oracle.
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    let envelope = spec.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 293.15;
    envelope.ambient_hi.value = 293.15;
}

/// Extend the existing tracked-reference generator pattern without rewriting
/// the original geometry, project or conductivity pack.
#[test]
#[ignore = "generator: writes only the radiating reference variant and emissivity pack"]
fn dump_radiation_reference_project_fixture() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    let mut spec = conduction_fixture_project(7, &tetra_stl());
    declare(&mut spec, &card, "air", 293.15);
    spec.metadata.as_mut().unwrap().name = "cooling-radiation-reference".to_string();
    spec.metadata.as_mut().unwrap().context_of_use =
        "manufactured cooling example with a declared gray surface; no experimental validation"
            .to_string();
    let pack = cards
        .materials()
        .iter()
        .find(|pack| pack.card().to_hex() == card)
        .unwrap();
    let dir = reference_project_dir();
    std::fs::write(
        dir.join("cooling-radiation.fsim"),
        print_sexpr(&spec).unwrap(),
    )
    .unwrap();
    std::fs::write(dir.join("gray-surface.fsmcdpk"), pack.bytes()).unwrap();
    println!(
        "wrote cooling-radiation.fsim and gray-surface.fsmcdpk in {}",
        dir.display()
    );
}

fn n(row: &J, field: &str) -> f64 {
    row.f64_field(field).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected} (tolerance {tolerance})"
    );
}

fn slab_stl() -> Vec<u8> {
    // Dyadic lengths survive the importer's f32 STL coordinates exactly.
    // The independent oracle therefore uses the actual admitted geometry:
    // L=1/16 m, W=H=1/8 m, exterior area=1/16 m2.
    let points = [
        [0., 0., 0.],
        [0.0625, 0., 0.],
        [0.0625, 0.125, 0.],
        [0., 0.125, 0.],
        [0., 0., 0.125],
        [0.0625, 0., 0.125],
        [0.0625, 0.125, 0.125],
        [0., 0.125, 0.125],
    ];
    let faces = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [3, 7, 6],
        [3, 6, 2],
        [0, 4, 7],
        [0, 7, 3],
        [1, 2, 6],
        [1, 6, 5],
    ];
    let mut stl = String::from("solid slab\n");
    for [a, b, c] in faces {
        stl.push_str(&facet(points[a], points[b], points[c]));
    }
    stl.push_str("endsolid slab\n");
    stl.into_bytes()
}

/// Global balance fixes the area-mean wall temperature independently of the
/// solid's discretization: Q = h A (Tmean-Tair) + eps sigma A(Tmean^4-Trad^4).
/// Bisection here neither assembles FEM nor uses any production radiation API.
fn slab_oracle(reservoir: f64) -> f64 {
    let (mut lo, mut hi) = (200.0_f64, 450.0_f64);
    for _ in 0..100 {
        let t = 0.5 * (lo + hi);
        let out =
            10.0 * 0.0625 * (t - 293.15) + 0.85 * SIGMA * 0.0625 * (t.powi(4) - reservoir.powi(4));
        if out > 5.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_fsim_radiation_slab_matches_balance_and_reaches_retained_report() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    let mut maxima = Vec::new();
    for reservoir in [270.0, 293.15, 350.0] {
        let bytes = slab_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        spec.cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .regions[0]
            .seed = [0.03125, 0.0625, 0.0625]
            .map(|value| QtyAny::new(value, fs_project::spec::dims::LENGTH));
        declare(&mut spec, &card, "air", reservoir);
        if reservoir == 293.15 {
            let solver = spec.solver.as_mut().unwrap();
            solver.fidelity = "ladder".to_string();
            // The 1536-tet rung reached a recomputed linear residual floor
            // of 3.91e-13. The driver's 100x inner tightening makes this
            // explicit 1e-10 stop request a resolvable 1e-12 linear gate.
            // The independent temperature and watt checks stay unchanged.
            solver.tolerance_rel = 1e-10;
        }
        let source = print_sexpr(&spec).unwrap();
        assert!(
            source.contains(":radiation (radiation"),
            "the test must traverse the actual .fsim declaration"
        );
        let decoded = fs_project::parse_sexpr(&source).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (run, qoi, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
        let retained = J::parse(&conduction).unwrap();
        let radiation = retained.get("radiation").unwrap();
        let surface = &radiation.get("surfaces").unwrap().as_array().unwrap()[0];
        near(n(surface, "area_m2"), 0.0625, 1e-12);
        near(
            n(surface, "mean_temperature_k"),
            slab_oracle(reservoir),
            2e-5,
        );
        near(
            n(radiation, "radiative_out_w") + n(radiation, "convective_out_w"),
            5.0,
            2e-6,
        );
        near(
            n(surface, "applied_heat_w"),
            n(surface, "nonlinear_heat_w"),
            1e-7,
        );
        assert!(n(radiation, "max_temperature_change_k") <= 1e-8);
        assert_eq!(
            radiation
                .path(&["controls", "relative_heat_tolerance"])
                .unwrap()
                .as_f64(),
            Some(0.0)
        );
        if reservoir == 350.0 {
            assert!(n(radiation, "radiative_out_w") < 0.0);
        }
        if reservoir == 293.15 {
            assert_eq!(
                retained
                    .path(&["ladder", "rungs"])
                    .unwrap()
                    .as_array()
                    .unwrap()
                    .len(),
                3
            );
        }
        let qoi = J::parse(&qoi).unwrap();
        let sensitivity = qoi.path(&["model_form_sensitivity", "radiation"]).unwrap();
        near(
            n(sensitivity, "delta_on_minus_off_k"),
            n(sensitivity, "radiation_on_k") - n(sensitivity, "radiation_off_k"),
            1e-12,
        );
        maxima.push(n(sensitivity, "radiation_on_k"));
        let receipts = stage_receipt_hashes(&ledger, &run);
        let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
        let html = String::from_utf8(artifact_bytes(
            &ledger,
            &receipt_str_field(&report, "report_html"),
        ))
        .unwrap();
        assert!(
            html.contains("Radiation model-form sensitivity (Estimated)"),
            "{html}"
        );
        assert!(html.contains("radiative_out_w"));
        assert!(html.contains("not a bound on model-form error"));
        assert!(
            qoi.get("budget").unwrap().as_array().unwrap()[0]
                .get("terms")
                .unwrap()
                .as_array()
                .unwrap()
                .iter()
                .any(|term| term.str_field("kind") == Some("model-form")
                    && term.str_field("state") == Some("no-data")),
            "a sensitivity must not fill unknown model-form error"
        );
    }
    assert!(
        maxima.windows(2).all(|pair| pair[0] < pair[1]),
        "hotter surroundings must change the actual solution"
    );
}

#[test]
fn g1_fsim_radiation_airflow_receives_only_convective_watts() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    let bytes = tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    declare(&mut spec, &card, "air", 293.15);
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes);
    let (_, _, conduction, _) = run_conjugate_to_completion(&ledger, &decode(&spec), &cards);
    let retained = J::parse(&conduction).unwrap();
    let radiation = retained.get("radiation").unwrap();
    let air = retained.get("conjugate").unwrap();
    let qrad = n(radiation, "radiative_out_w");
    assert!(
        qrad > 0.01,
        "the fixture must expose leaking radiation into air"
    );
    near(
        n(air, "air_total_w"),
        n(radiation, "convective_out_w"),
        2e-6,
    );
    near(n(air, "air_total_w") + qrad, 5.0, 2e-6);
    near(
        n(air, "mass_flow_kg_s") * 1007.0 * (n(air, "outlet_k") - n(air, "inlet_k")),
        n(air, "air_total_w"),
        2e-6,
    );
    assert!(n(air, "decomposition_residual_w").abs() < 2e-6);
}

#[test]
fn g1_fsim_radiation_preserves_heterogeneous_contact_heat_path() {
    let base = contact_cards();
    let mut raw: Vec<RawCardPack> = base
        .iter()
        .map(|entry| raw_pack(entry.kind(), "contact-fixture-card", entry.bytes().to_vec()))
        .collect();
    raw.push(raw_pack(
        CardPackKind::Material,
        "nonlinear-hot-body",
        material_pack_bytes_with_property(
            "hot-nonlinear",
            "hot-nonlinear-radiation",
            200.0,
            450.0,
            fs_matdb::PropertyValue::Curve {
                abscissa: "T".to_string(),
                abscissa_dims: fs_project::spec::dims::TEMPERATURE,
                knots: vec![(200.0, 10.0), (450.0, 20.0)],
                dims: CONDUCTIVITY_DIMS,
            },
            fs_matdb::InterpolationPolicy::LinearInside,
        ),
    ));
    let base = CardPackSet::admit(raw).unwrap();
    let hot = base
        .materials()
        .iter()
        .find(|entry| entry.identity().starts_with("hot-nonlinear/"))
        .unwrap();
    let (cards, card) = emissivity_cards(&base, 0.7);
    let mut spec = multi_region_contact_project();
    let binding = spec
        .materials
        .as_mut()
        .unwrap()
        .iter_mut()
        .find(|binding| binding.region == "hot")
        .unwrap();
    binding.card = hot.card().to_hex();
    binding.state = hot.identity().to_string();
    spec.requirements.as_mut().unwrap()[0].region = "hot".to_string();
    declare(&mut spec, &card, "hot", 293.15);
    let ledger = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&ledger, &spec);
    let (_, _, conduction, _) = run_conjugate_to_completion(&ledger, &decode(&spec), &cards);
    let retained = J::parse(&conduction).unwrap();
    let fluxes = retained
        .path(&["interfaces", "fluxes"])
        .unwrap()
        .as_array()
        .unwrap();
    assert!(!fluxes.is_empty());
    assert!(n(&fluxes[0], "heat_rate_a_to_b_w").abs() > 0.01);
    let radiation = retained.get("radiation").unwrap();
    let energy = retained.get("energy").unwrap();
    near(
        n(radiation, "radiative_out_w") + n(radiation, "convective_out_w")
            - n(energy, "dirichlet_in_w"),
        5.0,
        2e-6,
    );
    assert!(n(radiation, "radiative_out_w") > 0.01);
}

#[test]
fn g0_fsim_radiation_refuses_adaptive_and_missing_emissivity_claim() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    let conductivity_only_card = fixture_card_identity().0;
    for adaptive in [false, true] {
        let bytes = tetra_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        declare(
            &mut spec,
            if adaptive {
                &card
            } else {
                &conductivity_only_card
            },
            "air",
            293.15,
        );
        if adaptive {
            spec.solver.as_mut().unwrap().fidelity = "adaptive".to_string();
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let error = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decode(&spec),
            &cards,
            &mut Vec::new(),
        )
        .unwrap_err();
        assert_eq!(
            error.code,
            if adaptive {
                "cli-solve-conduction-radiation-adaptive"
            } else {
                "cli-solve-conduction-radiation-card"
            }
        );
        assert_eq!(error.stage, Some("conduction"));
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4
        );
    }
}
