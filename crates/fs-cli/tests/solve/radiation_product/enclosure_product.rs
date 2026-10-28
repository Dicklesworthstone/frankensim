//! G1/G3/G5 native enclosure exchange: real import, material cards, fields and report.
use super::*;
use fs_project::{ConductionEnclosure, EnclosureRadiatingSurface, EnclosureViewFactorEvidence};

fn enclosure_project(card: &str) -> ProjectSpec {
    let mut spec = conduction_fixture_project(7, &tetra_stl());
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    for (name, normal) in [("wall-x", [1.0, 0.0, 0.0]), ("wall-z", [0.0, 0.0, 1.0])] {
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".into(),
            name: name.into(),
            display: "Manufactured radiating patch".into(),
            expect_id: None,
        });
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".into(),
            target: name.into(),
            length_unit: "m".into(),
            selector: MeshSelector::HalfSpace {
                normal,
                offset: 0.0,
                side: HalfSpaceSide::AtMost,
                tolerance: 0.0,
            },
            allow_overlap: true,
        });
    }
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries = [("wall-x", 290.0), ("wall-z", 310.0)]
        .into_iter()
        .map(|(name, reference)| ThermalBoundary {
            target: name.into(),
            condition: ThermalBoundaryCondition::Convection {
                coefficient: QtyAny::new(10.0, fs_project::spec::dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(reference, fs_project::spec::dims::TEMPERATURE),
            },
        })
        .collect();
    setup.adiabatic_remainder = true;
    setup.enclosure_radiation = Some(ConductionEnclosure {
        surfaces: ["wall-x", "wall-z"]
            .into_iter()
            .enumerate()
            .map(|(i, target)| EnclosureRadiatingSurface {
                name: format!("radiator-{i}"),
                target: target.into(),
                card: card.into(),
                claim: None,
                query_temperature: QtyAny::new(300.0, fs_project::spec::dims::TEMPERATURE),
            })
            .collect(),
        view_factors: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        // This is an exact algebraic exchange oracle. The test does not claim
        // that two tetrahedral faces form an actual opaque geometric enclosure.
        evidence: EnclosureViewFactorEvidence::Analytic {
            geometry: "manufactured-equal-area-two-patch-closure; not finite-geometry view factors"
                .into(),
        },
        closure_tolerance: 1e-12,
        reciprocity_relative_tolerance: 1e-12,
        max_iterations: 128,
        temperature_tolerance: QtyAny::new(1e-9, fs_project::spec::dims::TEMPERATURE),
        heat_tolerance: QtyAny::new(1e-6, fs_project::spec::dims::POWER),
        relaxation: 0.5,
    });
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    spec.envelope.as_mut().unwrap().ambient_lo.value = 300.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 300.0;
    spec
}

fn run_enclosure(spec: &ProjectSpec, cards: &CardPackSet) -> (J, J, String) {
    let source = print_sexpr(spec).unwrap();
    assert!(source.contains(":enclosure-radiation (enclosure-radiation"));
    let decoded = fs_project::parse_sexpr(&source).unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, tetra_stl());
    let (run, qoi, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    let hashes = stage_receipt_hashes(&ledger, &run);
    let report = String::from_utf8(artifact_bytes(&ledger, &hashes[6])).unwrap();
    let html = String::from_utf8(artifact_bytes(
        &ledger,
        &receipt_str_field(&report, "report_html"),
    ))
    .unwrap();
    let qoi = J::parse(&qoi).unwrap();
    assert!(qoi.path(&["model_form_sensitivity", "radiation"]).is_some());
    (J::parse(&receipt).unwrap(), J::parse(&field).unwrap(), html)
}

fn check_enclosure(receipt: &J, epsilon: f64) -> f64 {
    let radiation = receipt.get("radiation").unwrap();
    assert_eq!(
        radiation.str_field("model"),
        Some("closed-gray-diffuse-enclosure")
    );
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    assert_eq!(surfaces.len(), 2);
    assert_eq!(surfaces[0].str_field("name"), Some("radiator-0"));
    assert_eq!(surfaces[0].str_field("target"), Some("wall-x"));
    near(n(&surfaces[0], "area_m2"), 0.5, 1e-12);
    near(n(&surfaces[1], "area_m2"), 0.5, 1e-12);
    let a = n(&surfaces[0], "mean_temperature_k");
    let b = n(&surfaces[1], "mean_temperature_k");
    // Independent two-gray-surface formula, with stable fourth-power difference.
    let expected = SIGMA * 0.5 * (a - b) * (a + b) * (a * a + b * b) / (2.0 / epsilon - 1.0);
    assert!(
        expected.abs() > 1e-5,
        "the fixture must exercise mutual exchange"
    );
    near(n(&surfaces[0], "nonlinear_heat_w"), expected, 1e-6);
    near(n(&surfaces[1], "nonlinear_heat_w"), -expected, 1e-6);
    near(n(radiation, "radiative_out_w"), 0.0, 1e-6);
    near(n(radiation, "convective_out_w"), 5.0, 1e-5);
    near(n(radiation, "decomposition_residual_w"), 0.0, 1e-6);
    assert!(n(radiation, "max_heat_mismatch_w") <= 1e-6);
    assert!(n(radiation, "nonlinear_residual_w") <= n(radiation, "nonlinear_residual_threshold_w"));
    assert!(n(radiation, "solid_solves") > n(radiation, "surface_iterations"));
    expected.abs()
}

#[test]
fn g1_g3_native_enclosure_preserves_reflected_exchange_and_original_convection() {
    let mut transfer = Vec::new();
    for epsilon in [0.2, 0.85] {
        let (cards, card) = emissivity_cards(&fixture_cards(), epsilon);
        let spec = enclosure_project(&card);
        let (receipt, _, html) = run_enclosure(&spec, &cards);
        transfer.push(check_enclosure(&receipt, epsilon));
        assert!(html.contains("closed-gray-diffuse-enclosure"), "{html}");
        // A frozen enclosure load must not trigger a linear goal correction.
        let control = receipt.get("solver_control").unwrap();
        assert_ne!(control.str_field("mode"), Some("physical-goal-correction"));
    }
    assert!(
        transfer[1] > transfer[0],
        "the physical emissivity card must change exchange"
    );
}

#[test]
fn g3_g5_native_enclosure_rebinds_refined_traces_and_replays_exact_fields() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.7);
    let mut spec = enclosure_project(&card);
    spec.solver.as_mut().unwrap().fidelity = "ladder".into();
    let (a, field_a, _) = run_enclosure(&spec, &cards);
    let (b, field_b, _) = run_enclosure(&spec, &cards);
    assert_eq!(a, b);
    assert_eq!(field_a, field_b);
    check_enclosure(&a, 0.7);
    let rungs = a.path(&["ladder", "rungs"]).unwrap().as_array().unwrap();
    assert!(
        rungs.len() >= 2,
        "the test must exercise actual spatial refinement"
    );
    assert!(n(rungs.last().unwrap(), "tets") > n(&rungs[0], "tets"));
}

#[test]
fn g0_native_enclosure_refuses_unsupported_bounds_and_bad_mesh_reciprocity() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    for kind in ["mean-bound", "adaptive", "unequal-area"] {
        let mut spec = enclosure_project(&card);
        match kind {
            "mean-bound" => spec.outputs.as_mut().unwrap().push(OutputRequest {
                name: "temperature-volume-mean-bound".into(),
                kind: "report".into(),
                region: None,
            }),
            "adaptive" => spec.solver.as_mut().unwrap().fidelity = "adaptive".into(),
            _ => {
                // The sloping tetrahedral face has area sqrt(3)/2, so an
                // equal-area F12=F21=1 declaration must fail at mesh binding.
                spec.assignments
                    .as_mut()
                    .unwrap()
                    .iter_mut()
                    .find(|row| row.target == "wall-z")
                    .unwrap()
                    .selector = MeshSelector::HalfSpace {
                    normal: [1.0, 1.0, 1.0],
                    offset: 1.0,
                    side: HalfSpaceSide::AtLeast,
                    tolerance: 0.0,
                };
            }
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, tetra_stl());
        let error = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decode(&spec),
            &cards,
            &mut Vec::new(),
        )
        .unwrap_err();
        assert!(
            error.what.contains("enclosure") || error.what.contains("reciprocity"),
            "{kind}: {error:?}"
        );
    }
}
