//! Actual .fsim product regressions, not isolated lab or frozen-Robin solves.
use super::*;

fn adaptive_controls(spec: &mut ProjectSpec) {
    let solver = spec.solver.as_mut().unwrap();
    solver.fidelity = "adaptive".to_string();
    // The driver tightens this to 1e-10 for the complete primal/dual residual.
    // Keep a resolvable gate for independently re-marched air references.
    solver.tolerance_rel = 1e-8;
    spec.budgets.as_mut().unwrap().accuracy_rel = 0.1;
    let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .radiation.as_mut().unwrap();
    radiation.temperature_tolerance.value = 1e-9;
    radiation.heat_tolerance.value = 1e-8;
}

fn check_comparison(conduction: &str, air_feedback: bool) {
    let retained = J::parse(conduction).unwrap();
    let history = retained.path(&["adaptive", "history"]).unwrap().as_array().unwrap();
    assert!(!history.is_empty(), "adaptive radiation must execute an enriched comparison");
    assert!(conduction.contains(&format!("\"air_reference_feedback\":{air_feedback}")), "{conduction}");
    assert!(conduction.contains("\"continuum_error_bound\":false"), "{conduction}");
    for row in history {
        let balance = row.get("algebraic_balance").unwrap();
        assert_eq!(balance.str_field("comparison_model"), Some("area-mean-radiation-full-transpose"));
        assert_eq!(balance.str_field("status"), Some("unsupported-model"),
            "nonlinear radiation has no certified linear algebraic balance");
        let linear = n(row, "signed_linear_change_k");
        let nonlinear = n(row, "linearization_remainder_k");
        let relocation = n(row, "maximum_remainder_k");
        near(linear + nonlinear + relocation, n(row, "measured_change_k"), 1e-9);
        near(linear.abs() + nonlinear.abs() + relocation.abs(), n(row, "estimated_change_k"), 1e-9);
        assert_eq!(n(row, "probe_tets"), 8.0*n(row, "tets"));
        assert!(n(row, "primal_residual") < 1e-10);
        assert!(n(row, "dual_residual") < 1e-10);
    }
}

#[test]
fn g1_fsim_radiation_adaptive_slab_keeps_balance_and_hot_reservoir_sign() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    for reservoir in [293.15, 350.0] {
        let bytes = slab_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().regions[0].seed =
            [0.03125, 0.0625, 0.0625].map(|v| QtyAny::new(v, fs_project::spec::dims::LENGTH));
        declare(&mut spec, &card, "air", reservoir);
        adaptive_controls(&mut spec);
        let source = print_sexpr(&spec).unwrap();
        let decoded = fs_project::parse_sexpr(&source).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (run, qoi, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
        check_comparison(&conduction, false);
        let retained = J::parse(&conduction).unwrap();
        let radiation = retained.get("radiation").unwrap();
        let surface = &radiation.get("surfaces").unwrap().as_array().unwrap()[0];
        near(n(surface, "mean_temperature_k"), slab_oracle(reservoir), 1e-3);
        near(n(radiation, "radiative_out_w") + n(radiation, "convective_out_w"), 5.0, 2e-4);
        near(n(surface, "applied_heat_w"), n(surface, "nonlinear_heat_w"), 1e-8);
        assert_eq!(n(radiation, "radiative_out_w") < 0.0, reservoir == 350.0);
        let receipts = stage_receipt_hashes(&ledger, &run);
        assert_eq!(receipts.len(), 7);
        let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
        let html = String::from_utf8(artifact_bytes(&ledger, &receipt_str_field(&report, "report_html"))).unwrap();
        assert!(html.contains("radiative_out_w"));
        let qoi = J::parse(&qoi).unwrap();
        assert!(qoi.get("budget").unwrap().as_array().unwrap()[0]
            .get("terms").unwrap().as_array().unwrap().iter().any(|term|
                term.str_field("kind") == Some("model-form") && term.str_field("state") == Some("no-data")),
            "adaptive convergence cannot fill an unknown model-form error");
    }
}

#[test]
fn g1_fsim_radiation_adaptive_air_uses_only_convection_and_replays() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    let bytes = tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    declare(&mut spec, &card, "air", 293.15);
    adaptive_controls(&mut spec);
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes);
    let decoded = decode(&spec);
    let (_, qoi, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
    check_comparison(&conduction, true);
    let retained = J::parse(&conduction).unwrap();
    let radiation = retained.get("radiation").unwrap();
    let air = retained.get("conjugate").unwrap();
    let qrad = n(radiation, "radiative_out_w");
    assert!(qrad > 0.01, "fixture must distinguish radiative from convective heat");
    near(n(air, "air_total_w"), n(radiation, "convective_out_w"), 2e-4);
    near(n(air, "air_total_w") + qrad, 5.0, 2e-4);
    near(n(air, "mass_flow_kg_s") * 1007.0 * (n(air, "outlet_k") - n(air, "inlet_k")),
        n(air, "air_total_w"), 2e-4);
    let (_, replay_qoi, replay_conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
    assert_eq!(replay_conduction, conduction);
    assert_eq!(replay_qoi, qoi);
}

#[test]
fn g1_fsim_radiation_adaptive_retains_nonlinear_material_and_contact() {
    let base = contact_cards();
    let mut raw: Vec<RawCardPack> = base.iter()
        .map(|entry| raw_pack(entry.kind(), "contact-fixture-card", entry.bytes().to_vec())).collect();
    raw.push(raw_pack(CardPackKind::Material, "nonlinear-hot-body",
        material_pack_bytes_with_property("hot-nonlinear", "hot-nonlinear-radiation", 200.0, 450.0,
            fs_matdb::PropertyValue::Curve {
                abscissa: "T".to_string(),
                abscissa_dims: fs_project::spec::dims::TEMPERATURE,
                knots: vec![(200.0, 10.0), (450.0, 20.0)],
                dims: CONDUCTIVITY_DIMS,
            }, fs_matdb::InterpolationPolicy::LinearInside)));
    let base = CardPackSet::admit(raw).unwrap();
    let hot = base.materials().iter().find(|entry| entry.identity().starts_with("hot-nonlinear/")).unwrap();
    let (cards, card) = emissivity_cards(&base, 0.7);
    let mut spec = multi_region_contact_project();
    let binding = spec.materials.as_mut().unwrap().iter_mut().find(|b| b.region == "hot").unwrap();
    binding.card = hot.card().to_hex();
    binding.state = hot.identity().to_string();
    spec.requirements.as_mut().unwrap()[0].region = "hot".to_string();
    declare(&mut spec, &card, "hot", 293.15);
    adaptive_controls(&mut spec);
    let ledger = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&ledger, &spec);
    let (_, _, conduction, _) = run_conjugate_to_completion(&ledger, &decode(&spec), &cards);
    check_comparison(&conduction, false);
    assert!(conduction.contains("\"uses_nonlinear_jacobian\":true"), "{conduction}");
    let retained = J::parse(&conduction).unwrap();
    let fluxes = retained.path(&["interfaces", "fluxes"]).unwrap().as_array().unwrap();
    assert!(!fluxes.is_empty());
    assert!(n(&fluxes[0], "heat_rate_a_to_b_w").abs() > 0.01);
    let radiation = retained.get("radiation").unwrap();
    let energy = retained.get("energy").unwrap();
    near(n(radiation, "radiative_out_w") + n(radiation, "convective_out_w")
        - n(energy, "dirichlet_in_w"), 5.0, 2e-4);
}
