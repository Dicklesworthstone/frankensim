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
    nonlinear_contact_case(false);
}

#[test]
fn g1_fsim_natural_radiation_adaptive_retains_nonlinear_material_and_contact() {
    nonlinear_contact_case(true);
}

fn nonlinear_contact_case(natural: bool) {
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
    if natural { declare_natural(&mut spec, "hot", 293.15); }
    let ledger = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&ledger, &spec);
    let (_, _, conduction, _) = run_conjugate_to_completion(&ledger, &decode(&spec), &cards);
    check_comparison(&conduction, false);
    assert!(conduction.contains("\"uses_nonlinear_jacobian\":true"), "{conduction}");
    if natural { assert!(conduction.contains("\"natural_coefficient_feedback\":true"), "{conduction}"); }
    let retained = J::parse(&conduction).unwrap();
    let fluxes = retained.path(&["interfaces", "fluxes"]).unwrap().as_array().unwrap();
    assert!(!fluxes.is_empty());
    assert!(n(&fluxes[0], "heat_rate_a_to_b_w").abs() > 0.01);
    let radiation = retained.get("radiation").unwrap();
    let energy = retained.get("energy").unwrap();
    near(n(radiation, "radiative_out_w") + n(radiation, "convective_out_w")
        - n(energy, "dirichlet_in_w"), 5.0, 2e-4);
}

fn declare_natural(spec: &mut ProjectSpec, target: &str, ambient: f64) {
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let row = setup.boundaries.iter_mut().find(|row| row.target == target).unwrap();
    row.condition = fs_project::ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.0625, fs_project::spec::dims::LENGTH),
        ambient_temperature: QtyAny::new(ambient, fs_project::spec::dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".to_string(),
    };
    let envelope = spec.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = ambient;
    envelope.ambient_hi.value = ambient;
}

/// Independent scalar heat balance, with no FEM, production card evaluation
/// or radiative secant helper. All faces share one mean-temperature law, so
/// conservation fixes that mean even when the interior field is nonuniform.
fn passive_slab_oracle(reservoir: f64, pressure: f64) -> f64 {
    let ambient = 300.0;
    let length: f64 = 0.0625;
    let area = 0.0625;
    let (mut lo, mut hi) = (200.0, 450.0);
    for _ in 0..100 {
        let mean = 0.5*(lo+hi);
        let film = 0.5*(mean+ambient);
        let density = pressure/(287.05*film);
        let nu = 1.846e-5/density;
        let ra = 9.80665*(mean-ambient).abs()*length.powi(3)*0.707/(film*nu*nu);
        let nu_number = (0.825 + 0.387*ra.powf(1.0/6.0)
            / (1.0+(0.492_f64/0.707).powf(9.0/16.0)).powf(8.0/27.0)).powi(2);
        let convection = (0.0263/length)*nu_number*area*(mean-ambient);
        let radiation = 0.85*SIGMA*area*(mean.powi(4)-reservoir.powi(4));
        if convection+radiation > 5.0 { hi = mean; } else { lo = mean; }
    }
    0.5*(lo+hi)
}

#[test]
fn g1_fsim_natural_radiation_adaptive_slab_keeps_both_heat_signs_and_replays() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    for (reservoir, pressure) in [(270.0, 80000.0), (350.0, 101325.0)] {
        let bytes = slab_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().regions[0].seed =
            [0.03125, 0.0625, 0.0625].map(|v| QtyAny::new(v, fs_project::spec::dims::LENGTH));
        declare(&mut spec, &card, "air", reservoir);
        declare_natural(&mut spec, "air", 300.0);
        spec.envelope.as_mut().unwrap().pressure.value = pressure;
        adaptive_controls(&mut spec);
        let source = print_sexpr(&spec).unwrap();
        assert!(source.contains("natural-convection") && source.contains(":radiation"));
        let decoded = fs_project::parse_sexpr(&source).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (run, qoi, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
        check_comparison(&conduction, false);
        assert!(conduction.contains("\"natural_coefficient_feedback\":true"), "{conduction}");
        let retained = J::parse(&conduction).unwrap();
        let radiation = retained.get("radiation").unwrap();
        let surface = &radiation.get("surfaces").unwrap().as_array().unwrap()[0];
        let mean = n(surface, "mean_temperature_k");
        near(mean, passive_slab_oracle(reservoir, pressure), 1e-3);
        near(n(surface, "area_m2"), 0.0625, 1e-12);
        near(n(radiation, "radiative_out_w")+n(radiation, "convective_out_w"), 5.0, 2e-4);
        assert_eq!(mean < 300.0, reservoir < 300.0);
        assert_eq!(n(radiation, "convective_out_w") < 0.0, reservoir < 300.0);
        assert_eq!(n(radiation, "radiative_out_w") < 0.0, reservoir > 300.0);
        let receipts = stage_receipt_hashes(&ledger, &run);
        assert_eq!(receipts.len(), 7);
        let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
        let html_hash = receipt_str_field(&report, "report_html");
        let html = artifact_bytes(&ledger, &html_hash);
        let (_, replay_qoi, replay_conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
        assert_eq!(conduction, replay_conduction);
        assert_eq!(qoi, replay_qoi);
        assert_eq!(receipts, stage_receipt_hashes(&ledger, &run));
        assert_eq!(html, artifact_bytes(&ledger, &html_hash));
    }
}

#[test]
fn g0_tracked_natural_radiation_adaptive_project_is_admitted() {
    let source = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../data/reference-project/cooling-natural-radiation.fsim"));
    let decoded = fs_project::parse_sexpr(source).unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    assert!(source.contains("natural-convection") && source.contains(":radiation"));
    assert!(source.contains(":fidelity \"adaptive\""));
}
