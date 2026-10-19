//! Native vacuum-facing traces: real project admission, meshing, solve and reports.
//! Reuse the existing radiation/card fixtures and independent balance equations.

use super::*;
use fs_project::spec::dims;

fn vacuum(spec: &mut ProjectSpec) {
    let cooling = spec.cooling.as_mut().unwrap();
    cooling.fans.clear();
    cooling.vents.clear();
    cooling.fan_system = None;
    cooling.airflow_leakage = None;
    cooling.conduction.as_mut().unwrap().boundaries[0].condition =
        ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
}

fn setup(spec: &mut ProjectSpec) -> &mut ConductionSetup {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
}

#[test]
fn g1_native_vacuum_steady_balance_and_no_singular_off_comparison() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    for reservoir in [270.0_f64, 293.15, 350.0] {
        let bytes = slab_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        setup(&mut spec).regions[0].seed = [0.03125, 0.0625, 0.0625]
            .map(|value| QtyAny::new(value, dims::LENGTH));
        declare(&mut spec, &card, "air", reservoir);
        vacuum(&mut spec);
        let decoded = fs_project::parse_sexpr(&print_sexpr(&spec).unwrap()).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (run, qoi, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &cards);
        let retained = J::parse(&conduction).unwrap();
        let radiation = retained.get("radiation").unwrap();
        let patch = &radiation.get("surfaces").unwrap().as_array().unwrap()[0];
        // Global area-mean balance, independent of FEM stiffness or radiation code.
        let expected = (reservoir.powi(4) + 5.0 / (0.85 * SIGMA * 0.0625)).sqrt().sqrt();
        near(n(patch, "area_m2"), 0.0625, 1e-12);
        near(n(patch, "mean_temperature_k"), expected, 2e-5);
        near(n(radiation, "radiative_out_w"), 5.0, 2e-6);
        assert_eq!(n(radiation, "convective_out_w"), 0.0);
        assert_eq!(n(radiation, "radiative_magnitude_share"), 1.0);
        near(n(patch, "applied_heat_w"), n(patch, "nonlinear_heat_w"), 1e-7);
        let qoi = J::parse(&qoi).unwrap();
        let comparison = qoi.path(&["model_form_sensitivity", "radiation"]).unwrap();
        assert_eq!(comparison.str_field("state"), Some("no-data"));
        assert!(comparison.str_field("reason").unwrap().contains("no temperature anchor"));
        assert!(comparison.get("radiation_off_k").is_none());
        let receipts = stage_receipt_hashes(&ledger, &run);
        assert_eq!(receipts.len(), 7);
        let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
        let html = String::from_utf8(artifact_bytes(&ledger,
            &receipt_str_field(&report, "report_html"))).unwrap();
        assert!(html.contains("radiative_out_w"));
    }
}

#[test]
fn g0_zero_flux_needs_its_matched_sourced_radiation_and_cannot_erase_other_laws() {
    let (_, card) = emissivity_cards(&fixture_cards(), 0.85);
    let mut original = conduction_fixture_project(7, &tetra_stl());
    declare(&mut original, &card, "air", 293.15);
    vacuum(&mut original);
    assert!(original.validate().is_empty());
    let has = |spec: &ProjectSpec, code| spec.validate().iter().any(|v| v.code == code);

    let mut spec = original.clone();
    setup(&mut spec).radiation = None;
    assert!(has(&spec, "project-conduction-boundary-unanchored"));
    let mut spec = original.clone();
    setup(&mut spec).radiation.as_mut().unwrap().surfaces.clear();
    assert!(has(&spec, "project-radiation-surfaces-empty"));
    assert!(has(&spec, "project-conduction-boundary-unanchored"));
    let mut spec = original.clone();
    setup(&mut spec).radiation.as_mut().unwrap().surfaces[0].target = "missing".into();
    assert!(has(&spec, "project-radiation-target-boundary"));
    assert!(has(&spec, "project-conduction-boundary-unanchored"));

    for value in [-1.0, 1.0, f64::NAN] {
        let mut spec = original.clone();
        setup(&mut spec).boundaries[0].condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(value, dims::HEAT_FLUX),
        };
        assert!(has(&spec, "project-radiation-target-boundary"));
        assert!(has(&spec, "project-conduction-boundary-unanchored"));
    }
    let mut spec = original.clone();
    setup(&mut spec).boundaries[0].condition = ThermalBoundaryCondition::HeatFlux {
        outward_flux: QtyAny::new(0.0, fs_qty::Dims::NONE),
    };
    assert!(has(&spec, "project-conduction-boundary-dims"));
    assert!(has(&spec, "project-radiation-target-boundary"));
    let mut spec = original.clone();
    setup(&mut spec).boundaries[0].condition = ThermalBoundaryCondition::FixedTemperature {
        temperature: QtyAny::new(300.0, dims::TEMPERATURE),
    };
    assert!(has(&spec, "project-radiation-target-boundary"));
    let mut spec = original;
    setup(&mut spec).boundaries[0].condition = ThermalBoundaryCondition::HeatFlux {
        outward_flux: QtyAny::new(-0.0, dims::HEAT_FLUX),
    };
    assert!(spec.validate().is_empty(), "signed zero is still explicitly zero exchange");
}

#[test]
fn g1_native_vacuum_transient_matches_each_implicit_radiative_endpoint() {
    let (cards, card) = emissivity_cards(&fixture_cards(), 0.85);
    let area = 8.0 * 3.0_f64.sqrt();
    let capacity = 100.0 * 8.0 / 3.0;
    for reservoir in [270.0_f64, 350.0] {
        let bytes = super::super::transient_radiation::regular_tetra_stl();
        let mut spec = conduction_fixture_project(7, &bytes);
        super::super::transient_product::declaration(&mut spec, 100.0, 1.0);
        setup(&mut spec).regions[0].seed = [QtyAny::new(0.0, dims::LENGTH); 3];
        setup(&mut spec).transient.as_mut().unwrap().initial_temperature =
            fs_project::TransientInitialTemperature::Uniform(QtyAny::new(330.0, dims::TEMPERATURE));
        declare(&mut spec, &card, "air", reservoir);
        vacuum(&mut spec);
        let decoded = fs_project::parse_sexpr(&print_sexpr(&spec).unwrap()).unwrap();
        assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, bytes);
        let (_, qoi, conduction, field) = run_conjugate_to_completion(&ledger, &decoded, &cards);
        let retained = J::parse(&conduction).unwrap();
        let time = retained.get("transient").unwrap();
        assert_eq!(n(time, "total_steps"), 6.0);
        for grid in ["coarse", "fine"] {
            let mut old = 330.0;
            let mut old_actual = 330.0;
            for row in time.get(grid).unwrap().as_array().unwrap() {
                let dt = n(row, "dt_s");
                let (mut low, mut high) = (200.0_f64, 450.0_f64);
                for _ in 0..100 {
                    let t = 0.5 * (low + high);
                    let residual = capacity * (t - old) / dt
                        + 0.85 * SIGMA * area * (t.powi(4) - reservoir.powi(4)) - 5.0;
                    if residual > 0.0 { high = t; } else { low = t; }
                }
                let expected = 0.5 * (low + high);
                let actual = n(row, "final_region_max_k");
                near(actual, expected, 2e-8);
                let radiation = row.get("radiation").unwrap();
                assert_eq!(n(radiation, "convective_out_w"), 0.0);
                let heat = 0.85 * SIGMA * area * (actual.powi(4) - reservoir.powi(4));
                near(n(radiation, "nonlinear_radiative_out_w"), heat, 2e-6);
                near(n(row, "stored_energy_change_j"), capacity * (actual - old_actual), 2e-6);
                near(capacity * (actual - old_actual) - dt * (5.0 - heat), 0.0, 2e-6);
                assert!(n(radiation, "physical_residual_norm_j")
                    <= n(radiation, "physical_residual_tolerance_j"));
                assert!(n(radiation, "physical_energy_residual_j").abs() <= 1e-6);
                assert_eq!(heat > 0.0, reservoir < 330.0);
                old = expected;
                old_actual = actual;
            }
        }
        let last = time.get("fine").unwrap().as_array().unwrap().last().unwrap();
        let final_k = n(last, "final_region_max_k");
        let solution = J::parse(&field).unwrap();
        let temperatures = solution.get("temperature").unwrap().as_array().unwrap();
        assert_eq!(temperatures.len(), 4, "regular tetra uniform-mode oracle");
        for value in temperatures { near(value.as_f64().unwrap(), final_k, 2e-8); }
        let qoi = J::parse(&qoi).unwrap();
        let comparison = qoi.path(&["model_form_sensitivity", "radiation"]).unwrap();
        assert_eq!(comparison.str_field("state"), Some("measured"),
            "without radiation, finite-time capacity storage still defines a unique endpoint");
        near(n(comparison, "radiation_off_k"), 330.0 + 5.0 / capacity, 2e-8);
    }
}
