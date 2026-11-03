//! Genuine two-stage time accuracy, source support and physical bath accounting.
use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA};

fn input(problem: ConductionProblem<'_>) -> EnthalpySdirk2Input<'_> {
    EnthalpySdirk2Input { problem, interfaces: None, radiation: EnthalpyRadiation::None, contacts: None }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e} (tol {tolerance:e})");
}

#[test]
fn sdirk2_spatial_eigenmode_has_second_order_accuracy_and_original_energy() {
    let mesh = tetrahedron(); let boundary = adiabatic(&mesh); let curve = sensible_curve();
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        // For this P1 tetra, [0,1,-1,0] is an exact discrete mode with
        // lambda=4*k/(rho0*Cp)=1.2/s, independently of the production assembly.
        let lambda = 1.2; let duration = 1.0; let exact = 20.0 * (-lambda * duration as f64).exp();
        let p = input(problem(&mesh, &boundary, &conductivity, &source));
        let mut errors = Vec::new();
        for count in [4, 8, 16] {
            let dt = duration / f64::from(count); let z = -lambda * dt;
            let factor = (1.0 + (1.0 - 2.0 * SDIRK2_GAMMA) * z) / (1.0 - SDIRK2_GAMMA * z).powi(2);
            let mut h = vec![350.0, 450.0, 250.0, 350.0];
            for step in 1..=count {
                let report = engine.advance_sdirk2(cx, [p; 2], &h, dt, config()).unwrap();
                close(report.first_stage.temperature[1] - 350.0,
                    (h[1] / 5.0 + 280.0 - 350.0) / (1.0 - SDIRK2_GAMMA * z), 2e-9);
                close(report.second_stage.temperature[1], 350.0 + 20.0 * factor.powi(step), 2e-9);
                close(report.second_stage.temperature[2], 350.0 - 20.0 * factor.powi(step), 2e-9);
                close(report.stored_energy_change_j, 0.0, 1e-9);
                close(report.source_j, 0.0, 0.0); close(report.dirichlet_in_j, 0.0, 0.0);
                close(report.energy_residual_j, 0.0, config().energy_tolerance_j);
                assert!(report.nodal_residual_j.iter().all(|r| r.abs() <= config().energy_tolerance_j));
                h = report.second_stage.specific_enthalpy_j_kg;
            }
            errors.push((h[1] / 5.0 + 280.0 - 350.0 - exact).abs());
        }
        for pair in errors.windows(2) { assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.3, "{errors:?}"); }
    });
}

#[test]
fn sdirk2_carry_can_leave_chart_support_while_both_real_stages_remain_admissible() {
    let mesh = tetrahedron(); let boundary = prescribed_wall(&mesh, 285.0);
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let zero = ScalarField::Uniform(0.0); let impossible = ScalarField::Uniform(-1e6);
    for datum in [0.0, 1_048_576.0] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0x20; 32]),
            SolidLiquidPhase::Solid, [(0.0, 280.0), (1100.0, 500.0)].map(|(h, t)|
                EnthalpyPhaseKnot { specific_enthalpy_j_kg: datum + h, temperature_k: t,
                    liquid_mass_fraction: 0.0, bulk_density_kg_m3: REFERENCE_DENSITY }).to_vec()).unwrap();
        with_cx(|cx| {
            let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
            // A 2^20 J/kg reference has an h ULP near 2e-10 J/kg. The absolute
            // Newton target must resolve physical joules, not demand an h
            // correction smaller than that represented state can retain.
            let mut policy = config();
            policy.newton.absolute_tolerance = 1e-10;
            policy.newton.relative_tolerance = 0.0;
            let old = [datum + 25.0, datum + 100.0, datum + 25.0, datum + 25.0];
            let p = input(problem(&mesh, &boundary, &conductivity, &zero));
            let dt = 5.0 / 1.2;
            let step = engine.advance_sdirk2(cx, [p; 2], &old, dt, policy).unwrap();
            let first = step.first_stage.specific_enthalpy_j_kg[1] - datum;
            let final_h = step.second_stage.specific_enthalpy_j_kg[1] - datum;
            let false_history = 100.0 + (1.0 - SDIRK2_GAMMA) / SDIRK2_GAMMA * (first - 100.0);
            assert!(false_history < 0.0);
            assert!(first > 0.0 && final_h > 0.0);
            let z = -5.0;
            let factor = (1.0 + (1.0 - 2.0 * SDIRK2_GAMMA) * z) / (1.0 - SDIRK2_GAMMA * z).powi(2);
            close(first, 25.0 + 75.0 / (1.0 - SDIRK2_GAMMA * z), 2e-8);
            close(final_h, 25.0 + 75.0 * factor, 2e-8);
            close(step.dirichlet_in_j, engine.reference_nodal_masses_kg()[1] * (final_h - 100.0), 2e-9);
            // L-stability can overshoot equilibrium, so no monotonicity claim
            // is inferred; it still must evaluate only real in-domain stages.
            assert!(final_h < 25.0);
            let bad = input(problem(&mesh, &boundary, &conductivity, &impossible));
            assert!(engine.advance_sdirk2(cx, [p, bad], &old, dt, policy).is_err());
            with_cancelled_cx(|cancelled| assert!(matches!(
                engine.advance_sdirk2(cancelled, [p; 2], &old, dt, policy),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. })))));
            let replay = engine.advance_sdirk2(cx, [p; 2], &old, dt, policy).unwrap();
            assert_eq!(bits(&step.second_stage.specific_enthalpy_j_kg), bits(&replay.second_stage.specific_enthalpy_j_kg));
        });
    }
}

#[test]
fn sdirk2_prescribed_stage_storage_and_sources_are_integrated_once() {
    let mesh = tetrahedron(); let curve = sensible_curve();
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let boundary = [340.0, 360.0].map(|t| ThermalBoundaryBuilder::new(&mesh)
        .region("bath", |_| true, ThermalBc::dirichlet(t).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap());
    let sources = [ScalarField::Uniform(100.0), ScalarField::Uniform(-50.0)];
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        let stages = [0, 1].map(|i| input(problem(&mesh, &boundary[i], &conductivity, &sources[i])));
        let dt = 0.2;
        let result = engine.advance_sdirk2(cx, stages, &[150.0; 4], dt, config()).unwrap();
        assert_eq!(result.first_stage.temperature, [340.0; 4]);
        assert_eq!(result.second_stage.temperature, [360.0; 4]);
        let energy = REFERENCE_DENSITY * VOLUME * 5.0 * (360.0 - 310.0);
        let source = VOLUME * dt * ((1.0 - SDIRK2_GAMMA) * 100.0 - SDIRK2_GAMMA * 50.0);
        close(result.stored_energy_change_j, energy, 1e-11);
        close(result.source_j, source, 1e-11);
        close(result.dirichlet_in_j, energy - source, 1e-11);
        close(result.energy_residual_j, 0.0, 1e-11);
        assert_eq!(result.first_stage.source_w, 100.0 * VOLUME);
        assert_eq!(result.second_stage.source_w, -50.0 * VOLUME);
    });
}
