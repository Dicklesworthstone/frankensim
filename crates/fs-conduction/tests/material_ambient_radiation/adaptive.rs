//! Independent continuous-time latent radiation, with a changing emissivity.

use super::*;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus,
};
use std::ops::ControlFlow;

#[test]
fn phase_error_controls_adaptive_emissivity_at_constant_temperature() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let card = card(Some("melt"));
    let context = context();
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &conductivity,
        source: &source,
        element_materials: None,
    };
    for datum in [0.0, 8192.0] {
        let curve = curve(&card, datum);
        for bath in [300.0_f64, 500.0] {
            let coefficient = AREA * STEFAN_BOLTZMANN_W_M2_K4
                * (400.0_f64.powi(4) - bath.powi(4));
            // Three radiating nodes have one kg each. All T stay at 400 K,
            // so m du/dt = -C (0.2 + 0.6 u/L) has this exact exponential.
            let rate = coefficient * 0.6 / LATENT / 3.0;
            let exact = 500.0 + (500.0 + 0.2 * LATENT / 0.6) * (-rate).exp_m1();
            let mut errors = Vec::new();
            let mut counts = Vec::new();
            with_cx(|cx| {
                let radiation = MaterialPointwiseAmbientRadiation::new(
                    cx, &mesh, &boundary, &[patch(&card, &curve, &context, bath)],
                ).unwrap();
                let engine = EnthalpyBackwardEuler::uniform(
                    cx, &mesh, &curve, RHO, budget(),
                ).unwrap();
                for phase_tolerance in [1e-3, 1e-5] {
                    let controls = AdaptiveEnthalpyConfig {
                        enthalpy_absolute_tolerance_j_kg: LATENT,
                        temperature_absolute_tolerance_k: 100.0,
                        liquid_fraction_absolute_tolerance: phase_tolerance,
                        relative_tolerance: 0.0,
                        minimum_step_s: 1e-8,
                        maximum_step_s: 0.5,
                        endpoint: config(),
                    };
                    let initial = [datum + 500.0; 4];
                    let mut state = engine.adaptive_state(
                        cx, &initial, 0.0, 0.5, controls,
                    ).unwrap();
                    let mut integrated = 0.0;
                    let mut previous_time = 0.0_f64;
                    let mut visits = 0;
                    let report = state.advance_to_with_material_ambient_radiation_observed(
                        cx, problem, None, &radiation, 1.0, 4096, |trial| {
                            assert_eq!(trial.start_time_s.to_bits(), previous_time.to_bits());
                            assert!(trial.error.normalized <= 1.0);
                            assert_eq!(trial.error.maximum_temperature_difference_k, 0.0);
                            close(trial.error.normalized,
                                trial.error.maximum_liquid_fraction_difference / phase_tolerance,
                                2e-10);
                            for (half, dt) in [
                                (trial.first_half, trial.midpoint_time_s - trial.start_time_s),
                                (trial.second_half, trial.end_time_s - trial.midpoint_time_s),
                            ] {
                                assert_eq!(half.temperature, vec![400.0; 4]);
                                let radiation = half.pointwise_radiation.as_ref().unwrap();
                                assert_eq!(radiation.patches[0].material_states.len(), 3);
                                integrated -= dt * radiation.outward_heat_w;
                            }
                            visits += 1;
                            previous_time = trial.end_time_s;
                            ControlFlow::Continue(())
                        },
                    ).unwrap();
                    assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
                    assert_eq!(state.time_s().to_bits(), 1.0_f64.to_bits());
                    assert_eq!(visits, report.accepted);
                    assert_eq!(report.endpoint_solves, 3 * (report.accepted + report.rejected));
                    close(report.integrated_net_input_j, integrated, 2e-9);
                    close(report.stored_energy_change_j, integrated, 2e-7);
                    close(state.specific_enthalpy_j_kg()[3], initial[3], 2e-9);
                    let endpoint = state.specific_enthalpy_j_kg()[0] - datum;
                    errors.push((endpoint - exact).abs());
                    counts.push(report.accepted);
                }
            });
            assert!(errors[1] < 0.75, "fine trajectory error: {errors:?}");
            assert!(errors[1] < 0.25 * errors[0], "trajectory refinement: {errors:?}");
            assert!(counts[1] > counts[0]);
        }
    }
}
