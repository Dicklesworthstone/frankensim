//! Continuous latent exchange with changing reflected irradiation.
use super::*;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus,
};
use std::ops::ControlFlow;

#[test]
fn phase_error_refines_reflected_exchange_while_both_material_temperatures_stay_fixed() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let context = context();
    let cards = [material_card("left"), material_card("right")];
    let enclosure = closed(&mesh, &cards, &context);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, source: &source, element_materials: None };
    let driving = AREA * SIGMA * (400.0_f64.powi(4) - 350.0_f64.powi(4));
    let slope = 0.6 / LATENT;
    // Three kg radiate on each side. Let h_left=500-x, h_right=500+x.
    // Then eps_left=.5-k*x, eps_right=.5+k*x, dx/dt=Q(x)/3. Integrating
    // the reciprocal gray resistance gives this independent exact time law.
    let elapsed = |x: f64| 3.0 / driving
        * (2.0 / slope * (2.0 * slope * x).atanh() - x);
    let (mut lo, mut hi) = (0.0, 499.0);
    for _ in 0..80 {
        let x = (lo + hi) / 2.0;
        if elapsed(x) > 5.0 { hi = x; } else { lo = x; }
    }
    let exact = (lo + hi) / 2.0;
    for datum in [0.0, 8192.0] {
        let curves = [phase_curve(&cards[0], 400.0, datum),
            phase_curve(&cards[1], 350.0, datum)];
        let declarations = bindings(&cards, &curves, &context);
        let mut errors = Vec::new();
        let mut counts = Vec::new();
        with_cx(|cx| {
            let base = PointwiseEnclosureRadiation::new(cx, &mesh, &boundary,
                &enclosure, PointwiseEnclosureTolerance::default()).unwrap();
            let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
            let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
                &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: RHO },
                    ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: RHO }],
                &[0, 1], budget()).unwrap();
            for phase_tolerance in [1e-3, 1e-5] {
                let controls = AdaptiveEnthalpyConfig {
                    enthalpy_absolute_tolerance_j_kg: LATENT,
                    temperature_absolute_tolerance_k: 100.0,
                    liquid_fraction_absolute_tolerance: phase_tolerance,
                    relative_tolerance: 0.0, minimum_step_s: 1e-8,
                    maximum_step_s: 5.0, endpoint: config(),
                };
                let initial = [datum + 500.0; 8];
                let mut state = engine.adaptive_state(cx, &initial, 0.0, 5.0, controls).unwrap();
                let mut left_out_j = 0.0;
                let mut visits = 0;
                let mut previous_time = 0.0_f64;
                let report = state.advance_to_with_material_enclosure_radiation_observed(
                    cx, problem, None, &owner, 5.0, 4096, |trial| {
                        assert_eq!(trial.start_time_s.to_bits(), previous_time.to_bits());
                        assert_eq!(trial.error.maximum_temperature_difference_k, 0.0);
                        close(trial.error.normalized,
                            trial.error.maximum_liquid_fraction_difference / phase_tolerance, 2e-9);
                        for (half, dt) in [
                            (trial.first_half, trial.midpoint_time_s - trial.start_time_s),
                            (trial.second_half, trial.end_time_s - trial.midpoint_time_s),
                        ] {
                            assert_eq!(&half.temperature[..4], &[400.0; 4]);
                            assert_eq!(&half.temperature[4..], &[350.0; 4]);
                            let radiation = half.pointwise_enclosure_radiation.as_ref().unwrap();
                            assert_eq!(radiation.patches[0].material_states.len(), 3);
                            assert_eq!(radiation.patches[1].material_states.len(), 3);
                            assert!(radiation.patches[0].outward_heat_w > 0.0);
                            close(radiation.outward_heat_w, 0.0, 1e-8);
                            left_out_j += dt * radiation.patches[0].outward_heat_w;
                        }
                        visits += 1;
                        previous_time = trial.end_time_s;
                        ControlFlow::Continue(())
                    },
                ).unwrap();
                assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
                assert_eq!(state.time_s(), 5.0);
                assert_eq!(visits, report.accepted);
                close(report.integrated_net_input_j, 0.0, 1e-8);
                close(report.stored_energy_change_j, 0.0, 2e-7);
                let transfer = initial[0] - state.specific_enthalpy_j_kg()[0];
                close(left_out_j, 3.0 * transfer, 2e-7);
                close(state.specific_enthalpy_j_kg()[4] - initial[4], transfer, 2e-7);
                close(state.specific_enthalpy_j_kg()[3], initial[3], 2e-9);
                close(state.specific_enthalpy_j_kg()[7], initial[7], 2e-9);
                errors.push((transfer - exact).abs());
                counts.push(report.accepted);
            }
        });
        assert!(errors[1] < 0.3, "fine trajectory error: {errors:?}");
        assert!(errors[1] < 0.25 * errors[0], "reflection refinement: {errors:?}");
        assert!(counts[1] > counts[0]);
    }
}
