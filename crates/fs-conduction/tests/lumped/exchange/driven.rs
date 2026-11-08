//! Original sensible and latent inventory under held convection and signed power.

use super::*;
use fs_conduction::lumped::solve_lumped_enthalpy_exact;
use fs_conduction::lumped::exchange::{
    LumpedDrivenExchange, LumpedDrivenExchangeConfig, solve_lumped_driven_exchange,
    solve_lumped_driven_exchange_response,
};

fn driven_config(
    initial: f64,
    bath: f64,
    power: f64,
    duration: f64,
) -> LumpedDrivenExchangeConfig {
    LumpedDrivenExchangeConfig {
        initial_specific_enthalpy_j_kg: initial,
        reservoir_temperature_k: bath,
        applied_power_w: power,
        duration_s: duration,
        energy_tolerance_j: 1e-7,
    }
}

fn balanced(result: &LumpedDrivenExchange, tolerance: f64) {
    near(
        result.source_heat_j + result.reservoir_heat_j,
        result.heat_into_body_j,
        tolerance,
    );
    near(
        result.source_heat_j + result.reservoir_heat_j,
        result.stored_enthalpy_change_j,
        tolerance,
    );
    assert!(result.heat_residual_j.abs() <= tolerance);
    assert!(result.heat_partition_residual_j.abs() <= tolerance);
    assert!(result.energy_balance_residual_j.abs() <= tolerance);
}

#[test]
fn g1_signed_power_crosses_original_fusion_and_keeps_both_heat_accounts() {
    let curve = phase_curve();
    let body = body(&curve, 100.0, 0.0); // G=1 W/K, mass=2 kg.
    // Integrate Cp/(Teq-T) separately on each sensible branch and L/(Teq-Tm)
    // on the original latent plateau. Teq is used only for this analytic oracle.
    let heating =
        100.0 * (500.0_f64 / 300.0).ln() + 25_000.0 / 300.0 + 200.0 * (300.0_f64 / 200.0).ln();
    let cooling =
        200.0 * (1_000.0_f64 / 800.0).ln() + 25_000.0 / 800.0 + 100.0 * (800.0_f64 / 600.0).ln();
    with_cx(|cx| {
        for (h0, bath, power, exposure, target_h, target_t, liquid) in [
            (10_000.0, 300.0, 600.0, heating, 75_000.0, 700.0, 1.0),
            (95_000.0, 800.0, -1_000.0, cooling, 10_000.0, 400.0, 0.0),
        ] {
            let cfg = driven_config(h0, bath, power, 2.0 * exposure);
            let result =
                solve_lumped_driven_exchange(cx, &body, BiotGate::corpus_default(), cfg).unwrap();
            let expected_heat = 2.0 * (target_h - h0);
            let expected_mean = bath + power - expected_heat / cfg.duration_s;
            assert_eq!(result.body_identity, body.identity());
            near(result.caloric.state.temperature_k(), target_t, 2e-10);
            near(result.caloric.state.specific_enthalpy_j_kg(), target_h, 2e-8);
            assert_eq!(result.caloric.state.liquid_mass_fraction(), liquid);
            near(result.caloric.mean_temperature_k, expected_mean, 2e-10);
            near(result.heat_into_body_j, expected_heat, 1e-7);
            assert_eq!(result.source_heat_j, power * cfg.duration_s);
            near(
                result.reservoir_heat_j,
                cfg.duration_s * (bath - expected_mean),
                1e-7,
            );
            assert!(result.source_heat_j * result.reservoir_heat_j < 0.0);
            assert!(result.maximum_biot < LUMPED_BIOT_CEILING);
            balanced(&result, cfg.energy_tolerance_j);

            // Both half windows join inside fusion. The second retains the
            // actual latent inventory instead of restarting from Tm alone.
            let first = solve_lumped_driven_exchange(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedDrivenExchangeConfig { duration_s: cfg.duration_s * 0.5, ..cfg },
            )
            .unwrap();
            assert_eq!(first.caloric.state.temperature_k(), 600.0);
            assert!(first.caloric.state.liquid_mass_fraction() > 0.0);
            assert!(first.caloric.state.liquid_mass_fraction() < 1.0);
            let second = solve_lumped_driven_exchange(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedDrivenExchangeConfig {
                    initial_specific_enthalpy_j_kg: first.caloric.state.specific_enthalpy_j_kg(),
                    duration_s: cfg.duration_s * 0.5,
                    ..cfg
                },
            )
            .unwrap();
            near(second.caloric.state.temperature_k(), target_t, 2e-10);
            near(first.heat_into_body_j + second.heat_into_body_j, expected_heat, 1e-7);
            near(
                first.reservoir_heat_j + second.reservoir_heat_j,
                result.reservoir_heat_j,
                1e-7,
            );
        }
    });
}

#[test]
fn g1_isolated_signed_source_is_exact_enthalpy_calorimetry_through_fusion() {
    let curve = phase_curve();
    let body = body(&curve, 0.0, 0.0);
    with_cx(|cx| {
        for (h0, target_h, target_t, power, duration, expected_mean) in [
            (
                10_000.0, 75_000.0, 700.0, 1_000.0, 130.0,
                (20_000.0 * 500.0 + 25_000.0 * 600.0 + 20_000.0 * 650.0) / 65_000.0,
            ),
            (
                95_000.0, 10_000.0, 400.0, -1_000.0, 170.0,
                (40_000.0 * 700.0 + 25_000.0 * 600.0 + 20_000.0 * 500.0) / 85_000.0,
            ),
        ] {
            // dh/dt=P/m makes the exact time mean equal to integral(T dh)/Delta h.
            let cfg = driven_config(h0, 300.0, power, duration);
            let response = solve_lumped_driven_exchange_response(
                cx, &body, BiotGate::corpus_default(), cfg,
            )
            .unwrap();
            let result = response.exchange;
            near(result.caloric.state.specific_enthalpy_j_kg(), target_h, 2e-8);
            near(result.caloric.state.temperature_k(), target_t, 2e-10);
            near(result.caloric.mean_temperature_k, expected_mean, 2e-10);
            assert_eq!(result.heat_into_body_j, power * duration);
            assert_eq!(result.source_heat_j, result.heat_into_body_j);
            assert_eq!(result.reservoir_heat_j, 0.0);
            assert_eq!(response.heat_initial_enthalpy_derivative_kg, 0.0);
            assert_eq!(response.heat_reservoir_derivative_j_k, 0.0);
            assert_eq!(response.heat_source_power_derivative_s, duration);
            assert_eq!(response.heat_duration_derivative_w, power);
            balanced(&result, cfg.energy_tolerance_j);
            assert_eq!(
                result,
                solve_lumped_driven_exchange(cx, &body, BiotGate::corpus_default(), cfg).unwrap(),
            );
        }
    });
}

#[test]
fn g1_latent_inventory_and_stationary_through_flow_preserve_signed_source_work() {
    let curve = phase_curve();
    let body = body(&curve, 100.0, 0.0);
    with_cx(|cx| {
        for (initial, target, power, fraction) in [
            (35_000.0, 45_000.0, 200.0, 0.6),
            (45_000.0, 35_000.0, -200.0, 0.2),
        ] {
            let cfg = driven_config(initial, 600.0, power, 100.0);
            let response = solve_lumped_driven_exchange_response(
                cx, &body, BiotGate::corpus_default(), cfg,
            )
            .unwrap();
            let result = response.exchange;
            assert_eq!(result.caloric.state.temperature_k(), 600.0);
            assert_eq!(result.caloric.mean_temperature_k, 600.0);
            assert_eq!(result.caloric.state.specific_enthalpy_j_kg(), target);
            near(result.caloric.state.liquid_mass_fraction(), fraction, 2e-15);
            assert_eq!(result.heat_into_body_j, power * 100.0);
            assert_eq!(result.reservoir_heat_j, 0.0);
            assert_eq!(response.heat_initial_enthalpy_derivative_kg, 0.0);
            assert_eq!(response.heat_reservoir_derivative_j_k, 100.0);
            assert_eq!(response.heat_source_power_derivative_s, 100.0);
            assert_eq!(response.heat_duration_derivative_w, power);
            balanced(&result, cfg.energy_tolerance_j);
        }
        for (bath, power) in [(300.0, 300.0), (800.0, -200.0)] {
            let cfg = driven_config(40_000.0, bath, power, 100.0);
            let response = solve_lumped_driven_exchange_response(
                cx, &body, BiotGate::corpus_default(), cfg,
            )
            .unwrap();
            let result = response.exchange;
            assert_eq!(result.caloric.state, result.caloric.initial);
            assert_eq!(result.heat_into_body_j, 0.0);
            assert_eq!(result.source_heat_j, 100.0 * power);
            assert_eq!(result.reservoir_heat_j, -result.source_heat_j);
            assert_eq!(response.heat_source_power_derivative_s, 100.0);
            assert_eq!(response.heat_reservoir_derivative_j_k, 100.0);
            assert_eq!(response.heat_duration_derivative_w, 0.0);
            balanced(&result, cfg.energy_tolerance_j);
        }
        let zero = solve_lumped_driven_exchange_response(
            cx, &body, BiotGate::corpus_default(), driven_config(10_000.0, 300.0, 600.0, 0.0),
        )
        .unwrap();
        assert_eq!(zero.exchange.caloric.state, zero.exchange.caloric.initial);
        assert_eq!(zero.exchange.source_heat_j, 0.0);
        assert_eq!(zero.heat_source_power_derivative_s, 0.0);
        assert_eq!(zero.heat_reservoir_derivative_j_k, 0.0);
        assert_eq!(zero.heat_duration_derivative_w, 500.0);
    });
}

#[test]
fn g1_source_balanced_at_melting_never_spends_unavailable_latent_heat() {
    let curve = phase_curve();
    let body = body(&curve, 100.0, 0.0);
    with_cx(|cx| {
        for duration in [1_000.0, 100_000.0] {
            // P+G*(Tb-T)=600-T: fusion is approached from below. Even when
            // rounded h reaches its corner, no latent heat has been supplied.
            let cfg = driven_config(10_000.0, 500.0, 100.0, duration);
            let result =
                solve_lumped_driven_exchange(cx, &body, BiotGate::corpus_default(), cfg).unwrap();
            let remaining = (-duration / 200.0).exp();
            let expected_heat = -40_000.0 * (-duration / 200.0).exp_m1();
            near(result.caloric.state.temperature_k(), 600.0 - 200.0 * remaining, 2e-10);
            near(
                result.caloric.state.specific_enthalpy_j_kg(),
                30_000.0 - 20_000.0 * remaining,
                2e-8,
            );
            assert_eq!(result.caloric.state.liquid_mass_fraction(), 0.0);
            near(result.heat_into_body_j, expected_heat, 1e-7);
            near(result.reservoir_heat_j, expected_heat - 100.0 * duration, 1e-7);
            balanced(&result, cfg.energy_tolerance_j);
        }
    });
}

#[test]
fn g1_cross_phase_source_and_boundary_jacobians_match_independent_flow_oracles() {
    let curve = phase_curve();
    let body = body(&curve, 250.0, 0.0); // G=2.5 W/K checks physical unit scaling.
    let conductance = 2.5;
    let heating =
        100.0 * (500.0_f64 / 300.0).ln() + 25_000.0 / 300.0 + 200.0 * (300.0_f64 / 200.0).ln();
    let cooling =
        200.0 * (975.0_f64 / 800.0).ln() + 25_000.0 / 800.0 + 100.0 * (800.0_f64 / 600.0).ln();
    with_cx(|cx| {
        for (initial, bath, power, exposure, target_t, heating) in [
            (10_000.0, 300.0, 1_500.0, heating, 700.0, true),
            (90_000.0, 800.0, -2_500.0, cooling, 400.0, false),
        ] {
            let duration = 2.0 * exposure / conductance;
            let cfg = driven_config(initial, bath, power, duration);
            let response = solve_lumped_driven_exchange_response(
                cx, &body, BiotGate::corpus_default(), cfg,
            )
            .unwrap();
            let reference = |h: f64, tb: f64, p: f64, t: f64| {
                three_phase_heat(h, tb + p / conductance, conductance * t, heating)
            };
            near(response.exchange.heat_into_body_j, reference(initial, bath, power, duration), 1e-7);
            near(response.exchange.caloric.state.temperature_k(), target_t, 2e-10);
            let derivative_h = (reference(initial + 0.01, bath, power, duration)
                - reference(initial - 0.01, bath, power, duration)) / 0.02;
            let derivative_bath = (reference(initial, bath + 0.001, power, duration)
                - reference(initial, bath - 0.001, power, duration)) / 0.002;
            let derivative_power = (reference(initial, bath, power + 0.01, duration)
                - reference(initial, bath, power - 0.01, duration)) / 0.02;
            let derivative_duration = (reference(initial, bath, power, duration + 0.001)
                - reference(initial, bath, power, duration - 0.001)) / 0.002;
            for (actual, expected) in [
                (response.heat_initial_enthalpy_derivative_kg, derivative_h),
                (response.heat_reservoir_derivative_j_k, derivative_bath),
                (response.heat_source_power_derivative_s, derivative_power),
                (response.heat_duration_derivative_w, derivative_duration),
            ] {
                near(actual, expected, 2e-6 * expected.abs().max(1.0));
            }
            near(
                response.heat_reservoir_derivative_j_k,
                conductance * response.heat_source_power_derivative_s,
                1e-10,
            );
            near(
                response.heat_duration_derivative_w,
                power + conductance * (bath - target_t),
                2e-8,
            );
            balanced(&response.exchange, cfg.energy_tolerance_j);
        }
    });
}

#[test]
fn g0_zero_drive_phase_corner_has_a_value_but_no_two_sided_source_response() {
    let curve = phase_curve();
    let coupled = body(&curve, 100.0, 0.0);
    let isolated = body(&curve, 0.0, 0.0);
    with_cx(|cx| {
        for (body, power) in [(&coupled, 300.0), (&isolated, 0.0)] {
            let cfg = driven_config(30_000.0, 300.0, power, 100.0);
            let result =
                solve_lumped_driven_exchange(cx, body, BiotGate::corpus_default(), cfg).unwrap();
            assert_eq!(result.caloric.state, result.caloric.initial);
            assert_eq!(result.caloric.state.liquid_mass_fraction(), 0.0);
            assert_eq!(result.heat_into_body_j, 0.0);
            assert_eq!(result.source_heat_j, power * cfg.duration_s);
            assert_eq!(result.reservoir_heat_j, -result.source_heat_j);
            balanced(&result, cfg.energy_tolerance_j);
            // At h=30000, the one-sided T(h) slopes are 1/100 and zero.
            // Even with G=0, a source perturbation has two different mean-T
            // derivatives. Value-only calorimetry must still be accepted.
            assert!(
                solve_lumped_driven_exchange_response(cx, body, BiotGate::corpus_default(), cfg)
                    .is_err()
            );
        }
    });
}

#[test]
fn g0_driven_exchange_preserves_source_support_biot_and_native_energy_admission() {
    let curve = phase_curve();
    let good = body(&curve, 100.0, 0.0);
    let isolated = body(&curve, 0.0, 0.0);
    let radiative = body(&curve, 100.0, 0.2);
    let high_biot = LumpedEnthalpyBody::try_new(
        "not spatially uniform", 2.0, 0.01, 100.0, 0.0, 1.0, 1.0, &curve,
    )
    .unwrap();
    with_cx(|cx| {
        for rejected in [&radiative, &high_biot] {
            assert!(
                solve_lumped_driven_exchange(
                    cx, rejected, BiotGate::corpus_default(),
                    driven_config(10_000.0, 300.0, 600.0, 100.0),
                )
                .is_err()
            );
        }
        for (body, cfg) in [
            (&isolated, driven_config(95_000.0, 300.0, 1_000.0, 1.0)),
            (&isolated, driven_config(0.0, 300.0, -1_000.0, 1.0)),
            (&good, driven_config(10_000.0, 300.0, 600.0, 1e6)),
            (&good, driven_config(95_000.0, 800.0, -1_000.0, 1e6)),
        ] {
            assert!(
                solve_lumped_driven_exchange(cx, body, BiotGate::corpus_default(), cfg).is_err()
            );
        }
    });

    // Preserve the native phase owner's strict absolute-h rounding contract.
    let datum = 2.0_f64.powi(52);
    let shifted = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xb3; 32]),
        curve.knots().iter().map(|knot| {
            let mut knot = *knot;
            knot.specific_enthalpy_j_kg += datum;
            knot
        }).collect(),
    )
    .unwrap();
    let shifted_body = body(&shifted, 100.0, 0.0);
    with_cx(|cx| {
        let cfg = LumpedDrivenExchangeConfig {
            energy_tolerance_j: 1.0,
            ..driven_config(datum + 10_000.0, 300.0, 600.0, 0.007)
        };
        let accepted = solve_lumped_driven_exchange(
            cx, &shifted_body, BiotGate::corpus_default(), cfg,
        )
        .unwrap();
        assert!(accepted.heat_residual_j.abs() > 0.1);
        assert!(accepted.heat_residual_j.abs() < 1.0);
        near(
            accepted.heat_residual_j,
            accepted.heat_into_body_j - accepted.stored_enthalpy_change_j,
            1e-14,
        );
        balanced(&accepted, cfg.energy_tolerance_j);
        for refused in [
            LumpedDrivenExchangeConfig { energy_tolerance_j: 1e-6, ..cfg },
            LumpedDrivenExchangeConfig { duration_s: 1e-8, ..cfg },
        ] {
            assert!(
                solve_lumped_driven_exchange(
                    cx, &shifted_body, BiotGate::corpus_default(), refused,
                )
                .is_err()
            );
            assert!(
                solve_lumped_driven_exchange_response(
                    cx, &shifted_body, BiotGate::corpus_default(), refused,
                )
                .is_err()
            );
        }
    });
    support::with_cancelled_cx(|cx| {
        let cfg = driven_config(10_000.0, 300.0, 600.0, 100.0);
        assert!(matches!(
            solve_lumped_driven_exchange(cx, &good, BiotGate::corpus_default(), cfg),
            Err(ConductionError::Cancelled { .. }),
        ));
        assert!(matches!(
            solve_lumped_driven_exchange_response(cx, &good, BiotGate::corpus_default(), cfg),
            Err(ConductionError::Cancelled { .. }),
        ));
    });
}

// Independent piecewise solution for the original test chart, G=1 W/K and m=2 kg.
fn driven_history_enthalpy(initial_h: f64, equilibrium: f64, elapsed: f64, heating: bool) -> f64 {
    let x = elapsed / 2.0;
    if heating {
        let initial_t = 300.0 + initial_h / 100.0;
        let solid = 100.0 * ((equilibrium - initial_t) / (equilibrium - 600.0)).ln();
        if x <= solid {
            let t = equilibrium - (equilibrium - initial_t) * (-x / 100.0).exp();
            return 100.0 * (t - 300.0);
        }
        let latent = 25_000.0 / (equilibrium - 600.0);
        if x <= solid + latent {
            return 30_000.0 + (x - solid) * (equilibrium - 600.0);
        }
        let t = equilibrium - (equilibrium - 600.0) * (-(x - solid - latent) / 200.0).exp();
        55_000.0 + 200.0 * (t - 600.0)
    } else {
        let initial_t = 600.0 + (initial_h - 55_000.0) / 200.0;
        let liquid = 200.0 * ((initial_t - equilibrium) / (600.0 - equilibrium)).ln();
        if x <= liquid {
            let t = equilibrium + (initial_t - equilibrium) * (-x / 200.0).exp();
            return 55_000.0 + 200.0 * (t - 600.0);
        }
        let latent = 25_000.0 / (600.0 - equilibrium);
        if x <= liquid + latent {
            return 55_000.0 - (x - liquid) * (600.0 - equilibrium);
        }
        let t = equilibrium + (600.0 - equilibrium) * (-(x - liquid - latent) / 100.0).exp();
        100.0 * (t - 300.0)
    }
}

#[test]
fn g1_exact_history_resolves_heating_and_freezing_independently_of_output_cadence() {
    let curve = phase_curve();
    let body = body(&curve, 100.0, 0.0);
    let heating =
        100.0 * (500.0_f64 / 300.0).ln() + 25_000.0 / 300.0 + 200.0 * (300.0_f64 / 200.0).ln();
    let cooling =
        200.0 * (1_000.0_f64 / 800.0).ln() + 25_000.0 / 800.0 + 100.0 * (800.0_f64 / 600.0).ln();
    with_cx(|cx| {
        for (initial, bath, power, exposure, target_h, target_t, heating) in [
            (10_000.0, 300.0, 600.0, heating, 75_000.0, 700.0, true),
            (95_000.0, 800.0, -1_000.0, cooling, 10_000.0, 400.0, false),
        ] {
            let duration = 2.0 * exposure;
            for cadence in [duration, 20.0] {
                let cfg = LumpedEnthalpyMarchConfig {
                    initial_specific_enthalpy_j_kg: initial,
                    ambient_temperature_k: bath,
                    radiation_temperature_k: 300.0,
                    internal_power_w: power,
                    duration_s: duration,
                    maximum_step_s: cadence,
                    maximum_steps: 32,
                    enthalpy_tolerance_j_kg: 5e-8,
                };
                let march =
                    solve_lumped_enthalpy_exact(cx, &body, BiotGate::corpus_default(), cfg).unwrap();
                let samples = march.samples();
                if cadence == duration {
                    assert_eq!(samples.len(), 2);
                } else {
                    assert!(samples.len() > 2);
                    assert!(samples.iter().any(|sample| {
                        let fraction = sample.phase_state.liquid_mass_fraction();
                        fraction > 0.0 && fraction < 1.0
                    }));
                }
                assert_eq!(samples[0].time_s, 0.0);
                let last = samples.last().unwrap();
                assert_eq!(last.time_s, duration);
                near(last.phase_state.specific_enthalpy_j_kg(), target_h, 2e-7);
                near(last.phase_state.temperature_k(), target_t, 2e-9);
                for sample in samples {
                    let expected_h = driven_history_enthalpy(initial, bath + power, sample.time_s, heating);
                    near(sample.phase_state.specific_enthalpy_j_kg(), expected_h, 2e-7);
                    // Reported powers are instantaneous at the actual native
                    // state, while interval energy uses the exact time mean.
                    near(
                        sample.convection_into_body_w,
                        bath - sample.phase_state.temperature_k(),
                        1e-12,
                    );
                    near(
                        sample.net_power_into_body_w,
                        power + bath - sample.phase_state.temperature_k(),
                        1e-12,
                    );
                    assert_eq!(sample.radiation_into_body_w, 0.0);
                    assert!(sample.step_energy_residual_j.abs() <= 2.0 * cfg.enthalpy_tolerance_j_kg);
                }
                let budget = (samples.len() - 1) as f64 * 2.0 * cfg.enthalpy_tolerance_j_kg;
                assert!(march.cumulative_absolute_energy_residual_j() <= budget);
            }
        }
    });
}

#[test]
fn g0_exact_history_retains_energy_radiation_step_budget_and_cancellation_admission() {
    let curve = phase_curve();
    let good = body(&curve, 100.0, 0.0);
    let isolated = body(&curve, 0.0, 0.0);
    let radiative = body(&curve, 100.0, 0.2);
    let cfg = LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 10_000.0,
        ambient_temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        internal_power_w: 600.0,
        duration_s: 100.0,
        maximum_step_s: 100.0,
        maximum_steps: 1,
        enthalpy_tolerance_j_kg: 5e-8,
    };
    with_cx(|cx| {
        assert!(
            solve_lumped_enthalpy_exact(cx, &radiative, BiotGate::corpus_default(), cfg).is_err()
        );
        assert!(
            solve_lumped_enthalpy_exact(
                cx, &good, BiotGate::corpus_default(),
                LumpedEnthalpyMarchConfig { maximum_step_s: 20.0, ..cfg },
            )
            .is_err()
        );
        for refused in [
            LumpedEnthalpyMarchConfig {
                initial_specific_enthalpy_j_kg: 95_000.0,
                internal_power_w: 1_000.0,
                duration_s: 1.0,
                ..cfg
            },
            LumpedEnthalpyMarchConfig {
                initial_specific_enthalpy_j_kg: 20_000.0,
                internal_power_w: 1e-12,
                duration_s: 1.0,
                enthalpy_tolerance_j_kg: 1e-16,
                ..cfg
            },
        ] {
            assert!(
                solve_lumped_enthalpy_exact(cx, &isolated, BiotGate::corpus_default(), refused)
                    .is_err()
            );
        }
    });
    support::with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_lumped_enthalpy_exact(cx, &good, BiotGate::corpus_default(), cfg),
            Err(ConductionError::Cancelled { .. }),
        ));
    });
}
