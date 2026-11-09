//! Native uniform-body exchange through sensible branches and sharp fusion.

#[path = "exchange/driven.rs"]
mod driven;
#[path = "exchange/pair.rs"]
mod pair;

use super::*;
use fs_conduction::lumped::exchange::{
    LumpedLinearExchangeConfig, solve_lumped_linear_exchange, solve_lumped_linear_exchange_response,
};

fn body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    htc: f64,
    emissivity: f64,
) -> LumpedEnthalpyBody<'_> {
    LumpedEnthalpyBody::try_new(
        "exact phase bath",
        2.0,
        0.01,
        htc,
        emissivity,
        0.001,
        35.0,
        curve,
    )
    .unwrap()
}

fn config(initial: f64, bath: f64, duration: f64) -> LumpedLinearExchangeConfig {
    LumpedLinearExchangeConfig {
        initial_specific_enthalpy_j_kg: initial,
        reservoir_temperature_k: bath,
        duration_s: duration,
        energy_tolerance_j: 1e-7,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_exact_bath_crosses_fusion_with_original_mass_heat_and_convection_mean() {
    let curve = phase_curve();
    let body = body(&curve, 100.0, 0.0); // G=1 W/K, mass=2 kg.
    let heating =
        100.0 * (500.0_f64 / 300.0).ln() + 25_000.0 / 300.0 + 200.0 * (300.0_f64 / 200.0).ln();
    let cooling =
        200.0 * (600.0_f64 / 400.0).ln() + 25_000.0 / 400.0 + 100.0 * (400.0_f64 / 200.0).ln();
    for (h0, bath, exposure, target_h, target_t) in [
        (10_000.0, 900.0, heating, 75_000.0, 700.0),
        (95_000.0, 200.0, cooling, 10_000.0, 400.0),
    ] {
        with_cx(|cx| {
            let cfg = config(h0, bath, 2.0 * exposure);
            let result =
                solve_lumped_linear_exchange(cx, &body, BiotGate::corpus_default(), cfg).unwrap();
            assert_eq!(result.body_identity, body.identity());
            near(result.caloric.state.temperature_k(), target_t, 2e-10);
            near(
                result.caloric.state.specific_enthalpy_j_kg(),
                target_h,
                2e-8,
            );
            near(result.heat_into_body_j, 2.0 * (target_h - h0), 5e-8);
            near(
                result.stored_enthalpy_change_j,
                result.heat_into_body_j,
                1e-7,
            );
            near(
                cfg.duration_s * (bath - result.caloric.mean_temperature_k),
                result.heat_into_body_j,
                1e-7,
            );
            assert!(result.heat_residual_j.abs() <= cfg.energy_tolerance_j);
            assert!(result.maximum_biot < LUMPED_BIOT_CEILING);
            // Splitting a held bath creates no time-discretization model.
            let first = solve_lumped_linear_exchange(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedLinearExchangeConfig {
                    duration_s: cfg.duration_s * 0.5,
                    ..cfg
                },
            )
            .unwrap();
            let second = solve_lumped_linear_exchange(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedLinearExchangeConfig {
                    initial_specific_enthalpy_j_kg: first.caloric.state.specific_enthalpy_j_kg(),
                    duration_s: cfg.duration_s * 0.5,
                    ..cfg
                },
            )
            .unwrap();
            near(
                second.caloric.state.temperature_k(),
                result.caloric.state.temperature_k(),
                2e-10,
            );
            near(
                first.heat_into_body_j + second.heat_into_body_j,
                result.heat_into_body_j,
                1e-7,
            );
        });
    }
}

#[test]
fn g1_latent_heat_and_zero_transfer_retain_original_phase_inventory() {
    let curve = phase_curve();
    let coupled = body(&curve, 100.0, 0.0);
    with_cx(|cx| {
        let result = solve_lumped_linear_exchange(
            cx,
            &coupled,
            BiotGate::corpus_default(),
            config(35_000.0, 800.0, 100.0),
        )
        .unwrap();
        // Q=G*dt*(800-600)=20000 J; m=2 kg, so h rises by 10000 J/kg.
        assert_eq!(result.caloric.state.temperature_k(), 600.0);
        assert_eq!(result.caloric.mean_temperature_k, 600.0);
        assert_eq!(result.caloric.state.specific_enthalpy_j_kg(), 45_000.0);
        near(result.caloric.state.liquid_mass_fraction(), 0.6, 2e-15);
        assert_eq!(result.heat_into_body_j, 20_000.0);
        assert_eq!(result.heat_residual_j, 0.0);
        let isolated = body(&curve, 0.0, 0.0);
        for (body, cfg) in [
            (&isolated, config(35_000.0, 800.0, 1e10)),
            (&coupled, config(35_000.0, 600.0, 1e10)),
            (&coupled, config(35_000.0, 800.0, 0.0)),
        ] {
            let stopped =
                solve_lumped_linear_exchange(cx, body, BiotGate::corpus_default(), cfg).unwrap();
            assert_eq!(stopped.caloric.state, stopped.caloric.initial);
            assert_eq!(stopped.heat_into_body_j, 0.0);
        }
    });
}

#[test]
fn g0_original_biot_radiation_domain_and_cancellation_admission_remain_mandatory() {
    let curve = phase_curve();
    let good = body(&curve, 100.0, 0.0);
    let radiative = body(&curve, 100.0, 0.2);
    let high_biot = LumpedEnthalpyBody::try_new(
        "not spatially uniform",
        2.0,
        0.01,
        100.0,
        0.0,
        1.0,
        1.0,
        &curve,
    )
    .unwrap();
    with_cx(|cx| {
        for rejected in [&radiative, &high_biot] {
            assert!(
                solve_lumped_linear_exchange(
                    cx,
                    rejected,
                    BiotGate::corpus_default(),
                    config(10_000.0, 900.0, 100.0)
                )
                .is_err()
            );
        }
        for cfg in [
            config(10_000.0, 900.0, -1.0),
            config(10_000.0, 0.0, 1.0),
            config(10_000.0, 900.0, 1e10),
            LumpedLinearExchangeConfig {
                energy_tolerance_j: f64::NAN,
                ..config(10_000.0, 900.0, 1.0)
            },
        ] {
            assert!(
                solve_lumped_linear_exchange(cx, &good, BiotGate::corpus_default(), cfg).is_err()
            );
        }
    });
    support::with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_lumped_linear_exchange(
                cx,
                &good,
                BiotGate::corpus_default(),
                config(10_000.0, 900.0, 100.0)
            ),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn g0_energy_budget_measures_actual_state_quantization_instead_of_absolute_h() {
    let original = phase_curve();
    let datum = 2.0_f64.powi(52);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xb2; 32]),
        original
            .knots()
            .iter()
            .map(|knot| {
                let mut knot = *knot;
                knot.specific_enthalpy_j_kg += datum;
                knot
            })
            .collect(),
    )
    .unwrap();
    let body = body(&curve, 100.0, 0.0);
    with_cx(|cx| {
        let cfg = LumpedLinearExchangeConfig {
            energy_tolerance_j: 1.0,
            ..config(datum + 10_000.0, 900.0, 0.007)
        };
        let result =
            solve_lumped_linear_exchange(cx, &body, BiotGate::corpus_default(), cfg).unwrap();
        let response =
            solve_lumped_linear_exchange_response(cx, &body, BiotGate::corpus_default(), cfg)
                .unwrap();
        assert_eq!(response.exchange, result);
        assert!(
            solve_lumped_linear_exchange_response(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedLinearExchangeConfig {
                    energy_tolerance_j: 1e-6,
                    ..cfg
                },
            )
            .is_err()
        );
        assert!(result.heat_residual_j.abs() > 0.1);
        assert!(result.heat_residual_j.abs() < 1.0);
        near(
            result.heat_residual_j,
            result.heat_into_body_j - result.stored_enthalpy_change_j,
            1e-14,
        );
        assert!(
            solve_lumped_linear_exchange(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedLinearExchangeConfig {
                    energy_tolerance_j: 1e-6,
                    ..cfg
                }
            )
            .is_err()
        );
        assert!(
            solve_lumped_linear_exchange(
                cx,
                &body,
                BiotGate::corpus_default(),
                LumpedLinearExchangeConfig {
                    duration_s: 1e-8,
                    ..cfg
                }
            )
            .is_err()
        );
    });
}

// Independent closed-form endpoint after all three phases of the manufactured
// constant-Cp fixture. Its perturbed endpoints stay inside the same final phase.
fn three_phase_heat(initial_h: f64, bath: f64, duration: f64, heating: bool) -> f64 {
    let x = duration / 2.0; // G=1 W/K and m=2 kg.
    let final_h = if heating {
        let initial_t = 300.0 + initial_h / 100.0;
        let remaining =
            x - 100.0 * ((bath - initial_t) / (bath - 600.0)).ln() - 25_000.0 / (bath - 600.0);
        let final_t = bath - (bath - 600.0) * (-remaining / 200.0).exp();
        55_000.0 + 200.0 * (final_t - 600.0)
    } else {
        let initial_t = 600.0 + (initial_h - 55_000.0) / 200.0;
        let remaining =
            x - 200.0 * ((initial_t - bath) / (600.0 - bath)).ln() - 25_000.0 / (600.0 - bath);
        let final_t = bath + (600.0 - bath) * (-remaining / 100.0).exp();
        100.0 * (final_t - 300.0)
    };
    2.0 * (final_h - initial_h)
}

#[test]
fn g1_native_heat_jacobian_matches_independent_cross_phase_derivatives_and_units() {
    let curve = phase_curve();
    let body = body(&curve, 100.0, 0.0);
    let heating =
        100.0 * (500.0_f64 / 300.0).ln() + 25_000.0 / 300.0 + 200.0 * (300.0_f64 / 200.0).ln();
    let cooling =
        200.0 * (600.0_f64 / 400.0).ln() + 25_000.0 / 400.0 + 100.0 * (400.0_f64 / 200.0).ln();
    with_cx(|cx| {
        for (initial, bath, duration, is_heating) in [
            (10_000.0, 900.0, 2.0 * heating, true),
            (95_000.0, 200.0, 2.0 * cooling, false),
        ] {
            let response = solve_lumped_linear_exchange_response(
                cx,
                &body,
                BiotGate::corpus_default(),
                config(initial, bath, duration),
            )
            .unwrap();
            near(
                response.exchange.heat_into_body_j,
                three_phase_heat(initial, bath, duration, is_heating),
                1e-7,
            );
            let derivative_h = (three_phase_heat(initial + 0.01, bath, duration, is_heating)
                - three_phase_heat(initial - 0.01, bath, duration, is_heating))
                / 0.02;
            let derivative_bath = (three_phase_heat(initial, bath + 0.001, duration, is_heating)
                - three_phase_heat(initial, bath - 0.001, duration, is_heating))
                / 0.002;
            let derivative_duration =
                (three_phase_heat(initial, bath, duration + 0.001, is_heating)
                    - three_phase_heat(initial, bath, duration - 0.001, is_heating))
                    / 0.002;
            near(
                response.heat_initial_enthalpy_derivative_kg,
                derivative_h,
                1e-6,
            );
            near(
                response.heat_reservoir_derivative_j_k,
                derivative_bath,
                1e-4,
            );
            near(
                response.heat_duration_derivative_w,
                derivative_duration,
                1e-5,
            );
            assert_eq!(
                response.heat_initial_enthalpy_derivative_kg,
                -2.0 * response.caloric_jacobian.inlet_enthalpy_loss
            );
            assert_eq!(
                response.heat_reservoir_derivative_j_k,
                2.0 * response
                    .caloric_jacobian
                    .reservoir_enthalpy_derivative_j_kg_k
            );
            assert_eq!(
                response.heat_duration_derivative_w,
                response.caloric_jacobian.exposure_enthalpy_derivative_k
            ); // G=1.
            assert_eq!(
                response.exchange,
                solve_lumped_linear_exchange(
                    cx,
                    &body,
                    BiotGate::corpus_default(),
                    config(initial, bath, duration),
                )
                .unwrap()
            );
        }
    });
}

#[test]
fn g0_native_response_retains_latent_inventory_energy_budget_and_cancellation() {
    let curve = phase_curve();
    let coupled = body(&curve, 100.0, 0.0);
    with_cx(|cx| {
        let plateau = solve_lumped_linear_exchange_response(
            cx,
            &coupled,
            BiotGate::corpus_default(),
            config(35_000.0, 600.0, 100.0),
        )
        .unwrap();
        assert_eq!(plateau.exchange.heat_into_body_j, 0.0);
        assert_eq!(plateau.heat_initial_enthalpy_derivative_kg, 0.0);
        assert_eq!(plateau.heat_reservoir_derivative_j_k, 100.0); // G*dt while latent.
        assert_eq!(plateau.heat_duration_derivative_w, 0.0);
        assert!(
            solve_lumped_linear_exchange_response(
                cx,
                &coupled,
                BiotGate::corpus_default(),
                config(30_000.0, 600.0, 100.0),
            )
            .is_err()
        );
        let zero = solve_lumped_linear_exchange_response(
            cx,
            &coupled,
            BiotGate::corpus_default(),
            config(10_000.0, 900.0, 0.0),
        )
        .unwrap();
        assert_eq!(zero.heat_reservoir_derivative_j_k, 0.0);
        assert_eq!(zero.heat_duration_derivative_w, 500.0); // Initial physical G*(Tb-T).
    });
    support::with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_lumped_linear_exchange_response(
                cx,
                &coupled,
                BiotGate::corpus_default(),
                config(10_000.0, 900.0, 100.0),
            ),
            Err(ConductionError::Cancelled { .. }),
        ));
    });
}
