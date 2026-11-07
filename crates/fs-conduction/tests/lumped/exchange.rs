//! Native uniform-body exchange through sensible branches and sharp fusion.

use super::*;
use fs_conduction::lumped::exchange::{LumpedLinearExchangeConfig, solve_lumped_linear_exchange};

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
