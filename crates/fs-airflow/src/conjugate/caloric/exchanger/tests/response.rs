//! Analytic-response checks against independent steady balances and perturbations.

use super::*;

fn response_config() -> ExchangerResponseConfig {
    ExchangerResponseConfig::default()
}

#[test]
fn constant_cp_response_matches_independent_ntu_and_reciprocal_outlet_balances() {
    with_cx(|cx| {
        for arrangement in [FlowArrangement::CoCurrent, FlowArrangement::Counterflow] {
            let exchanger = CaloricHeatExchanger::new(
                stream(600.0, 1.0, affine(1000.0, 0.0, 400.0, 100_000.0)),
                stream(300.0, 1.5, affine(2000.0, 0.0, 650.0, 250_000.0)),
                arrangement,
                750.0,
            )
            .unwrap();
            let (state, response) = exchanger
                .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                .unwrap();
            assert_eq!(
                state,
                exchanger.solve(cx, ExchangerConfig::default()).unwrap()
            );
            let ntu: f64 = 0.75;
            let ratio: f64 = 1.0 / 3.0;
            let effectiveness = match arrangement {
                FlowArrangement::CoCurrent => -(-ntu * (1.0 + ratio)).exp_m1() / (1.0 + ratio),
                FlowArrangement::Counterflow => {
                    let decay = (-ntu * (1.0 - ratio)).exp();
                    (1.0 - decay) / (1.0 - ratio * decay)
                }
            };
            let heat = effectiveness * 300_000.0;
            let expected = [heat / (1000.0 * 300.0), -heat / (2000.0 * 300.0)];
            close(
                response.heat_from_first_inlet_enthalpy_kg_per_s,
                expected[0],
                2e-14,
            );
            close(
                response.heat_from_second_inlet_enthalpy_kg_per_s,
                expected[1],
                2e-14,
            );
            close(
                response.first_outlet_enthalpy_gain,
                1.0 - expected[0],
                2e-14,
            );
            close(
                response.second_outlet_enthalpy_gain,
                1.0 + expected[1] / 1.5,
                2e-14,
            );
            let end_gap = state.first_outlet_temperature_k
                - if arrangement == FlowArrangement::CoCurrent {
                    state.second_outlet_temperature_k
                } else {
                    300.0
                };
            let conductance = if arrangement == FlowArrangement::CoCurrent {
                end_gap
            } else {
                (600.0 - state.second_outlet_temperature_k) * end_gap / 300.0
            };
            close(response.heat_from_conductance_k, conductance, 3e-11);
            // One differentiated Q still credits and debits the two mass flows.
            let direction = [37.0, -19.0];
            let dq = response.heat_from_first_inlet_enthalpy_kg_per_s * direction[0]
                + response.heat_from_second_inlet_enthalpy_kg_per_s * direction[1];
            let first_out = direction[0] - dq;
            let second_out = direction[1] + dq / 1.5;
            close(
                (first_out - direction[0]) + 1.5 * (second_out - direction[1]),
                0.0,
                2e-14,
            );
        }
    });
}

#[test]
fn original_variable_cp_response_matches_independent_implicit_ua_moments() {
    with_cx(|cx| {
        // Independently integrate the manufactured quadratic inverse with
        // fixed Simpson panels; these are not values from the response owner.
        let cases = [
            (
                FlowArrangement::CoCurrent,
                919.5863752576657,
                0.31216291445439015,
                -0.3584898604415899,
                96.52480449312009,
            ),
            (
                FlowArrangement::Counterflow,
                816.2887925555515,
                0.306213977060121,
                -0.35602595128030173,
                140.87438559528732,
            ),
        ];
        for (arrangement, ua, first, second, conductance) in cases {
            for (first_datum, second_datum) in [(250.0, 700.0), (700.0, 250.0)] {
                let exchanger = CaloricHeatExchanger::new(
                    stream(620.0, 0.8, affine(600.0, 2.0, first_datum, 100_000.0)),
                    stream(280.0, 1.1, affine(1800.0, -1.0, second_datum, 300_000.0)),
                    arrangement,
                    ua,
                )
                .unwrap();
                let (state, response) = exchanger
                    .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                    .unwrap();
                close(state.heat_rate_first_to_second_w, 180_000.0, 4e-5);
                close(
                    response.heat_from_first_inlet_enthalpy_kg_per_s,
                    first,
                    5e-9,
                );
                close(
                    response.heat_from_second_inlet_enthalpy_kg_per_s,
                    second,
                    5e-9,
                );
                close(response.heat_from_conductance_k, conductance, 2e-6);
                close(response.first_outlet_enthalpy_gain, 1.0 - first / 0.8, 7e-9);
                close(
                    response.second_outlet_enthalpy_gain,
                    1.0 + second / 1.1,
                    7e-9,
                );
                audit(state, &exchanger);
                let reversed = CaloricHeatExchanger::new(
                    exchanger.second.clone(),
                    exchanger.first.clone(),
                    arrangement,
                    ua,
                )
                .unwrap();
                let (_, reverse) = reversed
                    .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                    .unwrap();
                assert_eq!(
                    reverse.heat_from_first_inlet_enthalpy_kg_per_s,
                    -response.heat_from_second_inlet_enthalpy_kg_per_s
                );
                assert_eq!(
                    reverse.heat_from_second_inlet_enthalpy_kg_per_s,
                    -response.heat_from_first_inlet_enthalpy_kg_per_s
                );
                assert_eq!(
                    reverse.heat_from_conductance_k,
                    -response.heat_from_conductance_k
                );
                assert_eq!(
                    reverse.first_outlet_enthalpy_gain,
                    response.second_outlet_enthalpy_gain
                );
            }
        }
    });
}

#[test]
fn equilibrium_uses_actual_local_cp_and_retains_extremely_weak_counterflow_anchor() {
    with_cx(|cx| {
        let law = affine(600.0, 2.0, 250.0, 100_000.0);
        let cp = 600.0 + 2.0 * 420.0;
        for difference in [0.0, 1e-12] {
            let exchanger = CaloricHeatExchanger::new(
                stream(420.0, 1.0, law.clone()),
                stream(420.0, 1.0 + difference, law.clone()),
                FlowArrangement::Counterflow,
                1700.0,
            )
            .unwrap();
            let (state, response) = exchanger
                .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                .unwrap();
            assert_eq!(state.heat_rate_first_to_second_w, 0.0);
            assert_eq!(response.heat_from_conductance_k, 0.0);
            close(
                response.heat_from_first_inlet_enthalpy_kg_per_s,
                1700.0 / (cp + 1700.0),
                1e-12,
            );
            close(
                response.first_outlet_enthalpy_gain,
                cp / (cp + 1700.0),
                1e-12,
            );
        }
        let exchanger = CaloricHeatExchanger::new(
            stream(420.0, 1.0, law.clone()),
            stream(420.0, 1.0, law),
            FlowArrangement::Counterflow,
            1e22,
        )
        .unwrap();
        let (_, response) = exchanger
            .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
            .unwrap();
        assert_eq!(response.heat_from_first_inlet_enthalpy_kg_per_s, 1.0);
        assert!(response.first_outlet_enthalpy_gain > 0.0);
        close(
            response.first_outlet_enthalpy_gain,
            1.0 / (1.0 + 1e22 / cp),
            1e-33,
        );
        assert_eq!(
            response.first_outlet_enthalpy_gain,
            response.second_outlet_enthalpy_gain
        );
    });
}

#[test]
fn weak_variable_source_exchange_retains_nonzero_heat_sensitivity() {
    with_cx(|cx| {
        for arrangement in [FlowArrangement::CoCurrent, FlowArrangement::Counterflow] {
            let exchanger = CaloricHeatExchanger::new(
                stream(500.0, 1.0, affine(600.0, 2.0, 700.0, 100_000.0)),
                stream(300.0, 1.0, affine(1800.0, -1.0, 250.0, 200_000.0)),
                arrangement,
                1e-12,
            )
            .unwrap();
            let (_, response) = exchanger
                .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                .unwrap();
            assert!(response.heat_from_first_inlet_enthalpy_kg_per_s > 0.0);
            close(
                response.heat_from_first_inlet_enthalpy_kg_per_s,
                1e-12 / 1600.0,
                1e-27,
            );
            close(
                response.heat_from_second_inlet_enthalpy_kg_per_s,
                -1e-12 / 1500.0,
                1e-27,
            );
            close(response.heat_from_conductance_k, 200.0, 2e-10);
        }
    });
}

#[test]
fn response_crosses_original_cp_kinks_and_matches_independent_perturbed_solves() {
    with_cx(|cx| {
        let first = stream(
            650.0,
            1.0,
            source(
                &[(250.0, 1000.0), (450.0, 1000.0), (700.0, 2000.0)],
                450.0,
                100_000.0,
            ),
        );
        let second = stream(
            280.0,
            1.3,
            source(
                &[(250.0, 1300.0), (350.0, 1500.0), (700.0, 1150.0)],
                350.0,
                200_000.0,
            ),
        );
        let ua = 1826.5130379126617;
        let exchanger = CaloricHeatExchanger::new(
            first.clone(),
            second.clone(),
            FlowArrangement::Counterflow,
            ua,
        )
        .unwrap();
        let (_, response) = exchanger
            .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
            .unwrap();
        let evaluate = |first_h: f64, second_h: f64, conductance: f64| {
            let first_t = first
                .caloric
                .temperature_after_specific_heat(650.0, first_h)
                .unwrap();
            let second_t = second
                .caloric
                .temperature_after_specific_heat(280.0, second_h)
                .unwrap();
            CaloricHeatExchanger::new(
                stream(first_t, 1.0, first.caloric.clone()),
                stream(second_t, 1.3, second.caloric.clone()),
                FlowArrangement::Counterflow,
                conductance,
            )
            .unwrap()
            .solve(cx, ExchangerConfig::default())
            .unwrap()
            .heat_rate_first_to_second_w
        };
        let dh = 200.0;
        close(
            response.heat_from_first_inlet_enthalpy_kg_per_s,
            (evaluate(dh, 0.0, ua) - evaluate(-dh, 0.0, ua)) / (2.0 * dh),
            3e-7,
        );
        close(
            response.heat_from_second_inlet_enthalpy_kg_per_s,
            (evaluate(0.0, dh, ua) - evaluate(0.0, -dh, ua)) / (2.0 * dh),
            3e-7,
        );
        close(
            response.heat_from_conductance_k,
            evaluate(0.0, 0.0, ua + 0.5) - evaluate(0.0, 0.0, ua - 0.5),
            2e-4,
        );
    });
}

#[test]
fn a_response_that_misses_unity_by_three_summed_allowances_is_still_refused() {
    with_cx(|cx| {
        for arrangement in [FlowArrangement::CoCurrent, FlowArrangement::Counterflow] {
            let exchanger = CaloricHeatExchanger::new(
                stream(620.0, 0.8, affine(600.0, 2.0, 250.0, 100_000.0)),
                stream(280.0, 1.1, affine(1800.0, -1.0, 700.0, 300_000.0)),
                arrangement,
                900.0,
            )
            .unwrap();
            let (_, response) = exchanger
                .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                .unwrap();
            // The response the solver returned passes its own gate.
            response.check(&exchanger, response_config()).unwrap();
            // loss + gain - 1 holds to the summed allowance (loss error + gain error); moving the
            // gain by three of them leaves the identity at least twice that allowance out.
            let mut first = response;
            first.first_outlet_enthalpy_gain += 3.0 * response.first_enthalpy_gain_error_allowance;
            assert!(matches!(
                first.check(&exchanger, response_config()),
                Err(ExchangerError::Unresolved("exchanger enthalpy-response gate"))
            ));
            let mut second = response;
            second.second_outlet_enthalpy_gain +=
                3.0 * response.second_enthalpy_gain_error_allowance;
            assert!(matches!(
                second.check(&exchanger, response_config()),
                Err(ExchangerError::Unresolved("exchanger enthalpy-response gate"))
            ));
        }
    });
}

#[test]
fn zero_ua_and_refusals_preserve_original_admission_and_cancellation() {
    let first = stream(620.0, 0.8, affine(600.0, 2.0, 250.0, 100_000.0));
    let second = stream(280.0, 1.1, affine(1800.0, -1.0, 700.0, 300_000.0));
    let zero = CaloricHeatExchanger::new(
        first.clone(),
        second.clone(),
        FlowArrangement::CoCurrent,
        0.0,
    )
    .unwrap();
    with_cx(|cx| {
        let (_, response) = zero
            .solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
            .unwrap();
        assert_eq!(response.heat_from_first_inlet_enthalpy_kg_per_s, 0.0);
        assert_eq!(response.first_outlet_enthalpy_gain, 1.0);
        assert_eq!(response.heat_from_conductance_k, 340.0);
        assert!(
            zero.solve_with_inlet_response(
                cx,
                ExchangerConfig::default(),
                ExchangerResponseConfig {
                    absolute_enthalpy_gain_tolerance: 0.0,
                    ..response_config()
                }
            )
            .is_err()
        );
        let exchanger = CaloricHeatExchanger::new(
            first.clone(),
            second.clone(),
            FlowArrangement::CoCurrent,
            350.0,
        )
        .unwrap();
        assert!(exchanger.solve(cx, ExchangerConfig::default()).is_ok());
        let refusal = exchanger.solve_with_inlet_response(
            cx,
            ExchangerConfig::default(),
            ExchangerResponseConfig {
                absolute_enthalpy_gain_tolerance: 1e-300,
                relative_enthalpy_gain_tolerance: 0.0,
            },
        );
        assert!(matches!(
            refusal,
            Err(ExchangerError::WorkLimit | ExchangerError::Unresolved(_))
        ));
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert_eq!(
            zero.solve_with_inlet_response(cx, ExchangerConfig::default(), response_config())
                .unwrap_err(),
            ExchangerError::Interrupted
        );
    });
}
