use super::*;
use fs_qty::Area;

const SIGMA: f64 = 5.670_374_419e-8;

fn gray_body<'a>(chart: &'a EquilibriumEnthalpyPhaseCurve, h: f64, g: f64) -> FiniteHeatBody<'a> {
    body(chart, h, g)
        .with_gray_radiation(Area::new(0.2), 0.8)
        .unwrap()
}

#[test]
fn finite_gray_source_cp_stage_closes_heat_and_exact_eliminated_tangent() {
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    for cp in [[10.0, 110.0], [110.0, 10.0]] {
        for datum in [0.0, 1_048_576.0, -1_048_576.0] {
            let chart = source_cp_curve(datum, cp);
            for (initial, gas_temperature) in [(350.0, 800.0), (900.0, 300.0)] {
                let h0 = source_cp_enthalpy(cp, initial);
                let original = gray_body(&chart, datum + h0, 40.0);
                let (volume, diagonal_step) = (2.0, 0.5);
                let mut fixed = Sum::default();
                fixed.add(250.0);
                let stage = HeatStage::implicit(
                    fixed,
                    Some(HeatExchange::Body(original)),
                    volume,
                    diagonal_step,
                    20.0,
                )
                .unwrap();
                let a = diagonal_step * 40.0 / volume;
                let b = diagonal_step * SIGMA * 0.16 / volume;
                // Independent original-Cp polynomial and fourth-power equation,
                // with no native stage, caloric inversion or root helper.
                let (mut lo, mut hi) = (200.0, 1200.0);
                for _ in 0..80 {
                    let t = 0.5 * (lo + hi);
                    let heat = 270.0
                        + a * (t - gas_temperature)
                        + b * (t.powi(4) - gas_temperature.powi(4));
                    let residual = source_cp_enthalpy(cp, t) - h0 + volume * heat;
                    if residual < 0.0 {
                        lo = t;
                    } else {
                        hi = t;
                    }
                }
                let expected = 0.5 * (lo + hi);
                let (heat, feedback) = stage
                    .heat(gas_temperature, controls, &mut || false)
                    .unwrap();
                let actual_h = original.state().specific_enthalpy_j_kg() - heat.value() * volume;
                let actual = chart.temperature_at_specific_enthalpy(actual_h).unwrap();
                assert!((actual - expected).abs() < 3e-9, "{actual} != {expected}");
                let flux = 270.0
                    + a * (actual - gas_temperature)
                    + b * (actual.powi(4) - gas_temperature.powi(4));
                assert!((heat.value() - flux).abs() < 2e-8);
                let expected_feedback = (a + 4.0 * b * gas_temperature.powi(3))
                    / (1.0 + volume * (a + 4.0 * b * actual.powi(3)) / source_cp_at(cp, actual));
                assert!((feedback - expected_feedback).abs() < 1e-9);
                let epsilon = 0.001;
                let plus = stage
                    .heat(gas_temperature + epsilon, controls, &mut || false)
                    .unwrap()
                    .0
                    .value();
                let minus = stage
                    .heat(gas_temperature - epsilon, controls, &mut || false)
                    .unwrap()
                    .0
                    .value();
                assert!(
                    (feedback + (plus - minus) / (2.0 * epsilon)).abs() < 3e-7 * feedback.max(1.0)
                );
                assert_eq!(original.state().specific_enthalpy_j_kg(), datum + h0);
            }
        }
    }
}

fn rk4(initial: f64, duration: f64, steps: usize, derivative: impl Fn(f64, f64) -> f64) -> f64 {
    let mut value = initial;
    let dt = duration / steps as f64;
    for index in 0..steps {
        let time = index as f64 * dt;
        let a = derivative(time, value);
        let b = derivative(time + 0.5 * dt, value + 0.5 * dt * a);
        let c = derivative(time + 0.5 * dt, value + 0.5 * dt * b);
        let d = derivative(time + dt, value + dt * c);
        value += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    value
}

#[test]
fn finite_gray_body_changes_temperature_during_reaction_and_frozen_exchange() {
    for pressure in [false, true] {
        let gas_capacity = if pressure { 4.0 } else { 3.0 } * R;
        for (cp, gas_initial, body_initial, binding, reacting) in [
            ([10.0, 110.0], 800.0, 300.0, 512.0, true),
            ([110.0, 10.0], 300.0, 900.0, 0.0, false),
        ] {
            let h0 = source_cp_enthalpy(cp, body_initial);
            let gas_temperature = |time: f64, body_t: f64| {
                gas_initial
                    + (binding * R * (1.0 - (-2.0 * time).exp())
                        - (source_cp_enthalpy(cp, body_t) - h0))
                        / gas_capacity
            };
            let derivative = |time: f64, body_t: f64| {
                let gas_t = gas_temperature(time, body_t);
                (40.0 * (gas_t - body_t) + SIGMA * 0.16 * (gas_t.powi(4) - body_t.powi(4)))
                    / source_cp_at(cp, body_t)
            };
            let expected_body = rk4(body_initial, 0.3, 4096, derivative);
            assert!((expected_body - rk4(body_initial, 0.3, 8192, derivative)).abs() < 1e-9);
            let expected_gas = gas_temperature(0.3, expected_body);
            let chart = source_cp_curve(0.0, cp);
            let original = gray_body(&chart, h0, 40.0);
            for method in [
                IsochoricKineticsMethod::BackwardEuler,
                IsochoricKineticsMethod::Sdirk2,
            ] {
                let result = run(
                    pressure,
                    method,
                    gas_initial,
                    0.3,
                    original,
                    binding,
                    reacting,
                );
                let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                    0.2
                } else {
                    0.003
                };
                assert!(
                    (result.gas_t - expected_gas).abs() < tolerance,
                    "{method:?}, p={pressure}: gas {} != {expected_gas}",
                    result.gas_t
                );
                assert!((result.body_t - expected_body).abs() < tolerance);
                assert!((result.body_t - body_initial).abs() > 50.0);
                if reacting {
                    assert!((result.a - (-0.6_f64).exp()).abs() < 3e-4);
                }
                balances(result, 3.0 * R * gas_initial);
            }
        }
    }
}

#[test]
fn finite_gray_body_crosses_sensible_and_latent_storage_without_extra_heat() {
    for pressure in [false, true] {
        let gas_capacity = if pressure { 4.0 } else { 3.0 } * R;
        for (h0, gas_initial, binding, reacting) in [
            (19_900.0, 1000.0, 512.0, true),
            (40_010.0, 300.0, 0.0, false),
        ] {
            let gas_temperature = |time: f64, h: f64| {
                gas_initial + (binding * R * (1.0 - (-2.0 * time).exp()) - (h - h0)) / gas_capacity
            };
            // The manufactured original chart has Cp=100 and a 20 kJ/kg
            // plateau at 400 K. This independent enthalpy ODE includes the
            // initial sensible heating/cooling and ensuing phase inventory.
            let derivative = |time: f64, h: f64| {
                let body_t = if h < 20_000.0 {
                    200.0 + h / 100.0
                } else if h <= 40_000.0 {
                    400.0
                } else {
                    400.0 + (h - 40_000.0) / 100.0
                };
                let gas_t = gas_temperature(time, h);
                SIGMA * 0.16 * (gas_t.powi(4) - body_t.powi(4))
            };
            let expected_h = rk4(h0, 0.3, 8192, derivative);
            assert!((expected_h - rk4(h0, 0.3, 16384, derivative)).abs() < 1e-4);
            let expected_gas = gas_temperature(0.3, expected_h);
            for method in [
                IsochoricKineticsMethod::BackwardEuler,
                IsochoricKineticsMethod::Sdirk2,
            ] {
                let chart = curve(0.0, true);
                let result = run(
                    pressure,
                    method,
                    gas_initial,
                    0.3,
                    gray_body(&chart, h0, 0.0),
                    binding,
                    reacting,
                );
                let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                    0.15
                } else {
                    0.004
                };
                assert!(
                    (result.gas_t - expected_gas).abs() < tolerance,
                    "{method:?}, p={pressure}: {} != {expected_gas}",
                    result.gas_t
                );
                assert!((result.body_h - expected_h).abs() < gas_capacity * tolerance + 0.5);
                assert_eq!(result.body_t, 400.0);
                assert!(result.body_f > 0.0 && result.body_f < 1.0);
                assert!((result.body_f - (result.body_h - 20_000.0) / 20_000.0).abs() < 1e-14);
                assert!(result.body_work.abs() > 0.0);
                balances(result, 3.0 * R * gas_initial);
            }
        }
    }
}

#[test]
fn finite_gray_source_brackets_zero_mode_and_cancelled_roots_remain_explicit() {
    let chart = curve(0.0, false);
    let controls = thermal_controls(IsochoricKineticsMethod::Sdirk2);
    let original = body(&chart, 0.0, 0.0);
    let radiating = original.with_gray_radiation(Area::new(1.0), 1.0).unwrap();
    assert_eq!(
        radiating.with_gray_radiation(Area::new(1.0), 0.0).unwrap(),
        original
    );
    assert_eq!(
        radiating
            .with_conductance(ThermalConductance::new(2.0))
            .unwrap()
            .effective_radiating_area(),
        Area::new(1.0)
    );
    let stage = HeatStage::implicit(
        Sum::default(),
        Some(HeatExchange::Body(radiating)),
        1.0,
        0.5,
        0.0,
    )
    .unwrap();
    let (lower, upper) = stage
        .temperature_interval(200.0, 2000.0, controls, &mut || false)
        .unwrap();
    let expected_upper = (1200.0_f64.powi(4) + 20_000.0 / (0.5 * SIGMA))
        .sqrt()
        .sqrt();
    assert_eq!(lower, 200.0);
    assert!((upper - expected_upper).abs() < 1e-9);
    assert!(
        stage
            .heat(lower, controls, &mut || false)
            .unwrap()
            .0
            .value()
            .abs()
            < 1e-8
    );
    assert!(stage.heat(upper, controls, &mut || false).is_ok());
    assert!(stage.heat(upper + 1e-5, controls, &mut || false).is_err());
    assert!(matches!(
        stage.heat(700.0, controls, &mut || true),
        Err(IsochoricKineticsError::Cancelled)
    ));
    let limited = IsochoricKineticsSolve {
        maximum_temperature_iterations: 1,
        ..controls
    };
    assert!(matches!(
        stage.heat(1100.0, limited, &mut || false),
        Err(IsochoricKineticsError::TemperatureNoConvergence { .. })
    ));
    assert_eq!(radiating.state(), original.state());
}
