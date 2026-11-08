use super::*;
use crate::chemistry::heat::{HeatExchange, HeatStage};
use fs_qty::Area;

const SIGMA: f64 = 5.670_374_419e-8;

// Exact separated integral of dT/dt = k (Tb^4-T^4)/C. Its inverse uses only
// this analytic expression, independently of the native caloric/time solvers.
fn frozen_endpoint(initial: f64, reservoir: f64, capacity: f64, area: f64, time: f64) -> f64 {
    let primitive = |x: f64| 0.25 * ((1.0 + x) / (1.0 - x)).abs().ln() + 0.5 * x.atan();
    let target =
        primitive(initial / reservoir) + SIGMA * area * reservoir.powi(3) * time / capacity;
    let mut lo = initial.min(reservoir);
    let mut hi = initial.max(reservoir);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if (primitive(mid / reservoir) < target) == (initial < reservoir) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn gray_bath_heats_and_cools_against_exact_fourth_power_integral() {
    for pressure in [false, true] {
        let capacity = if pressure { 4.0 } else { 3.0 } * R;
        for method in [
            IsochoricKineticsMethod::BackwardEuler,
            IsochoricKineticsMethod::Sdirk2,
        ] {
            for (initial, reservoir) in [(1400.0, 350.0), (350.0, 900.0)] {
                let radiation = bath(0.0, reservoir)
                    .with_gray_radiation(Area::new(0.2), 0.8)
                    .unwrap();
                let result = advance(
                    thermal_mechanism(0.0, 2.0, 0.0),
                    pressure,
                    [Amount::new(0.0), Amount::new(1.0)],
                    initial,
                    0.4,
                    radiation,
                    thermal_controls(method),
                );
                let expected = frozen_endpoint(initial, reservoir, capacity, 0.16, 0.4);
                let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                    0.2
                } else {
                    0.005
                };
                assert!(
                    (result.temperature - expected).abs() < tolerance,
                    "{method:?}, p={pressure}: {} != {expected}",
                    result.temperature
                );
                assert_eq!(result.a, 0.0);
                assert!(
                    (result.balance.heat_input_j - capacity * (result.temperature - initial)).abs()
                        < 3e-8
                );
                assert!(
                    (result.balance.entropy_change_j_k
                        - capacity * (result.temperature / initial).ln())
                    .abs()
                        < 2e-11
                );
                check_balance(&result, 3.0 * R * initial, reservoir);
            }
        }
    }
}

// Independent scalar ODE: the first-order species population is known
// analytically. The fourth-power loss and chemical formation-energy release
// therefore share one explicit temperature equation, with no native stage,
// reaction solver, source inversion or heat helper used by this oracle.
fn reacting_endpoint(capacity: f64, steps: usize) -> f64 {
    let derivative = |time: f64, temperature: f64| {
        (512.0 * R * 2.0 * (-2.0 * time).exp()
            + 2.0 * (400.0 - temperature)
            + SIGMA * 0.12 * (400.0_f64.powi(4) - temperature.powi(4)))
            / capacity
    };
    let mut temperature = 1000.0;
    let dt = 0.3 / steps as f64;
    for index in 0..steps {
        let time = index as f64 * dt;
        let a = derivative(time, temperature);
        let b = derivative(time + 0.5 * dt, temperature + 0.5 * dt * a);
        let c = derivative(time + 0.5 * dt, temperature + 0.5 * dt * b);
        let d = derivative(time + dt, temperature + dt * c);
        temperature += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    temperature
}

#[test]
fn gray_reaction_couples_formation_energy_heat_and_pressure_work_once() {
    let remaining = (-0.6_f64).exp();
    for pressure in [false, true] {
        let capacity = if pressure { 4.0 } else { 3.0 } * R;
        let expected = reacting_endpoint(capacity, 2048);
        assert!((expected - reacting_endpoint(capacity, 4096)).abs() < 1e-9);
        for method in [
            IsochoricKineticsMethod::BackwardEuler,
            IsochoricKineticsMethod::Sdirk2,
        ] {
            let result = advance(
                thermal_mechanism(512.0, 2.0, 0.0),
                pressure,
                [Amount::new(1.0), Amount::new(0.0)],
                1000.0,
                0.3,
                bath(2.0, 400.0)
                    .with_gray_radiation(Area::new(0.15), 0.8)
                    .unwrap(),
                thermal_controls(method),
            );
            let tolerance = if method == IsochoricKineticsMethod::BackwardEuler {
                0.12
            } else {
                0.003
            };
            assert!(
                (result.temperature - expected).abs() < tolerance,
                "{method:?}, p={pressure}: {} != {expected}",
                result.temperature
            );
            assert!((result.a - remaining).abs() < 3e-4);
            let heat = capacity * (result.temperature - 1000.0) - 512.0 * R * (1.0 - result.a);
            assert!((result.balance.heat_input_j - heat).abs() < 3e-8);
            assert!(result.balance.heat_input_j < 0.0);
            check_balance(&result, 3.0 * R * 1000.0, 400.0);
        }
    }
}

#[test]
fn gray_stage_retains_kelvin_separation_original_linear_arithmetic_and_admission() {
    let reservoir = 1000.0_f64;
    let temperature = f64::from_bits(reservoir.to_bits() - 1);
    let radiation = bath(0.0, reservoir)
        .with_gray_radiation(Area::new(2.0), 0.5)
        .unwrap();
    assert_eq!(radiation.effective_radiating_area().value(), 1.0);
    let stage = HeatStage::implicit(
        Sum::default(),
        Some(HeatExchange::Bath(&radiation)),
        0.25,
        0.125,
        0.0,
    )
    .unwrap();
    let (heat, tangent) = stage.heat(temperature).unwrap();
    let expected_tangent = 0.5 * 4.0 * SIGMA * reservoir.powi(3);
    let expected_heat = expected_tangent * (reservoir - temperature);
    assert!((heat.value() / expected_heat - 1.0).abs() < 8e-15);
    assert!((tangent / expected_tangent - 1.0).abs() < 8e-15);
    assert_eq!(
        heat.value(),
        radiation
            .increment_density(0.125, temperature, 0.25)
            .unwrap()
    );
    let (equilibrium_heat, equilibrium_tangent) = stage.heat(reservoir).unwrap();
    assert_eq!(equilibrium_heat.value(), 0.0);
    assert!(equilibrium_tangent > 0.0);

    let linear = bath(37.0, 600.0);
    for disabled in [
        linear.with_gray_radiation(Area::new(1.0), 0.0).unwrap(),
        linear.with_gray_radiation(Area::new(0.0), 1.0).unwrap(),
    ] {
        assert_eq!(disabled, linear);
        let original = linear.increment_density(0.003, 500.0, 0.17).unwrap();
        assert_eq!(
            disabled
                .increment_density(0.003, 500.0, 0.17)
                .unwrap()
                .to_bits(),
            original.to_bits()
        );
    }
    for area in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(linear.with_gray_radiation(Area::new(area), 0.5).is_err());
    }
    for emissivity in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(
            linear
                .with_gray_radiation(Area::new(1.0), emissivity)
                .is_err()
        );
    }
    assert!(
        bath(0.0, 1e160)
            .with_gray_radiation(Area::new(1.0), 1.0)
            .unwrap()
            .increment_density(1.0, 300.0, 1.0)
            .is_err()
    );
}
