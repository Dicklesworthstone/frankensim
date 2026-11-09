//! Thermal loading uses the original physical-coordinate ODE as its oracle.

use super::*;
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_material::visco_tensor::{TemperatureRamp, TemperatureRampControls};

fn calibration() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 100.0,
        },
        300.0,
        [250.0, 350.0],
    )
    .unwrap()
}

fn ramp(start_k: f64, end_k: f64, duration_s: f64) -> TemperatureRamp {
    TemperatureRamp {
        start_k,
        end_k,
        duration_s,
    }
}

fn controls(limit: f64) -> TemperatureRampControls {
    TemperatureRampControls {
        maximum_log_shift_change: limit,
        maximum_substeps: 10_000,
    }
}

// Independent scalar WLF evaluation, without the production calibration owner.
fn factor(temperature: f64) -> f64 {
    10.0_f64.powf(-2.0 * (temperature - 300.0) / (100.0 + temperature - 300.0))
}

#[test]
fn g1_constant_temperature_matches_shifted_full_law_and_reference_bits() {
    let (law, _, _, _) = coupled();
    let calibration = calibration();
    for temperature in [275.0, 300.0, 335.0] {
        let shifted = calibration.shift_tensor(&law, temperature).unwrap();
        let mut expected = shifted.law().history();
        let reference = Controller::new(&expected, Budget::default(), || false).unwrap();
        let mut history = law.history();
        let controller = Controller::new(&history, Budget::default(), || false).unwrap();
        for (stress, dt) in [
            ([1.5, -0.8, 0.6, 0.7, -0.5, 0.4], 0.0),
            ([1.5, -0.8, 0.6, 0.7, -0.5, 0.4], 0.42),
            ([-0.7, 1.1, -0.3, -0.4, 0.5, 0.2], 0.61),
            ([0.0; 6], 0.35),
        ] {
            let direct = reference
                .advance(&mut expected, stress, dt, || false)
                .unwrap();
            let actual = controller
                .advance_temperature_ramp(
                    &mut history,
                    stress,
                    ramp(temperature, temperature, dt),
                    &calibration,
                    controls(0.05),
                    || false,
                )
                .unwrap();
            assert_eq!(actual.substeps, 1);
            assert_eq!(actual.maximum_log_shift_change, 0.0);
            assert_eq!(actual.start.temperature_k(), temperature);
            assert_eq!(actual.end.temperature_k(), temperature);
            for i in 0..6 {
                close(actual.response.strain[i], direct.strain[i], 2e-13);
                close(actual.response.strain_rate[i], direct.strain_rate[i], 3e-12);
            }
            close(actual.response.work, direct.work, 3e-13);
            close(actual.response.dissipated, direct.dissipated, 3e-13);
            close(
                history.branch_stored_energy(),
                expected.branch_stored_energy(),
                3e-13,
            );
            if temperature == 300.0 {
                assert_eq!(actual.response, direct);
                assert_eq!(history, expected);
            }
        }
    }
}

#[test]
fn g1_endpoint_rate_uses_end_temperature_with_midpoint_history() {
    let (e, viscosity, stress, dt) = (400.0, 75.0, 0.8, 0.2);
    let law = kelvin(e, viscosity);
    let mut history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    let result = controller
        .advance_temperature_ramp(
            &mut history,
            [stress, 0.0, 0.0, 0.0, 0.0, 0.0],
            ramp(300.0, 340.0, dt),
            &calibration(),
            controls(10.0),
            || false,
        )
        .unwrap();
    assert_eq!(result.substeps, 1);
    let r = dt * e / (viscosity * factor(320.0));
    let strain = stress / e * (1.0 + (-r).exp_m1() / r);
    let end_rate = (stress - e * strain) / (viscosity * factor(340.0));
    relative(result.response.strain[0], strain, 3e-13);
    relative(result.response.strain_rate[0], end_rate, 3e-13);
    let midpoint_rate = (stress - e * strain) / (viscosity * factor(320.0));
    assert!((result.response.strain_rate[0] - midpoint_rate).abs() > 1e-3);
    assert!(result.response.dissipated > 0.0);
    close(result.response.balance_residual, 0.0, 2e-17);
}

#[allow(clippy::too_many_arguments)] // Independent original equations, physical temperature and loading.
fn integrate_temperature(
    y: &mut [f64; 20],
    start: Voigt,
    target: Voigt,
    thermal: TemperatureRamp,
    count: usize,
    eq: &Tangent6,
    viscosity: &Tangent6,
    branches: &[Tangent6; 2],
) -> [Voigt; 2] {
    let h = thermal.duration_s / count as f64;
    let mut bounds = [core::array::from_fn(|i| y[i]); 2];
    for step in 0..count {
        let rhs = |state: &[f64; 20], offset: f64| {
            let fraction = (step as f64 + offset) / count as f64;
            let temperature = thermal.start_k + fraction * (thermal.end_k - thermal.start_k);
            let stress = core::array::from_fn(|i| start[i] + fraction * (target[i] - start[i]));
            // Scaling BOTH physical viscosities and Maxwell times divides the
            // original strain/memory/work/heat derivatives by a_T. This uses
            // direct partial-pivot elimination, not the prepared mode kernel.
            derivative(state, stress, eq, viscosity, branches).map(|v| v / factor(temperature))
        };
        let k1 = rhs(y, 0.0);
        let k2 = rhs(&core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]), 0.5);
        let k3 = rhs(&core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]), 0.5);
        let k4 = rhs(&core::array::from_fn(|i| y[i] + h * k3[i]), 1.0);
        for i in 0..20 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        for i in 0..6 {
            bounds[0][i] = bounds[0][i].min(y[i]);
            bounds[1][i] = bounds[1][i].max(y[i]);
        }
    }
    bounds
}

#[test]
fn g1_noncommuting_thermal_reversals_converge_to_physical_ode_and_keep_memory() {
    let (law, eq, viscosity, branches) = coupled();
    let calibration = calibration();
    let cases = [
        ([1.5, -0.8, 0.6, 0.7, -0.5, 0.4], ramp(300.0, 335.0, 0.42)),
        ([-0.7, 1.1, -0.3, -0.4, 0.5, 0.2], ramp(335.0, 280.0, 0.61)),
        ([0.0; 6], ramp(280.0, 300.0, 0.35)),
    ];
    let mut oracle = [0.0; 20];
    let mut observations = Vec::new();
    let mut start = [0.0; 6];
    for (target, thermal) in cases {
        integrate_temperature(
            &mut oracle,
            start,
            target,
            thermal,
            2048,
            &eq,
            &viscosity,
            &branches,
        );
        observations.push(oracle);
        start = target;
    }
    let mut errors = [0.0_f64; 3];
    for (rung, limit) in [0.1, 0.05, 0.025].into_iter().enumerate() {
        let mut history = law.history();
        let controller = Controller::new(&history, Budget::default(), || false).unwrap();
        for ((target, thermal), expected) in cases.into_iter().zip(&observations) {
            let result = controller
                .advance_temperature_ramp(
                    &mut history,
                    tensor(target),
                    thermal,
                    &calibration,
                    controls(limit),
                    || false,
                )
                .unwrap();
            assert!(result.substeps > 1);
            assert!(result.maximum_log_shift_change <= limit);
            assert!(result.response.dissipated >= 0.0);
            close(result.response.balance_residual, 0.0, 2e-15);
            for i in 0..6 {
                errors[rung] = errors[rung].max((mandel(history.strain())[i] - expected[i]).abs());
            }
            if rung == 2 {
                close(history.work(), expected[18], 2e-7);
                close(history.dissipated(), expected[19], 2e-7);
                let rate = derivative(expected, target, &eq, &viscosity, &branches);
                for i in 0..6 {
                    close(
                        mandel(result.response.strain_rate)[i],
                        rate[i] / factor(thermal.end_k),
                        3e-6,
                    );
                }
            }
        }
        if rung == 2 {
            // The thermal cycle returns to the original reference law. Its
            // ORIGINAL memory must remain usable by prescribed-strain control.
            let next_strain = [0.0007, 0.0001, -0.0003, 0.0004, -0.0005, 0.0002];
            let h = 0.31;
            let next_mandel = mandel(next_strain);
            let rate = core::array::from_fn(|i| (next_mandel[i] - oracle[i]) / h);
            let equilibrium = apply(&eq, next_mandel);
            let kelvin = apply(&viscosity, rate);
            let mut stress: Voigt = core::array::from_fn(|i| equilibrium[i] + kelvin[i]);
            for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
                let decay = (-h / tau).exp();
                let q = core::array::from_fn(|i| {
                    decay * oracle[6 + 6 * j + i] + tau * (1.0 - decay) * rate[i]
                });
                let branch = apply(&branches[j], q);
                for i in 0..6 {
                    stress[i] += branch[i];
                }
            }
            let next = history.advance(next_strain, h, || false).unwrap();
            for i in 0..6 {
                close(mandel(next.stress)[i], stress[i], 2e-4);
            }
        }
    }
    assert!(
        errors[1] < errors[0] / 2.5 && errors[2] < errors[1] / 2.5,
        "{errors:?}"
    );
    assert!(errors[2] < 2e-7, "{errors:?}");
}

#[test]
fn g0_thermal_budget_binding_and_late_cancellation_are_atomic() {
    let (law, _, _, _) = coupled();
    let mut history = law.history();
    history
        .advance(
            [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005],
            0.4,
            || false,
        )
        .unwrap();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    let before = history.clone();
    let calibration = calibration();
    let target = [-0.5, 0.7, -0.2, 0.3, 0.4, -0.6];
    let thermal = ramp(300.0, 335.0, 0.3);
    let mut calls = 0;
    let mut complete = history.clone();
    controller
        .advance_temperature_ramp(
            &mut complete,
            target,
            thermal,
            &calibration,
            controls(0.05),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 50);
    for stop in [1, calls / 2, calls - 1, calls] {
        let mut count = 0;
        let result = controller.advance_temperature_ramp(
            &mut history,
            target,
            thermal,
            &calibration,
            controls(0.05),
            || {
                count += 1;
                count == stop
            },
        );
        assert_eq!(result, Err(TensorPronyError::Cancelled));
        assert_eq!(history, before);
    }
    for (thermal, control) in [
        (
            thermal,
            TemperatureRampControls {
                maximum_substeps: 1,
                ..controls(0.05)
            },
        ),
        (thermal, controls(0.0)),
        (ramp(300.0, 360.0, 0.3), controls(0.05)),
        (ramp(300.0, 335.0, 0.0), controls(0.05)),
        (ramp(300.0, 300.0, f64::INFINITY), controls(0.05)),
    ] {
        assert!(
            controller
                .advance_temperature_ramp(
                    &mut history,
                    target,
                    thermal,
                    &calibration,
                    control,
                    || false,
                )
                .is_err()
        );
        assert_eq!(history, before);
    }
    let mut unrelated = law.history();
    let saved = unrelated.clone();
    assert!(
        controller
            .advance_temperature_ramp(
                &mut unrelated,
                target,
                thermal,
                &calibration,
                controls(0.05),
                || false,
            )
            .is_err()
    );
    assert_eq!(unrelated, saved);
}

#[test]
fn g0_physical_endpoint_rate_overflow_refuses_after_admissible_reference_jump() {
    let law = kelvin(1.0, 1.0);
    let mut history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    let before = history.clone();
    let calibration = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 300.0,
            c2_k: 100.0,
        },
        300.0,
        [300.0, 350.0],
    )
    .unwrap();
    let target = [1e300, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut reference = history.clone();
    let direct = controller
        .advance(&mut reference, target, 0.0, || false)
        .unwrap();
    assert!(direct.strain_rate[0].is_finite());
    assert!((direct.strain_rate[0] / calibration.evaluate(350.0).unwrap().factor()).is_infinite());
    assert_eq!(
        controller.advance_temperature_ramp(
            &mut history,
            target,
            ramp(350.0, 350.0, 0.0),
            &calibration,
            controls(0.05),
            || false,
        ),
        Err(TensorPronyError::Invalid(
            "unrepresentable physical Kelvin endpoint strain rate"
        ))
    );
    assert_eq!(history, before);
}
