//! Original differential/algebraic equations, without an invented viscosity.

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

fn controls(limit: f64) -> TemperatureRampControls {
    TemperatureRampControls {
        maximum_log_shift_change: limit,
        maximum_substeps: 10_000,
    }
}

#[test]
fn g1_directional_kelvin_jumps_ramps_and_shear_keep_exact_elastic_nullspace() {
    let e = [400.0, 600.0, 800.0, 300.0, 500.0, 700.0];
    let v = [75.0, 0.0, 0.0, 25.0, 0.0, 0.0];
    let law = TensorProny::new_mandel(diagonal(e), vec![], 0.1)
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal(v))
        .unwrap();
    let mut history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    assert_eq!(controller.viscosity_rank(), 2);
    assert_eq!(controller.state_coordinates(), 6);
    let stress = [1.4, -0.9, 0.8, 1.3, -0.7, 0.6];
    let jump = controller
        .advance(&mut history, tensor(stress), 0.0, || false)
        .unwrap();
    let mut jump_work = 0.0;
    for i in 0..6 {
        let expected = if v[i] == 0.0 { stress[i] / e[i] } else { 0.0 };
        close(mandel(jump.strain)[i], expected, 2e-17);
        close(
            mandel(jump.strain_rate)[i],
            if v[i] == 0.0 { 0.0 } else { stress[i] / v[i] },
            2e-16,
        );
        jump_work += 0.5 * stress[i] * expected;
    }
    close(jump.work, jump_work, 2e-17);
    close(jump.stored_energy_change, jump_work, 2e-17);
    assert_eq!(jump.dissipated, 0.0);
    let dt = 0.4;
    let hold = controller
        .advance(&mut history, tensor(stress), dt, || false)
        .unwrap();
    let mut heat = 0.0;
    let mut work = 0.0;
    for i in 0..6 {
        let expected = if v[i] == 0.0 {
            stress[i] / e[i]
        } else {
            let r = dt * e[i] / v[i];
            heat += stress[i] * stress[i] / (2.0 * e[i]) * -(-2.0 * r).exp_m1();
            work += stress[i] * stress[i] / e[i] * -(-r).exp_m1();
            stress[i] / e[i] * -(-r).exp_m1()
        };
        close(mandel(hold.strain)[i], expected, 2e-16);
    }
    close(hold.dissipated, heat, 3e-16);
    close(hold.work, work, 3e-16);
    let target = [-0.5, 0.7, -0.3, -0.8, 0.5, -0.4];
    let dt = 0.3;
    let old = mandel(history.strain());
    let result = controller
        .advance(&mut history, tensor(target), dt, || false)
        .unwrap();
    let mut heat = 0.0;
    for i in 0..6 {
        let (expected, rate) = if v[i] == 0.0 {
            (target[i] / e[i], (target[i] - stress[i]) / (e[i] * dt))
        } else {
            let lambda = e[i] / v[i];
            let r = lambda * dt;
            let decay = (-r).exp();
            let ramp = -(-r).exp_m1() / r;
            let expected = decay * old[i]
                + stress[i] / e[i] * (1.0 - decay)
                + (target[i] - stress[i]) / e[i] * (1.0 - ramp);
            let constant = (target[i] - stress[i]) / (e[i] * dt);
            let transient = (stress[i] / e[i] - old[i]) * lambda - constant;
            heat += v[i]
                * (constant * constant * dt
                    + 2.0 * constant * transient * -(-r).exp_m1() / lambda
                    + transient * transient * -(-2.0 * r).exp_m1() / (2.0 * lambda));
            (expected, (target[i] - e[i] * expected) / v[i])
        };
        close(mandel(result.strain)[i], expected, 4e-16);
        close(mandel(result.strain_rate)[i], rate, 4e-15);
    }
    close(result.dissipated, heat, 3e-16);
    close(result.balance_residual, 0.0, 3e-16);
}

#[test]
fn g1_shear_only_kelvin_retains_weak_positive_channel_without_a_floor() {
    for shear_viscosity in [1e-120, 35.0, 1e120] {
        for axial_viscosity in [0.0, 75.0] {
            let mut viscosity = [0.0; 6];
            viscosity[0] = axial_viscosity;
            viscosity[3] = shear_viscosity;
            let law = TensorProny::new_mandel(diagonal([400.0; 6]), vec![], 0.1)
                .unwrap()
                .with_kelvin_viscosity_mandel(diagonal(viscosity))
                .unwrap();
            let mut history = law.history();
            let controller = Controller::new(&history, Budget::default(), || false).unwrap();
            assert_eq!(
                controller.viscosity_rank(),
                if axial_viscosity == 0.0 { 1 } else { 2 }
            );
            let mut stress = [0.0; 6];
            stress[3] = 0.8;
            let jump = controller
                .advance(&mut history, tensor(stress), 0.0, || false)
                .unwrap();
            assert_eq!(jump.strain, [0.0; 6]);
            assert_eq!(jump.work, 0.0);
            let result = controller
                .advance(
                    &mut history,
                    tensor(stress),
                    shear_viscosity / 400.0,
                    || false,
                )
                .unwrap();
            relative(
                mandel(result.strain)[3],
                0.8 / 400.0 * -(-1.0_f64).exp_m1(),
                4e-13,
            );
            relative(
                result.dissipated,
                0.8 * 0.8 / 800.0 * -(-2.0_f64).exp_m1(),
                4e-13,
            );
            for i in [0, 1, 2, 4, 5] {
                assert_eq!(result.strain[i], 0.0);
                assert_eq!(result.strain_rate[i], 0.0);
            }
        }
    }
}

fn directional() -> (TensorProny, Tangent6, Tangent6, [Tangent6; 2]) {
    let (law, eq, _, branches) = coupled();
    let viscosity = diagonal([90.0, 0.0, 0.0, 55.0, 0.0, 0.0]);
    let law = TensorProny::new_mandel(eq, law.branches().to_vec(), 0.1)
        .unwrap()
        .with_kelvin_viscosity_mandel(viscosity)
        .unwrap();
    (law, eq, viscosity, branches)
}

// An independent index-one DAE solve in the ORIGINAL strain and branch-strain
// coordinates. Viscous rows use stress balance; exactly nonviscous rows use
// its derivative. No production null basis, compliance, SVD or mode is used.
fn dae_derivative(
    y: &[f64; 20],
    stress: Voigt,
    slope: Voigt,
    scale: f64,
    eq: &Tangent6,
    viscosity: &Tangent6,
    branches: &[Tangent6; 2],
) -> [f64; 20] {
    let strain = core::array::from_fn(|i| y[i]);
    let elastic = apply(eq, strain);
    let mut rhs = core::array::from_fn(|i| stress[i] - elastic[i]);
    let mut released = [0.0; 6];
    let mut heat = 0.0;
    for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
        let q = core::array::from_fn(|i| y[6 + 6 * j + i]);
        let sigma = apply(&branches[j], q);
        for i in 0..6 {
            rhs[i] -= sigma[i];
            released[i] += sigma[i] / (tau * scale);
        }
        heat += dot(q, sigma) / (tau * scale);
    }
    let matrix = core::array::from_fn(|i| {
        if viscosity[i][i] != 0.0 {
            viscosity[i].map(|v| v * scale)
        } else {
            rhs[i] = slope[i] + released[i];
            core::array::from_fn(|j| eq[i][j] + branches[0][i][j] + branches[1][i][j])
        }
    });
    let rate = solve(matrix, rhs);
    let mut result = [0.0; 20];
    result[..6].copy_from_slice(&rate);
    for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
        for i in 0..6 {
            result[6 + 6 * j + i] = rate[i] - y[6 + 6 * j + i] / (tau * scale);
        }
    }
    result[18] = dot(stress, rate);
    result[19] = heat + scale * dot(rate, apply(viscosity, rate));
    result
}

fn shift(temperature: f64) -> f64 {
    10.0_f64.powf(-2.0 * (temperature - 300.0) / (100.0 + temperature - 300.0))
}

#[allow(clippy::too_many_arguments)]
fn integrate_dae(
    y: &mut [f64; 20],
    start: Voigt,
    target: Voigt,
    dt: f64,
    count: usize,
    temperatures: [f64; 2],
    eq: &Tangent6,
    viscosity: &Tangent6,
    branches: &[Tangent6; 2],
) {
    let h = dt / count as f64;
    let slope = core::array::from_fn(|i| (target[i] - start[i]) / dt);
    for step in 0..count {
        let rhs = |state: &[f64; 20], offset: f64| {
            let t = (step as f64 + offset) / count as f64;
            let stress = core::array::from_fn(|i| start[i] + t * (target[i] - start[i]));
            let scale = shift(temperatures[0] + t * (temperatures[1] - temperatures[0]));
            dae_derivative(state, stress, slope, scale, eq, viscosity, branches)
        };
        let k1 = rhs(y, 0.0);
        let k2 = rhs(&core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]), 0.5);
        let k3 = rhs(&core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]), 0.5);
        let k4 = rhs(&core::array::from_fn(|i| y[i] + h * k3[i]), 1.0);
        for i in 0..20 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
}

#[test]
fn g1_noncommuting_semidefinite_maxwell_kelvin_matches_physical_dae_refinement() {
    let (law, eq, viscosity, branches) = directional();
    let mut history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    assert_eq!(controller.state_coordinates(), 18);
    assert_eq!(controller.viscosity_rank(), 2);
    let cases = [
        ([1.5, -0.8, 0.6, 0.7, -0.5, 0.4], 0.42),
        ([-0.7, 1.1, -0.3, -0.4, 0.5, 0.2], 0.61),
        ([0.0; 6], 0.35),
    ];
    let mut reports = Vec::new();
    for (target, dt) in cases {
        reports.push(
            controller
                .advance(&mut history, tensor(target), dt, || false)
                .unwrap(),
        );
    }
    let mut errors = [0.0_f64; 3];
    for (rung, count) in [32, 64, 128].into_iter().enumerate() {
        let mut y = [0.0; 20];
        let mut start = [0.0; 6];
        for ((target, dt), report) in cases.into_iter().zip(&reports) {
            integrate_dae(
                &mut y, start, target, dt, count, [300.0; 2], &eq, &viscosity, &branches,
            );
            for i in 0..6 {
                errors[rung] = errors[rung].max((y[i] - mandel(report.strain)[i]).abs());
            }
            if rung == 2 {
                let slope = core::array::from_fn(|i| (target[i] - start[i]) / dt);
                let rate = dae_derivative(&y, target, slope, 1.0, &eq, &viscosity, &branches);
                for i in 0..6 {
                    close(rate[i], mandel(report.strain_rate)[i], 2e-9);
                }
                assert!(report.dissipated >= 0.0);
                close(report.balance_residual, 0.0, 5e-15);
            }
            start = target;
        }
        if rung == 2 {
            close(y[18], history.work(), 3e-10);
            close(y[19], history.dissipated(), 3e-10);
            let branch_energy = (0..2)
                .map(|j| {
                    let q = core::array::from_fn(|i| y[6 + 6 * j + i]);
                    0.5 * dot(q, apply(&branches[j], q))
                })
                .sum::<f64>();
            close(branch_energy, history.branch_stored_energy(), 3e-12);
        }
    }
    assert!(
        errors[1] < errors[0] / 8.0 && errors[2] < errors[1] / 8.0,
        "{errors:?}"
    );
}

#[test]
fn g3_semidefinite_rotation_subdivision_jumps_and_strain_switch_preserve_memory() {
    let (law, _, _, _) = directional();
    let rotation = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    let rotate = |v| fs_material::tensor::rotate(&v, &rotation);
    let mut history = law.history();
    let mut rotated = law.rotated(rotation).unwrap().history();
    let preload = [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005];
    history.advance(preload, 0.4, || false).unwrap();
    rotated.advance(rotate(preload), 0.4, || false).unwrap();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    let rotated_controller = Controller::new(&rotated, Budget::default(), || false).unwrap();
    assert_eq!(rotated_controller.viscosity_rank(), 2);
    let old_stress = mandel(history.stress());
    let old_strain = history.strain();
    let target = [-0.5, 0.7, -0.2, 0.3, 0.4, -0.6];
    let jump = controller
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let rotated_jump = rotated_controller
        .advance(&mut rotated, rotate(target), 0.0, || false)
        .unwrap();
    assert_eq!(jump.dissipated, 0.0);
    close(
        jump.work,
        0.5 * dot(
            core::array::from_fn(|i| old_stress[i] + mandel(target)[i]),
            mandel(jump.strain_increment),
        ),
        3e-15,
    );
    close(jump.work, jump.stored_energy_change, 3e-15);
    close(jump.strain[0], old_strain[0], 1e-17);
    close(jump.strain[3], old_strain[3], 1e-17);
    close(rotated_jump.work, jump.work, 3e-13);
    let mut split = history.clone();
    let target = [0.8, -0.5, 0.3, -0.2, 0.1, 0.7];
    let middle = core::array::from_fn(|i| 0.5 * (history.stress()[i] + target[i]));
    controller
        .advance(&mut history, target, 0.73, || false)
        .unwrap();
    controller
        .advance(&mut split, middle, 0.365, || false)
        .unwrap();
    controller
        .advance(&mut split, target, 0.365, || false)
        .unwrap();
    rotated_controller
        .advance(&mut rotated, rotate(target), 0.73, || false)
        .unwrap();
    for i in 0..6 {
        close(history.strain()[i], split.strain()[i], 3e-13);
        close(rotated.strain()[i], rotate(history.strain())[i], 3e-13);
    }
    close(history.work(), split.work(), 3e-13);
    close(history.dissipated(), split.dissipated(), 3e-13);
    close(history.work(), rotated.work(), 3e-13);
    close(history.dissipated(), rotated.dissipated(), 3e-13);
    let strain = [0.0007, 0.0001, -0.0003, 0.0004, -0.0005, 0.0002];
    let a = history.advance(strain, 0.31, || false).unwrap();
    let b = split.advance(strain, 0.31, || false).unwrap();
    let c = rotated.advance(rotate(strain), 0.31, || false).unwrap();
    for i in 0..6 {
        close(a.stress[i], b.stress[i], 3e-10);
        close(c.stress[i], rotate(a.stress)[i], 3e-10);
    }
}

#[test]
fn g4_semidefinite_jump_and_late_cancellation_are_atomic() {
    let (law, _, _, _) = directional();
    let history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    for (target, dt) in [([1e6; 6], 0.0), ([f64::NAN; 6], 0.1), ([1.0; 6], -1.0)] {
        let mut trial = history.clone();
        assert!(
            controller
                .advance(&mut trial, target, dt, || false)
                .is_err()
        );
        assert_eq!(trial, history);
    }
    for dt in [0.0, 0.4] {
        let mut accepted = history.clone();
        let mut polls = 0;
        controller
            .advance(&mut accepted, [0.2; 6], dt, || {
                polls += 1;
                false
            })
            .unwrap();
        let mut trial = history.clone();
        let mut actual = 0;
        let error = controller
            .advance(&mut trial, [0.2; 6], dt, || {
                actual += 1;
                actual == polls
            })
            .unwrap_err();
        assert_eq!(error, TensorPronyError::Cancelled);
        assert_eq!(trial, history);
    }
    let mut unrelated = law.history();
    assert!(
        controller
            .advance(&mut unrelated, [0.2; 6], 0.0, || false)
            .is_err()
    );
    assert_eq!(unrelated, history);
}

#[test]
fn g1_semidefinite_temperature_preserves_elastic_rates_and_shifted_law() {
    let (law, _, _, _) = directional();
    let calibration = calibration();
    for temperature in [275.0, 300.0, 335.0] {
        let shifted = calibration.shift_tensor(&law, temperature).unwrap();
        let mut direct = shifted.law().history();
        let direct_controller = Controller::new(&direct, Budget::default(), || false).unwrap();
        let mut history = law.history();
        let controller = Controller::new(&history, Budget::default(), || false).unwrap();
        for (stress, dt) in [
            ([0.8, -0.7, 0.6, -0.5, 0.4, -0.3], 0.0),
            ([-0.2, 0.3, -0.4, 0.5, -0.6, 0.7], 0.4),
        ] {
            let expected = direct_controller
                .advance(&mut direct, stress, dt, || false)
                .unwrap();
            let actual = controller
                .advance_temperature_ramp(
                    &mut history,
                    stress,
                    TemperatureRamp {
                        start_k: temperature,
                        end_k: temperature,
                        duration_s: dt,
                    },
                    &calibration,
                    controls(0.05),
                    || false,
                )
                .unwrap()
                .response;
            for i in 0..6 {
                close(actual.strain[i], expected.strain[i], 3e-13);
                close(actual.strain_rate[i], expected.strain_rate[i], 3e-11);
            }
            close(actual.work, expected.work, 3e-13);
            close(actual.dissipated, expected.dissipated, 3e-13);
            if temperature == 300.0 {
                assert_eq!(actual, expected);
                assert_eq!(history, direct);
            }
        }
    }
    let law = TensorProny::new_mandel(diagonal([400.0; 6]), vec![], 0.1)
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal([75.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        .unwrap();
    let mut history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    let result = controller
        .advance_temperature_ramp(
            &mut history,
            [0.8, 1.2, 0.0, 0.0, 0.0, 0.0],
            TemperatureRamp {
                start_k: 300.0,
                end_k: 340.0,
                duration_s: 0.2,
            },
            &calibration,
            controls(10.0),
            || false,
        )
        .unwrap();
    assert_eq!(result.substeps, 1);
    close(result.response.strain[1], 1.2 / 400.0, 2e-17);
    close(result.response.strain_rate[1], 1.2 / (400.0 * 0.2), 2e-16);
    close(
        result.response.strain_rate[0],
        (0.8 - 400.0 * result.response.strain[0]) / (75.0 * shift(340.0)),
        2e-15,
    );
}

#[test]
fn g1_noncommuting_semidefinite_thermal_reversals_converge_to_physical_dae() {
    let (law, eq, viscosity, branches) = directional();
    let cases = [
        ([1.5, -0.8, 0.6, 0.7, -0.5, 0.4], 0.42, [300.0, 335.0]),
        ([-0.7, 1.1, -0.3, -0.4, 0.5, 0.2], 0.61, [335.0, 280.0]),
    ];
    let mut oracle = [0.0; 20];
    let mut start = [0.0; 6];
    for (target, dt, temperatures) in cases {
        integrate_dae(
            &mut oracle,
            start,
            target,
            dt,
            4096,
            temperatures,
            &eq,
            &viscosity,
            &branches,
        );
        start = target;
    }
    let calibration = calibration();
    let mut errors = [0.0_f64; 3];
    for (rung, limit) in [0.12, 0.06, 0.03].into_iter().enumerate() {
        let mut history = law.history();
        let controller = Controller::new(&history, Budget::default(), || false).unwrap();
        let mut last = None;
        for (target, dt, temperatures) in cases {
            let report = controller
                .advance_temperature_ramp(
                    &mut history,
                    tensor(target),
                    TemperatureRamp {
                        start_k: temperatures[0],
                        end_k: temperatures[1],
                        duration_s: dt,
                    },
                    &calibration,
                    controls(limit),
                    || false,
                )
                .unwrap();
            assert!(report.response.dissipated >= 0.0);
            close(report.response.balance_residual, 0.0, 2e-14);
            last = Some(report.response);
        }
        for i in 0..6 {
            errors[rung] = errors[rung].max((mandel(history.strain())[i] - oracle[i]).abs());
        }
        if rung == 2 {
            close(history.work(), oracle[18], 2e-7);
            close(history.dissipated(), oracle[19], 2e-7);
            let (target, dt, temperatures) = cases[1];
            let slope = core::array::from_fn(|i| (target[i] - cases[0].0[i]) / dt);
            let rate = dae_derivative(
                &oracle,
                target,
                slope,
                shift(temperatures[1]),
                &eq,
                &viscosity,
                &branches,
            );
            let actual = mandel(last.unwrap().strain_rate);
            for i in 0..6 {
                close(actual[i], rate[i], 2e-7);
            }
        }
    }
    assert!(
        errors[1] < errors[0] / 2.5 && errors[2] < errors[1] / 2.5,
        "{errors:?}"
    );
}

#[test]
fn g4_semidefinite_physical_rate_overflow_does_not_commit_a_stress_jump() {
    let calibration = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 300.0,
            c2_k: 100.0,
        },
        300.0,
        [300.0, 350.0],
    )
    .unwrap();
    let law = TensorProny::new_mandel(diagonal([1.0; 6]), vec![], 0.1)
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        .unwrap();
    let history = law.history();
    let controller = Controller::new(&history, Budget::default(), || false).unwrap();
    let stress = [1e300, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut reference = history.clone();
    assert!(
        controller
            .advance(&mut reference, stress, 0.0, || false)
            .unwrap()
            .strain_rate[0]
            .is_finite()
    );
    let mut trial = history.clone();
    assert!(
        controller
            .advance_temperature_ramp(
                &mut trial,
                stress,
                TemperatureRamp {
                    start_k: 350.0,
                    end_k: 350.0,
                    duration_s: 0.0
                },
                &calibration,
                controls(0.1),
                || false
            )
            .is_err()
    );
    assert_eq!(trial, history);
}
