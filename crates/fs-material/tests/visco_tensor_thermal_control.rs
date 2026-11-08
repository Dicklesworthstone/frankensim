//! Thermal creep and mixed loading against independent original-source equations.
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_material::visco_tensor::stress_control::{
    PreparedTensorMixedControl as Mixed, PreparedTensorStressControl as Stress,
    TensorStressControlBudget as Budget,
};
use fs_material::visco_tensor::{
    TemperatureRamp as Ramp, TemperatureRampControls as Controls, TensorProny, TensorPronyBranch,
    TensorPronyError,
};
use fs_material::{Tangent6, Voigt};

const PLANE: [bool; 6] = [false, false, true, false, true, true];

fn diagonal(values: Voigt) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 80.0,
        },
        300.0,
        [280.0, 340.0],
    )
    .unwrap()
}
fn ramp(start_k: f64, end_k: f64, duration_s: f64) -> Ramp {
    Ramp {
        start_k,
        end_k,
        duration_s,
    }
}
fn controls(limit: f64) -> Controls {
    Controls {
        maximum_log_shift_change: limit,
        maximum_substeps: 4096,
    }
}
fn sls(limit: f64) -> TensorProny {
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = 40.0;
    TensorProny::new_mandel(
        diagonal([400.0; 6]),
        vec![TensorPronyBranch::from_gram_rows(rows, 0.6).unwrap()],
        limit,
    )
    .unwrap()
}
// Independent empirical formula; neither the calibration evaluator nor the
// production midpoint partition participates in this continuous-time oracle.
fn inverse_shift(temperature: f64) -> f64 {
    10.0_f64.powf(2.0 * (temperature - 300.0) / (80.0 + temperature - 300.0))
}
fn exposure(ramp: Ramp) -> f64 {
    let n = 8192;
    let sum: f64 = (0..=n)
        .map(|j| {
            let weight = if j == 0 || j == n {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * inverse_shift(ramp.start_k + (ramp.end_k - ramp.start_k) * j as f64 / n as f64)
        })
        .sum();
    ramp.duration_s * sum / (3.0 * n as f64)
}

#[test]
fn reference_temperature_preserves_original_stress_and_mixed_response_bits() {
    let mut initial = sls(0.1).history();
    initial
        .advance([0.001, -0.0003, 0.0002, 0.0004, 0.0, 0.0], 0.23, || false)
        .unwrap();
    for mask in [[true; 6], PLANE, [false; 6]] {
        let control = Mixed::new(&initial, mask, Budget::default(), || false).unwrap();
        let target = core::array::from_fn(|i| if mask[i] { 0.7 } else { -0.0002 });
        for duration in [0.0, 0.31] {
            let mut native = initial.clone();
            let mut thermal = initial.clone();
            let expected = control
                .advance(&mut native, target, duration, || false)
                .unwrap();
            let actual = control
                .advance_temperature_ramp(
                    &mut thermal,
                    target,
                    ramp(300.0, 300.0, duration),
                    &shift(),
                    controls(0.1),
                    || false,
                )
                .unwrap();
            assert_eq!(actual.response, expected);
            assert_eq!(thermal, native);
            assert_eq!(actual.substeps, 1);
            assert_eq!(actual.maximum_log_shift_change, 0.0);
        }
    }
}

#[test]
fn isothermal_zener_creep_and_recovery_keep_physical_work_and_heat() {
    let mut history = sls(0.1).history();
    let control = Stress::new(&history, Budget::default(), || false).unwrap();
    let target = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let jump = control
        .advance_temperature_ramp(
            &mut history,
            target,
            ramp(320.0, 320.0, 0.0),
            &shift(),
            controls(0.1),
            || false,
        )
        .unwrap()
        .response;
    close(jump.strain[0], 0.001, 2e-15);
    close(jump.work, 0.001, 2e-15);
    assert_eq!(jump.dissipated, 0.0);
    let duration = 0.7;
    let decay = (-duration * inverse_shift(320.0) / 3.0).exp();
    let increase = 0.004 * (1.0 - decay);
    let creep = control
        .advance_temperature_ramp(
            &mut history,
            target,
            ramp(320.0, 320.0, duration),
            &shift(),
            controls(0.1),
            || false,
        )
        .unwrap()
        .response;
    close(creep.strain[0], 0.001 + increase, 2e-14);
    close(creep.strain_increment[0], increase, 2e-14);
    close(creep.work, 2.0 * increase, 2e-14);
    close(creep.dissipated, 0.004 * (1.0 - decay * decay), 2e-14);
    close(creep.balance_residual, 0.0, 2e-14);
    control
        .advance(&mut history, [0.0; 6], 0.0, || false)
        .unwrap();
    close(history.strain()[0], increase, 2e-14);
    let recovery = control
        .advance_temperature_ramp(
            &mut history,
            [0.0; 6],
            ramp(290.0, 290.0, 0.8),
            &shift(),
            controls(0.1),
            || false,
        )
        .unwrap()
        .response;
    let recovery_decay = (-0.8 * inverse_shift(290.0) / 3.0).exp();
    close(recovery.strain[0], increase * recovery_decay, 2e-14);
    close(recovery.work, 0.0, 2e-15);
    close(
        recovery.dissipated,
        250.0 * increase * increase * (1.0 - recovery_decay * recovery_decay),
        2e-14,
    );
    close(
        history.work(),
        history.stored_energy() + history.dissipated(),
        4e-14,
    );
}

#[test]
fn thermal_creep_converges_to_the_continuous_reduced_time_closed_form() {
    let mut initial = sls(0.1).history();
    let control = Stress::new(&initial, Budget::default(), || false).unwrap();
    let target = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut initial, target, 0.0, || false)
        .unwrap();
    let thermal = ramp(290.0, 330.0, 1.2);
    let decay = (-exposure(thermal) / 3.0).exp();
    let expected_strain = 0.005 - 0.004 * decay;
    let expected_heat = 0.004 * (1.0 - decay * decay);
    let mut errors = [[0.0; 2]; 3];
    for (i, limit) in [0.3, 0.15, 0.075].into_iter().enumerate() {
        let mut history = initial.clone();
        let report = control
            .advance_temperature_ramp(
                &mut history,
                target,
                thermal,
                &shift(),
                controls(limit),
                || false,
            )
            .unwrap();
        let response = &report.response;
        errors[i] = [
            (response.strain[0] - expected_strain).abs(),
            (response.dissipated - expected_heat).abs(),
        ];
        close(response.work, 2.0 * response.strain_increment[0], 3e-14);
        close(response.balance_residual, 0.0, 3e-14);
        assert!(response.strain_bounds[0][0] <= initial.strain()[0]);
        assert!(response.strain_bounds[1][0] >= response.strain[0]);
        assert!(report.maximum_log_shift_change <= limit);
        assert_eq!(report.end.temperature_k(), thermal.end_k);
        let mut replay = initial.clone();
        let repeated = control
            .advance_temperature_ramp(
                &mut replay,
                target,
                thermal,
                &shift(),
                controls(limit),
                || false,
            )
            .unwrap();
        assert_eq!(repeated, report);
        assert_eq!(replay, history);
    }
    for component in 0..2 {
        assert!(errors[1][component] < 0.30 * errors[0][component]);
        assert!(errors[2][component] < 0.30 * errors[1][component]);
        assert!(errors[2][component] < 6e-8);
    }
}

// Original rank-one tensor equations with a free zz surface:
// q' = zeta*(30*ex' + sqrt(2)*11*exy') - zeta*q/(a_T*tau).
// Work uses both tensor-shear ports; heat independently integrates q²/(a_T*tau).
fn mixed_ode(y: [f64; 3], time: f64, thermal: Ramp, start: [f64; 2], end: [f64; 2]) -> [f64; 3] {
    let rate = core::array::from_fn::<_, 2, _>(|i| (end[i] - start[i]) / thermal.duration_s);
    let strain = core::array::from_fn::<_, 2, _>(|i| start[i] + time * rate[i]);
    let zeta = 600.0 / (600.0 + 24.0 * 24.0);
    let drive = 30.0 * rate[0] + std::f64::consts::SQRT_2 * 11.0 * rate[1];
    let kinetics = inverse_shift(
        thermal.start_k + (thermal.end_k - thermal.start_k) * time / thermal.duration_s,
    );
    [
        zeta * drive - zeta * kinetics * y[0] / 0.4,
        1200.0 * strain[0] * rate[0] + 1600.0 * strain[1] * rate[1] + drive * y[0],
        y[0] * y[0] * kinetics / 0.4,
    ]
}
fn integrate_mixed(y: &mut [f64; 3], thermal: Ramp, start: [f64; 2], end: [f64; 2]) {
    let h = thermal.duration_s / 8192.0;
    for j in 0..8192 {
        let time = j as f64 * h;
        let k1 = mixed_ode(*y, time, thermal, start, end);
        let k2 = mixed_ode(
            core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]),
            time + 0.5 * h,
            thermal,
            start,
            end,
        );
        let k3 = mixed_ode(
            core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]),
            time + 0.5 * h,
            thermal,
            start,
            end,
        );
        let k4 = mixed_ode(
            core::array::from_fn(|i| y[i] + h * k3[i]),
            time + h,
            thermal,
            start,
            end,
        );
        for i in 0..3 {
            y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
}

#[test]
fn mixed_thermal_reversals_match_original_physical_time_equations() {
    let mut rows = [[0.0; 6]; 6];
    rows[0] = [30.0, 0.0, 24.0, 11.0, 0.0, 0.0];
    let law = TensorProny::new_mandel(
        diagonal([1200.0, 900.0, 600.0, 800.0, 700.0, 1100.0]),
        vec![TensorPronyBranch::from_gram_rows(rows, 0.4).unwrap()],
        0.1,
    )
    .unwrap();
    let mut history = law.history();
    let control = Mixed::new(&history, PLANE, Budget::default(), || false).unwrap();
    let mut oracle = [0.0; 3];
    let mut start = [0.0; 2];
    let mut temperature = 300.0;
    for (end, next_temperature, duration) in [
        ([0.0015, 0.0004], 325.0, 0.37),
        ([-0.0008, -0.0003], 290.0, 0.53),
        ([0.0005, -0.0001], 310.0, 0.29),
        ([0.0, 0.0], 300.0, 0.41),
    ] {
        let thermal = ramp(temperature, next_temperature, duration);
        let old = history.clone();
        let response = control
            .advance_temperature_ramp(
                &mut history,
                [end[0], 0.0, 0.0, end[1], 0.0, 0.0],
                thermal,
                &shift(),
                controls(0.005),
                || false,
            )
            .unwrap()
            .response;
        integrate_mixed(&mut oracle, thermal, start, end);
        close(response.stress[0], 1200.0 * end[0] + 30.0 * oracle[0], 2e-6);
        close(
            response.stress[3],
            800.0 * end[1] + 11.0 * oracle[0] / std::f64::consts::SQRT_2,
            1e-6,
        );
        close(response.strain[2], -24.0 * oracle[0] / 600.0, 3e-9);
        close(response.strain[0], end[0], 2e-17);
        close(response.strain[3], end[1], 2e-17);
        for i in [2, 4, 5] {
            assert_eq!(response.stress[i], 0.0);
        }
        close(history.work(), oracle[1], 1e-8);
        close(history.dissipated(), oracle[2], 1e-8);
        close(response.balance_residual, 0.0, 2e-13);
        assert!(response.dissipated >= 0.0);
        for i in 0..6 {
            assert!(response.strain_bounds[0][i] <= old.strain()[i]);
            assert!(response.strain_bounds[1][i] >= old.strain()[i]);
        }
        start = end;
        temperature = next_temperature;
    }
    assert!(history.branch_stored_energy() > 0.0);
}

#[test]
fn thermal_control_refusals_preserve_the_complete_accepted_history() {
    let mut history = sls(0.1).history();
    let control = Stress::new(&history, Budget::default(), || false).unwrap();
    let target = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let before = history.clone();
    for (thermal, budget) in [
        (
            ramp(290.0, 330.0, 1.0),
            Controls {
                maximum_substeps: 1,
                ..controls(0.01)
            },
        ),
        (ramp(300.0, 350.0, 1.0), controls(0.1)),
        (ramp(300.0, 310.0, 0.0), controls(0.1)),
    ] {
        assert!(control.advance_temperature_ramp(
            &mut history, target, thermal, &shift(), budget, || false,
        ).is_err());
        assert_eq!(history, before);
    }
    let mut polls = 0;
    let cancelled = control.advance_temperature_ramp(
        &mut history,
        target,
        ramp(290.0, 330.0, 1.0),
        &shift(),
        controls(0.01),
        || {
            polls += 1;
            polls == 60
        },
    );
    assert_eq!(cancelled, Err(TensorPronyError::Cancelled));
    assert_eq!(history, before);

    let mut limited = sls(0.0015).history();
    let limited_control = Stress::new(&limited, Budget::default(), || false).unwrap();
    limited_control
        .advance(&mut limited, target, 0.0, || false)
        .unwrap();
    let original = limited.clone();
    assert!(
        limited_control
            .advance_temperature_ramp(
                &mut limited,
                target,
                ramp(300.0, 330.0, 5.0),
                &shift(),
                controls(0.03),
                || false,
            )
            .is_err()
    );
    assert_eq!(limited, original);

    let extremes = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 200.0,
            c2_k: 10.0,
        },
        300.0,
        [295.0, 305.0],
    )
    .unwrap();
    for (temperature, duration) in [(295.0, 1e-200), (305.0, 1e300)] {
        let error = control.advance_temperature_ramp(
            &mut history,
            target,
            ramp(temperature, temperature, duration),
            &extremes,
            controls(0.1),
            || false,
        );
        assert_eq!(
            error,
            Err(TensorPronyError::Invalid(
                "unrepresentable reduced-time constitutive interval",
            ))
        );
        assert_eq!(history, before);
    }
}
