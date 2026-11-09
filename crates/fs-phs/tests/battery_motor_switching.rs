//! Native finite-source PWM drive: independent SI dynamics, all ten physical
//! heats, chemical regeneration, edge/source continuity and atomic refusal.
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadMode, RegulatedPowerLoad,
};
use fs_phs::circuit::converter::switching::battery_motor::{
    BatteryMotorSwitchCommand, PreparedBatteryMotorSwitchCycle,
};
use fs_phs::circuit::converter::switching::{
    BridgeSwitchState, SwitchCycleControls, SwitchInterval, centered_pwm,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PreparedStepError, SolveBudget};

const CAPACITY: f64 = 80.0;
const LEAD: f64 = 0.08;
const BATTERY: BatteryParameters = BatteryParameters {
    series_resistance_ohm: 0.1,
    polarization_resistance_ohm: 0.2,
    polarization_capacitance_f: 0.5,
    bus_capacitance_f: 0.03,
};
const CONVERTER: BuckBoostParameters = BuckBoostParameters {
    inductance_h: 0.04,
    bus_capacitance_f: [0.04, 0.05],
    winding_resistance_ohm: 0.07,
    input_on_resistance_ohm: [0.02, 0.03],
    output_on_resistance_ohm: [0.02, 0.03],
};
const MOTOR: DcMotorParameters = DcMotorParameters {
    inductance_h: 0.02,
    winding_resistance_ohm: 0.2,
    rotor_inertia_kg_m2: 0.08,
    stator_inertia_kg_m2: 0.4,
    coupling_nm_per_a: 0.4,
    bearing_drag_nm_s_per_rad: 0.01,
};
const MODULATION: BridgeModulation = BridgeModulation {
    input_high_fraction: 0.55,
    output_high_fraction: 0.85,
};
const HIGH: BridgeSwitchState = BridgeSwitchState {
    input_high: true,
    output_high: true,
};
const COMMAND: BatteryMotorSwitchCommand = BatteryMotorSwitchCommand {
    shaft_torque_nm: [-0.2, 0.03],
    auxiliary_load: None,
};

fn initial() -> [f64; 9] {
    [
        CAPACITY * 0.65,
        0.5 * 0.1,
        0.03 * 12.1,
        0.04 * 0.5,
        0.04 * 11.8,
        0.05 * 6.4,
        0.02 * 0.1,
        0.08 * 4.0,
        0.4 * -0.4,
    ]
}
fn model_with_curve(curve: BatteryOcvCurve, modulation: BridgeModulation) -> BatteryMotorDrive {
    AveragedBuckBoost::new(CONVERTER, modulation)
        .unwrap()
        .with_battery_source(BatterySource::new(curve, BATTERY).unwrap(), LEAD)
        .unwrap()
        .with_dc_motor(DcMotor::new(MOTOR).unwrap())
        .unwrap()
}
fn model(modulation: BridgeModulation) -> BatteryMotorDrive {
    model_with_curve(
        BatteryOcvCurve::new(CAPACITY, &[(0.0, 10.0), (1.0, 14.0)]).unwrap(),
        modulation,
    )
}
fn electrical() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-9,
        relative_energy_tolerance: 1e-11,
    }
}
fn cycle_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        absolute_energy_tolerance_j: 2e-9,
        relative_energy_tolerance: 2e-11,
        ..SwitchCycleControls::default()
    }
}
fn prepare(intervals: &[SwitchInterval], subdivisions: usize) -> PreparedBatteryMotorSwitchCycle {
    PreparedBatteryMotorSwitchCycle::new_controlled(
        model(MODULATION),
        intervals,
        subdivisions,
        electrical(),
        cycle_controls(),
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn one_static_gate_is_the_original_nine_state_owner() {
    let interval = SwitchInterval {
        state: HIGH,
        duration_s: 0.002,
    };
    let mut cycle = prepare(&[interval], 1);
    let mut original = model(BridgeModulation {
        input_high_fraction: 1.0,
        output_high_fraction: 1.0,
    })
    .prepare_controlled(interval.duration_s, electrical(), || false)
    .unwrap();
    let (mut expected, mut actual) = ([0.0; 9], [0.0; 9]);
    let direct = original
        .step_controlled(
            &initial(),
            &COMMAND.shaft_torque_nm,
            None,
            &mut expected,
            || false,
        )
        .unwrap();
    let step = cycle
        .step_controlled(&initial(), COMMAND, &mut actual, || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(step.heat_joules(), direct.heat_joules());
    assert_eq!(
        step.mean_bus_voltage_v,
        direct.electrical.discrete_bus_voltage_v
    );
    assert_eq!(step.mean_motor_current_a, direct.discrete_motor_current_a);
    assert_eq!(
        step.mean_shaft_speed_rad_s,
        direct.discrete_shaft_speed_rad_s
    );
    assert_eq!(step.shaft_work_j, direct.shaft_work_j);
    assert_eq!(
        step.converter_terminal_work_j,
        direct.electrical.converter_terminal_work_j
    );
    assert_eq!(step.native.newton_iters, direct.native.newton_iters);
    assert_eq!(cycle.cached_interval_count(), 1);
    assert_eq!(step.integration_steps, 1);
    close(step.balance_residual_j(), 0.0, 1e-9);
}

// G1 oracle: explicit KCL, KVL, shaft momentum, and SI power products.
// It uses no production port matrices, native steps or loss/account helpers.
fn derivative(x: &[f64; 27], gates: BridgeSwitchState) -> [f64; 27] {
    let ocv = 10.0 + 4.0 * x[0] / CAPACITY;
    let pol = x[1] / BATTERY.polarization_capacitance_f;
    let vb = x[2] / BATTERY.bus_capacitance_f;
    let ic = x[3] / CONVERTER.inductance_h;
    let vi = x[4] / CONVERTER.bus_capacitance_f[0];
    let vo = x[5] / CONVERTER.bus_capacitance_f[1];
    let im = x[6] / MOTOR.inductance_h;
    let wr = x[7] / MOTOR.rotor_inertia_kg_m2;
    let ws = x[8] / MOTOR.stator_inertia_kg_m2;
    let relative = wr - ws;
    let cell = (ocv - pol - vb) / BATTERY.series_resistance_ohm;
    let lead = (vb - vi) / LEAD;
    let a = if gates.input_high { 1.0 } else { 0.0 };
    let b = if gates.output_high { 1.0 } else { 0.0 };
    let resistance = [
        CONVERTER.winding_resistance_ohm,
        a * CONVERTER.input_on_resistance_ohm[0],
        (1.0 - a) * CONVERTER.input_on_resistance_ohm[1],
        b * CONVERTER.output_on_resistance_ohm[0],
        (1.0 - b) * CONVERTER.output_on_resistance_ohm[1],
    ];
    let internal_torque = MOTOR.coupling_nm_per_a * im - MOTOR.bearing_drag_nm_s_per_rad * relative;
    let mut out = [0.0; 27];
    out[..9].copy_from_slice(&[
        -cell,
        cell - pol / BATTERY.polarization_resistance_ohm,
        cell - lead,
        a * vi - b * vo - resistance.iter().sum::<f64>() * ic,
        lead - a * ic,
        b * ic - im,
        vo - MOTOR.winding_resistance_ohm * im - MOTOR.coupling_nm_per_a * relative,
        internal_torque + COMMAND.shaft_torque_nm[0],
        -internal_torque + COMMAND.shaft_torque_nm[1],
    ]);
    for i in 0..5 {
        out[9 + i] = resistance[i] * ic * ic;
    }
    out[14] = BATTERY.series_resistance_ohm * cell * cell;
    out[15] = pol * pol / BATTERY.polarization_resistance_ohm;
    out[16] = LEAD * lead * lead;
    out[17] = MOTOR.winding_resistance_ohm * im * im;
    out[18] = MOTOR.bearing_drag_nm_s_per_rad * relative * relative;
    out[19..].copy_from_slice(&[
        -vb * lead,
        vi * lead,
        -vo * im,
        vo * im,
        MOTOR.coupling_nm_per_a * relative * im,
        COMMAND.shaft_torque_nm[0] * wr,
        COMMAND.shaft_torque_nm[1] * ws,
        0.0,
    ]);
    out
}
fn reference(intervals: &[SwitchInterval], cycles: usize, subdivisions: usize) -> [f64; 27] {
    let mut x = [0.0; 27];
    x[..9].copy_from_slice(&initial());
    for _ in 0..cycles {
        for interval in intervals {
            let h = interval.duration_s / subdivisions as f64;
            for _ in 0..subdivisions {
                let a = derivative(&x, interval.state);
                let b = derivative(
                    &core::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
                    interval.state,
                );
                let c = derivative(
                    &core::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
                    interval.state,
                );
                let d = derivative(&core::array::from_fn(|i| x[i] + h * c[i]), interval.state);
                x = core::array::from_fn(|i| {
                    x[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
                });
            }
        }
    }
    x
}

#[test]
fn g1_pwm_battery_winding_and_both_shafts_match_independent_ode_and_ten_heats() {
    let intervals = centered_pwm(0.006, MODULATION).unwrap();
    let oracle = reference(&intervals, 8, 512);
    let second_oracle = reference(&intervals, 8, 256);
    for i in 0..27 {
        close(oracle[i], second_oracle[i], 1e-9);
    }
    let mut errors = Vec::new();
    for subdivisions in [4, 8, 16] {
        let mut cycle = prepare(&intervals, subdivisions);
        let mut x = initial();
        let mut heat = [0.0; 10];
        let mut work = [0.0; 8];
        for _ in 0..8 {
            let mut next = [0.0; 9];
            let step = cycle
                .step_controlled(&x, COMMAND, &mut next, || false)
                .unwrap();
            close(step.balance_residual_j(), 0.0, 2e-9);
            close(step.angular_momentum_residual_nms, 0.0, 2e-11);
            for (sum, q) in heat.iter_mut().zip(step.heat_joules()) {
                assert!(q >= 0.0);
                *sum += q;
            }
            for (sum, w) in work.iter_mut().zip([
                step.battery_terminal_work_j,
                step.converter_terminal_work_j[0],
                step.converter_terminal_work_j[1],
                step.motor_electrical_work_j,
                step.electromagnetic_work_j,
                step.shaft_work_j[0],
                step.shaft_work_j[1],
                step.load_work_j,
            ]) {
                *sum += w;
            }
            x = next;
        }
        // Fixed 1 SI unit scales for charge, flux and momentum, shared by all refinements.
        errors.push(
            x.iter()
                .zip(&oracle[..9])
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
        if subdivisions == 16 {
            for i in 0..10 {
                close(heat[i], oracle[9 + i], 2e-5);
            }
            for i in 0..8 {
                close(work[i], oracle[19 + i], 2e-5);
            }
            assert!(x[0] < initial()[0]);
            assert!(heat.iter().all(|q| *q > 0.0));
        }
    }
    assert!(errors[2] < 3e-5, "{errors:?}");
    for pair in errors.windows(2) {
        assert!(
            (3.2..4.8).contains(&(pair[0] / pair[1])),
            "second order: {errors:?}"
        );
    }
}

#[test]
fn g1_constant_regeneration_increases_finite_chemical_energy_with_exact_shaft_and_heat_work() {
    // Independently solved steady charging state: I_cell=I_converter=I_motor=-2 A.
    // Both shaft reactions are explicitly authored; neither inertia is pinned.
    let current = -2.0;
    let pol = BATTERY.polarization_resistance_ohm * current;
    let vb = 12.0 - pol - BATTERY.series_resistance_ohm * current;
    let vi = vb - LEAD * current;
    let converter_r = 0.07 + 0.02 + 0.02;
    let vo = vi - converter_r * current;
    let relative = (vo - MOTOR.winding_resistance_ohm * current) / MOTOR.coupling_nm_per_a;
    let reaction = -MOTOR.coupling_nm_per_a * current + MOTOR.bearing_drag_nm_s_per_rad * relative;
    let x = [
        0.5 * CAPACITY,
        pol * BATTERY.polarization_capacitance_f,
        vb * BATTERY.bus_capacitance_f,
        current * CONVERTER.inductance_h,
        vi * CONVERTER.bus_capacitance_f[0],
        vo * CONVERTER.bus_capacitance_f[1],
        current * MOTOR.inductance_h,
        relative * MOTOR.rotor_inertia_kg_m2,
        0.0,
    ];
    let intervals = [SwitchInterval {
        state: HIGH,
        duration_s: 0.01,
    }; 3];
    let model = model_with_curve(
        BatteryOcvCurve::new(CAPACITY, &[(0.0, 12.0), (1.0, 12.0)]).unwrap(),
        MODULATION,
    );
    let mut cycle = PreparedBatteryMotorSwitchCycle::new_controlled(
        model,
        &intervals,
        4,
        electrical(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(cycle.cached_interval_count(), 1);
    let dt = cycle.period_s();
    let mut next = [0.0; 9];
    let step = cycle
        .step_controlled(
            &x,
            BatteryMotorSwitchCommand {
                shaft_torque_nm: [reaction, -reaction],
                auxiliary_load: None,
            },
            &mut next,
            || false,
        )
        .unwrap();
    close(next[0], x[0] + 2.0 * dt, 1e-11);
    for i in 1..9 {
        close(next[i], x[i], 1e-11);
    }
    close(step.storage_change_j[0], 24.0 * dt, 2e-10);
    close(step.motor_electrical_work_j, vo * current * dt, 2e-10);
    close(step.battery_terminal_work_j, -vb * current * dt, 2e-10);
    close(step.shaft_work_j[0], reaction * relative * dt, 2e-10);
    close(step.shaft_work_j[1], 0.0, 2e-10);
    let expected_heat = [
        0.07 * 4.0,
        0.02 * 4.0,
        0.0,
        0.02 * 4.0,
        0.0,
        0.1 * 4.0,
        pol * pol / 0.2,
        LEAD * 4.0,
        0.2 * 4.0,
        0.01 * relative * relative,
    ];
    for (q, p) in step.heat_joules().into_iter().zip(expected_heat) {
        close(q, p * dt, 2e-10);
    }
    assert!(step.motor_electrical_work_j < 0.0 && step.storage_change_j[0] > 0.0);
    close(step.balance_residual_j(), 0.0, 2e-9);
}

#[test]
fn authored_commands_keep_their_edges_even_when_native_modes_share_one_cache() {
    let interval = SwitchInterval {
        state: HIGH,
        duration_s: 0.001,
    };
    let commands = [
        BatteryMotorSwitchCommand {
            shaft_torque_nm: [-0.2, 0.03],
            auxiliary_load: Some(RegulatedPowerLoad::try_new(2.0, 10.0, 0.2, 100.0).unwrap()),
        },
        BatteryMotorSwitchCommand {
            shaft_torque_nm: [0.6, -0.1],
            auxiliary_load: Some(RegulatedPowerLoad::try_new(8.0, 10.0, 0.2, 100.0).unwrap()),
        },
    ];
    let mut cycle = prepare(&[interval; 2], 1);
    assert_eq!(cycle.cached_interval_count(), 1);
    let mut native = model(BridgeModulation {
        input_high_fraction: 1.0,
        output_high_fraction: 1.0,
    })
    .prepare_controlled(interval.duration_s, electrical(), || false)
    .unwrap();
    let mut x = initial();
    let mut shaft_work = [0.0; 2];
    for command in commands {
        let mut next = [0.0; 9];
        let step = native
            .step_controlled(
                &x,
                &command.shaft_torque_nm,
                command.auxiliary_load,
                &mut next,
                || false,
            )
            .unwrap();
        for i in 0..2 {
            shaft_work[i] += step.shaft_work_j[i];
        }
        x = next;
    }
    let mut next = [0.0; 9];
    let step = cycle
        .step_interval_commands_controlled(&initial(), &commands, &mut next, || false)
        .unwrap();
    assert_eq!(next, x);
    close(step.load_work_j, 0.001 * (2.0 + 8.0), 2e-12);
    assert_eq!(step.last_auxiliary_load.mode, PowerLoadMode::ConstantPower);
    assert_eq!(step.last_auxiliary_load.requested_power_w, 8.0);
    for i in 0..2 {
        close(step.shaft_work_j[i], shaft_work[i], 2e-12);
    }
    close(
        step.native.supplied,
        shaft_work.iter().sum::<f64>() - step.load_work_j,
        2e-10,
    );
}

#[test]
fn late_source_refusal_and_cancellation_leave_all_nine_states_unchanged_then_retry() {
    let interval = SwitchInterval {
        state: HIGH,
        duration_s: 0.001,
    };
    let mut cycle = prepare(&[interval; 2], 4);
    let unsupported = BatteryMotorSwitchCommand {
        auxiliary_load: Some(RegulatedPowerLoad::try_new(1.0, 1.0, 0.5, 1.0).unwrap()),
        ..COMMAND
    };
    let mut next = [91.0; 9];
    assert!(
        cycle
            .step_interval_commands_controlled(
                &initial(),
                &[COMMAND, unsupported],
                &mut next,
                || false,
            )
            .is_err()
    );
    assert_eq!(next, [91.0; 9]);
    let mut count = 0;
    let first = cycle
        .step_controlled(&initial(), COMMAND, &mut next, || {
            count += 1;
            false
        })
        .unwrap();
    let expected = next;
    next = [91.0; 9];
    let mut polls = 0;
    assert!(matches!(
        cycle.step_controlled(&initial(), COMMAND, &mut next, || {
            polls += 1;
            polls >= count / 2
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [91.0; 9]);
    let retry = cycle
        .step_controlled(&initial(), COMMAND, &mut next, || false)
        .unwrap();
    assert_eq!(next, expected);
    assert_eq!(retry.native, first.native);
    assert_eq!(retry.heat_joules(), first.heat_joules());
    let mut depleted = initial();
    depleted[0] = 0.0;
    depleted[2] = 9.0 * BATTERY.bus_capacitance_f;
    next = [91.0; 9];
    assert!(
        cycle
            .step_controlled(&depleted, COMMAND, &mut next, || false)
            .is_err()
    );
    assert_eq!(next, [91.0; 9]);
}

#[test]
fn physical_edges_and_numerical_work_are_bounded_before_preparation() {
    let interval = SwitchInterval {
        state: HIGH,
        duration_s: 0.001,
    };
    for subdivisions in [0, 257] {
        assert!(
            PreparedBatteryMotorSwitchCycle::new_controlled(
                model(MODULATION),
                &[interval],
                subdivisions,
                electrical(),
                cycle_controls(),
                || false,
            )
            .is_err()
        );
    }
    for intervals in [
        Vec::new(),
        vec![interval; 65],
        vec![SwitchInterval {
            duration_s: 0.0,
            ..interval
        }],
    ] {
        assert!(
            PreparedBatteryMotorSwitchCycle::new_controlled(
                model(MODULATION),
                &intervals,
                1,
                electrical(),
                cycle_controls(),
                || false,
            )
            .is_err()
        );
    }
    let mut cycle = prepare(&[interval], 1);
    let mut next = [91.0; 9];
    assert!(
        cycle
            .step_interval_commands_controlled(&initial(), &[], &mut next, || false)
            .is_err()
    );
    assert_eq!(next, [91.0; 9]);
}
