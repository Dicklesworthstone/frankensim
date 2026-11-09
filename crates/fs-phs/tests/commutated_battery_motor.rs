//! G1/G3/G4: finite chemical energy, physical commutation and two moving shafts.
//! The reference solves scalar device KCL by bisection and integrates SI laws;
//! it does not call the production half-bridge, port, step or loss owners.

use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::battery_motor::{
    CommutatedBatteryMotorCycleStep, PreparedCommutatedBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationInterval, DeviceGateState, DeviceOperatingBounds, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PreparedStepError, SolveBudget};
use std::cell::Cell;

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
    // Deliberately different from DEVICE: the new owner must replace these
    // averaged channel laws, while retaining the actual storage and winding.
    input_on_resistance_ohm: [0.031, 0.047],
    output_on_resistance_ohm: [0.053, 0.071],
};
const MOTOR: DcMotorParameters = DcMotorParameters {
    inductance_h: 0.02,
    winding_resistance_ohm: 0.2,
    rotor_inertia_kg_m2: 0.08,
    stator_inertia_kg_m2: 0.4,
    coupling_nm_per_a: 0.4,
    bearing_drag_nm_s_per_rad: 0.01,
};
const DEVICE: SwitchDevice = SwitchDevice {
    on_resistance_ohm: 0.08,
    off_resistance_ohm: 500.0,
    body_diode: BodyDiode {
        forward_voltage_v: 0.45,
        differential_resistance_ohm: 0.12,
    },
    bounds: DeviceOperatingBounds {
        maximum_reverse_voltage_v: 100.0,
        maximum_channel_current_a: 200.0,
        maximum_diode_current_a: 100.0,
    },
};
const HIGH: [bool; 4] = [true, false, true, false];

fn model(constant_ocv: bool) -> BatteryMotorDrive {
    let knots = if constant_ocv {
        [(0.0, 12.0), (1.0, 12.0)]
    } else {
        [(0.0, 10.0), (1.0, 14.0)]
    };
    AveragedBuckBoost::new(
        CONVERTER,
        BridgeModulation {
            input_high_fraction: 0.2,
            output_high_fraction: 0.8,
        },
    )
    .unwrap()
    .with_battery_source(
        BatterySource::new(BatteryOcvCurve::new(CAPACITY, &knots).unwrap(), BATTERY).unwrap(),
        LEAD,
    )
    .unwrap()
    .with_dc_motor(DcMotor::new(MOTOR).unwrap())
    .unwrap()
}
fn initial() -> [f64; 9] {
    [
        0.65 * CAPACITY,
        0.1 * BATTERY.polarization_capacitance_f,
        12.1 * BATTERY.bus_capacitance_f,
        0.5 * CONVERTER.inductance_h,
        11.8 * CONVERTER.bus_capacitance_f[0],
        6.4 * CONVERTER.bus_capacitance_f[1],
        0.1 * MOTOR.inductance_h,
        4.0 * MOTOR.rotor_inertia_kg_m2,
        -0.4 * MOTOR.stator_inertia_kg_m2,
    ]
}
fn electrical() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 10,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-11,
    }
}
fn cycle_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        absolute_energy_tolerance_j: 4e-9,
        relative_energy_tolerance: 2e-11,
        ..SwitchCycleControls::default()
    }
}
fn interval(on: [bool; 4], duration_s: f64) -> CommutationInterval {
    CommutationInterval {
        gates: DeviceGateState { on },
        duration_s,
    }
}
fn prepare(
    intervals: &[CommutationInterval],
    subdivisions: usize,
) -> PreparedCommutatedBatteryMotorCycle {
    PreparedCommutatedBatteryMotorCycle::new_controlled(
        model(false),
        [DEVICE; 4],
        intervals,
        subdivisions,
        electrical(),
        cycle_controls(),
        || false,
    )
    .unwrap()
}
fn command(power: f64) -> BatteryMotorSwitchCommand {
    BatteryMotorSwitchCommand {
        shaft_torque_nm: [-0.2, 0.03],
        auxiliary_load: (power > 0.0)
            .then(|| RegulatedPowerLoad::try_new(power, 3.0, 0.5, 50.0).unwrap()),
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

struct Leg {
    node: f64,
    vds: [f64; 2],
    channel: [f64; 2],
    diode: [f64; 2],
}
impl Leg {
    fn high_current(&self) -> f64 {
        self.channel[0] - self.diode[0]
    }
}

// Independent monotone scalar KCL, not the native affine active-set elimination.
fn leg(bus: f64, draw: f64, on: [bool; 2]) -> Leg {
    let resistance = on.map(|gate| {
        if gate {
            DEVICE.on_resistance_ohm
        } else {
            DEVICE.off_resistance_ohm
        }
    });
    let current = |voltage: f64, index: usize| {
        voltage / resistance[index]
            - ((-voltage - DEVICE.body_diode.forward_voltage_v)
                / DEVICE.body_diode.differential_resistance_ohm)
                .max(0.0)
    };
    let residual = |node: f64| current(bus - node, 0) - current(node, 1) - draw;
    let span = bus.abs() + (draw.abs() + 1.0) * DEVICE.off_resistance_ohm + 1.0;
    let (mut low, mut high) = (-span, span);
    assert!(residual(low) > 0.0 && residual(high) < 0.0);
    for _ in 0..80 {
        let middle = f64::midpoint(low, high);
        if residual(middle) > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let node = f64::midpoint(low, high);
    let vds = [bus - node, node];
    let channel = core::array::from_fn(|i| vds[i] / resistance[i]);
    let diode = vds.map(|v| {
        ((-v - DEVICE.body_diode.forward_voltage_v) / DEVICE.body_diode.differential_resistance_ohm)
            .max(0.0)
    });
    Leg {
        node,
        vds,
        channel,
        diode,
    }
}

fn load_current(load: Option<RegulatedPowerLoad>, voltage: f64) -> f64 {
    let Some(load) = load else {
        return 0.0;
    };
    assert!(voltage.abs() <= load.maximum_voltage_v());
    let magnitude = if voltage.abs() <= load.foldback_voltage_v() {
        load.requested_power_w() * voltage.abs() / load.foldback_voltage_v().powi(2)
    } else {
        load.requested_power_w() / voltage.abs()
    };
    magnitude.min(load.maximum_current_a()).copysign(voltage)
}

// Nine original storage coordinates, fourteen independent heat integrals and
// eight signed work integrals. All powers are explicit SI products.
fn derivative(x: &[f64; 31], on: [bool; 4], cmd: BatteryMotorSwitchCommand) -> [f64; 31] {
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
    let input = leg(vi, ic, [on[0], on[1]]);
    let output = leg(vo, -ic, [on[2], on[3]]);
    let load = load_current(cmd.auxiliary_load, vo);
    let torque = MOTOR.coupling_nm_per_a * im - MOTOR.bearing_drag_nm_s_per_rad * relative;
    let mut out = [0.0; 31];
    out[..9].copy_from_slice(&[
        -cell,
        cell - pol / BATTERY.polarization_resistance_ohm,
        cell - lead,
        input.node - output.node - CONVERTER.winding_resistance_ohm * ic,
        lead - input.high_current(),
        -output.high_current() - im - load,
        vo - MOTOR.winding_resistance_ohm * im - MOTOR.coupling_nm_per_a * relative,
        torque + cmd.shaft_torque_nm[0],
        -torque + cmd.shaft_torque_nm[1],
    ]);
    out[9] = CONVERTER.winding_resistance_ohm * ic * ic;
    for (p, leg) in [input, output].iter().enumerate() {
        for j in 0..2 {
            out[10 + 2 * p + j] = leg.vds[j] * leg.channel[j];
            out[19 + 2 * p + j] = -leg.vds[j] * leg.diode[j];
        }
    }
    out[14] = BATTERY.series_resistance_ohm * cell * cell;
    out[15] = pol * pol / BATTERY.polarization_resistance_ohm;
    out[16] = LEAD * lead * lead;
    out[17] = MOTOR.winding_resistance_ohm * im * im;
    out[18] = MOTOR.bearing_drag_nm_s_per_rad * relative * relative;
    out[23..].copy_from_slice(&[
        -vb * lead,
        vi * lead,
        -vo * (im + load),
        vo * im,
        MOTOR.coupling_nm_per_a * relative * im,
        cmd.shaft_torque_nm[0] * wr,
        cmd.shaft_torque_nm[1] * ws,
        vo * load,
    ]);
    out
}
fn reference(
    intervals: &[CommutationInterval],
    cycles: usize,
    subdivisions: usize,
    cmd: BatteryMotorSwitchCommand,
) -> [f64; 31] {
    let mut x = [0.0; 31];
    x[..9].copy_from_slice(&initial());
    for _ in 0..cycles {
        for interval in intervals {
            let h = interval.duration_s / subdivisions as f64;
            for _ in 0..subdivisions {
                let a = derivative(&x, interval.gates.on, cmd);
                let b = derivative(
                    &core::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
                    interval.gates.on,
                    cmd,
                );
                let c = derivative(
                    &core::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
                    interval.gates.on,
                    cmd,
                );
                let d = derivative(
                    &core::array::from_fn(|i| x[i] + h * c[i]),
                    interval.gates.on,
                    cmd,
                );
                x = core::array::from_fn(|i| {
                    x[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
                });
            }
        }
    }
    x
}

fn storage(x: &[f64; 9], constant_ocv: bool) -> [f64; 7] {
    [
        if constant_ocv {
            12.0 * x[0]
        } else {
            10.0 * x[0] + 2.0 * x[0] * x[0] / CAPACITY
        },
        0.5 * x[1] * x[1] / BATTERY.polarization_capacitance_f,
        0.5 * x[2] * x[2] / BATTERY.bus_capacitance_f,
        0.5 * x[3] * x[3] / CONVERTER.inductance_h
            + 0.5 * x[4] * x[4] / CONVERTER.bus_capacitance_f[0]
            + 0.5 * x[5] * x[5] / CONVERTER.bus_capacitance_f[1],
        0.5 * x[6] * x[6] / MOTOR.inductance_h,
        0.5 * x[7] * x[7] / MOTOR.rotor_inertia_kg_m2,
        0.5 * x[8] * x[8] / MOTOR.stator_inertia_kg_m2,
    ]
}
fn works(step: &CommutatedBatteryMotorCycleStep) -> [f64; 8] {
    [
        step.battery_terminal_work_j,
        step.converter_terminal_work_j[0],
        step.converter_terminal_work_j[1],
        step.motor_electrical_work_j,
        step.electromagnetic_work_j,
        step.shaft_work_j[0],
        step.shaft_work_j[1],
        step.load_work_j,
    ]
}
fn account(
    before: &[f64; 9],
    after: &[f64; 9],
    step: &CommutatedBatteryMotorCycleStep,
    constant_ocv: bool,
) {
    let a = storage(before, constant_ocv);
    let b = storage(after, constant_ocv);
    let change: [f64; 7] = core::array::from_fn(|i| b[i] - a[i]);
    let q = step.heat_joules();
    let w = works(step);
    for i in 0..7 {
        close(step.storage_change_j[i], change[i], 3e-9);
    }
    assert!(q.iter().all(|value| value.is_finite() && *value >= 0.0));
    close(change[..3].iter().sum::<f64>() + q[5] + q[6], w[0], 3e-9);
    close(w[0] + q[7] + w[1], 0.0, 3e-9);
    close(
        change[3] + q[..5].iter().sum::<f64>() + q[10..].iter().sum::<f64>(),
        w[1] + w[2],
        3e-9,
    );
    close(w[2] + w[3] + w[7], 0.0, 3e-9);
    close(change[4] + q[8] + w[4], w[3], 3e-9);
    close(change[5] + change[6] + q[9], w[4] + w[5] + w[6], 3e-9);
    close(step.native.delta_h, change.iter().sum(), 3e-9);
    close(step.native.dissipated, q.iter().sum(), 3e-9);
    close(step.native.supplied, w[5] + w[6] - w[7], 3e-9);
    close(
        change.iter().sum::<f64>() + q.iter().sum::<f64>(),
        w[5] + w[6] - w[7],
        3e-9,
    );
    close(step.balance_residual_j(), 0.0, 3e-9);
    close(step.angular_momentum_residual_nms, 0.0, 3e-11);
}

#[test]
fn g1_nine_states_fourteen_heats_and_each_internal_work_refine_to_independent_kcl_ode() {
    let intervals = [
        interval(HIGH, 0.0015),
        interval([false; 4], 0.0002),
        interval([false, true, false, true], 0.0007),
        interval([false, true, true, false], 0.0004),
        interval([true; 4], 0.0001),
        interval(HIGH, 0.0011),
    ];
    let cmd = command(2.0);
    let oracle = reference(&intervals, 6, 512, cmd);
    let refined = reference(&intervals, 6, 256, cmd);
    for i in 0..31 {
        close(oracle[i], refined[i], 1e-9);
    }
    let mut errors = Vec::new();
    for subdivisions in [4, 8, 16] {
        let mut native = prepare(&intervals, subdivisions);
        let mut x = initial();
        let mut heat = [0.0; 14];
        let mut work = [0.0; 8];
        for _ in 0..6 {
            let mut next = [0.0; 9];
            let step = native
                .step_controlled(&x, cmd, &mut next, || false)
                .unwrap();
            account(&x, &next, &step, false);
            assert_eq!(step.integration_steps, intervals.len() * subdivisions);
            for (total, value) in heat.iter_mut().zip(step.heat_joules()) {
                *total += value;
            }
            for (total, value) in work.iter_mut().zip(works(&step)) {
                *total += value;
            }
            x = next;
        }
        let state_error = x
            .iter()
            .zip(&oracle[..9])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let power_error = heat
            .iter()
            .chain(&work)
            .zip(&oracle[9..])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push([state_error, power_error]);
        if subdivisions == 16 {
            assert!(state_error < 2e-5 && power_error < 8e-5, "{errors:?}");
            assert!(x[0] < initial()[0]);
            assert!(heat[..10].iter().all(|q| *q > 0.0));
            assert!(heat[11] > 0.0 && heat[12] > 0.0);
        }
    }
    for metric in 0..2 {
        for pair in errors.windows(2) {
            let order = pair[0][metric] / pair[1][metric];
            assert!(
                (3.2..4.8).contains(&order),
                "second-order state/work refinement: {errors:?}"
            );
        }
    }
}

#[test]
fn g1_regeneration_increases_chemical_charge_with_finite_leakage_and_two_authored_shafts() {
    // Constant OCV admits an exact charging trajectory: only chemical charge
    // changes. Leakage makes I_cell, I_converter and I_motor DIFFERENT.
    let cell = -2.0;
    let pol = BATTERY.polarization_resistance_ohm * cell;
    let vb = 12.0 - pol - BATTERY.series_resistance_ohm * cell;
    let vi = vb - LEAD * cell;
    let left = vi - DEVICE.on_resistance_ohm * cell;
    let ic = cell - left / DEVICE.off_resistance_ohm;
    let right = left - CONVERTER.winding_resistance_ohm * ic;
    let output_high = -ic + right / DEVICE.off_resistance_ohm;
    let vo = right + DEVICE.on_resistance_ohm * output_high;
    let im = -output_high;
    let relative = (vo - MOTOR.winding_resistance_ohm * im) / MOTOR.coupling_nm_per_a;
    let reaction = -MOTOR.coupling_nm_per_a * im + MOTOR.bearing_drag_nm_s_per_rad * relative;
    let ws = -2.0;
    let wr = ws + relative;
    let before = [
        CAPACITY / 2.0,
        pol * BATTERY.polarization_capacitance_f,
        vb * BATTERY.bus_capacitance_f,
        ic * CONVERTER.inductance_h,
        vi * CONVERTER.bus_capacitance_f[0],
        vo * CONVERTER.bus_capacitance_f[1],
        im * MOTOR.inductance_h,
        wr * MOTOR.rotor_inertia_kg_m2,
        ws * MOTOR.stator_inertia_kg_m2,
    ];
    let intervals = [interval(HIGH, 1.0 / 256.0); 3];
    let mut native = PreparedCommutatedBatteryMotorCycle::new_controlled(
        model(true),
        [DEVICE; 4],
        &intervals,
        4,
        electrical(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    let dt = native.period_s();
    let mut next = [0.0; 9];
    let step = native
        .step_controlled(
            &before,
            BatteryMotorSwitchCommand {
                shaft_torque_nm: [reaction, -reaction],
                auxiliary_load: None,
            },
            &mut next,
            || false,
        )
        .unwrap();
    close(next[0], before[0] - cell * dt, 2e-11);
    for i in 1..9 {
        close(next[i], before[i], 2e-11);
    }
    account(&before, &next, &step, true);
    close(step.storage_change_j[0], -12.0 * cell * dt, 2e-10);
    close(step.battery_terminal_work_j, -vb * cell * dt, 2e-10);
    close(step.converter_terminal_work_j[0], vi * cell * dt, 2e-10);
    close(step.motor_electrical_work_j, vo * im * dt, 2e-10);
    close(step.shaft_work_j[0], reaction * wr * dt, 2e-10);
    close(step.shaft_work_j[1], -reaction * ws * dt, 2e-10);
    let expected = [
        CONVERTER.winding_resistance_ohm * ic * ic,
        DEVICE.on_resistance_ohm * cell * cell,
        left * left / DEVICE.off_resistance_ohm,
        DEVICE.on_resistance_ohm * output_high * output_high,
        right * right / DEVICE.off_resistance_ohm,
        BATTERY.series_resistance_ohm * cell * cell,
        pol * pol / BATTERY.polarization_resistance_ohm,
        LEAD * cell * cell,
        MOTOR.winding_resistance_ohm * im * im,
        MOTOR.bearing_drag_nm_s_per_rad * relative * relative,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    for (actual, power) in step.heat_joules().into_iter().zip(expected) {
        close(actual, dt * power, 2e-10);
    }
    assert!(step.motor_electrical_work_j < 0.0 && step.storage_change_j[0] > 0.0);
    assert!(step.shaft_work_j.iter().all(|w| *w > 0.0));
}

#[test]
fn dead_time_both_current_signs_follow_each_physical_diode_at_the_same_discrete_effort() {
    let dt = 1.0 / 65536.0;
    for current in [2.0, -2.0] {
        let mut before = initial();
        before[3] = current * CONVERTER.inductance_h;
        let mut native = prepare(&[interval([false; 4], dt)], 1);
        let mut next = [0.0; 9];
        let cmd = command(0.0);
        let step = native
            .step_controlled(&before, cmd, &mut next, || false)
            .unwrap();
        account(&before, &next, &step, false);
        let mut midpoint = [0.0; 31];
        for i in 0..9 {
            midpoint[i] = f64::midpoint(before[i], next[i]);
        }
        let reference = derivative(&midpoint, [false; 4], cmd);
        for (q, power) in step.heat_joules().iter().zip(&reference[9..23]) {
            close(*q, dt * power, 2e-10);
        }
        let input = leg(
            midpoint[4] / CONVERTER.bus_capacitance_f[0],
            midpoint[3] / CONVERTER.inductance_h,
            [false; 2],
        );
        let output = leg(
            midpoint[5] / CONVERTER.bus_capacitance_f[1],
            -midpoint[3] / CONVERTER.inductance_h,
            [false; 2],
        );
        for (p, expected) in [input, output].iter().enumerate() {
            for j in 0..2 {
                let actual = step.last_devices[2 * p + j];
                close(actual.drain_source_voltage_v, expected.vds[j], 2e-10);
                close(actual.channel_current_a, expected.channel[j], 2e-10);
                close(actual.diode_current_a, expected.diode[j], 2e-10);
            }
        }
        let active = if current > 0.0 { [1, 2] } else { [0, 3] };
        for i in 0..4 {
            if active.contains(&i) {
                assert!(step.heat_joules()[10 + i] > 0.0);
            } else {
                assert_eq!(step.heat_joules()[10 + i], 0.0);
            }
        }
    }
}

#[test]
fn physical_gate_overlap_draws_real_finite_bus_energy_without_inventing_a_source() {
    let dt = 1.0 / 4096.0;
    let before = initial();
    let mut native = prepare(&[interval([true; 4], dt)], 1);
    let mut next = [0.0; 9];
    let cmd = BatteryMotorSwitchCommand {
        shaft_torque_nm: [0.0; 2],
        auxiliary_load: None,
    };
    let step = native
        .step_controlled(&before, cmd, &mut next, || false)
        .unwrap();
    account(&before, &next, &step, false);
    let mut midpoint = [0.0; 31];
    for i in 0..9 {
        midpoint[i] = f64::midpoint(before[i], next[i]);
    }
    let reference = derivative(&midpoint, [true; 4], cmd);
    for (q, power) in step.heat_joules().iter().zip(&reference[9..23]) {
        close(*q, dt * power, 2e-10);
    }
    assert!(next[4] < before[4] && next[5] < before[5]);
    assert!(step.heat_joules()[1..5].iter().all(|q| *q > 0.01));
    assert_eq!(step.heat_joules()[10..], [0.0; 4]);
    close(step.native.supplied, 0.0, 1e-12);
    close(step.native.delta_h, -step.native.dissipated, 3e-9);
}

#[test]
fn authored_commands_continue_at_their_gate_edges_including_actual_load_curtailment() {
    let dt = 1.0 / 1024.0;
    let deck = [interval(HIGH, dt); 2];
    let commands = [
        command(2.0),
        BatteryMotorSwitchCommand {
            shaft_torque_nm: [0.7, -0.1],
            auxiliary_load: Some(RegulatedPowerLoad::try_new(10.0, 1.0, 0.5, 50.0).unwrap()),
        },
    ];
    let mut together = prepare(&deck, 3);
    let mut split = prepare(&deck[..1], 3);
    let mut complete = [0.0; 9];
    let before = initial();
    let receipt = together
        .step_interval_commands_controlled(&before, &commands, &mut complete, || false)
        .unwrap();
    let mut x = before;
    let mut heat = [0.0; 14];
    let mut work = [0.0; 8];
    for command in commands {
        let mut next = [0.0; 9];
        let step = split
            .step_controlled(&x, command, &mut next, || false)
            .unwrap();
        account(&x, &next, &step, false);
        for (q, value) in heat.iter_mut().zip(step.heat_joules()) {
            *q += value;
        }
        for (w, value) in work.iter_mut().zip(works(&step)) {
            *w += value;
        }
        x = next;
    }
    assert_eq!(complete.map(f64::to_bits), x.map(f64::to_bits));
    for (a, b) in receipt.heat_joules().into_iter().zip(heat) {
        close(a, b, 2e-12);
    }
    for (a, b) in works(&receipt).into_iter().zip(work) {
        close(a, b, 2e-12);
    }
    close(
        receipt.load_work_j + receipt.curtailed_load_energy_j,
        12.0 * dt,
        2e-12,
    );
    assert!(receipt.curtailed_load_energy_j > 0.0);
    assert_eq!(receipt.integration_steps, 6);
    account(&before, &complete, &receipt, false);
}

#[test]
fn late_device_load_soc_and_cancellation_refusals_publish_none_of_the_nine_states() {
    let dt = 1.0 / 1024.0;
    let before = initial();
    let sentinel = [77.0; 9];
    let deck = [interval(HIGH, dt); 2];
    let cmd = command(0.0);
    let mut native = prepare(&deck, 3);
    let mut expected = [0.0; 9];
    let polls = Cell::new(0usize);
    let good = native
        .step_controlled(&before, cmd, &mut expected, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let stop = polls.get() - 1;
    polls.set(0);
    let mut next = sentinel;
    assert!(matches!(
        native.step_controlled(&before, cmd, &mut next, || {
            polls.set(polls.get() + 1);
            polls.get() >= stop
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, sentinel);
    let invalid_load = BatteryMotorSwitchCommand {
        shaft_torque_nm: [0.0; 2],
        auxiliary_load: Some(RegulatedPowerLoad::try_new(1.0, 1.0, 0.25, 0.5).unwrap()),
    };
    assert!(
        native
            .step_interval_commands_controlled(&before, &[cmd, invalid_load], &mut next, || false)
            .is_err()
    );
    assert_eq!(next, sentinel);
    let mut empty = [0.0; 9];
    empty[0] = 1e-8;
    assert!(
        native
            .step_controlled(&empty, cmd, &mut next, || false)
            .is_err()
    );
    assert_eq!(next, sentinel);
    let retry = native
        .step_controlled(&before, cmd, &mut next, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(retry.heat_joules(), good.heat_joules());
    let mut limited = [DEVICE; 4];
    for d in &mut limited {
        d.bounds.maximum_channel_current_a = 10.0;
    }
    let mut hidden = PreparedCommutatedBatteryMotorCycle::new_controlled(
        model(false),
        limited,
        &[interval(HIGH, dt), interval([true; 4], dt)],
        1,
        electrical(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    next = sentinel;
    assert!(
        hidden
            .step_controlled(&before, cmd, &mut next, || false)
            .is_err()
    );
    assert_eq!(next, sentinel);
    for subdivisions in [0, 257] {
        assert!(
            PreparedCommutatedBatteryMotorCycle::new_controlled(
                model(false),
                [DEVICE; 4],
                &deck,
                subdivisions,
                electrical(),
                cycle_controls(),
                || false
            )
            .is_err()
        );
    }
    assert_eq!(before, initial());
}
