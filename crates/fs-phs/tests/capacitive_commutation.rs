//! G1/G3/G4: physical switch-node charge, exact finite-capacitor RC laws,
//! independent branch work, and transactional gate-edge continuation.
//! All capacitances below are explicitly constant synthetic fixture values.

use fs_phs::circuit::converter::commutation::capacitance::{
    CapacitiveCommutatedStep, CapacitiveCommutationParameters,
    PreparedCapacitiveCommutatedSwitchCycle,
};
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationInterval, CommutationParameters, DeviceGateState, DeviceOperatingBounds,
    HalfBridgeConduction, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::power_load::PowerSwitchCommand;
use fs_phs::{PreparedStepError, SolveBudget};
use std::cell::Cell;

const DEVICE: SwitchDevice = SwitchDevice {
    on_resistance_ohm: 0.2,
    off_resistance_ohm: 40.0,
    body_diode: BodyDiode {
        forward_voltage_v: 0.5,
        differential_resistance_ohm: 0.1,
    },
    bounds: DeviceOperatingBounds {
        maximum_reverse_voltage_v: 100.0,
        maximum_channel_current_a: 100.0,
        maximum_diode_current_a: 100.0,
    },
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
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
fn controls() -> SwitchCycleControls {
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
fn parameters() -> CapacitiveCommutationParameters {
    CapacitiveCommutationParameters {
        converter: CommutationParameters {
            inductance_h: 0.4,
            bus_capacitance_f: [0.3, 0.4],
            winding_resistance_ohm: 0.1,
            devices: [DEVICE; 4],
        },
        output_capacitance_f: [0.2, 0.1, 0.3, 0.2],
    }
}
fn prepare(
    p: CapacitiveCommutationParameters,
    resistance: [f64; 2],
    deck: &[CommutationInterval],
    subdivisions: usize,
) -> PreparedCapacitiveCommutatedSwitchCycle {
    PreparedCapacitiveCommutatedSwitchCycle::new_controlled(
        p,
        resistance,
        deck,
        subdivisions,
        electrical(),
        controls(),
        || false,
    )
    .unwrap()
}

// Invert the explicitly authored two-electrode capacitance matrix ourselves.
// q_bus is TOTAL bus-electrode charge; it is not C_bus * V_bus.
fn physical(p: CapacitiveCommutationParameters, x: &[f64; 5]) -> [f64; 5] {
    let mut y = [0.0; 5];
    y[0] = x[0] / p.converter.inductance_h;
    for leg in 0..2 {
        let cb = p.converter.bus_capacitance_f[leg];
        let ch = p.output_capacitance_f[2 * leg];
        let cl = p.output_capacitance_f[2 * leg + 1];
        let det = cb * (ch + cl) + ch * cl;
        y[1 + leg] = ((ch + cl) * x[1 + leg] + ch * x[3 + leg]) / det;
        y[3 + leg] = (ch * x[1 + leg] + (cb + ch) * x[3 + leg]) / det;
    }
    y
}
fn energies(p: CapacitiveCommutationParameters, x: &[f64; 5]) -> [f64; 7] {
    let y = physical(p, x);
    [
        0.5 * p.converter.inductance_h * y[0].powi(2),
        0.5 * p.converter.bus_capacitance_f[0] * y[1].powi(2),
        0.5 * p.converter.bus_capacitance_f[1] * y[2].powi(2),
        0.5 * p.output_capacitance_f[0] * (y[1] - y[3]).powi(2),
        0.5 * p.output_capacitance_f[1] * y[3].powi(2),
        0.5 * p.output_capacitance_f[2] * (y[2] - y[4]).powi(2),
        0.5 * p.output_capacitance_f[3] * y[4].powi(2),
    ]
}
fn heat(step: &CapacitiveCommutatedStep) -> [f64; 11] {
    let c = &step.converter;
    [
        c.winding_heat_j,
        c.channel_heat_j[0],
        c.channel_heat_j[1],
        c.channel_heat_j[2],
        c.channel_heat_j[3],
        c.diode_heat_j[0],
        c.diode_heat_j[1],
        c.diode_heat_j[2],
        c.diode_heat_j[3],
        step.terminal_resistor_heat_j[0],
        step.terminal_resistor_heat_j[1],
    ]
}
fn account(
    p: CapacitiveCommutationParameters,
    before: &[f64; 5],
    after: &[f64; 5],
    step: &CapacitiveCommutatedStep,
) {
    let a = energies(p, before);
    let b = energies(p, after);
    let change: [f64; 7] = core::array::from_fn(|i| b[i] - a[i]);
    let c = &step.converter;
    close(c.inductor_energy_change_j, change[0], 2e-9);
    for i in 0..2 {
        close(c.bus_capacitor_energy_change_j[i], change[1 + i], 2e-9);
        close(
            step.source_work_j[i] - step.terminal_resistor_heat_j[i] - step.load_work_j[i],
            c.terminal_work_j[i],
            2e-9,
        );
    }
    for i in 0..4 {
        close(c.device_capacitor_energy_change_j[i], change[3 + i], 2e-9);
    }
    let q = heat(step);
    assert!(q.iter().all(|v| v.is_finite() && *v >= 0.0));
    close(c.storage_change_j, change.iter().sum(), 2e-9);
    close(step.native.delta_h, change.iter().sum(), 2e-9);
    close(step.native.dissipated, q[..9].iter().sum(), 2e-9);
    close(step.native.supplied, c.terminal_work_j.iter().sum(), 2e-9);
    close(
        change.iter().sum::<f64>() + q.iter().sum::<f64>() + step.load_work_j.iter().sum::<f64>(),
        step.source_work_j.iter().sum(),
        2e-9,
    );
}

#[test]
fn stored_node_conduction_keeps_actual_charge_injection_and_exact_active_slopes() {
    let leg = HalfBridgeConduction::new(DEVICE, DEVICE, [false, true]).unwrap();
    let point = leg.evaluate_at_node(6.0, -1.0).unwrap();
    close(point.switch_node_voltage_v, -1.0, 0.0);
    close(point.high.terminal_current_a, 7.0 / 40.0, 1e-14);
    close(point.low.channel_current_a, -5.0, 1e-14);
    close(point.low.diode_current_a, 5.0, 1e-14);
    close(
        point.high.terminal_current_a - point.low.terminal_current_a,
        10.175,
        1e-13,
    );
    close(
        point.d_high_current_d_drain_source_voltage_s,
        1.0 / 40.0,
        0.0,
    );
    close(point.d_low_current_d_drain_source_voltage_s, 15.0, 0.0);
    let eps = 1e-5;
    let high_slope = (leg
        .evaluate_at_node(6.0 + eps, -1.0)
        .unwrap()
        .high
        .terminal_current_a
        - leg
            .evaluate_at_node(6.0 - eps, -1.0)
            .unwrap()
            .high
            .terminal_current_a)
        / (2.0 * eps);
    let low_slope = (leg
        .evaluate_at_node(6.0, -1.0 + eps)
        .unwrap()
        .low
        .terminal_current_a
        - leg
            .evaluate_at_node(6.0, -1.0 - eps)
            .unwrap()
            .low
            .terminal_current_a)
        / (2.0 * eps);
    close(
        high_slope,
        point.d_high_current_d_drain_source_voltage_s,
        1e-10,
    );
    close(
        low_slope,
        point.d_low_current_d_drain_source_voltage_s,
        1e-9,
    );
    close(point.low.channel_power_w, 5.0, 1e-14);
    close(point.low.diode_power_w, 5.0, 1e-14);
    assert!(leg.evaluate_at_node(101.0, 0.0).is_err());
    assert!(leg.evaluate_at_node(6.0, -11.0).is_err());
    assert!(leg.evaluate_at_node(f64::NAN, 0.0).is_err());
}

#[test]
fn finite_capacitor_rc_mode_has_exact_midpoint_state_and_seven_separate_storage_parts() {
    let c = 0.25;
    let r = 2.0;
    let mut p = parameters();
    p.converter.bus_capacitance_f = [c; 2];
    p.output_capacitance_f = [c; 4];
    for d in &mut p.converter.devices {
        d.on_resistance_ohm = r;
        d.off_resistance_ohm = r;
    }
    let dt = 0.125;
    let before = p.state(0.0, [2.0; 2], [1.0; 2]).unwrap();
    // Bus total charge includes the high-device capacitor; each node also
    // includes that capacitor with the opposite electrode sign.
    assert_eq!(before, [0.0, 3.0 * c, 3.0 * c, 0.0, 0.0]);
    let mut native = prepare(p, [r; 2], &[interval([true, false, true, false], dt)], 1);
    let (mut next, mut mean) = ([0.0; 5], [0.0; 2]);
    let step = native
        .step_controlled(&before, &[0.0; 2], &[None; 2], &mut next, &mut mean, || {
            false
        })
        .unwrap();
    let ratio = (1.0 - dt / (2.0 * r * c)) / (1.0 + dt / (2.0 * r * c));
    for (actual, expected) in next.into_iter().zip(before.map(|x| ratio * x)) {
        close(actual, expected, 2e-12);
    }
    let released = 1.0 - ratio * ratio;
    for port in 0..2 {
        close(mean[port], 1.0 + ratio, 2e-12);
        close(
            step.terminal_resistor_heat_j[port],
            2.0 * c * released,
            2e-12,
        );
        close(
            step.converter.bus_capacitor_energy_change_j[port],
            -2.0 * c * released,
            2e-12,
        );
        assert_eq!(step.source_work_j[port], 0.0);
    }
    for device in 0..4 {
        close(
            step.converter.channel_heat_j[device],
            0.5 * c * released,
            2e-12,
        );
        close(
            step.converter.device_capacitor_energy_change_j[device],
            -0.5 * c * released,
            2e-12,
        );
        assert_eq!(step.converter.diode_heat_j[device], 0.0);
    }
    close(step.converter.winding_heat_j, 0.0, 1e-25);
    account(p, &before, &next, &step);
}

#[test]
fn finite_charge_redistribution_and_each_resistor_heat_refine_to_two_exponential_solution() {
    charge_redistribution(1.0);
    // Exact SI similarity: 0.2/0.3 nF output capacitors over 0.5 ns. The
    // physical voltage error and scale-normalized heat error must still refine;
    // a tiny-joule defect cannot hide beneath a fixed absolute energy allowance.
    charge_redistribution(1e-9);
}

fn charge_redistribution(scale: f64) {
    // Two identical legs keep the inductor at zero. Each has the exact physical
    // M=[[Cb+Ch,-Ch],[-Ch,Ch+Cl]] and G=[[Gs+Gh,-Gh],[-Gh,Gh+Gl]].
    // Unlike an algebraic switch, the initially empty node charges continuously.
    let mut p = parameters();
    let (cb, ch, cl, gh, gl, gs) = (0.5 * scale, 0.2 * scale, 0.3 * scale, 1.0, 0.01, 0.01);
    p.converter.bus_capacitance_f = [cb; 2];
    p.output_capacitance_f = [ch, cl, ch, cl];
    for d in &mut p.converter.devices {
        d.on_resistance_ohm = 1.0 / gh;
        d.off_resistance_ohm = 1.0 / gl;
    }
    let det = cb * (ch + cl) + ch * cl;
    let a: [[f64; 2]; 2] = [
        [
            ((ch + cl) * (gs + gh) - ch * gh) / det,
            (-(ch + cl) * gh + ch * (gh + gl)) / det,
        ],
        [
            (ch * (gs + gh) - (cb + ch) * gh) / det,
            (-ch * gh + (cb + ch) * (gh + gl)) / det,
        ],
    ];
    let trace = a[0][0] + a[1][1];
    let gap = ((a[0][0] - a[1][1]).powi(2) + 4.0 * a[0][1] * a[1][0]).sqrt();
    let slow = (trace - gap) / 2.0;
    let fast = (trace + gap) / 2.0;
    let u = [(2.0 * fast - 2.0 * a[0][0]) / gap, -2.0 * a[1][0] / gap];
    let v = [2.0 - u[0], -u[1]];
    let dt = 0.5 * scale;
    let expected =
        core::array::from_fn::<_, 2, _>(|i| u[i] * (-slow * dt).exp() + v[i] * (-fast * dt).exp());
    let integral = |x: f64, y: f64| {
        let decay = |rate: f64| -(-rate * dt).exp_m1() / rate;
        x * x * decay(2.0 * slow) + 2.0 * x * y * decay(slow + fast) + y * y * decay(2.0 * fast)
    };
    let expected_heat = [
        gh * integral(u[0] - u[1], v[0] - v[1]),
        gl * integral(u[1], v[1]),
        gs * integral(u[0], v[0]),
    ];
    let before = p.state(0.0, [2.0; 2], [0.0; 2]).unwrap();
    let mut errors = Vec::new();
    for n in [4, 8, 16] {
        let mut native = prepare(
            p,
            [1.0 / gs; 2],
            &[interval([true, false, true, false], dt)],
            n,
        );
        let (mut next, mut mean) = ([0.0; 5], [0.0; 2]);
        let step = native
            .step_controlled(&before, &[0.0; 2], &[None; 2], &mut next, &mut mean, || {
                false
            })
            .unwrap();
        account(p, &before, &next, &step);
        let actual = physical(p, &next);
        let mut error: f64 = 0.0;
        for leg in 0..2 {
            error = error
                .max((actual[1 + leg] - expected[0]).abs())
                .max((actual[3 + leg] - expected[1]).abs());
            for (got, want) in [
                step.converter.channel_heat_j[2 * leg],
                step.converter.channel_heat_j[2 * leg + 1],
                step.terminal_resistor_heat_j[leg],
            ]
            .into_iter()
            .zip(expected_heat)
            {
                error = error.max((got - want).abs() / scale);
            }
        }
        assert!(actual[3] > 0.8 && actual[3] < actual[1]);
        assert!(step.converter.device_capacitor_energy_change_j[1] > 0.0);
        assert!(step.converter.device_capacitor_energy_change_j[0] < 0.0);
        assert_eq!(step.converter.diode_heat_j, [0.0; 4]);
        close(actual[0], 0.0, 1e-12);
        let energy_before: f64 = energies(p, &before).iter().sum();
        let energy_after: f64 = energies(p, &next).iter().sum();
        close(
            (energy_after - energy_before + heat(&step).iter().sum::<f64>()) / scale,
            0.0,
            2e-9,
        );
        errors.push(error);
    }
    assert!(errors[2] < 3e-4, "scale {scale:e}: {errors:?}");
    for pair in errors.windows(2) {
        assert!(
            (3.8..4.2).contains(&(pair[0] / pair[1])),
            "scale {scale:e}: {errors:?}"
        );
    }
}

#[test]
fn gate_edges_continue_all_five_charges_and_preserve_source_load_correlation() {
    let p = parameters();
    let dt = 1.0 / 1024.0;
    let deck = [
        interval([true, false, true, false], dt),
        interval([false; 4], dt),
        interval([true; 4], dt),
    ];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [8.0, 3.0],
            loads: [
                None,
                Some(RegulatedPowerLoad::try_new(2.0, 1.0, 0.5, 30.0).unwrap()),
            ],
        },
        PowerSwitchCommand {
            source_voltage_v: [7.0, 4.0],
            loads: [None; 2],
        },
        PowerSwitchCommand {
            source_voltage_v: [6.0, 2.0],
            loads: [None; 2],
        },
    ];
    let resistance = [0.5, 0.7];
    let before = p.state(1.0, [7.0, 3.0], [2.0, 1.0]).unwrap();
    let mut together = prepare(p, resistance, &deck, 3);
    let (mut all, mut mean) = ([0.0; 5], [0.0; 2]);
    let total = together
        .step_interval_commands_controlled(&before, &commands, &mut all, &mut mean, || false)
        .unwrap();
    let mut x = before;
    let (mut q, mut source, mut load) = ([0.0; 11], [0.0; 2], [0.0; 2]);
    for (interval, cmd) in deck.into_iter().zip(commands) {
        let mut split = prepare(p, resistance, &[interval], 3);
        let (mut next, mut voltage) = ([0.0; 5], [0.0; 2]);
        let step = split
            .step_controlled(
                &x,
                &cmd.source_voltage_v,
                &cmd.loads,
                &mut next,
                &mut voltage,
                || false,
            )
            .unwrap();
        account(p, &x, &next, &step);
        // Finite capacitance keeps the next interval's node near its supplied
        // physical history, even across dead time and overlap.
        let a = physical(p, &x);
        let b = physical(p, &next);
        assert!((a[3] - b[3]).abs() < 0.3 && (a[4] - b[4]).abs() < 0.3);
        for (q, got) in q.iter_mut().zip(heat(&step)) {
            *q += got;
        }
        for i in 0..2 {
            source[i] += step.source_work_j[i];
            load[i] += step.load_work_j[i];
        }
        x = next;
    }
    assert_eq!(all.map(f64::to_bits), x.map(f64::to_bits));
    for (got, want) in heat(&total).into_iter().zip(q) {
        close(got, want, 2e-12);
    }
    for i in 0..2 {
        close(total.source_work_j[i], source[i], 2e-12);
        close(total.load_work_j[i], load[i], 2e-12);
    }
    assert_eq!(total.integration_steps, 9);
    account(p, &before, &all, &total);
}

#[test]
fn unsupported_capacitance_and_late_commands_or_cancellation_publish_no_state_or_voltage() {
    let p = parameters();
    let dt = 1.0 / 1024.0;
    let deck = [interval([true, false, true, false], dt); 2];
    let before = p.state(1.0, [7.0, 3.0], [2.0, 1.0]).unwrap();
    let mut native = prepare(p, [0.5; 2], &deck, 3);
    let (mut expected, mut mean) = ([0.0; 5], [0.0; 2]);
    let polls = Cell::new(0usize);
    let good = native
        .step_controlled(
            &before,
            &[8.0, 3.0],
            &[None; 2],
            &mut expected,
            &mut mean,
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let stop = polls.get() - 1;
    polls.set(0);
    let (mut next, mut voltage) = ([77.0; 5], [88.0; 2]);
    assert!(matches!(
        native.step_controlled(
            &before,
            &[8.0, 3.0],
            &[None; 2],
            &mut next,
            &mut voltage,
            || {
                polls.set(polls.get() + 1);
                polls.get() >= stop
            }
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [77.0; 5]);
    assert_eq!(voltage, [88.0; 2]);
    let mut commands = [PowerSwitchCommand {
        source_voltage_v: [8.0, 3.0],
        loads: [None; 2],
    }; 2];
    commands[1].loads[1] = Some(RegulatedPowerLoad::try_new(1.0, 1.0, 0.25, 0.5).unwrap());
    assert!(
        native
            .step_interval_commands_controlled(&before, &commands, &mut next, &mut voltage, || {
                false
            })
            .is_err()
    );
    assert_eq!(next, [77.0; 5]);
    assert_eq!(voltage, [88.0; 2]);
    let retry = native
        .step_controlled(
            &before,
            &[8.0, 3.0],
            &[None; 2],
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(voltage.map(f64::to_bits), mean.map(f64::to_bits));
    assert_eq!(heat(&retry), heat(&good));
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut bad = p;
        bad.output_capacitance_f[2] = value;
        assert!(bad.state(1.0, [7.0, 3.0], [2.0, 1.0]).is_err());
        assert!(
            PreparedCapacitiveCommutatedSwitchCycle::new_controlled(
                bad,
                [0.5; 2],
                &deck,
                1,
                electrical(),
                controls(),
                || false
            )
            .is_err()
        );
    }
}
