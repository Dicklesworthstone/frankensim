//! Synthetic, explicitly authored conduction data and independent physical oracles.
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutatedStep, CommutationInterval, CommutationParameters, DeviceGateState,
    DeviceOperatingBounds, HalfBridgeConduction, HalfBridgePoint, PreparedCommutatedSwitchCycle,
    SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::power_load::PowerSwitchCommand;
use fs_phs::{PreparedStepError, SolveBudget};

fn device() -> SwitchDevice {
    SwitchDevice {
        on_resistance_ohm: 0.2,
        off_resistance_ohm: 1000.0,
        body_diode: BodyDiode {
            forward_voltage_v: 0.7,
            differential_resistance_ohm: 0.1,
        },
        bounds: DeviceOperatingBounds {
            maximum_reverse_voltage_v: 100.0,
            maximum_channel_current_a: 100.0,
            maximum_diode_current_a: 100.0,
        },
    }
}
fn parameters() -> CommutationParameters {
    CommutationParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        winding_resistance_ohm: 0.05,
        devices: [device(); 4],
    }
}
fn controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-12,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
    }
}
fn cycle_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        max_intervals: 64,
        absolute_energy_tolerance_j: 3e-9,
        relative_energy_tolerance: 1e-11,
    }
}
fn interval(on: [bool; 4], duration_s: f64) -> CommutationInterval {
    CommutationInterval {
        gates: DeviceGateState { on },
        duration_s,
    }
}
fn prepare(
    p: CommutationParameters,
    resistance: [f64; 2],
    schedule: &[CommutationInterval],
    subdivisions: usize,
) -> PreparedCommutatedSwitchCycle {
    PreparedCommutatedSwitchCycle::new_controlled(
        p,
        resistance,
        schedule,
        subdivisions,
        controls(),
        cycle_controls(),
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}
fn leg_balance(point: HalfBridgePoint, bus: f64, draw: f64) {
    close(
        point.high.terminal_current_a - point.low.terminal_current_a,
        draw,
        1e-12,
    );
    let heat = point.high.channel_power_w
        + point.high.diode_power_w
        + point.low.channel_power_w
        + point.low.diode_power_w;
    close(
        bus * point.high.terminal_current_a - point.switch_node_voltage_v * draw,
        heat,
        1e-11,
    );
    for d in [point.high, point.low] {
        assert!(d.channel_power_w >= 0.0 && d.diode_power_w >= 0.0);
        close(
            d.channel_power_w,
            d.drain_source_voltage_v * d.channel_current_a,
            1e-12,
        );
        close(
            d.diode_power_w,
            -d.drain_source_voltage_v * d.diode_current_a,
            1e-12,
        );
    }
}

#[test]
fn original_kcl_closed_forms_cover_both_diodes_zero_current_and_gate_overlap() {
    let d = device();
    let off = HalfBridgeConduction::new(d, d, [false; 2]).unwrap();
    let zero = off.evaluate(12.0, 0.0).unwrap();
    close(zero.switch_node_voltage_v, 6.0, 1e-13);
    close(zero.high.terminal_current_a, 0.006, 1e-14);
    close(
        zero.high.channel_power_w + zero.low.channel_power_w,
        0.072,
        1e-14,
    );
    assert_eq!(
        (zero.high.diode_current_a, zero.low.diode_current_a),
        (0.0, 0.0)
    );
    leg_balance(zero, 12.0, 0.0);

    let positive = off.evaluate(12.0, 2.0).unwrap();
    close(
        positive.switch_node_voltage_v,
        (0.001 * 12.0 - 10.0 * 0.7 - 2.0) / 10.002,
        1e-13,
    );
    assert_eq!(positive.high.diode_current_a, 0.0);
    assert!(positive.low.diode_current_a > 1.9);
    leg_balance(positive, 12.0, 2.0);
    let negative = off.evaluate(12.0, -2.0).unwrap();
    close(
        negative.switch_node_voltage_v,
        ((10.001 * 12.0) + 7.0 + 2.0) / 10.002,
        1e-13,
    );
    close(
        positive.switch_node_voltage_v + negative.switch_node_voltage_v,
        12.0,
        1e-13,
    );
    assert!(negative.high.diode_current_a > 1.9);
    assert_eq!(negative.low.diode_current_a, 0.0);
    leg_balance(negative, 12.0, -2.0);

    // A negative supported bus can forward-bias BOTH diodes simultaneously.
    let both_diodes = off.evaluate(-4.0, 0.0).unwrap();
    close(both_diodes.switch_node_voltage_v, -2.0, 1e-13);
    close(both_diodes.high.diode_current_a, 13.0, 1e-12);
    close(both_diodes.low.diode_current_a, 13.0, 1e-12);
    close(-4.0 * both_diodes.high.terminal_current_a, 52.008, 1e-11);
    leg_balance(both_diodes, -4.0, 0.0);

    // An ON channel remains bidirectional and shares reverse current with its diode.
    let on = HalfBridgeConduction::new(d, d, [true, false]).unwrap();
    let forward = on.evaluate(12.0, 2.0).unwrap();
    close(forward.switch_node_voltage_v, (60.0 - 2.0) / 5.001, 1e-13);
    assert_eq!(forward.high.diode_current_a, 0.0);
    leg_balance(forward, 12.0, 2.0);
    let reverse = on.evaluate(12.0, -20.0).unwrap();
    close(reverse.switch_node_voltage_v, 207.0 / 15.001, 1e-13);
    assert!(reverse.high.channel_current_a < -8.0);
    assert!(reverse.high.diode_current_a > 10.0);
    leg_balance(reverse, 12.0, -20.0);

    let shoot_through = HalfBridgeConduction::new(d, d, [true; 2])
        .unwrap()
        .evaluate(12.0, 0.0)
        .unwrap();
    close(shoot_through.high.channel_current_a, 30.0, 1e-12);
    close(shoot_through.low.channel_current_a, 30.0, 1e-12);
    close(
        shoot_through.high.channel_power_w + shoot_through.low.channel_power_w,
        360.0,
        1e-11,
    );
    leg_balance(shoot_through, 12.0, 0.0);
}

#[test]
fn eliminated_node_tangent_is_the_original_differential_response() {
    for (gates, bus, draw) in [
        ([false; 2], 12.0, 0.0),
        ([false; 2], 12.0, 2.0),
        ([false; 2], 12.0, -2.0),
        ([false; 2], -4.0, 0.0),
        ([true, false], 12.0, -20.0),
        ([true; 2], 12.0, 3.0),
    ] {
        let owner = HalfBridgeConduction::new(device(), device(), gates).unwrap();
        let p = owner.evaluate(bus, draw).unwrap();
        let h = 1e-5;
        let vm = owner.evaluate(bus - h, draw).unwrap();
        let vp = owner.evaluate(bus + h, draw).unwrap();
        let im = owner.evaluate(bus, draw - h).unwrap();
        let ip = owner.evaluate(bus, draw + h).unwrap();
        close(
            (vp.switch_node_voltage_v - vm.switch_node_voltage_v) / (2.0 * h),
            p.d_node_voltage_d_bus_voltage,
            2e-9,
        );
        close(
            (ip.switch_node_voltage_v - im.switch_node_voltage_v) / (2.0 * h),
            p.d_node_voltage_d_draw_current_ohm,
            2e-9,
        );
        close(
            (vp.high.terminal_current_a - vm.high.terminal_current_a) / (2.0 * h),
            p.d_high_current_d_bus_voltage_s,
            2e-8,
        );
        close(
            (ip.high.terminal_current_a - im.high.terminal_current_a) / (2.0 * h),
            p.d_high_current_d_draw_current,
            2e-8,
        );
    }
}

fn physical_energy(p: CommutationParameters, state: [f64; 3]) -> f64 {
    0.5 * state[0] * state[0] / p.inductance_h
        + 0.5 * state[1] * state[1] / p.bus_capacitance_f[0]
        + 0.5 * state[2] * state[2] / p.bus_capacitance_f[1]
}
fn account(p: CommutationParameters, initial: [f64; 3], next: [f64; 3], r: &CommutatedStep) {
    let heat = r.winding_heat_j
        + r.channel_heat_j.iter().sum::<f64>()
        + r.diode_heat_j.iter().sum::<f64>();
    close(
        r.native.delta_h,
        physical_energy(p, next) - physical_energy(p, initial),
        2e-10,
    );
    close(r.native.dissipated, heat, 2e-10);
    close(r.native.delta_h + heat, r.native.supplied, 2e-10);
    close(r.native.supplied, r.terminal_work_j.iter().sum(), 2e-10);
    close(
        r.native.delta_h
            + heat
            + r.terminal_resistor_heat_j.iter().sum::<f64>()
            + r.load_work_j.iter().sum::<f64>(),
        r.source_work_j.iter().sum(),
        2e-10,
    );
    for k in 0..2 {
        close(
            r.source_work_j[k] - r.terminal_resistor_heat_j[k] - r.load_work_j[k],
            r.terminal_work_j[k],
            2e-10,
        );
    }
    assert!(
        r.channel_heat_j
            .iter()
            .chain(&r.diode_heat_j)
            .chain(&r.terminal_resistor_heat_j)
            .chain(&r.load_work_j)
            .all(|x| *x >= 0.0)
    );
}

#[test]
fn dead_time_commutates_both_current_directions_into_separately_owned_diodes() {
    let p = parameters();
    let dt = 1e-5;
    for sign in [-1.0, 1.0] {
        let initial = p.state(sign * 2.0, [12.0, 8.0]).unwrap();
        let mut cycle = prepare(p, [2.0; 2], &[interval([false; 4], dt)], 1);
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let r = cycle
            .step_controlled(
                &initial,
                &[12.0, 8.0],
                &[None; 2],
                &mut next,
                &mut voltage,
                || false,
            )
            .unwrap();
        let conducting = if sign > 0.0 { [1, 2] } else { [0, 3] };
        for d in 0..4 {
            assert_eq!(r.diode_heat_j[d] > 0.0, conducting.contains(&d));
            close(
                r.diode_heat_j[d],
                dt * (-r.last_devices[d].drain_source_voltage_v)
                    * r.last_devices[d].diode_current_a,
                1e-15,
            );
            close(
                r.channel_heat_j[d],
                dt * r.last_devices[d].drain_source_voltage_v.powi(2)
                    / p.devices[d].off_resistance_ohm,
                1e-15,
            );
        }
        assert!(r.winding_heat_j > 0.0);
        assert_eq!(r.integration_steps, 1);
        for k in 0..2 {
            close(
                voltage[k],
                0.5 * (initial[k + 1] + next[k + 1]) / p.bus_capacitance_f[k],
                1e-12,
            );
        }
        account(p, initial, next, &r);
    }
}

#[test]
fn deliberate_overlap_draws_finite_shoot_through_energy_from_both_capacitors() {
    let mut p = parameters();
    p.bus_capacitance_f = [0.2; 2];
    let initial = p.state(0.0, [12.0; 2]).unwrap();
    let dt = 0.005;
    let mut cycle = prepare(p, [2.0; 2], &[interval([true; 4], dt)], 1);
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let r = cycle
        .step_controlled(
            &initial,
            &[0.0; 2],
            &[None; 2],
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    // Independent parallel RC equation: each bridge is two 0.2-ohm channels
    // in series across its capacitor, with a separate 2-ohm source resistor.
    let conductance = 1.0 / 0.4 + 1.0 / 2.0;
    let expected = 12.0 * (0.2 - dt * conductance / 2.0) / (0.2 + dt * conductance / 2.0);
    close(next[0], 0.0, 1e-12);
    for k in 0..2 {
        close(next[k + 1] / 0.2, expected, 1e-11);
        close(
            r.terminal_resistor_heat_j[k],
            dt * voltage[k].powi(2) / 2.0,
            1e-12,
        );
    }
    for d in 0..4 {
        close(
            r.channel_heat_j[d],
            dt * (voltage[d / 2] / 2.0).powi(2) / 0.2,
            1e-12,
        );
    }
    assert_eq!(r.diode_heat_j, [0.0; 4]);
    assert_eq!(r.source_work_j, [0.0; 2]);
    close(r.winding_heat_j, 0.0, 1e-20);
    account(p, initial, next, &r);
}

#[test]
fn manufactured_steady_diode_transfer_retains_original_flow_and_power_conjugacy() {
    let mut p = parameters();
    p.winding_resistance_ohm = 0.1;
    for d in &mut p.devices {
        d.on_resistance_ohm = 0.1;
        d.off_resistance_ohm = 100.0;
        d.body_diode = BodyDiode {
            forward_voltage_v: 0.5,
            differential_resistance_ohm: 0.2,
        };
    }
    let current = 1.0;
    let output = 2.0;
    let right_node = ((0.01 + 5.0) * output + 5.0 * 0.5 + current) / 5.02;
    let left_node = right_node + p.winding_resistance_ohm * current;
    let input = left_node + (current + 0.01 * left_node) * 0.1;
    let high_left = current + 0.01 * left_node;
    let high_right = 0.01 * (output - right_node) - (right_node - output - 0.5) / 0.2;
    let source = [input + high_left, output + high_right + 0.5 / output];
    let load = RegulatedPowerLoad::try_new(0.5, 5.0, 0.5, 100.0).unwrap();
    let initial = p.state(current, [input, output]).unwrap();
    let dt = 0.01;
    let mut cycle = prepare(p, [1.0; 2], &[interval([true, false, false, false], dt)], 1);
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let r = cycle
        .step_controlled(
            &initial,
            &source,
            &[None, Some(load)],
            &mut next,
            &mut mean,
            || false,
        )
        .unwrap();
    for k in 0..3 {
        close(next[k], initial[k], 2e-12);
    }
    close(
        r.last_devices[2].diode_current_a,
        (right_node - output - 0.5) / 0.2,
        2e-12,
    );
    close(
        r.diode_heat_j[2],
        dt * (right_node - output) * (right_node - output - 0.5) / 0.2,
        2e-12,
    );
    close(r.load_work_j[1], dt * 0.5, 2e-12);
    account(p, initial, next, &r);
}

// Independent test-only root solve of the original scalar KCL, deliberately
// unrelated to the production four-affine-candidate elimination or tangent.
fn reference_leg(
    devices: [SwitchDevice; 2],
    on: [bool; 2],
    bus: f64,
    draw: f64,
) -> (f64, f64, [f64; 2], [f64; 2]) {
    let g: [f64; 2] = core::array::from_fn(|k| {
        1.0 / if on[k] {
            devices[k].on_resistance_ohm
        } else {
            devices[k].off_resistance_ohm
        }
    });
    let currents = |node: f64| {
        let v = [bus - node, node];
        let diode: [f64; 2] = core::array::from_fn(|k| {
            ((-v[k] - devices[k].body_diode.forward_voltage_v)
                / devices[k].body_diode.differential_resistance_ohm)
                .max(0.0)
        });
        (v, diode, [g[0] * v[0] - diode[0], g[1] * v[1] - diode[1]])
    };
    let (mut lo, mut hi) = (-1000.0, 1000.0);
    for _ in 0..64 {
        let mid = (lo + hi) / 2.0;
        let (_, _, i) = currents(mid);
        if i[0] - i[1] > draw {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let node = (lo + hi) / 2.0;
    let (v, diode, i) = currents(node);
    (
        node,
        i[0],
        core::array::from_fn(|k| g[k] * v[k] * v[k]),
        core::array::from_fn(|k| -v[k] * diode[k]),
    )
}
fn reference_rate(
    p: CommutationParameters,
    resistance: [f64; 2],
    gates: [bool; 4],
    command: PowerSwitchCommand,
    y: &[f64; 18],
) -> [f64; 18] {
    let left = reference_leg(
        [p.devices[0], p.devices[1]],
        [gates[0], gates[1]],
        y[1],
        y[0],
    );
    let right = reference_leg(
        [p.devices[2], p.devices[3]],
        [gates[2], gates[3]],
        y[2],
        -y[0],
    );
    let source: [f64; 2] =
        core::array::from_fn(|k| (command.source_voltage_v[k] - y[k + 1]) / resistance[k]);
    let loads: [f64; 2] = core::array::from_fn(|k| {
        command.loads[k].map_or(0.0, |load| {
            let v = y[k + 1];
            let unlimited = if v.abs() <= load.foldback_voltage_v() {
                load.requested_power_w() * v.abs() / load.foldback_voltage_v().powi(2)
            } else {
                load.requested_power_w() / v.abs()
            };
            unlimited.min(load.maximum_current_a()).copysign(v)
        })
    });
    let mut d = [0.0; 18];
    d[0] = (left.0 - right.0 - p.winding_resistance_ohm * y[0]) / p.inductance_h;
    d[1] = (source[0] - left.1 - loads[0]) / p.bus_capacitance_f[0];
    d[2] = (source[1] - right.1 - loads[1]) / p.bus_capacitance_f[1];
    d[3] = p.winding_resistance_ohm * y[0] * y[0];
    d[4..8].copy_from_slice(&[left.2[0], left.2[1], right.2[0], right.2[1]]);
    d[8..12].copy_from_slice(&[left.3[0], left.3[1], right.3[0], right.3[1]]);
    for k in 0..2 {
        d[12 + k] = resistance[k] * source[k] * source[k];
        d[14 + k] = command.source_voltage_v[k] * source[k];
        d[16 + k] = y[k + 1] * loads[k];
    }
    d
}
fn reference_cycle(
    p: CommutationParameters,
    resistance: [f64; 2],
    schedule: &[CommutationInterval],
    commands: &[PowerSwitchCommand],
    physical: [f64; 3],
) -> [f64; 18] {
    let mut y = [0.0; 18];
    y[..3].copy_from_slice(&physical);
    for (interval, command) in schedule.iter().zip(commands) {
        let dt = interval.duration_s / 2048.0;
        for _ in 0..2048 {
            let f = |s: &[f64; 18]| reference_rate(p, resistance, interval.gates.on, *command, s);
            let a = f(&y);
            let b = f(&core::array::from_fn(|k| y[k] + dt * a[k] / 2.0));
            let c = f(&core::array::from_fn(|k| y[k] + dt * b[k] / 2.0));
            let d = f(&core::array::from_fn(|k| y[k] + dt * c[k]));
            for k in 0..18 {
                y[k] += dt * (a[k] + 2.0 * b[k] + 2.0 * c[k] + d[k]) / 6.0;
            }
        }
    }
    y
}
fn packed(p: CommutationParameters, state: [f64; 3], r: &CommutatedStep) -> [f64; 18] {
    let mut y = [0.0; 18];
    y[..3].copy_from_slice(&[
        state[0] / p.inductance_h,
        state[1] / p.bus_capacitance_f[0],
        state[2] / p.bus_capacitance_f[1],
    ]);
    y[3] = r.winding_heat_j;
    y[4..8].copy_from_slice(&r.channel_heat_j);
    y[8..12].copy_from_slice(&r.diode_heat_j);
    y[12..14].copy_from_slice(&r.terminal_resistor_heat_j);
    y[14..16].copy_from_slice(&r.source_work_j);
    y[16..18].copy_from_slice(&r.load_work_j);
    y
}

#[test]
fn refined_gate_cycle_converges_to_independent_kcl_rk4_and_retains_every_heat_owner() {
    let p = parameters();
    let resistance = [0.4, 1.2];
    let schedule = [
        interval([true, false, true, false], 0.0008),
        interval([false; 4], 0.001),
        interval([true, false, false, true], 0.0008),
        interval([true, true, true, false], 0.0003),
    ];
    let demand = |p| Some(RegulatedPowerLoad::try_new(p, 10.0, 1.0, 100.0).unwrap());
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [13.0, 7.5],
            loads: [demand(6.0), demand(10.0)],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 9.0],
            loads: [None, demand(16.0)],
        },
        PowerSwitchCommand {
            source_voltage_v: [14.0, 6.0],
            loads: [demand(4.0), demand(8.0)],
        },
        PowerSwitchCommand {
            source_voltage_v: [13.0, 7.5],
            loads: [demand(6.0), demand(10.0)],
        },
    ];
    let initial = p.state(2.0, [12.0, 8.0]).unwrap();
    let reference = reference_cycle(p, resistance, &schedule, &commands, [2.0, 12.0, 8.0]);
    let mut errors = [0.0_f64; 2];
    for (index, subdivisions) in [4, 8].into_iter().enumerate() {
        let mut cycle = prepare(p, resistance, &schedule, subdivisions);
        assert_eq!(cycle.integration_steps(), schedule.len() * subdivisions);
        close(cycle.period_s(), 0.0029, 1e-18);
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        let r = cycle
            .step_interval_commands_controlled(&initial, &commands, &mut next, &mut mean, || false)
            .unwrap();
        assert_eq!(r.integration_steps, schedule.len() * subdivisions);
        let result = packed(p, next, &r);
        for k in 0..18 {
            errors[index] = errors[index].max((result[k] - reference[k]).abs());
        }
        for k in 3..18 {
            close(result[k], reference[k], 3e-6);
        }
        assert!(r.diode_heat_j[1] > 0.0 && r.diode_heat_j[2] > 0.0);
        assert!(r.channel_heat_j.iter().all(|h| *h > 0.0));
        account(p, initial, next, &r);
    }
    assert!(
        errors[1] < errors[0] * 0.30,
        "second-order refinement: {errors:?}"
    );
    assert!(errors[1] < 3e-6, "fine trajectory error: {errors:?}");
}

#[test]
fn source_limits_are_checked_at_endpoints_and_hidden_gate_intervals_without_clipping() {
    let mut p = parameters();
    p.bus_capacitance_f = [0.2; 2];
    for d in &mut p.devices {
        d.bounds.maximum_reverse_voltage_v = 15.0;
        d.bounds.maximum_channel_current_a = 1000.0;
    }
    let initial = p.state(0.0, [10.0; 2]).unwrap();
    let mut cycle = prepare(p, [0.1; 2], &[interval([true; 4], 0.01)], 1);
    // Independent RC endpoint=43 1/3 V, midpoint=26 2/3 V. The actual
    // channel reverse voltage at the endpoint leaves 15 V support while its
    // initial and accepted midpoint values would both remain supported.
    let leg = HalfBridgeConduction::new(p.devices[0], p.devices[1], [true; 2]).unwrap();
    assert!(leg.evaluate(10.0, 0.0).is_ok());
    assert!(leg.evaluate(80.0 / 3.0, 0.0).is_ok());
    assert!(leg.evaluate(130.0 / 3.0, 0.0).is_err());
    let (mut next, mut mean) = ([37.0; 3], [41.0; 2]);
    assert!(
        cycle
            .step_controlled(
                &initial,
                &[100.0; 2],
                &[None; 2],
                &mut next,
                &mut mean,
                || false
            )
            .is_err()
    );
    assert_eq!((next, mean), ([37.0; 3], [41.0; 2]));

    let mut p = parameters();
    for d in &mut p.devices {
        d.bounds.maximum_channel_current_a = 20.0;
    }
    let initial = p.state(0.0, [12.0; 2]).unwrap();
    let mut cycle = prepare(
        p,
        [2.0; 2],
        &[
            interval([false; 4], 1e-6),
            interval([true; 4], 1e-6),
            interval([false; 4], 1e-6),
        ],
        2,
    );
    assert!(
        cycle
            .step_controlled(
                &initial,
                &[12.0; 2],
                &[None; 2],
                &mut next,
                &mut mean,
                || false
            )
            .is_err()
    );
    assert_eq!((next, mean), ([37.0; 3], [41.0; 2]));
    let mut d = device();
    for resistance in [0.0, f64::INFINITY, f64::NAN] {
        d.off_resistance_ohm = resistance;
        assert!(d.validate().is_err());
    }
    assert!(
        PreparedCommutatedSwitchCycle::new_controlled(
            parameters(),
            [1.0; 2],
            &[interval([false; 4], 0.01)],
            0,
            controls(),
            cycle_controls(),
            || false
        )
        .is_err()
    );
    assert!(
        PreparedCommutatedSwitchCycle::new_controlled(
            parameters(),
            [1.0; 2],
            &[interval([false; 4], f64::from_bits(1))],
            2,
            controls(),
            cycle_controls(),
            || false
        )
        .is_err()
    );
}

#[test]
fn later_command_failure_and_cancellation_leave_outputs_untouched_and_retry_identical() {
    let p = parameters();
    let schedule = [
        interval([true, false, true, false], 0.0007),
        interval([false; 4], 0.0004),
        interval([true, false, true, false], 0.0007),
    ];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [13.0, 7.0],
            loads: [None; 2],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 9.0],
            loads: [None; 2],
        },
        PowerSwitchCommand {
            source_voltage_v: [11.0, 8.0],
            loads: [None; 2],
        },
    ];
    let initial = p.state(2.0, [12.0, 8.0]).unwrap();
    let mut cycle = prepare(p, [2.0; 2], &schedule, 3);
    assert!(cycle.cached_interval_count() < 2 * schedule.len());
    let (mut expected, mut expected_mean, mut polls) = ([0.0; 3], [0.0; 2], 0usize);
    let expected_report = cycle
        .step_interval_commands_controlled(
            &initial,
            &commands,
            &mut expected,
            &mut expected_mean,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut bad = commands;
    bad[2].source_voltage_v[0] = f64::NAN;
    let (mut next, mut mean) = ([37.0; 3], [41.0; 2]);
    assert!(
        cycle
            .step_interval_commands_controlled(&initial, &bad, &mut next, &mut mean, || false)
            .is_err()
    );
    assert_eq!((next, mean), ([37.0; 3], [41.0; 2]));
    for stop in [polls / 2, polls] {
        let mut count = 0;
        let result = cycle.step_interval_commands_controlled(
            &initial,
            &commands,
            &mut next,
            &mut mean,
            || {
                count += 1;
                count >= stop
            },
        );
        assert!(matches!(result, Err(PreparedStepError::Cancelled)));
        assert_eq!((next, mean), ([37.0; 3], [41.0; 2]));
    }
    let retry = cycle
        .step_interval_commands_controlled(&initial, &commands, &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(mean.map(f64::to_bits), expected_mean.map(f64::to_bits));
    assert_eq!(
        retry.channel_heat_j.map(f64::to_bits),
        expected_report.channel_heat_j.map(f64::to_bits)
    );
    assert_eq!(
        retry.diode_heat_j.map(f64::to_bits),
        expected_report.diode_heat_j.map(f64::to_bits)
    );
    assert_eq!(
        retry.source_work_j.map(f64::to_bits),
        expected_report.source_work_j.map(f64::to_bits)
    );
    account(p, initial, next, &retry);
}
