//! Regulated PWM source/demand correlation, ripple heat, refinement and rollback.
use super::*;
use fs_phs::circuit::converter::switching::power_load::{
    MAX_POWER_SUBSTEPS_PER_INTERVAL, PowerLoadedSwitchCycleStep, PowerSwitchCommand,
    PreparedPowerLoadedSwitchCycle,
};
use fs_phs::circuit::converter::switching::{
    BridgeSwitchState, SwitchCycleControls, SwitchInterval, centered_pwm,
};

fn cycle_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        max_intervals: 64,
        absolute_energy_tolerance_j: 3e-9,
        relative_energy_tolerance: 1e-11,
    }
}
fn interval(a: bool, b: bool, duration_s: f64) -> SwitchInterval {
    SwitchInterval {
        state: BridgeSwitchState {
            input_high: a,
            output_high: b,
        },
        duration_s,
    }
}
fn prepared(
    p: BuckBoostParameters,
    r: [f64; 2],
    intervals: &[SwitchInterval],
    subdivisions: usize,
) -> PreparedPowerLoadedSwitchCycle {
    PreparedPowerLoadedSwitchCycle::new_controlled(
        p,
        r,
        intervals,
        subdivisions,
        controls(),
        cycle_controls(),
        || false,
    )
    .unwrap()
}
fn balance(
    p: BuckBoostParameters,
    before: [f64; 3],
    after: [f64; 3],
    r: &PowerLoadedSwitchCycleStep,
) {
    let energy = |x: [f64; 3]| {
        x[0] * x[0] / (2.0 * p.inductance_h)
            + x[1] * x[1] / (2.0 * p.bus_capacitance_f[0])
            + x[2] * x[2] / (2.0 * p.bus_capacitance_f[1])
    };
    close(r.converter.native.delta_h, energy(after) - energy(before));
    close(
        energy(after) - energy(before)
            + r.converter.conduction.joules.iter().sum::<f64>()
            + r.terminal_resistor_heat_j.iter().sum::<f64>()
            + r.load_work_j.iter().sum::<f64>(),
        r.source_work_j.iter().sum::<f64>(),
    );
    for port in 0..2 {
        close(
            r.source_work_j[port] - r.terminal_resistor_heat_j[port] - r.load_work_j[port],
            r.converter.terminal_work_j[port],
        );
    }
}

#[test]
fn a_single_physical_interval_preserves_the_original_native_result() {
    let p = parameters();
    let intervals = [interval(true, false, 0.003)];
    let sources = [24.0, 0.0];
    let terminal_r = [0.5, 10.0];
    let demands = [None, Some(load(12.0, 10.0, 2.0))];
    let native = model(1.0, 0.0);
    let initial = native.state(1.0, [20.0, 10.0]).unwrap();
    let mut native = native
        .prepare_power_loaded_controlled(0.003, controls(), || false)
        .unwrap();
    let mut cycle = prepared(p, terminal_r, &intervals, 1);
    let (mut expected, mut expected_v) = ([0.0; 3], [0.0; 2]);
    let expected_step = native
        .step_voltage_driven_controlled(
            &initial,
            &sources,
            &terminal_r,
            &demands,
            &mut expected,
            &mut expected_v,
            || false,
        )
        .unwrap();
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let step = cycle
        .step_controlled(
            &initial,
            &sources,
            &demands,
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(voltage.map(f64::to_bits), expected_v.map(f64::to_bits));
    assert_eq!(step.source_work_j, expected_step.source_work_j);
    assert_eq!(
        step.terminal_resistor_heat_j,
        expected_step.terminal_resistor_heat_j
    );
    assert_eq!(step.load_work_j, expected_step.load_work_j);
    assert_eq!(step.mean_source_current_a, expected_step.source_current_a);
    assert_eq!(step.last_loads, expected_step.loads);
    assert_eq!(
        step.converter.native.newton_iters,
        expected_step.converter.native.newton_iters
    );
    assert_eq!(step.integration_steps, 1);
    assert_eq!(cycle.integration_steps(), 1);
    assert_eq!(cycle.period_s(), 0.003);
    balance(p, initial, next, &step);
}

// Independent continuous SI equations, including every physical power owner.
// The chosen reference trajectory stays in constant-power operation; no native
// callback, discrete gradient, power/load report or matrix assembly is reused.
fn power_ode(
    p: BuckBoostParameters,
    r: [f64; 2],
    gates: BridgeSwitchState,
    source: [f64; 2],
    demand_w: [f64; 2],
    y: [f64; 15],
) -> [f64; 15] {
    let (a, b) = (
        if gates.input_high { 1.0 } else { 0.0 },
        if gates.output_high { 1.0 } else { 0.0 },
    );
    let channel = [
        p.winding_resistance_ohm,
        a * p.input_on_resistance_ohm[0],
        (1.0 - a) * p.input_on_resistance_ohm[1],
        b * p.output_on_resistance_ohm[0],
        (1.0 - b) * p.output_on_resistance_ohm[1],
    ];
    let source_i = [(source[0] - y[1]) / r[0], (source[1] - y[2]) / r[1]];
    assert!(y[1] > 2.0 && y[2] > 2.0);
    let load_i = [demand_w[0] / y[1], demand_w[1] / y[2]];
    assert!(load_i.iter().all(|i| *i < 10.0));
    let mut dy = [0.0; 15];
    dy[0] = (a * y[1] - b * y[2] - channel.iter().sum::<f64>() * y[0]) / p.inductance_h;
    dy[1] = (source_i[0] - a * y[0] - load_i[0]) / p.bus_capacitance_f[0];
    dy[2] = (source_i[1] + b * y[0] - load_i[1]) / p.bus_capacitance_f[1];
    for port in 0..2 {
        dy[3 + port] = source[port] * source_i[port];
        dy[10 + port] = r[port] * source_i[port] * source_i[port];
        dy[12 + port] = demand_w[port];
    }
    for channel_index in 0..5 {
        dy[5 + channel_index] = channel[channel_index] * y[0] * y[0];
    }
    dy[14] = y[0];
    dy
}
fn reference(p: BuckBoostParameters, r: [f64; 2], intervals: &[SwitchInterval]) -> [f64; 15] {
    let mut y = [0.0; 15];
    y[..3].copy_from_slice(&[2.0, 15.0, 8.0]);
    for interval in intervals {
        let dt = interval.duration_s / 2048.0;
        for _ in 0..2048 {
            let rhs = |state| power_ode(p, r, interval.state, [18.0, 0.0], [2.0, 15.0], state);
            let k1 = rhs(y);
            let k2 = rhs(core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]));
            let k3 = rhs(core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]));
            let k4 = rhs(core::array::from_fn(|i| y[i] + dt * k3[i]));
            y = core::array::from_fn(|i| {
                y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            });
        }
    }
    y
}

#[test]
fn resolved_pwm_refines_to_an_independent_ode_and_ripple_heat() {
    let p = BuckBoostParameters {
        inductance_h: 0.002,
        bus_capacitance_f: [0.003, 0.002],
        ..parameters()
    };
    let r = [0.4, 30.0];
    let intervals = centered_pwm(
        0.002,
        BridgeModulation {
            input_high_fraction: 0.4,
            output_high_fraction: 0.75,
        },
    )
    .unwrap();
    let initial = [
        p.inductance_h * 2.0,
        p.bus_capacitance_f[0] * 15.0,
        p.bus_capacitance_f[1] * 8.0,
    ];
    let truth = reference(p, r, &intervals);
    let run = |subdivisions| {
        let mut cycle = prepared(p, r, &intervals, subdivisions);
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let record = cycle
            .step_controlled(
                &initial,
                &[18.0, 0.0],
                &[Some(load(2.0, 10.0, 2.0)), Some(load(15.0, 10.0, 2.0))],
                &mut next,
                &mut voltage,
                || false,
            )
            .unwrap();
        balance(p, initial, next, &record);
        assert_eq!(record.integration_steps, intervals.len() * subdivisions);
        let state = [
            next[0] / p.inductance_h,
            next[1] / p.bus_capacitance_f[0],
            next[2] / p.bus_capacitance_f[1],
        ];
        let error = (0..3)
            .map(|i| (state[i] - truth[i]).abs())
            .fold(0.0_f64, f64::max);
        (error, record)
    };
    let (coarse, _) = run(2);
    let (middle, middle_record) = run(4);
    let (fine, fine_record) = run(8);
    assert!(middle < 0.35 * coarse, "{coarse:e} -> {middle:e}");
    assert!(fine < 0.35 * middle, "{middle:e} -> {fine:e}");
    let mut middle_heat_error = 0.0_f64;
    let mut fine_heat_error = 0.0_f64;
    for channel in 0..7 {
        let (middle_heat, fine_heat, exact) = if channel < 5 {
            (
                middle_record.converter.conduction.joules[channel],
                fine_record.converter.conduction.joules[channel],
                truth[5 + channel],
            )
        } else {
            (
                middle_record.terminal_resistor_heat_j[channel - 5],
                fine_record.terminal_resistor_heat_j[channel - 5],
                truth[10 + channel - 5],
            )
        };
        middle_heat_error = middle_heat_error.max((middle_heat - exact).abs());
        fine_heat_error = fine_heat_error.max((fine_heat - exact).abs());
    }
    assert!(fine_heat_error < 0.35 * middle_heat_error);
    close(fine_record.load_work_j[0], 2.0 * 0.002);
    close(fine_record.load_work_j[1], 15.0 * 0.002);
    // Resolved positive variance adds real winding heat above the square of
    // the mean current. A cycle-averaged-current heat shortcut misses it.
    let mean_i = truth[14] / 0.002;
    let mean_square_heat = 0.002 * p.winding_resistance_ohm * mean_i * mean_i;
    assert!(truth[5] > mean_square_heat * 1.01);
    assert!(fine_record.converter.conduction.joules[0] > mean_square_heat * 1.005);
}

#[test]
fn repeated_gate_caches_keep_interval_source_and_load_commands_correlated() {
    let intervals = [interval(false, false, 0.01); 3];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [0.0, 20.0],
            loads: [None, Some(load(20.0, 10.0, 2.0))],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 5.0],
            loads: [None, None],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 15.0],
            loads: [None, Some(load(10.0, 10.0, 2.0))],
        },
    ];
    let mut cycle = prepared(parameters(), [1.0, 2.0], &intervals, 3);
    assert!(cycle.cached_interval_count() <= 2);
    assert_eq!(cycle.integration_steps(), 9);
    let initial = [0.0, 0.0, 20.0];
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let record = cycle
        .step_interval_commands_controlled(&initial, &commands, &mut next, &mut voltage, || false)
        .unwrap();
    // Independently solve each capacitor midpoint quadratic, preserving the
    // physical command order. Three equal steps agree to representable roundoff
    // with the owner's compensated last remainder.
    let (mut v, mut work, mut heat, mut charge, mut integral_v) = (10.0, 0.0, 0.0, 0.0, 0.0);
    for command in commands {
        let power = command.loads[1].map_or(0.0, |load| load.requested_power_w());
        for _ in 0..3 {
            let h = 0.01 / 3.0;
            let aa = 2.0 * 2.0 / h + 0.5;
            let bb = 2.0 * 2.0 * v / h + command.source_voltage_v[1] * 0.5;
            let midpoint = (bb + (bb * bb - 4.0 * aa * power).sqrt()) / (2.0 * aa);
            let current = (command.source_voltage_v[1] - midpoint) / 2.0;
            work += h * command.source_voltage_v[1] * current;
            heat += h * 2.0 * current * current;
            charge += h * current;
            integral_v += h * midpoint;
            v = 2.0 * midpoint - v;
        }
    }
    close(next[2] / 2.0, v);
    close(record.source_work_j[1], work);
    close(record.terminal_resistor_heat_j[1], heat);
    close(record.mean_source_current_a[1], charge / 0.03);
    close(voltage[1], integral_v / 0.03);
    close(record.load_work_j[1], 0.3);
    assert_eq!(record.last_loads[1].requested_power_w, 10.0);
    let averaged_source_v = (20.0 + 5.0 + 15.0) / 3.0;
    assert!(
        (record.source_work_j[1] - averaged_source_v * record.mean_source_current_a[1] * 0.03)
            .abs()
            > 0.1
    );
    balance(parameters(), initial, next, &record);
}

#[test]
fn cycle_curtailment_uses_each_actual_load_branch_and_retains_only_the_final_point() {
    let intervals = [interval(false, false, 0.02); 3];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [0.0, 6.0],
            loads: [None, Some(load(100.0, 0.2, 1.0))],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 6.0],
            loads: [None, Some(load(4.0, 10.0, 10.0))],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 6.0],
            loads: [None, None],
        },
    ];
    let mut cycle = prepared(parameters(), [1.0, 2.0], &intervals, 2);
    let initial = [0.0, 0.0, 8.0];
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let record = cycle
        .step_interval_commands_controlled(&initial, &commands, &mut next, &mut voltage, || false)
        .unwrap();
    let (mut v, mut delivered) = (4.0, 0.0);
    for (conductance, fixed_i) in [(0.0, 0.2), (0.04, 0.0), (0.0, 0.0)] {
        for _ in 0..2 {
            let midpoint = (2.0 * 2.0 * v / 0.01 + 6.0 / 2.0 - fixed_i)
                / (2.0 * 2.0 / 0.01 + 0.5 + conductance);
            delivered += 0.01 * midpoint * (fixed_i + conductance * midpoint);
            v = 2.0 * midpoint - v;
        }
    }
    close(next[2] / 2.0, v);
    close(record.load_work_j[1], delivered);
    close(
        record.curtailed_load_energy_j[1],
        0.02 * (100.0 + 4.0) - delivered,
    );
    assert!(delivered < 0.1);
    assert_eq!(record.last_loads[1].mode, PowerLoadMode::Disabled);
    balance(parameters(), initial, next, &record);
}

#[test]
fn subdivision_limits_late_commands_and_cancellation_preserve_outputs_and_retry() {
    let intervals = [interval(false, false, 0.01); 2];
    for subdivisions in [0, MAX_POWER_SUBSTEPS_PER_INTERVAL + 1] {
        assert!(
            PreparedPowerLoadedSwitchCycle::new_controlled(
                parameters(),
                [1.0; 2],
                &intervals,
                subdivisions,
                controls(),
                cycle_controls(),
                || false,
            )
            .is_err()
        );
    }
    assert!(
        PreparedPowerLoadedSwitchCycle::new_controlled(
            parameters(),
            [1.0; 2],
            &[interval(false, false, f64::from_bits(1))],
            2,
            controls(),
            cycle_controls(),
            || false,
        )
        .is_err()
    );
    let max = prepared(
        parameters(),
        [1.0; 2],
        &[interval(false, false, 0.01); 64],
        MAX_POWER_SUBSTEPS_PER_INTERVAL,
    );
    assert_eq!(
        max.integration_steps(),
        64 * MAX_POWER_SUBSTEPS_PER_INTERVAL
    );
    let mut cycle = prepared(parameters(), [1.0, 2.0], &intervals, 3);
    let initial = [0.0, 0.0, 20.0];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [0.0, 12.0],
            loads: [None, Some(load(12.0, 10.0, 2.0))],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 15.0],
            loads: [None, Some(load(24.0, 10.0, 2.0))],
        },
    ];
    let (mut expected, mut expected_v) = ([0.0; 3], [0.0; 2]);
    let mut polls = 0;
    cycle
        .step_interval_commands_controlled(
            &initial,
            &commands,
            &mut expected,
            &mut expected_v,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut next, mut voltage) = ([91.0; 3], [73.0; 2]);
    for stop in [polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            cycle.step_interval_commands_controlled(
                &initial,
                &commands,
                &mut next,
                &mut voltage,
                || {
                    seen += 1;
                    seen == stop
                },
            ),
            Err(PreparedStepError::Cancelled)
        ));
        assert_eq!(next, [91.0; 3]);
        assert_eq!(voltage, [73.0; 2]);
    }
    let mut invalid = commands;
    invalid[1].source_voltage_v[1] = f64::NAN;
    assert!(
        cycle
            .step_interval_commands_controlled(&initial, &invalid, &mut next, &mut voltage, || {
                false
            })
            .is_err()
    );
    assert!(
        cycle
            .step_interval_commands_controlled(
                &initial,
                &commands[..1],
                &mut next,
                &mut voltage,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    cycle
        .step_interval_commands_controlled(&initial, &commands, &mut next, &mut voltage, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(voltage.map(f64::to_bits), expected_v.map(f64::to_bits));
}
