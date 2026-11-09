//! Actual switched load work enters the original finite thermal receivers.
use super::*;
use fs_phs::circuit::converter::switching::power_load::PowerSwitchCommand;

fn switch_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-11,
        ..SwitchCycleControls::default()
    }
}

#[test]
fn interval_commands_keep_independent_capacitor_work_and_exact_shared_calorimetry() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let initial = native.state(0.0, [0.0, 10.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let intervals = [
        SwitchInterval {
            state: BridgeSwitchState {
                input_high: false,
                output_high: false,
            },
            duration_s: 1.0 / 128.0,
        },
        SwitchInterval {
            state: BridgeSwitchState {
                input_high: true,
                output_high: false,
            },
            duration_s: 1.0 / 64.0,
        },
    ];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [0.0, 12.0],
            loads: [None, Some(load(24.0))],
        },
        PowerSwitchCommand {
            source_voltage_v: [0.0, 8.0],
            loads: [None, Some(load(8.0))],
        },
    ];
    let subdivisions = 4;
    let period = intervals.iter().map(|i| i.duration_s).sum::<f64>();
    let run = |receiver| {
        network
            .step_power_loaded_switching_candidate(
                &thermal,
                period,
                &converter,
                &intervals,
                &commands,
                subdivisions,
                [None, receiver],
                power_controls(),
                switch_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .unwrap()
    };
    let exported = run(None);
    let heated = run(Some(ALL_IDS[0]));

    // The input stays at zero, so either input gate leaves the winding at rest.
    // Each output bus interval independently solves the capacitor/CPL quadratic.
    let c = parameters().bus_capacitance_f[1];
    let (mut voltage, mut supplied, mut resistor_heat, mut delivered) = (10.0, 0.0, 0.0, 0.0);
    for (index, interval) in intervals.iter().enumerate() {
        let dt = interval.duration_s / subdivisions as f64;
        let source_v = commands[index].source_voltage_v[1];
        let power = [24.0, 8.0][index];
        for _ in 0..subdivisions {
            let a = 2.0 * c / dt + 1.0;
            let b = 2.0 * c * voltage / dt + source_v;
            let effort = (b + (b * b - 4.0 * a * power).sqrt()) / (2.0 * a);
            let current = source_v - effort;
            voltage = 2.0 * effort - voltage;
            supplied += dt * source_v * current;
            resistor_heat += dt * current * current;
            delivered += dt * power;
        }
    }
    let report = heated.report();
    near(heated.converter().state()[2], c * voltage);
    near(
        report.electrical.converter.native.delta_h,
        0.5 * c * (voltage * voltage - 100.0),
    );
    near(report.electrical.source_work_j[1], supplied);
    near(report.electrical.terminal_resistor_heat_j[1], resistor_heat);
    near(report.electrical.load_work_j[1], delivered);
    near(delivered, 0.3125);
    assert_eq!(report.electrical.integration_steps, 8);
    assert_eq!(report.circuit_trials % 8, 0);
    near(report.deposited_joule_j, resistor_heat + delivered);
    near(exported.report().deposited_joule_j, resistor_heat);
    near(
        heated.thermal().temperature_k(ALL_IDS[0]).unwrap() - 300.0,
        delivered,
    );
    near(
        heated.thermal().temperature_k(ALL_IDS[6]).unwrap() - 300.0,
        resistor_heat,
    );
    near(report.total_energy_residual_j, 0.0);
    near(
        report.electrical.converter.native.delta_h + heated.thermal().energy_j()
            - thermal.energy_j(),
        supplied,
    );
    assert_eq!(heated.converter().state(), exported.converter().state());
    assert_eq!(report.source.load_heat[1].unwrap().cell, ALL_IDS[0]);
    near(report.source.load_heat[1].unwrap().energy_j, delivered);
}

#[test]
fn switching_ripple_and_load_heat_share_actual_temperature_dependent_resistance() {
    let cp = [0.001; 7];
    let alpha = 0.0039;
    let (cards, curves) = sources(all_r(), cp, alpha);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let resistance_key = key();
    let channels = paths(ALL_IDS, &resistance_key);
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.5, 1.0)).unwrap();
    let initial = native.state(1.0, [12.0, 6.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let period = 1.0 / 1024.0;
    let intervals = [
        SwitchInterval {
            state: BridgeSwitchState {
                input_high: true,
                output_high: true,
            },
            duration_s: period / 2.0,
        },
        SwitchInterval {
            state: BridgeSwitchState {
                input_high: false,
                output_high: true,
            },
            duration_s: period / 2.0,
        },
    ];
    let commands = [PowerSwitchCommand {
        source_voltage_v: [14.0, 0.0],
        loads: [None, Some(load(2.0))],
    }; 2];
    let run = |receiver| {
        network
            .step_power_loaded_switching_candidate(
                &thermal,
                period,
                &converter,
                &intervals,
                &commands,
                4,
                [None, receiver],
                power_controls(),
                switch_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .unwrap()
    };
    let exported = run(None);
    let heated = run(Some(ALL_IDS[6]));
    let report = heated.report();
    let theta = heated.thermal().temperature_k(ALL_IDS[6]).unwrap() - 300.0;
    let r = 1.0 + alpha * theta;
    near(heated.converter().terminal_resistance_ohm()[1], r);
    near(
        cp[6] * theta,
        report.electrical.terminal_resistor_heat_j[1] + report.electrical.load_work_j[1],
    );
    near(report.electrical.load_work_j[1], 2.0 * period);
    near(report.total_energy_residual_j, 0.0);
    assert!(
        heated.converter().terminal_resistance_ohm()[1]
            > exported.converter().terminal_resistance_ohm()[1]
    );
    assert!(heated.converter().state()[2] > exported.converter().state()[2] + 1e-9);
    // Heat observes each substep's actual current. Squaring the numerical mean
    // would omit the positive ripple variance even though the ideal source is0V.
    let squared_mean_heat = period * r * report.electrical.mean_source_current_a[1].powi(2);
    assert!(report.electrical.terminal_resistor_heat_j[1] > squared_mean_heat + 1e-10);
    near(report.electrical.source_work_j[1], 0.0);
}

#[test]
fn complete_cycle_clock_work_budget_and_cancellation_preserve_both_histories() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.5, 1.0)).unwrap();
    let initial = native.state(1.0, [12.0, 6.0]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let period = 1.0 / 1024.0;
    let intervals = [
        SwitchInterval {
            state: BridgeSwitchState {
                input_high: true,
                output_high: true,
            },
            duration_s: period / 2.0,
        },
        SwitchInterval {
            state: BridgeSwitchState {
                input_high: false,
                output_high: true,
            },
            duration_s: period / 2.0,
        },
    ];
    let commands = [PowerSwitchCommand {
        source_voltage_v: [14.0, 0.0],
        loads: [None, Some(load(2.0))],
    }; 2];
    let energy = thermal.energy_j();
    let run = |target, solve, receivers, cancelled: &mut dyn FnMut() -> bool| {
        network.step_power_loaded_switching_candidate(
            &thermal,
            target,
            &converter,
            &intervals,
            &commands,
            4,
            receivers,
            power_controls(),
            switch_controls(),
            solve,
            &[0.0; 7],
            cancelled,
        )
    };
    let mut too_small = controls();
    too_small.maximum_circuit_trials = 7; // One complete proposal needs8 native steps.
    assert!(run(period, too_small, [None; 2], &mut || false).is_err());
    assert!(run(2.0 * period, controls(), [None; 2], &mut || false).is_err());
    assert!(
        run(period, controls(), [Some(CellId(9999)), None], &mut || {
            false
        })
        .is_err()
    );
    let calls = Cell::new(0usize);
    let expected = run(period, controls(), [None, Some(ALL_IDS[6])], &mut || {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let stop = calls.get() / 2;
    assert!(stop > 10);
    calls.set(0);
    assert!(
        run(period, controls(), [None, Some(ALL_IDS[6])], &mut || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        })
        .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.energy_j(), energy);
    assert_eq!(converter.state(), &initial);
    let retried = run(period, controls(), [None, Some(ALL_IDS[6])], &mut || false).unwrap();
    assert_eq!(retried.converter().state(), expected.converter().state());
    assert_eq!(retried.thermal().energy_j(), expected.thermal().energy_j());
}

#[test]
fn large_absolute_clock_requires_the_same_represented_cycle_endpoint() {
    let cp = [1.0; 7];
    let (cards, curves) = sources(all_r(), cp, 0.0);
    let context = point(300.0);
    let cells = bindings(&cards, &curves, &context);
    let channels = fixed_paths();
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let start = 1e16;
    let thermal = body(&curves, cp)
        .advance_to(
            start,
            StepConfig {
                max_dt_s: start,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let initial = native.state(0.0, [0.0; 2]).unwrap();
    let converter = network
        .bind(&thermal, &native, [1.0; 2], initial, || false)
        .unwrap();
    let intervals = [SwitchInterval {
        state: BridgeSwitchState {
            input_high: false,
            output_high: false,
        },
        duration_s: 2.0,
    }];
    let commands = [PowerSwitchCommand {
        source_voltage_v: [0.0; 2],
        loads: [None; 2],
    }];
    let run = |target| {
        network.step_power_loaded_switching_candidate(
            &thermal,
            target,
            &converter,
            &intervals,
            &commands,
            1,
            [None; 2],
            power_controls(),
            switch_controls(),
            controls(),
            &[0.0; 7],
            || false,
        )
    };
    // At this epoch adjacent represented clocks differ by two seconds.
    // A loose clock-spacing tolerance must not accept a different endpoint.
    assert!(run(start + 4.0).is_err());
    let accepted = run(start + 2.0).unwrap();
    assert_eq!(
        accepted.thermal().time_s().to_bits(),
        (start + 2.0).to_bits()
    );
    assert_eq!(accepted.converter().state(), &initial);
    assert_eq!(thermal.time_s().to_bits(), start.to_bits());
}
