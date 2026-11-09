//! Regulated demand, independent capacitor energy and converter heat ownership.
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadMode, RegulatedPowerLoad,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::{PreparedStepError, SolveBudget};

#[path = "circuit_converter_power_load/voltage.rs"]
mod voltage;

fn parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [2.0, 2.0],
        winding_resistance_ohm: 0.03,
        input_on_resistance_ohm: [0.02, 0.04],
        output_on_resistance_ohm: [0.015, 0.025],
    }
}
fn model(a: f64, b: f64) -> AveragedBuckBoost {
    AveragedBuckBoost::new(
        parameters(),
        BridgeModulation {
            input_high_fraction: a,
            output_high_fraction: b,
        },
    )
    .unwrap()
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
fn load(power: f64, current: f64, foldback: f64) -> RegulatedPowerLoad {
    RegulatedPowerLoad::try_new(power, current, foldback, 100.0).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-10, "{a:.17e} != {b:.17e}");
}

#[test]
fn load_law_is_continuous_bipolar_and_has_the_actual_incremental_tangent() {
    let demand = load(100.0, 5.0, 4.0);
    for sign in [-1.0, 1.0] {
        let low = demand.evaluate(sign * 0.5).unwrap();
        assert_eq!(low.mode, PowerLoadMode::VoltageFoldback);
        close(low.current_a, sign * 3.125);
        close(low.delivered_power_w, 1.5625);
        close(low.incremental_conductance_s, 6.25);
        let capped = demand.evaluate(sign * 10.0).unwrap();
        assert_eq!(capped.mode, PowerLoadMode::CurrentLimited);
        close(capped.current_a, sign * 5.0);
        close(capped.delivered_power_w, 50.0);
        close(capped.curtailed_power_w, 50.0);
        close(capped.incremental_conductance_s, 0.0);
        let regulated = demand.evaluate(sign * 40.0).unwrap();
        assert_eq!(regulated.mode, PowerLoadMode::ConstantPower);
        close(regulated.current_a, sign * 2.5);
        close(regulated.delivered_power_w, 100.0);
        close(regulated.incremental_conductance_s, -0.0625);
    }
    for voltage in [-40.0, -10.0, -0.5, 0.0, 0.5, 10.0, 40.0] {
        let derivative = (demand.evaluate(voltage + 1e-5).unwrap().current_a
            - demand.evaluate(voltage - 1e-5).unwrap().current_a)
            / (2e-5);
        close(
            derivative,
            demand.evaluate(voltage).unwrap().incremental_conductance_s,
        );
    }
    let no_cap = load(100.0, 100.0, 4.0);
    close(
        no_cap.evaluate(4.0 - 1e-12).unwrap().current_a,
        no_cap.evaluate(4.0 + 1e-12).unwrap().current_a,
    );
    close(
        demand.evaluate(20.0 - 1e-12).unwrap().current_a,
        demand.evaluate(20.0 + 1e-12).unwrap().current_a,
    );
    close(demand.evaluate(0.0).unwrap().delivered_power_w, 0.0);
    assert_eq!(
        load(0.0, 1.0, 1.0).evaluate(20.0).unwrap().mode,
        PowerLoadMode::Disabled
    );
    assert!(demand.evaluate(101.0).is_err());
    assert!(demand.evaluate(f64::NAN).is_err());
    for (p, i, vf, vm) in [
        (f64::NAN, 1.0, 1.0, 10.0),
        (-1.0, 1.0, 1.0, 10.0),
        (1.0, 0.0, 1.0, 10.0),
        (1.0, 1.0, 0.0, 10.0),
        (1.0, 1.0, 20.0, 10.0),
    ] {
        assert!(RegulatedPowerLoad::try_new(p, i, vf, vm).is_err());
    }
}

#[test]
fn capacitor_supplies_constant_requested_energy_without_converter_heating() {
    let owner = model(0.0, 0.0);
    let mut state = owner.state(0.0, [10.0, 12.0]).unwrap();
    let mut prepared = owner
        .prepare_power_loaded_controlled(0.1, controls(), || false)
        .unwrap();
    let demand = [None, Some(load(24.0, 10.0, 2.0))];
    for k in 1..=5 {
        let before = state;
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let result = prepared
            .step_controlled(&state, &[0.0; 2], &demand, &mut next, &mut voltage, || {
                false
            })
            .unwrap();
        // C=2: V(t)^2=144-24*t. This independent capacitor energy equation
        // fixes the positive continuous solution as well as its discrete endpoint.
        let expected = (144.0 - 2.4 * f64::from(k)).sqrt();
        close(next[2] / 2.0, expected);
        close(voltage[1], 0.5 * (before[2] + next[2]) / 2.0);
        close(result.load_work_j[1], 2.4);
        close(result.loads[1].current_a, 24.0 / voltage[1]);
        assert_eq!(result.loads[1].mode, PowerLoadMode::ConstantPower);
        close(result.converter.native.delta_h, -2.4);
        close(result.converter.native.supplied, -2.4);
        close(result.converter.native.dissipated, 0.0);
        assert_eq!(result.converter.conduction.joules, [0.0; 5]);
        assert_eq!(result.source_work_j, [0.0; 2]);
        close(result.converter.terminal_work_j[1], -2.4);
        state = next;
    }
}

#[test]
fn low_voltage_foldback_and_current_limits_follow_independent_capacitor_equations() {
    for sign in [-1.0, 1.0] {
        let owner = model(0.0, 0.0);
        let initial = owner.state(0.0, [0.0, sign]).unwrap();
        let mut prepared = owner
            .prepare_power_loaded_controlled(0.2, controls(), || false)
            .unwrap();
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        let r = prepared
            .step_controlled(
                &initial,
                &[0.0; 2],
                &[None, Some(load(20.0, 100.0, 10.0))],
                &mut next,
                &mut mean,
                || false,
            )
            .unwrap();
        // C=2 and foldback conductance=0.2 S, implicit midpoint RC decay.
        let expected = sign * (2.0 - 0.2 * 0.2 / 2.0) / (2.0 + 0.2 * 0.2 / 2.0);
        close(next[2] / 2.0, expected);
        assert_eq!(r.loads[1].mode, PowerLoadMode::VoltageFoldback);
        close(r.load_work_j[1], 1.0 - expected * expected);
        close(r.load_work_j[1], 0.2 * 0.2 * mean[1] * mean[1]);
        close(r.curtailed_load_energy_j[1], 4.0 - r.load_work_j[1]);

        let owner = model(0.0, 0.0);
        let initial = owner.state(0.0, [0.0, 10.0 * sign]).unwrap();
        let mut prepared = owner
            .prepare_power_loaded_controlled(0.1, controls(), || false)
            .unwrap();
        let r = prepared
            .step_controlled(
                &initial,
                &[0.0; 2],
                &[None, Some(load(100.0, 2.0, 1.0))],
                &mut next,
                &mut mean,
                || false,
            )
            .unwrap();
        close(next[2] / 2.0, 9.9 * sign);
        close(r.loads[1].current_a, 2.0 * sign);
        close(r.load_work_j[1], 1.99);
        close(r.curtailed_load_energy_j[1], 8.01);
        assert_eq!(r.loads[1].mode, PowerLoadMode::CurrentLimited);
    }
}

#[test]
fn winding_heat_and_exported_power_are_separate_on_the_same_transient() {
    let owner = model(0.0, 0.0);
    let initial = owner.state(1.0, [10.0, 12.0]).unwrap();
    let dt = 0.05;
    let mut prepared = owner
        .prepare_power_loaded_controlled(dt, controls(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let r = prepared
        .step_controlled(
            &initial,
            &[0.0; 2],
            &[None, Some(load(24.0, 10.0, 2.0))],
            &mut next,
            &mut mean,
            || false,
        )
        .unwrap();
    let resistance = 0.03 + 0.04 + 0.025;
    let current = (2.0 * 0.02 - dt * resistance) / (2.0 * 0.02 + dt * resistance);
    let heat = 0.5 * 0.02 * (1.0 - current * current);
    close(next[0] / 0.02, current);
    close(r.load_work_j[1], 1.2);
    close(r.converter.native.dissipated, heat);
    close(r.converter.native.delta_h, -heat - 1.2);
    for (actual, resistance) in r
        .converter
        .conduction
        .joules
        .into_iter()
        .zip([0.03, 0.0, 0.04, 0.0, 0.025])
    {
        close(actual, dt * resistance * (0.5 + 0.5 * current).powi(2));
    }
}

#[test]
fn both_conversion_directions_close_real_source_work_load_work_and_five_losses() {
    let (a, b, dt) = (0.5, 1.0, 0.03);
    let channel_r = [0.03, 0.01, 0.02, 0.015, 0.0];
    let resistance: f64 = channel_r.iter().sum();
    for (current, input_v) in [(3.0, 24.0), (-3.0, 12.0)] {
        let output_v = (a * input_v - resistance * current) / b;
        let (sources, loads) = if current > 0.0 {
            (
                [a * current, 0.0],
                [None, Some(load(output_v * b * current, 10.0, 2.0))],
            )
        } else {
            (
                [0.0, -b * current],
                [Some(load(-input_v * a * current, 10.0, 2.0)), None],
            )
        };
        let owner = model(a, b);
        let initial = owner.state(current, [input_v, output_v]).unwrap();
        let mut prepared = owner
            .prepare_power_loaded_controlled(dt, controls(), || false)
            .unwrap();
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        let r = prepared
            .step_controlled(&initial, &sources, &loads, &mut next, &mut mean, || false)
            .unwrap();
        for i in 0..3 {
            close(next[i], initial[i]);
        }
        for p in 0..2 {
            close(r.source_work_j[p], dt * mean[p] * sources[p]);
        }
        for p in 0..5 {
            close(
                r.converter.conduction.joules[p],
                dt * channel_r[p] * current * current,
            );
        }
        close(r.converter.native.delta_h, 0.0);
        close(
            r.source_work_j.iter().sum::<f64>() - r.load_work_j.iter().sum::<f64>(),
            dt * resistance * current * current,
        );
    }
}

#[test]
fn late_cancellation_source_refusal_and_solver_budget_leave_outputs_unchanged() {
    let owner = model(0.0, 0.0);
    let initial = owner.state(0.0, [10.0, 12.0]).unwrap();
    let mut prepared = owner
        .prepare_power_loaded_controlled(0.1, controls(), || false)
        .unwrap();
    let loads = [None, Some(load(24.0, 10.0, 2.0))];
    let (mut expected, mut expected_v) = ([0.0; 3], [0.0; 2]);
    let mut polls = 0;
    prepared
        .step_controlled(
            &initial,
            &[0.0; 2],
            &loads,
            &mut expected,
            &mut expected_v,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut next, mut voltage) = ([91.0; 3], [73.0; 2]);
    let mut seen = 0;
    assert!(matches!(
        prepared.step_controlled(&initial, &[0.0; 2], &loads, &mut next, &mut voltage, || {
            seen += 1;
            seen == polls
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    let invalid = [
        None,
        Some(RegulatedPowerLoad::try_new(24.0, 10.0, 2.0, 5.0).unwrap()),
    ];
    assert!(
        prepared
            .step_controlled(
                &initial,
                &[0.0; 2],
                &invalid,
                &mut next,
                &mut voltage,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    let mut limited = controls();
    limited.solve.max_iterations = 0;
    let mut no_updates = model(0.0, 0.0)
        .prepare_power_loaded_controlled(0.1, limited, || false)
        .unwrap();
    assert!(
        no_updates
            .step_controlled(&initial, &[0.0; 2], &loads, &mut next, &mut voltage, || {
                false
            })
            .is_err()
    );
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    prepared
        .step_controlled(&initial, &[0.0; 2], &loads, &mut next, &mut voltage, || {
            false
        })
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(voltage.map(f64::to_bits), expected_v.map(f64::to_bits));
}
