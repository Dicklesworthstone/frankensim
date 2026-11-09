//! Thevenin supply and regulated demand share one original implicit trajectory.
use super::*;
use fs_phs::circuit::converter::power_load::VoltagePowerLoadedStep;

fn total_balance(step: &VoltagePowerLoadedStep, initial: [f64; 3], next: [f64; 3]) {
    // Independent original L/C energy; no report total is used as the reference.
    let energy = |x: [f64; 3]| {
        x[0] * x[0] / (2.0 * 0.02) + x[1] * x[1] / (2.0 * 2.0) + x[2] * x[2] / (2.0 * 2.0)
    };
    let delta = energy(next) - energy(initial);
    close(step.converter.native.delta_h, delta);
    close(
        delta
            + step.converter.conduction.joules.iter().sum::<f64>()
            + step.terminal_resistor_heat_j.iter().sum::<f64>()
            + step.load_work_j.iter().sum::<f64>(),
        step.source_work_j.iter().sum::<f64>(),
    );
}

#[test]
fn finite_sources_support_both_conversion_directions_and_all_seven_loss_owners() {
    let (a, b, dt) = (0.5, 1.0, 0.03);
    let channel_r = [0.03, 0.01, 0.02, 0.015, 0.0];
    let winding_path_r: f64 = channel_r.iter().sum();
    for (current, input_v) in [(3.0, 24.0), (-3.0, 12.0)] {
        let output_v = (a * input_v - winding_path_r * current) / b;
        let (source_v, source_r, loads, expected_source_i) = if current > 0.0 {
            let resistor_current = output_v / 20.0;
            (
                [input_v + 0.4 * a * current, 0.0],
                [0.4, 20.0],
                [
                    None,
                    Some(load(output_v * (b * current - resistor_current), 10.0, 2.0)),
                ],
                [a * current, -resistor_current],
            )
        } else {
            let resistor_current = input_v / 20.0;
            (
                [0.0, output_v - 0.4 * b * current],
                [20.0, 0.4],
                [
                    Some(load(input_v * (-a * current - resistor_current), 10.0, 2.0)),
                    None,
                ],
                [-resistor_current, -b * current],
            )
        };
        let owner = model(a, b);
        let initial = owner.state(current, [input_v, output_v]).unwrap();
        let mut prepared = owner
            .prepare_power_loaded_controlled(dt, controls(), || false)
            .unwrap();
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let result = prepared
            .step_voltage_driven_controlled(
                &initial,
                &source_v,
                &source_r,
                &loads,
                &mut next,
                &mut voltage,
                || false,
            )
            .unwrap();
        for i in 0..3 {
            close(next[i], initial[i]);
        }
        for p in 0..2 {
            close(result.source_current_a[p], expected_source_i[p]);
            close(
                result.source_work_j[p],
                dt * source_v[p] * expected_source_i[p],
            );
            close(
                result.terminal_resistor_heat_j[p],
                dt * source_r[p] * expected_source_i[p].powi(2),
            );
            close(
                result.converter.terminal_work_j[p],
                dt * voltage[p] * expected_source_i[p] - result.load_work_j[p],
            );
        }
        for p in 0..5 {
            close(
                result.converter.conduction.joules[p],
                dt * channel_r[p] * current.powi(2),
            );
        }
        total_balance(&result, initial, next);
    }
}

#[test]
fn nearly_balanced_source_and_load_preserve_the_actual_net_flow_power_identity() {
    let owner = model(0.0, 0.0);
    let initial = owner.state(0.0, [0.0, 10.0]).unwrap();
    let mut prepared = owner
        .prepare_power_loaded_controlled(0.2, controls(), || false)
        .unwrap();
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let result = prepared
        .step_voltage_driven_controlled(
            &initial,
            &[0.0, 14.8],
            &[1.0, 2.0],
            &[None, Some(load(24.0, 10.0, 2.0))],
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    // Decimal source/branch commands differ in their last-bit arithmetic:
    // ((14.8-10)/2)*10-24 is NOT 10*((14.8-10)/2-24/10).
    // This is a near equilibrium, not grounds to invent or reject finite power.
    close(next[2], initial[2]);
    close(result.source_current_a[1], 2.4);
    close(result.source_work_j[1], 0.2 * 14.8 * 2.4);
    close(
        result.terminal_resistor_heat_j[1],
        0.2 * 2.0 * 2.4_f64.powi(2),
    );
    close(result.load_work_j[1], 0.2 * 24.0);
    total_balance(&result, initial, next);
}

#[test]
fn voltage_sag_and_constant_power_follow_an_independent_capacitor_quadratic() {
    let dt = 0.2;
    let owner = model(0.0, 0.0);
    let initial = owner.state(0.0, [10.0, 10.0]).unwrap();
    let mut prepared = owner
        .prepare_power_loaded_controlled(dt, controls(), || false)
        .unwrap();
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let result = prepared
        .step_voltage_driven_controlled(
            &initial,
            &[10.0, 12.0],
            &[1.0, 2.0],
            &[None, Some(load(24.0, 10.0, 2.0))],
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    // C(V1-V0)/dt=(Vs-Vm)/R-P/Vm, Vm=(V0+V1)/2.
    // The physical root near V0 is the larger root of A Vm²-B Vm+P=0.
    let aa: f64 = 2.0 * 2.0 / dt + 1.0 / 2.0;
    let bb: f64 = 2.0 * 2.0 * 10.0 / dt + 12.0 / 2.0;
    let midpoint = (bb + (bb * bb - 4.0 * aa * 24.0).sqrt()) / (2.0 * aa);
    let source_current = (12.0 - midpoint) / 2.0;
    close(voltage[1], midpoint);
    close(next[2] / 2.0, 2.0 * midpoint - 10.0);
    close(result.source_current_a[1], source_current);
    assert!(source_current > (12.0 - 10.0) / 2.0);
    assert_eq!(result.loads[1].mode, PowerLoadMode::ConstantPower);
    close(
        result.loads[1].incremental_conductance_s,
        -24.0 / midpoint.powi(2),
    );
    close(result.load_work_j[1], dt * 24.0);
    close(
        result.terminal_resistor_heat_j[1],
        dt * 2.0 * source_current.powi(2),
    );
    close(result.source_work_j[1], dt * 12.0 * source_current);
    assert_eq!(result.converter.conduction.joules, [0.0; 5]);
    total_balance(&result, initial, next);
}

#[test]
fn finite_impedance_startup_foldback_and_current_limit_use_actual_accepted_power() {
    let dt = 0.2;
    for sign in [-1.0, 1.0] {
        for (initial_v, source_v, demand, expected_mode, conductance, load_current) in [
            (
                0.0,
                6.0,
                load(20.0, 100.0, 10.0),
                PowerLoadMode::VoltageFoldback,
                0.2,
                0.0,
            ),
            (
                10.0,
                12.0,
                load(100.0, 2.0, 1.0),
                PowerLoadMode::CurrentLimited,
                0.0,
                2.0,
            ),
        ] {
            let owner = model(0.0, 0.0);
            let initial = owner.state(0.0, [0.0, sign * initial_v]).unwrap();
            let mut prepared = owner
                .prepare_power_loaded_controlled(dt, controls(), || false)
                .unwrap();
            let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
            let result = prepared
                .step_voltage_driven_controlled(
                    &initial,
                    &[0.0, sign * source_v],
                    &[1.0, 2.0],
                    &[None, Some(demand)],
                    &mut next,
                    &mut voltage,
                    || false,
                )
                .unwrap();
            let midpoint = sign * (2.0 * 2.0 * initial_v / dt + source_v / 2.0 - load_current)
                / (2.0 * 2.0 / dt + 1.0 / 2.0 + conductance);
            let actual_load_i = conductance * midpoint + sign * load_current;
            close(voltage[1], midpoint);
            close(next[2] / 2.0, 2.0 * midpoint - sign * initial_v);
            close(result.loads[1].current_a, actual_load_i);
            assert_eq!(result.loads[1].mode, expected_mode);
            close(result.load_work_j[1], dt * midpoint * actual_load_i);
            close(
                result.curtailed_load_energy_j[1],
                dt * demand.requested_power_w() - result.load_work_j[1],
            );
            assert!(result.load_work_j[1] < dt * demand.requested_power_w());
            total_balance(&result, initial, next);
        }
    }
}

#[test]
fn a_charged_bus_regenerates_into_its_source_and_heats_only_the_physical_series_resistor() {
    let dt = 0.4;
    for sign in [-1.0, 1.0] {
        let owner = model(0.0, 0.0);
        let initial = owner.state(0.0, [0.0, sign * 12.0]).unwrap();
        let mut prepared = owner
            .prepare_power_loaded_controlled(dt, controls(), || false)
            .unwrap();
        let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
        let result = prepared
            .step_voltage_driven_controlled(
                &initial,
                &[0.0, sign * 4.0],
                &[1.0, 2.0],
                &[None, None],
                &mut next,
                &mut voltage,
                || false,
            )
            .unwrap();
        let midpoint = sign * (2.0 * 2.0 * 12.0 / dt + 4.0 / 2.0) / (2.0 * 2.0 / dt + 1.0 / 2.0);
        let source_current = (sign * 4.0 - midpoint) / 2.0;
        close(voltage[1], midpoint);
        close(next[2] / 2.0, 2.0 * midpoint - sign * 12.0);
        assert!(result.source_work_j[1] < 0.0);
        close(result.source_work_j[1], dt * sign * 4.0 * source_current);
        close(
            result.terminal_resistor_heat_j[1],
            dt * 2.0 * source_current.powi(2),
        );
        assert!((result.terminal_resistor_heat_j[1] - dt * midpoint.powi(2) / 2.0).abs() > 1.0);
        assert_eq!(result.load_work_j, [0.0; 2]);
        total_balance(&result, initial, next);
    }
}

#[test]
fn invalid_sources_final_load_domain_budget_and_late_cancellation_are_transactional() {
    let owner = model(0.0, 0.0);
    let initial = owner.state(0.0, [10.0, 10.0]).unwrap();
    let mut prepared = owner
        .prepare_power_loaded_controlled(0.2, controls(), || false)
        .unwrap();
    let sources = [10.0, 12.0];
    let resistances = [1.0, 2.0];
    let loads = [None, Some(load(24.0, 10.0, 2.0))];
    let (mut expected, mut expected_v) = ([0.0; 3], [0.0; 2]);
    let mut polls = 0;
    prepared
        .step_voltage_driven_controlled(
            &initial,
            &sources,
            &resistances,
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
        prepared.step_voltage_driven_controlled(
            &initial,
            &sources,
            &resistances,
            &loads,
            &mut next,
            &mut voltage,
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    for invalid_r in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-320] {
        assert!(
            prepared
                .step_voltage_driven_controlled(
                    &initial,
                    &sources,
                    &[1.0, invalid_r],
                    &loads,
                    &mut next,
                    &mut voltage,
                    || false,
                )
                .is_err()
        );
        assert_eq!(next, [91.0; 3]);
        assert_eq!(voltage, [73.0; 2]);
    }
    assert!(
        prepared
            .step_voltage_driven_controlled(
                &initial,
                &[10.0, f64::NAN],
                &resistances,
                &loads,
                &mut next,
                &mut voltage,
                || false,
            )
            .is_err()
    );
    // The accepted midpoint would stay below 5 V, but the final bus exceeds
    // the caller's load support. Endpoint applicability still gates publication.
    let bounded = [
        None,
        Some(RegulatedPowerLoad::try_new(1.0, 10.0, 1.0, 5.0).unwrap()),
    ];
    assert!(
        prepared
            .step_voltage_driven_controlled(
                &[0.0, 20.0, 8.0],
                &[10.0, 20.0],
                &[1.0; 2],
                &bounded,
                &mut next,
                &mut voltage,
                || false,
            )
            .is_err()
    );
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    let mut limited = controls();
    limited.solve.max_iterations = 0;
    let mut exhausted = model(0.0, 0.0)
        .prepare_power_loaded_controlled(0.2, limited, || false)
        .unwrap();
    assert!(
        exhausted
            .step_voltage_driven_controlled(
                &initial,
                &sources,
                &resistances,
                &loads,
                &mut next,
                &mut voltage,
                || false,
            )
            .is_err()
    );
    assert_eq!(next, [91.0; 3]);
    assert_eq!(voltage, [73.0; 2]);
    prepared
        .step_voltage_driven_controlled(
            &initial,
            &sources,
            &resistances,
            &loads,
            &mut next,
            &mut voltage,
            || false,
        )
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(voltage.map(f64::to_bits), expected_v.map(f64::to_bits));
}
