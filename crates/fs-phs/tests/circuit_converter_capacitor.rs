//! Finite terminal storage, independent resistor work, and original converter dynamics.
use fs_phs::PreparedStepError;
use fs_phs::circuit::converter::capacitor_terminals::{CapacitorFedBuckBoost, CapacitorTerminals};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};

fn parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 1.0,
        bus_capacitance_f: [2.0, 1.0],
        winding_resistance_ohm: 0.3,
        input_on_resistance_ohm: [0.1, 0.2],
        output_on_resistance_ohm: [0.15, 0.25],
    }
}

fn terminals() -> CapacitorTerminals {
    CapacitorTerminals {
        capacitance_f: [3.0, 4.0],
        resistance_ohm: [4.0, 4.0],
    }
}

fn model(
    p: BuckBoostParameters,
    terminal: CapacitorTerminals,
    a: f64,
    b: f64,
) -> CapacitorFedBuckBoost {
    AveragedBuckBoost::new(
        p,
        BridgeModulation {
            input_high_fraction: a,
            output_high_fraction: b,
        },
    )
    .unwrap()
    .with_capacitor_terminals(terminal)
    .unwrap()
}

fn linear_budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 5,
        max_omega_dt: 0.04,
    }
}

fn observation_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        // Five states, two held currents and a constant require eight flow
        // coordinates; positive source-square observers require ten.
        max_dimension: 10,
        max_dense_bytes: 1_000_000,
        max_multiply_terms: 100_000_000,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "{actual:.17e} != {expected:.17e} (scale {scale:e})"
    );
}

#[test]
fn two_capacitor_charge_sharing_has_independent_forward_and_reverse_work() {
    let p = parameters();
    let terminal = terminals();
    let duration = 4.8 * 2.0_f64.ln();
    for (bus0, source0, bus1, source1, bus_work, cap_work) in [
        (2.0, 12.0, 5.0, 10.0, 21.0, 66.0),
        (12.0, 2.0, 9.0, 4.0, -63.0, -18.0),
    ] {
        let owner = model(p, terminal, 0.0, 0.0);
        let initial = owner.state(0.0, [bus0, 4.0], [source0, 4.0]).unwrap();
        let mut prepared = owner
            .prepare_controlled(duration, linear_budget(), observation_budget(), || false)
            .unwrap();
        let (mut next, mut mean) = ([0.0; 5], [0.0; 4]);
        let record = prepared
            .step_controlled(&initial, &[0.0; 2], &mut next, &mut mean, || false)
            .unwrap();
        close(next[0], 0.0, 1.0);
        close(next[1] / p.bus_capacitance_f[0], bus1, 12.0);
        close(next[3] / terminal.capacitance_f[0], source1, 12.0);
        close(next[2], initial[2], 20.0);
        close(next[4], initial[4], 20.0);
        // Total charge is unchanged. Ceq=1.2 F, tau=4.8 s and exp(-t/tau)=1/2.
        close(next[1] + next[3], initial[1] + initial[3], 40.0);
        let weighted_voltage = (2.0 * bus0 + 3.0 * source0) / 5.0;
        let mean_difference = (source0 - bus0) * 2.4 / duration;
        close(mean[0], weighted_voltage - 0.6 * mean_difference, 12.0);
        close(mean[2], weighted_voltage + 0.4 * mean_difference, 12.0);
        close(mean[1], 4.0, 4.0);
        close(mean[3], 4.0, 4.0);
        close(record.converter.terminal_work_j[0], bus_work, 100.0);
        close(record.terminal_outflow_work_j[0], cap_work, 100.0);
        close(record.terminal_resistor_heat_j[0], 45.0, 100.0);
        close(record.converter.native.delta_h, bus_work, 100.0);
        close(record.terminal_storage_change_j[0], -cap_work, 100.0);
        assert_eq!(record.external_work_j, [0.0; 2]);
        assert_eq!(record.converter.conduction.joules, [0.0; 5]);
        close(record.terminal_resistor_heat_j[1], 0.0, 1.0);
        close(record.converter.terminal_work_j[1], 0.0, 1.0);
    }
}

#[test]
fn held_charger_and_load_currents_use_actual_capacitor_voltage_integrals() {
    let p = BuckBoostParameters {
        bus_capacitance_f: [1.0, 2.0],
        ..parameters()
    };
    let terminal = CapacitorTerminals {
        capacitance_f: [1.0, 4.0],
        resistance_ohm: [2.0, 4.0],
    };
    let duration = 2.0_f64.ln();
    for current in [-2.0, 2.0] {
        let owner = model(p, terminal, 0.0, 0.0);
        let initial = owner.state(0.0, [3.0, 4.0], [3.0, 4.0]).unwrap();
        let mut prepared = owner
            .prepare_controlled(duration, linear_budget(), observation_budget(), || false)
            .unwrap();
        let (mut next, mut mean) = ([0.0; 5], [0.0; 4]);
        let record = prepared
            .step_controlled(&initial, &[current, 0.0], &mut next, &mut mean, || false)
            .unwrap();
        // Cb=Cs=1, R=2: tau=1, dss=I, d(t)=I*(1-exp(-t)).
        let expected_bus = 3.0 + 0.5 * current * duration - 0.25 * current;
        let expected_cap = 3.0 + 0.5 * current * duration + 0.25 * current;
        let integrated_cap_voltage =
            3.0 * duration + current * duration * duration / 4.0 + current * (duration - 0.5) / 2.0;
        let expected_work = current * integrated_cap_voltage;
        let expected_heat = current * current / 2.0 * (duration - 0.625);
        close(next[1], expected_bus, 5.0);
        close(next[3], expected_cap, 5.0);
        close(mean[2], integrated_cap_voltage / duration, 5.0);
        close(record.external_work_j[0], expected_work, 10.0);
        close(record.terminal_resistor_heat_j[0], expected_heat, 1.0);
        close(next[1] + next[3], 6.0 + current * duration, 10.0);
        close(
            record.converter.terminal_work_j[0],
            0.5 * (expected_bus * expected_bus - 9.0),
            10.0,
        );
        close(
            record.terminal_storage_change_j[0],
            0.5 * (expected_cap * expected_cap - 9.0),
            10.0,
        );
        assert_eq!(record.external_work_j[0].is_sign_positive(), current > 0.0);
        assert!((expected_work - current * expected_cap * duration).abs() > 0.1);
    }
}

#[test]
fn low_side_freewheel_retains_copper_and_each_conducting_device_heat() {
    let p = parameters();
    let owner = model(
        p,
        CapacitorTerminals {
            capacitance_f: [4.0; 2],
            resistance_ohm: [4.0; 2],
        },
        0.0,
        0.0,
    );
    let current = -2.0;
    let duration = 0.3;
    let initial = owner.state(current, [4.0, 4.0], [4.0, 4.0]).unwrap();
    let total_r =
        p.winding_resistance_ohm + p.input_on_resistance_ohm[1] + p.output_on_resistance_ohm[1];
    let decay = (-total_r * duration / p.inductance_h).exp();
    let integrated_current_square =
        current * current * p.inductance_h / (2.0 * total_r) * (1.0 - decay * decay);
    let mut prepared = owner
        .prepare_controlled(duration, linear_budget(), observation_budget(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 5], [0.0; 4]);
    let record = prepared
        .step_controlled(&initial, &[0.0; 2], &mut next, &mut mean, || false)
        .unwrap();
    close(next[0] / p.inductance_h, current * decay, 2.0);
    for (heat, resistance) in record.converter.conduction.joules.iter().zip([
        p.winding_resistance_ohm,
        0.0,
        p.input_on_resistance_ohm[1],
        0.0,
        p.output_on_resistance_ohm[1],
    ]) {
        close(*heat, resistance * integrated_current_square, 2.0);
    }
    close(
        record.converter.native.delta_h,
        -total_r * integrated_current_square,
        2.0,
    );
    close(
        record.converter.native.dissipated,
        total_r * integrated_current_square,
        2.0,
    );
    for heat in record.terminal_resistor_heat_j {
        close(heat, 0.0, 1.0);
    }
    assert_eq!(record.external_work_j, [0.0; 2]);
}

fn rk4<const N: usize>(state: &mut [f64; N], dt: f64, rhs: impl Fn(&[f64; N]) -> [f64; N]) {
    let k1 = rhs(state);
    let k2 = rhs(&core::array::from_fn(|i| state[i] + 0.5 * dt * k1[i]));
    let k3 = rhs(&core::array::from_fn(|i| state[i] + 0.5 * dt * k2[i]));
    let k4 = rhs(&core::array::from_fn(|i| state[i] + dt * k3[i]));
    for i in 0..N {
        state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn active_reverse_conversion_matches_independent_current_voltage_and_power_odes() {
    let p = BuckBoostParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        ..parameters()
    };
    let terminal = CapacitorTerminals {
        capacitance_f: [0.5, 0.7],
        resistance_ohm: [0.4, 0.6],
    };
    let (a, b) = (0.6, 0.8);
    let external = [0.3, -0.2];
    let duration = 0.04;
    let owner = model(p, terminal, a, b);
    let initial = owner.state(-1.0, [12.0, 24.0], [8.0, 30.0]).unwrap();
    let mut prepared = owner
        .prepare_controlled(duration, linear_budget(), observation_budget(), || false)
        .unwrap();
    let resistance = [
        p.winding_resistance_ohm,
        a * p.input_on_resistance_ohm[0],
        (1.0 - a) * p.input_on_resistance_ohm[1],
        b * p.output_on_resistance_ohm[0],
        (1.0 - b) * p.output_on_resistance_ohm[1],
    ];
    let total_r: f64 = resistance.iter().sum();
    // Physical i,Vb0,Vb1,Vs0,Vs1, followed by independent bus work,
    // capacitor-outflow work, charger work, ESR heat, converter heat, and
    // four voltage-time integrals. No production matrices supply this oracle.
    let mut oracle = [0.0; 18];
    oracle[..5].copy_from_slice(&[-1.0, 12.0, 24.0, 8.0, 30.0]);
    for _ in 0..512 {
        rk4(&mut oracle, duration / 512.0, |x| {
            let flow = [
                (x[3] - x[1]) / terminal.resistance_ohm[0],
                (x[4] - x[2]) / terminal.resistance_ohm[1],
            ];
            [
                (a * x[1] - b * x[2] - total_r * x[0]) / p.inductance_h,
                (flow[0] - a * x[0]) / p.bus_capacitance_f[0],
                (flow[1] + b * x[0]) / p.bus_capacitance_f[1],
                (external[0] - flow[0]) / terminal.capacitance_f[0],
                (external[1] - flow[1]) / terminal.capacitance_f[1],
                x[1] * flow[0],
                x[2] * flow[1],
                x[3] * flow[0],
                x[4] * flow[1],
                x[3] * external[0],
                x[4] * external[1],
                terminal.resistance_ohm[0] * flow[0] * flow[0],
                terminal.resistance_ohm[1] * flow[1] * flow[1],
                total_r * x[0] * x[0],
                x[1],
                x[2],
                x[3],
                x[4],
            ]
        });
    }
    let (mut next, mut mean) = ([0.0; 5], [0.0; 4]);
    let record = prepared
        .step_controlled(&initial, &external, &mut next, &mut mean, || false)
        .unwrap();
    for (index, coefficient) in [
        p.inductance_h,
        p.bus_capacitance_f[0],
        p.bus_capacitance_f[1],
        terminal.capacitance_f[0],
        terminal.capacitance_f[1],
    ]
    .into_iter()
    .enumerate()
    {
        close(next[index] / coefficient, oracle[index], 40.0);
    }
    for port in 0..2 {
        close(
            record.converter.terminal_work_j[port],
            oracle[5 + port],
            30.0,
        );
        close(record.terminal_outflow_work_j[port], oracle[7 + port], 30.0);
        close(record.external_work_j[port], oracle[9 + port], 1.0);
        close(
            record.terminal_resistor_heat_j[port],
            oracle[11 + port],
            3.0,
        );
    }
    for index in 0..4 {
        close(mean[index], oracle[14 + index] / duration, 40.0);
    }
    close(record.converter.native.dissipated, oracle[13], 3.0);
    for (heat, r) in record.converter.conduction.joules.iter().zip(resistance) {
        close(*heat, oracle[13] * r / total_r, 3.0);
    }
    close(
        record.total.delta_h,
        oracle[9] + oracle[10] - oracle[11] - oracle[12] - oracle[13],
        10.0,
    );
    assert!(next[3] > initial[3]);
    assert!(next[4] < initial[4]);
    assert!(record.converter.terminal_work_j[0] < 0.0);
    assert!(record.converter.terminal_work_j[1] > 0.0);
}

#[test]
fn exact_equilibrium_preserves_a_large_unchanged_energy_datum() {
    let p = BuckBoostParameters {
        bus_capacitance_f: [2.0, 4.0],
        ..parameters()
    };
    let terminal = CapacitorTerminals {
        capacitance_f: [8.0, 16.0],
        resistance_ohm: [1.0, 4.0],
    };
    let voltage = [2.0_f64.powi(80), -2.0_f64.powi(60)];
    let owner = model(p, terminal, 0.0, 0.0);
    let initial = owner.state(0.0, voltage, voltage).unwrap();
    let mut prepared = owner
        .prepare_controlled(0.2, linear_budget(), observation_budget(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 5], [0.0; 4]);
    let record = prepared
        .step_controlled(&initial, &[0.0; 2], &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), initial.map(f64::to_bits));
    assert_eq!(mean, [voltage[0], voltage[1], voltage[0], voltage[1]]);
    assert_eq!(record.converter.native.delta_h, 0.0);
    assert_eq!(record.converter.terminal_work_j, [0.0; 2]);
    assert_eq!(record.terminal_outflow_work_j, [0.0; 2]);
    assert_eq!(record.terminal_resistor_heat_j, [0.0; 2]);
}

#[test]
fn late_cancellation_and_invalid_sources_leave_both_outputs_unchanged() {
    let owner = model(parameters(), terminals(), 0.6, 0.8);
    let initial = owner.state(1.0, [12.0, 8.0], [15.0, 7.0]).unwrap();
    let mut prepared = owner
        .prepare_controlled(0.02, linear_budget(), observation_budget(), || false)
        .unwrap();
    let (mut expected, mut expected_mean) = ([0.0; 5], [0.0; 4]);
    let mut polls = 0;
    let expected_record = prepared
        .step_controlled(
            &initial,
            &[0.3, -0.2],
            &mut expected,
            &mut expected_mean,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 1);
    let mut seen = 0;
    let (mut next, mut mean) = ([123.0; 5], [456.0; 4]);
    assert!(matches!(
        prepared.step_controlled(&initial, &[0.3, -0.2], &mut next, &mut mean, || {
            seen += 1;
            seen == polls
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [123.0; 5]);
    assert_eq!(mean, [456.0; 4]);
    assert!(
        prepared
            .step_controlled(&initial, &[f64::NAN, 0.0], &mut next, &mut mean, || false)
            .is_err()
    );
    assert_eq!(next, [123.0; 5]);
    assert_eq!(mean, [456.0; 4]);
    let retried = prepared
        .step_controlled(&initial, &[0.3, -0.2], &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(mean.map(f64::to_bits), expected_mean.map(f64::to_bits));
    assert_eq!(
        retried.converter.terminal_work_j.map(f64::to_bits),
        expected_record.converter.terminal_work_j.map(f64::to_bits)
    );
}

#[test]
fn physical_terminal_data_and_complete_observation_budgets_are_required() {
    let make_base = || {
        AveragedBuckBoost::new(
            parameters(),
            BridgeModulation {
                input_high_fraction: 0.0,
                output_high_fraction: 0.0,
            },
        )
        .unwrap()
    };
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN, f64::from_bits(1)] {
        let mut terminal = terminals();
        terminal.capacitance_f[0] = value;
        assert!(make_base().with_capacitor_terminals(terminal).is_err());
        terminal = terminals();
        terminal.resistance_ohm[1] = value;
        assert!(make_base().with_capacitor_terminals(terminal).is_err());
    }
    let make_owner = || model(parameters(), terminals(), 0.0, 0.0);
    let mut linear = linear_budget();
    linear.max_modes = 4;
    assert!(
        make_owner()
            .prepare_controlled(0.1, linear, observation_budget(), || false)
            .is_err()
    );
    for dense in [
        DenseExponentialBudget {
            max_dimension: 9,
            ..observation_budget()
        },
        DenseExponentialBudget {
            max_dense_bytes: 1,
            ..observation_budget()
        },
        DenseExponentialBudget {
            max_multiply_terms: 1,
            ..observation_budget()
        },
    ] {
        assert!(
            make_owner()
                .prepare_controlled(0.1, linear_budget(), dense, || false)
                .is_err()
        );
    }
    assert!(matches!(
        make_owner().prepare_controlled(0.1, linear_budget(), observation_budget(), || true),
        Err(PreparedStepError::Cancelled)
    ));
}
