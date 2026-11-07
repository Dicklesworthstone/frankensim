//! G1/G3: physical source/load response and heat across commanded PWM edges.
use super::{budget, close, modulation, parameters, resistances};
use fs_phs::PreparedStepError;
use fs_phs::circuit::converter::switching::voltage_driven::PreparedVoltageDrivenSwitchCycle;
use fs_phs::circuit::converter::switching::{
    BridgeSwitchState, SwitchCycleControls, SwitchInterval, centered_pwm,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BuckBoostParameters};
use fs_phs::linear_exact::DenseExponentialBudget;

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 6,
        max_dense_bytes: 16_384,
        max_multiply_terms: 500_000,
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
    intervals: &[SwitchInterval],
    resistance: [f64; 2],
) -> PreparedVoltageDrivenSwitchCycle {
    PreparedVoltageDrivenSwitchCycle::new_controlled(
        parameters(),
        resistance,
        intervals,
        budget(),
        heat_budget(),
        SwitchCycleControls::default(),
        || false,
    )
    .unwrap()
}

// Independent SI Kirchhoff ODE with explicit physical current, five conducting
// heat owners, two source resistors, two ideal-source works and voltage-time
// integrals. No production matrix, state gradient or accounting result enters.
fn rhs(
    x: &[f64; 14],
    p: BuckBoostParameters,
    gate: BridgeSwitchState,
    resistance: [f64; 2],
    source: [f64; 2],
) -> [f64; 14] {
    let a = if gate.input_high { 1.0 } else { 0.0 };
    let b = if gate.output_high { 1.0 } else { 0.0 };
    let r = resistances(p, a, b);
    let current = [
        (source[0] - x[1]) / resistance[0],
        (source[1] - x[2]) / resistance[1],
    ];
    [
        (a * x[1] - b * x[2] - r.iter().sum::<f64>() * x[0]) / p.inductance_h,
        (current[0] - a * x[0]) / p.bus_capacitance_f[0],
        (current[1] + b * x[0]) / p.bus_capacitance_f[1],
        r[0] * x[0] * x[0],
        r[1] * x[0] * x[0],
        r[2] * x[0] * x[0],
        r[3] * x[0] * x[0],
        r[4] * x[0] * x[0],
        resistance[0] * current[0] * current[0],
        resistance[1] * current[1] * current[1],
        source[0] * current[0],
        source[1] * current[1],
        x[1],
        x[2],
    ]
}
fn reference(
    x: &mut [f64; 14],
    p: BuckBoostParameters,
    gate: SwitchInterval,
    resistance: [f64; 2],
    source: [f64; 2],
) {
    let dt = gate.duration_s / 512.0;
    for _ in 0..512 {
        let k1 = rhs(x, p, gate.state, resistance, source);
        let k2 = rhs(
            &core::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]),
            p,
            gate.state,
            resistance,
            source,
        );
        let k3 = rhs(
            &core::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]),
            p,
            gate.state,
            resistance,
            source,
        );
        let k4 = rhs(
            &core::array::from_fn(|i| x[i] + dt * k3[i]),
            p,
            gate.state,
            resistance,
            source,
        );
        for i in 0..14 {
            x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
}

#[test]
fn switched_sources_and_loads_match_independent_current_voltage_and_power_odes() {
    let p = parameters();
    let resistance = [0.2, 0.7];
    let gates = [
        interval(false, true, 0.004),
        interval(true, true, 0.007),
        interval(true, false, 0.003),
        interval(false, false, 0.002),
    ];
    let sources = [[24.0, 0.0], [18.0, 12.0], [0.0, 24.0], [28.0, 0.0]];
    let mut oracle = [
        -1.2, 24.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    let model = AveragedBuckBoost::new(p, modulation(0.0, 0.0)).unwrap();
    let initial = model.state(oracle[0], [oracle[1], oracle[2]]).unwrap();
    let mut owner = prepared(&gates, resistance);
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let record = owner
        .step_interval_sources_controlled(&initial, &sources, &mut next, &mut mean, || false)
        .unwrap();
    for (&gate, &source) in gates.iter().zip(&sources) {
        reference(&mut oracle, p, gate, resistance, source);
    }
    close(next[0] / p.inductance_h, oracle[0], 100.0);
    for port in 0..2 {
        close(
            next[port + 1] / p.bus_capacitance_f[port],
            oracle[port + 1],
            100.0,
        );
        close(mean[port] * owner.period_s(), oracle[12 + port], 1.0);
        close(
            record.terminal_resistor_heat_j[port],
            oracle[8 + port],
            100.0,
        );
        close(record.source_work_j[port], oracle[10 + port], 100.0);
        close(
            record.converter.terminal_work_j[port],
            oracle[10 + port] - oracle[8 + port],
            100.0,
        );
    }
    for channel in 0..5 {
        close(
            record.converter.conduction.joules[channel],
            oracle[3 + channel],
            1.0,
        );
    }
    close(
        record.converter.native.delta_h
            + record.converter.native.dissipated
            + record.terminal_resistor_heat_j.iter().sum::<f64>()
            - record.source_work_j.iter().sum::<f64>(),
        0.0,
        100.0,
    );
    let (mut retry, mut retry_mean) = ([0.0; 3], [0.0; 2]);
    let repeated = owner
        .step_interval_sources_controlled(&initial, &sources, &mut retry, &mut retry_mean, || false)
        .unwrap();
    assert_eq!(next, retry);
    assert_eq!(mean, retry_mean);
    assert_eq!(record.source_work_j, repeated.source_work_j);
}

#[test]
fn static_connection_reuses_native_behavior_and_duplicate_intervals_share_cache() {
    let gates = [interval(true, false, 0.01)];
    let r = [0.2, 6.0];
    let model = AveragedBuckBoost::new(parameters(), modulation(1.0, 0.0)).unwrap();
    let initial = model.state(-1.3, [24.0, 12.0]).unwrap();
    let mut native = model
        .prepare_voltage_driven_controlled(r, 0.01, budget(), heat_budget(), || false)
        .unwrap();
    let mut cycle = prepared(&gates, r);
    let (mut a, mut ya, mut b, mut yb) = ([0.0; 3], [0.0; 2], [0.0; 3], [0.0; 2]);
    let expected = native
        .step_controlled(&initial, &[18.0, 0.0], &mut a, &mut ya, || false)
        .unwrap();
    let actual = cycle
        .step_controlled(&initial, &[18.0, 0.0], &mut b, &mut yb, || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(ya, yb);
    assert_eq!(expected.source_work_j, actual.source_work_j);
    assert_eq!(
        expected.terminal_resistor_heat_j,
        actual.terminal_resistor_heat_j
    );
    assert_eq!(expected.converter.conduction, actual.converter.conduction);
    assert_eq!(expected.converter.native, actual.converter.native);
    let cached = prepared(&[gates[0]; 8], r);
    assert_eq!(cached.cached_interval_count(), 1);
    assert_eq!(cached.parameters(), parameters());
    assert_eq!(cached.terminal_resistance_ohm(), r);
    assert_eq!(cached.intervals().len(), 8);
    assert_eq!(cached.averaged_modulation(), modulation(1.0, 0.0));
}

#[test]
fn loaded_pwm_retains_ripple_heat_and_regeneration() {
    let p = parameters();
    let r = [0.2, 0.7];
    let gates = centered_pwm(0.02, modulation(0.4, 0.75)).unwrap();
    let model = AveragedBuckBoost::new(p, modulation(0.4, 0.75)).unwrap();
    let initial = model.state(-2.0, [28.0, 8.0]).unwrap();
    let mut owner = prepared(&gates, r);
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let source = [24.0, 0.0];
    let record = owner
        .step_controlled(&initial, &source, &mut next, &mut mean, || false)
        .unwrap();
    assert!(
        record.source_work_j[0] < 0.0,
        "overcharged input regenerates into its source"
    );
    assert_eq!(record.source_work_j[1], 0.0, "passive load is not a source");
    assert!(record.converter.terminal_work_j[1] < 0.0);
    for port in 0..2 {
        let mean_current = (source[port] - mean[port]) / r[port];
        let square_of_mean = owner.period_s() * r[port] * mean_current * mean_current;
        assert!(
            record.terminal_resistor_heat_j[port] > square_of_mean + 1e-6,
            "bus ripple must increase RMS resistor heat"
        );
    }
    assert!(record.converter.conduction.joules.iter().all(|j| *j > 0.0));
}

#[test]
fn late_invalid_sources_and_every_cancellation_checkpoint_publish_nothing() {
    let gates = centered_pwm(0.02, modulation(0.4, 0.75)).unwrap();
    let mut owner = prepared(&gates, [0.2, 0.7]);
    let initial = AveragedBuckBoost::new(parameters(), modulation(0.4, 0.75))
        .unwrap()
        .state(-2.0, [28.0, 8.0])
        .unwrap();
    let (mut expected, mut expected_mean) = ([0.0; 3], [0.0; 2]);
    let mut polls = 0;
    owner
        .step_controlled(
            &initial,
            &[24.0, 0.0],
            &mut expected,
            &mut expected_mean,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in 1..=polls {
        let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
        let mut calls = 0;
        let result = owner.step_controlled(&initial, &[24.0, 0.0], &mut next, &mut mean, || {
            calls += 1;
            calls == stop
        });
        assert!(matches!(result, Err(PreparedStepError::Cancelled)));
        assert_eq!(next, [123.0; 3]);
        assert_eq!(mean, [456.0; 2]);
    }
    let mut sources = vec![[24.0, 0.0]; gates.len()];
    sources.last_mut().unwrap()[0] = f64::NAN;
    let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
    assert!(
        owner
            .step_interval_sources_controlled(&initial, &sources, &mut next, &mut mean, || false,)
            .is_err()
    );
    assert_eq!(next, [123.0; 3]);
    assert_eq!(mean, [456.0; 2]);
    assert!(
        owner
            .step_interval_sources_controlled(
                &initial,
                &sources[..sources.len() - 1],
                &mut next,
                &mut mean,
                || false,
            )
            .is_err()
    );
    owner
        .step_controlled(&initial, &[24.0, 0.0], &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!(next, expected);
    assert_eq!(mean, expected_mean);
}

#[test]
fn voltage_cycle_refuses_invalid_physical_parameters_and_event_limits() {
    let gates = [interval(true, false, 0.01)];
    for r in [[0.0, 1.0], [-1.0, 1.0], [f64::NAN, 1.0]] {
        assert!(
            PreparedVoltageDrivenSwitchCycle::new_controlled(
                parameters(),
                r,
                &gates,
                budget(),
                heat_budget(),
                SwitchCycleControls::default(),
                || false,
            )
            .is_err()
        );
    }
    for invalid in [vec![], vec![gates[0]; 65], vec![interval(true, false, 0.0)]] {
        assert!(
            PreparedVoltageDrivenSwitchCycle::new_controlled(
                parameters(),
                [0.2, 6.0],
                &invalid,
                budget(),
                heat_budget(),
                SwitchCycleControls::default(),
                || false,
            )
            .is_err()
        );
    }
}
