//! G1/G2/G3: original gate intervals resolve ripple, channel heat and port work.
use super::*;
use fs_phs::circuit::converter::switching::{
    BridgeSwitchState, MAX_SWITCH_INTERVALS, PreparedSwitchCycle, SwitchCycleControls,
    SwitchInterval, centered_pwm,
};

fn interval(a: bool, b: bool, duration_s: f64) -> SwitchInterval {
    SwitchInterval {
        state: BridgeSwitchState {
            input_high: a,
            output_high: b,
        },
        duration_s,
    }
}

fn prepared(p: BuckBoostParameters, intervals: &[SwitchInterval]) -> PreparedSwitchCycle {
    PreparedSwitchCycle::new_controlled(
        p,
        intervals,
        budget(),
        SwitchCycleControls::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn centered_edges_retain_both_leg_duties_and_simultaneous_or_static_connections() {
    let cycle = centered_pwm(0.008, modulation(0.25, 0.75)).unwrap();
    assert_eq!(cycle.len(), 5);
    assert_eq!(
        cycle.iter().map(|step| step.state).collect::<Vec<_>>(),
        vec![
            interval(false, false, 1.0).state,
            interval(false, true, 1.0).state,
            interval(true, true, 1.0).state,
            interval(false, true, 1.0).state,
            interval(false, false, 1.0).state
        ]
    );
    for (step, duration) in cycle.iter().zip([0.001, 0.002, 0.002, 0.002, 0.001]) {
        close(step.duration_s, duration, 0.01);
    }
    let owner = prepared(parameters(), &cycle);
    close(owner.period_s(), 0.008, 0.01);
    close(owner.averaged_modulation().input_high_fraction, 0.25, 1.0);
    close(owner.averaged_modulation().output_high_fraction, 0.75, 1.0);
    assert_eq!(centered_pwm(0.008, modulation(0.5, 0.5)).unwrap().len(), 3);
    for (a, b) in [(false, false), (false, true), (true, false), (true, true)] {
        let m = modulation(if a { 1.0 } else { 0.0 }, if b { 1.0 } else { 0.0 });
        assert_eq!(centered_pwm(0.008, m).unwrap(), vec![interval(a, b, 0.008)]);
    }
    // A requested nonzero pulse must not silently disappear when its two
    // command times round onto the same representable clock value.
    assert!(centered_pwm(1.0, modulation(f64::MIN_POSITIVE, 0.5)).is_err());
    assert!(centered_pwm(1.0, modulation(f64::from_bits(1.0_f64.to_bits() - 1), 0.5)).is_err());
}

#[test]
fn static_gate_cycle_is_the_same_original_propagator_and_repeated_intervals_share_cache() {
    let p = parameters();
    let model = AveragedBuckBoost::new(p, modulation(1.0, 0.0)).unwrap();
    let state = model.state(-1.3, [24.0, 12.0]).unwrap();
    let mut native = model.prepare_controlled(0.01, budget(), || false).unwrap();
    let mut cycle = prepared(p, &[interval(true, false, 0.01)]);
    let (mut a, mut ya, mut b, mut yb) = ([0.0; 3], [0.0; 2], [0.0; 3], [0.0; 2]);
    let reference = native
        .step_controlled(&state, &[0.4, -0.3], &mut a, &mut ya, || false)
        .unwrap();
    let actual = cycle
        .step_controlled(&state, &[0.4, -0.3], &mut b, &mut yb, || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(ya, yb);
    assert_eq!(reference.conduction, actual.conduction);
    assert_eq!(reference.terminal_work_j, actual.terminal_work_j);
    assert_eq!(reference.native, actual.native);
    let intervals = [interval(true, false, 0.001); 8];
    let cached = prepared(p, &intervals);
    assert_eq!(cached.intervals(), &intervals[..]);
    assert_eq!(cached.cached_interval_count(), 1);
}

// Independent physical-current/voltage integration with five separate channel
// heat states and two actual terminal-work states. It never reads production
// matrices, gradients, means or heat records.
fn rhs(
    x: &[f64; 10],
    p: BuckBoostParameters,
    gate: BridgeSwitchState,
    current: [f64; 2],
) -> [f64; 10] {
    let a = if gate.input_high { 1.0 } else { 0.0 };
    let b = if gate.output_high { 1.0 } else { 0.0 };
    let rs = resistances(p, a, b);
    let r: f64 = rs.iter().sum();
    [
        (a * x[1] - b * x[2] - r * x[0]) / p.inductance_h,
        (current[0] - a * x[0]) / p.bus_capacitance_f[0],
        (current[1] + b * x[0]) / p.bus_capacitance_f[1],
        rs[0] * x[0] * x[0],
        rs[1] * x[0] * x[0],
        rs[2] * x[0] * x[0],
        rs[3] * x[0] * x[0],
        rs[4] * x[0] * x[0],
        current[0] * x[1],
        current[1] * x[2],
    ]
}
fn reference(x: &mut [f64; 10], p: BuckBoostParameters, step: SwitchInterval, current: [f64; 2]) {
    let dt = step.duration_s / 256.0;
    for _ in 0..256 {
        let k1 = rhs(x, p, step.state, current);
        let k2 = rhs(
            &core::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]),
            p,
            step.state,
            current,
        );
        let k3 = rhs(
            &core::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]),
            p,
            step.state,
            current,
        );
        let k4 = rhs(
            &core::array::from_fn(|i| x[i] + dt * k3[i]),
            p,
            step.state,
            current,
        );
        for i in 0..10 {
            x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
}

#[test]
fn arbitrary_gate_order_resolves_channel_ripple_and_each_sources_own_voltage_work() {
    let p = parameters();
    let intervals = [
        interval(false, true, 0.004),
        interval(true, true, 0.007),
        interval(true, false, 0.003),
        interval(false, false, 0.002),
    ];
    let current = [[0.4, -0.7], [-0.9, 0.3], [0.8, -0.4], [-0.1, 0.6]];
    let mut oracle = [-1.2, 24.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let state = AveragedBuckBoost::new(p, modulation(0.0, 0.0))
        .unwrap()
        .state(oracle[0], [oracle[1], oracle[2]])
        .unwrap();
    let mut owner = prepared(p, &intervals);
    let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
    let step = owner
        .step_interval_sources_controlled(&state, &current, &mut next, &mut voltage, || false)
        .unwrap();
    for (gate, input) in intervals.into_iter().zip(current) {
        reference(&mut oracle, p, gate, input);
    }
    close(next[0] / p.inductance_h, oracle[0], 30.0);
    close(next[1] / p.bus_capacitance_f[0], oracle[1], 30.0);
    close(next[2] / p.bus_capacitance_f[1], oracle[2], 30.0);
    for i in 0..5 {
        close(step.conduction.joules[i], oracle[3 + i], 1.0);
    }
    for i in 0..2 {
        close(step.terminal_work_j[i], oracle[8 + i], 1.0);
    }
    close(
        step.native.delta_h + step.native.dissipated - step.native.supplied,
        0.0,
        100.0,
    );
    let mean_current: [f64; 2] = core::array::from_fn(|port| {
        intervals
            .iter()
            .zip(current)
            .map(|(gate, input)| gate.duration_s * input[port])
            .sum::<f64>()
            / owner.period_s()
    });
    assert!(
        (mean_current[0] * voltage[0] * owner.period_s() - step.terminal_work_j[0]).abs() > 1e-5
    );
    // Reusing the prepared caches does not leak an earlier accepted state.
    let (mut again, mut again_voltage) = ([0.0; 3], [0.0; 2]);
    let retry = owner
        .step_interval_sources_controlled(&state, &current, &mut again, &mut again_voltage, || {
            false
        })
        .unwrap();
    assert_eq!(again, next);
    assert_eq!(again_voltage, voltage);
    assert_eq!(retry.conduction, step.conduction);
}

fn ripple_parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 0.001,
        bus_capacitance_f: [0.1; 2],
        winding_resistance_ohm: 0.05,
        input_on_resistance_ohm: [0.02; 2],
        output_on_resistance_ohm: [0.03, 0.07],
    }
}

#[test]
fn nonzero_pwm_ripple_heat_is_resolved_when_the_averaged_current_is_exactly_zero() {
    let p = ripple_parameters();
    let t = 1e-4;
    let model = AveragedBuckBoost::new(p, modulation(0.5, 1.0)).unwrap();
    let initial = model.state(0.0, [48.0, 24.0]).unwrap();
    let mut averaged = model.prepare_controlled(t, budget(), || false).unwrap();
    let intervals = centered_pwm(t, modulation(0.5, 1.0)).unwrap();
    let mut switched = prepared(p, &intervals);
    let (mut state, mut voltage) = ([0.0; 3], [0.0; 2]);
    let coarse = averaged
        .step_controlled(&initial, &[0.0; 2], &mut state, &mut voltage, || false)
        .unwrap();
    assert!(coarse.native.dissipated.abs() < 1e-15);
    let exact = switched
        .step_controlled(&initial, &[0.0; 2], &mut state, &mut voltage, || false)
        .unwrap();
    let mut oracle = [0.0, 48.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for interval in intervals {
        reference(&mut oracle, p, interval, [0.0; 2]);
    }
    for channel in 0..5 {
        assert!((exact.conduction.joules[channel] - oracle[3 + channel]).abs() < 2e-13);
    }
    assert!(exact.native.dissipated > 1.19e-6 && exact.native.dissipated < 1.21e-6);
    assert_eq!(exact.conduction.joules[4], 0.0);
    close(state[0] / p.inductance_h, oracle[0], 1.0);
    close(state[1] / p.bus_capacitance_f[0], oracle[1], 50.0);
    close(state[2] / p.bus_capacitance_f[1], oracle[2], 50.0);
    // Independent stiff-bus leading-order triangular ripple: Δi_pp=1.2 A,
    // <i²>=Δi_pp²/12. Finite bus capacitances explain its small correction.
    let triangular = 0.1 * t * 1.2_f64.powi(2) / 12.0;
    assert!((exact.native.dissipated / triangular - 1.0).abs() < 3e-5);
}

#[test]
fn centered_switching_refines_to_averaged_dynamics_over_the_same_physical_horizon() {
    let p = ripple_parameters();
    let initial = AveragedBuckBoost::new(p, modulation(0.5, 1.0))
        .unwrap()
        .state(0.0, [48.0, 24.0])
        .unwrap();
    let mut errors = Vec::new();
    let mut heats = Vec::new();
    for cycles in [20, 40, 80] {
        let intervals = centered_pwm(0.002 / cycles as f64, modulation(0.5, 1.0)).unwrap();
        let mut owner = prepared(p, &intervals);
        let mut state = initial;
        let mut heat = 0.0;
        for _ in 0..cycles {
            let (mut next, mut voltage) = ([0.0; 3], [0.0; 2]);
            let step = owner
                .step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || false)
                .unwrap();
            state = next;
            heat += step.native.dissipated;
        }
        errors.push(
            (state[0] / p.inductance_h)
                .abs()
                .max((state[1] / p.bus_capacitance_f[0] - 48.0).abs())
                .max((state[2] / p.bus_capacitance_f[1] - 24.0).abs()),
        );
        heats.push(heat);
    }
    for series in [errors, heats] {
        assert!(series.iter().all(|v| *v > 0.0));
        for pair in series.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{series:?}"
            );
        }
    }
}

#[test]
fn invalid_edges_work_limits_sources_and_last_checkpoint_refuse_atomically() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(centered_pwm(bad, modulation(0.5, 0.5)).is_err());
        assert!(
            PreparedSwitchCycle::new_controlled(
                parameters(),
                &[interval(true, true, bad)],
                budget(),
                SwitchCycleControls::default(),
                || false
            )
            .is_err()
        );
    }
    for bad in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(centered_pwm(0.01, modulation(bad, 0.5)).is_err());
        assert!(centered_pwm(0.01, modulation(0.5, bad)).is_err());
    }
    for intervals in [
        Vec::new(),
        vec![interval(true, true, 0.001); MAX_SWITCH_INTERVALS + 1],
    ] {
        assert!(
            PreparedSwitchCycle::new_controlled(
                parameters(),
                &intervals,
                budget(),
                SwitchCycleControls::default(),
                || false
            )
            .is_err()
        );
    }
    let controls = SwitchCycleControls {
        max_intervals: 1,
        ..SwitchCycleControls::default()
    };
    assert!(
        PreparedSwitchCycle::new_controlled(
            parameters(),
            &[interval(true, true, 0.001); 2],
            budget(),
            controls,
            || false
        )
        .is_err()
    );
    let intervals = centered_pwm(0.01, modulation(0.5, 1.0)).unwrap();
    let mut owner = prepared(parameters(), &intervals);
    let state = AveragedBuckBoost::new(parameters(), modulation(0.5, 1.0))
        .unwrap()
        .state(1.0, [24.0, 12.0])
        .unwrap();
    let (mut next, mut voltage) = ([123.0; 3], [456.0; 2]);
    let mut polls = 0;
    owner
        .step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || {
            polls += 1;
            false
        })
        .unwrap();
    for boundary in [1, polls / 2, polls] {
        next = [123.0; 3];
        voltage = [456.0; 2];
        let mut seen = 0;
        assert!(matches!(
            owner.step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || {
                seen += 1;
                seen == boundary
            }),
            Err(PreparedStepError::Cancelled)
        ));
        assert_eq!(next, [123.0; 3]);
        assert_eq!(voltage, [456.0; 2]);
    }
    next = [123.0; 3];
    voltage = [456.0; 2];
    let mut sources = vec![[0.0; 2]; intervals.len()];
    sources.last_mut().unwrap()[1] = f64::NAN;
    assert!(
        owner
            .step_interval_sources_controlled(&state, &sources, &mut next, &mut voltage, || false)
            .is_err()
    );
    assert_eq!(next, [123.0; 3]);
    assert_eq!(voltage, [456.0; 2]);
    assert!(
        owner
            .step_interval_sources_controlled(&state, &[], &mut next, &mut voltage, || false)
            .is_err()
    );
    owner
        .step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || false)
        .unwrap();
}
