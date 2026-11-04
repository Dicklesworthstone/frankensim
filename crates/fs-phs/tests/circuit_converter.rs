//! G1/G2/G3: original electrical dynamics, reciprocal power and physical heat.
use fs_phs::PreparedStepError;

#[path = "circuit_converter/switching.rs"]
mod switching;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};

fn parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        winding_resistance_ohm: 0.03,
        input_on_resistance_ohm: [0.02, 0.04],
        output_on_resistance_ohm: [0.015, 0.025],
    }
}
fn ideal() -> BuckBoostParameters {
    BuckBoostParameters {
        winding_resistance_ohm: 0.0,
        input_on_resistance_ohm: [0.0; 2],
        output_on_resistance_ohm: [0.0; 2],
        ..parameters()
    }
}
fn modulation(a: f64, b: f64) -> BridgeModulation {
    BridgeModulation {
        input_high_fraction: a,
        output_high_fraction: b,
    }
}
fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 3,
        max_omega_dt: 0.04,
    }
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * scale,
        "{actual:.17e} != {expected:.17e} (scale={scale:e})"
    );
}
fn resistances(p: BuckBoostParameters, a: f64, b: f64) -> [f64; 5] {
    [
        p.winding_resistance_ohm,
        a * p.input_on_resistance_ohm[0],
        (1.0 - a) * p.input_on_resistance_ohm[1],
        b * p.output_on_resistance_ohm[0],
        (1.0 - b) * p.output_on_resistance_ohm[1],
    ]
}

#[test]
fn buck_boost_and_reverse_conversion_keep_both_terminal_work_and_every_heat_owner() {
    let p = parameters();
    let dt = 0.037;
    for (a, b) in [(0.4, 1.0), (1.0, 0.4), (0.7, 0.6)] {
        let resistance = resistances(p, a, b);
        let total: f64 = resistance.iter().sum();
        for i in [-3.0, 3.0] {
            let voltage = [24.0, (a * 24.0 - total * i) / b];
            let current = [a * i, -b * i];
            let model = AveragedBuckBoost::new(p, modulation(a, b)).unwrap();
            let state = model.state(i, voltage).unwrap();
            let mut prepared = model.prepare_controlled(dt, budget(), || false).unwrap();
            let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
            let step = prepared
                .step_controlled(&state, &current, &mut next, &mut mean, || false)
                .unwrap();
            for k in 0..3 {
                close(next[k], state[k], state[k].abs().max(1.0));
            }
            for k in 0..2 {
                close(mean[k], voltage[k], voltage[k]);
                close(step.terminal_work_j[k], dt * voltage[k] * current[k], 10.0);
            }
            for (heat, r) in step.conduction.joules.iter().zip(resistance) {
                close(*heat, dt * i * i * r, 1.0);
            }
            assert_eq!(step.terminal_work_j[0].is_sign_positive(), i > 0.0);
            assert_eq!(step.terminal_work_j[1].is_sign_positive(), i < 0.0);
            assert!(step.native.dissipated > 0.0);
            close(step.native.dissipated, dt * i * i * total, 1.0);
            close(step.native.delta_h, 0.0, 10.0);
            close(
                step.terminal_work_j.iter().sum(),
                step.native.dissipated,
                10.0,
            );
        }
    }
}

#[test]
fn lossless_transient_matches_independent_oscillator_and_trajectory_integrals() {
    let p = ideal();
    let (a, b) = (0.63, 0.48);
    let (i0, v0, current, dt) = (1.7, [24.0, 8.0], [0.3, -0.4], 0.13_f64);
    let model = AveragedBuckBoost::new(p, modulation(a, b)).unwrap();
    let state = model.state(i0, v0).unwrap();
    // Eliminate the two capacitor equations only for this analytic oracle:
    // L i'' + (a^2/C1+b^2/C2)i = a I1/C1-b I2/C2.
    let k = a * a / p.bus_capacitance_f[0] + b * b / p.bus_capacitance_f[1];
    let w = (k / p.inductance_h).sqrt();
    let forced_i =
        (a * current[0] / p.bus_capacitance_f[0] - b * current[1] / p.bus_capacitance_f[1]) / k;
    let drive = a * v0[0] - b * v0[1];
    let (s, c) = (w * dt).sin_cos();
    let end_i = forced_i + (i0 - forced_i) * c + drive / (p.inductance_h * w) * s;
    let integrated_i =
        forced_i * dt + (i0 - forced_i) * s / w + drive / (p.inductance_h * w * w) * (1.0 - c);
    let twice_integrated_i = forced_i * dt * dt / 2.0
        + (i0 - forced_i) * (1.0 - c) / (w * w)
        + drive / (p.inductance_h * w) * (dt / w - s / (w * w));
    let expected_v = [
        v0[0] + (current[0] * dt - a * integrated_i) / p.bus_capacitance_f[0],
        v0[1] + (current[1] * dt + b * integrated_i) / p.bus_capacitance_f[1],
    ];
    let mean_v = [
        v0[0] + (current[0] * dt / 2.0 - a * twice_integrated_i / dt) / p.bus_capacitance_f[0],
        v0[1] + (current[1] * dt / 2.0 + b * twice_integrated_i / dt) / p.bus_capacitance_f[1],
    ];
    let mut prepared = model.prepare_controlled(dt, budget(), || false).unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let step = prepared
        .step_controlled(&state, &current, &mut next, &mut mean, || false)
        .unwrap();
    close(next[0] / p.inductance_h, end_i, 100.0);
    for k in 0..2 {
        close(next[k + 1] / p.bus_capacitance_f[k], expected_v[k], 30.0);
        close(mean[k], mean_v[k], 30.0);
        close(step.terminal_work_j[k], current[k] * dt * mean_v[k], 10.0);
    }
    close(
        b * next[1] + a * next[2],
        b * state[1] + a * state[2] + dt * (b * current[0] + a * current[1]),
        10.0,
    );
    assert_eq!(step.conduction.joules, [0.0; 5]);
    close(
        step.native.delta_h,
        step.terminal_work_j.iter().sum(),
        100.0,
    );
    // A final-voltage work approximation is detectably different here.
    assert!((step.terminal_work_j[0] - current[0] * dt * expected_v[0]).abs() > 1e-3);
}

#[test]
fn low_side_freewheel_has_exact_resistive_decay_without_draining_either_bus() {
    let p = parameters();
    let dt = 0.41;
    let i = -2.3;
    let r = p.winding_resistance_ohm + p.input_on_resistance_ohm[1] + p.output_on_resistance_ohm[1];
    let decay = (-r * dt / p.inductance_h).exp();
    let integral_i2 = i * i * p.inductance_h / (2.0 * r) * (1.0 - decay * decay);
    let model = AveragedBuckBoost::new(p, modulation(0.0, 0.0)).unwrap();
    let state = model.state(i, [12.0, 25.0]).unwrap();
    let mut prepared = model.prepare_controlled(dt, budget(), || false).unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let step = prepared
        .step_controlled(&state, &[0.0; 2], &mut next, &mut mean, || false)
        .unwrap();
    close(next[0], state[0] * decay, 1.0);
    close(next[1], state[1], 10.0);
    close(next[2], state[2], 10.0);
    for (heat, resistance) in step.conduction.joules.iter().zip(resistances(p, 0.0, 0.0)) {
        close(*heat, resistance * integral_i2, 1.0);
    }
    assert_eq!(step.conduction.joules[1], 0.0);
    assert_eq!(step.conduction.joules[3], 0.0);
    close(step.native.dissipated, -step.native.delta_h, 1.0);
}

// Independent physical-current/voltage ODE; no production matrices or port
// outputs reconstruct this reference. The final two coordinates integrate
// winding/channel heat and actual electrical source work on the trajectory.
fn rk4(state: &mut [f64; 5], dt: f64, rhs: impl Fn(&[f64; 5]) -> [f64; 5]) {
    let k1 = rhs(state);
    let k2 = rhs(&core::array::from_fn(|i| state[i] + 0.5 * dt * k1[i]));
    let k3 = rhs(&core::array::from_fn(|i| state[i] + 0.5 * dt * k2[i]));
    let k4 = rhs(&core::array::from_fn(|i| state[i] + dt * k3[i]));
    for i in 0..5 {
        state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn modulation_changes_keep_original_flux_charge_and_match_independent_forced_dynamics() {
    let p = parameters();
    let dt = 0.012;
    let current = [0.7, -0.9];
    let mut oracle = [2.0, 24.0, 12.0, 0.0, 0.0];
    let first = AveragedBuckBoost::new(p, modulation(0.5, 1.0)).unwrap();
    let mut state = first.state(oracle[0], [oracle[1], oracle[2]]).unwrap();
    let initial_energy = first.system().hamiltonian(&state);
    let (mut work, mut heat) = (0.0, 0.0);
    for (a, b) in [
        (0.5, 1.0),
        (0.8, 1.0),
        (1.0, 0.7),
        (1.0, 0.4),
        (0.3, 0.8),
        (0.6, 0.5),
    ] {
        let model = AveragedBuckBoost::new(p, modulation(a, b)).unwrap();
        let r: f64 = resistances(p, a, b).iter().sum();
        let mut prepared = model.prepare_controlled(dt, budget(), || false).unwrap();
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        let step = prepared
            .step_controlled(&state, &current, &mut next, &mut mean, || false)
            .unwrap();
        for _ in 0..128 {
            rk4(&mut oracle, dt / 128.0, |x| {
                [
                    (a * x[1] - b * x[2] - r * x[0]) / p.inductance_h,
                    (current[0] - a * x[0]) / p.bus_capacitance_f[0],
                    (current[1] + b * x[0]) / p.bus_capacitance_f[1],
                    r * x[0] * x[0],
                    current[0] * x[1] + current[1] * x[2],
                ]
            });
        }
        state = next;
        heat += step.conduction.joules.iter().sum::<f64>();
        work += step.terminal_work_j.iter().sum::<f64>();
        close(state[0] / p.inductance_h, oracle[0], 50.0);
        close(state[1] / p.bus_capacitance_f[0], oracle[1], 30.0);
        close(state[2] / p.bus_capacitance_f[1], oracle[2], 30.0);
        close(heat, oracle[3], 10.0);
        close(work, oracle[4], 10.0);
        close(
            prepared.system().hamiltonian(&state) - initial_energy + heat - work,
            0.0,
            100.0,
        );
    }
}

#[test]
fn original_source_impedance_and_load_close_the_live_converter_ports() {
    let p = parameters();
    let (a, b) = (0.65, 0.8);
    let (source_v, source_r, load_r) = (24.0, 0.2, 6.0);
    let model = AveragedBuckBoost::new(p, modulation(a, b)).unwrap();
    let mut state = model.state(0.0, [24.0, 0.0]).unwrap();
    let system = model
        .into_system()
        .with_port_conductance(&[0, 1], &[1.0 / source_r, 0.0, 0.0, 1.0 / load_r])
        .unwrap();
    let dt = 0.2;
    let mut prepared = PreparedExactZoh::new_owned(system, dt, budget()).unwrap();
    let mut oracle = [0.0, 24.0, 0.0, 0.0, 0.0];
    let r: f64 = resistances(p, a, b).iter().sum();
    for _ in 0..64 {
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        prepared
            .step_into(&state, &[source_v / source_r, 0.0], &mut next, &mut mean)
            .unwrap();
        for _ in 0..512 {
            rk4(&mut oracle, dt / 512.0, |x| {
                [
                    (a * x[1] - b * x[2] - r * x[0]) / p.inductance_h,
                    ((source_v - x[1]) / source_r - a * x[0]) / p.bus_capacitance_f[0],
                    (b * x[0] - x[2] / load_r) / p.bus_capacitance_f[1],
                    0.0,
                    0.0,
                ]
            });
        }
        state = next;
        close(state[0] / p.inductance_h, oracle[0], 100.0);
        close(state[1] / p.bus_capacitance_f[0], oracle[1], 100.0);
        close(state[2] / p.bus_capacitance_f[1], oracle[2], 100.0);
    }
    let current = a * source_v / (r + source_r * a * a + load_r * b * b);
    close(state[0] / p.inductance_h, current, 100.0);
    close(
        state[1] / p.bus_capacitance_f[0],
        source_v - source_r * a * current,
        100.0,
    );
    close(
        state[2] / p.bus_capacitance_f[1],
        load_r * b * current,
        100.0,
    );
}

#[test]
fn malformed_data_budgets_and_cancelled_or_nonfinite_trials_refuse_without_publication() {
    for bad in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
        assert!(AveragedBuckBoost::new(parameters(), modulation(bad, 0.5)).is_err());
        assert!(AveragedBuckBoost::new(parameters(), modulation(0.5, bad)).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY, -1.0, 0.0, f64::from_bits(1)] {
        assert!(
            AveragedBuckBoost::new(
                BuckBoostParameters {
                    inductance_h: bad,
                    ..parameters()
                },
                modulation(0.5, 1.0)
            )
            .is_err()
        );
    }
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(
            AveragedBuckBoost::new(
                BuckBoostParameters {
                    winding_resistance_ohm: bad,
                    ..parameters()
                },
                modulation(0.5, 1.0)
            )
            .is_err()
        );
    }
    assert!(
        AveragedBuckBoost::new(
            BuckBoostParameters {
                winding_resistance_ohm: f64::MAX,
                input_on_resistance_ohm: [f64::MAX; 2],
                ..parameters()
            },
            modulation(0.5, 1.0)
        )
        .is_err()
    );
    let model = || AveragedBuckBoost::new(parameters(), modulation(0.5, 1.0)).unwrap();
    assert!(
        model()
            .prepare_controlled(
                0.01,
                LinearSplitBudget {
                    max_modes: 2,
                    ..budget()
                },
                || false
            )
            .is_err()
    );
    assert!(
        model()
            .prepare_controlled(-0.1, budget(), || false)
            .is_err()
    );
    assert!(matches!(
        model().prepare_controlled(0.01, budget(), || true),
        Err(PreparedStepError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        model().prepare_controlled(0.01, budget(), || {
            polls += 1;
            polls == 2
        }),
        Err(PreparedStepError::Cancelled)
    ));
    let initial = model().state(1.0, [24.0, 12.0]).unwrap();
    let mut prepared = model()
        .prepare_controlled(0.01, budget(), || false)
        .unwrap();
    for boundary in [1, 2] {
        let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
        let mut polls = 0;
        assert!(matches!(
            prepared.step_controlled(&initial, &[0.0; 2], &mut next, &mut mean, || {
                polls += 1;
                polls == boundary
            }),
            Err(PreparedStepError::Cancelled)
        ));
        assert_eq!(next, [123.0; 3]);
        assert_eq!(mean, [456.0; 2]);
    }
    let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
    assert!(
        prepared
            .step_controlled(&initial, &[f64::NAN, 0.0], &mut next, &mut mean, || false)
            .is_err()
    );
    assert_eq!(next, [123.0; 3]);
    assert_eq!(mean, [456.0; 2]);
    let retry = prepared
        .step_controlled(&initial, &[0.0; 2], &mut next, &mut mean, || false)
        .unwrap();
    assert!(retry.native.dissipated > 0.0);
    let mut zero = model().prepare_controlled(0.0, budget(), || false).unwrap();
    let rec = zero
        .step_controlled(&initial, &[1.0, -2.0], &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!(next, initial);
    assert_eq!(rec.terminal_work_j, [0.0; 2]);
    assert_eq!(rec.conduction.joules, [0.0; 5]);
}
