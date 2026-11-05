//! Physical source/load powers on the original converter's exact trajectory.
use super::{budget, close, modulation, parameters, resistances};
use fs_phs::PreparedStepError;
use fs_phs::circuit::converter::AveragedBuckBoost;
use fs_phs::linear_exact::DenseExponentialBudget;

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 6,
        max_dense_bytes: 16_384,
        max_multiply_terms: 500_000,
    }
}

#[test]
fn voltage_driven_buck_boost_and_regeneration_keep_every_physical_heat_owner() {
    let p = parameters();
    let terminal_r = [0.2, 6.0];
    let dt = 0.07;
    for (a, b) in [(0.4, 1.0), (1.0, 0.4)] {
        let channel_r = resistances(p, a, b);
        let r: f64 = channel_r.iter().sum();
        for source in [[24.0, 0.0], [0.0, 24.0], [24.0, 18.0]] {
            // Independent stationary Kirchhoff reduction in amperes/volts.
            let i = (a * source[0] - b * source[1])
                / (r + a * a * terminal_r[0] + b * b * terminal_r[1]);
            let terminal_i = [a * i, -b * i];
            let v = [
                source[0] - terminal_r[0] * terminal_i[0],
                source[1] - terminal_r[1] * terminal_i[1],
            ];
            let model = AveragedBuckBoost::new(p, modulation(a, b)).unwrap();
            let initial = model.state(i, v).unwrap();
            let mut prepared = model
                .prepare_voltage_driven_controlled(terminal_r, dt, budget(), heat_budget(), || {
                    false
                })
                .unwrap();
            let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
            let record = prepared
                .step_controlled(&initial, &source, &mut next, &mut mean, || false)
                .unwrap();
            for k in 0..3 {
                close(next[k], initial[k], 100.0);
            }
            for k in 0..2 {
                close(mean[k], v[k], 100.0);
                close(
                    record.source_work_j[k],
                    dt * source[k] * terminal_i[k],
                    100.0,
                );
                close(
                    record.terminal_resistor_heat_j[k],
                    dt * terminal_r[k] * terminal_i[k] * terminal_i[k],
                    100.0,
                );
                close(
                    record.converter.terminal_work_j[k],
                    dt * v[k] * terminal_i[k],
                    100.0,
                );
                if source[k] * terminal_i[k] != 0.0 {
                    assert_eq!(
                        record.source_work_j[k].is_sign_positive(),
                        (source[k] * terminal_i[k]).is_sign_positive()
                    );
                }
            }
            for (heat, resistance) in record.converter.conduction.joules.iter().zip(channel_r) {
                close(*heat, dt * resistance * i * i, 100.0);
            }
            close(record.converter.native.delta_h, 0.0, 100.0);
            close(record.converter.native.dissipated, dt * r * i * i, 100.0);
            close(record.converter.native.supplied, dt * r * i * i, 100.0);
        }
    }
}

#[test]
fn freewheel_rc_transients_distinguish_series_heat_from_norton_shunt_heat() {
    let p = parameters();
    let terminal_r = [0.4, 0.7];
    let (dt, i0, v0, source) = (0.41, -2.3, [24.000_000_1, 28.0], [24.0, 12.0]);
    let model = AveragedBuckBoost::new(p, modulation(0.0, 0.0)).unwrap();
    let initial = model.state(i0, v0).unwrap();
    let mut prepared = model
        .prepare_voltage_driven_controlled(terminal_r, dt, budget(), heat_budget(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    let record = prepared
        .step_controlled(&initial, &source, &mut next, &mut mean, || false)
        .unwrap();
    for k in 0..2 {
        let c = p.bus_capacitance_f[k];
        let rate = 1.0 / (terminal_r[k] * c);
        let decay = (-rate * dt).exp();
        let difference = v0[k] - source[k];
        let endpoint = source[k] + difference * decay;
        let mean_v = source[k] + difference * -(-rate * dt).exp_m1() / (rate * dt);
        let heat = 0.5 * c * difference * difference * -(-2.0 * rate * dt).exp_m1();
        close(next[k + 1] / c, endpoint, 100.0);
        close(mean[k], mean_v, 100.0);
        close(
            record.source_work_j[k],
            source[k] * c * (endpoint - v0[k]),
            100.0,
        );
        assert!(record.terminal_resistor_heat_j[k] > 0.0);
        assert!((record.terminal_resistor_heat_j[k] / heat - 1.0).abs() < 2e-6);
        close(
            record.converter.terminal_work_j[k],
            0.5 * c * (endpoint - v0[k]) * (endpoint + v0[k]),
            100.0,
        );
    }
    assert!(
        record.source_work_j[1] < 0.0,
        "overcharged bus must regenerate"
    );
    assert!(
        record.terminal_resistor_heat_j[0] < 1e-14,
        "near-matched physical source resistor must not receive Norton shunt heat"
    );
    let r: f64 = resistances(p, 0.0, 0.0).iter().sum();
    let inductor_decay = (-r * dt / p.inductance_h).exp();
    close(next[0], initial[0] * inductor_decay, 1.0);
    close(
        record.converter.native.dissipated,
        0.5 * p.inductance_h * i0 * i0 * (1.0 - inductor_decay * inductor_decay),
        1.0,
    );
}

fn rk4(state: &mut [f64; 8], dt: f64, rhs: impl Fn(&[f64; 8]) -> [f64; 8]) {
    let k1 = rhs(state);
    let k2 = rhs(&core::array::from_fn(|i| state[i] + 0.5 * dt * k1[i]));
    let k3 = rhs(&core::array::from_fn(|i| state[i] + 0.5 * dt * k2[i]));
    let k4 = rhs(&core::array::from_fn(|i| state[i] + dt * k3[i]));
    for i in 0..8 {
        state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn loaded_startup_and_changed_voltage_sources_match_independent_power_odes() {
    let p = parameters();
    let (a, b, dt) = (0.65, 0.8, 0.05);
    let terminal_r = [0.2, 6.0];
    let r: f64 = resistances(p, a, b).iter().sum();
    let model = AveragedBuckBoost::new(p, modulation(a, b)).unwrap();
    let mut initial = model.state(0.0, [24.0, 0.0]).unwrap();
    let mut prepared = model
        .prepare_voltage_driven_controlled(terminal_r, dt, budget(), heat_budget(), || false)
        .unwrap();
    let mut oracle = [0.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for source in [[24.0, 0.0], [24.0, 8.0], [8.0, 24.0], [-12.0, 0.0]] {
        let before = oracle;
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        let record = prepared
            .step_controlled(&initial, &source, &mut next, &mut mean, || false)
            .unwrap();
        for _ in 0..512 {
            rk4(&mut oracle, dt / 512.0, |x| {
                let current = [
                    (source[0] - x[1]) / terminal_r[0],
                    (source[1] - x[2]) / terminal_r[1],
                ];
                [
                    (a * x[1] - b * x[2] - r * x[0]) / p.inductance_h,
                    (current[0] - a * x[0]) / p.bus_capacitance_f[0],
                    (current[1] + b * x[0]) / p.bus_capacitance_f[1],
                    terminal_r[0] * current[0] * current[0],
                    terminal_r[1] * current[1] * current[1],
                    r * x[0] * x[0],
                    source[0] * current[0],
                    source[1] * current[1],
                ]
            });
        }
        close(next[0] / p.inductance_h, oracle[0], 100.0);
        close(next[1] / p.bus_capacitance_f[0], oracle[1], 100.0);
        close(next[2] / p.bus_capacitance_f[1], oracle[2], 100.0);
        for k in 0..2 {
            let heat = oracle[3 + k] - before[3 + k];
            let work = oracle[6 + k] - before[6 + k];
            close(record.terminal_resistor_heat_j[k], heat, 100.0);
            close(record.source_work_j[k], work, 100.0);
            close(record.converter.terminal_work_j[k], work - heat, 100.0);
            let mean_current = (source[k] - mean[k]) / terminal_r[k];
            assert!(
                record.terminal_resistor_heat_j[k]
                    > dt * terminal_r[k] * mean_current * mean_current,
                "transient RMS current must retain additional physical heat"
            );
        }
        close(
            record.converter.native.dissipated,
            oracle[5] - before[5],
            100.0,
        );
        initial = next;
    }
}

#[test]
fn voltage_terminal_admission_and_every_cancellation_checkpoint_preserve_outputs() {
    let model = || AveragedBuckBoost::new(parameters(), modulation(0.6, 0.8)).unwrap();
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(
            model()
                .prepare_voltage_driven_controlled(
                    [bad, 6.0],
                    0.02,
                    budget(),
                    heat_budget(),
                    || false,
                )
                .is_err()
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model()
                .prepare_voltage_driven_controlled([0.2, 6.0], bad, budget(), heat_budget(), || {
                    false
                },)
                .is_err()
        );
    }
    assert!(
        model()
            .prepare_voltage_driven_controlled(
                [0.2, 6.0],
                0.02,
                budget(),
                DenseExponentialBudget {
                    max_dimension: 5,
                    ..heat_budget()
                },
                || false,
            )
            .is_err()
    );
    assert!(matches!(
        model().prepare_voltage_driven_controlled(
            [0.2, 6.0],
            0.02,
            budget(),
            heat_budget(),
            || true,
        ),
        Err(PreparedStepError::Cancelled)
    ));
    let initial = model().state(1.0, [24.0, 12.0]).unwrap();
    let mut prepared = model()
        .prepare_voltage_driven_controlled([0.2, 6.0], 0.02, budget(), heat_budget(), || false)
        .unwrap();
    let (mut expected, mut expected_mean, mut polls) = ([0.0; 3], [0.0; 2], 0);
    prepared
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
        let (mut next, mut mean, mut at) = ([123.0; 3], [456.0; 2], 0);
        assert!(matches!(
            prepared.step_controlled(&initial, &[24.0, 0.0], &mut next, &mut mean, || {
                at += 1;
                at == stop
            },),
            Err(PreparedStepError::Cancelled)
        ));
        assert_eq!(next, [123.0; 3]);
        assert_eq!(mean, [456.0; 2]);
    }
    for bad in [f64::NAN, f64::INFINITY, f64::from_bits(1), f64::MAX] {
        let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
        assert!(
            prepared
                .step_controlled(&initial, &[bad, 0.0], &mut next, &mut mean, || false)
                .is_err()
        );
        assert_eq!(next, [123.0; 3]);
        assert_eq!(mean, [456.0; 2]);
    }
    let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
    prepared
        .step_controlled(&initial, &[24.0, 0.0], &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!(next, expected);
    assert_eq!(mean, expected_mean);
}

fn unit_freewheel() -> AveragedBuckBoost {
    let mut p = super::ideal();
    p.inductance_h = 1.0;
    p.bus_capacitance_f = [1.0; 2];
    AveragedBuckBoost::new(p, modulation(0.0, 0.0)).unwrap()
}

#[test]
fn tiny_source_current_is_integrated_without_subtracting_a_rounded_mean_voltage() {
    let dt = 0.3_f64;
    for offset in [1e-7, 1e-10] {
        let model = unit_freewheel();
        // The second bus has a resolved passive discharge, while the input
        // source independently receives a much smaller regenerative current.
        let initial = model.state(0.0, [24.0 + offset, 10.0]).unwrap();
        let difference = initial[1] - 24.0;
        let source_work = 24.0 * difference * (-dt).exp_m1();
        let heat = -0.5 * difference * difference * (-2.0 * dt).exp_m1();
        let mut prepared = model
            .prepare_voltage_driven_controlled([1.0; 2], dt, budget(), heat_budget(), || false)
            .unwrap();
        let (mut next, mut mean) = ([0.0; 3], [0.0; 2]);
        let record = prepared
            .step_controlled(&initial, &[24.0, 0.0], &mut next, &mut mean, || false)
            .unwrap();
        for (actual, expected) in [
            (record.source_work_j[0], source_work),
            (record.terminal_resistor_heat_j[0], heat),
            (record.converter.terminal_work_j[0], source_work - heat),
        ] {
            assert!(
                (actual - expected).abs() <= 3e-12 * expected.abs(),
                "{actual:e} != {expected:e}",
            );
        }
    }
}

#[test]
fn unresolved_small_storage_transfers_refuse_without_using_the_large_energy_datum() {
    for (voltage, offset) in [(24.0, 1e-7), (24.0, 1e-10), (1e16, 2.0)] {
        let model = unit_freewheel();
        let initial = model.state(0.0, [voltage + offset, 0.0]).unwrap();
        let mut prepared = model
            .prepare_voltage_driven_controlled([1.0; 2], 0.3, budget(), heat_budget(), || false)
            .unwrap();
        let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
        assert!(
            prepared
                .step_controlled(&initial, &[voltage, 0.0], &mut next, &mut mean, || false,)
                .is_err()
        );
        assert_eq!(next, [123.0; 3]);
        assert_eq!(mean, [456.0; 2]);
        // Refusal retains the prepared cache: a representable driven interval
        // with the same coefficients remains usable.
        prepared
            .step_controlled(&[0.0; 3], &[1.0, 0.0], &mut next, &mut mean, || false)
            .unwrap();
    }
}

#[test]
fn exact_unloaded_equilibrium_preserves_state_without_exponential_roundoff() {
    for voltage in [[24.0, -12.0], [1e16, 0.0]] {
        let model = unit_freewheel();
        let initial = model.state(0.0, voltage).unwrap();
        let mut prepared = model
            .prepare_voltage_driven_controlled([1.0; 2], 0.3, budget(), heat_budget(), || false)
            .unwrap();
        let (mut next, mut mean) = ([123.0; 3], [456.0; 2]);
        let record = prepared
            .step_controlled(&initial, &voltage, &mut next, &mut mean, || false)
            .unwrap();
        assert_eq!(next.map(f64::to_bits), initial.map(f64::to_bits));
        assert_eq!(mean, voltage);
        assert_eq!(record.source_work_j, [0.0; 2]);
        assert_eq!(record.terminal_resistor_heat_j, [0.0; 2]);
        assert_eq!(record.converter.terminal_work_j[0], 0.0);
        assert_eq!(record.converter.terminal_work_j[1], 0.0);
        assert_eq!(record.converter.native.delta_h, 0.0);
        assert_eq!(record.converter.native.dissipated, 0.0);
    }
}
