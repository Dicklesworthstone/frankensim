//! Prescribed switching with independently integrated finite source/storage power.
use fs_phs::PreparedStepError;
use fs_phs::circuit::converter::capacitor_terminals::CapacitorTerminals;
use fs_phs::circuit::converter::switching::capacitor_fed::PreparedCapacitorFedSwitchCycle;
use fs_phs::circuit::converter::switching::{
    BridgeSwitchState, SwitchCycleControls, SwitchInterval,
};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};

fn parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        winding_resistance_ohm: 0.03,
        input_on_resistance_ohm: [0.02, 0.04],
        output_on_resistance_ohm: [0.015, 0.025],
    }
}
fn terminals() -> CapacitorTerminals {
    CapacitorTerminals {
        capacitance_f: [0.5, 0.8],
        resistance_ohm: [0.2, 0.7],
    }
}
fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 5,
        max_omega_dt: 0.04,
    }
}
fn observations() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 10,
        max_dense_bytes: 1_000_000,
        max_multiply_terms: 100_000_000,
    }
}
fn interval(input_high: bool, output_high: bool, duration_s: f64) -> SwitchInterval {
    SwitchInterval {
        state: BridgeSwitchState {
            input_high,
            output_high,
        },
        duration_s,
    }
}
fn cycle(intervals: &[SwitchInterval]) -> PreparedCapacitorFedSwitchCycle {
    PreparedCapacitorFedSwitchCycle::new_controlled(
        parameters(),
        terminals(),
        intervals,
        budget(),
        observations(),
        SwitchCycleControls::default(),
        || false,
    )
    .unwrap()
}
fn initial() -> [f64; 5] {
    AveragedBuckBoost::new(
        parameters(),
        BridgeModulation {
            input_high_fraction: 0.0,
            output_high_fraction: 0.0,
        },
    )
    .unwrap()
    .with_capacitor_terminals(terminals())
    .unwrap()
    .state(1.7, [18.0, 7.0], [24.0, 3.0])
    .unwrap()
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 5e-11 * scale, "{a:.17e} != {b:.17e}");
}
fn rk4<const N: usize>(x: &mut [f64; N], h: f64, rhs: impl Fn(&[f64; N]) -> [f64; N]) {
    let k1 = rhs(x);
    let k2 = rhs(&core::array::from_fn(|i| x[i] + h * k1[i] * 0.5));
    let k3 = rhs(&core::array::from_fn(|i| x[i] + h * k2[i] * 0.5));
    let k4 = rhs(&core::array::from_fn(|i| x[i] + h * k3[i]));
    for i in 0..N {
        x[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn every_gate_keeps_finite_source_sag_ripple_heat_and_correlated_charger_work() {
    let intervals = [
        interval(true, true, 0.012),
        interval(false, false, 0.008),
        interval(true, false, 0.012),
        interval(false, false, 0.008),
    ];
    let external = [[0.4, -0.1], [0.0, 0.0], [0.1, 0.2], [-0.2, 0.1]];
    let mut native = cycle(&intervals);
    assert_eq!(native.cached_interval_count(), 3);
    let mut oracle = [0.0; 22];
    oracle[..5].copy_from_slice(&[1.7, 18.0, 7.0, 24.0, 3.0]);
    let p = parameters();
    let terminal = terminals();
    for (gate, source) in intervals.iter().zip(external) {
        let a = if gate.state.input_high { 1.0 } else { 0.0 };
        let b = if gate.state.output_high { 1.0 } else { 0.0 };
        let resistance = [
            p.winding_resistance_ohm,
            a * p.input_on_resistance_ohm[0],
            (1.0 - a) * p.input_on_resistance_ohm[1],
            b * p.output_on_resistance_ohm[0],
            (1.0 - b) * p.output_on_resistance_ohm[1],
        ];
        for _ in 0..256 {
            rk4(&mut oracle, gate.duration_s / 256.0, |x| {
                let flow = [
                    (x[3] - x[1]) / terminal.resistance_ohm[0],
                    (x[4] - x[2]) / terminal.resistance_ohm[1],
                ];
                let mut derivative = [0.0; 22];
                derivative[0] =
                    (a * x[1] - b * x[2] - resistance.iter().sum::<f64>() * x[0]) / p.inductance_h;
                derivative[1] = (flow[0] - a * x[0]) / p.bus_capacitance_f[0];
                derivative[2] = (flow[1] + b * x[0]) / p.bus_capacitance_f[1];
                derivative[3] = (source[0] - flow[0]) / terminal.capacitance_f[0];
                derivative[4] = (source[1] - flow[1]) / terminal.capacitance_f[1];
                for k in 0..5 {
                    derivative[5 + k] = resistance[k] * x[0] * x[0];
                }
                for k in 0..2 {
                    derivative[10 + k] = terminal.resistance_ohm[k] * flow[k] * flow[k];
                    derivative[12 + k] = x[1 + k] * flow[k];
                    derivative[14 + k] = x[3 + k] * flow[k];
                    derivative[16 + k] = x[3 + k] * source[k];
                }
                derivative[18..22].copy_from_slice(&x[1..5]);
                derivative
            });
        }
    }
    let (mut next, mut mean) = ([0.0; 5], [0.0; 4]);
    let result = native
        .step_interval_sources_controlled(&initial(), &external, &mut next, &mut mean, || false)
        .unwrap();
    for (k, coefficient) in [
        p.inductance_h,
        p.bus_capacitance_f[0],
        p.bus_capacitance_f[1],
        terminal.capacitance_f[0],
        terminal.capacitance_f[1],
    ]
    .into_iter()
    .enumerate()
    {
        close(next[k] / coefficient, oracle[k], 40.0);
    }
    for k in 0..5 {
        close(result.converter.conduction.joules[k], oracle[5 + k], 1.0);
    }
    for k in 0..2 {
        close(result.terminal_resistor_heat_j[k], oracle[10 + k], 10.0);
        close(result.converter.terminal_work_j[k], oracle[12 + k], 40.0);
        close(result.terminal_outflow_work_j[k], oracle[14 + k], 40.0);
        close(result.external_work_j[k], oracle[16 + k], 1.0);
        let initial_v = [24.0, 3.0][k];
        close(
            result.terminal_storage_change_j[k],
            0.5 * terminal.capacitance_f[k]
                * (oracle[3 + k] * oracle[3 + k] - initial_v * initial_v),
            50.0,
        );
    }
    for k in 0..4 {
        close(mean[k], oracle[18 + k] / native.period_s(), 40.0);
    }
    let heat: f64 = oracle[5..12].iter().sum();
    close(result.total.dissipated, heat, 10.0);
    close(result.total.delta_h, oracle[16] + oracle[17] - heat, 10.0);
    assert!(next[3] < initial()[3], "finite input supply must discharge");
    assert!(
        next[4] > initial()[4],
        "receiving terminal capacitor must charge"
    );
    assert_eq!(native.cached_interval_count(), 3);
}

#[test]
fn static_gate_has_exact_native_parity_and_cycle_refusal_is_atomic() {
    let gates = [interval(true, false, 0.03)];
    let mut switched = cycle(&gates);
    let mut held = AveragedBuckBoost::new(
        parameters(),
        BridgeModulation {
            input_high_fraction: 1.0,
            output_high_fraction: 0.0,
        },
    )
    .unwrap()
    .with_capacitor_terminals(terminals())
    .unwrap()
    .prepare_controlled(0.03, budget(), observations(), || false)
    .unwrap();
    let source = [0.2, -0.1];
    let (mut direct, mut direct_mean, mut actual, mut mean) =
        ([0.0; 5], [0.0; 4], [0.0; 5], [0.0; 4]);
    let expected = held
        .step_controlled(&initial(), &source, &mut direct, &mut direct_mean, || false)
        .unwrap();
    let mut polls = 0;
    let result = switched
        .step_controlled(&initial(), &source, &mut actual, &mut mean, || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(direct.map(f64::to_bits), actual.map(f64::to_bits));
    assert_eq!(direct_mean.map(f64::to_bits), mean.map(f64::to_bits));
    assert_eq!(
        expected.external_work_j.map(f64::to_bits),
        result.external_work_j.map(f64::to_bits)
    );
    let (mut next, mut voltage) = ([91.0; 5], [73.0; 4]);
    let mut seen = 0;
    assert!(matches!(
        switched.step_controlled(&initial(), &source, &mut next, &mut voltage, || {
            seen += 1;
            seen == polls
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [91.0; 5]);
    assert_eq!(voltage, [73.0; 4]);
    assert!(
        switched
            .step_interval_sources_controlled(&initial(), &[], &mut next, &mut voltage, || false)
            .is_err()
    );
    assert_eq!(next, [91.0; 5]);
    assert_eq!(voltage, [73.0; 4]);
    switched
        .step_controlled(&initial(), &source, &mut next, &mut voltage, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), direct.map(f64::to_bits));
    assert_eq!(voltage.map(f64::to_bits), direct_mean.map(f64::to_bits));
}
