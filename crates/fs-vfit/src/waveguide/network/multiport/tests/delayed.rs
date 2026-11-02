//! G1/G3: independent delayed-resistor network and shared-state transaction checks.
use super::*;
use crate::mimo::delayed::{DelayedMimoLoadLimits, PreparedDelayedMimoLoad};
use crate::mimo::DelayedMimoModel;
use crate::passivity::CertificateClass;

fn delayed_spec(delay_samples: f64, model: MimoModel) -> DelayedMimoTerminalSpec {
    DelayedMimoTerminalSpec {
        nodes: vec![2, 3],
        model: DelayedMimoModel {
            delays: vec![0.0, delay_samples * DT, delay_samples * DT, 0.0],
            model,
        },
        limits: DelayedMimoLoadLimits {
            source: MimoLoadLimits {
                band_rad_s: [1.0, core::f64::consts::PI / DT],
                sampled_passivity_points: 129,
                max_states: 32,
                max_memory_bytes: 1 << 20,
            },
            max_delayed_states: 128,
            max_history_values: 256,
        },
    }
}

fn delayed_snapshot(net: &WaveguideNetwork) -> Vec<u64> {
    let mut state = snapshot(net);
    let load = net.delayed_mimo_load().unwrap();
    state.push(load.history_cursor() as u64);
    state.extend(load.input_history().iter().map(|x| x.to_bits()));
    state
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-12 * (1.0 + a.abs() + b.abs()), "{a} != {b}");
}

#[test]
fn fractional_mutual_arrivals_match_an_independent_held_flow_network() {
    // Z_ii is instantaneous; Z_ij = 0.5 exp(-s tau). The exact interval
    // average of a delayed held flow has two overlap lengths, not a chosen
    // interpolation filter. Solve the independent 2x2 system at each interval.
    for delay in [0.0_f64, 0.25, 2.25, 7.0] {
        let mut net = graph().with_delayed_mimo_terminals(delayed_spec(delay, model())).unwrap();
        let whole = delay.floor() as usize;
        let fraction = delay - whole as f64;
        let current = if whole == 0 { 0.5 * (1.0 - fraction) } else { 0.0 };
        let determinant = 5.0 * 7.0 - current * current;
        let mut flows: Vec<[f64; 2]> = Vec::new();
        let mut inlet_work = 0.0;
        let mut exterior_work = 0.0;
        for step in 0..96 {
            let incoming: [f64; 2] = [2, 3].map(|node|
                net.arriving(net.ports[net.offsets[node]]));
            let past = |lag: usize, port: usize| {
                if lag == 0 || step < lag { 0.0 } else { flows[step - lag][port] }
            };
            let history: [f64; 2] = std::array::from_fn(|i| 0.5 * (
                (1.0 - fraction) * past(whole, 1 - i)
                    + fraction * past(whole + 1, 1 - i)));
            let rhs = [2.0 * incoming[0] - history[0], 2.0 * incoming[1] - history[1]];
            let expected = [(7.0 * rhs[0] - current * rhs[1]) / determinant,
                (5.0 * rhs[1] - current * rhs[0]) / determinant];
            let frame = net.step(if step == 0 { 1.0 } else { 0.0 }).unwrap();
            let ports = net.mimo_port_frames().unwrap();
            for i in 0..2 {
                let pressure = [2.0, 3.0][i] * expected[i]
                    + current * expected[1 - i] + history[i];
                close(ports[i].flow_m3_s, expected[i]);
                close(ports[i].pressure_pa, pressure);
                close(ports[i].work_j, pressure * expected[i] * DT);
            }
            if delay == 2.25 && step < 4 { assert_eq!(ports[1].pressure_pa, 0.0); }
            if delay == 2.25 && step == 4 { assert!(ports[1].pressure_pa > 0.0); }
            flows.push([ports[0].flow_m3_s, ports[1].flow_m3_s]);
            inlet_work += frame.inlet_work_j;
            exterior_work += frame.coupled_terminal_work_j;
            assert!(frame.balance_residual_j().abs() < 2e-17);
        }
        near(net.stored_energy_j() + exterior_work, inlet_work);
        let (continuous, sampled) = net.mimo_passivity().unwrap();
        assert!(continuous.passive && sampled.passive);
        assert_eq!(continuous.class, CertificateClass::GridOnly);
        assert_eq!(sampled.class, CertificateClass::GridOnly);
    }
}

#[test]
fn delayed_rational_ports_commit_the_same_poles_and_rings_as_the_original_owner() {
    let mut source = model();
    source.poles = vec![PoleTerm::Real { pole: -400.0, residue: 0.0 }];
    source.residues = vec![vec![C64::from_re(200.0), C64::from_re(50.0),
        C64::from_re(50.0), C64::from_re(180.0)]];
    let spec = delayed_spec(1.375, source);
    let mut reference = PreparedDelayedMimoLoad::new(&spec.model, DT, spec.limits).unwrap();
    let mut net = graph().with_delayed_mimo_terminals(spec).unwrap();
    for sample in 0..96 {
        let input = if sample < 37 { f64::from(sample % 7) * 0.1 - 0.3 } else { 0.0 };
        let before = delayed_snapshot(&net);
        let preview = net.preview_step(input).unwrap();
        assert_eq!(delayed_snapshot(&net), before);
        let mut other = graph();
        assert!(net.step_pair(&mut other, [input, f64::NAN]).is_err());
        assert_eq!(delayed_snapshot(&net), before);
        let accepted = net.step(input).unwrap();
        assert_eq!(accepted, preview);
        let ports = net.mimo_port_frames().unwrap();
        let mut pressure = [0.0; 2];
        reference.step_mean_into(&[ports[0].flow_m3_s, ports[1].flow_m3_s], &mut pressure).unwrap();
        close(pressure[0], ports[0].pressure_pa);
        close(pressure[1], ports[1].pressure_pa);
        let actual = net.delayed_mimo_load().unwrap();
        assert_eq!(actual.pole_states(), reference.pole_states());
        assert_eq!(actual.input_history(), reference.input_history());
        assert_eq!(actual.history_cursor(), reference.history_cursor());
        assert!(accepted.balance_residual_j().abs() < 2e-17);
    }
    let before = delayed_snapshot(&net);
    assert!(net.step(f64::MAX).is_err());
    assert_eq!(delayed_snapshot(&net), before);
}

#[test]
fn zero_delay_resistor_group_preserves_existing_network_results() {
    let mut original = graph().with_mimo_terminals(spec(model())).unwrap();
    let mut delayed = graph().with_delayed_mimo_terminals(delayed_spec(0.0, model())).unwrap();
    for sample in 0..64 {
        let input = if sample < 19 { f64::from(sample % 5) * 0.2 - 0.4 } else { 0.0 };
        assert_eq!(delayed.step(input).unwrap(), original.step(input).unwrap());
        assert_eq!(delayed.mimo_port_frames(), original.mimo_port_frames());
    }
}

#[test]
fn delayed_terminal_identity_limits_and_failure_atomic_attachment() {
    let source = delayed_spec(2.25, model());
    let mut net = graph();
    let mut duplicate = source.clone();
    duplicate.nodes[1] = duplicate.nodes[0];
    assert!(net.install_delayed_mimo_terminals(duplicate).is_err());
    let mut no_history = source.clone();
    no_history.limits.max_history_values = 0;
    assert!(net.install_delayed_mimo_terminals(no_history).is_err());
    let mut no_storage = source.clone();
    no_storage.limits.source.max_memory_bytes = 1;
    assert!(net.install_delayed_mimo_terminals(no_storage).is_err());
    assert!(net.mimo_port_frames().is_none());
    assert_eq!(net.stored_energy_j(), 0.0);
    net.install_delayed_mimo_terminals(source.clone()).unwrap();
    assert!(net.mimo_terminal_spec().is_none());
    let installed = net.delayed_mimo_terminal_spec().unwrap();
    assert_eq!(installed.nodes, source.nodes);
    assert_eq!(installed.model.delays, source.model.delays);
    assert_eq!(installed.model.model, source.model.model);
    assert_eq!(installed.limits.max_history_values, source.limits.max_history_values);
    assert!(net.delayed_mimo_load().unwrap().limits().source.max_memory_bytes
        < source.limits.source.max_memory_bytes);
    assert!(net.install_mimo_terminals(spec(model())).is_err());
    net.step(0.5).unwrap();
    let before = delayed_snapshot(&net);
    assert!(net.install_delayed_mimo_terminals(source).is_err());
    assert_eq!(delayed_snapshot(&net), before);
}
