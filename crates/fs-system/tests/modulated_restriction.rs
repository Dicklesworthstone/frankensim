//! G0/G3: a signal opens a native stream law without creating transported inventory.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortSpec, PortTangent, Proposal, Receipt, StreamQuantity};
use fs_qty::chemistry::SpeciesId;
use fs_system::elements::{GasContents, GasVolume, Restriction};
use fs_system::SystemBuilder;

const G: f64 = 2.0e-8;
const R: f64 = 287.0;
const DT: f64 = 0.125;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey { seed: 53, kernel_id: 11, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn gas() -> GasVolume {
    GasVolume::identified(
        R,
        1.4,
        ["N2", "O2"]
            .into_iter()
            .map(|name| StreamQuantity::SpeciesMass(SpeciesId::new(name).unwrap()))
            .collect(),
    )
    .unwrap()
}

fn point(reverse: bool, opening: f64) -> [f64; 9] {
    let (a, b) = if reverse { (1.0e5, 2.0e5) } else { (2.0e5, 1.0e5) };
    [a, 3.5e5, 0.2, 0.7, b, 4.2e5, 0.9, 0.1, opening]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-11 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

fn evaluate(valve: &Restriction, inputs: &[f64], cx: &Cx<'_>) -> (Receipt, Vec<f64>) {
    let mut outputs = vec![f64::NAN; valve.ports().iter().map(PortSpec::output_width).sum()];
    match valve.propose(&(), DT, inputs, &mut (), &mut outputs, cx).unwrap() {
        Proposal::Step(receipt) => (receipt, outputs),
        _ => panic!("a restriction does not shorten a step"),
    }
}

fn identity(valve: &Restriction) -> Vec<u64> {
    let mut words = Vec::new();
    valve.identity_words(&mut words);
    words
}

/// Differentiate the authored physical law independently of either native tangent.
fn reference(inputs: &[f64; 9], direction: &[f64; 9]) -> ([f64; 8], [f64; 8]) {
    let delta = inputs[0] - inputs[4];
    let upstream = if delta >= 0.0 { 0 } else { 4 };
    let rate = G * inputs[8] * delta;
    let d_rate = G * (inputs[8] * (direction[0] - direction[4]) + delta * direction[8]);
    let (mut values, mut slopes) = ([0.0; 8], [0.0; 8]);
    for coordinate in 0..4 {
        let (value, slope) = if coordinate == 0 {
            (rate, d_rate)
        } else {
            (
                rate * inputs[upstream + coordinate],
                d_rate * inputs[upstream + coordinate] + rate * direction[upstream + coordinate],
            )
        };
        values[coordinate] = value;
        values[4 + coordinate] = -value;
        slopes[coordinate] = slope;
        slopes[4 + coordinate] = -slope;
    }
    (values, slopes)
}

#[test]
fn partial_opening_preserves_the_donor_law_and_has_a_complete_port_tangent() {
    let valve = Restriction::for_gas(G, &gas()).unwrap().modulated();
    assert_eq!(valve.ports().len(), 3);
    assert_eq!(valve.ports()[2].id().as_str(), "opening");
    assert_eq!(valve.ports()[2].input_width(), 1);
    assert_eq!(valve.ports()[2].output_width(), 0);
    let direction = [100.0, -2300.0, 0.01, -0.03, -200.0, 1200.0, -0.01, 0.02, 0.2];
    let covector = [2.0, -0.5, 3.0, 4.0, -1.0, 0.25, -2.0, 0.75];
    with_cx(|cx| {
        for reverse in [false, true] {
            let inputs = point(reverse, 0.4);
            let case = Case { state: &(), stale: &(), dt_s: DT, inputs: &inputs };
            let Observed::Stepped { outputs, receipt, .. } =
                check_component(&valve, &case, cx).unwrap()
            else {
                panic!("restriction shortened");
            };
            assert_eq!(outputs.len(), 8);
            check_port_tangent(&valve, &case, 1.0e-6, 1.0e-9, cx).unwrap();
            let (expected, expected_push) = reference(&inputs, &direction);
            let (mut push, mut pull) = ([f64::NAN; 8], [f64::NAN; 9]);
            valve.jvp(&(), DT, &inputs, &direction, &mut push, cx).unwrap();
            valve.vjp(&(), DT, &inputs, &covector, &mut pull, cx).unwrap();
            for i in 0..8 {
                close(outputs[i], expected[i]);
                close(push[i], expected_push[i]);
            }
            let left: f64 = push.iter().zip(covector).map(|(a, b)| a * b).sum();
            let right: f64 = direction.iter().zip(pull).map(|(a, b)| a * b).sum();
            close(left, right);
            // An ideal-gas throttle keeps h constant; its entropy rise is R ln(p_up/p_down).
            let entropy = DT * expected[0].abs() * R * 2.0_f64.ln();
            close(receipt.entropy_generated_j_per_k, entropy);
            assert_eq!(receipt.stored_change_j, 0.0);
            assert_eq!(receipt.dissipated_j, 0.0);
        }
    });
}

#[test]
fn fully_open_matches_the_original_restriction_and_preserves_stream_identity() {
    let gas = gas();
    let plain = Restriction::for_gas(G, &gas).unwrap();
    let valve = plain.clone().modulated();
    let repeated = valve.clone().modulated();
    assert_eq!(repeated.ports().len(), 3);
    assert_eq!(identity(&repeated), identity(&valve));
    assert_ne!(identity(&plain), identity(&valve));
    for port in &valve.ports()[..2] {
        assert_eq!(port.stream_basis(), gas.ports()[0].stream_basis());
        assert_eq!(port.carried_quantities(), 2);
    }
    let mut words = Vec::new();
    assert!(valve.encode_state(&(), &mut words));
    assert!(words.is_empty());
    assert_eq!(valve.decode_state(&words), Some(()));
    assert_eq!(valve.decode_state(&[0]), None);
    with_cx(|cx| {
        for reverse in [false, true] {
            let inputs = point(reverse, 1.0);
            let (ordinary, original_outputs) = evaluate(&plain, &inputs[..8], cx);
            let (controlled, controlled_outputs) = evaluate(&valve, &inputs, cx);
            assert_eq!(
                original_outputs.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                controlled_outputs.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
            assert_eq!(
                ordinary.entropy_generated_j_per_k.to_bits(),
                controlled.entropy_generated_j_per_k.to_bits()
            );
            assert_eq!(ordinary.energy_tolerance_j.to_bits(), controlled.energy_tolerance_j.to_bits());
        }
    });
}

#[test]
fn closed_valve_has_zero_transport_and_the_correct_one_sided_opening_derivative() {
    let valve = Restriction::for_gas(G, &gas()).unwrap().modulated();
    let mut direction = [0.0; 9];
    direction[8] = 1.0;
    let covector = [2.0, -0.5, 3.0, 4.0, -1.0, 0.25, -2.0, 0.75];
    with_cx(|cx| {
        for reverse in [false, true] {
            let inputs = point(reverse, 0.0);
            let (receipt, outputs) = evaluate(&valve, &inputs, cx);
            assert!(outputs.iter().all(|value| *value == 0.0));
            assert_eq!(receipt.entropy_generated_j_per_k, 0.0);
            let (_, expected) = reference(&inputs, &direction);
            let (mut push, mut pull) = ([f64::NAN; 8], [f64::NAN; 9]);
            valve.jvp(&(), DT, &inputs, &direction, &mut push, cx).unwrap();
            valve.vjp(&(), DT, &inputs, &covector, &mut pull, cx).unwrap();
            let opened = point(reverse, 1.0e-6);
            let one_sided = evaluate(&valve, &opened, cx).1;
            for i in 0..8 {
                close(push[i], expected[i]);
                close(push[i], one_sided[i] / 1.0e-6);
            }
            assert!(pull[..8].iter().all(|value| *value == 0.0));
            close(pull[8], expected.iter().zip(covector).map(|(a, b)| a * b).sum());
        }
    });
}

#[test]
fn modulated_calls_refuse_invalid_commands_layouts_steps_and_cancellation() {
    let valve = Restriction::new(G, 2).unwrap().modulated();
    with_cx(|cx| {
        let inputs = point(false, 0.5);
        for opening in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            let mut bad = inputs;
            bad[8] = opening;
            assert!(valve.propose(&(), DT, &bad, &mut (), &mut [0.0; 8], cx).is_err());
        }
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(valve.propose(&(), dt, &inputs, &mut (), &mut [0.0; 8], cx).is_err());
        }
        for length in [0, 8, 10] {
            assert!(valve.propose(&(), DT, &vec![0.5; length], &mut (), &mut [0.0; 8], cx).is_err());
        }
        for length in [0, 7, 9] {
            assert!(valve.propose(&(), DT, &inputs, &mut (), &mut vec![0.0; length], cx).is_err());
        }
        let mut bad = inputs;
        bad[1] = f64::NAN;
        assert!(valve.propose(&(), DT, &bad, &mut (), &mut [0.0; 8], cx).is_err());
        assert!(valve.jvp(&(), DT, &inputs, &[0.0; 8], &mut [0.0; 8], cx).is_err());
        assert!(valve.jvp(&(), DT, &inputs, &[0.0; 9], &mut [0.0; 7], cx).is_err());
        assert!(valve.vjp(&(), DT, &inputs, &[0.0; 7], &mut [0.0; 9], cx).is_err());
        assert!(valve.vjp(&(), DT, &inputs, &[0.0; 8], &mut [0.0; 8], cx).is_err());
        let mut direction = [0.0; 9];
        direction[8] = f64::INFINITY;
        assert!(valve.jvp(&(), DT, &inputs, &direction, &mut [0.0; 8], cx).is_err());
        cx.cancel_gate().request();
        assert!(valve.propose(&(), DT, &inputs, &mut (), &mut [0.0; 8], cx).is_err());
        assert!(valve.jvp(&(), DT, &inputs, &[0.0; 9], &mut [0.0; 8], cx).is_err());
        assert!(valve.vjp(&(), DT, &inputs, &[0.0; 8], &mut [0.0; 9], cx).is_err());
    });
}

#[test]
fn scheduled_valve_conserves_two_vessels_and_replays_from_checkpoint() {
    let gas = gas();
    let start = [
        gas.contents(5.0e5, 350.0, 0.01, &[0.7, 0.3]).unwrap(),
        gas.contents(1.0e5, 290.0, 0.03, &[0.1, 0.9]).unwrap(),
    ];
    let total_mass: f64 = start.iter().map(|s| s.mass_kg).sum();
    let total_energy: f64 = start.iter().map(|s| s.energy_j).sum();
    let total_carried: [f64; 2] =
        std::array::from_fn(|i| start.iter().map(|s| s.carried_kg[i]).sum());
    let mut builder = SystemBuilder::new();
    let valve = builder
        .add_with_tangent("valve", Restriction::for_gas(G, &gas).unwrap().modulated(), ())
        .unwrap();
    let opening = builder.port(valve, "opening").unwrap();
    builder
        .scheduled_boundary(
            opening,
            Box::new(|time, _, values| {
                values[0] = if time < 0.25 { 0.0 } else if time < 0.5 { 1.0 } else { 0.3 };
            }),
            vec![0.25, 0.5],
        )
        .unwrap();
    let mut nodes = Vec::new();
    for (index, name) in ["a", "b"].into_iter().enumerate() {
        let node = builder.add_with_tangent(name, gas.clone(), start[index].clone()).unwrap();
        builder.bond(builder.port(node, "stream").unwrap(), builder.port(valve, name).unwrap()).unwrap();
        builder.boundary(builder.port(node, "heat").unwrap(), Box::new(|_, _, values| values[0] = 0.0)).unwrap();
        nodes.push(node);
    }
    let mut system = builder.seal().unwrap();
    let observe = |system: &fs_system::System| {
        (
            system.time_s().to_bits(),
            system.steps(),
            nodes.iter().map(|&node| system.state::<GasContents>(node).unwrap().clone()).collect::<Vec<_>>(),
            system.ledger().clone(),
        )
    };
    with_cx(|cx| {
        let mut checkpoint = None;
        for step in 1..=32 {
            let report = system.step(1.0 / 32.0, cx).unwrap();
            assert_eq!(report.dt_s, 1.0 / 32.0);
            let states: Vec<_> = nodes.iter().map(|&node| system.state::<GasContents>(node).unwrap()).collect();
            if step <= 8 {
                assert_eq!(states[0], &start[0]);
                assert_eq!(states[1], &start[1]);
                assert_eq!(system.last_output(system.port(valve, "a").unwrap())[0], 0.0);
            } else {
                assert!(system.last_output(system.port(valve, "a").unwrap())[0] > 0.0);
            }
            let mass: f64 = states.iter().map(|s| s.mass_kg).sum();
            let energy: f64 = states.iter().map(|s| s.energy_j).sum();
            assert!((mass - total_mass).abs() < 2.0e-14);
            assert!((energy - total_energy).abs() < 1.0e-8);
            for i in 0..2 {
                let held: f64 = states.iter().map(|s| s.carried_kg[i]).sum();
                assert!((held - total_carried[i]).abs() < 2.0e-14);
                assert!(system.ledger().carried_residual(i).unwrap().abs() < 2.0e-14);
            }
            assert!(system.ledger().energy_residual_j().abs() < 1.0e-8);
            if step == 13 {
                checkpoint = Some(system.checkpoint().unwrap());
            }
        }
        assert!(system.state::<GasContents>(nodes[0]).unwrap().mass_kg < start[0].mass_kg);
        assert!(system.state::<GasContents>(nodes[1]).unwrap().mass_kg > start[1].mass_kg);
        let expected = observe(&system);
        system.restore(&checkpoint.unwrap()).unwrap();
        for _ in 13..32 {
            system.step(1.0 / 32.0, cx).unwrap();
        }
        assert_eq!(observe(&system), expected);
    });
}
