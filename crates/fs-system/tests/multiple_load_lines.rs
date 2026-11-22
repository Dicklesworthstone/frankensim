//! G0/G4/G5: several non-smooth owners in one coupled block, with the
//! mechanical response and each owner's dissipated heat independently checked.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::elements::{
    Force, HeatCapacity, Mass, Motion, Slider, Slip, SourceAccount, Spring, Waveform,
};
use fs_system::{NodeId, PortRef, SolveConfig, System, SystemBuilder, SystemError};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn port(builder: &SystemBuilder, node: NodeId, name: &str) -> PortRef {
    builder.port(node, name).unwrap()
}

fn grounded(builder: &mut SystemBuilder, name: &str, force: f64, velocity: f64) -> NodeId {
    let initial = if velocity == 0.0 {
        Slip::stuck()
    } else {
        Slip::sliding(-velocity)
    };
    let node = builder
        .add_on_line(name, Slider::new(force, 1.0, 1.0).unwrap(), initial)
        .unwrap();
    builder
        .boundary(port(builder, node, "b"), Box::new(|_, _, out| out[0] = 0.0))
        .unwrap();
    node
}

/// Two independent friction laws load the same body and deliver their losses
/// to one thermal store. Their forces are 2 N and 1.5 N while sliding.
fn two_brakes(velocity: f64, drive: Option<f64>) -> (System, NodeId, [NodeId; 2], NodeId) {
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent(
            "body",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: velocity,
            },
        )
        .unwrap();
    let brakes = [
        grounded(&mut builder, "front", 2.0, velocity),
        grounded(&mut builder, "rear", 1.5, velocity),
    ];
    let heat = builder
        .add_with_tangent("heat", HeatCapacity::new(100.0).unwrap(), 300.0)
        .unwrap();
    let mut loads: Vec<_> = brakes
        .iter()
        .map(|&brake| port(&builder, brake, "a"))
        .collect();
    if let Some(force) = drive {
        let source = builder
            .add_with_tangent(
                "drive",
                Force::new(Waveform::constant(force).unwrap()),
                SourceAccount::default(),
            )
            .unwrap();
        loads.push(port(&builder, source, "p"));
    }
    builder
        .common_flow(port(&builder, body, "p"), &loads)
        .unwrap();
    let warmth: Vec<_> = brakes
        .iter()
        .map(|&brake| port(&builder, brake, "heat"))
        .collect();
    builder
        .common_effort(port(&builder, heat, "heat"), &warmth)
        .unwrap();
    let mut system = builder.seal().unwrap();
    for heat_port in warmth {
        assert!(system.guess(heat_port, &[300.0]));
    }
    (system, body, brakes, heat)
}

fn close(got: f64, want: f64, tolerance: f64) {
    assert!(
        (got - want).abs() <= tolerance,
        "{got:.17e} != {want:.17e} within {tolerance:.1e}"
    );
}

#[test]
fn two_brakes_on_one_body_recover_force_motion_and_heat() {
    let (mut system, body, brakes, heat) = two_brakes(3.0, None);
    assert_eq!(system.plan_shape().1, 1, "both laws must be in one loop");
    with_cx(|cx, _| {
        for step in 1..=50 {
            system.step(0.01, cx).unwrap();
            let time = f64::from(step) * 0.01;
            let expected_velocity = 3.0 - 3.5 / 2.0 * time;
            let expected_position = 3.0 * time - 0.5 * 3.5 / 2.0 * time * time;
            let state = system.state::<Motion>(body).unwrap();
            close(state.velocity_m_s, expected_velocity, 3e-12);
            close(state.position_m, expected_position, 3e-12);
            for (&brake, force) in brakes.iter().zip([2.0, 1.5]) {
                close(
                    system.last_output(system.port(brake, "a").unwrap())[0],
                    -force,
                    2e-12,
                );
                close(
                    system.state::<Slip>(brake).unwrap().slip_m_s(),
                    -expected_velocity,
                    3e-12,
                );
            }
            close(
                100.0 * (system.state::<f64>(heat).unwrap() - 300.0),
                3.5 * expected_position,
                5e-10,
            );
            close(system.ledger().energy_residual_j(), 0.0, 2e-10);
        }
    });
}

#[test]
fn connected_bodies_solve_one_sticking_and_one_sliding_law_together() {
    // The left brake holds against the spring. The right body is therefore
    // a forced oscillator with known midpoint update, while both brake laws
    // and their mutual compliance are found by the generic kernel.
    let mut builder = SystemBuilder::new();
    let left = builder
        .add_with_tangent(
            "left",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let right = builder
        .add_with_tangent(
            "right",
            Mass::new(3.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 1.0,
            },
        )
        .unwrap();
    let brakes = [
        grounded(&mut builder, "holding", 2.0, 0.0),
        grounded(&mut builder, "sliding", 0.1, 1.0),
    ];
    let spring = builder
        .add_with_tangent("spring", Spring::new(4.0, 0.0).unwrap(), 0.0)
        .unwrap();
    for (mass, brake, end) in [(left, brakes[0], "a"), (right, brakes[1], "b")] {
        builder
            .common_flow(
                port(&builder, mass, "p"),
                &[port(&builder, brake, "a"), port(&builder, spring, end)],
            )
            .unwrap();
        builder
            .boundary(
                port(&builder, brake, "heat"),
                Box::new(|_, _, out| out[0] = 300.0),
            )
            .unwrap();
    }
    let mut system = builder.seal().unwrap();
    assert_eq!(system.plan_shape().1, 1);
    let (mut position, mut velocity) = (0.0, 1.0);
    let dt = 0.01;
    with_cx(|cx, _| {
        for _ in 0..20 {
            let next_position = ((1.0 - dt * dt / 3.0) * position + dt * velocity
                - dt * dt * 0.1 / 6.0)
                / (1.0 + dt * dt / 3.0);
            velocity -= dt / 3.0 * (2.0 * (position + next_position) + 0.1);
            position = next_position;
            system.step(dt, cx).unwrap();
            assert!(system.state::<Slip>(brakes[0]).unwrap().is_stuck());
            assert!(!system.state::<Slip>(brakes[1]).unwrap().is_stuck());
            close(
                system.state::<Motion>(left).unwrap().velocity_m_s,
                0.0,
                2e-12,
            );
            close(
                system.state::<Motion>(right).unwrap().position_m,
                position,
                3e-12,
            );
            close(
                system.state::<Motion>(right).unwrap().velocity_m_s,
                velocity,
                3e-12,
            );
            close(system.ledger().energy_residual_j(), 0.0, 2e-10);
        }
    });
}

/// Accepted observable words, including every wire of the brake laws. A
/// discarded trial may remake matrices, but must not move any of these.
fn accepted_words(system: &System, body: NodeId, brakes: [NodeId; 2], heat: NodeId) -> Vec<u64> {
    let body_state = system.state::<Motion>(body).unwrap();
    let mut words = vec![
        system.time_s().to_bits(),
        system.steps(),
        body_state.position_m.to_bits(),
        body_state.velocity_m_s.to_bits(),
        system.state::<f64>(heat).unwrap().to_bits(),
    ];
    for brake in brakes {
        let state = system.state::<Slip>(brake).unwrap();
        words.extend([state.slip_m_s().to_bits(), u64::from(state.is_stuck())]);
        for name in ["a", "b", "heat"] {
            let port = system.port(brake, name).unwrap();
            words.extend(
                system
                    .last_input(port)
                    .iter()
                    .chain(system.last_output(port))
                    .map(|value| value.to_bits()),
            );
        }
    }
    words
}

#[test]
fn failed_iteration_cancellation_and_checkpoint_replay_leave_every_owner_atomic() {
    let (mut system, body, brakes, heat) = two_brakes(3.0, None);
    with_cx(|cx, _| {
        for _ in 0..5 {
            system.step(0.01, cx).unwrap();
        }
    });
    let checkpoint = system.checkpoint().unwrap();
    let before = accepted_words(&system, body, brakes, heat);
    let accounts = system.ledger().clone();
    system.set_config(SolveConfig {
        max_iterations: 0,
        ..SolveConfig::default()
    });
    with_cx(|cx, _| {
        assert!(matches!(
            system.step(0.02, cx),
            Err(SystemError::NotConverged { .. })
        ))
    });
    assert_eq!(accepted_words(&system, body, brakes, heat), before);
    assert_eq!(system.ledger(), &accounts);
    system.restore(&checkpoint).unwrap();
    with_cx(|cx, gate| {
        gate.request();
        assert!(matches!(system.step(0.02, cx), Err(SystemError::Cancelled)));
    });
    assert_eq!(accepted_words(&system, body, brakes, heat), before);
    assert_eq!(system.ledger(), &accounts);
    let steps = [0.02, 0.01, 0.005, 0.015];
    with_cx(|cx, _| {
        for dt in steps {
            system.step(dt, cx).unwrap();
        }
    });
    let uninterrupted = accepted_words(&system, body, brakes, heat);
    let accounts = system.ledger().clone();
    system.restore(&checkpoint).unwrap();
    with_cx(|cx, _| {
        for dt in steps {
            system.step(dt, cx).unwrap();
        }
    });
    assert_eq!(accepted_words(&system, body, brakes, heat), uninterrupted);
    assert_eq!(system.ledger(), &accounts);
}

#[test]
fn redundant_loaded_stick_constraints_refuse_without_publishing() {
    // Two holds constrain the same velocity, leaving their force split
    // undetermined. This driver does not invent a force-sharing law.
    let (mut system, body, brakes, heat) = two_brakes(0.0, Some(1.0));
    let before = accepted_words(&system, body, brakes, heat);
    let accounts = system.ledger().clone();
    with_cx(|cx, _| {
        let error = system.step(0.125, cx).unwrap_err();
        match error {
            SystemError::Singular { nodes } => {
                assert!(nodes.iter().any(|name| name == "front"));
                assert!(nodes.iter().any(|name| name == "rear"));
            }
            other => panic!("expected a named singular block: {other:?}"),
        }
    });
    assert_eq!(accepted_words(&system, body, brakes, heat), before);
    assert_eq!(system.ledger(), &accounts);
}
