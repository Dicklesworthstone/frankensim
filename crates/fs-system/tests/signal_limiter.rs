//! G0/G2/G3: bounded controls through the ordinary signal and physical kernel.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortRole, PortTangent};
use fs_system::elements::{
    AffineSignal, ControlledSource, ElementError, Mass, Motion, SignalLimiter, SourceAccount,
};
use fs_system::{NodeId, System, SystemBuilder};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn vector_limits_conform_and_the_interior_tangent_is_its_transpose() {
    let limiter = SignalLimiter::new(vec![0.0, -2.0, 3.0], vec![1.0, 2.0, 3.0]).unwrap();
    assert_eq!(limiter.class(), Class::NonSmooth);
    assert_eq!(limiter.ports().len(), 2);
    assert_eq!(limiter.ports()[0].id().as_str(), "input");
    assert_eq!(limiter.ports()[0].role(), PortRole::SignalIn);
    assert_eq!(limiter.ports()[1].id().as_str(), "output");
    assert_eq!(limiter.ports()[1].role(), PortRole::SignalOut);
    assert!(limiter.ports().iter().all(|port| port.width() == 3));
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (inputs, expected, diagonal) in [
            ([-0.25, 0.0, 5.0], [0.0, 0.0, 3.0], [0.0, 1.0, 0.0]),
            ([0.5, 5.0, 2.0], [0.5, 2.0, 3.0], [1.0, 0.0, 0.0]),
            ([0.25, -1.0, 3.0], [0.25, -1.0, 3.0], [1.0, 1.0, 0.0]),
        ] {
            let case = Case {
                state: &(),
                stale: &(),
                dt_s: 0.125,
                inputs: &inputs,
            };
            let Observed::Stepped {
                outputs, receipt, ..
            } = check_component(&limiter, &case, cx).unwrap()
            else {
                panic!("a limiter does not shorten a step")
            };
            assert_eq!(outputs, expected);
            assert_eq!(receipt.stored_change_j, 0.0);
            assert_eq!(receipt.dissipated_j, 0.0);
            assert_eq!(receipt.energy_tolerance_j, 0.0);
            check_port_tangent(&limiter, &case, 1.0e-6, 1.0e-8, cx).unwrap();
            let (direction, covector) = ([2.0, -3.0, 7.0], [-4.0, 5.0, 9.0]);
            let (mut push, mut pull) = ([f64::NAN; 3], [f64::NAN; 3]);
            limiter
                .jvp(&(), case.dt_s, &inputs, &direction, &mut push, cx)
                .unwrap();
            limiter
                .vjp(&(), case.dt_s, &inputs, &covector, &mut pull, cx)
                .unwrap();
            for index in 0..3 {
                assert_eq!(push[index], diagonal[index] * direction[index]);
                assert_eq!(pull[index], diagonal[index] * covector[index]);
            }
            near(
                push.iter().zip(covector).map(|(a, b)| a * b).sum(),
                pull.iter().zip(direction).map(|(a, b)| a * b).sum(),
                0.0,
            );
        }
        // At a kink there is no classical derivative: the documented
        // generalized Jacobian selects the held branch, in both directions.
        let mut tangent = [f64::NAN; 3];
        for inputs in [[0.0, -2.0, 3.0], [1.0, 2.0, 3.0]] {
            limiter
                .jvp(&(), 1.0, &inputs, &[1.0; 3], &mut tangent, cx)
                .unwrap();
            assert_eq!(tangent, [0.0; 3]);
            limiter
                .vjp(&(), 1.0, &inputs, &[-1.0; 3], &mut tangent, cx)
                .unwrap();
            assert_eq!(tangent, [0.0; 3]);
        }
    });
}

#[test]
fn bound_identity_and_empty_state_codec_are_explicit() {
    let limiter = SignalLimiter::new(vec![-1.0, 0.0], vec![1.0, 2.0]).unwrap();
    let mut state = Vec::new();
    assert!(limiter.encode_state(&(), &mut state));
    assert_eq!(limiter.decode_state(&state), Some(()));
    assert_eq!(limiter.decode_state(&[0]), None);
    assert!(SignalLimiter::STATE_CHARTER.is_some());
    let mut identity = Vec::new();
    limiter.identity_words(&mut identity);
    for changed in [
        SignalLimiter::new(vec![-2.0, 0.0], vec![1.0, 2.0]).unwrap(),
        SignalLimiter::new(vec![-1.0, 0.0], vec![1.0, 3.0]).unwrap(),
        SignalLimiter::new(vec![-1.0], vec![1.0]).unwrap(),
    ] {
        let mut other = Vec::new();
        changed.identity_words(&mut other);
        assert_ne!(identity, other);
    }
}

#[test]
fn invalid_vectors_and_cancelled_calls_refuse_before_writing() {
    for (lower, upper) in [
        (vec![], vec![]),
        (vec![0.0], vec![1.0, 2.0]),
        (vec![2.0], vec![1.0]),
        (vec![f64::NAN], vec![1.0]),
        (vec![0.0], vec![f64::INFINITY]),
        (vec![f64::NEG_INFINITY], vec![1.0]),
    ] {
        assert!(SignalLimiter::new(lower, upper).is_err());
    }
    let limiter = SignalLimiter::new(vec![-1.0; 2], vec![1.0; 2]).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (dt, inputs) in [
            (0.0, &[0.0, 0.0][..]),
            (-1.0, &[0.0, 0.0][..]),
            (f64::NAN, &[0.0, 0.0][..]),
            (f64::INFINITY, &[0.0, 0.0][..]),
            (1.0, &[0.0][..]),
            (1.0, &[0.0, f64::NAN][..]),
            (1.0, &[0.0, f64::INFINITY][..]),
        ] {
            let mut output = [91.0; 2];
            assert!(
                limiter
                    .propose(&(), dt, inputs, &mut (), &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [91.0; 2]);
        }
        assert!(
            limiter
                .propose(&(), 1.0, &[0.0; 2], &mut (), &mut [91.0], cx)
                .is_err()
        );
        for seeds in [&[0.0][..], &[0.0, f64::NAN][..], &[0.0, f64::INFINITY][..]] {
            let mut output = [91.0; 2];
            assert!(
                limiter
                    .jvp(&(), 1.0, &[0.0; 2], seeds, &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [91.0; 2]);
            assert!(
                limiter
                    .vjp(&(), 1.0, &[0.0; 2], seeds, &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [91.0; 2]);
        }
        // Finite inputs may be arbitrarily outside the interval; no
        // subtraction or multiplication is needed to bound them.
        let mut output = [91.0; 2];
        limiter
            .propose(&(), 1.0, &[f64::MIN, f64::MAX], &mut (), &mut output, cx)
            .unwrap();
        assert_eq!(output, [-1.0, 1.0]);
        cx.cancel_gate().request();
        assert_eq!(
            limiter.propose(&(), 1.0, &[0.0; 2], &mut (), &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            limiter.jvp(&(), 1.0, &[0.0; 2], &[1.0; 2], &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            limiter.vjp(&(), 1.0, &[0.0; 2], &[1.0; 2], &mut output, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(output, [-1.0, 1.0]);
    });
}

fn reference(time: f64) -> f64 {
    if time < 1.0 {
        3.0
    } else if time < 2.0 {
        0.3
    } else {
        -3.0
    }
}

fn feedback_graph() -> (System, NodeId, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent(
            "body",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let source = builder
        .add_with_tangent(
            "actuator",
            ControlledSource::effort(PortKind::MechanicalForceVelocity),
            SourceAccount::default(),
        )
        .unwrap();
    let control = builder
        .add_with_tangent(
            "control",
            AffineSignal::new(2, 1, vec![1.0, -0.5], vec![0.0]).unwrap(),
            (),
        )
        .unwrap();
    let limiter = builder
        .add_with_tangent(
            "limiter",
            SignalLimiter::new(vec![-1.0], vec![1.0]).unwrap(),
            (),
        )
        .unwrap();
    builder
        .bond(
            builder.port(body, "p").unwrap(),
            builder.port(source, "p").unwrap(),
        )
        .unwrap();
    for (from, output, to, input) in [
        (control, "output", limiter, "input"),
        (limiter, "output", control, "input1"),
        (limiter, "output", source, "command"),
    ] {
        builder
            .signal(
                builder.port(from, output).unwrap(),
                builder.port(to, input).unwrap(),
            )
            .unwrap();
    }
    builder
        .scheduled_boundary(
            builder.port(control, "input0").unwrap(),
            Box::new(|start, _, values| values[0] = reference(start)),
            vec![1.0, 2.0],
        )
        .unwrap();
    (builder.seal().unwrap(), body, source, limiter)
}

#[test]
fn saturated_algebraic_feedback_drives_real_work_and_replays_across_active_sets() {
    let (mut system, body, source, limiter) = feedback_graph();
    let output = system.port(limiter, "output").unwrap();
    let mut position = 0.0;
    let mut velocity = 0.0;
    let mut work = 0.0;
    let mut saved = None;
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for step in 0..256 {
            let start = system.time_s();
            // Independently solve y = clamp(reference - y/2, -1, 1).
            let force = if start < 1.0 {
                1.0
            } else if start < 2.0 {
                0.2
            } else {
                -1.0
            };
            let dt = 1.0 / 64.0;
            let next_velocity = velocity + dt * force / 2.0;
            let displacement = 0.5 * dt * (velocity + next_velocity);
            position += displacement;
            velocity = next_velocity;
            work += force * displacement;
            system.step(dt, cx).unwrap();
            near(system.last_output(output)[0], force, 1.0e-12);
            let actual = system.state::<Motion>(body).unwrap();
            near(actual.position_m, position, 1.0e-11);
            near(actual.velocity_m_s, velocity, 1.0e-11);
            near(
                system.state::<SourceAccount>(source).unwrap().energy_in_j,
                -work,
                1.0e-11,
            );
            near(system.ledger().energy_residual_j(), 0.0, 1.0e-11);
            if step == 95 {
                saved = Some(system.checkpoint().unwrap());
            }
        }
        let final_motion = *system.state::<Motion>(body).unwrap();
        let final_source = *system.state::<SourceAccount>(source).unwrap();
        let final_ledger = system.ledger().clone();
        let checkpoint = saved.unwrap();
        system.restore(&checkpoint).unwrap();
        for _ in 96..256 {
            system.step(1.0 / 64.0, cx).unwrap();
        }
        assert_eq!(system.state::<Motion>(body), Some(&final_motion));
        assert_eq!(system.state::<SourceAccount>(source), Some(&final_source));
        assert_eq!(system.ledger(), &final_ledger);
        near(system.time_s(), 4.0, 0.0);
        cx.cancel_gate().request();
        assert!(system.step(1.0 / 64.0, cx).is_err());
        assert_eq!(system.state::<Motion>(body), Some(&final_motion));
        assert_eq!(system.state::<SourceAccount>(source), Some(&final_source));
        assert_eq!(system.ledger(), &final_ledger);
    });
    near(work, velocity * velocity, 1.0e-11);
}
