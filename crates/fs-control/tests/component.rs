//! Sampled policy, actual kernel storage and an independent digital recurrence.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_control::ControlError;
use fs_control::component::{SampledSpeedController, SampledSpeedState};
use fs_control::speed::{ShaftPowerLimits, SpeedControllerParameters};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, check_component, check_port_tangent};
use fs_ports::component::{Component, Proposal};
use fs_system::elements::{ControlledSource, Force, Mass, Motion, SourceAccount, Waveform};
use fs_system::{AdaptiveConfig, AdaptiveError, NodeId, System, SystemBuilder};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn parameters() -> SpeedControllerParameters {
    SpeedControllerParameters {
        sample_period_s: 0.125,
        proportional_gain_nm_s_per_rad: 1.0,
        integral_gain_nm_per_rad: 2.0,
        anti_windup_time_s: 0.25,
    }
}

fn limits() -> ShaftPowerLimits {
    ShaftPowerLimits {
        positive_torque_limit_nm: 3.0,
        negative_torque_limit_nm: 2.0,
        motoring_power_w: 100.0,
        regeneration_power_w: 100.0,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn controller_and_sampled_measurement_pass_component_and_tangent_conformance() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let controller = SampledSpeedController::new(
            parameters(),
            ShaftPowerLimits {
                motoring_power_w: 0.6,
                regeneration_power_w: 0.3,
                ..limits()
            },
        )
        .unwrap();
        let state = SampledSpeedState::default();
        let stale = SampledSpeedState::new(4.0).unwrap();
        // Interior, constant torque limits, and power limits in both directions.
        for inputs in [
            [0.1, 0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0, 0.0],
            [-100.0, 0.0, 0.0, 0.0],
            [100.0, 0.0, 2.0, 0.0],
            [-100.0, 0.0, 2.0, 0.0],
            [100.0, 0.0, -2.0, 0.0],
            [-100.0, 0.0, -2.0, 0.0],
        ] {
            let case = Case {
                state: &state,
                stale: &stale,
                dt_s: 0.0625,
                inputs: &inputs,
            };
            check_component(&controller, &case, cx).unwrap();
            check_port_tangent(&controller, &case, 1e-6, 1e-8, cx).unwrap();
        }
        let mass = Mass::in_domain(PortKind::RotationalTorqueAngularVelocity, 2.0)
            .unwrap()
            .publishing()
            .publishing();
        assert_eq!(mass.ports().len(), 3);
        let old = Motion {
            position_m: 0.3,
            velocity_m_s: 1.2,
        };
        let stale = Motion {
            position_m: -3.0,
            velocity_m_s: -4.0,
        };
        let case = Case {
            state: &old,
            stale: &stale,
            dt_s: 0.125,
            inputs: &[1.5],
        };
        check_component(&mass, &case, cx).unwrap();
        check_port_tangent(&mass, &case, 1e-7, 1e-9, cx).unwrap();
        let mut next = stale;
        let mut outputs = [f64::NAN; 3];
        mass.propose(&old, 0.125, &[1.5], &mut next, &mut outputs, cx)
            .unwrap();
        assert_eq!(&outputs[1..], &[0.3, 1.2]);
        assert_ne!(
            outputs[0], outputs[2],
            "power mean must differ from sampled speed"
        );
    });
}

#[test]
fn plant_substeps_hold_one_sample_and_checkpoint_codec_retains_the_phase() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let controller = SampledSpeedController::new(parameters(), limits()).unwrap();
        let initial = SampledSpeedState::default();
        let mut held = initial;
        let mut output = [f64::NAN];
        controller
            .propose(
                &initial,
                0.03125,
                &[2.0, 0.0, 0.0, 0.0],
                &mut held,
                &mut output,
                cx,
            )
            .unwrap();
        assert_eq!(output, [2.5]);
        assert_eq!(held.integral_torque_nm(), 0.5);
        assert_eq!(held.samples(), 1);
        let mut words = Vec::new();
        assert!(controller.encode_state(&held, &mut words));
        assert_eq!(controller.decode_state(&words), Some(held));
        let case = Case {
            state: &held,
            stale: &initial,
            dt_s: 0.03125,
            inputs: &[-100.0, 1.0, 90.0, -20.0],
        };
        check_component(&controller, &case, cx).unwrap();
        check_port_tangent(&controller, &case, 1e-7, 1e-9, cx).unwrap();
        for _ in 0..3 {
            let old = held;
            controller
                .propose(&old, 0.03125, case.inputs, &mut held, &mut output, cx)
                .unwrap();
            assert_eq!(output, [2.5]);
            assert_eq!(held.samples(), 1);
            assert_eq!(held.integral_torque_nm(), 0.5);
        }
        assert_eq!(held.phase_s(), 0.0);
        let old = held;
        controller
            .propose(&old, 0.125, &[0.0; 4], &mut held, &mut output, cx)
            .unwrap();
        assert_eq!(held.samples(), 2);
        assert_eq!(output, [0.5]);
        let mut untouched = initial;
        output = [91.0];
        assert_eq!(
            controller
                .propose(&initial, 0.5, &[0.0; 4], &mut untouched, &mut output, cx)
                .unwrap(),
            Proposal::Shorten { dt_s: 0.125 }
        );
        assert_eq!(untouched, initial);
        assert_eq!(output, [91.0]);
        words[3] = f64::NAN.to_bits();
        assert!(controller.decode_state(&words).is_none());
    });
}

#[test]
fn decimal_sample_periods_respect_the_cap_without_clock_rounding_microsteps() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let controller = SampledSpeedController::new(
            SpeedControllerParameters {
                sample_period_s: 0.1,
                ..parameters()
            },
            limits(),
        )
        .unwrap();
        let mut state = SampledSpeedState::default();
        for sample in 1..=1000 {
            let mut steps = 0;
            loop {
                let dt = controller.limit_dt_s(&state, 0.01).unwrap();
                assert!(
                    dt > 0.004 && dt <= 0.01,
                    "unresolvable or over-cap interval: {dt}"
                );
                let old = state;
                let mut output = [f64::NAN];
                assert!(matches!(
                    controller
                        .propose(&old, dt, &[0.0; 4], &mut state, &mut output, cx)
                        .unwrap(),
                    Proposal::Step(_)
                ));
                steps += 1;
                assert!(steps <= 11);
                if state.phase_s() == 0.0 {
                    break;
                }
            }
            assert_eq!(state.samples(), sample);
        }
    });
}

fn graph(reference: Rc<Cell<f64>>) -> (System, NodeId, NodeId, NodeId) {
    let kind = PortKind::RotationalTorqueAngularVelocity;
    let mut builder = SystemBuilder::new();
    let rotor = builder
        .add_with_tangent(
            "rotor",
            Mass::in_domain(kind, 1.0).unwrap().publishing(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let control = builder
        .add_with_tangent(
            "control",
            SampledSpeedController::new(parameters(), limits()).unwrap(),
            SampledSpeedState::default(),
        )
        .unwrap();
    let source = builder
        .add_with_tangent(
            "drive",
            ControlledSource::effort(kind),
            SourceAccount::default(),
        )
        .unwrap();
    let load = builder
        .add_with_tangent(
            "load",
            Force::in_domain(kind, Waveform::constant(-0.25).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .common_flow(
            builder.port(rotor, "p").unwrap(),
            &[
                builder.port(source, "p").unwrap(),
                builder.port(load, "p").unwrap(),
            ],
        )
        .unwrap();
    builder
        .signal(
            builder.port(rotor, "velocity").unwrap(),
            builder.port(control, "rotor_speed").unwrap(),
        )
        .unwrap();
    builder
        .signal(
            builder.port(control, "torque").unwrap(),
            builder.port(source, "command").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(control, "reference").unwrap(),
            Box::new(move |_, _, values| values[0] = reference.get()),
        )
        .unwrap();
    for name in ["feedforward", "stator_speed"] {
        builder
            .boundary(
                builder.port(control, name).unwrap(),
                Box::new(|_, _, values| values[0] = 0.0),
            )
            .unwrap();
    }
    (builder.seal().unwrap(), rotor, control, source)
}

#[test]
fn closed_loop_matches_digital_pi_and_physical_inertia_for_each_sample() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for substeps in [1, 4] {
            let reference = Rc::new(Cell::new(10.0));
            let (mut system, rotor, control, source) = graph(reference.clone());
            let (mut speed, mut position, mut integral) = (0.0_f64, 0.0, 0.0);
            for sample in 0..120 {
                let target = if sample < 16 {
                    10.0
                } else if sample < 64 {
                    2.0
                } else {
                    -1.0
                };
                reference.set(target);
                let error = target - speed;
                let predicted = integral + 0.25 * error;
                let unlimited = error + predicted;
                let torque = unlimited.clamp(-2.0, 3.0);
                integral = predicted + (torque - unlimited) / 3.0;
                let next_speed = speed + 0.125 * (torque - 0.25);
                position += 0.125 * 0.5 * (speed + next_speed);
                speed = next_speed;
                for _ in 0..substeps {
                    // Even a caller asking for several periods gets one held command.
                    let cap = if substeps == 1 {
                        1.0
                    } else {
                        0.125 / f64::from(substeps)
                    };
                    system.step(cap, cx).unwrap();
                }
                let state = system.state::<Motion>(rotor).unwrap();
                let policy = system.state::<SampledSpeedState>(control).unwrap();
                near(state.velocity_m_s, speed, 2e-11);
                near(state.position_m, position, 2e-11);
                near(policy.integral_torque_nm(), integral, 2e-11);
                near(policy.held_torque_nm(), torque, 2e-11);
                assert_eq!(policy.samples(), sample + 1);
                assert_eq!(policy.phase_s(), 0.0);
                // Work comes from the physical source and load, never from the PI state.
                let drive = system.state::<SourceAccount>(source).unwrap().energy_in_j;
                near(0.5 * speed * speed + drive + 0.25 * position, 0.0, 5e-10);
            }
        }
    });
}

#[test]
fn rejected_and_cancelled_intervals_do_not_resample_and_rewind_replays_exactly() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let reference = Rc::new(Cell::new(2.0));
        let (mut system, rotor, control, _) = graph(reference.clone());
        system.step(0.03125, cx).unwrap();
        let accepted = *system.state::<SampledSpeedState>(control).unwrap();
        let checkpoint = system.checkpoint().unwrap();
        let ledger = format!("{:?}", system.ledger());
        // Refuse the error estimate after the coarse and refined trials. This
        // proves rollback of sampled policy state, not an adaptive error chart.
        let refused = system.step_adaptive(1.0, AdaptiveConfig::default(), |_, _, _| f64::NAN, cx);
        assert!(matches!(
            refused,
            Err(AdaptiveError::InvalidEstimate { .. })
        ));
        assert_eq!(system.time_s(), checkpoint.time_s());
        assert_eq!(system.steps(), checkpoint.steps());
        assert_eq!(
            system.state::<SampledSpeedState>(control),
            checkpoint.state::<SampledSpeedState>(control)
        );
        assert_eq!(
            system.state::<Motion>(rotor),
            checkpoint.state::<Motion>(rotor)
        );
        assert_eq!(system.ledger(), checkpoint.ledger());
        reference.set(f64::NAN);
        assert!(system.step(0.03125, cx).is_err());
        assert_eq!(
            *system.state::<SampledSpeedState>(control).unwrap(),
            accepted
        );
        assert_eq!(format!("{:?}", system.ledger()), ledger);
        reference.set(2.0);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        let cancelled_cx = Cx::new(
            &cancelled,
            cx.arena(),
            cx.stream_key(),
            cx.budget(),
            cx.mode(),
        );
        assert!(system.step(0.03125, &cancelled_cx).is_err());
        assert_eq!(
            *system.state::<SampledSpeedState>(control).unwrap(),
            accepted
        );
        for _ in 0..39 {
            system.step(0.03125, cx).unwrap();
        }
        let expected_motion = *system.state::<Motion>(rotor).unwrap();
        let expected_control = *system.state::<SampledSpeedState>(control).unwrap();
        let expected_ledger = format!("{:?}", system.ledger());
        system.restore(&checkpoint).unwrap();
        for _ in 0..39 {
            system.step(0.03125, cx).unwrap();
        }
        assert_eq!(*system.state::<Motion>(rotor).unwrap(), expected_motion);
        assert_eq!(
            *system.state::<SampledSpeedState>(control).unwrap(),
            expected_control
        );
        assert_eq!(format!("{:?}", system.ledger()), expected_ledger);
    });
}

#[test]
fn invalid_parameters_and_cancelled_proposals_are_atomic() {
    assert!(
        SampledSpeedController::new(
            SpeedControllerParameters {
                sample_period_s: 0.0,
                ..parameters()
            },
            limits()
        )
        .is_err()
    );
    assert!(
        SampledSpeedController::new(
            parameters(),
            ShaftPowerLimits {
                motoring_power_w: -1.0,
                ..limits()
            }
        )
        .is_err()
    );
    let controller = SampledSpeedController::new(parameters(), limits()).unwrap();
    let old = SampledSpeedState::default();
    let stale = SampledSpeedState::new(7.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for inputs in [vec![], vec![0.0; 3], vec![f64::NAN; 4]] {
            let mut next = stale;
            let mut output = [42.0];
            assert!(
                controller
                    .propose(&old, 0.125, &inputs, &mut next, &mut output, cx)
                    .is_err()
            );
            assert_eq!(next, stale);
            assert_eq!(output, [42.0]);
        }
    });
    gate.request();
    with_cx(&gate, |cx| {
        let mut next = stale;
        let mut output = [42.0];
        assert_eq!(
            controller.propose(&old, 0.125, &[0.0; 4], &mut next, &mut output, cx),
            Err(ControlError::Cancelled)
        );
        assert_eq!(next, stale);
        assert_eq!(output, [42.0]);
    });
}
