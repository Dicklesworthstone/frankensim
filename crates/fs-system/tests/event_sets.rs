//! Earliest competing thresholds share one search and one accepted coupled step.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::elements::{Mass, Motion, Spring};
use fs_system::{
    EventConfig, EventDirection, EventSetConfig, EventSetError, EventSetReport, EventSetStatus,
    EventStatus, Ledger, LocatedEvent, NodeId, PortRef, SolveConfig, System, SystemBuilder,
    SystemError,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 47,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn config() -> EventSetConfig {
    EventSetConfig {
        time_tolerance_s: 1.0e-10,
        max_bisections: 64,
    }
}

fn located(report: &EventSetReport) -> &[LocatedEvent] {
    let EventSetStatus::Located(events) = &report.status else {
        panic!("{report:?}")
    };
    events
}

struct Device {
    mass: NodeId,
    spring: Option<NodeId>,
    ports: Vec<PortRef>,
}

fn mechanical(oscillator: bool) -> (System, Device) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(1.0).unwrap(),
            if oscillator {
                Motion {
                    position_m: 1.0,
                    velocity_m_s: 0.0,
                }
            } else {
                Motion {
                    position_m: 1.37,
                    velocity_m_s: -0.24,
                }
            },
        )
        .unwrap();
    let p = builder.port(mass, "p").unwrap();
    let (spring, ports) = if oscillator {
        let spring = builder
            .add_with_tangent("spring", Spring::new(4.0, 0.0).unwrap(), -1.0)
            .unwrap();
        let a = builder.port(spring, "a").unwrap();
        let b = builder.port(spring, "b").unwrap();
        builder.bond(p, a).unwrap();
        builder
            .boundary(b, Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
        (Some(spring), vec![p, a, b])
    } else {
        builder
            .boundary(p, Box::new(|_, _, values| values[0] = -9.81))
            .unwrap();
        (None, vec![p])
    };
    (
        builder.seal().unwrap(),
        Device {
            mass,
            spring,
            ports,
        },
    )
}

#[derive(Debug, PartialEq)]
struct Observed {
    time_s: f64,
    steps: u64,
    motion: [f64; 2],
    spring: Option<f64>,
    wires: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)>,
    ledger: Ledger,
    config: SolveConfig,
}

fn observe(system: &System, device: &Device) -> Observed {
    let motion = system.state::<Motion>(device.mass).unwrap();
    Observed {
        time_s: system.time_s(),
        steps: system.steps(),
        motion: [motion.position_m, motion.velocity_m_s],
        spring: device
            .spring
            .map(|node| *system.state::<f64>(node).unwrap()),
        wires: device
            .ports
            .iter()
            .map(|&port| {
                (
                    system.last_input(port).to_vec(),
                    system.last_output(port).to_vec(),
                    system.last_allowance(port).to_vec(),
                )
            })
            .collect(),
        ledger: system.ledger().clone(),
        config: system.checkpoint().unwrap().config(),
    }
}

#[test]
fn earliest_physical_threshold_is_independent_of_declaration_order() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut previous = None;
        for order in [[0, 1, 2], [2, 1, 0], [1, 0, 2]] {
            let (mut system, device) = mechanical(false);
            let initial = system.checkpoint().unwrap();
            let directions = order.map(|i| {
                if i == 2 {
                    EventDirection::Rising // The earlier falling height is deliberately excluded.
                } else {
                    EventDirection::Falling
                }
            });
            let report = system
                .step_to_events(
                    1.0,
                    config(),
                    &directions,
                    |s, values| {
                        let motion = s.state::<Motion>(device.mass).unwrap();
                        let guards = [
                            motion.position_m,
                            motion.velocity_m_s + 2.0,
                            motion.position_m - 1.3,
                        ];
                        for (out, &index) in values.iter_mut().zip(&order) {
                            *out = guards[index];
                        }
                    },
                    cx,
                )
                .unwrap();
            let events = located(&report);
            assert_eq!(events.len(), 1);
            assert_eq!(order[events[0].index], 1);
            let bracket = events[0].bracket;
            let expected = (2.0 - 0.24) / 9.81;
            assert!(bracket.lower_dt_s <= expected && expected <= bracket.upper_dt_s);
            assert!(bracket.upper_dt_s - bracket.lower_dt_s <= config().time_tolerance_s);
            assert!(bracket.lower_value >= 0.0 && bracket.upper_value <= 0.0);
            assert_eq!(report.values[events[0].index], bracket.upper_value);
            assert_eq!(system.steps(), 1);
            assert_eq!(system.ledger().steps, 1);
            let accepted = observe(&system, &device);
            if let Some((state, step, trials)) = &previous {
                assert_eq!(&accepted, state);
                assert_eq!(&report.step, step);
                assert_eq!(&report.trial_steps, trials);
            }
            system.restore(&initial).unwrap();
            system.step(report.step.unwrap().dt_s, cx).unwrap();
            assert_eq!(observe(&system, &device), accepted);
            previous = Some((accepted, report.step, report.trial_steps));
        }
    });
}

#[test]
fn exact_later_roots_never_hide_an_earlier_crossing() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        // The later root is exact at either the full trial or first midpoint.
        for roots in [[0.25, 1.0], [1.0, 0.25], [0.125, 0.5], [0.5, 0.125]] {
            let mut system = SystemBuilder::new().seal().unwrap();
            let report = system
                .step_to_events(
                    1.0,
                    config(),
                    &[EventDirection::Rising; 2],
                    |s, values| {
                        for (value, root) in values.iter_mut().zip(roots) {
                            *value = s.time_s() - root;
                        }
                    },
                    cx,
                )
                .unwrap();
            let expected_index = usize::from(roots[1] < roots[0]);
            let events = located(&report);
            assert_eq!(events.len(), 1);
            assert_eq!(events[0].index, expected_index);
            let expected = roots[expected_index];
            assert_eq!(report.time_s, expected);
            assert_eq!(events[0].bracket.lower_dt_s, expected);
            assert_eq!(events[0].bracket.upper_dt_s, expected);
            assert_eq!(events[0].bracket.lower_value, 0.0);
            assert_eq!(events[0].bracket.upper_value, 0.0);
            assert_eq!(report.trial_steps, if expected == 0.25 { 3 } else { 4 });
            assert_eq!(system.steps(), 1);
        }
    });
}

#[test]
fn unresolved_crossings_share_a_bracket_and_have_stable_index_order() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for roots in [[0.502, 0.501, 0.508], [0.5, 0.5, 0.508]] {
            let mut system = SystemBuilder::new().seal().unwrap();
            let limits = EventSetConfig {
                time_tolerance_s: 0.01,
                ..config()
            };
            let report = system
                .step_to_events(
                    1.0,
                    limits,
                    &[EventDirection::Rising; 3],
                    |s, values| {
                        for (value, root) in values.iter_mut().zip(roots) {
                            *value = s.time_s() - root;
                        }
                    },
                    cx,
                )
                .unwrap();
            let events = located(&report);
            assert_eq!(
                events.iter().map(|event| event.index).collect::<Vec<_>>(),
                [0, 1]
            );
            let first = events[0].bracket;
            assert!(first.upper_dt_s - first.lower_dt_s <= limits.time_tolerance_s);
            assert_eq!(first.upper_dt_s, report.step.unwrap().dt_s);
            for event in events {
                assert_eq!(event.bracket.lower_dt_s, first.lower_dt_s);
                assert_eq!(event.bracket.upper_dt_s, first.upper_dt_s);
                assert!(event.bracket.lower_value < 0.0);
                assert!(event.bracket.upper_value >= 0.0);
            }
            // This later root is within the tolerance of the first root, but
            // has not yet been reached and is therefore not reported.
            assert!(report.values[2] < 0.0);
            assert_eq!(system.steps(), 1);
        }
    });
}

#[test]
fn entry_zeros_are_collected_only_after_validating_every_output() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut system = SystemBuilder::new().seal().unwrap();
        let directions = [
            EventDirection::Rising,
            EventDirection::Falling,
            EventDirection::Either,
        ];
        let report = system
            .step_to_events(
                1.0,
                config(),
                &directions,
                |_, values| {
                    values.copy_from_slice(&[-0.0, 1.0, 0.0]);
                },
                cx,
            )
            .unwrap();
        assert_eq!(report.status, EventSetStatus::AtStart(vec![0, 2]));
        assert_eq!(report.values[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(report.step, None);
        assert_eq!(report.trial_steps, 0);
        for invalid in [Some(f64::NAN), Some(f64::INFINITY), None] {
            let result = system.step_to_events(
                1.0,
                config(),
                &directions,
                |_, values| {
                    values[0] = 0.0;
                    values[1] = 1.0;
                    if let Some(value) = invalid {
                        values[2] = value;
                    }
                },
                cx,
            );
            assert!(matches!(
                result,
                Err(EventSetError::InvalidGuard { index: 2, .. })
            ));
        }
        let result = system.step_to_events(
            1.0,
            config(),
            &[],
            |_, _| panic!("empty guard set called"),
            cx,
        );
        assert!(matches!(result, Err(EventSetError::EmptySet)));
        assert_eq!(system.time_s(), 0.0);
        assert_eq!(system.steps(), 0);
        assert_eq!(system.ledger().steps, 0);
    });
}

#[test]
fn singleton_sets_preserve_scalar_status_values_steps_and_trial_counts() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for direction in [
            EventDirection::Either,
            EventDirection::Rising,
            EventDirection::Falling,
        ] {
            for root in [0.0, 0.37, 0.5, 1.0, 2.0] {
                let mut single = SystemBuilder::new().seal().unwrap();
                let mut set = SystemBuilder::new().seal().unwrap();
                // Negation deliberately produces -0.0 at an exact root.
                let value = |s: &System| -(s.time_s() - root);
                let scalar = single
                    .step_to_event(
                        1.0,
                        EventConfig {
                            direction,
                            time_tolerance_s: config().time_tolerance_s,
                            max_bisections: config().max_bisections,
                        },
                        value,
                        cx,
                    )
                    .unwrap();
                let batch = set
                    .step_to_events(1.0, config(), &[direction], |s, out| out[0] = value(s), cx)
                    .unwrap();
                match scalar.status {
                    EventStatus::AtStart => {
                        assert_eq!(batch.status, EventSetStatus::AtStart(vec![0]))
                    }
                    EventStatus::Advanced => assert_eq!(batch.status, EventSetStatus::Advanced),
                    EventStatus::Located(bracket) => {
                        let events = located(&batch);
                        assert_eq!(events, &[LocatedEvent { index: 0, bracket }]);
                        if bracket.upper_value == 0.0 {
                            assert_eq!(bracket.lower_value.to_bits(), 0.0_f64.to_bits());
                            assert_eq!(bracket.upper_value.to_bits(), (-0.0_f64).to_bits());
                        }
                    }
                }
                assert_eq!(scalar.value.to_bits(), batch.values[0].to_bits());
                assert_eq!(scalar.time_s.to_bits(), batch.time_s.to_bits());
                assert_eq!(scalar.step, batch.step);
                assert_eq!(scalar.trial_steps, batch.trial_steps);
                assert_eq!(single.steps(), set.steps());
                assert_eq!(single.ledger(), set.ledger());
            }
        }
    });
}

#[test]
fn batch_failures_restore_states_warm_starts_accounts_and_allowances() {
    for failure in 0..6 {
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let (mut system, device) = mechanical(true);
            system.step(0.01, cx).unwrap();
            assert!(
                device
                    .ports
                    .iter()
                    .any(|&p| system.last_allowance(p).iter().any(|&a| a > 0.0))
            );
            assert!(system.guess(device.ports[0], &[2.75]));
            let before = observe(&system, &device);
            let mut limits = config();
            if failure == 1 {
                limits.max_bisections = 1;
            }
            let mut calls = 0;
            let result = system.step_to_events(
                2.0,
                limits,
                &[
                    EventDirection::Falling,
                    EventDirection::Rising,
                    EventDirection::Either,
                ],
                |s, values| {
                    calls += 1;
                    values[0] = s.state::<Motion>(device.mass).unwrap().position_m - 0.23;
                    if failure != 5 || calls != 3 {
                        values[1] = s.time_s() - 1.65;
                    }
                    values[2] = if failure == 0 && calls == 3 {
                        f64::NAN
                    } else {
                        1.0
                    };
                    if (failure == 2 && calls == 3) || (failure == 3 && calls == 1) {
                        gate.request();
                    }
                    if failure == 3 {
                        values[0] = 0.0;
                    }
                    if failure == 4 {
                        values.fill(1.0);
                        if calls == 2 {
                            gate.request();
                        }
                    }
                },
                cx,
            );
            match failure {
                0 => assert!(matches!(
                    result,
                    Err(EventSetError::InvalidGuard { index: 2, .. })
                )),
                1 => assert!(matches!(
                    result,
                    Err(EventSetError::Bisections { count: 1, .. })
                )),
                5 => assert!(matches!(
                    result,
                    Err(EventSetError::InvalidGuard { index: 1, .. })
                )),
                _ => assert!(matches!(
                    result,
                    Err(EventSetError::Step(SystemError::Cancelled))
                )),
            }
            assert_eq!(observe(&system, &device), before);
        });
    }
}

#[test]
fn selected_coupled_endpoint_retains_its_own_allowances_and_ledger() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut system, device) = mechanical(true);
        system.step(0.01, cx).unwrap();
        let initial = system.checkpoint().unwrap();
        let report = system
            .step_to_events(
                2.0,
                config(),
                &[EventDirection::Falling, EventDirection::Rising],
                |s, values| {
                    values[0] = s.state::<Motion>(device.mass).unwrap().position_m - 0.23;
                    values[1] = s.time_s() - 1.65;
                },
                cx,
            )
            .unwrap();
        assert_eq!(located(&report)[0].index, 0);
        assert_eq!(system.steps(), initial.steps() + 1);
        assert_eq!(system.ledger().steps, initial.ledger().steps + 1);
        assert!(
            device
                .ports
                .iter()
                .any(|&p| system.last_allowance(p).iter().any(|&a| a > 0.0))
        );
        let accepted = observe(&system, &device);
        system.restore(&initial).unwrap();
        system.step(report.step.unwrap().dt_s, cx).unwrap();
        assert_eq!(observe(&system, &device), accepted);
    });
}

#[test]
fn many_guards_share_the_same_physics_trials_and_callback_evaluations() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut expected = None;
        for count in [1, 3, 32, 128] {
            let boundary_calls = Rc::new(Cell::new(0_u64));
            let observed_calls = Rc::clone(&boundary_calls);
            let mut builder = SystemBuilder::new();
            let mass = builder
                .add_with_tangent(
                    "mass",
                    Mass::new(1.0).unwrap(),
                    Motion {
                        position_m: 0.0,
                        velocity_m_s: 1.0,
                    },
                )
                .unwrap();
            let port = builder.port(mass, "p").unwrap();
            builder
                .boundary(
                    port,
                    Box::new(move |_, _, values| {
                        observed_calls.set(observed_calls.get() + 1);
                        values[0] = 0.0;
                    }),
                )
                .unwrap();
            let mut system = builder.seal().unwrap();
            let mut guard_calls = 0_u64;
            let report = system
                .step_to_events(
                    1.0,
                    config(),
                    &vec![EventDirection::Rising; count],
                    |s, values| {
                        guard_calls += 1;
                        values.fill(s.state::<Motion>(mass).unwrap().position_m - 0.37);
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(located(&report).len(), count);
            assert_eq!(guard_calls, report.trial_steps + 1);
            assert!(boundary_calls.get() >= report.trial_steps);
            let observed = (report.time_s, report.trial_steps, boundary_calls.get());
            if let Some(expected) = expected {
                assert_eq!(observed, expected);
            }
            expected = Some(observed);
            assert_eq!(system.steps(), 1);
            assert_eq!(system.ledger().steps, 1);
        }
    });
}
