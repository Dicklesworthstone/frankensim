//! State-triggered thresholds: independent event-time oracles and trial rollback.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::{Class, Component, PortSpec, Proposal, Receipt};
use fs_system::elements::{HeatCapacity, Mass, Motion, Spring};
use fs_system::{
    EventConfig, EventDirection, EventError, EventStatus, Ledger, NodeId, PortRef, SolveConfig,
    System, SystemBuilder, SystemError,
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

fn config(direction: EventDirection) -> EventConfig {
    EventConfig {
        direction,
        time_tolerance_s: 1.0e-10,
        max_bisections: 64,
    }
}

struct Device {
    mass: NodeId,
    spring: Option<NodeId>,
    ports: Vec<PortRef>,
}

fn falling() -> (System, Device) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 1.37,
                velocity_m_s: -0.24,
            },
        )
        .unwrap();
    let port = builder.port(mass, "p").unwrap();
    builder
        .boundary(port, Box::new(|_, _, value| value[0] = -2.0 * 9.81))
        .unwrap();
    (
        builder.seal().unwrap(),
        Device {
            mass,
            spring: None,
            ports: vec![port],
        },
    )
}

fn oscillator() -> (System, Device) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 1.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(4.0, 0.0).unwrap(), -1.0)
        .unwrap();
    let ports = [
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "a").unwrap(),
        builder.port(spring, "b").unwrap(),
    ];
    builder.bond(ports[0], ports[1]).unwrap();
    builder
        .boundary(ports[2], Box::new(|_, _, value| value[0] = 0.0))
        .unwrap();
    (
        builder.seal().unwrap(),
        Device {
            mass,
            spring: Some(spring),
            ports: ports.to_vec(),
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
fn falling_mass_locates_the_quadratic_root_and_retains_one_audited_step() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut system, device) = falling();
        let initial = system.checkpoint().unwrap();
        let expected_s = (-0.24 + (0.24_f64.powi(2) + 2.0 * 9.81 * 1.37).sqrt()) / 9.81;
        let report = system
            .step_to_event(
                1.0,
                config(EventDirection::Falling),
                |s| s.state::<Motion>(device.mass).unwrap().position_m,
                cx,
            )
            .unwrap();
        let EventStatus::Located(bracket) = report.status else {
            panic!("{report:?}")
        };
        assert!(bracket.lower_dt_s <= expected_s && expected_s <= bracket.upper_dt_s);
        assert!(bracket.upper_dt_s - bracket.lower_dt_s <= 1.0e-10);
        assert!(bracket.lower_value >= 0.0 && bracket.upper_value <= 0.0);
        assert!((report.time_s - expected_s).abs() <= 1.0e-10);
        assert!(report.trial_steps > 2);
        assert_eq!(system.steps(), 1);
        assert_eq!(system.ledger().steps, 1);
        assert_eq!(report.step.unwrap().dt_s, bracket.upper_dt_s);
        let accepted = observe(&system, &device);

        system.restore(&initial).unwrap();
        system.step(bracket.upper_dt_s, cx).unwrap();
        assert_eq!(observe(&system, &device), accepted);
        system.restore(&initial).unwrap();
        let replay = system
            .step_to_event(
                1.0,
                config(EventDirection::Falling),
                |s| s.state::<Motion>(device.mass).unwrap().position_m,
                cx,
            )
            .unwrap();
        assert_eq!(replay, report);
        assert_eq!(observe(&system, &device), accepted);
    });
}

#[test]
fn heating_threshold_respects_an_earlier_scheduled_boundary() {
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(10.0).unwrap(), 300.0)
        .unwrap();
    let heat = builder.port(body, "heat").unwrap();
    builder
        .scheduled_boundary(heat, Box::new(|_, _, value| value[0] = 7.0), vec![0.2])
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let first = system
            .step_to_event(
                2.0,
                config(EventDirection::Rising),
                |s| *s.state::<f64>(body).unwrap() - 301.0,
                cx,
            )
            .unwrap();
        assert_eq!(first.status, EventStatus::Advanced);
        assert_eq!(first.time_s, 0.2);
        assert_eq!(first.trial_steps, 1);
        let second = system
            .step_to_event(
                2.0,
                config(EventDirection::Rising),
                |s| *s.state::<f64>(body).unwrap() - 301.0,
                cx,
            )
            .unwrap();
        assert!(matches!(second.status, EventStatus::Located(_)));
        // C * (T_target - T_initial) / prescribed power.
        assert!((second.time_s - 10.0 / 7.0).abs() < 1.0e-9);
        assert_eq!(system.steps(), 2);
        assert!((system.ledger().boundary_energy_in_j() - 10.0).abs() < 1.0e-8);
        assert!(system.ledger().energy_residual_j().abs() < 1.0e-10);
    });
}

#[test]
fn no_crossing_and_wrong_direction_are_the_ordinary_step_to_the_bit() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (cap, direction) in [(0.1, EventDirection::Either), (1.0, EventDirection::Rising)] {
            let (mut system, device) = falling();
            let initial = system.checkpoint().unwrap();
            let report = system
                .step_to_event(
                    cap,
                    config(direction),
                    |s| s.state::<Motion>(device.mass).unwrap().position_m,
                    cx,
                )
                .unwrap();
            assert_eq!(report.status, EventStatus::Advanced);
            assert_eq!(report.trial_steps, 1);
            let expected = observe(&system, &device);
            system.restore(&initial).unwrap();
            assert_eq!(system.step(cap, cx).unwrap(), report.step.unwrap());
            assert_eq!(observe(&system, &device), expected);
        }
    });
}

#[test]
fn an_initial_zero_is_reported_without_advancing_for_each_direction() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut system, device) = falling();
        let initial = observe(&system, &device);
        for direction in [
            EventDirection::Either,
            EventDirection::Rising,
            EventDirection::Falling,
        ] {
            let report = system
                .step_to_event(1.0, config(direction), |_| -0.0, cx)
                .unwrap();
            assert_eq!(report.status, EventStatus::AtStart);
            assert_eq!(report.step, None);
            assert_eq!(report.trial_steps, 0);
            assert_eq!(observe(&system, &device), initial);
        }
    });
}

#[test]
fn invalid_guards_exhausted_bisections_and_cancellation_restore_complete_entry_state() {
    for failure in 0..5 {
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let (mut system, device) = oscillator();
            system.step(0.01, cx).unwrap();
            // A loop gives a nonzero last-step allowance; public restore alone
            // would clear it, so it must be preserved by the event transaction.
            assert!(
                device
                    .ports
                    .iter()
                    .any(|&p| system.last_allowance(p).iter().any(|&a| a > 0.0))
            );
            assert!(system.guess(device.ports[0], &[2.75]));
            let before = observe(&system, &device);
            let mut limits = config(EventDirection::Falling);
            if failure == 1 {
                limits.max_bisections = 1;
            }
            let mut calls = 0;
            let result = system.step_to_event(
                2.0,
                limits,
                |s| {
                    calls += 1;
                    if failure == 0 && calls == 3 {
                        return f64::NAN;
                    }
                    if (failure == 2 && calls == 3) || (failure == 3 && calls == 1) {
                        gate.request();
                    }
                    if failure == 4 {
                        if calls == 2 {
                            gate.request();
                        }
                        return 1.0; // Cancellation on the no-crossing return path.
                    }
                    if failure == 3 {
                        return 0.0;
                    } // Cancellation on AtStart.
                    s.state::<Motion>(device.mass).unwrap().position_m
                },
                cx,
            );
            match failure {
                0 => assert!(matches!(result, Err(EventError::InvalidGuard { .. }))),
                1 => assert!(matches!(result, Err(EventError::Bisections { .. }))),
                _ => assert!(matches!(
                    result,
                    Err(EventError::Step(SystemError::Cancelled))
                )),
            }
            assert_eq!(observe(&system, &device), before);
        });
    }
}

#[test]
fn a_located_loop_retains_its_own_interface_allowance_and_accounts() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut system, device) = oscillator();
        system.step(0.01, cx).unwrap();
        let initial = system.checkpoint().unwrap();
        let report = system
            .step_to_event(
                2.0,
                config(EventDirection::Falling),
                |s| s.state::<Motion>(device.mass).unwrap().position_m - 0.23,
                cx,
            )
            .unwrap();
        assert!(matches!(report.status, EventStatus::Located(_)));
        assert!(
            device
                .ports
                .iter()
                .any(|&p| system.last_allowance(p).iter().any(|&a| a > 0.0))
        );
        let accepted = observe(&system, &device);
        assert_eq!(system.steps(), initial.steps() + 1);
        system.restore(&initial).unwrap();
        system.step(report.step.unwrap().dt_s, cx).unwrap();
        assert_eq!(observe(&system, &device), accepted);
    });
}

#[derive(Debug)]
struct Refused;

/// A small owner that limits just one trial interval; all other steps are
/// normal. This exercises negotiation without a second integration scheme.
struct TrialLimit {
    duration_s: f64,
    shorten_to_s: Option<f64>,
}

impl Component for TrialLimit {
    type State = ();
    type Error = Refused;
    fn ports(&self) -> &[PortSpec] {
        &[]
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &()) -> Result<f64, Refused> {
        Ok(0.0)
    }
    fn propose(
        &self,
        _: &(),
        dt_s: f64,
        _: &[f64],
        _: &mut (),
        _: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<Proposal, Refused> {
        if dt_s == self.duration_s {
            match self.shorten_to_s {
                Some(dt_s) => Ok(Proposal::Shorten { dt_s }),
                None => Err(Refused),
            }
        } else {
            Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
        }
    }
}

fn limited(limit: TrialLimit) -> System {
    let mut builder = SystemBuilder::new();
    builder.add("trial_limit", limit, ()).unwrap();
    builder.seal().unwrap()
}

#[test]
fn shortened_trials_use_actual_duration_and_refuse_a_stalled_bracket() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut system = limited(TrialLimit {
            duration_s: 0.5,
            shorten_to_s: Some(0.25),
        });
        let report = system
            .step_to_event(
                1.0,
                config(EventDirection::Rising),
                |s| s.time_s() - 0.73,
                cx,
            )
            .unwrap();
        assert!(matches!(report.status, EventStatus::Located(_)));
        assert!((report.time_s - 0.73).abs() <= 1.0e-10);
        assert_eq!(system.steps(), 1);

        for (duration_s, shorten_to_s) in [(0.75, Some(0.25)), (0.5, None)] {
            let mut system = limited(TrialLimit {
                duration_s,
                shorten_to_s,
            });
            let initial = system.checkpoint().unwrap();
            let result = system.step_to_event(
                1.0,
                config(EventDirection::Rising),
                |s| s.time_s() - 0.73,
                cx,
            );
            if shorten_to_s.is_some() {
                assert!(matches!(result, Err(EventError::Resolution { .. })));
            } else {
                assert!(matches!(result, Err(EventError::Step(_))));
            }
            assert_eq!(system.time_s(), initial.time_s());
            assert_eq!(system.steps(), initial.steps());
            assert_eq!(system.ledger(), initial.ledger());
        }
    });
}

#[test]
fn exact_endpoint_zero_and_clock_resolution_are_explicit() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut system = SystemBuilder::new().seal().unwrap();
        let report = system
            .step_to_event(
                0.5,
                config(EventDirection::Rising),
                |s| s.time_s() - 0.5,
                cx,
            )
            .unwrap();
        let EventStatus::Located(bracket) = report.status else {
            panic!("{report:?}")
        };
        assert_eq!(bracket.lower_dt_s, 0.5);
        assert_eq!(bracket.upper_dt_s, 0.5);
        assert_eq!(report.trial_steps, 1);

        system.step(2.0_f64.powi(54), cx).unwrap();
        let initial = system.checkpoint().unwrap();
        let result = system.step_to_event(0.25, config(EventDirection::Either), |_| 1.0, cx);
        assert!(matches!(result, Err(EventError::Resolution { .. })));
        assert_eq!(system.time_s(), initial.time_s());
        assert_eq!(system.steps(), initial.steps());
        assert_eq!(system.ledger(), initial.ledger());
    });
}

#[test]
fn invalid_configuration_refuses_before_calling_the_guard() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (cap, tolerance, bisections) in [
            (1.0, 0.0, 64),
            (1.0, f64::NAN, 64),
            (1.0, f64::INFINITY, 64),
            (1.0, 1.0e-9, 0),
            (0.0, 1.0e-9, 64),
            (f64::NAN, 1.0e-9, 64),
        ] {
            let (mut system, device) = falling();
            let initial = observe(&system, &device);
            let limits = EventConfig {
                direction: EventDirection::Either,
                time_tolerance_s: tolerance,
                max_bisections: bisections,
            };
            assert!(
                system
                    .step_to_event(cap, limits, |_| panic!("invalid call ran guard"), cx)
                    .is_err()
            );
            assert_eq!(observe(&system, &device), initial);
        }
    });
}
