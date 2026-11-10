//! Step-doubling accuracy, retained physical accounts, and atomic refusals.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::{Class, Component, PortSpec, Proposal, Receipt};
use fs_system::elements::{Damper, Force, Mass, Motion, SourceAccount, Spring, Waveform};
use fs_system::{
    AdaptiveConfig, AdaptiveError, Checkpoint, ErrorNorm, ErrorTolerance, Ledger, NodeId, PortRef,
    SolveConfig, System, SystemBuilder, SystemError,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 43,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn controller() -> AdaptiveConfig {
    AdaptiveConfig {
        order: 2,
        min_step_s: 1.0e-12,
        max_attempts: 16,
    }
}

struct Device {
    mass: NodeId,
    spring: NodeId,
    source: NodeId,
    damper: Option<NodeId>,
    ports: Vec<PortRef>,
}

fn oscillator(driven: bool) -> (System, Device) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 1.0,
                velocity_m_s: 0.25,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(4.0, 0.0).unwrap(), -1.0)
        .unwrap();
    let wave = if driven {
        Waveform::sine(0.8, 0.7, 0.2, 0.3).unwrap()
    } else {
        Waveform::constant(0.0).unwrap()
    };
    let source = builder
        .add_with_tangent("source", Force::new(wave), SourceAccount::default())
        .unwrap();
    let mut ports = vec![
        builder.port(mass, "p").unwrap(),
        builder.port(spring, "a").unwrap(),
        builder.port(spring, "b").unwrap(),
        builder.port(source, "p").unwrap(),
    ];
    let mut members = vec![ports[1], ports[3]];
    let damper = if driven {
        let node = builder
            .add_with_tangent("damper", Damper::new(0.25).unwrap(), ())
            .unwrap();
        let a = builder.port(node, "a").unwrap();
        let b = builder.port(node, "b").unwrap();
        let heat = builder.port(node, "heat").unwrap();
        members.push(a);
        builder
            .boundary(b, Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
        builder
            .boundary(heat, Box::new(|_, _, values| values[0] = 300.0))
            .unwrap();
        ports.extend([a, b, heat]);
        Some(node)
    } else {
        None
    };
    builder.common_flow(ports[0], &members).unwrap();
    builder
        .boundary(ports[2], Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    (
        builder.seal().unwrap(),
        Device {
            mass,
            spring,
            source,
            damper,
            ports,
        },
    )
}

#[derive(Debug, PartialEq)]
struct Observed {
    time_bits: u64,
    steps: u64,
    motion_bits: [u64; 2],
    spring_bits: u64,
    source: SourceAccount,
    wires: Vec<(Vec<u64>, Vec<u64>)>,
    ledger: Ledger,
    config: SolveConfig,
}

fn observe(system: &System, device: &Device) -> Observed {
    let motion = system.state::<Motion>(device.mass).unwrap();
    Observed {
        time_bits: system.time_s().to_bits(),
        steps: system.steps(),
        motion_bits: [motion.position_m.to_bits(), motion.velocity_m_s.to_bits()],
        spring_bits: system.state::<f64>(device.spring).unwrap().to_bits(),
        source: *system.state::<SourceAccount>(device.source).unwrap(),
        wires: device
            .ports
            .iter()
            .map(|&port| {
                (
                    system.last_input(port).iter().map(|x| x.to_bits()).collect(),
                    system.last_output(port).iter().map(|x| x.to_bits()).collect(),
                )
            })
            .collect(),
        ledger: system.ledger().clone(),
        config: system.checkpoint().unwrap().config(),
    }
}

fn prime(system: &mut System, device: &Device) {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for dt in [0.003, 0.0017, 0.0023] {
            system.step(dt, cx).unwrap();
        }
    });
    system.set_config(SolveConfig {
        relative_tolerance: 2.0e-13,
        absolute_tolerance: 1.0e-14,
        max_iterations: 8,
        max_shortenings: 9,
        interface_energy_tolerance: 2.0e-10,
    });
    assert!(system.guess(device.ports[0], &[8.75]));
}

fn motion_norm(
    initial: &Checkpoint,
    coarse: &Checkpoint,
    refined: &Checkpoint,
    mass: NodeId,
    absolute: f64,
) -> f64 {
    let a = initial.state::<Motion>(mass).unwrap();
    let b = coarse.state::<Motion>(mass).unwrap();
    let c = refined.state::<Motion>(mass).unwrap();
    let mut norm = ErrorNorm::new();
    norm.include(
        a.position_m,
        b.position_m,
        c.position_m,
        ErrorTolerance {
            absolute,
            relative: 0.0,
        },
    );
    // The oscillator has angular frequency 2 rad/s: these two floors give
    // equal weight to displacement and velocity expressed as displacement.
    norm.include(
        a.velocity_m_s,
        b.velocity_m_s,
        c.velocity_m_s,
        ErrorTolerance {
            absolute: 2.0 * absolute,
            relative: 0.0,
        },
    );
    norm.value()
}

#[test]
fn tighter_local_tolerances_reduce_oscillator_error_against_the_analytic_trajectory() {
    let run = |absolute| {
        let (mut system, device) = oscillator(false);
        let mut cap = 0.5_f64;
        let mut accepted = 0_u64;
        let mut rejected = 0_u64;
        with_cx(&CancelGate::new_clock_free(), |cx| {
            while system.time_s() < 2.0 {
                assert!(accepted < 2_000, "adaptive trajectory must make progress");
                let requested = cap.min(0.5).min(2.0 - system.time_s());
                let report = system
                    .step_adaptive(
                        requested,
                        controller(),
                        |a, b, c| motion_norm(a, b, c, device.mass, absolute),
                        cx,
                    )
                    .unwrap();
                assert!(report.error_norm >= 0.0 && report.error_norm <= 1.0);
                assert!(report.dt_s > 0.0 && report.dt_s <= requested);
                assert!(report.suggested_dt_s.is_finite() && report.suggested_dt_s > 0.0);
                assert_eq!(report.time_s.to_bits(), system.time_s().to_bits());
                assert_eq!(report.trial_steps, 3 * u64::from(report.attempts));
                accepted += 1;
                rejected += u64::from(report.attempts - 1);
                assert_eq!(system.steps(), 2 * accepted);
                assert_eq!(system.ledger().steps, system.steps());
                cap = report.suggested_dt_s;
            }
        });
        assert_eq!(system.time_s(), 2.0);
        let (sin, cos) = (2.0 * system.time_s()).sin_cos();
        let expected_position = cos + 0.125 * sin;
        let expected_velocity = 0.25 * cos - 2.0 * sin;
        let motion = system.state::<Motion>(device.mass).unwrap();
        let error = (motion.position_m - expected_position)
            .hypot(0.5 * (motion.velocity_m_s - expected_velocity));
        let spring = *system.state::<f64>(device.spring).unwrap();
        let energy = 0.5 * motion.velocity_m_s.powi(2) + 2.0 * spring.powi(2);
        assert!((energy - 2.03125).abs() < 2.0e-10);
        assert!(system.ledger().energy_residual_j().abs() < 2.0e-10);
        (error, accepted, rejected)
    };
    let loose = run(1.0e-3);
    let tight = run(1.0e-6);
    assert!(loose.0 > 1.0e-4, "loose run must expose phase error: {loose:?}");
    assert!(tight.0 < 0.1 * loose.0, "loose {loose:?}, tight {tight:?}");
    assert!(tight.0 < 5.0e-4, "tight trajectory error: {tight:?}");
    assert!(tight.1 > 4 * loose.1, "tighter control must refine the trajectory");
    assert!(loose.2 > 0 && tight.2 > 0, "both initial coarse trials must be rejected");
}

#[test]
fn rejected_trials_leave_exactly_the_two_retained_half_steps_in_every_account() {
    let (mut system, device) = oscillator(true);
    prime(&mut system, &device);
    let before = observe(&system, &device);
    let initial = system.checkpoint().unwrap();
    let mut comparisons = 0;
    let report = with_cx(&CancelGate::new_clock_free(), |cx| {
        system
            .step_adaptive(
                0.2,
                controller(),
                |a, b, c| {
                    comparisons += 1;
                    assert_eq!(a.time_s().to_bits(), initial.time_s().to_bits());
                    assert_eq!(b.time_s().to_bits(), c.time_s().to_bits());
                    // A deterministic metric makes the large interval fail
                    // and accepts the smaller one with a nonzero estimate.
                    if b.time_s() - a.time_s() > 0.05 {
                        300.0
                    } else {
                        1.5
                    }
                },
                cx,
            )
            .unwrap()
    });
    assert_eq!(comparisons, 2);
    assert_eq!(report.attempts, 2);
    assert_eq!(report.trial_steps, 6);
    assert_eq!(report.error_norm, 0.5, "apply the order-two divisor exactly once");
    assert_eq!(report.shortenings, 0);
    let retained = observe(&system, &device);
    assert_eq!(retained.steps, before.steps + 2);
    assert_eq!(retained.ledger.steps, before.ledger.steps + 2);
    assert_eq!(retained.source.time_s().to_bits(), retained.time_bits);
    assert!((retained.source.energy_in_j - before.source.energy_in_j).abs() > 1.0e-8);
    let damper = device.damper.unwrap().index();
    assert!(
        retained.ledger.nodes[damper].dissipated_j.value()
            > before.ledger.nodes[damper].dissipated_j.value()
    );

    system.restore(&initial).unwrap();
    let iterations = with_cx(&CancelGate::new_clock_free(), |cx| {
        let first = system.step(0.5 * report.dt_s, cx).unwrap();
        let second = system.step(0.5 * report.dt_s, cx).unwrap();
        u64::from(first.iterations) + u64::from(second.iterations)
    });
    assert_eq!(retained, observe(&system, &device));
    assert!(report.iterations >= iterations && report.iterations > 0);
}

/// The first half changes the accepted state so only the second half asks
/// for a shorter interval. All state and receipts obey the ordinary contract.
struct LateShortening;

impl Component for LateShortening {
    type State = Vec<f64>;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        _inputs: &[f64],
        new: &mut Self::State,
        _outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        if old.len() == 1 && dt_s > 0.125 {
            return Ok(Proposal::Shorten { dt_s: 0.125 });
        }
        new.clone_from(old);
        new.push(dt_s);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn a_shortened_second_half_restarts_at_entry_and_only_equal_endpoints_are_compared() {
    let mut builder = SystemBuilder::new();
    let node = builder.add("late", LateShortening, Vec::new()).unwrap();
    let mut system = builder.seal().unwrap();
    let mut comparisons = 0;
    let report = with_cx(&CancelGate::new_clock_free(), |cx| {
        system
            .step_adaptive(
                1.0,
                controller(),
                |initial, coarse, refined| {
                    comparisons += 1;
                    assert_eq!(initial.steps(), 0);
                    assert!(initial.state::<Vec<f64>>(node).unwrap().is_empty());
                    assert_eq!(coarse.time_s(), 0.25);
                    assert_eq!(coarse.time_s().to_bits(), refined.time_s().to_bits());
                    assert_eq!(coarse.state::<Vec<f64>>(node).unwrap(), &[0.25]);
                    assert_eq!(refined.state::<Vec<f64>>(node).unwrap(), &[0.125, 0.125]);
                    0.0
                },
                cx,
            )
            .unwrap()
    });
    assert_eq!(comparisons, 1);
    assert_eq!(report.dt_s, 0.25);
    assert_eq!(report.time_s, 0.25);
    assert_eq!(report.attempts, 2);
    assert_eq!(report.trial_steps, 6);
    assert_eq!(report.shortenings, 1);
    assert_eq!(report.iterations, 0);
    assert_eq!(system.state::<Vec<f64>>(node).unwrap(), &[0.125, 0.125]);
    assert_eq!(system.steps(), 2);
    assert_eq!(system.ledger().steps, 2);
    assert_eq!(system.ledger().duration_s.value(), 0.25);
}

#[derive(Clone, Copy, Debug)]
enum Failure {
    InvalidNorm,
    Attempts,
    Cancel,
}

#[test]
fn invalid_estimates_exhausted_attempts_and_late_cancellation_restore_and_retry() {
    for failure in [Failure::InvalidNorm, Failure::Attempts, Failure::Cancel] {
        let (mut system, device) = oscillator(true);
        prime(&mut system, &device);
        let before = observe(&system, &device);
        let initial = system.checkpoint().unwrap();
        with_cx(&CancelGate::new_clock_free(), |cx| {
            system.step(0.03125, cx).unwrap();
            system.step(0.03125, cx).unwrap();
        });
        let expected_retry = observe(&system, &device);
        system.restore(&initial).unwrap();

        let mut comparisons = 0;
        let gate = CancelGate::new_clock_free();
        let error = with_cx(&gate, |cx| {
            system.step_adaptive(
                0.0625,
                AdaptiveConfig {
                    max_attempts: 2,
                    ..controller()
                },
                |a, b, c| {
                    comparisons += 1;
                    assert_eq!(a.time_s().to_bits(), before.time_bits);
                    assert_eq!(b.time_s().to_bits(), c.time_s().to_bits());
                    match failure {
                        Failure::InvalidNorm => ErrorNorm::new().value(),
                        Failure::Attempts => 300.0,
                        Failure::Cancel => {
                            // Both fine steps have already completed. The
                            // final cancellation poll must rewind them too.
                            cx.cancel_gate().request();
                            0.0
                        }
                    }
                },
                cx,
            )
        })
        .expect_err("injected adaptive refusal");
        match (failure, error) {
            (Failure::InvalidNorm, AdaptiveError::InvalidEstimate { value }) => {
                assert!(value.is_nan());
                assert_eq!(comparisons, 1);
            }
            (Failure::Attempts, AdaptiveError::Attempts { attempts, dt_s }) => {
                assert_eq!(attempts, 2);
                assert_eq!(comparisons, 2);
                assert!(dt_s > 0.0 && dt_s < 0.0625);
            }
            (Failure::Cancel, AdaptiveError::Step(SystemError::Cancelled)) => {
                assert_eq!(comparisons, 1);
            }
            (_, error) => panic!("{failure:?} returned the wrong error: {error:?}"),
        }
        assert_eq!(observe(&system, &device), before, "rollback after {failure:?}");

        let retry = with_cx(&CancelGate::new_clock_free(), |cx| {
            system.step_adaptive(0.0625, controller(), |_, _, _| 0.0, cx)
        })
        .unwrap();
        assert_eq!(retry.attempts, 1);
        assert_eq!(retry.trial_steps, 3);
        assert_eq!(observe(&system, &device), expected_retry, "retry after {failure:?}");
    }
}

#[test]
fn the_error_norm_uses_initial_and_refined_scales_and_the_worst_coordinate() {
    for (initial, coarse, refined) in [(-4.0, 2.0, 1.0), (0.25, 5.0, 4.0)] {
        let mut norm = ErrorNorm::new();
        assert!(norm.value().is_nan());
        norm.include(
            initial,
            coarse,
            refined,
            ErrorTolerance {
                absolute: 0.5,
                relative: 0.125,
            },
        );
        assert_eq!(norm.value(), 1.0);
        norm.include(
            0.0,
            3.0,
            1.0,
            ErrorTolerance {
                absolute: 0.25,
                relative: 0.0,
            },
        );
        assert_eq!(norm.value(), 8.0);
    }
}

#[test]
fn invalid_norm_coordinates_and_tolerances_remain_poisoned_after_valid_data() {
    let valid = ErrorTolerance {
        absolute: 1.0,
        relative: 0.0,
    };
    let invalid_tolerances = [
        ErrorTolerance { absolute: 0.0, ..valid },
        ErrorTolerance { absolute: -1.0, ..valid },
        ErrorTolerance { absolute: f64::NAN, ..valid },
        ErrorTolerance { absolute: f64::INFINITY, ..valid },
        ErrorTolerance { relative: -0.5, ..valid },
        ErrorTolerance { relative: f64::NAN, ..valid },
        ErrorTolerance { relative: f64::INFINITY, ..valid },
    ];
    for tolerance in invalid_tolerances {
        let mut norm = ErrorNorm::new();
        norm.include(0.0, 2.0, 1.0, valid);
        norm.include(0.0, 2.0, 1.0, tolerance);
        norm.include(0.0, 3.0, 1.0, valid);
        assert!(norm.value().is_nan(), "{tolerance:?}");
    }
    for (initial, coarse, refined) in [
        (f64::NAN, 1.0, 0.0),
        (0.0, f64::INFINITY, 1.0),
        (0.0, 1.0, f64::NEG_INFINITY),
        (0.0, f64::MAX, -f64::MAX),
    ] {
        let mut norm = ErrorNorm::new();
        norm.include(initial, coarse, refined, valid);
        norm.include(0.0, 1.0, 0.0, valid);
        assert!(norm.value().is_nan());
    }
}
