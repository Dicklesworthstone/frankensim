//! In-memory checkpoints restore the accepted physics, accounting and solver
//! state needed to replay a driven device, including after a failed trial.

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::{Class, Component, PortSpec, Proposal, Receipt};
use fs_system::elements::{
    Damper, Force, HeatCapacity, Mass, Motion, SourceAccount, Spring, Waveform,
};
use fs_system::{
    CheckpointError, Ledger, NodeId, Phase, PortRef, RouteTemperature, SolveConfig, StepReport,
    System, SystemBuilder, SystemError, heat_sink,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 13,
            kernel_id: 37,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Interruption {
    None,
    Refuse,
    Cancel,
}

/// A final component whose owned vector records accepted step sizes. The
/// external switch injects failures; it is deliberately not checkpointed.
struct TrialControl(Rc<Cell<Interruption>>);

#[derive(Debug, PartialEq, Eq)]
struct Refused;

impl Component for TrialControl {
    type State = Vec<u64>;
    type Error = Refused;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &Self::State) -> Result<f64, Refused> {
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
    ) -> Result<Proposal, Refused> {
        new.clone_from(old);
        new.push(dt_s.to_bits());
        if self.0.get() == Interruption::Cancel {
            cx.cancel_gate().request();
        }
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }

    fn admit(
        &self,
        _old: &Self::State,
        _new: &Self::State,
        _receipt: &Receipt,
    ) -> Result<(), Refused> {
        if self.0.get() == Interruption::Refuse {
            Err(Refused)
        } else {
            Ok(())
        }
    }
}

struct Device {
    mass: NodeId,
    spring: NodeId,
    damper: NodeId,
    source: NodeId,
    body: NodeId,
    control: NodeId,
    ports: Vec<PortRef>,
    interruption: Rc<Cell<Interruption>>,
}

fn device() -> (System, Device) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(0.8).unwrap(),
            Motion {
                position_m: 0.07,
                velocity_m_s: 0.34,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(60.0, 0.0).unwrap(), -0.07)
        .unwrap();
    let damper = builder
        .add_with_tangent("damper", Damper::new(1.1).unwrap(), ())
        .unwrap();
    let source = builder
        .add_with_tangent(
            "source",
            Force::new(Waveform::sine(1.2, 1.7, 0.3, 0.2).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(4.0).unwrap(), 300.0)
        .unwrap();
    let interruption = Rc::new(Cell::new(Interruption::None));
    let control = builder
        .add(
            "control",
            TrialControl(Rc::clone(&interruption)),
            Vec::new(),
        )
        .unwrap();
    let ports = [
        (mass, "p"),
        (spring, "a"),
        (spring, "b"),
        (damper, "a"),
        (damper, "b"),
        (damper, "heat"),
        (source, "p"),
        (body, "heat"),
    ]
    .map(|(node, name)| builder.port(node, name).unwrap());
    builder
        .common_flow(ports[0], &[ports[1], ports[3], ports[6]])
        .unwrap();
    // A moving support gives the system a nonconstant external input as
    // well as the source's internally clocked sinusoid. This is the exact
    // mean of a velocity that is linear in time.
    builder
        .boundary(
            ports[2],
            Box::new(|t0, dt, value| value[0] = 0.05 + 0.02 * (t0 + 0.5 * dt)),
        )
        .unwrap();
    builder
        .boundary(ports[4], Box::new(|_, _, value| value[0] = 0.0))
        .unwrap();
    builder
        .heat_route(
            ports[5],
            &[(ports[7], 0.8)],
            Some(heat_sink("ambient", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    (
        builder.seal().unwrap(),
        Device {
            mass,
            spring,
            damper,
            source,
            body,
            control,
            ports: ports.to_vec(),
            interruption,
        },
    )
}

/// Observe through the ordinary System API, independently of checkpoint
/// accessors. Float bits include signed zeros; ledger equality compares the
/// private compensation of every Sum as well as its displayed value.
#[derive(Debug, PartialEq)]
struct Observed {
    time_bits: u64,
    steps: u64,
    motion_bits: [u64; 2],
    spring_bits: u64,
    temperature_bits: u64,
    source: SourceAccount,
    source_energy_bits: u64,
    accepted_sizes: Vec<u64>,
    wires: Vec<(Vec<u64>, Vec<u64>)>,
    ledger: Ledger,
}

fn observe(system: &System, device: &Device) -> Observed {
    let motion = system.state::<Motion>(device.mass).unwrap();
    let source = *system.state::<SourceAccount>(device.source).unwrap();
    Observed {
        time_bits: system.time_s().to_bits(),
        steps: system.steps(),
        motion_bits: [motion.position_m.to_bits(), motion.velocity_m_s.to_bits()],
        spring_bits: system.state::<f64>(device.spring).unwrap().to_bits(),
        temperature_bits: system.state::<f64>(device.body).unwrap().to_bits(),
        source,
        source_energy_bits: source.energy_in_j.to_bits(),
        accepted_sizes: system.state::<Vec<u64>>(device.control).unwrap().clone(),
        wires: device
            .ports
            .iter()
            .map(|&port| {
                (
                    system
                        .last_input(port)
                        .iter()
                        .map(|v| v.to_bits())
                        .collect(),
                    system
                        .last_output(port)
                        .iter()
                        .map(|v| v.to_bits())
                        .collect(),
                )
            })
            .collect(),
        ledger: system.ledger().clone(),
    }
}

fn dt(index: usize) -> f64 {
    [0.003, 0.0017, 0.0023, 0.0009][index % 4]
}

fn run(
    system: &mut System,
    device: &Device,
    start: usize,
    count: usize,
) -> Vec<(StepReport, Observed)> {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        (start..start + count)
            .map(|index| {
                let report = system.step(dt(index), cx).expect("accepted step");
                (report, observe(system, device))
            })
            .collect()
    })
}

#[test]
fn a_driven_dissipative_device_replays_every_step_from_one_reusable_checkpoint() {
    let (mut system, device) = device();
    run(&mut system, &device, 0, 73);
    let before = observe(&system, &device);
    let checkpoint = system.checkpoint().unwrap();
    assert_eq!(checkpoint.time_s().to_bits(), before.time_bits);
    assert_eq!(checkpoint.steps(), before.steps);
    assert_eq!(checkpoint.ledger(), &before.ledger);
    assert_eq!(
        checkpoint
            .state::<Motion>(device.mass)
            .unwrap()
            .position_m
            .to_bits(),
        before.motion_bits[0]
    );
    assert_eq!(
        checkpoint.state::<SourceAccount>(device.source),
        Some(&before.source)
    );
    assert!(checkpoint.state::<Vec<u64>>(device.mass).is_none());
    // The compared accounts carry actual source work, moving-boundary work,
    // dissipation, and delivered heat rather than just empty ledger rows.
    assert!(before.source.energy_in_j.abs() > 1.0e-8);
    assert!(before.ledger.boundaries[0].energy_in_j.value().abs() > 1.0e-8);
    assert!(
        before.ledger.nodes[device.damper.index()]
            .dissipated_j
            .value()
            > 0.0
    );
    assert!(before.ledger.sinks[0].heat_j.value() > 0.0);

    let expected = run(&mut system, &device, 73, 257);
    assert_eq!(
        checkpoint.state::<Vec<u64>>(device.control),
        Some(&before.accepted_sizes),
        "advancing must not mutate the checkpoint's owned vector"
    );
    for _ in 0..2 {
        system.restore(&checkpoint).unwrap();
        assert_eq!(observe(&system, &device), before);
        assert_eq!(run(&mut system, &device, 73, 257), expected);
    }
    assert_eq!(
        system
            .state::<SourceAccount>(device.source)
            .unwrap()
            .time_s()
            .to_bits(),
        system.time_s().to_bits()
    );
}

fn replay_after_interruption(interruption: Interruption) {
    let (mut system, device) = device();
    run(&mut system, &device, 0, 29);
    let before = observe(&system, &device);
    let checkpoint = system.checkpoint().unwrap();
    let expected = run(&mut system, &device, 29, 41);
    let after = observe(&system, &device);

    device.interruption.set(interruption);
    let gate = CancelGate::new_clock_free();
    let error = with_cx(&gate, |cx| system.step(0.0041, cx)).expect_err("injected failure");
    match interruption {
        Interruption::Refuse => {
            let SystemError::Component { node, phase, error } = error else {
                panic!("expected the final component's admission refusal");
            };
            assert_eq!(node, "control");
            assert_eq!(phase, Phase::Admit);
            assert_eq!(error.downcast_ref::<Refused>(), Some(&Refused));
        }
        Interruption::Cancel => assert!(matches!(error, SystemError::Cancelled)),
        Interruption::None => unreachable!(),
    }
    assert_eq!(observe(&system, &device), after);
    // Restore is independent of the cancelled execution context and may
    // rewind past successful steps even while candidate scratch is dirty.
    system.restore(&checkpoint).unwrap();
    assert_eq!(observe(&system, &device), before);
    device.interruption.set(Interruption::None);
    assert_eq!(run(&mut system, &device, 29, 41), expected);
}

#[test]
fn restore_after_a_late_refusal_replays_the_accepted_trajectory() {
    replay_after_interruption(Interruption::Refuse);
}

#[test]
fn restore_after_mid_trial_cancellation_replays_with_a_fresh_context() {
    replay_after_interruption(Interruption::Cancel);
}

#[test]
fn a_foreign_checkpoint_is_refused_before_mutation_even_for_identical_graphs() {
    let (mut origin, origin_device) = device();
    run(&mut origin, &origin_device, 0, 13);
    let checkpoint = origin.checkpoint().unwrap();
    let (mut other, other_device) = device();
    run(&mut other, &other_device, 0, 37);
    let before = observe(&other, &other_device);
    let config = other.checkpoint().unwrap().config();
    assert!(matches!(
        other.restore(&checkpoint),
        Err(CheckpointError::ForeignSystem)
    ));
    assert_eq!(observe(&other, &other_device), before);
    assert_eq!(other.checkpoint().unwrap().config(), config);

    // A foreign first state's different Rust type must not be downcast or
    // copied before the identity check, either.
    let mut builder = SystemBuilder::new();
    let body = builder
        .add("body", HeatCapacity::new(2.0).unwrap(), 290.0)
        .unwrap();
    let port = builder.port(body, "heat").unwrap();
    builder
        .boundary(port, Box::new(|_, _, v| v[0] = 0.0))
        .unwrap();
    let mut different = builder.seal().unwrap();
    let ledger = different.ledger().clone();
    assert!(matches!(
        different.restore(&checkpoint),
        Err(CheckpointError::ForeignSystem)
    ));
    assert_eq!(
        different.state::<f64>(body).unwrap().to_bits(),
        290.0_f64.to_bits()
    );
    assert_eq!(different.steps(), 0);
    assert_eq!(different.time_s().to_bits(), 0.0_f64.to_bits());
    assert_eq!(different.ledger(), &ledger);
}

#[test]
fn restoring_recovers_configuration_and_an_explicit_interface_guess() {
    let (mut system, device) = device();
    let config = SolveConfig {
        relative_tolerance: 2.0e-12,
        absolute_tolerance: 1.0e-14,
        max_iterations: 12,
        max_shortenings: 17,
        interface_energy_tolerance: 2.0e-9,
    };
    system.set_config(config);
    run(&mut system, &device, 0, 7);
    assert!(system.guess(device.ports[0], &[17.3]));
    let before = observe(&system, &device);
    let checkpoint = system.checkpoint().unwrap();
    assert_eq!(checkpoint.config(), config);

    run(&mut system, &device, 7, 11);
    assert!(system.guess(device.ports[0], &[-8.5]));
    let changed = SolveConfig {
        relative_tolerance: 0.0,
        absolute_tolerance: 0.0,
        max_iterations: 0,
        max_shortenings: 0,
        interface_energy_tolerance: 0.0,
    };
    system.set_config(changed);
    assert_eq!(system.checkpoint().unwrap().config(), changed);
    system.restore(&checkpoint).unwrap();
    assert_eq!(system.checkpoint().unwrap().config(), config);
    assert_eq!(observe(&system, &device), before);
    // Restore leaves the checkpoint reusable, including its configuration.
    let expected = run(&mut system, &device, 7, 11);
    system.restore(&checkpoint).unwrap();
    assert_eq!(run(&mut system, &device, 7, 11), expected);
}
