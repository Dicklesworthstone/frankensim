//! Sample clocks and discrete recurrence on a physical kernel device.

#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};
use fs_system::elements::{
    ControlledSource, ControllerState, ElementError, SampledController, SourceAccount, Storage,
};
use fs_system::{NodeId, System, SystemBuilder};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 23,
                kernel_id: 41,
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

fn scalar(period: f64) -> SampledController {
    SampledController::new(
        period,
        1,
        1,
        1,
        vec![0.5],
        vec![0.25],
        vec![2.0],
        vec![0.125],
    )
    .unwrap()
}

fn words(controller: &SampledController, state: &ControllerState) -> Vec<u64> {
    let mut words = Vec::new();
    assert!(controller.encode_state(state, &mut words));
    words
}

#[test]
fn rectangular_law_uses_preupdate_memory_and_conforms_at_sample_and_hold() {
    let controller = SampledController::new(
        0.5,
        2,
        2,
        3,
        vec![1.0, 0.25, -0.5, 0.75],
        vec![0.5, 0.0, 0.0, 0.25],
        vec![1.0, 0.0, 0.0, 1.0, 1.0, -1.0],
        vec![0.25, 0.5, -0.5, 0.25, 1.0, -1.0],
    )
    .unwrap();
    assert_eq!(controller.class(), Class::Smooth);
    assert_eq!(
        (
            controller.state_dimension(),
            controller.input_dimension(),
            controller.output_dimension()
        ),
        (2, 2, 3)
    );
    assert_eq!(controller.sample_period_s(), 0.5);
    assert_eq!(controller.ports()[0].id().as_str(), "input0");
    assert_eq!(controller.ports()[1].id().as_str(), "input1");
    assert_eq!(controller.ports()[2].id().as_str(), "output");
    let initial = ControllerState::new(&controller, vec![1.0, -2.0]).unwrap();
    let stale = ControllerState::new(&controller, vec![20.0, 30.0]).unwrap();
    let mut accepted = initial.clone();
    let mut out = [0.0; 3];
    let before = format!("{controller:?}");
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let case = Case {
            state: &initial,
            stale: &stale,
            dt_s: 0.125,
            inputs: &[2.0, -1.0],
        };
        check_component(&controller, &case, cx).unwrap();
        check_port_tangent(&controller, &case, 1.0e-7, 1.0e-8, cx).unwrap();
        controller
            .propose(&initial, 0.125, case.inputs, &mut accepted, &mut out, cx)
            .unwrap();
        assert_eq!(accepted.coordinates(), &[1.5, -2.25]);
        assert_eq!(out, [1.0, -3.25, 6.0]);
        assert_eq!(accepted.samples(), 1);
        let case = Case {
            state: &accepted,
            stale: &stale,
            dt_s: 0.125,
            inputs: &[200.0, 400.0],
        };
        check_component(&controller, &case, cx).unwrap();
        check_port_tangent(&controller, &case, 1.0e-7, 1.0e-8, cx).unwrap();
        let mut held = stale.clone();
        controller
            .propose(&accepted, 0.125, case.inputs, &mut held, &mut out, cx)
            .unwrap();
        assert_eq!(held.coordinates(), accepted.coordinates());
        assert_eq!(held.held_output(), accepted.held_output());
        assert_eq!(held.samples(), 1);
        assert_eq!(held.time_s(), 0.25);
    });
    assert_eq!(
        format!("{controller:?}"),
        before,
        "scratch is not observable state"
    );
}

fn prescribed_input(time_s: f64) -> f64 {
    if time_s < 1.0 { 2.0 } else { -1.0 }
}

struct Guard(Rc<Cell<bool>>);

impl Component for Guard {
    type State = ();
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &()) -> Result<f64, Self::Error> {
        Ok(0.0)
    }
    fn propose(
        &self,
        _: &(),
        _: f64,
        _: &[f64],
        _: &mut (),
        _: &mut [f64],
        _: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
    fn admit(&self, _: &(), _: &(), _: &Receipt) -> Result<(), Self::Error> {
        if self.0.get() {
            Err("deliberate final refusal")
        } else {
            Ok(())
        }
    }
}

struct Device {
    controller: NodeId,
    capacitor: NodeId,
    source: NodeId,
    refuse: Rc<Cell<bool>>,
}

fn device(period: f64) -> (System, Device) {
    let controller = scalar(period);
    let state = ControllerState::new(&controller, vec![-1.0]).unwrap();
    let mut builder = SystemBuilder::new();
    let control = builder
        .add_with_tangent("controller", controller, state)
        .unwrap();
    let kind = PortKind::ElectricalVoltageCurrent;
    let capacitor = builder
        .add_with_tangent("capacitor", Storage::in_domain(kind, 0.5).unwrap(), 0.75)
        .unwrap();
    let source = builder
        .add_with_tangent(
            "supply",
            ControlledSource::flow(kind),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .bond(
            builder.port(capacitor, "p").unwrap(),
            builder.port(source, "p").unwrap(),
        )
        .unwrap();
    builder
        .signal(
            builder.port(control, "output").unwrap(),
            builder.port(source, "command").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(control, "input0").unwrap(),
            Box::new(|t, _, out| out[0] = prescribed_input(t)),
        )
        .unwrap();
    let refuse = Rc::new(Cell::new(false));
    builder.add("guard", Guard(Rc::clone(&refuse)), ()).unwrap();
    (
        builder.seal().unwrap(),
        Device {
            controller: control,
            capacitor,
            source,
            refuse,
        },
    )
}

#[test]
fn physics_subcycles_without_resampling_and_the_held_source_closes_energy() {
    let (mut system, device) = device(0.25);
    let (mut memory, mut command, mut charge, mut source_work) = (-1.0, 0.0, 0.75, 0.0);
    let mut next_sample = 0_u64;
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        while system.time_s() < 2.0 {
            let t0 = system.time_s();
            if t0 == next_sample as f64 * 0.25 {
                let input = prescribed_input(t0);
                command = 2.0 * memory + 0.125 * input;
                memory = 0.5 * memory + 0.25 * input;
                next_sample += 1;
            }
            let report = system.step(0.07_f64.min(2.0 - t0), cx).unwrap();
            let next_charge = charge - report.dt_s * command;
            source_work += report.dt_s * command * (charge + next_charge);
            charge = next_charge;
            let state = system.state::<ControllerState>(device.controller).unwrap();
            assert_eq!(state.samples(), next_sample);
            assert_eq!(state.coordinates(), &[memory]);
            assert_eq!(state.held_output(), &[command]);
            assert_eq!(state.time_s().to_bits(), system.time_s().to_bits());
            near(
                *system.state::<f64>(device.capacitor).unwrap(),
                charge,
                1.0e-12,
            );
            near(
                system
                    .state::<SourceAccount>(device.source)
                    .unwrap()
                    .energy_in_j,
                source_work,
                1.0e-12,
            );
            near(system.ledger().energy_residual_j(), 0.0, 1.0e-11);
        }
    });
    assert_eq!(next_sample, 8);
    assert!(system.steps() > next_sample * 2);
    near(charge * charge + source_work, 0.75 * 0.75, 1.0e-11);
}

#[test]
fn nonbinary_periods_reach_indexed_boundaries_without_cumulative_sample_drift() {
    let (mut system, device) = device(0.1);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for index in 1..=40_u64 {
            let target = index as f64 * 0.1;
            while system.time_s() < target {
                let report = system.step(0.03, cx).unwrap();
                assert!(report.time_s <= target, "sample boundary was crossed");
            }
            assert_eq!(system.time_s().to_bits(), target.to_bits());
            assert_eq!(
                system
                    .state::<ControllerState>(device.controller)
                    .unwrap()
                    .samples(),
                index
            );
        }
    });
}

fn observation(
    system: &System,
    device: &Device,
) -> (ControllerState, u64, SourceAccount, fs_system::Ledger) {
    (
        system
            .state::<ControllerState>(device.controller)
            .unwrap()
            .clone(),
        system.state::<f64>(device.capacitor).unwrap().to_bits(),
        *system.state::<SourceAccount>(device.source).unwrap(),
        system.ledger().clone(),
    )
}

#[test]
fn final_refusal_and_checkpoint_replay_restore_sample_count_hold_and_physics_together() {
    let (mut system, device) = device(0.25);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        system.step(0.125, cx).unwrap();
        let checkpoint = system.checkpoint().unwrap();
        let start = observation(&system, &device);
        let mut expected = Vec::new();
        for _ in 0..5 {
            system.step(0.125, cx).unwrap();
            expected.push(observation(&system, &device));
        }
        system.restore(&checkpoint).unwrap();
        assert_eq!(observation(&system, &device), start);
        for wanted in expected {
            let before = observation(&system, &device);
            device.refuse.set(true);
            assert!(system.step(0.125, cx).is_err());
            assert_eq!(observation(&system, &device), before);
            device.refuse.set(false);
            system.step(0.125, cx).unwrap();
            assert_eq!(observation(&system, &device), wanted);
        }
    });
}

#[test]
fn sampled_gain_needs_no_memory_and_proposals_shorten_before_the_next_sample() {
    let controller =
        SampledController::new(0.25, 0, 1, 1, vec![], vec![], vec![], vec![3.0]).unwrap();
    let old = ControllerState::new(&controller, vec![]).unwrap();
    let mut candidate = old.clone();
    let mut output = [99.0];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        assert_eq!(controller.limit_dt_s(&old, 2.0).unwrap(), 0.25);
        assert_eq!(controller.next_event_s(&old), Some(0.25));
        assert_eq!(
            controller
                .propose(&old, 2.0, &[2.0], &mut candidate, &mut output, cx)
                .unwrap(),
            Proposal::Shorten { dt_s: 0.25 }
        );
        assert_eq!(candidate, old);
        assert_eq!(output, [99.0]);
        controller
            .propose(&old, 0.125, &[2.0], &mut candidate, &mut output, cx)
            .unwrap();
        assert_eq!(output, [6.0]);
        assert!(candidate.coordinates().is_empty());
        let sampled = candidate.clone();
        controller
            .propose(&sampled, 0.125, &[-3.0], &mut candidate, &mut output, cx)
            .unwrap();
        assert_eq!(output, [6.0]);
        assert_eq!(candidate.samples(), 1);
        let boundary = candidate.clone();
        controller
            .propose(&boundary, 0.125, &[-3.0], &mut candidate, &mut output, cx)
            .unwrap();
        assert_eq!(output, [-9.0]);
        assert_eq!(candidate.samples(), 2);
    });
}

#[test]
fn construction_rejects_nonfinite_or_wrong_dimensions_before_any_sample() {
    for period in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            SampledController::new(period, 1, 1, 1, vec![1.0], vec![1.0], vec![1.0], vec![0.0])
                .is_err()
        );
    }
    for (n, m, p) in [(1, 0, 1), (1, 1, 0), (usize::MAX, 1, 1), (0, usize::MAX, 2)] {
        assert!(SampledController::new(0.1, n, m, p, vec![], vec![], vec![], vec![]).is_err());
    }
    for which in 0..4 {
        let mut matrices = [vec![1.0], vec![1.0], vec![1.0], vec![1.0]];
        matrices[which] = vec![f64::NAN];
        let [a, b, c, d] = matrices;
        assert!(SampledController::new(0.1, 1, 1, 1, a, b, c, d).is_err());
        let mut matrices = [vec![1.0], vec![1.0], vec![1.0], vec![1.0]];
        matrices[which] = vec![];
        let [a, b, c, d] = matrices;
        assert!(SampledController::new(0.1, 1, 1, 1, a, b, c, d).is_err());
    }
    let controller = scalar(0.25);
    assert!(ControllerState::new(&controller, vec![]).is_err());
    assert!(ControllerState::new(&controller, vec![f64::INFINITY]).is_err());
}

#[test]
fn overflow_invalid_calls_and_cancellation_publish_neither_memory_nor_output() {
    let controller = SampledController::new(
        0.25,
        2,
        2,
        2,
        vec![1.0, 0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0, 1.0],
        vec![1.0, 0.0, 1.0, f64::MAX],
    )
    .unwrap();
    let old = ControllerState::new(&controller, vec![1.0, 2.0]).unwrap();
    let stale = ControllerState::new(&controller, vec![30.0, 40.0]).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for (dt, input) in [
            (0.125, vec![1.0, 2.0]),
            (0.0, vec![1.0, 0.0]),
            (f64::NAN, vec![1.0, 0.0]),
            (0.125, vec![1.0]),
            (0.125, vec![f64::NAN, 0.0]),
        ] {
            let mut candidate = stale.clone();
            let mut out = [33.0, 44.0];
            assert!(
                controller
                    .propose(&old, dt, &input, &mut candidate, &mut out, cx)
                    .is_err()
            );
            assert_eq!(words(&controller, &candidate), words(&controller, &stale));
            assert_eq!(out, [33.0, 44.0]);
        }
        let mut out = [33.0, 44.0];
        assert!(
            controller
                .jvp(&old, 0.125, &[1.0, 0.0], &[1.0, 2.0], &mut out, cx)
                .is_err()
        );
        assert_eq!(out, [33.0, 44.0]);
        assert!(
            controller
                .vjp(&old, 0.125, &[1.0, 0.0], &[1.0, 2.0], &mut out, cx)
                .is_err()
        );
        assert_eq!(out, [33.0, 44.0]);
    });
    gate.request();
    with_cx(&gate, |cx| {
        let mut candidate = stale.clone();
        let mut out = [33.0, 44.0];
        assert_eq!(
            controller.propose(&old, 0.125, &[1.0, 0.0], &mut candidate, &mut out, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(candidate, stale);
        assert_eq!(out, [33.0, 44.0]);
        assert_eq!(
            controller.jvp(&old, 0.125, &[1.0, 0.0], &[1.0, 2.0], &mut out, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(
            controller.vjp(&old, 0.125, &[1.0, 0.0], &[1.0, 2.0], &mut out, cx),
            Err(ElementError::Cancelled)
        );
        assert_eq!(out, [33.0, 44.0]);
    });
}

#[test]
fn native_codec_retains_compensated_clock_memory_and_hold_and_rejects_bad_state() {
    let controller = scalar(0.1);
    let mut accepted = ControllerState::new(&controller, vec![-1.0]).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..3 {
            let old = accepted.clone();
            controller
                .propose(&old, 0.03, &[2.0], &mut accepted, &mut [0.0], cx)
                .unwrap();
        }
    });
    let encoded = words(&controller, &accepted);
    let decoded = controller.decode_state(&encoded).unwrap();
    assert_eq!(words(&controller, &decoded), encoded);
    assert_eq!(
        controller.next_event_s(&decoded),
        controller.next_event_s(&accepted)
    );
    assert!(
        controller
            .decode_state(&encoded[..encoded.len() - 1])
            .is_none()
    );
    let mut extra = encoded.clone();
    extra.push(0);
    assert!(controller.decode_state(&extra).is_none());
    for (index, value) in [
        (0, u64::MAX),
        (0, 0),
        (1, f64::NAN.to_bits()),
        (3, f64::INFINITY.to_bits()),
        (4, f64::NAN.to_bits()),
    ] {
        let mut broken = encoded.clone();
        broken[index] = value;
        assert!(controller.decode_state(&broken).is_none());
    }
    let mut identity = Vec::new();
    controller.identity_words(&mut identity);
    let mut changed = Vec::new();
    scalar(0.2).identity_words(&mut changed);
    assert_ne!(identity, changed);
    let changed_gain =
        SampledController::new(0.1, 1, 1, 1, vec![0.5], vec![0.25], vec![2.0], vec![0.25]).unwrap();
    changed.clear();
    changed_gain.identity_words(&mut changed);
    assert_ne!(identity, changed);
}
