//! G0/G3 integration of authored physical frequencies into scenario histories.

use fs_qty::{Dims, QtyAny};
use fs_scenario::frame::{FrameKinematicsBudget, FrameKinematicsError};
use fs_scenario::ir::{parse_ir, write_ir};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Compat, DEFAULT_MAXIMUM_CYCLES_PER_STEP, Environment,
    Frame, FrameId, FrameMotion, FrameTree, LinearChirp, Physics, Scenario, ScenarioError,
    TimeSignal, WORLD,
};
use std::ops::ControlFlow;

const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MASS_FLOW: Dims = Dims([0, 1, -1, 0, 0, 0]);

fn chirped_frame(id: u32, parent: FrameId, hz: f64, amplitude: f64, tilt: bool) -> Frame {
    let dims = if tilt {
        Dims::NONE
    } else {
        Dims([1, 0, 0, 0, 0, 0])
    };
    let signal = TimeSignal::Chirp(
        LinearChirp::new(
            0.0,
            1.0,
            hz,
            hz,
            0.0,
            QtyAny::new(0.0, dims),
            QtyAny::new(amplitude, dims),
        )
        .unwrap(),
    );
    Frame {
        id: FrameId(id),
        name: format!("frequency-frame-{id}"),
        parent,
        motion: if tilt {
            FrameMotion::Tilt {
                axis: [0.0, 0.0, 1.0],
                center: Default::default(),
                angle: signal,
            }
        } else {
            FrameMotion::Translating {
                direction: [1.0, 0.0, 0.0],
                displacement: signal,
            }
        },
    }
}

#[test]
fn prepared_source_resolution_visits_only_the_selected_ancestor_chain() {
    let mut tree = FrameTree::new();
    tree.add(chirped_frame(1, WORLD, 1.0, 0.001, false));
    tree.add(chirped_frame(2, WORLD, 1000.0, 0.001, false));
    tree.add(chirped_frame(3, FrameId(2), 1.0, 0.001, false));
    let mut checkpoint = || ControlFlow::Continue(());
    let slow = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut checkpoint,
        )
        .unwrap();
    assert!(
        slow.admits_sampled_step(
            0.25,
            0.001,
            DEFAULT_MAXIMUM_CYCLES_PER_STEP,
            &mut checkpoint
        )
        .unwrap()
    );
    let fast_parent = tree
        .prepare_kinematics(
            FrameId(3),
            FrameKinematicsBudget::default(),
            &mut checkpoint,
        )
        .unwrap();
    assert!(
        !fast_parent
            .admits_sampled_step(
                0.25,
                0.001,
                DEFAULT_MAXIMUM_CYCLES_PER_STEP,
                &mut checkpoint
            )
            .unwrap()
    );
    assert!(
        fast_parent
            .admits_sampled_step(
                0.25,
                0.001 / 32.0,
                DEFAULT_MAXIMUM_CYCLES_PER_STEP,
                &mut checkpoint
            )
            .unwrap()
    );
}

#[test]
fn prepared_tilt_resolution_bounds_angular_travel_and_preserves_cancellation() {
    let mut tree = FrameTree::new();
    tree.add(chirped_frame(1, WORLD, 1.0, -8.0, true));
    tree.add(chirped_frame(2, WORLD, 1.0, 0.1, true));
    let mut checkpoint = || ControlFlow::Continue(());
    let large = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut checkpoint,
        )
        .unwrap();
    let small = tree
        .prepare_kinematics(
            FrameId(2),
            FrameKinematicsBudget::default(),
            &mut checkpoint,
        )
        .unwrap();
    assert!(
        small
            .admits_sampled_step(0.25, 0.01, DEFAULT_MAXIMUM_CYCLES_PER_STEP, &mut checkpoint)
            .unwrap()
    );
    assert!(
        !large
            .admits_sampled_step(0.25, 0.01, DEFAULT_MAXIMUM_CYCLES_PER_STEP, &mut checkpoint)
            .unwrap()
    );
    assert!(
        large
            .admits_sampled_step(
                0.25,
                0.001,
                DEFAULT_MAXIMUM_CYCLES_PER_STEP,
                &mut checkpoint
            )
            .unwrap()
    );
    let mut polls = 0;
    let cancelled =
        small.admits_sampled_step(0.25, 0.001, DEFAULT_MAXIMUM_CYCLES_PER_STEP, &mut || {
            polls += 1;
            if polls == 2 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        });
    assert!(matches!(cancelled, Err(FrameKinematicsError::Cancelled)));
    let world = tree
        .prepare_kinematics(WORLD, FrameKinematicsBudget::default(), &mut checkpoint)
        .unwrap();
    assert!(
        world
            .admits_sampled_step(0.25, 0.01, 0.0, &mut checkpoint)
            .is_err()
    );
}

fn source(dims: Dims) -> LinearChirp {
    LinearChirp::new(
        0.0,
        1.0,
        64.0,
        128.0,
        0.0,
        QtyAny::new(3.0, dims),
        QtyAny::new(2.0, dims),
    )
    .unwrap()
}

fn pressure_scenario() -> Scenario {
    let mut scenario = Scenario::new("frequency-sweep", 17, Environment::earth_lab());
    scenario.base_bcs.push(BoundaryCondition {
        region: "loaded-face".into(),
        physics: Physics::Elasticity,
        kind: BcKind::Traction,
        value: Some(BcValue::Signal(TimeSignal::Chirp(source(PRESSURE)))),
        compatibility: None,
        frame: 0,
    });
    scenario
}

#[test]
fn native_frequency_source_round_trips_through_canonical_scenario_ir() {
    let scenario = pressure_scenario();
    assert!(scenario.validate().is_empty(), "{:?}", scenario.validate());
    let wire = write_ir(&scenario);
    assert!(wire.contains("(chirp 0 1 64 128 0 "));
    let decoded = parse_ir(&wire).unwrap();
    assert_eq!(write_ir(decoded.scenario()), wire);
    assert!(decoded.scenario().validate().is_empty());
    let Some(BcValue::Signal(signal)) = &decoded.scenario().base_bcs[0].value else {
        panic!("frequency source lost on decode");
    };
    assert_eq!(signal.dims(), PRESSURE);
    for time in [-0.25, 0.0, 0.003, 0.1, 0.5, 1.0, 1.25, 1000.125] {
        assert_eq!(
            signal.eval(time).unwrap(),
            source(PRESSURE).eval(time).unwrap()
        );
    }
}

#[test]
fn invalid_wire_frequency_is_refused_at_the_original_signal_form() {
    let wire = write_ir(&pressure_scenario());
    let invalid = wire.replacen("(chirp 0 1 64 128 0 ", "(chirp 0 1 -64 128 0 ", 1);
    assert_ne!(wire, invalid);
    let error = parse_ir(&invalid).unwrap_err();
    assert!(matches!(error, ScenarioError::Parse { .. }), "{error}");
    assert!(error.to_string().contains("frequency"), "{error}");
}

#[test]
fn continuous_frequency_evaluation_preserves_cancellation_and_dimensions() {
    let signal = TimeSignal::Chirp(source(PRESSURE));
    for time in [-0.125, 0.0, 0.25, 1.0, 1.25] {
        let sample = signal
            .eval_with_checkpoint(time, &mut |_| Ok::<(), ()>(()))
            .unwrap()
            .unwrap();
        assert_eq!(sample, signal.eval(time).unwrap());
        assert_eq!(sample.dims, PRESSURE);
    }
    let refused = signal.eval_with_checkpoint(0.25, &mut |phase| {
        if phase == "signal evaluation publication" {
            Err("cancelled")
        } else {
            Ok(())
        }
    });
    assert!(matches!(refused, Err("cancelled")));
}

#[test]
fn oscillatory_flux_cannot_pass_an_all_time_check_by_aliasing_its_knots() {
    // cos(2*pi*32*t)-1 vanishes at both declared knots and at the scenario's
    // t=0 checkpoint, while it is -2 halfway through its first cycle.
    let chirp = LinearChirp::new(
        0.0,
        1.0,
        32.0,
        32.0,
        0.0,
        QtyAny::new(-1.0, MASS_FLOW),
        QtyAny::new(1.0, MASS_FLOW),
    )
    .unwrap();
    assert_eq!(chirp.eval(0.0).unwrap().value, 0.0);
    assert_eq!(chirp.eval(1.0).unwrap().value, 0.0);
    assert!((chirp.eval(1.0 / 64.0).unwrap().value + 2.0).abs() < 1e-14);
    let mut scenario = Scenario::new("oscillatory-flux", 7, Environment::earth_lab());
    scenario.base_bcs.push(BoundaryCondition {
        region: "inlet".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Signal(TimeSignal::Chirp(chirp))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    });
    let findings = scenario.validate();
    assert!(
        findings
            .iter()
            .any(|finding| finding.code == "flux-certification-unavailable"),
        "aliased checkpoints must not certify the source: {findings:?}"
    );
}
