//! G0/G3: dimensioned translation histories, canonical IR and frame semantics.

use fs_ga::{Point, Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::ir::{parse_ir, write_ir};
use fs_scenario::payload::{
    OrientationParity, Payload, PayloadId, PayloadMeta, QuantityContract, ReferenceSemantics,
    SampleSource, VectorPayload,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Environment, Frame, FrameId, FrameMotion, FrameMotionKind,
    FrameTree, Interp, Physics, Scenario, TimeSignal, ValidationBudget, WORLD,
};
use std::ops::ControlFlow;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);

fn translating() -> Frame {
    Frame {
        id: FrameId(2),
        name: "linear-drive".into(),
        parent: FrameId(1),
        motion: FrameMotion::Translating {
            direction: [1.0, 0.0, 0.0],
            displacement: TimeSignal::Table {
                times: vec![0.0, 1.0, 2.0],
                values: vec![-0.4, 0.2, 0.5],
                dims: LENGTH,
                interp: Interp::Linear,
            },
        },
    }
}

fn scenario() -> Scenario {
    let mut result = Scenario::new("translating-frame", 17, Environment::earth_lab());
    result.frames.add(Frame {
        id: FrameId(1),
        name: "drive-mount".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: Vec3::new(1.0, -2.0, 3.0),
        },
    });
    result.frames.add(translating());
    result
}

#[test]
fn dimensioned_translation_round_trips_and_composes_in_parent_axes() {
    let source = scenario();
    assert!(source.validate().is_empty());
    assert_eq!(
        source.frames.frames[1].motion.kind(),
        FrameMotionKind::Translating
    );
    let text = write_ir(&source);
    assert!(text.contains("(translating (vec 1 0 0)"));
    let decoded = parse_ir(&text).unwrap();
    assert_eq!(decoded.scenario(), &source);
    assert_eq!(write_ir(decoded.scenario()), text);
    for (time, displacement) in [(-1.0, -0.4), (0.5, -0.1), (1.5, 0.35), (3.0, 0.5)] {
        let pose = source.frames.world_pose(FrameId(2), time).unwrap();
        let point = pose.transform_point(Point::new(0.25, 0.0, 0.0)).unwrap();
        assert!((point.x - 1.0).abs() < 2e-14);
        assert!((point.y - (-2.0 + displacement + 0.25)).abs() < 2e-14);
        assert!((point.z - 3.0).abs() < 2e-14);
    }
}

#[test]
fn translation_does_not_add_carrier_speed_to_an_absolute_velocity_payload() {
    let source = scenario();
    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(VELOCITY),
        PayloadId::new("basis/frame-cartesian").unwrap(),
        FrameId(2),
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .unwrap();
    let value = VectorPayload::new(
        meta,
        SampleSource::Fixed([0.4, 0.0, 0.0].map(|v| QtyAny::new(v, VELOCITY)).to_vec()),
    )
    .unwrap();
    let boundary = BoundaryCondition {
        region: "inlet".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Typed(Payload::Vector(value))),
        compatibility: None,
        frame: 2,
    };
    // The mount rotates x to y. Its displacement rate is intentionally not
    // added: this existing boundary contract describes an absolute velocity.
    for time in [0.5, 1.5, 3.0] {
        let velocity = boundary
            .velocity_world_at(&source.frames, time, &mut || ControlFlow::Continue(()))
            .unwrap();
        assert!(velocity[0].abs() < 1e-14);
        assert!((velocity[1] - 0.4).abs() < 1e-14);
        assert!(velocity[2].abs() < 1e-14);
    }
}

#[test]
fn translation_rejects_bad_direction_units_and_nonfinite_displacement() {
    for (direction, distance) in [
        ([2.0, 0.0, 0.0], QtyAny::new(1.0, LENGTH)),
        ([f64::NAN, 0.0, 0.0], QtyAny::new(1.0, LENGTH)),
        ([1.0, 0.0, 0.0], QtyAny::new(1.0, VELOCITY)),
        ([1.0, 0.0, 0.0], QtyAny::new(f64::INFINITY, LENGTH)),
    ] {
        let mut source = scenario();
        source.frames.frames[1].motion = FrameMotion::Translating {
            direction,
            displacement: TimeSignal::Constant(distance),
        };
        assert!(!source.validate().is_empty());
        assert!(FrameTree::local_pose(&source.frames.frames[1], 0.5).is_err());
    }
}

#[test]
fn translation_history_scalars_are_included_in_admission_budget() {
    let source = scenario();
    let plan = source.validation_plan(ValidationBudget::default()).unwrap();
    assert_eq!(plan.signal_scalars, 6);
    assert!(
        source
            .validation_plan(ValidationBudget {
                max_signal_scalars: 5,
                ..ValidationBudget::default()
            })
            .is_err()
    );
}
