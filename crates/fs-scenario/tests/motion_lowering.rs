//! G0/G3 checks for the one-way `fs-scenario` -> `fs-motion` frame adapter.
//!
//! Pointwise sampling is a falsifier, not a certificate. The production
//! result remains the `fs-motion::CertifiedMotorTube`; these checks compare an
//! independent `FrameTree::world_pose` evaluation with its point-evaluation
//! view to catch composition, frame, phase, and double-cover mistakes.

use fs_ga::{Point, Quat, Vec3};
use fs_ivl::Interval;
use fs_motion::LowerToMotorTube;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{
    Frame, FrameId, FrameMotion, FrameMotionKind, FrameMotionLoweringError, FrameTree, TimeSignal,
    WORLD,
};

const RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn fixed(id: u32, parent: FrameId, translation: Vec3) -> Frame {
    fixed_pose(id, parent, Quat::identity(), translation)
}

fn fixed_pose(id: u32, parent: FrameId, orientation: Quat, translation: Vec3) -> Frame {
    Frame {
        id: FrameId(id),
        name: format!("fixed-{id}"),
        parent,
        motion: FrameMotion::Fixed {
            orientation,
            translation,
        },
    }
}

fn rotating(id: u32, parent: FrameId) -> Frame {
    Frame {
        id: FrameId(id),
        name: format!("rotating-{id}"),
        parent,
        motion: FrameMotion::Rotating {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::new(0.25, -0.4, 0.1),
            rate: QtyAny::new(2.75, RATE_DIMS),
        },
    }
}

fn point_distance(left: Point, right: Point) -> f64 {
    ((left.x - right.x).powi(2) + (left.y - right.y).powi(2) + (left.z - right.z).powi(2)).sqrt()
}

#[test]
fn rotating_frame_with_fixed_ancestors_agrees_with_lowered_motor_path() {
    let mut tree = FrameTree::new();
    tree.add(fixed_pose(
        1,
        WORLD,
        Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.35),
        Vec3::new(1.2, -0.6, 0.3),
    ));
    tree.add(fixed(2, FrameId(1), Vec3::new(-0.1, 0.35, 0.2)));
    tree.add(rotating(3, FrameId(2)));

    let builder = tree
        .rotating_motor_path(FrameId(3))
        .expect("fixed-ancestor rotating path admits");
    assert_eq!(builder.source_frame(), FrameId(3));

    let tube = builder
        .lower_to_motor_tube(Interval::new(-0.4, 1.3), 10, 6)
        .expect("certified screw tube builds");
    let path = tube.path();
    let probes = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.7, -0.2, 1.1),
        Point::new(-0.5, 0.9, -0.3),
    ];
    let mut max_deviation = 0.0f64;
    let mut max_enclosure_width = 0.0f64;
    for time in [-0.4, -0.125, 0.0, 0.375, 0.9, 1.3] {
        let expected = tree
            .world_pose(FrameId(3), time)
            .expect("scenario world pose");
        let (actual, sample) = path.motor_at(time).expect("tube point evaluation");
        max_enclosure_width = max_enclosure_width.max(sample.max_enclosure_width);
        for probe in probes {
            let expected_point = expected
                .transform_point(probe)
                .expect("finite scenario action");
            let actual_point = actual.transform_point(probe).expect("finite tube action");
            max_deviation = max_deviation.max(point_distance(expected_point, actual_point));
        }
    }

    println!(
        "{{\"suite\":\"fs-scenario/motion-lowering\",\"case\":\"rotating-fixed-chain\",\"max_deviation\":{max_deviation:.17e},\"max_enclosure_width\":{max_enclosure_width:.17e},\"tube_defect\":{:.17e}}}",
        tube.defect()
    );
    assert!(
        max_deviation <= 1.0e-9,
        "lowered point path disagrees with FrameTree world_pose: {max_deviation}"
    );
    assert!(
        tube.defect() <= 1.0e-8,
        "unit-axis lowering produced an excessive versor defect: {}",
        tube.defect()
    );
    assert!(
        max_enclosure_width <= 1.0e-8,
        "point-evaluation enclosure is unexpectedly loose: {max_enclosure_width}"
    );
}

#[test]
fn lowering_refuses_nonrotating_targets_and_dynamic_ancestors() {
    let mut fixed_target = FrameTree::new();
    fixed_target.add(fixed(1, WORLD, Vec3::new(0.0, 0.0, 0.0)));
    assert!(matches!(
        fixed_target.rotating_motor_path(FrameId(1)),
        Err(FrameMotionLoweringError::TargetNotRotating {
            frame: FrameId(1),
            actual: FrameMotionKind::Fixed,
        })
    ));

    let mut nested = FrameTree::new();
    nested.add(rotating(1, WORLD));
    nested.add(rotating(2, FrameId(1)));
    assert!(matches!(
        nested.rotating_motor_path(FrameId(2)),
        Err(FrameMotionLoweringError::DynamicAncestor {
            target: FrameId(2),
            ancestor: FrameId(1),
            actual: FrameMotionKind::Rotating,
        })
    ));

    let mut tilted = FrameTree::new();
    tilted.add(Frame {
        id: FrameId(1),
        name: "tilted-1".to_string(),
        parent: WORLD,
        motion: FrameMotion::Tilt {
            axis: [0.0, 1.0, 0.0],
            center: Vec3::new(0.0, 0.0, 0.0),
            angle: TimeSignal::Constant(QtyAny::dimensionless(0.2)),
        },
    });
    tilted.add(rotating(2, FrameId(1)));
    assert!(matches!(
        tilted.rotating_motor_path(FrameId(2)),
        Err(FrameMotionLoweringError::DynamicAncestor {
            target: FrameId(2),
            ancestor: FrameId(1),
            actual: FrameMotionKind::Tilt,
        })
    ));
}

#[test]
fn lowering_reconstruction_is_bit_stable() {
    let mut tree = FrameTree::new();
    tree.add(fixed(1, WORLD, Vec3::new(0.4, 0.2, -0.3)));
    tree.add(rotating(2, FrameId(1)));
    let builder = tree.rotating_motor_path(FrameId(2)).expect("path admits");
    let domain = Interval::new(-0.2, 0.8);
    let first = builder
        .lower_to_motor_tube(domain, 8, 4)
        .expect("first tube");
    let second = builder
        .lower_to_motor_tube(domain, 8, 4)
        .expect("second tube");

    assert_eq!(first.defect().to_bits(), second.defect().to_bits());
    for (left, right) in first.segments().iter().zip(second.segments()) {
        let left_components = left
            .components_over(left.domain())
            .expect("left components");
        let right_components = right
            .components_over(right.domain())
            .expect("right components");
        for (left, right) in left_components.iter().zip(right_components) {
            assert_eq!(left.lo().to_bits(), right.lo().to_bits());
            assert_eq!(left.hi().to_bits(), right.hi().to_bits());
        }
    }
}
