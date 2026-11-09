//! Closed-form falsifiers for rigorously enclosed accelerating trajectories.
//! Taylor remainders and the native complete time cover provide the proof;
//! the point samples below are independent falsifiers, never certificates.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Motor;
use fs_geom::Point3;
use fs_geom::fixtures::SphereChart;
use fs_ivl::Interval;
use fs_motion::events::{
    CrossingDirection, EventScanConfig, ScanVerdict, plane_crossing_guard, scan_events,
};
use fs_motion::{
    AcceleratedRigidParams, EnclosureClass, MotionError, SpacetimeChart, SweepDecision, SweptChart,
    SweptConfig, accelerated_rigid_tube,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xACC_E1,
                kernel_id: 0xAC,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn translation(reference_time_s: f64) -> AcceleratedRigidParams {
    AcceleratedRigidParams {
        axis: [0.0, 0.0, 1.0],
        center: [0.0; 3],
        angular_velocity: 0.0,
        angular_acceleration: 0.0,
        linear_velocity: [-4.0, 0.0, 0.0],
        linear_acceleration: [4.0, 0.0, 0.0],
        reference_time_s,
        base_pose: Motor::translator(1.0, 0.0, 0.0),
    }
}

fn encloses_oracle(actual: [Interval; 3], expected: [f64; 3]) {
    for k in 0..3 {
        // Constant pose products are rounded; the independent real-number
        // rotation formula is a cross-implementation falsifier with tolerance.
        assert!(
            actual[k].lo() <= expected[k] + 1e-9 && actual[k].hi() >= expected[k] - 1e-9,
            "axis {k}: {:?} missed {}",
            actual[k],
            expected[k]
        );
        assert!(
            actual[k].width() < 1e-5,
            "vacuous pointwise band: {:?}",
            actual[k]
        );
    }
}

#[test]
fn translation_keeps_reference_time_and_turning_point_acceleration() {
    with_cx(|cx, _| {
        for reference in [0.0, 16.0, -16.0] {
            let tube = accelerated_rigid_tube(
                &translation(reference),
                Interval::new(reference, reference + 2.0),
                8,
                4,
            )
            .unwrap();
            for tau in [0.0, 0.375, 1.0, 2.0] {
                let point = tube
                    .point_motion_over(
                        Point3::new(0.0, 2.0, 3.0),
                        Interval::point(reference + tau),
                        cx,
                    )
                    .unwrap();
                assert_eq!(point.class, EnclosureClass::Certified);
                assert!(point.position[0].contains(1.0 - 4.0 * tau + 2.0 * tau * tau));
                assert!(point.velocity[0].contains(-4.0 + 4.0 * tau));
                assert!(point.acceleration[0].contains(4.0));
                assert!(point.acceleration[0].width() < 1e-9);
                assert!(point.position[1].contains(2.0));
                assert!(point.position[2].contains(3.0));
            }
            let span = tube
                .point_motion_over(
                    Point3::new(0.0, 0.0, 0.0),
                    Interval::new(reference + 0.5, reference + 1.5),
                    cx,
                )
                .unwrap();
            assert!(span.velocity[0].contains(-2.0) && span.velocity[0].contains(2.0));
            assert!(span.position[0].contains(-1.0));
            assert!(span.acceleration[0].contains(4.0));
        }
    });
}

#[test]
fn nonaxial_translation_combines_tangential_centripetal_and_world_frame_motion() {
    with_cx(|cx, _| {
        let mut p = translation(8.0);
        p.center = [0.25, -0.5, 0.75];
        p.angular_velocity = 0.75;
        p.angular_acceleration = 1.25;
        p.linear_velocity = [0.25, -0.5, 0.125];
        p.linear_acceleration = [1.0, 0.5, -0.25];
        p.base_pose = Motor::translator(0.5, 1.0, 2.0)
            .compose(&Motor::rotor([1.0, 0.0, 0.0], std::f64::consts::FRAC_PI_2));
        let tube = accelerated_rigid_tube(&p, Interval::new(8.0, 9.0), 12, 8).unwrap();
        for tau in [0.0, 0.375, 1.0] {
            let angle: f64 = 0.75 * tau + 0.625 * tau * tau;
            let (s, c) = angle.sin_cos();
            let omega = 0.75 + 1.25 * tau;
            let pos: [f64; 3] = std::array::from_fn(|k| {
                p.center[k]
                    + [c, s, 0.0][k]
                    + p.linear_velocity[k] * tau
                    + 0.5 * p.linear_acceleration[k] * tau * tau
            });
            let vel: [f64; 3] = std::array::from_fn(|k| {
                p.linear_velocity[k]
                    + p.linear_acceleration[k] * tau
                    + [-omega * s, omega * c, 0.0][k]
            });
            let acc = [
                1.0 - 1.25 * s - omega * omega * c,
                0.5 + 1.25 * c - omega * omega * s,
                -0.25,
            ];
            let result = tube
                .point_motion_over(
                    Point3::new(1.25, -0.5, 0.75),
                    Interval::point(8.0 + tau),
                    cx,
                )
                .unwrap();
            encloses_oracle(result.position, [pos[0] + 0.5, -pos[2] + 1.0, pos[1] + 2.0]);
            encloses_oracle(result.velocity, [vel[0], -vel[2], vel[1]]);
            encloses_oracle(result.acceleration, [acc[0], -acc[2], acc[1]]);
        }
    });
}

#[test]
fn cover_sign_and_uniform_motor_scale_leave_physical_derivatives_invariant() {
    with_cx(|cx, _| {
        let mut p = translation(0.0);
        p.angular_velocity = -1.0;
        p.angular_acceleration = 1.0;
        p.linear_velocity = [0.0; 3];
        p.linear_acceleration = [0.0; 3];
        p.base_pose = Motor::identity();
        let positive = accelerated_rigid_tube(&p, Interval::new(0.0, 4.0), 12, 16).unwrap();
        p.base_pose = Motor(p.base_pose.0.scale(-2.0));
        let negative = accelerated_rigid_tube(&p, Interval::new(0.0, 4.0), 12, 16).unwrap();
        // Includes zero angular velocity at t=1 and the scalar sign change
        // when the accumulated angle crosses pi near the end of the path.
        for t in [0.0, 1.0, 2.0, 3.75, 4.0] {
            let angle: f64 = -t + 0.5 * t * t;
            let (s, c) = angle.sin_cos();
            let w = t - 1.0;
            for tube in [&positive, &negative] {
                let result = tube
                    .point_motion_over(Point3::new(1.0, 0.0, 0.0), Interval::point(t), cx)
                    .unwrap();
                encloses_oracle(result.position, [c, s, 0.0]);
                encloses_oracle(result.velocity, [-w * s, w * c, 0.0]);
                encloses_oracle(result.acceleration, [-s - w * w * c, c - w * w * s, 0.0]);
            }
        }
        assert!(negative.motion().defect() >= 3.0);
    });
}

#[test]
fn event_scan_finds_both_interior_crossings_when_endpoints_share_safe_side() {
    with_cx(|cx, _| {
        let tube =
            accelerated_rigid_tube(&translation(16.0), Interval::new(16.0, 18.0), 8, 4).unwrap();
        let point = Point3::new(0.0, 0.0, 0.0);
        for t in [16.0, 18.0] {
            assert!(
                tube.point_motion_over(point, Interval::point(t), cx)
                    .unwrap()
                    .position[0]
                    .lo()
                    > 0.9
            );
        }
        let guard =
            plane_crossing_guard(tube.motion(), tube.rate(), point, [1.0, 0.0, 0.0], 0.0).unwrap();
        let scan = scan_events(
            tube.motion(),
            &guard,
            tube.motion().domain(),
            &EventScanConfig::default(),
            cx,
        )
        .unwrap();
        assert_eq!(scan.count.verdict, ScanVerdict::Complete);
        assert_eq!(scan.count.confirmed, 2);
        assert_eq!(scan.count.possible_windows, 0);
        assert!(
            scan.certified[0]
                .window
                .contains(17.0 - std::f64::consts::FRAC_1_SQRT_2)
        );
        assert!(
            scan.certified[1]
                .window
                .contains(17.0 + std::f64::consts::FRAC_1_SQRT_2)
        );
        assert_eq!(scan.certified[0].direction, CrossingDirection::Falling);
        assert_eq!(scan.certified[1].direction, CrossingDirection::Rising);
    });
}

#[test]
fn turning_point_contact_remains_an_explicit_unresolved_event() {
    with_cx(|cx, _| {
        let mut p = translation(16.0);
        p.linear_velocity[0] = -2.0;
        p.linear_acceleration[0] = 2.0; // x=(t-17)^2: touching, no transverse crossing.
        let tube = accelerated_rigid_tube(&p, Interval::new(16.0, 18.0), 8, 4).unwrap();
        let guard = plane_crossing_guard(
            tube.motion(),
            tube.rate(),
            Point3::new(0.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            0.0,
        )
        .unwrap();
        let scan = scan_events(
            tube.motion(),
            &guard,
            tube.motion().domain(),
            &EventScanConfig::default(),
            cx,
        )
        .unwrap();
        assert_eq!(scan.count.confirmed, 0);
        assert_ne!(scan.count.verdict, ScanVerdict::Complete);
        assert!(
            scan.possible
                .iter()
                .any(|event| event.window.contains(17.0))
        );
    });
}

#[test]
fn native_sweep_contains_interior_excursion_missed_by_endpoint_poses() {
    with_cx(|cx, _| {
        let tube =
            accelerated_rigid_tube(&translation(16.0), Interval::new(16.0, 18.0), 8, 4).unwrap();
        let sweep = SweptChart::new(
            SpacetimeChart::new(
                SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 0.125,
                },
                tube.into_motion(),
            ),
            SweptConfig {
                value_tolerance: 0.01,
                time_tolerance: 1e-6,
                max_subdivisions: 256,
            },
            cx,
        )
        .unwrap();
        let receipt = sweep.evaluate(Point3::new(-1.0, 0.0, 0.0), cx).unwrap();
        assert_eq!(receipt.decision, SweepDecision::Enclosure);
        assert!(receipt.implicit.contains(-0.125));
        assert!(receipt.implicit.hi() < 0.0);
        assert_eq!(receipt.evaluated_cells, 1 + 2 * receipt.subdivisions);
    });
}

#[test]
fn invalid_sources_and_queries_refuse_without_changing_the_tube() {
    let p = translation(16.0);
    let domain = Interval::new(16.0, 18.0);
    let mut bad = p;
    bad.axis = [0.0; 3];
    assert!(matches!(
        accelerated_rigid_tube(&bad, domain, 8, 4),
        Err(MotionError::InvalidConfiguration { .. })
    ));
    bad = p;
    bad.linear_acceleration[2] = f64::NAN;
    assert!(matches!(
        accelerated_rigid_tube(&bad, domain, 8, 4),
        Err(MotionError::NonFiniteInput { .. })
    ));
    assert!(matches!(
        accelerated_rigid_tube(&p, domain, 8, 0),
        Err(MotionError::InvalidSegments)
    ));
    let tube = accelerated_rigid_tube(&p, domain, 8, 4).unwrap();
    with_cx(|cx, gate| {
        let point = Point3::new(0.0, 0.0, 0.0);
        let before = tube
            .point_motion_over(point, Interval::point(17.0), cx)
            .unwrap();
        assert!(matches!(
            tube.point_motion_over(point, Interval::new(15.0, 17.0), cx),
            Err(MotionError::OutOfDomain { .. })
        ));
        let after = tube
            .point_motion_over(point, Interval::point(17.0), cx)
            .unwrap();
        assert_eq!(before.acceleration, after.acceleration);
        gate.request();
        assert!(matches!(
            tube.point_motion_over(point, domain, cx),
            Err(MotionError::Cancelled)
        ));
    });
}
