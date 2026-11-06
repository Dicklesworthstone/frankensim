//! G0/G1/G3/G4: analytic carrier velocity through the original scenario frames.

use fs_cheb::Cheb1;
use fs_ga::{Point, Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::frame::{FrameKinematicsBudget, FrameKinematicsError, RATE_DIMS};
use fs_scenario::signal::{ChebProfile, Interp};
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, TimeSignal, WORLD};
use std::ops::ControlFlow;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn keep_going() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

fn near(actual: Vec3, expected: Vec3, tolerance: f64) {
    for (a, b) in [actual.x, actual.y, actual.z]
        .into_iter()
        .zip([expected.x, expected.y, expected.z])
    {
        assert!((a - b).abs() <= tolerance, "{actual:?} != {expected:?}");
    }
}

fn frame(id: u32, parent: FrameId, motion: FrameMotion) -> Frame {
    Frame {
        id: FrameId(id),
        name: format!("frame-{id}"),
        parent,
        motion,
    }
}

fn rotating(axis: [f64; 3], center: Vec3, rate: f64) -> FrameMotion {
    FrameMotion::Rotating {
        axis,
        center,
        rate: QtyAny::new(rate, RATE_DIMS),
    }
}

fn rotate_z(vector: Vec3, angle: f64) -> Vec3 {
    let (s, c) = angle.sin_cos();
    Vec3::new(
        c * vector.x - s * vector.y,
        s * vector.x + c * vector.y,
        vector.z,
    )
}

fn rotate_y(vector: Vec3, angle: f64) -> Vec3 {
    let (s, c) = angle.sin_cos();
    Vec3::new(
        c * vector.x + s * vector.z,
        vector.y,
        -s * vector.x + c * vector.z,
    )
}

fn tilt_tree(angle: TimeSignal) -> FrameTree {
    let mut tree = FrameTree::new();
    tree.add(frame(
        1,
        WORLD,
        FrameMotion::Tilt {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::default(),
            angle,
        },
    ));
    tree
}

#[test]
fn chebyshev_translation_and_tilt_expose_exact_second_derivatives() {
    let mut tree = FrameTree::new();
    tree.add(frame(
        1,
        WORLD,
        FrameMotion::Tilt {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::default(),
            angle: TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![0.5, 0.0, 0.25]),
                dims: Dims::NONE,
            }),
        },
    ));
    tree.add(frame(
        2,
        FrameId(1),
        FrameMotion::Translating {
            direction: [1.0, 0.0, 0.0],
            displacement: TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 0.0, 0.5]),
                dims: LENGTH,
            }),
        },
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(2),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let time = 0.25;
    let sample = path.sample(time, &mut keep_going).unwrap();
    let theta = 0.5 * time * time;
    let omega = time;
    let expected_origin = rotate_z(Vec3::new(time * time, 0.0, 0.0), theta);
    let rotated_v = rotate_z(Vec3::new(2.0 * time, 0.0, 0.0), theta);
    let rotated_a = rotate_z(Vec3::new(2.0, 0.0, 0.0), theta);
    let spin = Vec3::new(0.0, 0.0, omega);
    let expected_acceleration = rotated_a
        + Vec3::new(0.0, 0.0, 1.0).cross(expected_origin)
        + spin.scale(2.0).cross(rotated_v)
        + spin.cross(spin.cross(expected_origin));
    near(sample.origin_world(), expected_origin, 3e-13);
    near(sample.origin_acceleration_world(), expected_acceleration, 5e-13);
    near(
        sample.angular_acceleration_world(),
        Vec3::new(0.0, 0.0, 1.0),
        3e-13,
    );
}

#[test]
fn point_acceleration_contains_centripetal_coriolis_and_relative_terms() {
    let mut tree = FrameTree::new();
    tree.add(frame(
        1,
        WORLD,
        rotating([0.0, 0.0, 1.0], Vec3::default(), 2.0),
    ));
    let sample = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap()
        .sample(0.0, &mut keep_going)
        .unwrap();
    let point = Vec3::new(1.0, 0.0, 0.0);
    near(
        sample
            .material_acceleration_at_world_point(point)
            .unwrap(),
        Vec3::new(-4.0, 0.0, 0.0),
        1e-14,
    );
    let motion = sample
        .point_motion_in_world(
            point,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 3.0, 0.0),
        )
        .unwrap();
    near(motion.acceleration_world, Vec3::new(-4.0, 7.0, 0.0), 1e-14);
}

#[test]
fn acceleration_discontinuities_refuse_instead_of_selecting_a_temporal_side() {
    let tree = tilt_tree(TimeSignal::Ramp {
        t_start: -1.0,
        t_end: 1.0,
        from: QtyAny::new(-1.0, Dims::NONE),
        to: QtyAny::new(1.0, Dims::NONE),
    });
    let path = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    assert_eq!(
        path.sample(0.0, &mut keep_going)
            .unwrap()
            .angular_acceleration_world(),
        Vec3::default()
    );
    for time in [-1.0, 1.0] {
        assert!(matches!(
            path.sample(time, &mut keep_going),
            Err(FrameKinematicsError::Nondifferentiable { .. })
        ));
    }

    let mut endpoint = FrameTree::new();
    endpoint.add(frame(
        1,
        WORLD,
        FrameMotion::Tilt {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::default(),
            // theta=t² on [0,1]: zero endpoint rate at t=0 but nonzero
            // acceleration, so clamping makes the second derivative discontinuous.
            angle: TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(0.0, 1.0, vec![0.75, 0.5, 0.125]),
                dims: Dims::NONE,
            }),
        },
    ));
    let path = endpoint
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    assert!(matches!(
        path.sample(0.0, &mut keep_going),
        Err(FrameKinematicsError::Nondifferentiable { .. })
    ));
}

#[test]
fn fixed_coordinates_do_not_create_a_carrier_velocity() {
    let world = FrameTree::new();
    let path = world
        .prepare_kinematics(WORLD, FrameKinematicsBudget::default(), &mut keep_going)
        .unwrap();
    let point = Vec3::new(2.0, -1.0, 0.5);
    let relative = Vec3::new(3.0, 4.0, 5.0);
    let state = path
        .sample(9.0, &mut keep_going)
        .unwrap()
        .point_in_world(point, relative)
        .unwrap();
    assert_eq!(state.position_world, point);
    assert_eq!(state.velocity_world, relative);

    let shift = Vec3::new(1e12, -2e12, 3e12);
    let mut tree = FrameTree::new();
    tree.add(frame(
        1,
        WORLD,
        FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: shift,
        },
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let sample = path.sample(2.0, &mut keep_going).unwrap();
    let state = sample.point_in_world(point, relative).unwrap();
    near(state.position_world, shift + Vec3::new(1.0, 2.0, 0.5), 1e-3);
    near(state.velocity_world, Vec3::new(-4.0, 3.0, 5.0), 1e-13);
    assert_eq!(sample.origin_velocity_world(), Vec3::default());
    assert_eq!(sample.angular_velocity_world(), Vec3::default());
    assert_eq!(
        sample
            .velocity_at_world_point(state.position_world)
            .unwrap(),
        Vec3::default()
    );
}

#[test]
fn off_center_rotation_carries_points_under_a_fixed_parent() {
    let center = Vec3::new(1.0, -2.0, 0.5);
    let shift = Vec3::new(4.0, -3.0, 2.0);
    let rate = 2.5;
    let time_s = 0.2;
    let mut tree = FrameTree::new();
    tree.add(frame(
        1,
        WORLD,
        FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 1.0, 0.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: shift,
        },
    ));
    tree.add(frame(
        2,
        FrameId(1),
        rotating([0.0, 0.0, 1.0], center, rate),
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(2),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let sample = path.sample(time_s, &mut keep_going).unwrap();
    let point = Vec3::new(2.0, 1.0, -0.2);
    let relative = Vec3::new(0.1, 0.2, 0.3);
    let offset = rotate_z(point - center, rate * time_s);
    let parent_position = center + offset;
    let parent_velocity =
        rotate_z(relative, rate * time_s) + Vec3::new(-rate * offset.y, rate * offset.x, 0.0);
    let to_world = |p: Vec3| Vec3::new(p.z, p.y, -p.x);
    let state = sample.point_in_world(point, relative).unwrap();
    near(
        state.position_world,
        shift + to_world(parent_position),
        2e-13,
    );
    near(state.velocity_world, to_world(parent_velocity), 2e-13);
    near(
        sample.angular_velocity_world(),
        Vec3::new(rate, 0.0, 0.0),
        2e-14,
    );

    let carrier = sample.point_in_world(point, Vec3::default()).unwrap();
    near(
        sample
            .velocity_at_world_point(carrier.position_world)
            .unwrap(),
        carrier.velocity_world,
        2e-13,
    );
    // A fixed WORLD query is not the point obtained by advecting local coordinates.
    let fixed_world_point = Vec3::new(5.0, -4.0, 3.0);
    let world_center = shift + to_world(center);
    let displacement = fixed_world_point - world_center;
    near(
        sample.velocity_at_world_point(fixed_world_point).unwrap(),
        Vec3::new(0.0, -rate * displacement.z, rate * displacement.y),
        2e-13,
    );
}

#[test]
fn translating_child_and_noncommuting_ancestors_compose_analytic_rates() {
    let mut tree = FrameTree::new();
    let a = 0.7;
    let b = -0.9;
    tree.add(frame(
        1,
        WORLD,
        rotating([0.0, 0.0, 1.0], Vec3::default(), a),
    ));
    tree.add(frame(
        2,
        FrameId(1),
        rotating([0.0, 1.0, 0.0], Vec3::default(), b),
    ));
    tree.add(frame(
        3,
        FrameId(2),
        FrameMotion::Translating {
            direction: [1.0, 0.0, 0.0],
            displacement: TimeSignal::Ramp {
                t_start: -2.0,
                t_end: 2.0,
                from: QtyAny::new(-2.0, LENGTH),
                to: QtyAny::new(2.0, LENGTH),
            },
        },
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(3),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let point = Vec3::new(0.3, 0.5, -0.2);
    let relative = Vec3::new(0.4, -0.2, 0.3);
    for time_s in [-0.8, 0.0, 0.4, 1.2] {
        let sample = path.sample(time_s, &mut keep_going).unwrap();
        let inner = rotate_y(point + Vec3::new(time_s, 0.0, 0.0), b * time_s);
        let inner_rate = rotate_y(relative + Vec3::new(1.0, 0.0, 0.0), b * time_s)
            + Vec3::new(b * inner.z, 0.0, -b * inner.x);
        let expected_point = rotate_z(inner, a * time_s);
        let expected_velocity = rotate_z(inner_rate, a * time_s)
            + Vec3::new(-a * expected_point.y, a * expected_point.x, 0.0);
        let state = sample.point_in_world(point, relative).unwrap();
        near(state.position_world, expected_point, 3e-13);
        near(state.velocity_world, expected_velocity, 3e-13);
        near(
            sample.angular_velocity_world(),
            Vec3::new(-b * (a * time_s).sin(), b * (a * time_s).cos(), a),
            3e-13,
        );
        let original = tree
            .world_pose(FrameId(3), time_s)
            .unwrap()
            .transform_point(Point::new(point.x, point.y, point.z))
            .unwrap();
        near(
            state.position_world,
            Vec3::new(original.x, original.y, original.z),
            3e-13,
        );
        assert_eq!(sample, path.sample(time_s, &mut keep_going).unwrap());
    }
}

#[test]
fn inverse_motion_recovers_inertial_straight_lines_with_coriolis_signs() {
    let initial_world = Vec3::new(1.0, -0.2, 0.7);
    let velocity_world = Vec3::new(0.3, -0.4, 0.2);
    for omega in [-0.8, 0.0, 0.8] {
        let mut tree = FrameTree::new();
        tree.add(frame(
            1,
            WORLD,
            rotating([0.0, 0.0, 1.0], Vec3::default(), omega),
        ));
        let path = tree
            .prepare_kinematics(
                FrameId(1),
                FrameKinematicsBudget::default(),
                &mut keep_going,
            )
            .unwrap();
        for time_s in [-0.4, 0.0, 0.6] {
            let sample = path.sample(time_s, &mut keep_going).unwrap();
            let position_world = initial_world + velocity_world.scale(time_s);
            // Differentiate Rz(-omega*t) * (initial_world + velocity_world*t)
            // directly, independently of the forward transport implementation.
            let position = rotate_z(position_world, -omega * time_s);
            let rotated_velocity = rotate_z(velocity_world, -omega * time_s);
            let velocity = Vec3::new(
                rotated_velocity.x + omega * position.y,
                rotated_velocity.y - omega * position.x,
                rotated_velocity.z,
            );
            let acceleration = Vec3::new(
                2.0 * omega * rotated_velocity.y - omega * omega * position.x,
                -2.0 * omega * rotated_velocity.x - omega * omega * position.y,
                0.0,
            );
            let relative = sample
                .point_motion_in_frame(position_world, velocity_world, Vec3::default())
                .unwrap();
            near(relative.position_frame, position, 3e-13);
            near(relative.velocity_relative, velocity, 3e-13);
            near(relative.acceleration_relative, acceleration, 3e-13);
            let recovered = sample
                .point_motion_in_world(
                    relative.position_frame,
                    relative.velocity_relative,
                    relative.acceleration_relative,
                )
                .unwrap();
            near(recovered.position_world, position_world, 3e-13);
            near(recovered.velocity_world, velocity_world, 3e-13);
            near(recovered.acceleration_world, Vec3::default(), 3e-13);
        }
    }
}

#[test]
fn inverse_motion_removes_nested_euler_translation_and_coriolis_transport() {
    // Independent composition: Rz(t²/2) Ry(b*t) (q + t² e_x).
    let mut tree = tilt_tree(TimeSignal::Chebfun(ChebProfile {
        cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![0.5, 0.0, 0.25]),
        dims: Dims::NONE,
    }));
    let b = -0.9;
    tree.add(frame(
        2,
        FrameId(1),
        rotating([0.0, 1.0, 0.0], Vec3::default(), b),
    ));
    tree.add(frame(
        3,
        FrameId(2),
        FrameMotion::Translating {
            direction: [1.0, 0.0, 0.0],
            displacement: TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 0.0, 0.5]),
                dims: LENGTH,
            }),
        },
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(3),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let point = Vec3::new(0.3, 0.5, -0.2);
    let velocity = Vec3::new(0.4, -0.2, 0.3);
    let acceleration = Vec3::new(0.2, 0.1, -0.4);
    for time_s in [-0.5, 0.25, 0.7] {
        let sample = path.sample(time_s, &mut keep_going).unwrap();
        let inner = rotate_y(point + Vec3::new(time_s * time_s, 0.0, 0.0), b * time_s);
        let inner_rotated_velocity =
            rotate_y(velocity + Vec3::new(2.0 * time_s, 0.0, 0.0), b * time_s);
        let inner_velocity =
            inner_rotated_velocity + Vec3::new(b * inner.z, 0.0, -b * inner.x);
        let inner_acceleration =
            rotate_y(acceleration + Vec3::new(2.0, 0.0, 0.0), b * time_s)
                + Vec3::new(
                    2.0 * b * inner_rotated_velocity.z - b * b * inner.x,
                    0.0,
                    -2.0 * b * inner_rotated_velocity.x - b * b * inner.z,
                );
        let angle = 0.5 * time_s * time_s;
        let position_world = rotate_z(inner, angle);
        let rotated_velocity = rotate_z(inner_velocity, angle);
        let velocity_world = rotated_velocity
            + Vec3::new(-time_s * position_world.y, time_s * position_world.x, 0.0);
        let acceleration_world = rotate_z(inner_acceleration, angle)
            + Vec3::new(
                -position_world.y
                    - 2.0 * time_s * rotated_velocity.y
                    - time_s * time_s * position_world.x,
                position_world.x
                    + 2.0 * time_s * rotated_velocity.x
                    - time_s * time_s * position_world.y,
                0.0,
            );
        let relative = sample
            .point_motion_in_frame(position_world, velocity_world, acceleration_world)
            .unwrap();
        near(relative.position_frame, point, 5e-13);
        near(relative.velocity_relative, velocity, 5e-13);
        near(relative.acceleration_relative, acceleration, 5e-13);
    }
}

#[test]
fn inverse_motion_refuses_nonfinite_inputs_and_unrepresentable_transport() {
    let mut tree = FrameTree::new();
    tree.add(frame(
        1,
        WORLD,
        rotating([0.0, 0.0, 1.0], Vec3::default(), 2.0),
    ));
    let sample = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap()
        .sample(0.0, &mut keep_going)
        .unwrap();
    let zero = Vec3::default();
    for bad in [
        Vec3::new(f64::NAN, 0.0, 0.0),
        Vec3::new(0.0, f64::INFINITY, 0.0),
    ] {
        for (point, velocity, acceleration) in
            [(bad, zero, zero), (zero, bad, zero), (zero, zero, bad)]
        {
            assert!(matches!(
                sample.point_motion_in_frame(point, velocity, acceleration),
                Err(FrameKinematicsError::NonFinite(_))
            ));
        }
    }
    assert!(matches!(
        sample.point_motion_in_frame(Vec3::new(f64::MAX, 0.0, 0.0), zero, zero),
        Err(FrameKinematicsError::NonFinite(_))
    ));
}

#[test]
fn clamped_histories_use_analytic_rates_and_refuse_corners_or_jumps() {
    let linear = tilt_tree(TimeSignal::Table {
        times: vec![-1.0, 0.0, 1.0],
        values: vec![0.0, 1.0, 2.0],
        dims: Dims::NONE,
        interp: Interp::Linear,
    });
    let path = linear
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    for t in [-0.5, 0.0, 0.5] {
        near(
            path.sample(t, &mut keep_going)
                .unwrap()
                .angular_velocity_world(),
            Vec3::new(0.0, 0.0, 1.0),
            1e-14,
        );
    }
    for t in [-1.0, 1.0] {
        assert!(matches!(
            path.sample(t, &mut keep_going),
            Err(FrameKinematicsError::Nondifferentiable { .. })
        ));
    }
    assert_eq!(
        path.sample(2.0, &mut keep_going)
            .unwrap()
            .angular_velocity_world(),
        Vec3::default()
    );

    let held = tilt_tree(TimeSignal::Table {
        times: vec![-1.0, 0.0, 1.0],
        values: vec![0.0, 1.0, 1.0],
        dims: Dims::NONE,
        interp: Interp::Hold,
    });
    let path = held
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    for t in [-0.0, 0.0] {
        assert!(matches!(
            path.sample(t, &mut keep_going),
            Err(FrameKinematicsError::Nondifferentiable { .. })
        ));
    }
    assert_eq!(
        path.sample(1.0, &mut keep_going)
            .unwrap()
            .angular_velocity_world(),
        Vec3::default()
    );

    // theta(t) = t^2 on [-1,1], represented directly in native Chebyshev coefficients.
    let smooth = tilt_tree(TimeSignal::Chebfun(ChebProfile {
        cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 0.0, 0.5]),
        dims: Dims::NONE,
    }));
    let path = smooth
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let sample = path.sample(0.25, &mut keep_going).unwrap();
    near(
        sample.angular_velocity_world(),
        Vec3::new(0.0, 0.0, 0.5),
        1e-14,
    );
    let point = sample
        .point_in_world(Vec3::new(1.0, 0.0, 0.0), Vec3::default())
        .unwrap();
    near(
        point.velocity_world,
        Vec3::new(-0.5 * 0.0625_f64.sin(), 0.5 * 0.0625_f64.cos(), 0.0),
        1e-14,
    );
    assert!(matches!(
        path.sample(1.0, &mut keep_going),
        Err(FrameKinematicsError::Nondifferentiable { .. })
    ));
    assert_eq!(
        path.sample(2.0, &mut keep_going)
            .unwrap()
            .angular_velocity_world(),
        Vec3::default()
    );

    let mut translated = FrameTree::new();
    translated.add(frame(
        1,
        WORLD,
        FrameMotion::Translating {
            direction: [0.0, 1.0, 0.0],
            displacement: TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 0.0, 0.5]),
                dims: LENGTH,
            }),
        },
    ));
    let path = translated
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let sample = path.sample(0.5, &mut keep_going).unwrap();
    near(sample.origin_world(), Vec3::new(0.0, 0.25, 0.0), 1e-14);
    near(
        sample.origin_velocity_world(),
        Vec3::new(0.0, 1.0, 0.0),
        1e-14,
    );
}

#[test]
fn malformed_paths_rates_and_explicit_limits_refuse_before_publication() {
    let tree = tilt_tree(TimeSignal::Table {
        times: vec![0.0, 1.0],
        values: vec![0.0, 1.0],
        dims: Dims::NONE,
        interp: Interp::Linear,
    });
    assert!(matches!(
        tree.prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget {
                max_frames: 0,
                max_history_scalars: 16
            },
            &mut keep_going
        ),
        Err(FrameKinematicsError::Budget {
            resource: "frame rows",
            ..
        })
    ));
    assert!(matches!(
        tree.prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget {
                max_frames: 1,
                max_history_scalars: 3
            },
            &mut keep_going
        ),
        Err(FrameKinematicsError::Budget {
            resource: "history scalars",
            ..
        })
    ));
    assert!(
        tree.prepare_kinematics(
            FrameId(77),
            FrameKinematicsBudget::default(),
            &mut keep_going
        )
        .is_err()
    );
    let mut duplicate = tree.clone();
    duplicate.add(tree.frames[0].clone());
    assert!(
        duplicate
            .prepare_kinematics(
                FrameId(1),
                FrameKinematicsBudget::default(),
                &mut keep_going
            )
            .is_err()
    );
    let mut cycle = tree.clone();
    cycle.frames[0].parent = FrameId(1);
    assert!(
        cycle
            .prepare_kinematics(
                FrameId(1),
                FrameKinematicsBudget::default(),
                &mut keep_going
            )
            .is_err()
    );
    assert!(
        tilt_tree(TimeSignal::Constant(QtyAny::new(1.0, LENGTH)))
            .prepare_kinematics(
                FrameId(1),
                FrameKinematicsBudget::default(),
                &mut keep_going
            )
            .is_err()
    );

    let overflow = tilt_tree(TimeSignal::Chebfun(ChebProfile {
        cheb: Cheb1::from_coeffs(0.0, f64::from_bits(1), vec![0.0, 1.0]),
        dims: Dims::NONE,
    }));
    assert!(matches!(
        overflow.prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going
        ),
        Err(FrameKinematicsError::Derivative(_))
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    assert!(matches!(
        path.sample(f64::NAN, &mut keep_going),
        Err(FrameKinematicsError::NonFinite(_))
    ));
    let sample = path.sample(0.5, &mut keep_going).unwrap();
    assert!(
        sample
            .velocity_at_world_point(Vec3::new(f64::INFINITY, 0.0, 0.0))
            .is_err()
    );
    assert!(
        sample
            .point_in_world(Vec3::default(), Vec3::new(0.0, f64::NAN, 0.0))
            .is_err()
    );
}

#[test]
fn admission_and_time_sampling_are_cancellable_without_revalidating_tables() {
    let tree = tilt_tree(TimeSignal::Table {
        times: (0..128).map(|i| f64::from(i)).collect(),
        values: (0..128).map(|i| f64::from(i) * 0.01).collect(),
        dims: Dims::NONE,
        interp: Interp::Linear,
    });
    let mut polls = 0;
    assert!(matches!(
        tree.prepare_kinematics(FrameId(1), FrameKinematicsBudget::default(), &mut || {
            polls += 1;
            if polls == 5 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        }),
        Err(FrameKinematicsError::Cancelled)
    ));
    let path = tree
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let mut query_polls = 0;
    path.sample(2.5, &mut || {
        query_polls += 1;
        ControlFlow::Continue(())
    })
    .unwrap();
    assert!(query_polls < 8, "query rescanned admitted signal data");
    assert!(matches!(
        path.sample(2.5, &mut || ControlFlow::Break(())),
        Err(FrameKinematicsError::Cancelled)
    ));

    let smooth = tilt_tree(TimeSignal::Chebfun(ChebProfile {
        cheb: Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 0.1, 0.05, 0.025, 0.01, 0.005]),
        dims: Dims::NONE,
    }));
    let path = smooth
        .prepare_kinematics(
            FrameId(1),
            FrameKinematicsBudget::default(),
            &mut keep_going,
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        path.sample(0.2, &mut || {
            polls += 1;
            if polls == 4 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        }),
        Err(FrameKinematicsError::Cancelled)
    ));
}
