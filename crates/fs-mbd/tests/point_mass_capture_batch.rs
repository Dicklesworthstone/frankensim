//! Simultaneous capture: shared recoil, spatial momentum, and projection loss.

use fs_mbd::assembly::capture::{
    PointMassCapture, PointMassCaptureResult, capture_point_mass, capture_point_masses,
};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls(count: usize) -> AssemblyControls {
    AssemblyControls {
        maximum_parts: count + 1,
        maximum_impulses: count,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-8,
        energy_tolerance_j: 1e-9,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn inertia(a: MassProperties, b: MassProperties, tolerance: f64) {
    close(a.mass(), b.mass(), tolerance);
    for (x, y) in a
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(b.inertia_tensor_body().m)
    {
        close(x, y, tolerance);
    }
}
fn rotation(axis: Vec3, angle: f64) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(axis, angle).unwrap()
}
fn compose(q: UnitQuaternion, r: UnitQuaternion) -> UnitQuaternion {
    let [a, b, c, d] = q.components();
    let [e, f, g, h] = r.components();
    UnitQuaternion::new(
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )
    .unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn parcel(mass: f64, arm: Vec3, velocity: Vec3) -> PointMassCapture {
    PointMassCapture {
        mass_kg: mass,
        attachment_arm_body_m: arm,
        velocity_world_m_s: velocity,
    }
}
fn capture(
    body: &RigidBodyState,
    properties: MassProperties,
    incoming: &[PointMassCapture],
) -> PointMassCaptureResult {
    capture_point_masses(body, properties, incoming, controls(incoming.len()), || {
        false
    })
    .unwrap()
}
// Reconstruct physical totals at a separate fixed observer using native endpoint
// state and input parcel data, independently of the event's balance receipt.
fn conservation(
    before: RigidBodyState,
    properties: MassProperties,
    incoming: &[PointMassCapture],
    out: PointMassCaptureResult,
) {
    let origin = Vec3::new(-0.7, 1.1, 0.4);
    let mut linear = before.linear_momentum_world();
    let mut angular = before
        .pose()
        .orientation()
        .rotate_body_to_world(before.angular_momentum_body())
        .add(before.pose().position_world().sub(origin).cross(linear));
    let mut kinetic = before.kinetic_energy(properties).unwrap();
    for p in incoming {
        let momentum = p.velocity_world_m_s.scale(p.mass_kg);
        linear = linear.add(momentum);
        angular = angular.add(
            before
                .pose()
                .point_world_from_body(p.attachment_arm_body_m)
                .unwrap()
                .sub(origin)
                .cross(momentum),
        );
        kinetic += 0.5 * p.mass_kg * p.velocity_world_m_s.dot(p.velocity_world_m_s);
    }
    vector(out.state.linear_momentum_world(), linear, 2e-11);
    let after_angular = out
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(out.state.angular_momentum_body())
        .add(
            out.state
                .pose()
                .position_world()
                .sub(origin)
                .cross(out.state.linear_momentum_world()),
        );
    vector(after_angular, angular, 4e-10);
    close(
        out.state.kinetic_energy(out.properties).unwrap() + out.removed_relative_energy_j,
        kinetic,
        5e-9,
    );
}

#[test]
fn simultaneous_central_streams_include_shared_translation_cross_terms() {
    let properties = MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(1.0, 1.5, 2.0),
        rotation(Vec3::new(1.0, -0.3, 0.5), 0.4),
    )
    .unwrap();
    let body = state(
        properties,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::new(0.3, -0.5, 0.7),
    );
    for (sign, expected_loss, expected_p, expected_impulse) in [
        (-1.0, 4.0, Vec3::ZERO, Vec3::ZERO),
        (1.0, 2.0, Vec3::new(4.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)),
    ] {
        let incoming = [
            parcel(1.0, Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0)),
            parcel(1.0, Vec3::ZERO, Vec3::new(2.0 * sign, 0.0, 0.0)),
        ];
        let out = capture(&body, properties, &incoming);
        assert_eq!(out.properties.mass(), 4.0);
        assert_eq!(out.center_shift_body_m, Vec3::ZERO);
        assert_eq!(
            out.properties.inertia_tensor_body(),
            properties.inertia_tensor_body()
        );
        assert_eq!(
            out.state.angular_momentum_body(),
            body.angular_momentum_body()
        );
        vector(out.state.linear_momentum_world(), expected_p, 0.0);
        vector(out.captured_impulse_world_ns, expected_impulse, 0.0);
        close(out.removed_relative_energy_j, expected_loss, 2e-14);
        let independently_added: f64 = incoming
            .iter()
            .map(|p| {
                capture_point_mass(&body, properties, *p, controls(1), || false)
                    .unwrap()
                    .removed_relative_energy_j
            })
            .sum();
        assert!((out.removed_relative_energy_j - independently_added).abs() > 0.5);
        conservation(body, properties, &incoming, out);
    }
}

#[test]
fn off_center_capture_updates_com_inertia_and_spin_from_one_original_snapshot() {
    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 2.0)).unwrap();
    let body = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let incoming = [
        parcel(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 2.0, 0.0)),
        parcel(1.0, Vec3::new(3.0, 0.0, 0.0), Vec3::ZERO),
    ];
    let out = capture(&body, p, &incoming);
    inertia(
        out.properties,
        MassProperties::new(4.0, Vec3::ZERO, Vec3::new(1.0, 8.0, 8.0)).unwrap(),
        3e-14,
    );
    vector(out.center_shift_body_m, Vec3::new(1.0, 0.0, 0.0), 0.0);
    vector(out.state.angular_momentum_body(), Vec3::ZERO, 0.0);
    vector(out.captured_impulse_world_ns, Vec3::new(0.0, 1.0, 0.0), 0.0);
    close(out.removed_relative_energy_j, 1.5, 2e-14);
    for parcel in incoming {
        vector(
            out.state
                .pose()
                .point_world_from_body(parcel.attachment_arm_body_m.sub(out.center_shift_body_m))
                .unwrap(),
            body.pose()
                .point_world_from_body(parcel.attachment_arm_body_m)
                .unwrap(),
            0.0,
        );
    }
    conservation(body, p, &incoming, out);

    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let velocity = Vec3::new(1.0, -2.0, 0.0);
    let body = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        velocity,
        Vec3::new(0.0, 0.0, 1.0),
    );
    let incoming = [
        parcel(
            1.0,
            Vec3::new(1.0, 0.0, 0.0),
            velocity.add(Vec3::new(0.0, 3.0, 0.0)),
        ),
        parcel(
            1.0,
            Vec3::new(-1.0, 0.0, 0.0),
            velocity.add(Vec3::new(0.0, -1.0, 0.0)),
        ),
    ];
    let out = capture(&body, p, &incoming);
    inertia(
        out.properties,
        MassProperties::new(4.0, Vec3::ZERO, Vec3::new(2.0, 5.0, 6.0)).unwrap(),
        3e-14,
    );
    vector(
        out.state.linear_momentum_world(),
        Vec3::new(4.0, -6.0, 0.0),
        0.0,
    );
    vector(
        out.state.angular_momentum_body(),
        Vec3::new(0.0, 0.0, 8.0),
        0.0,
    );
    close(out.removed_relative_energy_j, 7.0 / 6.0, 3e-14);
    conservation(body, p, &incoming, out);
    let corotating = incoming.map(|item| PointMassCapture {
        velocity_world_m_s: body
            .point_kinematics(p, item.attachment_arm_body_m)
            .unwrap()
            .point_velocity_world,
        ..item
    });
    let locked = capture(&body, p, &corotating);
    assert_eq!(locked.removed_relative_energy_j, 0.0);
    vector(
        locked
            .properties
            .angular_velocity_body_checked(locked.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(0.0, 0.0, 1.0),
        0.0,
    );
    conservation(body, p, &corotating, locked);
}

fn tensor_fixture() -> (MassProperties, RigidBodyState, [PointMassCapture; 3]) {
    let p = MassProperties::from_principal_frame(
        3.5,
        Vec3::ZERO,
        Vec3::new(2.0, 3.0, 4.0),
        rotation(Vec3::new(1.0, 0.3, -0.7), 0.43),
    )
    .unwrap();
    let body = state(
        p,
        Vec3::new(2.0, -1.0, 0.4),
        rotation(Vec3::new(-0.3, 0.8, 1.0), -0.4),
        Vec3::new(1.0, -0.5, 0.25),
        Vec3::new(0.4, -0.3, 0.6),
    );
    let incoming = [
        parcel(0.5, Vec3::new(0.5, 0.2, -0.3), Vec3::new(-0.2, 0.7, 1.2)),
        parcel(1.25, Vec3::new(-0.3, 0.6, 0.4), Vec3::new(1.5, -0.4, 0.2)),
        parcel(0.75, Vec3::new(0.4, -0.4, 0.2), Vec3::new(-0.6, -1.1, 0.4)),
    ];
    (p, body, incoming)
}

#[test]
fn full_tensor_batch_is_covariant_and_permutation_equivalent() {
    let (p, body, incoming) = tensor_fixture();
    let out = capture(&body, p, &incoming);
    assert!(out.properties.inertia_tensor_body().m[1].abs() > 0.01);
    conservation(body, p, &incoming, out);
    for order in [[2, 1, 0], [1, 2, 0]] {
        let reordered = order.map(|i| incoming[i]);
        let other = capture(&body, p, &reordered);
        inertia(other.properties, out.properties, 2e-13);
        vector(
            other.state.linear_momentum_world(),
            out.state.linear_momentum_world(),
            3e-14,
        );
        vector(
            other.state.angular_momentum_body(),
            out.state.angular_momentum_body(),
            2e-13,
        );
        close(
            other.removed_relative_energy_j,
            out.removed_relative_energy_j,
            3e-13,
        );
    }
    let turn = rotation(Vec3::new(0.4, -0.2, 0.9), 1.1);
    let offset = Vec3::new(4.0, -2.0, 3.0);
    let boost = Vec3::new(64.0, -32.0, 16.0);
    let moved = RigidBodyState::new(
        Pose::new(
            turn.rotate_body_to_world(body.pose().position_world())
                .add(offset),
            compose(turn, body.pose().orientation()),
        )
        .unwrap(),
        turn.rotate_body_to_world(body.linear_momentum_world())
            .add(boost.scale(p.mass())),
        body.angular_momentum_body(),
    )
    .unwrap();
    let moved_parcels = incoming.map(|item| PointMassCapture {
        velocity_world_m_s: turn
            .rotate_body_to_world(item.velocity_world_m_s)
            .add(boost),
        ..item
    });
    let transformed = capture(&moved, p, &moved_parcels);
    assert_eq!(transformed.properties, out.properties);
    vector(
        transformed.center_shift_body_m,
        out.center_shift_body_m,
        0.0,
    );
    vector(
        transformed.state.pose().position_world(),
        turn.rotate_body_to_world(out.state.pose().position_world())
            .add(offset),
        4e-14,
    );
    vector(
        transformed.state.linear_momentum_world(),
        turn.rotate_body_to_world(out.state.linear_momentum_world())
            .add(boost.scale(out.properties.mass())),
        2e-13,
    );
    vector(
        transformed.state.angular_momentum_body(),
        out.state.angular_momentum_body(),
        3e-13,
    );
    close(
        transformed.removed_relative_energy_j,
        out.removed_relative_energy_j,
        3e-13,
    );
    conservation(moved, p, &moved_parcels, transformed);
}

#[test]
fn batch_matches_rebased_sequential_capture_and_preserves_singleton_receipt() {
    let (p, body, incoming) = tensor_fixture();
    let out = capture(&body, p, &incoming);
    let mut current = body;
    let mut properties = p;
    let mut shift = Vec3::ZERO;
    let mut loss = 0.0;
    for item in incoming {
        let single = capture_point_mass(
            &current,
            properties,
            PointMassCapture {
                attachment_arm_body_m: item.attachment_arm_body_m.sub(shift),
                ..item
            },
            controls(1),
            || false,
        )
        .unwrap();
        current = single.state;
        properties = single.properties;
        shift = shift.add(single.center_shift_body_m);
        loss += single.removed_relative_energy_j;
    }
    inertia(out.properties, properties, 2e-12);
    vector(
        out.state.pose().position_world(),
        current.pose().position_world(),
        4e-14,
    );
    vector(
        out.state.linear_momentum_world(),
        current.linear_momentum_world(),
        3e-14,
    );
    vector(
        out.state.angular_momentum_body(),
        current.angular_momentum_body(),
        5e-13,
    );
    close(out.removed_relative_energy_j, loss, 2e-12);
    for item in incoming {
        assert_eq!(
            capture_point_masses(&body, p, &[item], controls(1), || false),
            capture_point_mass(&body, p, item, controls(1), || false),
        );
    }
}

#[test]
fn batch_loss_survives_common_boost_and_lost_com_geometry_is_rejected() {
    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 2.0)).unwrap();
    let body = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let incoming = [
        parcel(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 2.0, 0.0)),
        parcel(1.0, Vec3::new(3.0, 0.0, 0.0), Vec3::ZERO),
    ];
    let plain = capture(&body, p, &incoming);
    let boost = Vec3::new(1_099_511_627_776.0, 0.0, 0.0);
    let boosted = state(p, Vec3::ZERO, UnitQuaternion::IDENTITY, boost, Vec3::ZERO);
    let shifted = incoming.map(|parcel| PointMassCapture {
        velocity_world_m_s: parcel.velocity_world_m_s.add(boost),
        ..parcel
    });
    let out = capture(&boosted, p, &shifted);
    assert_eq!(
        out.removed_relative_energy_j,
        plain.removed_relative_energy_j
    );
    assert_eq!(
        out.state.angular_momentum_body(),
        plain.state.angular_momentum_body()
    );
    assert_eq!(out.balance, plain.balance);
    let lab_energy = out.state.kinetic_energy(out.properties).unwrap();
    assert_eq!(lab_energy, lab_energy + out.removed_relative_energy_j);
    let distant = state(
        p,
        Vec3::new(36_028_797_018_963_968.0, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    assert!(matches!(
        capture_point_masses(&distant, p, &incoming, controls(2), || false),
        Err(AssemblyError::Balance(_)),
    ));
}

#[test]
fn aggregate_mass_budget_validation_and_late_cancellation_are_atomic() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let body = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let tiny = parcel(f64::EPSILON / 8.0, Vec3::ZERO, Vec3::ZERO);
    assert!(capture_point_mass(&body, p, tiny, controls(1), || false).is_err());
    let tiny_batch = [tiny; 8];
    let out = capture(&body, p, &tiny_batch);
    assert_eq!(out.properties.mass(), 1.0 + f64::EPSILON);
    assert_eq!(out.removed_relative_energy_j, 0.0);
    assert_eq!(out.state, body);
    let zero_pair = [parcel(0.0, Vec3::ZERO, Vec3::ZERO); 2];
    for incoming in [&[][..], &zero_pair[..]] {
        let identity = capture(&body, p, incoming);
        assert_eq!(identity.state, body);
        assert_eq!(identity.properties, p);
        assert_eq!(identity.removed_relative_energy_j, 0.0);
    }
    let (p, body, incoming) = tensor_fixture();
    let expected = capture(&body, p, &incoming);
    let before = (body, p, incoming);
    let mut limited = controls(3);
    limited.maximum_parts = 3;
    assert!(matches!(
        capture_point_masses(&body, p, &incoming, limited, || false),
        Err(AssemblyError::Budget { .. }),
    ));
    limited = controls(3);
    limited.maximum_impulses = 2;
    assert!(matches!(
        capture_point_masses(&body, p, &incoming, limited, || false),
        Err(AssemblyError::Budget { .. }),
    ));
    let zeros = [parcel(0.0, Vec3::ZERO, Vec3::ZERO); 3];
    assert!(matches!(
        capture_point_masses(&body, p, &zeros, controls(2), || false),
        Err(AssemblyError::Budget { .. }),
    ));
    for bad in [
        PointMassCapture {
            mass_kg: -1.0,
            ..incoming[0]
        },
        PointMassCapture {
            mass_kg: f64::NAN,
            ..incoming[0]
        },
        PointMassCapture {
            attachment_arm_body_m: Vec3::new(f64::INFINITY, 0.0, 0.0),
            ..incoming[0]
        },
        PointMassCapture {
            velocity_world_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
            ..incoming[0]
        },
    ] {
        let invalid = [incoming[0], bad];
        assert!(capture_point_masses(&body, p, &invalid, controls(2), || false).is_err());
    }
    let mut complete_polls = 0;
    assert_eq!(
        capture_point_masses(&body, p, &incoming, controls(3), || {
            complete_polls += 1;
            false
        })
        .unwrap(),
        expected,
    );
    let mut polls = 0;
    assert_eq!(
        capture_point_masses(&body, p, &incoming, controls(3), || {
            polls += 1;
            polls == complete_polls
        }),
        Err(AssemblyError::Cancelled),
    );
    assert_eq!(polls, complete_polls);
    assert_eq!((body, p, incoming), before);
    assert_eq!(capture(&body, p, &incoming), expected);
}
