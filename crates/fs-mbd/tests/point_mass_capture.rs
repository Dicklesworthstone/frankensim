use fs_mbd::assembly::capture::{PointMassCapture, PointMassCaptureResult, capture_point_mass};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn controls() -> AssemblyControls {
    AssemblyControls {
        maximum_parts: 2,
        maximum_impulses: 1,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-8,
        energy_tolerance_j: 1e-9,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.17e} != {b:.17e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    close(a.x, b.x, tol);
    close(a.y, b.y, tol);
    close(a.z, b.z, tol);
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
fn state(p: MassProperties, position: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn sphere() -> MassProperties {
    MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}
fn parcel() -> PointMassCapture {
    PointMassCapture {
        mass_kg: 1.0,
        attachment_arm_body_m: Vec3::new(1.0, 0.0, 0.0),
        velocity_world_m_s: Vec3::new(0.0, 3.0, 0.0),
    }
}
fn resting() -> RigidBodyState {
    state(
        sphere(),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    )
}
fn capture(
    b: &RigidBodyState,
    p: MassProperties,
    incoming: PointMassCapture,
) -> PointMassCaptureResult {
    capture_point_mass(b, p, incoming, controls(), || false).unwrap()
}
// Independent world-frame sums at a fixed observer; no returned balance is
// used to calculate the expected momentum or capture energy.
fn check_conservation(
    b: RigidBodyState,
    p: MassProperties,
    incoming: PointMassCapture,
    out: PointMassCaptureResult,
) {
    let observer = Vec3::new(-0.7, 1.1, 0.4);
    let momentum = incoming.velocity_world_m_s.scale(incoming.mass_kg);
    let point = b
        .pose()
        .point_world_from_body(incoming.attachment_arm_body_m)
        .unwrap();
    let before_l = b
        .pose()
        .orientation()
        .rotate_body_to_world(b.angular_momentum_body())
        .add(
            b.pose()
                .position_world()
                .sub(observer)
                .cross(b.linear_momentum_world()),
        )
        .add(point.sub(observer).cross(momentum));
    let after_l = out
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(out.state.angular_momentum_body())
        .add(
            out.state
                .pose()
                .position_world()
                .sub(observer)
                .cross(out.state.linear_momentum_world()),
        );
    vector(
        out.state.linear_momentum_world(),
        b.linear_momentum_world().add(momentum),
        2e-12,
    );
    vector(after_l, before_l, 3e-11);
    close(
        out.state.kinetic_energy(out.properties).unwrap() + out.removed_relative_energy_j,
        b.kinetic_energy(p).unwrap()
            + 0.5 * incoming.mass_kg * incoming.velocity_world_m_s.norm_squared(),
        2e-10,
    );
}

#[test]
fn central_capture_matches_reduced_mass_and_retains_spin() {
    let p = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let b = state(
        p,
        Vec3::new(2.0, -1.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(1.0, 2.0, -1.0),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let incoming = PointMassCapture {
        mass_kg: 1.0,
        attachment_arm_body_m: Vec3::ZERO,
        velocity_world_m_s: Vec3::new(5.0, -2.0, 1.0),
    };
    let out = capture(&b, p, incoming);
    assert_eq!(out.properties.mass(), 4.0);
    assert_eq!(out.state.pose(), b.pose());
    assert_eq!(out.state.angular_momentum_body(), b.angular_momentum_body());
    assert_eq!(
        out.properties.inertia_tensor_body(),
        p.inertia_tensor_body()
    );
    vector(
        out.state
            .center_of_mass_velocity_world(out.properties)
            .unwrap(),
        Vec3::new(2.0, 1.0, -0.5),
        0.0,
    );
    vector(
        out.captured_impulse_world_ns,
        Vec3::new(3.0, -3.0, 1.5),
        2e-15,
    );
    close(out.removed_relative_energy_j, 13.5, 1e-14);
    check_conservation(b, p, incoming, out);
}

#[test]
fn off_center_capture_matches_analytic_recoil_and_parallel_axis_inertia() {
    let b = resting();
    let incoming = parcel();
    let out = capture(&b, sphere(), incoming);
    vector(out.center_shift_body_m, Vec3::new(0.5, 0.0, 0.0), 0.0);
    vector(
        out.state
            .center_of_mass_velocity_world(out.properties)
            .unwrap(),
        Vec3::new(0.0, 1.5, 0.0),
        0.0,
    );
    vector(
        out.state.angular_momentum_body(),
        Vec3::new(0.0, 0.0, 1.5),
        5e-16,
    );
    vector(
        out.properties
            .angular_velocity_body_checked(out.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(0.0, 0.0, 1.0),
        5e-16,
    );
    let expected = [1.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.5];
    for (a, e) in out
        .properties
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(expected)
    {
        close(a, e, 5e-16);
    }
    vector(
        out.captured_impulse_world_ns,
        Vec3::new(0.0, 1.0, 0.0),
        5e-16,
    );
    close(out.removed_relative_energy_j, 1.5, 1e-15);
    let arm_after = incoming.attachment_arm_body_m.sub(out.center_shift_body_m);
    vector(
        out.state.pose().point_world_from_body(arm_after).unwrap(),
        b.pose()
            .point_world_from_body(incoming.attachment_arm_body_m)
            .unwrap(),
        0.0,
    );
    check_conservation(b, sphere(), incoming, out);
}

#[test]
fn full_tensor_capture_is_covariant_and_conserves_independent_world_totals() {
    let p = MassProperties::from_principal_frame(
        3.0,
        Vec3::ZERO,
        Vec3::new(2.0, 3.0, 4.0),
        rotation(Vec3::new(1.0, -0.4, 0.6), 0.7),
    )
    .unwrap();
    let q = rotation(Vec3::new(-0.3, 0.8, 1.0), -0.4);
    let v = Vec3::new(0.3, -0.8, 1.1);
    let b = state(
        p,
        Vec3::new(2.0, -1.0, 0.4),
        q,
        v,
        Vec3::new(0.5, -0.4, 0.7),
    );
    let incoming = PointMassCapture {
        mass_kg: 0.75,
        attachment_arm_body_m: Vec3::new(0.8, -0.4, 0.3),
        velocity_world_m_s: Vec3::new(-0.7, 1.2, 0.5),
    };
    let out = capture(&b, p, incoming);
    // The captured central tensor is far from diagonal in the body axes. By an
    // independent parallel-axis sum, I_xz = -0.517 and I_yz = -0.463; I_xy is
    // only 0.026 for this fixture, so it is not the entry to check.
    let tensor = out.properties.inertia_tensor_body().m;
    assert!(tensor[2].abs() > 0.4 && tensor[5].abs() > 0.4, "{tensor:?}");
    check_conservation(b, p, incoming, out);
    let turn = rotation(Vec3::new(0.4, -0.2, 0.9), 1.1);
    let offset = Vec3::new(4.0, -2.0, 3.0);
    let boost = Vec3::new(64.0, -32.0, 16.0);
    let moved = RigidBodyState::new(
        Pose::new(
            turn.rotate_body_to_world(b.pose().position_world())
                .add(offset),
            compose(turn, b.pose().orientation()),
        )
        .unwrap(),
        turn.rotate_body_to_world(b.linear_momentum_world())
            .add(boost.scale(p.mass())),
        b.angular_momentum_body(),
    )
    .unwrap();
    let moved_incoming = PointMassCapture {
        velocity_world_m_s: turn
            .rotate_body_to_world(incoming.velocity_world_m_s)
            .add(boost),
        ..incoming
    };
    let transformed = capture(&moved, p, moved_incoming);
    vector(
        transformed.state.pose().position_world(),
        turn.rotate_body_to_world(out.state.pose().position_world())
            .add(offset),
        3e-14,
    );
    vector(
        transformed.state.linear_momentum_world(),
        turn.rotate_body_to_world(out.state.linear_momentum_world())
            .add(boost.scale(out.properties.mass())),
        6e-14,
    );
    vector(
        transformed.state.angular_momentum_body(),
        out.state.angular_momentum_body(),
        3e-14,
    );
    vector(
        transformed.captured_impulse_world_ns,
        turn.rotate_body_to_world(out.captured_impulse_world_ns),
        3e-14,
    );
    close(
        transformed.removed_relative_energy_j,
        out.removed_relative_energy_j,
        3e-14,
    );
    assert_eq!(transformed.center_shift_body_m, out.center_shift_body_m);
    assert_eq!(transformed.properties, out.properties);
    check_conservation(moved, p, moved_incoming, transformed);
}

#[test]
fn capture_loss_survives_large_common_translation() {
    let plain = capture(&resting(), sphere(), parcel());
    let boost = Vec3::new(1_099_511_627_776.0, 0.0, 0.0);
    let b = state(
        sphere(),
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        boost,
        Vec3::ZERO,
    );
    let incoming = PointMassCapture {
        velocity_world_m_s: parcel().velocity_world_m_s.add(boost),
        ..parcel()
    };
    let out = capture(&b, sphere(), incoming);
    assert_eq!(
        out.removed_relative_energy_j,
        plain.removed_relative_energy_j
    );
    assert_eq!(
        out.captured_impulse_world_ns,
        plain.captured_impulse_world_ns
    );
    assert_eq!(
        out.state.angular_momentum_body(),
        plain.state.angular_momentum_body()
    );
    close(out.state.linear_momentum_world().x, 2.0 * boost.x, 0.0);
    close(out.removed_relative_energy_j, 1.5, 1e-15);
    // Subtracting total lab-frame energies cannot resolve this loss.
    assert_eq!(
        b.kinetic_energy(sphere()).unwrap(),
        b.kinetic_energy(sphere()).unwrap() + out.removed_relative_energy_j
    );
}

#[test]
fn native_flight_and_co_rotating_recapture_preserve_material_point_meaning() {
    let first = capture(&resting(), sphere(), parcel());
    let carried_arm = parcel()
        .attachment_arm_body_m
        .sub(first.center_shift_body_m);
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let elapsed = 0.2;
    let flight = integrator
        .step(first.state, first.properties, Wrench::ZERO, elapsed)
        .unwrap()
        .state_after;
    let point = flight
        .point_kinematics(first.properties, carried_arm)
        .unwrap();
    vector(
        point.point_world,
        Vec3::new(
            0.5 + 0.5 * elapsed.cos(),
            1.5 * elapsed + 0.5 * elapsed.sin(),
            0.0,
        ),
        2e-14,
    );
    let incoming = PointMassCapture {
        mass_kg: 1.0,
        attachment_arm_body_m: carried_arm,
        velocity_world_m_s: point.point_velocity_world,
    };
    let second = capture(&flight, first.properties, incoming);
    assert!(second.removed_relative_energy_j < 1e-28);
    vector(
        second
            .state
            .pose()
            .point_world_from_body(carried_arm.sub(second.center_shift_body_m))
            .unwrap(),
        point.point_world,
        2e-15,
    );
    vector(
        second
            .properties
            .angular_velocity_body_checked(second.state.angular_momentum_body())
            .unwrap(),
        first
            .properties
            .angular_velocity_body_checked(flight.angular_momentum_body())
            .unwrap(),
        2e-14,
    );
    check_conservation(flight, first.properties, incoming, second);
}

#[test]
fn zero_mass_is_identity_and_refusals_publish_no_changed_owner() {
    let b = resting();
    let p = sphere();
    let original = capture(&b, p, parcel());
    let mut one_part = controls();
    one_part.maximum_parts = 1;
    one_part.maximum_impulses = 0;
    let empty = capture_point_mass(
        &b,
        p,
        PointMassCapture {
            mass_kg: 0.0,
            ..parcel()
        },
        one_part,
        || false,
    )
    .unwrap();
    assert_eq!(empty.state, b);
    assert_eq!(empty.properties, p);
    assert_eq!(empty.center_shift_body_m, Vec3::ZERO);
    assert_eq!(empty.removed_relative_energy_j, 0.0);
    assert!(matches!(
        capture_point_mass(&b, p, parcel(), one_part, || false),
        Err(AssemblyError::Budget { .. })
    ));
    for mass_kg in [-1.0, f64::INFINITY, f64::EPSILON / 4.0] {
        assert!(
            capture_point_mass(
                &b,
                p,
                PointMassCapture {
                    mass_kg,
                    ..parcel()
                },
                controls(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        capture_point_mass(
            &b,
            p,
            PointMassCapture {
                velocity_world_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
                ..parcel()
            },
            controls(),
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        capture_point_mass(&b, p, parcel(), controls(), || {
            polls += 1;
            polls == 3
        }),
        Err(AssemblyError::Cancelled)
    );
    assert_eq!(polls, 3);
    assert_eq!(b, resting());
    assert_eq!(p, sphere());
    assert_eq!(capture(&b, p, parcel()), original);
}
