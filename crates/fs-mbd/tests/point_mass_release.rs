use fs_mbd::assembly::capture::{PointMassCapture, capture_point_mass};
use fs_mbd::assembly::release::{
    PointMassRelease, PointMassReleaseResult, ReleasedPointMass, release_point_mass,
};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls() -> AssemblyControls {
    AssemblyControls {
        maximum_parts: 2,
        maximum_impulses: 1,
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
fn spinning_fixture() -> (MassProperties, RigidBodyState, PointMassRelease) {
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(2.0, 6.0, 7.0)).unwrap();
    let b = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(5.0, -2.0, 1.0),
        Vec3::new(0.0, 0.0, 2.0),
    );
    let parcel = PointMassRelease {
        mass_kg: 1.0,
        attachment_arm_body_m: Vec3::new(1.5, 0.0, 0.0),
        impulse_on_parcel_body_ns: Vec3::new(0.0, 3.0, 0.0),
    };
    (p, b, parcel)
}
fn release(
    b: &RigidBodyState,
    p: MassProperties,
    parcel: PointMassRelease,
) -> PointMassReleaseResult {
    release_point_mass(b, p, parcel, controls(), || false).unwrap()
}
// Independently sum spin plus orbital momentum at a fixed world observer.
// No returned balance receipt or internal event formula supplies this oracle.
fn conservation(b: RigidBodyState, p: MassProperties, out: PointMassReleaseResult) {
    let parcel = out.parcel.unwrap();
    let observer = Vec3::new(-0.7, 1.1, 0.4);
    let before_l = b
        .pose()
        .orientation()
        .rotate_body_to_world(b.angular_momentum_body())
        .add(
            b.pose()
                .position_world()
                .sub(observer)
                .cross(b.linear_momentum_world()),
        );
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
        )
        .add(
            parcel
                .position_world_m
                .sub(observer)
                .cross(parcel.linear_momentum_world_ns),
        );
    vector(
        out.state
            .linear_momentum_world()
            .add(parcel.linear_momentum_world_ns),
        b.linear_momentum_world(),
        3e-12,
    );
    vector(after_l, before_l, 1e-10);
    let parcel_energy = 0.5
        * parcel
            .linear_momentum_world_ns
            .dot(parcel.velocity_world_m_s().unwrap());
    close(
        out.state.kinetic_energy(out.properties).unwrap() + parcel_energy,
        b.kinetic_energy(p).unwrap() + out.supplied_work_j,
        5e-10,
    );
}

#[test]
fn central_release_matches_analytic_recoil_and_keeps_full_spin() {
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let b = state(
        p,
        Vec3::new(1.0, -2.0, 0.5),
        UnitQuaternion::IDENTITY,
        Vec3::new(2.0, 1.0, 0.0),
        Vec3::new(0.5, 1.0, -0.25),
    );
    let outgoing = PointMassRelease {
        mass_kg: 1.0,
        attachment_arm_body_m: Vec3::ZERO,
        impulse_on_parcel_body_ns: Vec3::new(6.0, -3.0, 0.0),
    };
    let out = release(&b, p, outgoing);
    let parcel = out.parcel.unwrap();
    assert_eq!(out.properties.mass(), 3.0);
    assert_eq!(out.center_shift_body_m, Vec3::ZERO);
    assert_eq!(out.state.pose(), b.pose());
    assert_eq!(
        out.properties.inertia_tensor_body(),
        p.inertia_tensor_body()
    );
    vector(
        out.state.angular_momentum_body(),
        b.angular_momentum_body(),
        1e-15,
    );
    vector(
        out.state.linear_momentum_world(),
        Vec3::new(0.0, 6.0, 0.0),
        0.0,
    );
    vector(
        parcel.linear_momentum_world_ns,
        Vec3::new(8.0, -2.0, 0.0),
        0.0,
    );
    vector(
        parcel.velocity_world_m_s().unwrap(),
        Vec3::new(8.0, -2.0, 0.0),
        0.0,
    );
    assert_eq!(parcel.position_world_m, b.pose().position_world());
    close(out.supplied_work_j, 30.0, 2e-14);
    conservation(b, p, out);
}

#[test]
fn spinning_off_center_release_matches_exact_mass_inertia_and_work() {
    let (p, b, outgoing) = spinning_fixture();
    let out = release(&b, p, outgoing);
    let parcel = out.parcel.unwrap();
    let retained = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    inertia(out.properties, retained, 3e-14);
    vector(out.center_shift_body_m, Vec3::new(-0.5, 0.0, 0.0), 0.0);
    vector(
        out.state.linear_momentum_world(),
        Vec3::new(15.0, -12.0, 3.0),
        0.0,
    );
    vector(
        parcel.linear_momentum_world_ns,
        Vec3::new(5.0, 4.0, 1.0),
        0.0,
    );
    vector(
        out.state.angular_momentum_body(),
        Vec3::new(0.0, 0.0, 2.0),
        2e-14,
    );
    vector(
        out.properties
            .angular_velocity_body_checked(out.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(0.0, 0.0, 0.5),
        5e-15,
    );
    close(b.kinetic_energy(p).unwrap(), 74.0, 1e-14);
    close(out.supplied_work_j, 10.5, 3e-14);
    vector(
        out.state
            .pose()
            .point_world_from_body(outgoing.attachment_arm_body_m.sub(out.center_shift_body_m))
            .unwrap(),
        parcel.position_world_m,
        0.0,
    );
    conservation(b, p, out);
}

#[test]
fn zero_impulse_release_inherits_rigid_motion_and_reverses_corotating_capture() {
    let retained = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let original = state(
        retained,
        Vec3::new(-0.5, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(5.0, -3.0, 1.0),
        Vec3::new(0.0, 0.0, 2.0),
    );
    let incoming = PointMassCapture {
        mass_kg: 1.0,
        attachment_arm_body_m: Vec3::new(2.0, 0.0, 0.0),
        velocity_world_m_s: Vec3::new(5.0, 1.0, 1.0),
    };
    let captured = capture_point_mass(&original, retained, incoming, controls(), || false).unwrap();
    assert_eq!(captured.removed_relative_energy_j, 0.0);
    let outgoing = PointMassRelease {
        mass_kg: incoming.mass_kg,
        attachment_arm_body_m: incoming
            .attachment_arm_body_m
            .sub(captured.center_shift_body_m),
        impulse_on_parcel_body_ns: Vec3::ZERO,
    };
    let mut no_impulses = controls();
    no_impulses.maximum_impulses = 0;
    let out = release_point_mass(
        &captured.state,
        captured.properties,
        outgoing,
        no_impulses,
        || false,
    )
    .unwrap();
    assert_eq!(out.supplied_work_j, 0.0);
    inertia(out.properties, retained, 3e-14);
    vector(
        out.state.pose().position_world(),
        original.pose().position_world(),
        2e-15,
    );
    vector(
        out.state.linear_momentum_world(),
        original.linear_momentum_world(),
        3e-14,
    );
    vector(
        out.state.angular_momentum_body(),
        original.angular_momentum_body(),
        3e-14,
    );
    vector(
        out.parcel.unwrap().velocity_world_m_s().unwrap(),
        incoming.velocity_world_m_s,
        1e-14,
    );
    vector(
        out.parcel.unwrap().position_world_m,
        Vec3::new(1.5, 0.0, 0.0),
        2e-15,
    );
    conservation(captured.state, captured.properties, out);
}

#[test]
fn full_tensor_release_is_covariant_and_recapture_removes_exact_separation_work() {
    let p = MassProperties::from_principal_frame(
        5.0,
        Vec3::ZERO,
        Vec3::new(3.0, 4.0, 5.0),
        rotation(Vec3::new(1.0, -0.4, 0.6), 0.7),
    )
    .unwrap();
    let b = state(
        p,
        Vec3::new(2.0, -1.0, 0.4),
        rotation(Vec3::new(-0.3, 0.8, 1.0), -0.4),
        Vec3::new(0.3, -0.8, 1.1),
        Vec3::new(0.5, -0.4, 0.7),
    );
    let outgoing = PointMassRelease {
        mass_kg: 0.5,
        attachment_arm_body_m: Vec3::new(0.6, -0.3, 0.2),
        impulse_on_parcel_body_ns: Vec3::new(0.7, -1.2, 0.4),
    };
    let out = release(&b, p, outgoing);
    assert!(out.properties.inertia_tensor_body().m[1].abs() > 0.05);
    conservation(b, p, out);
    let parcel = out.parcel.unwrap();
    let captured = capture_point_mass(
        &out.state,
        out.properties,
        PointMassCapture {
            mass_kg: parcel.mass_kg,
            attachment_arm_body_m: outgoing.attachment_arm_body_m.sub(out.center_shift_body_m),
            velocity_world_m_s: parcel.velocity_world_m_s().unwrap(),
        },
        controls(),
        || false,
    )
    .unwrap();
    inertia(captured.properties, p, 5e-13);
    vector(
        captured.state.pose().position_world(),
        b.pose().position_world(),
        3e-14,
    );
    vector(
        captured.state.linear_momentum_world(),
        b.linear_momentum_world(),
        3e-14,
    );
    vector(
        captured.state.angular_momentum_body(),
        b.angular_momentum_body(),
        5e-13,
    );
    close(
        captured.removed_relative_energy_j,
        out.supplied_work_j,
        4e-13,
    );

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
    let transformed = release(&moved, p, outgoing);
    vector(
        transformed.state.pose().position_world(),
        turn.rotate_body_to_world(out.state.pose().position_world())
            .add(offset),
        3e-14,
    );
    vector(
        transformed.parcel.unwrap().position_world_m,
        turn.rotate_body_to_world(parcel.position_world_m)
            .add(offset),
        3e-14,
    );
    vector(
        transformed.state.linear_momentum_world(),
        turn.rotate_body_to_world(out.state.linear_momentum_world())
            .add(boost.scale(out.properties.mass())),
        9e-14,
    );
    vector(
        transformed.parcel.unwrap().linear_momentum_world_ns,
        turn.rotate_body_to_world(parcel.linear_momentum_world_ns)
            .add(boost.scale(parcel.mass_kg)),
        2e-14,
    );
    assert_eq!(
        transformed.state.angular_momentum_body(),
        out.state.angular_momentum_body()
    );
    assert_eq!(transformed.properties, out.properties);
    assert_eq!(transformed.supplied_work_j, out.supplied_work_j);
    conservation(moved, p, transformed);
}

#[test]
fn separation_work_survives_large_translation_and_lost_geometry_refuses() {
    let (p, b, outgoing) = spinning_fixture();
    let plain = release(&b, p, outgoing);
    let boost = Vec3::new(1_099_511_627_776.0, 0.0, 0.0);
    let boosted = RigidBodyState::new(
        b.pose(),
        b.linear_momentum_world().add(boost.scale(p.mass())),
        b.angular_momentum_body(),
    )
    .unwrap();
    let out = release(&boosted, p, outgoing);
    assert_eq!(out.supplied_work_j, plain.supplied_work_j);
    assert_eq!(
        out.state.angular_momentum_body(),
        plain.state.angular_momentum_body()
    );
    assert_eq!(out.balance, plain.balance);
    // Laboratory kinetic totals cannot resolve this physical work.
    assert_eq!(
        boosted.kinetic_energy(p).unwrap(),
        boosted.kinetic_energy(p).unwrap() + out.supplied_work_j,
    );
    let distant = RigidBodyState::new(
        Pose::new(
            Vec3::new(36_028_797_018_963_968.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        b.linear_momentum_world(),
        b.angular_momentum_body(),
    )
    .unwrap();
    assert!(matches!(
        release_point_mass(&distant, p, outgoing, controls(), || false),
        Err(AssemblyError::Balance(_)),
    ));
}

#[test]
fn identities_mass_inertia_work_budgets_and_late_cancellation_are_atomic() {
    let (p, b, outgoing) = spinning_fixture();
    let expected = release(&b, p, outgoing);
    let mut one_part = controls();
    one_part.maximum_parts = 1;
    one_part.maximum_impulses = 0;
    let identity = release_point_mass(
        &b,
        p,
        PointMassRelease {
            mass_kg: 0.0,
            impulse_on_parcel_body_ns: Vec3::ZERO,
            ..outgoing
        },
        one_part,
        || false,
    )
    .unwrap();
    assert_eq!(identity.state, b);
    assert_eq!(identity.properties, p);
    assert_eq!(identity.center_shift_body_m, Vec3::ZERO);
    assert_eq!(identity.parcel, None);
    assert_eq!(identity.supplied_work_j, 0.0);
    assert!(matches!(
        release_point_mass(&b, p, outgoing, one_part, || false),
        Err(AssemblyError::Budget { .. }),
    ));
    for mass_kg in [
        0.0,
        -1.0,
        p.mass(),
        p.mass() + 1.0,
        f64::NAN,
        f64::EPSILON / 8.0,
    ] {
        assert!(
            release_point_mass(
                &b,
                p,
                PointMassRelease {
                    mass_kg,
                    ..outgoing
                },
                controls(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        release_point_mass(
            &b,
            p,
            PointMassRelease {
                impulse_on_parcel_body_ns: Vec3::new(f64::INFINITY, 0.0, 0.0),
                ..outgoing
            },
            controls(),
            || false,
        )
        .is_err()
    );
    assert!(
        release_point_mass(
            &b,
            p,
            PointMassRelease {
                attachment_arm_body_m: Vec3::new(100.0, 0.0, 0.0),
                ..outgoing
            },
            controls(),
            || false,
        )
        .is_err()
    );
    // The retained tensor would remain SPD but violate the rigid inertia
    // triangle inequality. It must not be clipped or silently diagonalized.
    let triangle = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.9, 1.0, 1.0)).unwrap();
    let resting = state(
        triangle,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    assert!(
        release_point_mass(
            &resting,
            triangle,
            PointMassRelease {
                mass_kg: 1.0,
                attachment_arm_body_m: Vec3::new(0.5, 0.0, 0.0),
                impulse_on_parcel_body_ns: Vec3::ZERO,
            },
            controls(),
            || false,
        )
        .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        release_point_mass(&b, p, outgoing, controls(), || {
            polls += 1;
            polls == 3
        }),
        Err(AssemblyError::Cancelled),
    );
    assert_eq!(polls, 3);
    assert_eq!(release(&b, p, outgoing), expected);
    assert_eq!(b, spinning_fixture().1);
    assert_eq!(p, spinning_fixture().0);
    assert!(
        ReleasedPointMass {
            mass_kg: 0.0,
            position_world_m: Vec3::ZERO,
            linear_momentum_world_ns: Vec3::ZERO,
        }
        .velocity_world_m_s()
        .is_err()
    );
}
