use fs_mbd::assembly::{
    AssemblyControls, AssemblyError, AssemblyPart, RigidAssembly, SeparationImpulse,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls() -> AssemblyControls {
    AssemblyControls {
        maximum_parts: 8,
        maximum_impulses: 16,
        linear_momentum_tolerance_ns: 1e-10,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 1e-9,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; error={:.3e}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
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
fn parts() -> Vec<AssemblyPart> {
    let masses = [2.0, 3.0, 5.0];
    let centers = [
        Vec3::new(-0.8, 0.4, 0.7),
        Vec3::new(0.5, -0.2, 0.1),
        Vec3::new(0.6, 0.1, -0.3),
    ];
    (0..3)
        .map(|i| AssemblyPart {
            properties: MassProperties::from_principal_frame(
                masses[i],
                Vec3::ZERO,
                Vec3::new(0.4, 0.7, 0.9).scale(masses[i]),
                rotation(Vec3::new(0.4, 1.3, -0.6), 0.2 + 0.17 * i as f64),
            )
            .unwrap(),
            pose_in_reference: Pose::new(
                centers[i],
                rotation(Vec3::new(1.1, -0.5, 0.7), -0.4 + 0.3 * i as f64),
            )
            .unwrap(),
        })
        .collect()
}
fn state(
    p: MassProperties,
    position: Vec3,
    q: UnitQuaternion,
    velocity: Vec3,
    omega: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, q).unwrap(),
        velocity.scale(p.mass()),
        p.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap()
}
fn parent(a: &RigidAssembly) -> RigidBodyState {
    state(
        a.properties(),
        Vec3::new(3.0, -2.0, 1.0),
        rotation(Vec3::new(-1.0, 2.0, 0.3), 0.6),
        Vec3::new(1.4, -0.6, 0.2),
        Vec3::new(0.7, -0.8, 1.3),
    )
}
// Direct world sums about a fixed observer, independent of the assembly's
// translating-frame reduction and returned balance values.
fn totals(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    origin: Vec3,
) -> (Vec3, Vec3, f64) {
    let mut p = Vec3::ZERO;
    let mut l = Vec3::ZERO;
    let mut e = 0.0;
    for (b, m) in bodies.iter().zip(properties) {
        p = p.add(b.linear_momentum_world());
        l = l
            .add(
                b.pose()
                    .orientation()
                    .rotate_body_to_world(b.angular_momentum_body()),
            )
            .add(
                b.pose()
                    .position_world()
                    .sub(origin)
                    .cross(b.linear_momentum_world()),
            );
        e += b.kinetic_energy(*m).unwrap();
    }
    (p, l, e)
}

#[test]
fn composite_inertia_matches_independent_dense_parallel_axis_assembly() {
    let parts = parts();
    let a = RigidAssembly::new(&parts, controls(), || false).unwrap();
    let total: f64 = parts.iter().map(|p| p.properties.mass()).sum();
    let center = parts.iter().fold(Vec3::ZERO, |v, p| {
        v.add(
            p.pose_in_reference
                .position_world()
                .scale(p.properties.mass() / total),
        )
    });
    close(a.properties().mass(), total, 0.0);
    vector(a.center_of_mass_in_reference(), center, 2e-16);
    let basis = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let mut reference = [0.0; 9];
    for part in &parts {
        let r = part.pose_in_reference.position_world().sub(center);
        let q = part.pose_in_reference.orientation();
        for j in 0..3 {
            let column = q
                .rotate_body_to_world(
                    part.properties
                        .angular_momentum_body_checked(q.rotate_world_to_body(basis[j]))
                        .unwrap(),
                )
                .add(basis[j].scale(part.properties.mass() * r.dot(r)))
                .sub(r.scale(part.properties.mass() * r.dot(basis[j])));
            for i in 0..3 {
                reference[3 * i + j] += basis[i].dot(column);
            }
        }
    }
    assert!(reference[1].abs() + reference[2].abs() + reference[5].abs() > 0.1);
    for (actual, expected) in a
        .properties()
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(reference)
    {
        close(actual, expected, 2e-13);
    }
    let mut shifted = parts.clone();
    for p in &mut shifted {
        p.pose_in_reference = Pose::new(
            p.pose_in_reference
                .position_world()
                .add(Vec3::new(1024.0, -2048.0, 512.0)),
            p.pose_in_reference.orientation(),
        )
        .unwrap();
    }
    let shifted = RigidAssembly::new(&shifted, controls(), || false).unwrap();
    for (x, y) in a
        .properties()
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(shifted.properties().inertia_tensor_body().m)
    {
        close(x, y, 2e-12);
    }
}

#[test]
fn release_retains_every_component_pose_velocity_spin_and_orbital_momentum() {
    let a = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let p = parent(&a);
    let released = a.release(p, &[], controls(), || false).unwrap();
    let properties: Vec<_> = a.parts().iter().map(|p| p.properties).collect();
    for (part, body) in a.parts().iter().zip(&released.bodies) {
        let native = p
            .point_kinematics(a.properties(), part.pose_in_reference.position_world())
            .unwrap();
        vector(body.pose().position_world(), native.point_world, 2e-15);
        vector(
            body.center_of_mass_velocity_world(part.properties).unwrap(),
            native.point_velocity_world,
            2e-15,
        );
        vector(
            body.pose().orientation().rotate_body_to_world(
                part.properties
                    .angular_velocity_body_checked(body.angular_momentum_body())
                    .unwrap(),
            ),
            native.angular_velocity_world,
            3e-15,
        );
    }
    let origin = Vec3::new(-4.0, 2.0, 0.5);
    let before = totals(&[p], &[a.properties()], origin);
    let after = totals(&released.bodies, &properties, origin);
    vector(after.0, before.0, 2e-14);
    vector(after.1, before.1, 7e-14);
    close(after.2, before.2, 7e-14);
    assert_eq!(released.supplied_work_j, 0.0);
}

#[test]
fn shared_body_impulses_include_cross_work_and_match_native_sequential_application() {
    let a = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let p = parent(&a);
    let impulses = [
        SeparationImpulse {
            part_a: 0,
            part_b: 1,
            point_assembly_m: Vec3::new(0.1, 0.3, -0.2),
            impulse_on_a_assembly_ns: Vec3::new(0.3, -0.2, 0.5),
            couple_on_a_assembly_nms: Vec3::new(0.1, 0.04, -0.03),
        },
        SeparationImpulse {
            part_a: 0,
            part_b: 2,
            point_assembly_m: Vec3::new(-0.3, 0.2, 0.4),
            impulse_on_a_assembly_ns: Vec3::new(0.2, 0.6, -0.1),
            couple_on_a_assembly_nms: Vec3::new(-0.05, 0.07, 0.02),
        },
    ];
    let unforced = a.release(p, &[], controls(), || false).unwrap();
    let mut expected = unforced.bodies.clone();
    for input in impulses {
        let point = p
            .pose()
            .point_world_from_body(input.point_assembly_m)
            .unwrap();
        let j = p
            .pose()
            .orientation()
            .rotate_body_to_world(input.impulse_on_a_assembly_ns);
        let c = p
            .pose()
            .orientation()
            .rotate_body_to_world(input.couple_on_a_assembly_nms);
        for (i, sign) in [(input.part_a, 1.0), (input.part_b, -1.0)] {
            let arm = expected[i].pose().point_body_from_world(point).unwrap();
            let b = expected[i]
                .apply_impulse_at_body_point(a.parts()[i].properties, arm, j.scale(sign))
                .unwrap()
                .state_after;
            expected[i] = RigidBodyState::new(
                b.pose(),
                b.linear_momentum_world(),
                b.angular_momentum_body()
                    .add(b.pose().orientation().rotate_world_to_body(c.scale(sign))),
            )
            .unwrap();
        }
    }
    let result = a.release(p, &impulses, controls(), || false).unwrap();
    for (x, y) in result.bodies.iter().zip(&expected) {
        vector(x.linear_momentum_world(), y.linear_momentum_world(), 2e-15);
        vector(x.angular_momentum_body(), y.angular_momentum_body(), 2e-15);
    }
    let properties: Vec<_> = a.parts().iter().map(|p| p.properties).collect();
    let before = totals(&[p], &[a.properties()], Vec3::ZERO);
    let after = totals(&result.bodies, &properties, Vec3::ZERO);
    vector(after.0, before.0, 2e-14);
    vector(after.1, before.1, 3e-14);
    close(after.2 - before.2, result.supplied_work_j, 8e-14);
    let independent_pair_sum = a
        .release(p, &impulses[..1], controls(), || false)
        .unwrap()
        .supplied_work_j
        + a.release(p, &impulses[1..], controls(), || false)
            .unwrap()
            .supplied_work_j;
    assert!((independent_pair_sum - result.supplied_work_j).abs() > 0.01);
    let swapped = a
        .release(p, &[impulses[1], impulses[0]], controls(), || false)
        .unwrap();
    assert_eq!(result, swapped);
}

#[test]
fn perfectly_inelastic_lock_matches_closed_form_orbital_to_spin_transfer() {
    let mass = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 2.0)).unwrap();
    let bodies = [
        state(
            mass,
            Vec3::new(-2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(0.0, -4.0, 0.0),
            Vec3::ZERO,
        ),
        state(
            mass,
            Vec3::new(2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(0.0, 4.0, 0.0),
            Vec3::ZERO,
        ),
    ];
    let result = RigidAssembly::latch(
        &bodies,
        &[mass; 2],
        UnitQuaternion::IDENTITY,
        controls(),
        || false,
    )
    .unwrap();
    vector(result.state.linear_momentum_world(), Vec3::ZERO, 0.0);
    vector(
        result.state.angular_momentum_body(),
        Vec3::new(0.0, 0.0, 48.0),
        0.0,
    );
    close(
        result
            .assembly
            .properties()
            .angular_velocity_body_checked(result.state.angular_momentum_body())
            .unwrap()
            .z,
        12.0 / 7.0,
        3e-16,
    );
    close(result.removed_relative_energy_j, 48.0 / 7.0, 3e-14);
    let locked = result
        .assembly
        .release(result.state, &[], controls(), || false)
        .unwrap();
    for (i, b) in locked.bodies.iter().enumerate() {
        vector(
            b.pose().position_world(),
            bodies[i].pose().position_world(),
            0.0,
        );
        close(
            b.center_of_mass_velocity_world(mass).unwrap().y,
            if i == 0 { -24.0 / 7.0 } else { 24.0 / 7.0 },
            1e-15,
        );
        close(
            mass.angular_velocity_body_checked(b.angular_momentum_body())
                .unwrap()
                .z,
            12.0 / 7.0,
            3e-16,
        );
    }
    let again = RigidAssembly::latch(
        &locked.bodies,
        &[mass; 2],
        UnitQuaternion::IDENTITY,
        controls(),
        || false,
    )
    .unwrap();
    close(again.removed_relative_energy_j, 0.0, 1e-27);
}

#[test]
fn latch_and_release_are_world_frame_covariant_and_keep_arbitrary_component_frames() {
    let a = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let p = parent(&a);
    let mut bodies = a.release(p, &[], controls(), || false).unwrap().bodies;
    let properties: Vec<_> = a.parts().iter().map(|p| p.properties).collect();
    bodies[0] = bodies[0]
        .apply_impulse_at_body_point(
            properties[0],
            Vec3::new(0.2, -0.1, 0.3),
            Vec3::new(0.4, 0.5, -0.2),
        )
        .unwrap()
        .state_after;
    let q = rotation(Vec3::new(0.3, -0.7, 1.2), 0.8);
    let shift = Vec3::new(-7.0, 4.0, 2.0);
    let boost = Vec3::new(8.0, -3.0, 5.0);
    let frame = rotation(Vec3::new(0.2, 0.8, 1.3), 0.4);
    let first = RigidAssembly::latch(&bodies, &properties, frame, controls(), || false).unwrap();
    let transformed: Vec<_> = bodies
        .iter()
        .zip(&properties)
        .map(|(b, m)| {
            RigidBodyState::new(
                Pose::new(
                    q.rotate_body_to_world(b.pose().position_world()).add(shift),
                    compose(q, b.pose().orientation()),
                )
                .unwrap(),
                q.rotate_body_to_world(b.linear_momentum_world())
                    .add(boost.scale(m.mass())),
                b.angular_momentum_body(),
            )
            .unwrap()
        })
        .collect();
    let second = RigidAssembly::latch(
        &transformed,
        &properties,
        compose(q, frame),
        controls(),
        || false,
    )
    .unwrap();
    close(
        second.removed_relative_energy_j,
        first.removed_relative_energy_j,
        2e-14,
    );
    vector(
        second.state.linear_momentum_world(),
        q.rotate_body_to_world(first.state.linear_momentum_world())
            .add(boost.scale(a.properties().mass())),
        4e-14,
    );
    vector(
        second.state.angular_momentum_body(),
        first.state.angular_momentum_body(),
        2e-13,
    );
    let split = second
        .assembly
        .release(second.state, &[], controls(), || false)
        .unwrap();
    for (b, original) in split.bodies.iter().zip(&transformed) {
        vector(
            b.pose().position_world(),
            original.pose().position_world(),
            4e-15,
        );
        vector(
            b.pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(0.2, 0.3, 0.5)),
            original
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(0.2, 0.3, 0.5)),
            8e-16,
        );
    }
}

#[test]
fn small_latch_loss_survives_a_large_common_translational_energy() {
    let m = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let bodies = [
        state(
            m,
            Vec3::new(-1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(1e9 + 1.0, 0.0, 0.0),
            Vec3::ZERO,
        ),
        state(
            m,
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(1e9 - 1.0, 0.0, 0.0),
            Vec3::ZERO,
        ),
    ];
    let result = RigidAssembly::latch(
        &bodies,
        &[m; 2],
        UnitQuaternion::IDENTITY,
        controls(),
        || false,
    )
    .unwrap();
    close(result.removed_relative_energy_j, 2.0, 0.0);
    close(result.balance.energy_residual_j, 0.0, 0.0);
    let original_energy = bodies
        .iter()
        .map(|b| b.kinetic_energy(m).unwrap())
        .sum::<f64>();
    assert_eq!(
        original_energy,
        result
            .state
            .kinetic_energy(result.assembly.properties())
            .unwrap()
    );
}

#[test]
fn late_cancellation_caps_and_invalid_impulse_pairs_never_publish_partial_motion() {
    let components = parts();
    let a = RigidAssembly::new(&components, controls(), || false).unwrap();
    let p = parent(&a);
    assert!(matches!(
        RigidAssembly::new(&[], controls(), || false),
        Err(AssemblyError::Invalid(_))
    ));
    assert!(matches!(
        RigidAssembly::new(
            &components,
            AssemblyControls {
                maximum_parts: 2,
                ..controls()
            },
            || false
        ),
        Err(AssemblyError::Budget { .. })
    ));
    let bad = SeparationImpulse {
        part_a: 0,
        part_b: 0,
        point_assembly_m: Vec3::ZERO,
        impulse_on_a_assembly_ns: Vec3::ZERO,
        couple_on_a_assembly_nms: Vec3::ZERO,
    };
    assert!(matches!(
        a.release(p, &[bad], controls(), || false),
        Err(AssemblyError::Invalid(_))
    ));
    let mut polls = 0;
    let good = a
        .release(p, &[], controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let total = polls;
    polls = 0;
    assert!(matches!(
        a.release(p, &[], controls(), || {
            polls += 1;
            polls == total
        }),
        Err(AssemblyError::Cancelled)
    ));
    assert_eq!(a.release(p, &[], controls(), || false).unwrap(), good);
    let props: Vec<_> = a.parts().iter().map(|p| p.properties).collect();
    let saved = good.bodies.clone();
    let mut polls = 0;
    RigidAssembly::latch(&saved, &props, UnitQuaternion::IDENTITY, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let total = polls;
    polls = 0;
    assert!(matches!(
        RigidAssembly::latch(&saved, &props, UnitQuaternion::IDENTITY, controls(), || {
            polls += 1;
            polls == total
        }),
        Err(AssemblyError::Cancelled)
    ));
    assert_eq!(saved, good.bodies);
    assert!(matches!(
        a.release(
            p,
            &[],
            AssemblyControls {
                linear_momentum_tolerance_ns: f64::NAN,
                ..controls()
            },
            || false
        ),
        Err(AssemblyError::Invalid(_))
    ));
}
