use fs_mbd::assembly::partition::AssemblyPartition;
use fs_mbd::assembly::{
    AssemblyControls, AssemblyError, AssemblyPart, RigidAssembly, SeparationImpulse,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn controls() -> AssemblyControls {
    AssemblyControls {
        maximum_parts: 16,
        maximum_impulses: 16,
        linear_momentum_tolerance_ns: 1e-10,
        angular_momentum_tolerance_nms: 1e-9,
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
fn rotation(axis: Vec3, angle: f64) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(axis, angle).unwrap()
}
fn parts() -> Vec<AssemblyPart> {
    (0..6)
        .map(|i| {
            let x = i as f64;
            AssemblyPart {
                properties: MassProperties::from_principal_frame(
                    1.0 + x,
                    Vec3::ZERO,
                    Vec3::new(0.4, 0.7, 0.9).scale(1.0 + x),
                    rotation(Vec3::new(1.0, 2.0, -0.5), 0.1 + 0.17 * x),
                )
                .unwrap(),
                pose_in_reference: Pose::new(
                    Vec3::new(0.8 * x - 1.5, (x - 2.0).powi(2) * 0.1, -0.2 * x),
                    rotation(Vec3::new(-0.3, 0.7, 1.0), -0.4 + 0.3 * x),
                )
                .unwrap(),
            }
        })
        .collect()
}
fn parent(assembly: &RigidAssembly) -> RigidBodyState {
    let m = assembly.properties();
    RigidBodyState::new(
        Pose::new(
            Vec3::new(3.0, -2.0, 1.0),
            rotation(Vec3::new(1.0, -2.0, 0.3), 0.6),
        )
        .unwrap(),
        Vec3::new(1.4, -0.6, 0.2).scale(m.mass()),
        m.angular_momentum_body_checked(Vec3::new(0.7, -0.8, 1.3))
            .unwrap(),
    )
    .unwrap()
}
fn impulse(a: usize, b: usize, scale: f64) -> SeparationImpulse {
    SeparationImpulse {
        part_a: a,
        part_b: b,
        point_assembly_m: Vec3::new(0.2, -0.3, 0.4).scale(scale),
        impulse_on_a_assembly_ns: Vec3::new(0.6, 0.3, -0.5).scale(scale),
        couple_on_a_assembly_nms: Vec3::new(0.1, -0.2, 0.3).scale(scale),
    }
}
// Independent world sums, including spin and orbital angular momentum about
// one arbitrary observer. No returned balance receipt enters this calculation.
fn totals(bodies: &[RigidBodyState], properties: &[MassProperties]) -> (Vec3, Vec3, f64) {
    let origin = Vec3::new(-4.0, 2.0, -1.0);
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut energy = 0.0;
    for (body, &mass) in bodies.iter().zip(properties) {
        linear = linear.add(body.linear_momentum_world());
        angular = angular
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(body.angular_momentum_body()),
            )
            .add(
                body.pose()
                    .position_world()
                    .sub(origin)
                    .cross(body.linear_momentum_world()),
            );
        energy += body.kinetic_energy(mass).unwrap();
    }
    (linear, angular, energy)
}
fn same_motion(a: RigidBodyState, b: RigidBodyState, mass: MassProperties, tolerance: f64) {
    vector(
        a.pose().position_world(),
        b.pose().position_world(),
        tolerance,
    );
    vector(
        a.linear_momentum_world(),
        b.linear_momentum_world(),
        tolerance,
    );
    vector(
        a.pose()
            .orientation()
            .rotate_body_to_world(a.angular_momentum_body()),
        b.pose()
            .orientation()
            .rotate_body_to_world(b.angular_momentum_body()),
        tolerance,
    );
    let point = Vec3::new(0.4, -0.3, 0.2);
    let x = a.point_kinematics(mass, point).unwrap();
    let y = b.point_kinematics(mass, point).unwrap();
    vector(x.point_world, y.point_world, tolerance);
    vector(x.point_velocity_world, y.point_velocity_world, tolerance);
}

#[test]
fn interleaved_full_tensor_groups_preserve_original_material_motion() {
    let assembly = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let partition =
        AssemblyPartition::new(&assembly, &[0, 1, 0, 2, 1, 2], controls(), || false).unwrap();
    assert_eq!(partition.original_part_indices(0), Some([0, 2].as_slice()));
    assert_eq!(partition.original_part_indices(1), Some([1, 4].as_slice()));
    assert_eq!(partition.original_part_indices(2), Some([3, 5].as_slice()));
    assert_eq!(partition.original_part_indices(3), None);
    assert_eq!(partition.original_part_indices(usize::MAX), None);
    let p = parent(&assembly);
    let split = partition.release(p, &[], controls(), || false).unwrap();
    let expanded = partition
        .expand(&split.bodies, controls(), || false)
        .unwrap();
    let flat = assembly.release(p, &[], controls(), || false).unwrap();
    for ((a, b), part) in expanded.iter().zip(&flat.bodies).zip(assembly.parts()) {
        same_motion(*a, *b, part.properties, 2e-12);
    }
    let masses: Vec<_> = assembly.parts().iter().map(|x| x.properties).collect();
    let before = totals(&[p], &[assembly.properties()]);
    let after = totals(&expanded, &masses);
    vector(before.0, after.0, 2e-12);
    vector(before.1, after.1, 2e-12);
    close(before.2, after.2, 2e-12);
    assert_eq!(split.supplied_work_j, 0.0);
}

#[test]
fn off_center_separation_recoils_whole_stages_with_closed_form_work() {
    let mass = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.4, 0.4, 0.4)).unwrap();
    let parts: Vec<_> = [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
        .into_iter()
        .map(|(x, y)| AssemblyPart {
            properties: mass,
            pose_in_reference: Pose::new(Vec3::new(x, y, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        })
        .collect();
    let assembly = RigidAssembly::new(&parts, controls(), || false).unwrap();
    let partition = AssemblyPartition::new(&assembly, &[0, 0, 1, 1], controls(), || false).unwrap();
    let state = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap();
    let push = SeparationImpulse {
        part_a: 0,
        part_b: 2,
        point_assembly_m: Vec3::ZERO,
        impulse_on_a_assembly_ns: Vec3::new(0.0, 2.0, 0.0),
        couple_on_a_assembly_nms: Vec3::ZERO,
    };
    let split = partition
        .release(state, &[push], controls(), || false)
        .unwrap();
    close(split.supplied_work_j, 11.0 / 6.0, 2e-15);
    for (i, body) in split.bodies.iter().enumerate() {
        let properties = partition.groups()[i].properties();
        vector(body.angular_momentum_body(), Vec3::new(0.0, 0.0, 2.0), 0.0);
        vector(
            body.center_of_mass_velocity_world(properties).unwrap(),
            Vec3::new(0.0, if i == 0 { 0.5 } else { -0.5 }, 0.0),
            0.0,
        );
        close(
            properties
                .angular_velocity_body_checked(body.angular_momentum_body())
                .unwrap()
                .z,
            5.0 / 12.0,
            1e-15,
        );
    }
    let expanded = partition
        .expand(&split.bodies, controls(), || false)
        .unwrap();
    let result = totals(&expanded, &[mass; 4]);
    vector(result.0, Vec3::ZERO, 1e-15);
    vector(result.1, Vec3::ZERO, 3e-15);
    close(result.2, 11.0 / 6.0, 3e-15);
    // Detaching each leaf would allocate the kick to a single small inertia,
    // producing a different physical motion and spuriously larger work.
    let wrong_topology = assembly
        .release(state, &[push], controls(), || false)
        .unwrap();
    assert!(wrong_topology.supplied_work_j > split.supplied_work_j + 1.0);
}

#[test]
fn shared_group_impulses_match_native_full_inertia_updates_and_cross_work() {
    let assembly = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let membership = [0, 1, 0, 2, 1, 2];
    let partition = AssemblyPartition::new(&assembly, &membership, controls(), || false).unwrap();
    let p = parent(&assembly);
    let pushes = [impulse(0, 1, 1.0), impulse(2, 3, 0.7), impulse(4, 5, -0.3)];
    let base = partition.release(p, &[], controls(), || false).unwrap();
    let actual = partition.release(p, &pushes, controls(), || false).unwrap();
    let mut reference = base.bodies.clone();
    for push in &pushes {
        let point = p
            .pose()
            .point_world_from_body(push.point_assembly_m)
            .unwrap();
        for (part, sign) in [(push.part_a, 1.0), (push.part_b, -1.0)] {
            let group = membership[part];
            let old = reference[group];
            let arm = old
                .pose()
                .orientation()
                .rotate_world_to_body(point.sub(old.pose().position_world()));
            let j = p
                .pose()
                .orientation()
                .rotate_body_to_world(push.impulse_on_a_assembly_ns.scale(sign));
            let kick = old
                .apply_impulse_at_body_point(partition.groups()[group].properties(), arm, j)
                .unwrap()
                .state_after;
            let couple = p
                .pose()
                .orientation()
                .rotate_body_to_world(push.couple_on_a_assembly_nms.scale(sign));
            reference[group] = RigidBodyState::new(
                kick.pose(),
                kick.linear_momentum_world(),
                kick.angular_momentum_body()
                    .add(kick.pose().orientation().rotate_world_to_body(couple)),
            )
            .unwrap();
        }
    }
    let masses: Vec<_> = partition
        .groups()
        .iter()
        .map(RigidAssembly::properties)
        .collect();
    for ((a, b), &mass) in actual.bodies.iter().zip(&reference).zip(&masses) {
        same_motion(*a, *b, mass, 2e-12);
    }
    let before = totals(&base.bodies, &masses);
    let after = totals(&reference, &masses);
    close(after.2 - before.2, actual.supplied_work_j, 2e-12);
    vector(after.0, before.0, 2e-12);
    vector(after.1, before.1, 2e-12);
    let separate_pair_work: f64 = pushes
        .iter()
        .map(|push| {
            partition
                .release(p, &[*push], controls(), || false)
                .unwrap()
                .supplied_work_j
        })
        .sum();
    assert!((separate_pair_work - actual.supplied_work_j).abs() > 1e-3);
}

#[test]
fn independently_flying_stages_expand_and_separate_again_with_original_geometry() {
    let assembly = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let partition =
        AssemblyPartition::new(&assembly, &[0, 1, 0, 1, 0, 1], controls(), || false).unwrap();
    let split = partition
        .release(parent(&assembly), &[impulse(0, 1, 2.0)], controls(), || {
            false
        })
        .unwrap();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut states = split.bodies;
    for _ in 0..20 {
        for (state, group) in states.iter_mut().zip(partition.groups()) {
            *state = integrator
                .step(*state, group.properties(), Wrench::ZERO, 0.01)
                .unwrap()
                .state_after;
        }
    }
    let expanded = partition.expand(&states, controls(), || false).unwrap();
    for (group_id, (group, &state)) in partition.groups().iter().zip(&states).enumerate() {
        let expected = group.release(state, &[], controls(), || false).unwrap();
        for ((&original, local), part) in partition
            .original_part_indices(group_id)
            .unwrap()
            .iter()
            .zip(&expected.bodies)
            .zip(group.parts())
        {
            same_motion(expanded[original], *local, part.properties, 2e-12);
        }
    }
    let child = &partition.groups()[0];
    let nested = AssemblyPartition::new(child, &[0, 1, 1], controls(), || false).unwrap();
    let next = nested
        .release(states[0], &[impulse(0, 1, 0.8)], controls(), || false)
        .unwrap();
    let leaves = nested.expand(&next.bodies, controls(), || false).unwrap();
    let masses: Vec<_> = child.parts().iter().map(|p| p.properties).collect();
    let before = totals(&[states[0]], &[child.properties()]);
    let after = totals(&leaves, &masses);
    vector(before.0, after.0, 2e-12);
    vector(before.1, after.1, 2e-12);
    close(after.2 - before.2, next.supplied_work_j, 2e-12);
}

#[test]
fn single_group_and_singleton_groups_preserve_material_motion() {
    let assembly = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    let p = parent(&assembly);
    let expected = assembly.release(p, &[], controls(), || false).unwrap();
    for membership in [[0; 6], [0, 1, 2, 3, 4, 5]] {
        let partition =
            AssemblyPartition::new(&assembly, &membership, controls(), || false).unwrap();
        let groups = partition.release(p, &[], controls(), || false).unwrap();
        let bodies = partition
            .expand(&groups.bodies, controls(), || false)
            .unwrap();
        for ((a, b), part) in bodies.iter().zip(&expected.bodies).zip(assembly.parts()) {
            same_motion(*a, *b, part.properties, 3e-12);
        }
    }
}

#[test]
fn malformed_partitions_caps_and_late_cancel_refuse_without_changing_inputs() {
    let assembly = RigidAssembly::new(&parts(), controls(), || false).unwrap();
    for ids in [
        vec![],
        vec![0; 5],
        vec![1; 6],
        vec![0, 2, 0, 2, 0, 2],
        vec![usize::MAX; 6],
    ] {
        assert!(AssemblyPartition::new(&assembly, &ids, controls(), || false).is_err());
    }
    let membership = [0, 1, 0, 1, 0, 1];
    let partition = AssemblyPartition::new(&assembly, &membership, controls(), || false).unwrap();
    let saved = partition.clone();
    let p = parent(&assembly);
    assert!(
        partition
            .release(p, &[impulse(0, 2, 1.0)], controls(), || false)
            .is_err()
    );
    assert!(
        partition
            .release(p, &[impulse(0, 6, 1.0)], controls(), || false)
            .is_err()
    );
    assert!(
        partition
            .release(p, &[impulse(0, 1, f64::NAN)], controls(), || false)
            .is_err()
    );
    let small = AssemblyControls {
        maximum_parts: 5,
        ..controls()
    };
    assert!(matches!(
        partition.release(p, &[], small, || false),
        Err(AssemblyError::Budget { .. })
    ));
    let small = AssemblyControls {
        maximum_impulses: 0,
        ..controls()
    };
    assert!(matches!(
        partition.release(p, &[impulse(0, 1, 1.0)], small, || false),
        Err(AssemblyError::Budget { .. })
    ));
    let split = partition
        .release(p, &[impulse(0, 1, 1.0)], controls(), || false)
        .unwrap();
    let states = split.bodies.clone();
    let mut polls = 0;
    let expected = partition
        .expand(&states, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    let refused = partition.expand(&states, controls(), || {
        calls += 1;
        calls == polls - 1
    });
    assert_eq!(refused, Err(AssemblyError::Cancelled));
    assert_eq!(partition, saved);
    assert_eq!(states, split.bodies);
    assert_eq!(
        partition.expand(&states, controls(), || false).unwrap(),
        expected
    );
    assert!(
        partition
            .expand(&states[..1], controls(), || false)
            .is_err()
    );
}
