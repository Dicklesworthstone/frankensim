//! G0 independent impulse/work oracles, G3 frame/reference invariance, and
//! G4 source, actual-output and cancellation admission on native rigid states.

use fs_mbd::assembly::ejection::{
    RelativeSeparationControls, RelativeSeparationError, RelativeSeparationResult,
    RelativeSeparationTwist,
};
use fs_mbd::assembly::partition::AssemblyPartition;
use fs_mbd::assembly::{AssemblyControls, AssemblyError, AssemblyPart, RigidAssembly};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls() -> RelativeSeparationControls {
    RelativeSeparationControls {
        assembly: AssemblyControls {
            maximum_parts: 8,
            maximum_impulses: 8,
            linear_momentum_tolerance_ns: 1e-9,
            angular_momentum_tolerance_nms: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        relative_velocity_tolerance_m_s: 1e-9,
        relative_angular_velocity_tolerance_rad_s: 1e-9,
        maximum_supplied_work_j: 1e4,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    close(a.x, b.x, tol);
    close(a.y, b.y, tol);
    close(a.z, b.z, tol);
}
fn twist(v: Vec3, w: Vec3) -> RelativeSeparationTwist {
    RelativeSeparationTwist {
        relative_velocity_assembly_m_s: v,
        relative_angular_velocity_assembly_rad_s: w,
    }
}
fn sphere(m: f64, i: f64, r: Vec3) -> AssemblyPart {
    AssemblyPart {
        properties: MassProperties::new(m, Vec3::ZERO, Vec3::new(i, i, i)).unwrap(),
        pose_in_reference: Pose::new(r, UnitQuaternion::IDENTITY).unwrap(),
    }
}
fn assembly(parts: &[AssemblyPart]) -> RigidAssembly {
    RigidAssembly::new(parts, controls().assembly, || false).unwrap()
}
fn state(a: &RigidAssembly, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(a.properties().mass()),
        a.properties().angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn rest(a: &RigidAssembly) -> RigidBodyState {
    state(
        a,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    )
}
fn omega(b: RigidBodyState, m: MassProperties) -> Vec3 {
    b.pose().orientation().rotate_body_to_world(
        m.angular_velocity_body_checked(b.angular_momentum_body())
            .unwrap(),
    )
}
// Direct world sums about an independent observer, not returned balance values.
fn totals(states: &[RigidBodyState], masses: &[MassProperties]) -> (Vec3, Vec3, f64) {
    let origin = Vec3::new(-0.7, 1.3, -0.9);
    let (mut p, mut h, mut k) = (Vec3::ZERO, Vec3::ZERO, 0.0);
    for (&b, &m) in states.iter().zip(masses) {
        p = p.add(b.linear_momentum_world());
        h = h
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
        k += b.kinetic_energy(m).unwrap();
    }
    (p, h, k)
}
fn conservation(a: &RigidAssembly, parent: RigidBodyState, result: &RelativeSeparationResult) {
    let props: Vec<_> = a.parts().iter().map(|p| p.properties).collect();
    let before = totals(&[parent], &[a.properties()]);
    let after = totals(&result.release.bodies, &props);
    vector(after.0, before.0, 2e-8);
    vector(after.1, before.1, 2e-8);
    close(after.2 - before.2, result.release.supplied_work_j, 3e-8);
}
fn requested_motion(
    a: &RigidAssembly,
    parent: RigidBodyState,
    reference: usize,
    targets: &[RelativeSeparationTwist],
    result: &RelativeSeparationResult,
) {
    let b = result.release.bodies[reference];
    let v = b
        .center_of_mass_velocity_world(a.parts()[reference].properties)
        .unwrap();
    let w = omega(b, a.parts()[reference].properties);
    let q = parent.pose().orientation();
    for ((part, body), target) in a.parts().iter().zip(&result.release.bodies).zip(targets) {
        let point_velocity =
            v.add(w.cross(body.pose().position_world().sub(b.pose().position_world())));
        vector(
            body.center_of_mass_velocity_world(part.properties).unwrap(),
            point_velocity.add(q.rotate_body_to_world(target.relative_velocity_assembly_m_s)),
            2e-9,
        );
        vector(
            omega(*body, part.properties),
            w.add(q.rotate_body_to_world(target.relative_angular_velocity_assembly_rad_s)),
            2e-9,
        );
    }
}
fn two() -> (RigidAssembly, [RelativeSeparationTwist; 2]) {
    (
        assembly(&[sphere(2.0, 4.0, Vec3::ZERO), sphere(3.0, 6.0, Vec3::ZERO)]),
        [
            RelativeSeparationTwist::ZERO,
            twist(Vec3::new(5.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 2.0)),
        ],
    )
}

#[test]
fn g0_colocated_finite_components_recover_both_reduced_inertias() {
    let (a, targets) = two();
    let parent = rest(&a);
    let event = a
        .release_relative_twists(parent, 0, &targets, controls(), || false)
        .unwrap();
    // Linear reduced mass = 6/5; angular reduced inertia = 12/5.
    // J = (6/5)*5 = 6; K = (12/5)*2 = 24/5.
    vector(
        event.release.bodies[0].linear_momentum_world(),
        Vec3::new(-6.0, 0.0, 0.0),
        2e-14,
    );
    vector(
        event.release.bodies[1].linear_momentum_world(),
        Vec3::new(6.0, 0.0, 0.0),
        2e-14,
    );
    vector(
        omega(event.release.bodies[0], a.parts()[0].properties),
        Vec3::new(0.0, 0.0, -1.2),
        2e-14,
    );
    vector(
        omega(event.release.bodies[1], a.parts()[1].properties),
        Vec3::new(0.0, 0.0, 0.8),
        2e-14,
    );
    close(event.release.supplied_work_j, 19.8, 2e-13);
    assert_eq!(event.impulses.len(), 1);
    vector(
        event.impulses[0].couple_on_a_assembly_nms,
        Vec3::new(0.0, 0.0, 4.8),
        2e-14,
    );
    conservation(&a, parent, &event);
    requested_motion(&a, parent, 0, &targets, &event);
    // The solved equivalent wrench is accepted by the unchanged native owner.
    assert_eq!(
        event.release,
        a.release(parent, &event.impulses, controls().assembly, || false)
            .unwrap()
    );
}

#[test]
fn g0_transverse_separation_shares_spin_and_orbital_recoil() {
    let a = assembly(&[
        sphere(1.0, 2.0, Vec3::new(-1.0, 0.0, 0.0)),
        sphere(1.0, 2.0, Vec3::new(1.0, 0.0, 0.0)),
    ]);
    let target = [
        RelativeSeparationTwist::ZERO,
        twist(Vec3::new(0.0, 2.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    ];
    let parent = rest(&a);
    let event = a
        .release_relative_twists(parent, 0, &target, controls(), || false)
        .unwrap();
    // I_z=6, requested angular content=4, so common angular recoil=-2/3.
    // The orbital response reduces the linear impulse to 1/3, not the
    // independent translating-pair value 1.
    vector(
        event.release.bodies[0].linear_momentum_world(),
        Vec3::new(0.0, -1.0 / 3.0, 0.0),
        2e-14,
    );
    vector(
        event.release.bodies[1].linear_momentum_world(),
        Vec3::new(0.0, 1.0 / 3.0, 0.0),
        2e-14,
    );
    vector(
        event.release.bodies[0].angular_momentum_body(),
        Vec3::new(0.0, 0.0, -4.0 / 3.0),
        2e-14,
    );
    vector(
        event.release.bodies[1].angular_momentum_body(),
        Vec3::new(0.0, 0.0, 2.0 / 3.0),
        2e-14,
    );
    close(event.release.supplied_work_j, 2.0 / 3.0, 2e-14);
    conservation(&a, parent, &event);
    requested_motion(&a, parent, 0, &target, &event);
}

#[test]
fn g3_full_tensor_frames_reference_choice_and_component_order_preserve_motion() {
    let raw: Vec<_> = (0..3)
        .map(|i| {
            let z = i as f64;
            AssemblyPart {
                properties: MassProperties::from_principal_frame(
                    2.0 + z,
                    Vec3::ZERO,
                    Vec3::new(0.4, 0.7, 0.9).scale(2.0 + z),
                    UnitQuaternion::from_axis_angle(Vec3::new(0.4, 1.3, -0.6), 0.2 + 0.17 * z)
                        .unwrap(),
                )
                .unwrap(),
                pose_in_reference: Pose::new(
                    Vec3::new(0.8 * z - 1.0, 0.2 * z * z, -0.3 * z),
                    UnitQuaternion::from_axis_angle(Vec3::new(1.1, -0.5, 0.7), -0.4 + 0.3 * z)
                        .unwrap(),
                )
                .unwrap(),
            }
        })
        .collect();
    let a = assembly(&raw);
    let targets = [
        RelativeSeparationTwist::ZERO,
        twist(Vec3::new(1.2, -0.4, 0.7), Vec3::new(-0.2, 0.5, 0.3)),
        twist(Vec3::new(-0.8, 0.6, -0.1), Vec3::new(0.4, -0.1, 0.2)),
    ];
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.3), 0.6).unwrap();
    let parent = state(
        &a,
        Vec3::new(3.0, -2.0, 1.0),
        q,
        Vec3::new(1.4, -0.6, 0.2),
        Vec3::new(0.7, -0.8, 1.3),
    );
    let original = a
        .release_relative_twists(parent, 0, &targets, controls(), || false)
        .unwrap();
    conservation(&a, parent, &original);
    requested_motion(&a, parent, 0, &targets, &original);
    let rebase = |i: usize| {
        let arm = a.parts()[i]
            .pose_in_reference
            .position_world()
            .sub(a.parts()[2].pose_in_reference.position_world());
        twist(
            targets[i]
                .relative_velocity_assembly_m_s
                .sub(targets[2].relative_velocity_assembly_m_s)
                .sub(
                    targets[2]
                        .relative_angular_velocity_assembly_rad_s
                        .cross(arm),
                ),
            targets[i]
                .relative_angular_velocity_assembly_rad_s
                .sub(targets[2].relative_angular_velocity_assembly_rad_s),
        )
    };
    let changed = a
        .release_relative_twists(
            parent,
            2,
            &core::array::from_fn::<_, 3, _>(rebase),
            controls(),
            || false,
        )
        .unwrap();
    for (x, y) in original.release.bodies.iter().zip(&changed.release.bodies) {
        vector(x.linear_momentum_world(), y.linear_momentum_world(), 2e-12);
        vector(x.angular_momentum_body(), y.angular_momentum_body(), 2e-12);
    }
    close(
        original.release.supplied_work_j,
        changed.release.supplied_work_j,
        2e-12,
    );
    let order = [2, 0, 1];
    let reordered = assembly(&order.map(|i| raw[i]));
    let reordered_parent = state(
        &reordered,
        parent.pose().position_world(),
        q,
        parent
            .center_of_mass_velocity_world(a.properties())
            .unwrap(),
        Vec3::new(0.7, -0.8, 1.3),
    );
    let permuted = reordered
        .release_relative_twists(
            reordered_parent,
            1,
            &order.map(|i| targets[i]),
            controls(),
            || false,
        )
        .unwrap();
    for (index, &old) in order.iter().enumerate() {
        vector(
            permuted.release.bodies[index].linear_momentum_world(),
            original.release.bodies[old].linear_momentum_world(),
            2e-12,
        );
        vector(
            permuted.release.bodies[index].angular_momentum_body(),
            original.release.bodies[old].angular_momentum_body(),
            2e-12,
        );
    }
    let boost = Vec3::new(16.0, -8.0, 4.0);
    let boosted = state(
        &a,
        parent.pose().position_world(),
        q,
        parent
            .center_of_mass_velocity_world(a.properties())
            .unwrap()
            .add(boost),
        Vec3::new(0.7, -0.8, 1.3),
    );
    let moving = a
        .release_relative_twists(boosted, 0, &targets, controls(), || false)
        .unwrap();
    for ((body, old), part) in moving
        .release
        .bodies
        .iter()
        .zip(&original.release.bodies)
        .zip(a.parts())
    {
        vector(
            body.linear_momentum_world(),
            old.linear_momentum_world()
                .add(boost.scale(part.properties.mass())),
            2e-12,
        );
        vector(
            body.angular_momentum_body(),
            old.angular_momentum_body(),
            2e-12,
        );
    }
    close(
        moving.release.supplied_work_j,
        original.release.supplied_work_j,
        0.0,
    );
}

#[test]
fn g0_partial_groups_keep_internal_rigidity_and_their_complete_work() {
    let parts = [
        sphere(1.0, 0.4, Vec3::new(-1.5, -0.4, 0.2)),
        sphere(2.0, 0.7, Vec3::new(1.0, -0.3, 0.1)),
        sphere(3.0, 1.1, Vec3::new(-0.8, 0.6, -0.2)),
        sphere(4.0, 1.3, Vec3::new(1.4, 0.7, 0.3)),
    ];
    let a = assembly(&parts);
    let partition =
        AssemblyPartition::new(&a, &[0, 1, 0, 1], controls().assembly, || false).unwrap();
    let parent = state(
        &a,
        Vec3::new(2.0, -1.0, 0.5),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.4, -0.2, 0.3),
        Vec3::new(0.6, -0.5, 0.2),
    );
    let target = [
        RelativeSeparationTwist::ZERO,
        twist(Vec3::new(1.1, -0.3, 0.4), Vec3::new(0.2, -0.4, 0.7)),
    ];
    let event = partition
        .release_relative_twists(parent, 0, &target, controls(), || false)
        .unwrap();
    assert_eq!(event.release.bodies.len(), 2);
    let leaves = partition
        .expand(&event.release.bodies, controls().assembly, || false)
        .unwrap();
    let before = totals(&[parent], &[a.properties()]);
    let after = totals(&leaves, &parts.map(|part| part.properties));
    vector(after.0, before.0, 2e-8);
    vector(after.1, before.1, 2e-8);
    close(after.2 - before.2, event.release.supplied_work_j, 3e-8);
    for (group, assembly) in partition.groups().iter().enumerate() {
        let group_state = event.release.bodies[group];
        let group_spin = omega(group_state, assembly.properties());
        for &index in partition.original_part_indices(group).unwrap() {
            let arm = leaves[index]
                .pose()
                .position_world()
                .sub(group_state.pose().position_world());
            vector(
                leaves[index]
                    .center_of_mass_velocity_world(parts[index].properties)
                    .unwrap(),
                group_state
                    .center_of_mass_velocity_world(assembly.properties())
                    .unwrap()
                    .add(group_spin.cross(arm)),
                2e-12,
            );
            vector(
                omega(leaves[index], parts[index].properties),
                group_spin,
                2e-12,
            );
        }
    }
    let limited = RelativeSeparationControls {
        assembly: AssemblyControls {
            maximum_parts: 2,
            ..controls().assembly
        },
        ..controls()
    };
    assert!(matches!(
        partition.release_relative_twists(parent, 0, &target, limited, || false),
        Err(RelativeSeparationError::Assembly(AssemblyError::Budget {
            requested: 4,
            maximum: 2
        }))
    ));
}

#[test]
fn g0_dominant_mass_and_spin_inertia_keep_small_reference_recoil() {
    let a = assembly(&[sphere(1.0, 1.0, Vec3::ZERO), sphere(1e16, 1e16, Vec3::ZERO)]);
    let targets = [
        RelativeSeparationTwist::ZERO,
        twist(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    ];
    let result = a
        .release_relative_twists(rest(&a), 0, &targets, controls(), || false)
        .unwrap();
    // Directly subtracting the dominant target from its weighted mean rounds
    // both solved impulses to zero. The reduced quantities instead approach 1.
    vector(
        result.release.bodies[1].linear_momentum_world(),
        Vec3::new(1.0, 0.0, 0.0),
        3e-16,
    );
    vector(
        result.release.bodies[1].angular_momentum_body(),
        Vec3::new(0.0, 0.0, 1.0),
        3e-16,
    );
    vector(
        result.release.bodies[0].linear_momentum_world(),
        Vec3::new(-1.0, 0.0, 0.0),
        3e-16,
    );
    vector(
        result.release.bodies[0].angular_momentum_body(),
        Vec3::new(0.0, 0.0, -1.0),
        3e-16,
    );
    close(result.release.supplied_work_j, 1.0, 3e-16);
    requested_motion(&a, rest(&a), 0, &targets, &result);
}

#[test]
fn g4_budgets_sources_actual_endpoint_rates_and_terminal_cancel_are_atomic() {
    let (a, targets) = two();
    let saved = a.clone();
    let parent = rest(&a);
    let limited = RelativeSeparationControls {
        assembly: AssemblyControls {
            maximum_impulses: 0,
            ..controls().assembly
        },
        ..controls()
    };
    assert!(matches!(
        a.release_relative_twists(parent, 0, &targets, limited, || false),
        Err(RelativeSeparationError::Assembly(AssemblyError::Budget {
            requested: 1,
            maximum: 0
        }))
    ));
    let identity = a
        .release_relative_twists(
            parent,
            0,
            &[RelativeSeparationTwist::ZERO; 2],
            limited,
            || false,
        )
        .unwrap();
    assert!(identity.impulses.is_empty());
    close(identity.release.supplied_work_j, 0.0, 0.0);
    assert!(matches!(
        a.release_relative_twists(
            parent,
            0,
            &targets,
            RelativeSeparationControls {
                maximum_supplied_work_j: 19.0,
                ..controls()
            },
            || false
        ),
        Err(RelativeSeparationError::SourceWork { .. })
    ));
    assert!(
        a.release_relative_twists(parent, 2, &targets, controls(), || false)
            .is_err()
    );
    assert!(
        a.release_relative_twists(parent, 0, &targets[..1], controls(), || false)
            .is_err()
    );
    let mut invalid = targets;
    invalid[0] = targets[1];
    assert!(
        a.release_relative_twists(parent, 0, &invalid, controls(), || false)
            .is_err()
    );
    invalid = targets;
    invalid[1].relative_angular_velocity_assembly_rad_s.x = f64::NAN;
    assert!(
        a.release_relative_twists(parent, 0, &invalid, controls(), || false)
            .is_err()
    );
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            a.release_relative_twists(
                parent,
                0,
                &targets,
                RelativeSeparationControls {
                    relative_angular_velocity_tolerance_rad_s: bad,
                    ..controls()
                },
                || false
            )
            .is_err()
        );
    }
    // Relax native energy admission deliberately: the new rate gates must
    // still reject kicks erased in the ACTUAL canonical endpoint storage.
    let rate_gate = RelativeSeparationControls {
        assembly: AssemblyControls {
            energy_tolerance_j: 100.0,
            ..controls().assembly
        },
        ..controls()
    };
    let enormous = 2.0_f64.powi(60);
    let boosted = state(
        &a,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(enormous, 0.0, 0.0),
        Vec3::ZERO,
    );
    assert!(matches!(
        a.release_relative_twists(boosted, 0, &targets, rate_gate, || false),
        Err(RelativeSeparationError::RelativeVelocity { .. })
    ));
    let spinning = state(
        &a,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, enormous),
    );
    let spins = [
        RelativeSeparationTwist::ZERO,
        twist(Vec3::ZERO, Vec3::new(0.0, 0.0, 2.0)),
    ];
    assert!(matches!(
        a.release_relative_twists(spinning, 0, &spins, rate_gate, || false),
        Err(RelativeSeparationError::RelativeAngularVelocity { .. })
    ));
    let mut calls = 0;
    let completed = a
        .release_relative_twists(parent, 0, &targets, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        a.release_relative_twists(parent, 0, &targets, controls(), || {
            current += 1;
            current == calls
        }),
        Err(RelativeSeparationError::Assembly(AssemblyError::Cancelled))
    ));
    assert_eq!(current, calls);
    assert_eq!(a, saved);
    assert_eq!(
        a.release_relative_twists(parent, 0, &targets, controls(), || false)
            .unwrap(),
        completed
    );
}
