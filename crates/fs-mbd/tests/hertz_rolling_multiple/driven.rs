//! Powered shared-body assemblies retain physical coupled motion and work.

use super::*;
use fs_mbd::Wrench;

fn solve(f: &Fixture, external: &[Wrench]) -> RigidHertzRollingIslandCandidate {
    PreparedRigidHertzRollingIsland::new_driven(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        external,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.laws(), &CancelGate::new())
    .unwrap()
}

fn motor() -> [Wrench; 3] {
    [
        Wrench {
            force_world: Vec3::new(20.0, 10.0, 0.0),
            torque_body: Vec3::new(0.2, 0.3, 0.1),
        },
        Wrench {
            force_world: Vec3::new(0.4, -0.2, 0.0),
            torque_body: Vec3::new(0.002, -0.004, 0.003),
        },
        Wrench {
            force_world: Vec3::new(-0.2, 0.3, 0.0),
            torque_body: Vec3::new(-0.003, 0.002, -0.004),
        },
    ]
}

fn mixed() -> Fixture {
    let mut f = fixture(
        0.01,
        [0.12, 0.2],
        Vec3::new(0.012, 0.003, 0.0),
        Vec3::new(0.015, 0.025, 0.1),
    );
    for (i, principal, angle) in [
        (1, Vec3::new(0.08, 0.13, 0.18), 0.43),
        (2, Vec3::new(0.17, 0.21, 0.25), -0.31),
    ] {
        f.properties[i] = MassProperties::from_principal_frame(
            f.properties[i].mass(),
            Vec3::ZERO,
            principal,
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), angle).unwrap(),
        )
        .unwrap();
    }
    f
}

// Rebuild every native resultant and physical work from the common endpoints,
// without the production response matrix or any per-contact effective mass.
fn body_balance(
    before: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidHertzRollingIslandContact<'_>],
    external: &[Wrench],
    result: &RigidHertzRollingIslandCandidate,
) {
    let mut linear = vec![Vec3::ZERO; before.len()];
    let mut angular = vec![Vec3::ZERO; before.len()];
    let mut contact_work = 0.0;
    let after =
        |global: usize| result.bodies_after[result.body_indices.binary_search(&global).unwrap()];
    for report in &result.report.contacts {
        let original = contacts
            .iter()
            .find(|c| c.patch_id == report.patch_id)
            .unwrap()
            .contact;
        let binding = original.binding;
        let ids = [binding.body_a, binding.body_b];
        let arms = [binding.arm_a_body, binding.arm_b_body];
        let force = report.force_impulse_world_ns;
        let couple = report
            .spin_angular_impulse_world_nms
            .add(report.rolling_angular_impulse_world_nms);
        for (side, sign) in [(0, 1.0), (1, -1.0)] {
            let body = ids[side];
            let rotation = before[body].pose().orientation();
            linear[body] = linear[body].add(force.scale(sign));
            angular[body] = angular[body]
                .add(arms[side].cross(rotation.rotate_world_to_body(force.scale(sign))))
                .add(rotation.rotate_world_to_body(couple.scale(sign)));
            let old = before[body]
                .point_kinematics(properties[body], arms[side])
                .unwrap();
            let new = after(body)
                .point_kinematics(properties[body], arms[side])
                .unwrap();
            contact_work += force.scale(sign).dot(
                old.point_velocity_world
                    .add(new.point_velocity_world)
                    .scale(0.5),
            ) + couple.scale(sign).dot(
                old.angular_velocity_world
                    .add(new.angular_velocity_world)
                    .scale(0.5),
            );
        }
    }
    let (mut kinetic_change, mut external_work) = (0.0, 0.0);
    let (mut commanded_linear, mut commanded_angular) = (Vec3::ZERO, Vec3::ZERO);
    for &i in &result.body_indices {
        let body = after(i);
        let force = external[i].force_world.scale(DT);
        let torque = external[i].torque_body.scale(DT);
        close_vec(
            body.linear_momentum_world()
                .sub(before[i].linear_momentum_world()),
            linear[i].add(force),
            2e-12,
        );
        close_vec(
            body.angular_momentum_body()
                .sub(before[i].angular_momentum_body()),
            angular[i].add(torque),
            2e-12,
        );
        let mean_v = before[i]
            .center_of_mass_velocity_world(properties[i])
            .unwrap()
            .add(body.center_of_mass_velocity_world(properties[i]).unwrap())
            .scale(0.5);
        let mean_w = properties[i]
            .angular_velocity_body_checked(before[i].angular_momentum_body())
            .unwrap()
            .add(
                properties[i]
                    .angular_velocity_body_checked(body.angular_momentum_body())
                    .unwrap(),
            )
            .scale(0.5);
        external_work += force.dot(mean_v) + torque.dot(mean_w);
        kinetic_change += body.kinetic_energy(properties[i]).unwrap()
            - before[i].kinetic_energy(properties[i]).unwrap();
        commanded_linear = commanded_linear.add(force);
        commanded_angular = commanded_angular
            .add(before[i].pose().position_world().cross(force))
            .add(before[i].pose().orientation().rotate_body_to_world(torque));
    }
    close(result.report.external_work_j, external_work, 3e-13);
    close(result.report.midpoint_impulse_work_j, contact_work, 3e-13);
    close(kinetic_change - contact_work - external_work, 0.0, 3e-11);
    close(
        kinetic_change + result.report.surface_heat_j + result.report.bulk_heat_j - external_work,
        0.0,
        3e-11,
    );
    close_vec(
        result.report.external_linear_impulse_ns,
        commanded_linear,
        2e-12,
    );
    close_vec(
        result.report.external_angular_impulse_nms,
        commanded_angular,
        2e-12,
    );
    close_vec(
        result.report.linear_momentum_change_ns,
        commanded_linear,
        2e-12,
    );
    close_vec(
        result.report.angular_momentum_change_nms,
        commanded_angular,
        2e-12,
    );
}

#[test]
fn shared_sphere_motor_matches_independent_coupled_rolling_equations_from_rest() {
    let f = fixture(0.01, [0.0; 2], Vec3::ZERO, Vec3::ZERO);
    let torque = [0.8, -0.06, 0.04];
    let external = torque.map(|t| Wrench {
        force_world: Vec3::ZERO,
        torque_body: Vec3::new(0.0, t, 0.0),
    });
    let result = solve(&f, &external);
    body_balance(&f.bodies, &f.properties, &f.contacts(), &external, &result);
    assert_eq!(result.report.surface_heat_j, 0.0);
    assert!(result.report.bulk_heat_j > 0.0 && result.report.external_work_j > 0.0);
    let is = f.properties[0]
        .angular_momentum_body_checked(Vec3::new(0.0, 1.0, 0.0))
        .unwrap()
        .y;
    let ib = [0.13, 0.21];
    let c = f.patches.each_ref().map(|p| 0.01 * p.normal_load_n() * R);
    let free = [0, 1].map(|i| 0.5 * DT * (torque[0] / is - torque[i + 1] / ib[i]));
    let a = 1.0 + 0.5 * DT * c[0] * (1.0 / is + 1.0 / ib[0]);
    let b = 0.5 * DT * c[0] / is;
    let cc = 0.5 * DT * c[1] / is;
    let d = 1.0 + 0.5 * DT * c[1] * (1.0 / is + 1.0 / ib[1]);
    let rhs = [-c[0] * free[0], -c[1] * free[1]];
    let torques = [
        (rhs[0] * d - b * rhs[1]) / (a * d - b * cc),
        (a * rhs[1] - cc * rhs[0]) / (a * d - b * cc),
    ];
    close(
        result.bodies_after[0].angular_momentum_body().y,
        DT * (torque[0] + torques[0] + torques[1]),
        2e-16,
    );
    for i in 0..2 {
        close(
            result.bodies_after[i + 1].angular_momentum_body().y,
            DT * (torque[i + 1] - torques[i]),
            2e-16,
        );
        close(
            result.report.contacts[i]
                .rolling_angular_impulse_world_nms
                .y,
            DT * torques[i],
            2e-16,
        );
        close(
            result.report.contacts[i].bulk_heat_j,
            DT * torques[i].powi(2) / c[i],
            2e-18,
        );
    }

    // Two pads on one enclosing support give a dependent ten-coordinate
    // response. The SAME body motor still contributes once to the equations.
    let mut dependent = fixture(0.01, [0.0; 2], Vec3::ZERO, Vec3::ZERO);
    dependent.bodies[1] = RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::ZERO).unwrap();
    for i in 0..2 {
        dependent.bindings[i].body_b = 1;
        dependent.bindings[i].arm_b_body = dependent.bindings[i].arm_a_body;
    }
    let drive = [external[0], external[1], Wrench::ZERO];
    let reduced = solve(&dependent, &drive);
    assert_eq!(reduced.body_indices, vec![0, 1]);
    let midpoint = free[0] / (1.0 + 0.5 * DT * (c[0] + c[1]) * (1.0 / is + 1.0 / ib[0]));
    close(
        reduced.bodies_after[0].angular_momentum_body().y,
        DT * (torque[0] - (c[0] + c[1]) * midpoint),
        2e-16,
    );
    close(
        reduced.report.bulk_heat_j,
        DT * (c[0] + c[1]) * midpoint.powi(2),
        2e-18,
    );
    body_balance(
        &dependent.bodies,
        &dependent.properties,
        &dependent.contacts(),
        &drive,
        &reduced,
    );
}

#[test]
fn powered_mixed_shear_spin_and_rolling_use_full_inertia_and_one_body_work_port() {
    let f = mixed();
    let external = motor();
    let result = solve(&f, &external);
    body_balance(&f.bodies, &f.properties, &f.contacts(), &external, &result);
    assert!(result.report.kinetic_after_j > result.report.kinetic_before_j);
    assert!(result.report.contacts.iter().all(|c| c.surface_heat_j > 0.0
        && c.bulk_heat_j > 0.0
        && c.spin_angular_impulse_world_nms != Vec3::ZERO));
    assert!(
        result
            .report
            .contacts
            .iter()
            .all(|c| c.scaled_driving_acceleration_m_s2.iter().all(|a| *a != 0.0))
    );
    let prepared = PreparedRigidHertzRollingIsland::new_driven(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        &external,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let w = prepared.scaled_mobility_per_kg();
    assert!(w[5] != 0.0 && w[4 * 10 + 9] != 0.0);
    let common = Vec3::new(0.4, -0.2, 0.3);
    let boosted = [0, 1, 2].map(|i| Wrench {
        force_world: external[i]
            .force_world
            .add(common.scale(f.properties[i].mass())),
        torque_body: external[i].torque_body,
    });
    let accelerated = solve(&f, &boosted);
    for i in 0..3 {
        close_vec(
            accelerated.bodies_after[i].linear_momentum_world(),
            result.bodies_after[i]
                .linear_momentum_world()
                .add(common.scale(DT * f.properties[i].mass())),
            2e-12,
        );
        close_vec(
            accelerated.bodies_after[i].angular_momentum_body(),
            result.bodies_after[i].angular_momentum_body(),
            2e-12,
        );
    }
    close(
        accelerated.report.surface_heat_j,
        result.report.surface_heat_j,
        3e-14,
    );
    close(
        accelerated.report.bulk_heat_j,
        result.report.bulk_heat_j,
        3e-14,
    );
    body_balance(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        &boosted,
        &accelerated,
    );
}

#[test]
fn driven_patch_order_global_body_indices_singleton_and_zero_drive_are_preserved() {
    let f = mixed();
    let external = motor();
    let result = solve(&f, &external);
    let mut contacts = f.contacts();
    let mut laws = f.laws();
    contacts.reverse();
    laws.reverse();
    let reordered = PreparedRigidHertzRollingIsland::new_driven(
        &f.bodies,
        &f.properties,
        &contacts,
        &external,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&laws, &CancelGate::new())
    .unwrap();
    assert_eq!(result, reordered);
    assert_eq!(
        solve(&f, &[Wrench::ZERO; 3]),
        f.prepare().advance(&f.laws(), &CancelGate::new()).unwrap()
    );

    let spare =
        RigidBodyState::new(Pose::identity(), Vec3::new(4.0, 5.0, 6.0), Vec3::ZERO).unwrap();
    let bodies = [spare, f.bodies[2], f.bodies[0], f.bodies[1]];
    let properties = [
        f.properties[1],
        f.properties[2],
        f.properties[0],
        f.properties[1],
    ];
    let drives = [Wrench::ZERO, external[2], external[0], external[1]];
    let map = [2, 3, 1];
    let mut mapped = f.contacts();
    for c in &mut mapped {
        c.contact.binding.body_a = map[c.contact.binding.body_a];
        c.contact.binding.body_b = map[c.contact.binding.body_b];
    }
    let remapped = PreparedRigidHertzRollingIsland::new_driven(
        &bodies,
        &properties,
        &mapped,
        &drives,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&laws, &CancelGate::new())
    .unwrap();
    assert_eq!(remapped.body_indices, vec![1, 2, 3]);
    for i in 0..3 {
        assert_eq!(
            remapped.bodies_after[remapped.body_indices.binary_search(&map[i]).unwrap()],
            result.bodies_after[i]
        );
    }
    assert_eq!(bodies[0], spare);
    body_balance(&bodies, &properties, &mapped, &drives, &remapped);

    let singleton_drive = [external[0], external[1], Wrench::ZERO];
    let single = PreparedRigidHertzRollingIsland::new_driven(
        &f.bodies,
        &f.properties,
        &f.contacts()[..1],
        &singleton_drive,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.laws()[..1], &CancelGate::new())
    .unwrap();
    let pair = PreparedRigidHertzRollingPair::new_driven(
        &f.bodies,
        &f.properties,
        &f.histories[0],
        f.contact(0),
        [external[0], external[1]],
        DT,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.hertz[0], &CancelGate::new())
    .unwrap();
    for i in 0..2 {
        close_vec(
            single.bodies_after[i].linear_momentum_world(),
            pair.bodies_after[i].linear_momentum_world(),
            3e-12,
        );
        close_vec(
            single.bodies_after[i].angular_momentum_body(),
            pair.bodies_after[i].angular_momentum_body(),
            3e-12,
        );
    }
    close(
        single.report.surface_heat_j,
        pair.report.surface_heat_j,
        3e-12,
    );
    close(single.report.bulk_heat_j, pair.report.bulk_heat_j, 3e-12);
    close(
        single.report.external_work_j,
        pair.report.external_wrench_work_j,
        3e-12,
    );
}

#[test]
fn invalid_drive_unresolved_normal_motion_and_slip_events_preserve_retry() {
    let f = mixed();
    let bodies = f.bodies;
    let histories = f.histories.clone();
    let external = motor();
    let prepared = PreparedRigidHertzRollingIsland::new_driven(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        &external,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let accepted = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    assert!(
        PreparedRigidHertzRollingIsland::new_driven(
            &f.bodies,
            &f.properties,
            &f.contacts(),
            &external[..2],
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert!(
        PreparedRigidHertzRollingIsland::new_driven(
            &f.bodies,
            &f.properties,
            &f.contacts()[..1],
            &external,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut invalid = external;
    invalid[0].torque_body.x = f64::NAN;
    assert!(
        PreparedRigidHertzRollingIsland::new_driven(
            &f.bodies,
            &f.properties,
            &f.contacts(),
            &invalid,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut normal = external;
    normal[0].force_world.z = 1.0;
    let unsupported = PreparedRigidHertzRollingIsland::new_driven(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        &normal,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(unsupported.advance(&f.laws(), &CancelGate::new()).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(prepared.advance(&f.laws(), &cancelled).is_err());
    let braking = fixture(0.0, [0.1, 0.2], Vec3::new(0.012, 0.0, 0.0), Vec3::ZERO);
    let brakes = [
        Wrench {
            force_world: Vec3::new(-500.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
        Wrench::ZERO,
    ];
    let coarse = PreparedRigidHertzRollingIsland::new_driven(
        &braking.bodies,
        &braking.properties,
        &braking.contacts(),
        &brakes,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        coarse.advance(&braking.laws(), &CancelGate::new()),
        Err(RigidHertzRollingIslandError::UnresolvedSlipEvent { .. })
    ));
    assert_eq!(f.bodies, bodies);
    assert_eq!(f.histories, histories);
    assert_eq!(
        prepared.advance(&f.laws(), &CancelGate::new()).unwrap(),
        accepted
    );
}
