//! Actual body loads, contact timing and external work retain native balances.
use super::*;
use fs_mbd::Wrench;

fn axial_load(force: f64) -> [Wrench; 2] {
    [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, force),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::new(0.0, 0.0, -force),
            torque_body: Vec3::ZERO,
        },
    ]
}

fn check_loaded(
    f: &Fixture,
    c: &RigidJkrCandidate,
    wrenches: [Wrench; 2],
    dt: f64,
    normal: &SourceJkrSphereContact,
) {
    let limits = controls(normal);
    let arms = f.prepared(normal).reaction_arms_body();
    let mut kinetic = 0.0;
    let mut external_work = 0.0;
    let mut pair_linear = Vec3::ZERO;
    let mut pair_angular = Vec3::ZERO;
    let mut expected_linear = Vec3::ZERO;
    let mut expected_angular = Vec3::ZERO;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i];
        let after = c.bodies_after[i];
        let pose = before.pose();
        let contact_impulse = c.report.force_impulse_on_a_world_ns.scale(sign);
        let external_impulse = wrenches[i].force_world.scale(dt);
        let external_torque = wrenches[i].torque_body.scale(dt);
        let expected_torque = arms[i]
            .cross(pose.orientation().rotate_world_to_body(contact_impulse))
            .add(external_torque);
        let dp = after
            .linear_momentum_world()
            .sub(before.linear_momentum_world());
        let dl = after
            .angular_momentum_body()
            .sub(before.angular_momentum_body());
        close_vec(
            dp,
            contact_impulse.add(external_impulse),
            limits.linear_momentum_tolerance_ns,
        );
        close_vec(dl, expected_torque, limits.angular_momentum_tolerance_nms);
        assert_eq!(after.pose(), before.pose());
        kinetic += after.kinetic_energy(f.properties[i]).unwrap()
            - before.kinetic_energy(f.properties[i]).unwrap();
        external_work += wrenches[i]
            .force_world
            .dot(c.report.body_displacement_world_m[i])
            + wrenches[i]
                .torque_body
                .dot(c.report.angular_displacement_body_rad[i]);
        pair_linear = pair_linear.add(dp);
        pair_angular = pair_angular
            .add(pose.position_world().cross(dp))
            .add(pose.orientation().rotate_body_to_world(dl));
        expected_linear = expected_linear.add(external_impulse);
        expected_angular = expected_angular
            .add(pose.position_world().cross(external_impulse))
            .add(pose.orientation().rotate_body_to_world(external_torque));
    }
    assert_eq!(c.report.external_wrenches, wrenches);
    close(
        kinetic,
        c.report.kinetic_energy_change_j,
        limits.energy_tolerance_j,
    );
    close(
        external_work,
        c.report.external_work_j,
        limits.energy_tolerance_j,
    );
    close(
        kinetic + c.report.normal_storage_change_j + c.report.motion.released_energy_j,
        external_work,
        limits.energy_tolerance_j,
    );
    close(
        kinetic + c.report.motion.normal_work_j,
        external_work,
        limits.energy_tolerance_j,
    );
    close_vec(
        pair_linear,
        expected_linear,
        limits.linear_momentum_tolerance_ns,
    );
    close_vec(
        pair_angular,
        expected_angular,
        limits.angular_momentum_tolerance_nms,
    );
    close_vec(
        c.report.audit.linear_momentum_change_ns,
        expected_linear,
        limits.linear_momentum_tolerance_ns,
    );
    close_vec(
        c.report.audit.angular_momentum_change_nms,
        expected_angular,
        limits.angular_momentum_tolerance_nms,
    );
}

#[test]
fn zero_wrenches_preserve_original_api_and_common_gravity_preserves_relative_contact() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    let f = axial(0.0);
    let pair = f.prepared(&normal);
    let initial = normal.evaluate_attached(pair.approach_m()).unwrap().state();
    let dt = 0.2 * time;
    let unloaded = pair.advance(initial, dt, &CancelGate::new()).unwrap();
    assert_eq!(
        unloaded,
        pair.advance_with_wrenches(initial, [Wrench::ZERO; 2], dt, &CancelGate::new())
            .unwrap()
    );
    let gravity = Vec3::new(3.0, -7.0, -9.81);
    let wrenches = core::array::from_fn(|i| Wrench {
        force_world: gravity.scale(f.properties[i].mass()),
        torque_body: Vec3::ZERO,
    });
    let loaded = pair
        .advance_with_wrenches(initial, wrenches, dt, &CancelGate::new())
        .unwrap();
    close(loaded.report.applied_normal_force_n, 0.0, 1e-12 * FORCE);
    close(
        loaded.report.motion.normal_impulse_n_s,
        unloaded.report.motion.normal_impulse_n_s,
        controls(&normal).motion.impulse_absolute_tolerance_n_s,
    );
    close(
        loaded.state_after.indentation_m(),
        unloaded.state_after.indentation_m(),
        1e-12 * length,
    );
    for i in 0..2 {
        close_vec(
            loaded.bodies_after[i].linear_momentum_world(),
            unloaded.bodies_after[i]
                .linear_momentum_world()
                .add(wrenches[i].force_world.scale(dt)),
            1e-12 * MASS * speed,
        );
        close_vec(
            loaded.report.body_displacement_world_m[i],
            unloaded.report.body_displacement_world_m[i].add(gravity.scale(0.5 * dt * dt)),
            1e-12 * length,
        );
        assert_eq!(loaded.bodies_after[i].angular_momentum_body(), Vec3::ZERO);
    }
    close(
        loaded.report.external_work_j,
        0.5 * (f.properties[0].mass() + f.properties[1].mass()) * gravity.norm_squared() * dt * dt,
        1e-11 * FORCE * length,
    );
    check_loaded(&f, &loaded, wrenches, dt, &normal);
}

#[test]
fn actual_tensile_preload_balances_adhesion_and_stronger_loading_causes_original_pull_off() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    let f = axial(0.0);
    let pair = f.prepared(&normal);
    let response = normal.evaluate_attached(pair.approach_m()).unwrap();
    let initial = response.state();
    let hold = axial_load(response.normal_force_n);
    let held = pair
        .advance_with_wrenches(initial, hold, 0.5 * time, &CancelGate::new())
        .unwrap();
    close(
        held.state_after.indentation_m(),
        pair.approach_m(),
        1e-10 * length,
    );
    close(
        held.report.relative_normal_velocity_m_s[1],
        0.0,
        1e-10 * speed,
    );
    assert_eq!(held.report.motion.event_count(), 0);
    assert!(held.report.force_impulse_on_a_world_ns.z > 0.0);
    check_loaded(&f, &held, hold, 0.5 * time, &normal);

    let pulling = axial_load(-12.0 * FORCE);
    assert!(pulling[0].force_world.z < normal.force_controlled_pull_off_n());
    let dt = 3.0 * time;
    let pulled = pair
        .advance_with_wrenches(initial, pulling, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(pulled.state_after.branch(), JkrBranch::Detached);
    assert_eq!(pulled.report.motion.event_count(), 1);
    close(
        pulled
            .report
            .motion
            .events()
            .next()
            .unwrap()
            .before()
            .indentation_m(),
        normal.minimum_attached_indentation_m(),
        1e-12 * length,
    );
    close(
        pulled.report.motion.released_energy_j,
        2.4 * FORCE * length,
        1e-11 * FORCE * length,
    );
    assert!(pulled.report.relative_normal_velocity_m_s[1] < 0.0);
    assert!(pulled.report.external_work_j > 0.0);
    check_loaded(&f, &pulled, pulling, dt, &normal);
    assert_eq!(
        pulled,
        pair.advance_with_wrenches(initial, pulling, dt, &CancelGate::new())
            .unwrap()
    );
}

#[test]
fn off_axis_full_inertia_loads_use_actual_torque_and_first_moment_external_work() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    let orientations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.3).unwrap(),
    ];
    let arms = [
        Vec3::new(0.4, 0.2, 1.0).scale(ARM),
        Vec3::new(-0.3, 0.1, -1.0).scale(ARM),
    ];
    let mut f = fixture(0.0, arms, orientations);
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.6).unwrap();
    f.properties[0] = MassProperties::from_principal_frame(
        2.0 * MASS,
        Vec3::ZERO,
        Vec3::new(0.8, 1.1, 1.4).scale(INERTIA),
        principal,
    )
    .unwrap();
    f.set_motion(
        0,
        Vec3::new(0.3, -0.2, 1.0).scale(speed),
        Vec3::new(0.05, 0.02, -0.03).scale(speed / ARM),
    );
    let inverse = |i: usize, x: Vec3| {
        let rotation = if i == 0 {
            principal
        } else {
            UnitQuaternion::IDENTITY
        };
        let v = rotation.rotate_world_to_body(x);
        let d = f.properties[i].principal_inertia_body();
        rotation.rotate_body_to_world(Vec3::new(v.x / d.x, v.y / d.y, v.z / d.z))
    };
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.5, -0.4, 2.0).scale(FORCE),
            torque_body: Vec3::new(0.3, -0.4, 0.2).scale(FORCE * ARM),
        },
        Wrench {
            force_world: Vec3::new(-0.3, 0.2, 0.1).scale(FORCE),
            torque_body: Vec3::new(-0.2, 0.1, 0.4).scale(FORCE * ARM),
        },
    ];
    let pair = f.prepared(&normal);
    let n = Vec3::new(0.0, 0.0, 1.0);
    let mut mobility = 0.0;
    let mut normal_acceleration = 0.0;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let arm = pair.reaction_arms_body()[i];
        let covector = arm.cross(orientations[i].rotate_world_to_body(n));
        mobility += 1.0 / f.properties[i].mass() + covector.dot(inverse(i, covector));
        let acceleration = wrenches[i]
            .force_world
            .scale(f.properties[i].mass().recip())
            .add(
                orientations[i]
                    .rotate_body_to_world(inverse(i, wrenches[i].torque_body).cross(arm)),
            );
        normal_acceleration += sign * acceleration.dot(n);
    }
    let initial = normal.evaluate_attached(pair.approach_m()).unwrap().state();
    let dt = 0.5 * time;
    let c = pair
        .advance_with_wrenches(initial, wrenches, dt, &CancelGate::new())
        .unwrap();
    close(
        c.report.applied_normal_force_n,
        normal_acceleration / mobility,
        1e-12 * FORCE,
    );
    let q = c
        .report
        .force_impulse_on_a_world_ns
        .scale(dt)
        .sub(c.report.force_time_moment_on_a_world_ns2);
    let mut midpoint_external_work = 0.0;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i]
            .point_kinematics(f.properties[i], Vec3::ZERO)
            .unwrap();
        let after = c.bodies_after[i]
            .point_kinematics(f.properties[i], Vec3::ZERO)
            .unwrap();
        let expected_com = before.center_of_mass_velocity_world.scale(dt).add(
            q.scale(sign)
                .add(wrenches[i].force_world.scale(0.5 * dt * dt))
                .scale(f.properties[i].mass().recip()),
        );
        let torque_integral = pair.reaction_arms_body()[i]
            .cross(orientations[i].rotate_world_to_body(q.scale(sign)))
            .add(wrenches[i].torque_body.scale(0.5 * dt * dt));
        let expected_angle = before
            .angular_velocity_body
            .scale(dt)
            .add(inverse(i, torque_integral));
        close_vec(
            c.report.body_displacement_world_m[i],
            expected_com,
            1e-12 * length,
        );
        close_vec(
            c.report.angular_displacement_body_rad[i],
            expected_angle,
            1e-12 * length / ARM,
        );
        midpoint_external_work += 0.5
            * dt
            * (wrenches[i].force_world.dot(
                before
                    .center_of_mass_velocity_world
                    .add(after.center_of_mass_velocity_world),
            ) + wrenches[i].torque_body.dot(
                before
                    .angular_velocity_body
                    .add(after.angular_velocity_body),
            ));
    }
    // A midpoint load-work shortcut loses the curved contact response even
    // though the external loads themselves are constant during this phase.
    assert!(
        (midpoint_external_work - c.report.external_work_j).abs()
            > 100.0 * controls(&normal).energy_tolerance_j
    );
    check_loaded(&f, &c, wrenches, dt, &normal);
    let mut reversed_contact = f.contact(&normal);
    reversed_contact.binding.body_a = 1;
    reversed_contact.binding.body_b = 0;
    assert_eq!(
        c,
        prepare_rigid_jkr_contact_with_wrenches(
            &[f.bodies[1], f.bodies[0]],
            &[f.properties[1], f.properties[0]],
            reversed_contact,
            initial,
            wrenches,
            dt,
            controls(&normal),
            &CancelGate::new(),
        )
        .unwrap()
    );
}

#[test]
fn actual_load_domain_and_cancellation_refusals_leave_both_bodies_and_history_untouched() {
    let normal = source(1e7);
    let (length, _, time) = scales(&normal);
    let f = axial(0.0);
    let initial = normal.evaluate_attached(0.0).unwrap().state();
    let wrenches = axial_load(-12.0 * FORCE);
    let mut native = f.bodies;
    let mut history = initial;
    let mut nan_force = wrenches;
    nan_force[0].force_world.x = f64::NAN;
    let mut infinite_torque = wrenches;
    infinite_torque[1].torque_body.y = f64::INFINITY;
    for invalid in [nan_force, infinite_torque, axial_load(-1e6 * FORCE)] {
        assert!(
            apply_rigid_jkr_contact_with_wrenches(
                &mut native,
                &f.properties,
                f.contact(&normal),
                &mut history,
                invalid,
                time,
                controls(&normal),
                &CancelGate::new(),
            )
            .is_err()
        );
        assert_eq!(native, f.bodies);
        assert_eq!(history, initial);
    }
    let mut limited = controls(&normal);
    limited.motion.minimum_indentation_m = -3.1 * length;
    assert!(
        apply_rigid_jkr_contact_with_wrenches(
            &mut native,
            &f.properties,
            f.contact(&normal),
            &mut history,
            wrenches,
            3.0 * time,
            limited,
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_jkr_contact_with_wrenches(
            &mut native,
            &f.properties,
            f.contact(&normal),
            &mut history,
            wrenches,
            time,
            controls(&normal),
            &cancelled,
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let expected = f
        .prepared(&normal)
        .advance_with_wrenches(initial, wrenches, time, &CancelGate::new())
        .unwrap();
    let report = apply_rigid_jkr_contact_with_wrenches(
        &mut native,
        &f.properties,
        f.contact(&normal),
        &mut history,
        wrenches,
        time,
        controls(&normal),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(native, expected.bodies_after);
    assert_eq!(history, expected.state_after);
    assert_eq!(report, expected.report);
}
