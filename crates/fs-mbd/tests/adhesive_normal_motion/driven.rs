//! Original DMT source paths driven by actual COM forces and body torques.
use super::*;
use fs_mbd::Wrench;

// Independent principal-frame reconstruction, including products of inertia.
fn inverse_inertia(properties: MassProperties, torque: Vec3) -> Vec3 {
    let principal = properties.principal_axes_to_body();
    let t = principal.rotate_world_to_body(torque);
    let d = properties.principal_inertia_body();
    principal.rotate_body_to_world(Vec3::new(t.x / d.x, t.y / d.y, t.z / d.z))
}

fn external_point_acceleration(f: &Fixture, wrenches: [Wrench; 2]) -> Vec3 {
    let a: [Vec3; 2] = core::array::from_fn(|i| {
        wrenches[i]
            .force_world
            .scale(1.0 / f.properties[i].mass())
            .add(f.bodies[i].pose().orientation().rotate_body_to_world(
                inverse_inertia(f.properties[i], wrenches[i].torque_body).cross(f.reaction_arms[i]),
            ))
    });
    a[0].sub(a[1])
}

fn full_mobility(f: &Fixture) -> [[f64; 3]; 3] {
    let axes = [f.frame.normal(), f.frame.longitudinal(), f.frame.lateral()].map(vector);
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..2)
                .map(|body| {
                    let q = f.bodies[body].pose().orientation();
                    let a = f.reaction_arms[body].cross(q.rotate_world_to_body(axes[i]));
                    let b = f.reaction_arms[body].cross(q.rotate_world_to_body(axes[j]));
                    axes[i].dot(axes[j]) / f.properties[body].mass()
                        + a.dot(inverse_inertia(f.properties[body], b))
                })
                .sum()
        })
    })
}

// Reconstruct physical increments and signed work from native endpoints, not
// the source's stored energy or the production audit's residuals.
fn check_driven(f: &Fixture, c: &RigidDmtCandidate, wrenches: [Wrench; 2], dt: f64) {
    let mut external_work = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut expected_linear = Vec3::ZERO;
    let mut expected_angular = Vec3::ZERO;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i];
        let after = c.bodies_after[i];
        let orientation = before.pose().orientation();
        let contact_impulse = c.report.force_impulse_on_a_world_ns.scale(sign);
        let external_impulse = wrenches[i].force_world.scale(dt);
        let external_torque = wrenches[i].torque_body.scale(dt);
        let dp = after
            .linear_momentum_world()
            .sub(before.linear_momentum_world());
        let dl = after
            .angular_momentum_body()
            .sub(before.angular_momentum_body());
        close_vec(dp, contact_impulse.add(external_impulse), 1e-29);
        close_vec(
            dl,
            f.reaction_arms[i]
                .cross(orientation.rotate_world_to_body(contact_impulse))
                .add(external_torque),
            1e-34,
        );
        let travel = before
            .linear_momentum_world()
            .add(after.linear_momentum_world())
            .scale(0.5 * dt / f.properties[i].mass());
        let angular_travel = inverse_inertia(
            f.properties[i],
            before
                .angular_momentum_body()
                .add(after.angular_momentum_body()),
        )
        .scale(0.5 * dt);
        close_vec(c.report.body_displacement_world_m[i], travel, 2e-21);
        close_vec(
            c.report.angular_displacement_body_rad[i],
            angular_travel,
            1e-15,
        );
        external_work +=
            wrenches[i].force_world.dot(travel) + wrenches[i].torque_body.dot(angular_travel);
        linear = linear.add(dp);
        angular = angular
            .add(before.pose().position_world().cross(dp))
            .add(orientation.rotate_body_to_world(dl));
        expected_linear = expected_linear.add(external_impulse);
        expected_angular = expected_angular
            .add(before.pose().position_world().cross(external_impulse))
            .add(orientation.rotate_body_to_world(external_torque));
        assert_eq!(after.pose(), before.pose());
    }
    close_vec(linear, expected_linear, 1e-29);
    close_vec(angular, expected_angular, 1e-34);
    close_vec(c.report.audit.linear_momentum_change_ns, linear, 1e-29);
    close_vec(c.report.audit.angular_momentum_change_nms, angular, 1e-34);
    assert!(norm(c.report.audit.linear_momentum_residual_ns) < 1e-29);
    assert!(norm(c.report.audit.angular_momentum_residual_nms) < 1e-34);
    assert_eq!(c.report.external_wrenches, wrenches);
    close(c.report.external_work_j, external_work, 1e-28);
    let delta_k = kinetic(c.bodies_after, f.properties) - kinetic(f.bodies, f.properties);
    close(c.report.kinetic_energy_change_j, delta_k, 2e-28);
    close(
        delta_k,
        c.report.midpoint_impulse_work_j + external_work,
        2e-26,
    );
    close(
        delta_k + c.report.normal_storage_change_j + c.report.dissipated_energy_j,
        external_work,
        2e-26,
    );
    for residual in c.report.audit.body_work_energy_residual_j {
        assert!(residual.abs() <= controls().energy_tolerance_j);
    }
}

#[test]
fn zero_wrench_entry_points_preserve_original_gap_and_sliding_phases() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    for f in [
        axial(-0.5 * common::Z),
        crossing(&normal, -0.2 * common::Z, 0.2 * common::Z, dt).0,
    ] {
        let gate = CancelGate::new();
        let contact = f.contact(&normal, &shear);
        let mut prepared =
            PreparedRigidDmtPair::new(&f.bodies, &f.properties, contact, dt, controls(), &gate)
                .unwrap();
        let start = prepared.initialize([0.0; 2]).unwrap();
        let expected = prepared.advance(start, &gate).unwrap();
        assert_eq!(
            prepared
                .advance_with_wrenches(start, [Wrench::ZERO; 2], &gate)
                .unwrap(),
            expected
        );
        assert_eq!(
            prepare_rigid_dmt_contact_with_wrenches(
                &f.bodies,
                &f.properties,
                contact,
                start,
                [Wrench::ZERO; 2],
                dt,
                controls(),
                &gate,
            )
            .unwrap(),
            prepare_rigid_dmt_contact(
                &f.bodies,
                &f.properties,
                contact,
                start,
                dt,
                controls(),
                &gate,
            )
            .unwrap()
        );
        let mut legacy_bodies = f.bodies;
        let mut driven_bodies = f.bodies;
        let mut legacy_history = start;
        let mut driven_history = start;
        let legacy = apply_rigid_dmt_contact(
            &mut legacy_bodies,
            &f.properties,
            contact,
            &mut legacy_history,
            dt,
            controls(),
            &gate,
        )
        .unwrap();
        let driven = apply_rigid_dmt_contact_with_wrenches(
            &mut driven_bodies,
            &f.properties,
            contact,
            &mut driven_history,
            [Wrench::ZERO; 2],
            dt,
            controls(),
            &gate,
        )
        .unwrap();
        assert_eq!(driven_bodies, legacy_bodies);
        assert_eq!(driven_history, legacy_history);
        assert_eq!(driven, legacy);
        assert_eq!(driven.external_work_j, 0.0);
        check_driven(&f, &expected, [Wrench::ZERO; 2], dt);
    }
}

#[test]
fn driven_attractive_gap_matches_exact_tail_and_has_no_solid_shear_or_heat() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    for sign in [-1.0, 1.0] {
        let mut f = axial(-0.5 * common::Z);
        f.set_a_velocity(Vec3::new(0.0, 0.0, 0.02), 0.0);
        let wrenches = [
            Wrench {
                force_world: Vec3::new(2e-8, 0.0, sign * 1e-8),
                torque_body: Vec3::new(0.0, 0.0, 1e-18),
            },
            Wrench::ZERO,
        ];
        let gate = CancelGate::new();
        let mut prepared = PreparedRigidDmtPair::new(
            &f.bodies,
            &f.properties,
            f.contact(&normal, &shear),
            dt,
            controls(),
            &gate,
        )
        .unwrap();
        let start = prepared.initialize([0.0; 2]).unwrap();
        let candidate = prepared
            .advance_with_wrenches(start, wrenches, &gate)
            .unwrap();
        let x0 = prepared.approach_m();
        let w = 1.0 / f.properties[0].mass() + 1.0 / f.properties[1].mass();
        let compliance = 0.5 * dt * dt * w;
        let k = 2.0 * PI * normal.effective_radius_m() * ADHESION / common::Z;
        let acceleration = wrenches[0].force_world.z / f.properties[0].mass();
        let dx = (dt * 0.02 + 0.5 * dt * dt * acceleration + compliance * k * (common::Z + x0))
            / (1.0 - 0.5 * compliance * k);
        let force = -k * (common::Z + x0 + 0.5 * dx);
        close(candidate.state_after.position_m()[0], x0 + dx, 2e-21);
        close(candidate.report.applied_mean_force_n[0], force, 1e-18);
        close(
            candidate.report.relative_velocity_m_s[1][0],
            0.02 + dt * acceleration - dt * w * force,
            1e-12,
        );
        close(candidate.report.normal_storage_change_j, force * dx, 1e-28);
        assert!(force < 0.0);
        assert!(candidate.state_after.position_m()[0] < 0.0);
        assert_eq!(candidate.report.contact.path.mean_real_area_m2, 0.0);
        assert_eq!(candidate.report.applied_mean_force_n[1..], [0.0; 2]);
        assert_eq!(candidate.report.dissipated_energy_j, 0.0);
        assert_eq!(candidate.state_after.regime(), DmtCoupledRegime::Released);
        assert!(candidate.report.relative_spin_rad_s[1] > 0.0);
        assert!(sign * candidate.report.external_work_j > 0.0);
        check_driven(&f, &candidate, wrenches, dt);
    }
}

#[test]
fn driven_crossing_retains_all_cross_mobility_area_heat_and_original_body_order() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    let to = 0.2 * common::Z;
    let (mut f, expected_force, area) = crossing(&normal, -0.2 * common::Z, to, dt);
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, -0.5e-8, 0.7e-8),
            torque_body: Vec3::new(0.0, 2e12 * f.properties[0].principal_inertia_body().y, 0.0),
        },
        Wrench {
            force_world: Vec3::new(-0.3e-8, 0.2e-8, -0.6e-8),
            torque_body: Vec3::new(0.0, -1e12 * f.properties[1].principal_inertia_body().y, 0.0),
        },
    ];
    // Preserve the independent manufactured midpoint path, now with nonzero
    // external acceleration. Body torques cause normal AND tangential drive.
    let acceleration = external_point_acceleration(&f, wrenches);
    assert!(acceleration.x.abs() > 1e5 && acceleration.z.abs() > 1e5);
    let initial_velocity = f.bodies[0]
        .linear_momentum_world()
        .scale(1.0 / f.properties[0].mass());
    f.set_a_velocity(initial_velocity.sub(acceleration.scale(0.5 * dt)), 0.0);
    let contact = f.contact(&normal, &shear);
    let gate = CancelGate::new();
    let mut prepared =
        PreparedRigidDmtPair::new(&f.bodies, &f.properties, contact, dt, controls(), &gate)
            .unwrap();
    let start = prepared.initialize([0.0; 2]).unwrap();
    let candidate = prepared
        .advance_with_wrenches(start, wrenches, &gate)
        .unwrap();
    let w = mobility(&f);
    assert!(w[0][1].abs() > 1e12 && w[0][2].abs() > 1e12 && w[1][2].abs() > 1e12);
    let axes = [f.frame.normal(), f.frame.longitudinal(), f.frame.lateral()].map(vector);
    for i in 0..3 {
        close(
            candidate.report.applied_mean_force_n[i],
            expected_force[i],
            1e-18,
        );
        for j in 0..3 {
            close(candidate.report.mobility_per_kg[i][j], w[i][j], 2.0);
        }
        let response: f64 = (0..3).map(|j| w[i][j] * expected_force[j]).sum();
        close(
            candidate.report.relative_displacement_m[i],
            dt * candidate.report.relative_velocity_m_s[0][i]
                + 0.5 * dt * dt * (acceleration.dot(axes[i]) - response),
            2e-21,
        );
    }
    close(candidate.state_after.position_m()[0], to, 2e-21);
    close(candidate.report.contact.path.mean_real_area_m2, area, 1e-27);
    close(
        candidate.report.dissipated_energy_j,
        KINETIC * area * 0.1 * common::Z,
        1e-28,
    );
    assert!(candidate.report.dissipated_energy_j > 0.0);
    assert_eq!(candidate.report.contact.regime, DmtCoupledRegime::Sliding);
    close(candidate.report.relative_spin_rad_s[1], 0.0, 1e-6);
    check_driven(&f, &candidate, wrenches, dt);

    let mut bodies = [f.bodies[1], f.bodies[0]];
    let properties = [f.properties[1], f.properties[0]];
    let permuted = RigidDmtContact {
        binding: RigidContactBinding {
            body_a: 1,
            body_b: 0,
            ..contact.binding
        },
        ..contact
    };
    let mut history = start;
    let report = apply_rigid_dmt_contact_with_wrenches(
        &mut bodies,
        &properties,
        permuted,
        &mut history,
        wrenches,
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(
        bodies,
        [candidate.bodies_after[1], candidate.bodies_after[0]]
    );
    assert_eq!(history, candidate.state_after);
    assert_eq!(report, candidate.report);
}

#[test]
fn full_inertia_and_world_rotation_preserve_body_torque_work_and_native_gap_motion() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    let approach = -0.5 * common::Z;
    let mut f = fixture(
        approach,
        [
            Vec3::new(0.2, 0.1, 0.4).scale(LENGTH),
            Vec3::new(-0.12, 0.15, -0.3).scale(LENGTH),
        ],
    );
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.6).unwrap();
    f.properties[0] = MassProperties::from_principal_frame(
        f.properties[0].mass(),
        Vec3::ZERO,
        f.properties[0].principal_inertia_body(),
        principal,
    )
    .unwrap();
    let orientations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.3).unwrap(),
    ];
    for i in 0..2 {
        f.bodies[i] = RigidBodyState::new(
            Pose::new(
                orientations[i]
                    .rotate_body_to_world(f.reaction_arms[i])
                    .scale(-1.0),
                orientations[i],
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap();
    }
    f.binding.arm_a_body = f.reaction_arms[0].add(orientations[0].rotate_world_to_body(Vec3::new(
        0.0,
        0.0,
        0.5 * approach,
    )));
    f.binding.arm_b_body = f.reaction_arms[1].sub(orientations[1].rotate_world_to_body(Vec3::new(
        0.0,
        0.0,
        0.5 * approach,
    )));
    f.set_a_velocity(Vec3::new(0.01, 0.0, 0.02), 0.0);
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.2, -0.1, 0.3).scale(1e-8),
            torque_body: Vec3::new(0.4, 0.2, -0.3).scale(1e-15),
        },
        Wrench {
            force_world: Vec3::new(-0.3, 0.4, -0.1).scale(1e-8),
            torque_body: Vec3::new(-0.2, 0.1, 0.4).scale(1e-15),
        },
    ];
    let gate = CancelGate::new();
    let mut prepared = PreparedRigidDmtPair::new(
        &f.bodies,
        &f.properties,
        f.contact(&normal, &shear),
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    let start = prepared.initialize([0.0; 2]).unwrap();
    let candidate = prepared
        .advance_with_wrenches(start, wrenches, &gate)
        .unwrap();
    let w = full_mobility(&f);
    let acceleration = external_point_acceleration(&f, wrenches);
    let x0 = prepared.approach_m();
    let c = 0.5 * dt * dt * w[0][0];
    let k = 2.0 * PI * normal.effective_radius_m() * ADHESION / common::Z;
    let dx = (dt * 0.02 + 0.5 * dt * dt * acceleration.z + c * k * (common::Z + x0))
        / (1.0 - 0.5 * c * k);
    close(candidate.state_after.position_m()[0], x0 + dx, 2e-21);
    for i in 0..3 {
        for j in 0..3 {
            close(candidate.report.mobility_per_kg[i][j], w[i][j], 2.0);
        }
    }
    assert_eq!(candidate.report.contact.path.mean_real_area_m2, 0.0);
    check_driven(&f, &candidate, wrenches, dt);

    let rotation_vector = Vec3::new(1.0, 1.0, -1.0).scale(0.7 / 3.0_f64.sqrt());
    let rotation = UnitQuaternion::IDENTITY.right_exp(rotation_vector).unwrap();
    let mut rotated = fixture(approach, f.reaction_arms);
    rotated.properties = f.properties;
    rotated.binding = f.binding;
    for i in 0..2 {
        let pose = f.bodies[i].pose();
        let orientation = pose
            .orientation()
            .right_exp(pose.orientation().rotate_world_to_body(rotation_vector))
            .unwrap();
        rotated.bodies[i] = RigidBodyState::new(
            Pose::new(
                orientation
                    .rotate_body_to_world(f.reaction_arms[i])
                    .scale(-1.0),
                orientation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(f.bodies[i].linear_momentum_world()),
            f.bodies[i].angular_momentum_body(),
        )
        .unwrap();
    }
    let n = rotation.rotate_body_to_world(vector(f.frame.normal()));
    let t = rotation.rotate_body_to_world(vector(f.frame.longitudinal()));
    rotated.frame = TangentFrame::new([n.x, n.y, n.z], [t.x, t.y, t.z]).unwrap();
    let rotated_wrenches = wrenches.map(|w| Wrench {
        force_world: rotation.rotate_body_to_world(w.force_world),
        torque_body: w.torque_body,
    });
    let mut pair = PreparedRigidDmtPair::new(
        &rotated.bodies,
        &rotated.properties,
        rotated.contact(&normal, &shear),
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    let state = pair.initialize([0.0; 2]).unwrap();
    let actual = pair
        .advance_with_wrenches(state, rotated_wrenches, &gate)
        .unwrap();
    check_driven(&rotated, &actual, rotated_wrenches, dt);
    close(
        actual.report.external_work_j,
        candidate.report.external_work_j,
        1e-28,
    );
    close(
        actual.report.normal_storage_change_j,
        candidate.report.normal_storage_change_j,
        1e-28,
    );
    close_vec(
        actual.report.force_impulse_on_a_world_ns,
        rotation.rotate_body_to_world(candidate.report.force_impulse_on_a_world_ns),
        1e-27,
    );
    for i in 0..2 {
        close_vec(
            actual.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(candidate.bodies_after[i].linear_momentum_world()),
            1e-27,
        );
        close_vec(
            actual.bodies_after[i].angular_momentum_body(),
            candidate.bodies_after[i].angular_momentum_body(),
            1e-33,
        );
        close_vec(
            actual.report.angular_displacement_body_rad[i],
            candidate.report.angular_displacement_body_rad[i],
            1e-15,
        );
    }
}

#[test]
fn driven_cancellation_source_refusal_and_late_torsion_leave_owners_retryable() {
    let (normal, shear) = common::sources(ADHESION, 3e6, KINETIC);
    let dt = 5e-10;
    let f = axial(0.2 * common::Z);
    let contact = f.contact(&normal, &shear);
    let gate = CancelGate::new();
    let start = PreparedRigidDmtPair::new(&f.bodies, &f.properties, contact, dt, controls(), &gate)
        .unwrap()
        .initialize([0.0; 2])
        .unwrap();
    let valid = [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, 1e-8),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    let mut bodies = f.bodies;
    let mut history = start;
    for force in [f64::NAN, 1.0] {
        let mut invalid = valid;
        invalid[0].force_world.z = force;
        assert!(
            apply_rigid_dmt_contact_with_wrenches(
                &mut bodies,
                &f.properties,
                contact,
                &mut history,
                invalid,
                dt,
                controls(),
                &gate,
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
        assert_eq!(history, start);
    }
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_dmt_contact_with_wrenches(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            valid,
            dt,
            controls(),
            &cancelled,
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, start);
    let mut torsion = valid;
    torsion[0].torque_body.z = 1e-18;
    // This refusal follows the original source solve and actual native recoil.
    assert!(matches!(
        apply_rigid_dmt_contact_with_wrenches(
            &mut bodies,
            &f.properties,
            contact,
            &mut history,
            torsion,
            dt,
            controls(),
            &gate,
        ),
        Err(RigidDmtError::UnsupportedTorsion { .. })
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, start);
    let expected = prepare_rigid_dmt_contact_with_wrenches(
        &f.bodies,
        &f.properties,
        contact,
        start,
        valid,
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    let report = apply_rigid_dmt_contact_with_wrenches(
        &mut bodies,
        &f.properties,
        contact,
        &mut history,
        valid,
        dt,
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(history, expected.state_after);
    assert_eq!(report, expected.report);
    check_driven(&f, &expected, valid, dt);
}
