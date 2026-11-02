use super::*;
use fs_mbd::Wrench;

fn candidate(f: &Fixture, external: [Wrench; 2], dt: f64) -> RigidHertzRollingCandidate {
    PreparedRigidHertzRollingPair::new_driven(
        &f.bodies,
        &f.properties,
        &RigidRollingHistory::new(),
        f.contact(),
        external,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.hertz, &CancelGate::new())
    .unwrap()
}

fn drive() -> [Wrench; 2] {
    [
        Wrench {
            force_world: Vec3::new(0.7, -0.5, 0.0),
            torque_body: Vec3::new(0.003, -0.004, 0.08),
        },
        Wrench {
            force_world: Vec3::new(-0.2, 0.15, 0.0),
            torque_body: Vec3::new(-0.002, 0.001, -0.015),
        },
    ]
}

fn audit_drive(f: &Fixture, c: &RigidHertzRollingCandidate, external: [Wrench; 2], dt: f64) {
    let old = points(f, &f.bodies);
    let new = points(f, &c.bodies_after);
    let r = &c.report;
    let mut external_work = 0.0;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let point_impulse = r.force_impulse_world_ns.scale(sign);
        let body_couple = f.bodies[i].pose().orientation().rotate_world_to_body(
            r.spin_angular_impulse_world_nms
                .add(r.rolling_angular_impulse_world_nms)
                .scale(sign),
        );
        let arm = [f.binding.arm_a_body, f.binding.arm_b_body][i];
        let expected_linear = point_impulse.add(external[i].force_world.scale(dt));
        let expected_angular = arm
            .cross(
                f.bodies[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(point_impulse),
            )
            .add(body_couple)
            .add(external[i].torque_body.scale(dt));
        close_vec(
            c.bodies_after[i]
                .linear_momentum_world()
                .sub(f.bodies[i].linear_momentum_world()),
            expected_linear,
            3e-16,
        );
        close_vec(
            c.bodies_after[i]
                .angular_momentum_body()
                .sub(f.bodies[i].angular_momentum_body()),
            expected_angular,
            3e-17,
        );
        external_work += 0.5
            * dt
            * (external[i].force_world.dot(
                old[i]
                    .center_of_mass_velocity_world
                    .add(new[i].center_of_mass_velocity_world),
            ) + external[i].torque_body.dot(
                old[i]
                    .angular_velocity_body
                    .add(new[i].angular_velocity_body),
            ));
        assert_eq!(c.bodies_after[i].pose(), f.bodies[i].pose());
    }
    close(r.external_wrench_work_j, external_work, 3e-16);
    close(
        kinetic(&c.bodies_after, &f.properties) - kinetic(&f.bodies, &f.properties)
            + r.surface_heat_j
            + r.bulk_heat_j,
        external_work,
        3e-12,
    );
    close_vec(r.linear_momentum_residual_ns, Vec3::ZERO, 3e-14);
    close_vec(r.angular_momentum_residual_nms, Vec3::ZERO, 3e-14);
}

#[test]
fn motor_driven_mixed_shear_spin_and_rolling_match_pressure_integrals_and_native_work() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let c = candidate(&f, drive(), DT);
    audit_drive(&f, &c, drive(), DT);
    assert!(
        c.report
            .scaled_driving_acceleration_m_s2
            .iter()
            .all(|a| a.abs() > 1e-5)
    );
    assert!(c.report.kinetic_after_j > c.report.kinetic_before_j);
    let old = points(&f, &f.bodies);
    let new = points(&f, &c.bodies_after);
    let u = relative_point(old).add(relative_point(new)).scale(0.5);
    let omega = relative_omega(old).add(relative_omega(new)).scale(0.5);
    let (force, moment, heat) = surface_reference(
        [u.x, u.y],
        omega.z,
        f.patch.normal_load_n(),
        f.radius,
        0.2,
        DT,
        256,
    );
    close(c.report.hertz.force_n[0], force[0], 3e-8);
    close(c.report.hertz.force_n[1], force[1], 3e-8);
    close(c.report.hertz.spin_moment_nm, moment, 3e-10);
    close(c.report.surface_heat_j, heat, 3e-12);
    let damping = TAU * f.patch.normal_load_n() * R;
    close_vec(
        c.report.rolling_angular_impulse_world_nms,
        Vec3::new(omega.x, omega.y, 0.0).scale(-damping * DT),
        3e-17,
    );
    close(
        c.report.bulk_heat_j,
        damping * DT * (omega.x * omega.x + omega.y * omega.y),
        3e-17,
    );
    assert!(c.report.surface_heat_j > 0.0 && c.report.bulk_heat_j > 0.0);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let report = apply_rigid_hertz_rolling_driven(
        &mut bodies,
        &f.properties,
        &mut history,
        f.contact(),
        drive(),
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, c.bodies_after);
    assert_eq!(history, c.rolling_history);
    assert_eq!(report, c.report);
}

#[test]
fn offcentre_full_inertia_driven_response_is_covariant_in_world_and_tangent_frames() {
    let mut f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    f.properties[1] = MassProperties::from_principal_frame(
        7.0,
        Vec3::ZERO,
        Vec3::new(0.08, 0.13, 0.18),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 0.5), 0.63).unwrap(),
    )
    .unwrap();
    let expected = candidate(&f, drive(), DT);
    audit_drive(&f, &expected, drive(), DT);
    assert!(expected.report.scaled_mobility_per_kg[2][3].abs() > 1e-5);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.73).unwrap();
    let external = drive().map(|w| Wrench {
        force_world: rotation.rotate_body_to_world(w.force_world),
        ..w
    });
    f.bodies = f.bodies.map(|body| rotate_body(body, rotation));
    f.frame = TangentFrame::new(
        array(rotation.rotate_body_to_world(vector(f.frame.normal()))),
        array(rotation.rotate_body_to_world(vector(f.frame.longitudinal()))),
    )
    .unwrap()
    .rotated(0.41)
    .unwrap();
    let actual = candidate(&f, external, DT);
    audit_drive(&f, &actual, external, DT);
    for i in 0..2 {
        let transformed = rotate_body(expected.bodies_after[i], rotation);
        close_vec(
            actual.bodies_after[i].linear_momentum_world(),
            transformed.linear_momentum_world(),
            3e-11,
        );
        close_vec(
            actual.bodies_after[i].angular_momentum_body(),
            transformed.angular_momentum_body(),
            3e-12,
        );
    }
    close(
        actual.report.surface_heat_j,
        expected.report.surface_heat_j,
        3e-12,
    );
    close(
        actual.report.bulk_heat_j,
        expected.report.bulk_heat_j,
        3e-17,
    );
    close(
        actual.report.external_wrench_work_j,
        expected.report.external_wrench_work_j,
        3e-15,
    );
}

#[test]
fn zero_drive_is_exact_and_common_acceleration_keeps_relative_friction_and_heat() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let undriven = PreparedRigidHertzRollingPair::new(
        &f.bodies,
        &f.properties,
        &RigidRollingHistory::new(),
        f.contact(),
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.hertz, &CancelGate::new())
    .unwrap();
    assert_eq!(candidate(&f, [Wrench::ZERO; 2], DT), undriven);
    let acceleration = Vec3::new(0.3, -0.1, 0.2);
    let external = f.properties.map(|p| Wrench {
        force_world: acceleration.scale(p.mass()),
        torque_body: Vec3::ZERO,
    });
    let driven = candidate(&f, external, DT);
    audit_drive(&f, &driven, external, DT);
    close(
        driven.report.surface_heat_j,
        undriven.report.surface_heat_j,
        3e-16,
    );
    close(
        driven.report.bulk_heat_j,
        undriven.report.bulk_heat_j,
        3e-17,
    );
    for i in 0..2 {
        close_vec(
            driven.bodies_after[i]
                .linear_momentum_world()
                .sub(undriven.bodies_after[i].linear_momentum_world()),
            external[i].force_world.scale(DT),
            3e-16,
        );
        close_vec(
            driven.bodies_after[i].angular_momentum_body(),
            undriven.bodies_after[i].angular_momentum_body(),
            3e-17,
        );
    }
}

#[test]
fn driven_endpoint_reversal_normal_motion_and_cancellation_preserve_bodies_and_history() {
    let f = fixture(TAU, MIXED_SLIP, MIXED_ROLL, MIXED_SPIN, [0.3, 0.2]);
    let mut bodies = f.bodies;
    let mut history = RigidRollingHistory::new();
    let normal_drive = [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, 1.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    assert!(
        apply_rigid_hertz_rolling_driven(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            normal_drive,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, RigidRollingHistory::new());
    // Zero material viscosity keeps this deliberately coarse reversal inside
    // the original rolling applicability envelope, isolating endpoint work.
    let reversal = fixture(0.0, [0.01, 0.0], [0.0; 2], 0.0, [0.6, 0.2]);
    let mut reversing_bodies = reversal.bodies;
    let braking = [
        Wrench {
            force_world: Vec3::new(-100.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    assert!(matches!(
        apply_rigid_hertz_rolling_driven(
            &mut reversing_bodies,
            &reversal.properties,
            &mut history,
            reversal.contact(),
            braking,
            0.01,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidHertzRollingError::UnresolvedSlipEvent { .. })
    ));
    assert_eq!(reversing_bodies, reversal.bodies);
    assert_eq!(history, RigidRollingHistory::new());
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_hertz_rolling_driven(
            &mut bodies,
            &f.properties,
            &mut history,
            f.contact(),
            drive(),
            DT,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(history, RigidRollingHistory::new());
}
