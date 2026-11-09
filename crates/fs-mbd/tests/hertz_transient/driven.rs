//! Original-body powered shear/spin, common motion, work and atomic refusal.
use super::*;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::hertz::driven::*;

fn driven_controls() -> HertzDrivenControls {
    HertzDrivenControls {
        maximum_step_s: 2e-5,
        minimum_step_s: 1e-14,
        maximum_attempts: 20_000,
        maximum_evaluations: 200_000,
        scaled_rate_tolerance_m_s: 2e-11,
        rest_scaled_rate_tolerance_m_s: 2e-11,
    }
}
fn centred(common_velocity: f64) -> Fixture {
    let mut f = fixture();
    f.binding.arm_a_body = Vec3::ZERO;
    f.binding.arm_b_body = Vec3::ZERO;
    let pose = Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap();
    f.bodies = f.properties.map(|p| {
        RigidBodyState::new(
            pose,
            Vec3::new(p.mass() * common_velocity, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    });
    f
}
fn propose(
    f: &Fixture,
    law: &HertzWrenchContact,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    dt: f64,
) -> RigidHertzDrivenCandidate {
    PreparedRigidHertzDriven::new(
        &f.bodies,
        &f.properties,
        f.contact(law),
        external,
        mode,
        dt,
        controls(),
        driven_controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(law, &CancelGate::new())
    .unwrap()
}

#[test]
fn static_contact_transmits_drive_while_common_motion_receives_work_without_heat() {
    let f = centred(1.25);
    let external = [
        Wrench {
            force_world: Vec3::new(3.0, 0.0, 0.0),
            ..Wrench::ZERO
        },
        Wrench::ZERO,
    ];
    let dt = 0.4;
    let candidate = propose(&f, &law(), external, HertzDrivenMode::Resting, dt);
    let report = &candidate.report;
    assert_eq!(report.driven.mode, HertzDrivenMode::Resting);
    close(report.dissipated_energy_j, 0.0, 1e-15);
    close(report.driven.slip_increment_m[0], 0.0, 1e-15);
    let acceleration = 3.0 / 5.0;
    for (body, p) in candidate.bodies_after.iter().zip(f.properties) {
        close(
            body.linear_momentum_world().x / p.mass(),
            1.25 + acceleration * dt,
            2e-14,
        );
    }
    let work = 3.0 * (1.25 * dt + 0.5 * acceleration * dt * dt);
    close(report.external_work_j, work, 2e-14);
    close(report.kinetic_change_j, work, 2e-14);
    close(report.force_impulse_world_ns.x, -1.8 * dt, 2e-13);
    close(report.audit.total_energy_residual_j, 0.0, 3e-14);
    assert!(report.driven.accepted_steps <= 2);
}

#[test]
fn motor_breakaway_uses_original_hertz_spin_capacity_and_actual_both_body_work() {
    let f = centred(0.0);
    let torque = 0.3;
    let external = [
        Wrench {
            torque_body: Vec3::new(0.0, 0.0, torque),
            ..Wrench::ZERO
        },
        Wrench::ZERO,
    ];
    let dt = 0.1;
    let candidate = propose(&f, &law(), external, HertzDrivenMode::Resting, dt);
    let report = &candidate.report;
    let moment = 3.0 * core::f64::consts::PI / 16.0 * 4.0 * R;
    let alpha_a = (torque - moment) / f.properties[0].principal_inertia_body().z;
    let alpha_b = moment / f.properties[1].principal_inertia_body().z;
    assert_eq!(report.driven.mode, HertzDrivenMode::Sliding);
    close(report.spin_rad_s[1], (alpha_a - alpha_b) * dt, 3e-9);
    close(
        report.driven.spin_increment_rad,
        0.5 * (alpha_a - alpha_b) * dt * dt,
        3e-10,
    );
    close(
        report.dissipated_energy_j,
        moment * 0.5 * (alpha_a - alpha_b) * dt * dt,
        3e-11,
    );
    close(
        report.external_work_j,
        torque * 0.5 * alpha_a * dt * dt,
        3e-11,
    );
    for (i, alpha) in [alpha_a, alpha_b].into_iter().enumerate() {
        close(
            candidate.bodies_after[i].angular_momentum_body().z,
            f.properties[i].principal_inertia_body().z * alpha * dt,
            3e-12,
        );
    }
    close(report.audit.total_energy_residual_j, 0.0, 2e-10);
    let whole = propose(&f, &law(), external, HertzDrivenMode::Resting, 0.23);
    let mut continued = centred(0.0);
    continued.bodies = candidate.bodies_after;
    let remainder = propose(&continued, &law(), external, report.driven.mode, 0.13);
    for i in 0..2 {
        close_vec(
            remainder.bodies_after[i].angular_momentum_body(),
            whole.bodies_after[i].angular_momentum_body(),
            3e-12,
        );
    }
    close(
        report.external_work_j + remainder.report.external_work_j,
        whole.report.external_work_j,
        3e-10,
    );
    close(
        report.dissipated_energy_j + remainder.report.dissipated_energy_j,
        whole.report.dissipated_energy_j,
        3e-10,
    );
}

// Independent fixed-grid RK4 integrates each actual body's COM rate and spin,
// physical Hertz power, external power and relative travel. No reduced source-
// work formula or production contact-mobility solve enters this reference.
fn driven_reference(f: &Fixture, external: [Wrench; 2], dt: f64) -> [f64; 11] {
    let patch = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: R,
        friction_coefficient: 0.4,
    })
    .unwrap();
    let arms = [f.binding.arm_a_body, f.binding.arm_b_body];
    let rhs = |y: [f64; 11]| {
        let v = [Vec3::new(y[0], y[1], 0.0), Vec3::new(y[3], y[4], 0.0)];
        let omega = [Vec3::new(0.0, 0.0, y[2]), Vec3::new(0.0, 0.0, y[5])];
        let relative = v[0]
            .add(omega[0].cross(arms[0]))
            .sub(v[1].add(omega[1].cross(arms[1])));
        let spin = y[2] - y[5];
        let contact = patch.evaluate([relative.x, relative.y], spin).unwrap();
        let force = Vec3::new(contact.force_n[0], contact.force_n[1], 0.0);
        let mut result = [0.0; 11];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let acceleration = external[i]
                .force_world
                .add(force.scale(sign))
                .scale(f.properties[i].mass().recip());
            result[3 * i] = acceleration.x;
            result[3 * i + 1] = acceleration.y;
            result[3 * i + 2] = (external[i].torque_body.z
                + sign * (arms[i].cross(force).z + contact.spin_moment_nm))
                / f.properties[i].principal_inertia_body().z;
            result[7] += external[i].force_world.dot(v[i]) + external[i].torque_body.dot(omega[i]);
        }
        result[6] = contact.dissipated_power_w;
        result[8] = relative.x;
        result[9] = relative.y;
        result[10] = spin;
        result
    };
    let mut y = [0.0; 11];
    for i in 0..2 {
        y[3 * i] = f.bodies[i].linear_momentum_world().x / f.properties[i].mass();
        y[3 * i + 1] = f.bodies[i].linear_momentum_world().y / f.properties[i].mass();
        y[3 * i + 2] =
            f.bodies[i].angular_momentum_body().z / f.properties[i].principal_inertia_body().z;
    }
    let h = dt / 4096.0;
    for _ in 0..4096 {
        let a = rhs(y);
        let b = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn full_coupled_driven_motion_matches_independent_body_ode_and_external_power() {
    let f = fixture();
    let external = [
        Wrench {
            force_world: Vec3::new(1.2, 0.5, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.004),
        },
        Wrench {
            force_world: Vec3::new(-0.3, 0.7, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.002),
        },
    ];
    let dt = 0.006;
    let result = propose(&f, &law(), external, HertzDrivenMode::Sliding, dt);
    let expected = driven_reference(&f, external, dt);
    for i in 0..2 {
        let body = result.bodies_after[i];
        close(
            body.linear_momentum_world().x / f.properties[i].mass(),
            expected[3 * i],
            2e-5,
        );
        close(
            body.linear_momentum_world().y / f.properties[i].mass(),
            expected[3 * i + 1],
            2e-5,
        );
        close(
            R * body.angular_momentum_body().z / f.properties[i].principal_inertia_body().z,
            R * expected[3 * i + 2],
            2e-5,
        );
    }
    close(result.report.dissipated_energy_j, expected[6], 2e-6);
    close(result.report.external_work_j, expected[7], 2e-7);
    close(result.report.driven.slip_increment_m[0], expected[8], 2e-7);
    close(result.report.driven.slip_increment_m[1], expected[9], 2e-7);
    close(
        R * result.report.driven.spin_increment_rad,
        R * expected[10],
        2e-7,
    );
    close(result.report.audit.total_energy_residual_j, 0.0, 2e-9);

    // World forces rotate with the world frame; body torques retain their
    // original body coordinates. The driven response must transform together.
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.7).unwrap();
    let array = |v: Vec3| [v.x, v.y, v.z];
    let mut rotated = fixture();
    rotated.bodies = f.bodies.map(|body| {
        RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(body.pose().position_world()),
                rotation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(body.linear_momentum_world()),
            body.angular_momentum_body(),
        )
        .unwrap()
    });
    rotated.frame = TangentFrame::new(
        array(rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
        array(rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
    )
    .unwrap();
    let rotated_external = external.map(|w| Wrench {
        force_world: rotation.rotate_body_to_world(w.force_world),
        ..w
    });
    let transformed = propose(
        &rotated,
        &law(),
        rotated_external,
        HertzDrivenMode::Sliding,
        dt,
    );
    for i in 0..2 {
        close_vec(
            transformed.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(result.bodies_after[i].linear_momentum_world()),
            3e-10,
        );
        close_vec(
            transformed.bodies_after[i].angular_momentum_body(),
            result.bodies_after[i].angular_momentum_body(),
            3e-11,
        );
    }
    close(
        transformed.report.external_work_j,
        result.report.external_work_j,
        3e-11,
    );
    close(
        transformed.report.dissipated_energy_j,
        result.report.dissipated_energy_j,
        3e-11,
    );
}

#[test]
fn unrepresentable_common_drive_remains_external_work_and_refuses_without_publication() {
    let f = centred(2.0_f64.powi(54));
    let external = f.properties.map(|p| Wrench {
        force_world: Vec3::new(0.25 * p.mass(), 0.0, 0.0),
        ..Wrench::ZERO
    });
    let mut bodies = f.bodies;
    let mut mode = HertzDrivenMode::Resting;
    let error = apply_rigid_hertz_driven(
        &mut bodies,
        &f.properties,
        f.contact(&law()),
        external,
        &mut mode,
        0.5,
        controls(),
        driven_controls(),
        &CancelGate::new(),
    )
    .unwrap_err();
    let RigidHertzDrivenError::Audit(audit) = error else {
        panic!("{error:?}")
    };
    // The momenta cannot represent the commanded common acceleration. Its
    // independently integrated source work must still appear in the defect.
    assert!(audit.total_energy_residual_j.abs() > 1e15);
    assert_eq!(bodies, f.bodies);
    assert_eq!(mode, HertzDrivenMode::Resting);
}

#[test]
fn mode_and_bodies_publish_together_preserve_spectators_and_refuse_opening() {
    let f = centred(0.0);
    let law = law();
    let external = [
        Wrench {
            force_world: Vec3::new(20.0, 0.0, 0.0),
            ..Wrench::ZERO
        },
        Wrench::ZERO,
    ];
    let mut bodies = vec![f.bodies[0], f.bodies[1], f.bodies[1]];
    let properties = vec![f.properties[0], f.properties[1], f.properties[1]];
    let before = bodies.clone();
    let mut mode = HertzDrivenMode::Resting;
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_hertz_driven(
            &mut bodies,
            &properties,
            f.contact(&law),
            external,
            &mut mode,
            0.1,
            controls(),
            driven_controls(),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(bodies, before);
    assert_eq!(mode, HertzDrivenMode::Resting);
    let opening = [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, 1.0),
            ..Wrench::ZERO
        },
        Wrench::ZERO,
    ];
    assert!(
        apply_rigid_hertz_driven(
            &mut bodies,
            &properties,
            f.contact(&law),
            opening,
            &mut mode,
            0.1,
            controls(),
            driven_controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, before);
    assert_eq!(mode, HertzDrivenMode::Resting);
    let report = apply_rigid_hertz_driven(
        &mut bodies,
        &properties,
        f.contact(&law),
        external,
        &mut mode,
        0.1,
        controls(),
        driven_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, HertzDrivenMode::Sliding);
    assert_eq!(report.driven.mode, mode);
    assert_eq!(bodies[2], before[2]);
    assert_ne!(bodies[..2], before[..2]);
}
