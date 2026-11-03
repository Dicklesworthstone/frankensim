//! Actual contact motion must survive an unrelated inertial-frame translation.
use super::*;

fn translated_scene(
    boost: f64,
    origin: Vec3,
    normal_speed: f64,
) -> ([RigidBodyState; 1], RigidFrictionalContact) {
    let mut contact = plane(0.0);
    contact.friction_coefficient = 0.0;
    contact.other = RigidImpactCounterpart::PrescribedPlane {
        point_world: origin,
        velocity_world: Vec3::new(boost, 0.0, 0.0),
    };
    let state = RigidBodyState::new(
        Pose::new(
            block().pose().position_world().add(origin),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        Vec3::new(boost, 0.0, -normal_speed),
        Vec3::ZERO,
    )
    .unwrap();
    ([state], contact)
}

#[test]
fn plastic_contact_keeps_small_recoil_at_large_translation_and_world_origin() {
    // 300 m/s at an Earth-scale coordinate with gentle normal approach, plus
    // a deliberately extreme boost. No normal/tangent geometry is rounded away:
    // the shared moving support makes these the same physical normal impact.
    for (boost, origin, normal_speed) in [
        (0.0, Vec3::ZERO, 1e-3),
        (300.0, Vec3::new(0.37, 0.0, 6.4e6), 1e-3),
        (1e8, Vec3::ZERO, 1.0),
    ] {
        let (mut bodies, contact) = translated_scene(boost, origin, normal_speed);
        let initial = bodies[0];
        let dt = 1.0 / 32.0;
        let result = step_rigid_frictional_contacts(
            &mut bodies,
            &[property()],
            &[Wrench::ZERO],
            Gravity::ZERO,
            &[contact],
            dt,
            controls(),
            admit_box,
            &CancelGate::new(),
        )
        .unwrap();
        // Independent unit-mass perfectly inelastic normal-impact solution.
        close_vec(
            bodies[0].linear_momentum_world(),
            Vec3::new(boost, 0.0, 0.0),
            2e-12,
        );
        close_vec(bodies[0].angular_momentum_body(), Vec3::ZERO, 2e-12);
        close_vec(
            bodies[0].pose().position_world(),
            initial
                .pose()
                .position_world()
                .add(Vec3::new(boost * dt, 0.0, 0.0)),
            1e-10,
        );
        close(result.impact.contact.impulses[0][0], normal_speed, 2e-12);
        close(
            result.impact.contact.dissipated_energy,
            0.5 * normal_speed * normal_speed,
            2e-12,
        );
        close(
            result.impact.kinetic_change_j,
            -0.5 * normal_speed * normal_speed,
            2e-12,
        );
        close(result.total_energy_residual_j, 0.0, 2e-12);
        close_vec(result.linear_momentum_residual_ns, Vec3::ZERO, 1e-12);
        close_vec(result.angular_momentum_residual_nms, Vec3::ZERO, 1e-9);
        assert!(result.path_admission > 0.0);
    }
}

#[test]
fn supported_gravity_kick_retains_its_work_under_a_large_common_boost() {
    let boost = 1e8;
    let (mut bodies, contact) = translated_scene(boost, Vec3::ZERO, 0.0);
    let dt = 1.0 / 32.0;
    let result = step_rigid_frictional_contacts(
        &mut bodies,
        &[property()],
        &[Wrench::ZERO],
        gravity(0.0),
        &[contact],
        dt,
        controls(),
        admit_box,
        &CancelGate::new(),
    )
    .unwrap();
    // Kick/contact/drift: gravity creates -g*dt normal momentum and the
    // support removes it before drift. Its splitting loss stays explicit.
    let impulse = 8.0 * dt;
    close_vec(
        bodies[0].linear_momentum_world(),
        Vec3::new(boost, 0.0, 0.0),
        2e-12,
    );
    close(result.impact.contact.impulses[0][0], impulse, 2e-12);
    close(result.force_impulse_work_j, 0.5 * impulse * impulse, 2e-12);
    close(
        result.impact.contact.dissipated_energy,
        0.5 * impulse * impulse,
        2e-12,
    );
    close(result.force_displacement_work_j, 0.0, 2e-12);
    close(result.total_energy_residual_j, 0.0, 2e-12);
    close(result.kick_work_energy_residual_j, 0.0, 2e-12);
}

#[test]
fn an_unrepresentable_native_force_kick_still_refuses_without_publishing() {
    let initial = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(1e16, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    let mut bodies = [initial];
    let result = step_rigid_frictional_contacts(
        &mut bodies,
        &[property()],
        &[Wrench {
            force_world: Vec3::new(1.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        }],
        Gravity::ZERO,
        &[],
        1.0,
        controls(),
        |_| Ok::<_, &'static str>(()),
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(RigidContactStepError::Audit(_))));
    assert_eq!(bodies, [initial]);
}
