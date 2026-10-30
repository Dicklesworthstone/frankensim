//! G1/G3: original inertia, finite rotation, torque work and spatial impulse.
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(norm(a.sub(b)), 0.0, tolerance);
}

#[test]
fn g1_spherical_cayley_is_explicit_and_default_exact_spin_is_unchanged() {
    let mass = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.3, 0.3, 0.3)).unwrap();
    let before =
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, 0.6)).unwrap();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let native = integrator.step(before, mass, Wrench::ZERO, 0.4).unwrap();
    let cayley = integrator
        .step_cayley(before, mass, Wrench::ZERO, 0.4)
        .unwrap();
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let ordinary = native
        .state_after
        .pose()
        .orientation()
        .rotate_body_to_world(axis);
    close_vec(
        ordinary,
        Vec3::new(0.8_f64.cos(), 0.8_f64.sin(), 0.0),
        4e-16,
    );
    let rotated = cayley
        .state_after
        .pose()
        .orientation()
        .rotate_body_to_world(axis);
    close_vec(
        rotated,
        Vec3::new((1.0 - 0.16) / 1.16, 0.8 / 1.16, 0.0),
        4e-16,
    );
    assert!(norm(ordinary.sub(rotated)) > 0.03);
    assert_eq!(
        native.state_after.angular_momentum_body(),
        cayley.state_after.angular_momentum_body()
    );
    close(
        cayley.diagnostics_after.mechanical_energy,
        cayley.diagnostics_before.mechanical_energy,
        1e-16,
    );
}

#[test]
fn g1_cayley_torque_work_and_actual_spatial_impulse_use_arithmetic_rotation_average() {
    let principal_frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.4), 0.6).unwrap();
    for moments in [Vec3::new(0.3, 0.3, 0.3), Vec3::new(0.2, 0.3, 0.4)] {
        let mass = MassProperties::from_principal_frame(1.7, Vec3::ZERO, moments, principal_frame)
            .unwrap();
        let attitude = UnitQuaternion::from_axis_angle(Vec3::new(-1.0, 0.4, 2.0), 1.1).unwrap();
        let before = RigidBodyState::new(
            Pose::new(Vec3::new(0.2, -0.1, 0.4), attitude).unwrap(),
            Vec3::new(0.1, -0.2, 0.3),
            Vec3::new(0.5, 0.3, -0.2),
        )
        .unwrap();
        let wrench = Wrench {
            force_world: Vec3::new(0.3, -0.2, 0.1),
            torque_body: Vec3::new(0.2, -0.15, 0.05),
        };
        let h = 0.16;
        let receipt = RigidBodyIntegrator::new(Gravity::ZERO)
            .step_cayley(before, mass, wrench, h)
            .unwrap();
        let after = receipt.state_after;
        let old_q = before.pose().orientation();
        let new_q = after.pose().orientation();
        let mid_l = before
            .angular_momentum_body()
            .scale(0.5)
            .add(after.angular_momentum_body().scale(0.5));
        let omega = mass.angular_velocity_body(mid_l);
        let c = old_q.rotate_body_to_world(omega.scale(h));
        let arm = Vec3::new(0.3, -0.1, 0.2);
        let r0 = old_q.rotate_body_to_world(arm);
        let r1 = new_q.rotate_body_to_world(arm);
        close_vec(r1.sub(r0), c.cross(r0.scale(0.5).add(r1.scale(0.5))), 4e-15);
        let angular_impulse = old_q
            .rotate_body_to_world(wrench.torque_body)
            .scale(0.5 * h)
            .add(
                new_q
                    .rotate_body_to_world(wrench.torque_body)
                    .scale(0.5 * h),
            );
        close_vec(
            new_q
                .rotate_body_to_world(after.angular_momentum_body())
                .sub(old_q.rotate_body_to_world(before.angular_momentum_body())),
            angular_impulse,
            2e-14,
        );
        let work = wrench.force_world.dot(
            after
                .pose()
                .position_world()
                .sub(before.pose().position_world()),
        ) + h * wrench.torque_body.dot(omega);
        close(
            receipt.diagnostics_after.mechanical_energy
                - receipt.diagnostics_before.mechanical_energy,
            work,
            3e-14,
        );
    }
}

#[test]
fn g3_full_inertia_free_rotation_preserves_native_invariants_over_many_steps() {
    let principal_frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
    let properties = MassProperties::from_principal_frame(
        1.0,
        Vec3::ZERO,
        Vec3::new(0.4, 0.6, 0.7),
        principal_frame,
    )
    .unwrap();
    let initial =
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.3, 0.5, 0.9)).unwrap();
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    assert_eq!(
        integrator
            .step(initial, properties, Wrench::ZERO, 0.02)
            .unwrap(),
        integrator
            .step_cayley(initial, properties, Wrench::ZERO, 0.02)
            .unwrap()
    );
    let energy = initial.kinetic_energy(properties).unwrap();
    let mut state = initial;
    for _ in 0..500 {
        state = integrator
            .step_cayley(state, properties, Wrench::ZERO, 0.02)
            .unwrap()
            .state_after;
        close(state.kinetic_energy(properties).unwrap(), energy, 2e-12);
        close_vec(
            state
                .pose()
                .orientation()
                .rotate_body_to_world(state.angular_momentum_body()),
            initial.angular_momentum_body(),
            2e-12,
        );
    }
}
