use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn body(mass: f64, inertia: Vec3, linear: Vec3, angular: Vec3) -> (MassProperties, RigidBodyState) {
    (
        MassProperties::new(mass, Vec3::ZERO, inertia).unwrap(),
        RigidBodyState::new(Pose::identity(), linear, angular).unwrap(),
    )
}

#[test]
fn midpoint_forcing_and_gravity_close_the_original_work_balance() {
    let (properties, state) = body(
        2.0,
        Vec3::new(2.0, 2.0, 3.0),
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 6.0),
    );
    let integrator = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -10.0)).unwrap());
    let wrench = Wrench {
        force_world: Vec3::new(4.0, 0.0, 5.0),
        torque_body: Vec3::new(0.0, 0.0, 3.0),
    };
    let receipt = integrator.step(state, properties, wrench, 0.5).unwrap();
    // Independent constant-acceleration translation and single-axis rotation.
    let force_work = 4.0 * 0.75 + 5.0 * -0.9375;
    let torque_work = 3.0 * 1.125;
    let gravity_work = -20.0 * -0.9375;
    near(receipt.work.applied_force_work_j, force_work, 1e-14);
    near(receipt.work.applied_torque_work_j, torque_work, 1e-14);
    near(receipt.work.gravity_work_j, gravity_work, 1e-14);
    near(
        receipt.work.gravitational_potential_energy_change_j,
        -gravity_work,
        1e-14,
    );
    near(receipt.work.applied_work_j, force_work + torque_work, 1e-14);
    near(
        receipt.work.kinetic_energy_change_j,
        force_work + torque_work + gravity_work,
        1e-14,
    );
    near(
        receipt.work.mechanical_energy_change_j,
        force_work + torque_work,
        1e-14,
    );
    near(receipt.work.work_energy_residual_j, 0.0, 2e-14);
}

#[test]
fn braking_and_regeneration_keep_negative_work() {
    let (properties, state) = body(
        1.0,
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(4.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 3.0),
    );
    let wrench = Wrench {
        force_world: Vec3::new(-2.0, 0.0, 0.0),
        torque_body: Vec3::new(0.0, 0.0, -1.0),
    };
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    for receipt in [
        integrator.step(state, properties, wrench, 0.5).unwrap(),
        integrator
            .step_symplectic(state, properties, wrench, 0.5)
            .unwrap(),
    ] {
        near(receipt.work.applied_force_work_j, -3.5, 1e-14);
        near(receipt.work.applied_torque_work_j, -1.375, 1e-14);
        near(receipt.work.applied_work_j, -4.875, 1e-14);
        near(receipt.work.kinetic_energy_change_j, -4.875, 1e-14);
        near(receipt.work.work_energy_residual_j, 0.0, 1e-14);
    }
}

#[test]
fn torque_work_is_covariant_with_full_body_frame_inertia() {
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap();
    let moments = Vec3::new(2.0, 3.0, 4.0);
    let principal = MassProperties::new(1.0, Vec3::ZERO, moments).unwrap();
    let full = MassProperties::from_principal_frame(1.0, Vec3::ZERO, moments, q).unwrap();
    let momentum = Vec3::new(1.7, -2.1, 3.4);
    let torque = Vec3::new(-0.4, 0.9, 0.3);
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let first = integrator
        .step(
            RigidBodyState::new(Pose::identity(), Vec3::ZERO, momentum).unwrap(),
            principal,
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: torque,
            },
            0.02,
        )
        .unwrap();
    let second = integrator
        .step(
            RigidBodyState::new(
                Pose::identity(),
                Vec3::ZERO,
                q.rotate_body_to_world(momentum),
            )
            .unwrap(),
            full,
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: q.rotate_body_to_world(torque),
            },
            0.02,
        )
        .unwrap();
    near(
        second.work.applied_torque_work_j,
        first.work.applied_torque_work_j,
        1e-14,
    );
    near(
        second.work.kinetic_energy_change_j,
        first.work.kinetic_energy_change_j,
        3e-14,
    );
    near(second.work.work_energy_residual_j, 0.0, 3e-14);
    let cayley = integrator
        .step_cayley(
            first.state_before,
            principal,
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: torque,
            },
            0.02,
        )
        .unwrap();
    near(
        cayley.work.applied_torque_work_j,
        first.work.applied_torque_work_j,
        1e-14,
    );
}

#[test]
fn verlet_uses_both_original_endpoint_torque_kicks() {
    let (properties, state) = body(
        2.0,
        Vec3::new(2.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(1.5, -2.0, 3.0),
    );
    let first = Wrench {
        force_world: Vec3::new(2.0, 0.0, 0.0),
        torque_body: Vec3::new(0.3, -0.4, 0.5),
    };
    let last = Wrench {
        force_world: Vec3::new(-1.0, 0.0, 0.0),
        torque_body: Vec3::new(-0.2, 0.7, 0.4),
    };
    let dt = 0.3;
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let stage = integrator
        .begin_verlet(state, properties, first, dt)
        .unwrap();
    let drifted = stage.drifted_state();
    let initial_l = state.angular_momentum_body();
    let kicked_l = initial_l.add(first.torque_body.scale(0.5 * dt));
    let receipt = stage.finish(last).unwrap();
    let energy = |p: Vec3| 0.5 * p.dot(properties.angular_velocity_body(p));
    // Independent energy differences across each bounded torque kick only.
    let expected = energy(kicked_l) - energy(initial_l)
        + energy(receipt.state_after.angular_momentum_body())
        - energy(drifted.angular_momentum_body());
    near(receipt.work.applied_torque_work_j, expected, 3e-14);
    let dx = receipt
        .state_after
        .pose()
        .position_world()
        .sub(state.pose().position_world());
    near(
        receipt.work.applied_force_work_j,
        first
            .force_world
            .scale(0.5)
            .add(last.force_world.scale(0.5))
            .dot(dx),
        1e-14,
    );
    assert!(
        (receipt.work.applied_torque_work_j
            - integrator
                .step_symplectic(state, properties, first, dt)
                .unwrap()
                .work
                .applied_torque_work_j)
            .abs()
            > 0.01
    );
}

#[test]
fn free_rotation_error_is_not_invented_source_work_or_heat() {
    let (properties, state) = body(
        1.0,
        Vec3::new(2.0, 3.0, 4.0),
        Vec3::ZERO,
        Vec3::new(2.0, 3.0, 4.0),
    );
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let receipt = integrator
        .step_symplectic(state, properties, Wrench::ZERO, 0.3)
        .unwrap();
    assert_eq!(receipt.work.applied_work_j, 0.0);
    assert_eq!(receipt.work.gravity_work_j, 0.0);
    assert!(receipt.work.work_energy_residual_j.abs() > 1e-8);
    assert_eq!(
        receipt.work.work_energy_residual_j,
        receipt.work.kinetic_energy_change_j
    );
    let refined = integrator
        .step_symplectic(state, properties, Wrench::ZERO, 0.15)
        .unwrap();
    assert!(
        refined.work.work_energy_residual_j.abs() < receipt.work.work_energy_residual_j.abs() / 4.0
    );
    let fourth = integrator
        .step_free_fourth_order(state, properties, 0.3)
        .unwrap();
    assert_eq!(fourth.work.applied_work_j, 0.0);
    assert!(fourth.work.work_energy_residual_j.abs() < receipt.work.work_energy_residual_j.abs());
}

#[test]
fn kick_drift_reports_its_force_displacement_energy_defect() {
    let (properties, state) = body(2.0, Vec3::new(1.0, 1.0, 1.0), Vec3::ZERO, Vec3::ZERO);
    let receipt = RigidBodyIntegrator::new(Gravity::ZERO)
        .step_kick_drift(
            state,
            properties,
            Wrench {
                force_world: Vec3::new(4.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            0.5,
        )
        .unwrap();
    assert_eq!(receipt.work.applied_force_work_j, 2.0);
    assert_eq!(receipt.work.kinetic_energy_change_j, 1.0);
    assert_eq!(receipt.work.work_energy_residual_j, -1.0);
}

#[test]
fn small_work_survives_unchanged_dominant_stored_spin() {
    let (properties, state) = body(
        1.0,
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 1e12),
    );
    let receipt = RigidBodyIntegrator::new(Gravity::ZERO)
        .step(
            state,
            properties,
            Wrench {
                force_world: Vec3::new(1.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            1.0,
        )
        .unwrap();
    assert_eq!(
        receipt.diagnostics_after.mechanical_energy - receipt.diagnostics_before.mechanical_energy,
        0.0
    );
    assert_eq!(receipt.work.applied_work_j, 0.5);
    assert_eq!(receipt.work.kinetic_energy_change_j, 0.5);
    assert_eq!(receipt.work.work_energy_residual_j, 0.0);
}

#[test]
fn unrepresented_momentum_increment_is_a_visible_storage_defect() {
    let large = 2.0_f64.powi(60);
    let (properties, state) = body(
        large,
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(large, 0.0, 0.0),
        Vec3::ZERO,
    );
    let receipt = RigidBodyIntegrator::new(Gravity::ZERO)
        .step(
            state,
            properties,
            Wrench {
                force_world: Vec3::new(1.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            1.0,
        )
        .unwrap();
    assert_eq!(
        receipt.state_before.linear_momentum_world(),
        receipt.state_after.linear_momentum_world()
    );
    assert_eq!(receipt.work.kinetic_energy_change_j, 0.0);
    near(receipt.work.applied_work_j, 1.0, 1e-15);
    near(receipt.work.work_energy_residual_j, -1.0, 1e-15);
}

#[test]
fn torque_work_survives_an_underflowing_impulse_or_power() {
    // Both cases turn by 1e-10 rad and receive about 1e-230 J. In the first,
    // torque*dt underflows; in the second, torque*omega underflows. Choosing
    // either association unconditionally would hide a real source transfer.
    for (inertia, angular_momentum, dt) in [(1e-200, 1e-100, 1e-110), (1e100, 1e-10, 1e100)] {
        let torque = 1e-220;
        let (properties, state) = body(
            1.0,
            Vec3::new(inertia, inertia, inertia),
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, angular_momentum),
        );
        let wrench = Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(0.0, 0.0, torque),
        };
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let angular_rate = angular_momentum / inertia;
        let expected = torque * (angular_rate * dt);
        assert!(expected > 0.0);
        assert!(torque * dt == 0.0 || torque * angular_rate == 0.0);
        for receipt in [
            integrator.step(state, properties, wrench, dt).unwrap(),
            integrator
                .step_symplectic(state, properties, wrench, dt)
                .unwrap(),
            integrator
                .step_kick_drift(state, properties, wrench, dt)
                .unwrap(),
        ] {
            near(
                receipt.work.applied_torque_work_j,
                expected,
                1e-13 * expected,
            );
            assert_eq!(receipt.work.kinetic_energy_change_j, 0.0);
            near(
                receipt.work.work_energy_residual_j,
                -expected,
                1e-13 * expected,
            );
        }
    }
}

#[test]
fn unequal_subnormal_endpoint_forces_retain_their_displacement_work() {
    let unit = f64::from_bits(1);
    let (properties, state) = body(
        1.0,
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(1e150, 0.0, 0.0),
        Vec3::ZERO,
    );
    let first = Wrench {
        force_world: Vec3::new(unit, 0.0, 0.0),
        torque_body: Vec3::ZERO,
    };
    let last = Wrench {
        force_world: Vec3::new(2.0 * unit, 0.0, 0.0),
        torque_body: Vec3::ZERO,
    };
    let receipt = RigidBodyIntegrator::new(Gravity::ZERO)
        .begin_verlet(state, properties, first, 1e150)
        .unwrap()
        .finish(last)
        .unwrap();
    let displacement = receipt.state_after.pose().position_world().x;
    // Contract before halving; a rounded force average loses part of this work.
    let expected = 0.5 * (unit * displacement) + 0.5 * ((2.0 * unit) * displacement);
    assert!(expected > 0.0);
    near(
        receipt.work.applied_force_work_j,
        expected,
        1e-14 * expected,
    );
    assert_eq!(receipt.work.kinetic_energy_change_j, 0.0);
    near(
        receipt.work.work_energy_residual_j,
        -expected,
        1e-14 * expected,
    );
}
