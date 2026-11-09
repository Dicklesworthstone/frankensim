//! Independent held-ODE, integrated-work and composition checks for native drag.
#[path = "linear_drag/affine.rs"]
mod affine;

use fs_mbd::linear_drag::{IsotropicLinearDrag, LinearDragError};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn close_vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}
fn properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(3.0, 3.0, 3.0)).unwrap()
}
fn state(properties: MassProperties, velocity: Vec3, omega: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::identity(),
        velocity.scale(properties.mass()),
        properties.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap()
}
fn force(force_world: Vec3) -> Wrench {
    Wrench {
        force_world,
        torque_body: Vec3::ZERO,
    }
}

#[test]
fn held_force_velocity_displacement_spin_and_heat_match_independent_odes() {
    let p = properties();
    let ambient = Vec3::new(0.5, -1.0, 2.0);
    let gravity = Gravity::new(Vec3::new(0.1, 0.0, -0.2)).unwrap();
    let applied = force(Vec3::new(0.3, 0.4, 0.5));
    let velocity = Vec3::new(1.5, 0.4, -0.3);
    let omega = Vec3::new(0.3, -0.7, 1.1);
    let initial = state(p, velocity, omega);
    let duration = 1.7_f64;
    let translation_rate = 2.0_f64;
    let spin_rate = 2.0_f64;
    let e = (-translation_rate * duration).exp();
    let terminal_slip = applied
        .force_world
        .add(gravity.acceleration_world().scale(p.mass()))
        .scale(0.25);
    let deviation = velocity.sub(ambient).sub(terminal_slip);
    let expected_velocity = ambient.add(terminal_slip).add(deviation.scale(e));
    let integral = -(-translation_rate * duration).exp_m1() / translation_rate;
    let expected_displacement = ambient
        .add(terminal_slip)
        .scale(duration)
        .add(deviation.scale(integral));
    let heat = 4.0
        * (duration * terminal_slip.norm_squared()
            + 2.0 * integral * terminal_slip.dot(deviation)
            + deviation.norm_squared() * -(-2.0 * translation_rate * duration).exp_m1()
                / (2.0 * translation_rate));
    let spin_heat = 0.5 * 3.0 * omega.norm_squared() * -(-2.0 * spin_rate * duration).exp_m1();
    let step = IsotropicLinearDrag::new(p, gravity, ambient, 4.0, 6.0)
        .unwrap()
        .step(initial, applied, duration)
        .unwrap();
    close_vector(
        step.state_after.center_of_mass_velocity_world(p).unwrap(),
        expected_velocity,
        3e-14,
    );
    close_vector(step.displacement_world, expected_displacement, 3e-14);
    close_vector(
        step.state_after.pose().position_world(),
        expected_displacement,
        3e-14,
    );
    close_vector(
        p.angular_velocity_body_checked(step.state_after.angular_momentum_body())
            .unwrap(),
        omega.scale((-spin_rate * duration).exp()),
        3e-14,
    );
    close_vector(step.rotation_vector_body, omega.scale(integral), 3e-14);
    let expected_orientation =
        UnitQuaternion::from_axis_angle(omega, omega.norm_squared().sqrt() * integral).unwrap();
    close_vector(
        step.state_after
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        expected_orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        3e-14,
    );
    close(step.translation_heat_j, heat, 3e-14);
    close(step.rotation_heat_j, spin_heat, 3e-14);
    close(
        step.background_work_j,
        step.drag_impulse_world.dot(ambient),
        3e-14,
    );
    close(
        step.applied_work_j,
        applied.force_world.dot(expected_displacement),
        3e-14,
    );
    close(step.work_energy_residual_j, 0.0, 8e-14);
    close_vector(step.linear_momentum_residual_world, Vec3::ZERO, 3e-14);
    close_vector(step.angular_momentum_residual_world, Vec3::ZERO, 3e-14);
}

#[test]
fn stiff_entrainment_does_not_reverse_slip_and_accounts_for_bath_work() {
    let p = properties();
    let ambient = Vec3::new(3.0, -1.0, 0.5);
    let step = IsotropicLinearDrag::new(p, Gravity::ZERO, ambient, 4.0, 6.0)
        .unwrap()
        .step(state(p, Vec3::ZERO, Vec3::ZERO), Wrench::ZERO, 25.0)
        .unwrap();
    close_vector(
        step.state_after.center_of_mass_velocity_world(p).unwrap(),
        ambient,
        2e-14,
    );
    close_vector(step.displacement_world, ambient.scale(24.5), 2e-14);
    close(
        step.translation_heat_j,
        0.5 * p.mass() * ambient.norm_squared(),
        3e-14,
    );
    close(
        step.background_work_j,
        p.mass() * ambient.norm_squared(),
        3e-14,
    );
    close(step.work_energy_residual_j, 0.0, 3e-13);
}

#[test]
fn exact_path_heat_and_work_compose_across_unequal_intervals() {
    let p = properties();
    let flow = IsotropicLinearDrag::new(
        p,
        Gravity::new(Vec3::new(0.0, 0.0, -0.3)).unwrap(),
        Vec3::new(0.2, -0.1, 0.5),
        4.0,
        6.0,
    )
    .unwrap();
    let initial = state(p, Vec3::new(1.2, -0.6, 0.7), Vec3::new(0.4, 0.2, -0.3));
    let applied = force(Vec3::new(-0.1, 0.3, 0.7));
    let whole = flow.step(initial, applied, 20.0).unwrap();
    let mut current = initial;
    let mut heat = 0.0;
    let mut bath_work = 0.0;
    let mut drive_work = 0.0;
    for interval in [0.08, 0.04, 1.0, 18.88] {
        let part = flow.step(current, applied, interval).unwrap();
        current = part.state_after;
        heat += part.dissipated_energy_j;
        bath_work += part.background_work_j;
        drive_work += part.applied_work_j;
    }
    close_vector(
        current.linear_momentum_world(),
        whole.state_after.linear_momentum_world(),
        8e-14,
    );
    close_vector(
        current.pose().position_world(),
        whole.state_after.pose().position_world(),
        8e-14,
    );
    close_vector(
        current
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0)),
        whole
            .state_after
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0)),
        8e-14,
    );
    close(heat, whole.dissipated_energy_j, 8e-14);
    close(bath_work, whole.background_work_j, 8e-14);
    close(drive_work, whole.applied_work_j, 8e-14);
}

#[test]
fn exact_terminal_settling_persists_for_short_and_long_intervals() {
    let p = properties();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -2.0)).unwrap();
    let flow = IsotropicLinearDrag::new(p, gravity, Vec3::ZERO, 4.0, 6.0).unwrap();
    let initial = state(p, Vec3::new(0.0, 0.0, -1.0), Vec3::ZERO);
    for duration in [1e-12, 0.124, 0.126, 1.0, 1e6] {
        let step = flow.step(initial, Wrench::ZERO, duration).unwrap();
        close_vector(
            step.state_after.linear_momentum_world(),
            initial.linear_momentum_world(),
            2e-14,
        );
        close(step.displacement_world.z, -duration, 2e-14);
        close(step.translation_heat_j, 4.0 * duration, 2e-14);
        close(step.work_energy_residual_j / duration, 0.0, 2e-14);
    }
}

#[test]
fn tiny_unrepresented_momentum_loss_remains_physical_heat_and_a_visible_defect() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let initial = state(p, Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO);
    let step = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 1.0, 0.0)
        .unwrap()
        .step(initial, Wrench::ZERO, 1e-20)
        .unwrap();
    assert_eq!(
        step.state_after.linear_momentum_world(),
        initial.linear_momentum_world()
    );
    assert!(step.translation_heat_j > 0.0);
    close(step.translation_heat_j / 1e-20, 1.0, 2e-14);
    close(
        step.work_energy_residual_j / step.translation_heat_j,
        1.0,
        2e-14,
    );
    close(step.linear_momentum_residual_world.x / 1e-20, 1.0, 2e-14);
}

#[test]
fn zero_drag_is_exact_ballistic_motion_independent_of_the_bath() {
    let p = properties();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -2.0)).unwrap();
    let velocity = Vec3::new(0.3, 0.4, 0.5);
    let omega = Vec3::new(0.0, 0.0, 0.7);
    let initial = state(p, velocity, omega);
    let applied = force(Vec3::new(0.2, -0.4, 0.6));
    let a = gravity
        .acceleration_world()
        .add(applied.force_world.scale(0.5));
    for ambient in [Vec3::ZERO, Vec3::new(1e100, -1e200, 1e300)] {
        let step = IsotropicLinearDrag::new(p, gravity, ambient, 0.0, 0.0)
            .unwrap()
            .step(initial, applied, 2.0)
            .unwrap();
        close_vector(
            step.state_after.center_of_mass_velocity_world(p).unwrap(),
            velocity.add(a.scale(2.0)),
            3e-14,
        );
        close_vector(
            step.displacement_world,
            velocity.scale(2.0).add(a.scale(2.0)),
            3e-14,
        );
        assert_eq!(
            step.state_after.angular_momentum_body(),
            initial.angular_momentum_body()
        );
        assert_eq!(step.dissipated_energy_j, 0.0);
        assert_eq!(step.background_work_j, 0.0);
        close_vector(step.rotation_vector_body, omega.scale(2.0), 3e-14);
        close(step.work_energy_residual_j, 0.0, 3e-14);
    }
}

#[test]
fn galilean_boost_preserves_slip_heat_and_changes_bath_work_by_impulse_pairing() {
    let p = properties();
    let boost = Vec3::new(512.0, -1024.0, 256.0);
    let velocity = Vec3::new(0.5, -0.25, 0.75);
    let ambient = Vec3::new(-0.25, 0.5, 0.125);
    let applied = force(Vec3::new(0.2, -0.4, 0.3));
    let run = |offset: Vec3| {
        IsotropicLinearDrag::new(p, Gravity::ZERO, ambient.add(offset), 4.0, 6.0)
            .unwrap()
            .step(state(p, velocity.add(offset), Vec3::ZERO), applied, 0.12)
            .unwrap()
    };
    let original = run(Vec3::ZERO);
    let boosted = run(boost);
    close_vector(
        boosted
            .state_after
            .center_of_mass_velocity_world(p)
            .unwrap()
            .sub(boost),
        original
            .state_after
            .center_of_mass_velocity_world(p)
            .unwrap(),
        2e-13,
    );
    close_vector(
        boosted.displacement_world.sub(boost.scale(0.12)),
        original.displacement_world,
        2e-13,
    );
    close(
        boosted.translation_heat_j,
        original.translation_heat_j,
        2e-14,
    );
    close(
        boosted.background_work_j - original.background_work_j,
        original.drag_impulse_world.dot(boost),
        3e-14,
    );
}

#[test]
fn nonspherical_inertia_applied_torque_and_invalid_controls_refuse() {
    let p = properties();
    let asymmetric = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    assert!(matches!(
        IsotropicLinearDrag::new(asymmetric, Gravity::ZERO, Vec3::ZERO, 1.0, 1.0),
        Err(LinearDragError::Invalid(_))
    ));
    for resistance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, resistance, 1.0).is_err());
        assert!(IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 1.0, resistance).is_err());
    }
    let flow = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 1.0, 1.0).unwrap();
    let initial = state(p, Vec3::ZERO, Vec3::ZERO);
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(flow.step(initial, Wrench::ZERO, duration).is_err());
    }
    assert!(matches!(
        flow.step(
            initial,
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: Vec3::new(0.0, 0.0, 1.0)
            },
            1.0
        ),
        Err(LinearDragError::Invalid(_))
    ));
}
